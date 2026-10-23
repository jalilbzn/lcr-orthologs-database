import json
import logging
import math
import sys
import time
import numpy as np
import psycopg2
import psycopg2.extras
import psycopg2.extensions
from Bio.Align import PairwiseAligner, substitution_matrices
from ortholog_lcr_pipeline.inputs.pipeline.config import DB_CONFIG, LOGS_DIR

# Register numpy types with psycopg2
psycopg2.extensions.register_adapter(np.int64, lambda val: psycopg2.extensions.AsIs(int(val)))
psycopg2.extensions.register_adapter(np.int32, lambda val: psycopg2.extensions.AsIs(int(val)))
psycopg2.extensions.register_adapter(np.float64, lambda val: psycopg2.extensions.AsIs(float(val)))
psycopg2.extensions.register_adapter(np.float32, lambda val: psycopg2.extensions.AsIs(float(val)))

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
    handlers=[
        logging.FileHandler(LOGS_DIR / "anchor_map.log", mode="w"),
        logging.StreamHandler(sys.stdout),
    ],
)
log = logging.getLogger("anchor_map")

FLANK_LEN = 100
MIN_ANCHOR_COV = 0.50
MIN_INTERVAL_LEN = 1    # skip if candidate interval < 1 aa
MAX_INTERVAL_LEN = 500  # skip if candidate interval > 500 aa (likely misalignment)

TAXONOMY = {
    "Chimp": "Mammalia", "Mouse": "Mammalia", "Rat": "Mammalia",
    "Macaque": "Mammalia", "Cow": "Mammalia", "Pig": "Mammalia",
    "Dog": "Mammalia", "Opossum": "Mammalia",
    "Chicken": "Vertebrata", "Xenopus": "Vertebrata", "Zebrafish": "Vertebrata",
}

# ── Biophysics constants ──
KD_SCALE = {
    'A':  1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C':  2.5,
    'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I':  4.5,
    'L':  3.8, 'K': -3.9, 'M':  1.9, 'F':  2.8, 'P': -1.6,
    'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V':  4.2,
}
STANDARD_AA = set(KD_SCALE.keys())

# ── Smith-Waterman aligners ──
aligner_62 = PairwiseAligner()
aligner_62.mode = 'local'
aligner_62.substitution_matrix = substitution_matrices.load("BLOSUM62")
aligner_62.open_gap_score = -11
aligner_62.extend_gap_score = -1

aligner_45 = PairwiseAligner()
aligner_45.mode = 'local'
aligner_45.substitution_matrix = substitution_matrices.load("BLOSUM45")
aligner_45.open_gap_score = -11
aligner_45.extend_gap_score = -1


def sw_align(query, target, current_aligner):
    """
    Local-align query to target.
    Returns (score, t_start, t_end, identity, coverage) or None.
    t_start/t_end are 0-based.
    """
    if len(query) < 10 or len(target) < 10:
        return None
    try:
        alignments = current_aligner.align(query, target)
        if not alignments:
            return None
        top = alignments[0]
        aligned = top.aligned
        if len(aligned) != 2 or len(aligned[0]) == 0:
            return None

        t_blocks = aligned[1]
        t_start = t_blocks[0][0]
        t_end = t_blocks[-1][1] - 1

        q_aligned_len = sum(b[1] - b[0] for b in aligned[0])
        coverage = q_aligned_len / len(query)
        identity = min(1.0, top.score / (len(query) * 4.0))

        return (float(top.score), int(t_start), int(t_end), round(float(identity), 3), round(float(coverage), 3))
    except Exception:
        return None


def _py(val):
    """Convert numpy/non-standard types to native Python."""
    if val is None:
        return None
    if isinstance(val, (np.integer,)):
        return int(val)
    if isinstance(val, (np.floating,)):
        return float(val)
    return val


def sanitize_row(row):
    """Sanitize a tuple/list of values for psycopg2."""
    return tuple(_py(v) for v in row)


def compute_aa_composition(seq):
    n = len(seq)
    if n == 0:
        return {}
    counts = {}
    for aa in seq.upper():
        if aa in STANDARD_AA:
            counts[aa] = counts.get(aa, 0) + 1
    return {aa: round(c / n, 4) for aa, c in sorted(counts.items())}


def compute_biophys(seq):
    seq = seq.upper()
    n = len(seq)
    if n == 0:
        return None
    n_pos = sum(1 for aa in seq if aa in ('K', 'R'))
    n_neg = sum(1 for aa in seq if aa in ('D', 'E'))
    n_aro = sum(1 for aa in seq if aa in ('F', 'W', 'Y'))
    n_pro = sum(1 for aa in seq if aa == 'P')
    fcr = (n_pos + n_neg) / n
    ncpr = (n_pos - n_neg) / n
    hydro = sum(KD_SCALE.get(aa, 0.0) for aa in seq) / n
    # Kappa (blob method)
    kappa = None
    n_charged = n_pos + n_neg
    if n_charged >= 4 and n >= 6:
        blob_size = min(5, n)
        n_blobs = n - blob_size + 1
        sigma = 0.0
        for i in range(n_blobs):
            blob = seq[i:i + blob_size]
            bp = sum(1 for aa in blob if aa in ('K', 'R'))
            bn = sum(1 for aa in blob if aa in ('D', 'E'))
            sigma += ((bp - bn) / blob_size) ** 2
        sigma /= n_blobs
        sigma_max = fcr if fcr > 0 else 1
        kappa = round(min(sigma / sigma_max, 1.0), 4)
    return {
        'fcr': round(fcr, 4), 'ncpr': round(ncpr, 4), 'kappa': kappa,
        'mean_hydropathy': round(hydro, 4),
        'aromatic_frac': round(n_aro / n, 4), 'proline_frac': round(n_pro / n, 4),
    }


# ═══════════════════════════════════════════════════════════
# Database helpers
# ═══════════════════════════════════════════════════════════

def load_families_with_data(conn):
    """Load families + human protein + human LCRs + ortholog sequences in bulk."""
    with conn.cursor() as cur:
        # Families with members
        cur.execute("""
            SELECT f.family_id, f.source_gene, f.members,
                   p.protein_id, p.sequence
            FROM lcr_family f
            JOIN protein p ON p.uniprot_acc = f.source_gene
            WHERE f.members != '{}'::jsonb
            ORDER BY f.family_id
        """)
        families = cur.fetchall()
    return families


def load_human_lcrs_for_protein(conn, protein_id):
    with conn.cursor() as cur:
        cur.execute("""
            SELECT lcr_id, start_pos, end_pos, sequence
            FROM lcr_seq WHERE protein_id = %s ORDER BY start_pos
        """, (protein_id,))
        return cur.fetchall()


def get_protein_seq(conn, pid):
    with conn.cursor() as cur:
        cur.execute("SELECT sequence FROM protein WHERE protein_id = %s", (pid,))
        r = cur.fetchone()
        return r[0] if r else None


def family_already_done(conn, family_id):
    with conn.cursor() as cur:
        cur.execute("SELECT 1 FROM phylogeny_map WHERE family_id=%s LIMIT 1", (family_id,))
        return cur.fetchone() is not None


# ═══════════════════════════════════════════════════════════
# Core logic
# ═══════════════════════════════════════════════════════════

def process_family(conn, fam_id, members, human_pid, human_seq):
    """
    For each human LCR, align flanks to each ortholog protein.
    The candidate interval between anchors becomes the ortholog LCR.
    """
    if isinstance(members, str):
        members = json.loads(members)
    if not members:
        return 0

    human_lcrs = load_human_lcrs_for_protein(conn, human_pid)
    if not human_lcrs:
        return 0

    new_lcr_rows = []     # (protein_id, start, end, seq, comp_json, seg_mode)
    new_bio_rows = []     # biophys dicts
    new_phylo_rows = []   # phylogeny_map tuples (deferred until lcr_ids known)
    phylo_meta = []       # metadata for building phylo rows after lcr insert

    for species_name, orth_pid in members.items():
        orth_seq = get_protein_seq(conn, orth_pid)
        if not orth_seq:
            continue

        tax_level = TAXONOMY.get(species_name, "Vertebrata")
        current_aligner = aligner_62 if tax_level == "Mammalia" else aligner_45

        for lcr_id, lcr_start, lcr_end, lcr_seq_str in human_lcrs:
            # Extract flanks (0-based slicing, positions are 1-based)
            left_start = max(0, lcr_start - 1 - FLANK_LEN)
            left_flank = human_seq[left_start:lcr_start - 1]

            right_end_idx = min(len(human_seq), lcr_end + FLANK_LEN)
            right_flank = human_seq[lcr_end:right_end_idx]

            if len(left_flank) < 10 and len(right_flank) < 10:
                continue

            # Align flanks
            left_hit = sw_align(left_flank, orth_seq, current_aligner)
            right_hit = sw_align(right_flank, orth_seq, current_aligner)

            left_ok = left_hit and left_hit[4] >= MIN_ANCHOR_COV
            right_ok = right_hit and right_hit[4] >= MIN_ANCHOR_COV

            if not left_ok and not right_ok:
                continue

            # Define candidate interval
            if left_ok and right_ok:
                if left_hit[2] >= right_hit[1]:
                    continue  # anchors crossed
                cand_start = left_hit[2] + 1   # 0-based, after left anchor
                cand_end = right_hit[1] - 1     # 0-based, before right anchor
            elif left_ok:
                cand_start = left_hit[2] + 1
                cand_end = min(len(orth_seq) - 1,
                               cand_start + (lcr_end - lcr_start + 1) * 2)
            else:
                cand_end = right_hit[1] - 1
                cand_start = max(0, cand_end - (lcr_end - lcr_start + 1) * 2)

            if cand_start > cand_end:
                continue

            interval_len = cand_end - cand_start + 1
            if interval_len < MIN_INTERVAL_LEN or interval_len > MAX_INTERVAL_LEN:
                continue

            # The candidate interval IS the ortholog LCR
            orth_lcr_seq = orth_seq[cand_start:cand_end + 1]
            comp = compute_aa_composition(orth_lcr_seq)
            bio = compute_biophys(orth_lcr_seq)

            # 1-based positions for DB
            db_start = cand_start + 1
            db_end = cand_end + 1

            new_lcr_rows.append((
                orth_pid, db_start, db_end, orth_lcr_seq,
                json.dumps(comp), 'anchor'
            ))
            new_bio_rows.append(bio)

            # Synteny scores
            synteny_conf = min(
                left_hit[4] if left_ok else 1.0,
                right_hit[4] if right_ok else 1.0
            )
            if not left_ok or not right_ok:
                synteny_conf = left_hit[4] if left_ok else right_hit[4]

            anchor_score = 0
            n_anch = 0
            if left_ok:
                anchor_score += left_hit[0]; n_anch += 1
            if right_ok:
                anchor_score += right_hit[0]; n_anch += 1
            anchor_score /= max(n_anch, 1)

            phylo_meta.append({
                'source_lcr_id': lcr_id,    # Store the Human LCR ID
                'fam_id': fam_id,
                'tax': tax_level,
                'left': left_hit if left_ok else None,
                'right': right_hit if right_ok else None,
                'orth_pid': orth_pid,
                'synteny': round(synteny_conf, 2),
                'anchor': round(anchor_score, 2),
            })

    if not new_lcr_rows:
        return 0

    # Insert LCRs and get IDs
    with conn.cursor() as cur:
        sql_lcr = """
            INSERT INTO lcr_seq (protein_id, start_pos, end_pos, sequence,
                                 aa_composition_json, seg_mode)
            VALUES %s RETURNING lcr_id
        """
        result = psycopg2.extras.execute_values(
            cur, sql_lcr, [sanitize_row(r) for r in new_lcr_rows],
            page_size=200, fetch=True
        )
        new_ids = [r[0] for r in result]

        # Insert biophys
        bio_vals = []
        for i, lcr_id in enumerate(new_ids):
            b = new_bio_rows[i]
            if b:
                bio_vals.append((
                    lcr_id, b['fcr'], b['ncpr'], b['kappa'],
                    b['mean_hydropathy'], b['aromatic_frac'], b['proline_frac']
                ))
        if bio_vals:
            sql_bio = """
                INSERT INTO lcr_biophys (lcr_id, fcr, ncpr, kappa,
                                         mean_hydropathy, aromatic_frac, proline_frac)
                VALUES %s ON CONFLICT (lcr_id) DO NOTHING
            """
            psycopg2.extras.execute_values(cur, sql_bio, [sanitize_row(r) for r in bio_vals], page_size=200)

        # Insert phylogeny_map
        phylo_rows = []
        for i, lcr_id in enumerate(new_ids):
            m = phylo_meta[i]
            l = m['left']
            r = m['right']
            phylo_rows.append((
                lcr_id, m['fam_id'], m['tax'],
                m['orth_pid'] if l else None,
                l[1] + 1 if l else None, l[2] + 1 if l else None,
                l[3] if l else None, l[4] if l else None, None,
                m['orth_pid'] if r else None,
                r[1] + 1 if r else None, r[2] + 1 if r else None,
                r[3] if r else None, r[4] if r else None, None,
                m['synteny'], m['anchor'], m['source_lcr_id'],
            ))

        sql_phylo = """
            INSERT INTO phylogeny_map (
                lcr_id, family_id, taxonomy_level,
                left_protein_id, left_start, left_end, left_pid, left_cov, left_evalue,
                right_protein_id, right_start, right_end, right_pid, right_cov, right_evalue,
                synteny_conf, anchor_score, source_lcr_id
            ) VALUES %s ON CONFLICT (lcr_id, family_id) DO NOTHING
        """
        psycopg2.extras.execute_values(cur, sql_phylo, [sanitize_row(r) for r in phylo_rows], page_size=200)

    conn.commit()
    return len(new_ids)


# ═══════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════

def main():
    t_start = time.time()
    conn = psycopg2.connect(**DB_CONFIG)

    families = load_families_with_data(conn)
    total = len(families)
    log.info(f"Processing {total} families with members...")

    total_mappings = 0
    skipped = 0

    for i, (fam_id, src_gene, members, hpid, hseq) in enumerate(families):
        if family_already_done(conn, fam_id):
            skipped += 1
            continue

        try:
            n = process_family(conn, fam_id, members, hpid, hseq)
            total_mappings += n
        except Exception as e:
            import traceback
            log.warning(f"Family {fam_id} error: {e}\n{traceback.format_exc()}")
            conn.rollback()

        if (i + 1) % 200 == 0:
            elapsed = time.time() - t_start
            log.info(f"  {i+1}/{total} families, {total_mappings} mappings "
                     f"(skipped {skipped}) in {elapsed:.0f}s")

    elapsed = time.time() - t_start
    log.info(f"Done: {total_mappings} mappings from {total} families "
             f"(skipped {skipped}) in {elapsed:.0f}s")
    conn.close()


if __name__ == "__main__":
    main()
