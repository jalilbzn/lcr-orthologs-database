import gzip
import json
import logging
import os
import sys
import time
import psycopg2
import psycopg2.extras
from ortholog_lcr_pipeline.inputs.pipeline.config import DB_CONFIG, LOGS_DIR, OMA_DIR

# ── Logging (file + console) ──
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
    handlers=[
        logging.FileHandler(LOGS_DIR / "seed_bulk.log", mode='w'),
        logging.StreamHandler(sys.stdout),
    ],
)
log = logging.getLogger("seed_bulk")

FILE_UNIPROT = os.path.join(OMA_DIR, "oma-uniprot.txt.gz")
FILE_GROUPS  = os.path.join(OMA_DIR, "oma-groups.txt.gz")
FILE_SEQS    = os.path.join(OMA_DIR, "oma-seqs.fa.gz")

TARGET_SPECIES = {
    "MOUSE": ("Mouse",    10090),
    "RATNO": ("Rat",      10116),
    "BOVIN": ("Cow",       9913),
    "CANLF": ("Dog",       9615),
    "PIG":   ("Pig",       9823),
    "MACMU": ("Macaque",   9544),
    "PANTR": ("Chimp",     9598),
    "MONDO": ("Opossum",  13616),
    "CHICK": ("Chicken",   9031),
    "XENTR": ("Xenopus",   8364),
    "DANRE": ("Zebrafish", 7955),
}

# ── Phase 1 ──────────────────────────────────────────────────
def phase1_map_uniprot(human_accs):
    """Parse oma-uniprot -> {OMA_ID: UniProt_ACC} for our 16k humans."""
    log.info("Phase 1: Mapping UniProt -> OMA IDs...")
    t0 = time.time()
    mapping = {}

    with gzip.open(FILE_UNIPROT, 'rt') as f:
        for line in f:
            if line[0] == '#':
                continue
            parts = line.split('\t')
            if len(parts) < 2:
                continue
            oma_id = parts[0].strip()
            uni_acc = parts[1].strip()
            if oma_id[:5] == "HUMAN" and uni_acc in human_accs:
                mapping[oma_id] = uni_acc

    log.info(f"Phase 1 done: {len(mapping)} OMA IDs mapped in {time.time()-t0:.1f}s")
    return mapping


# ── Phase 2 ──────────────────────────────────────────────────
def phase2_find_orthologs(human_oma_map):
    """Scan oma-groups -> orthologs_needed + family_structure."""
    log.info("Phase 2: Scanning OMA Groups...")
    t0 = time.time()

    orthologs_needed = {}   # {OMA_ID: (human_acc, species_code)}
    family_structure = {}   # {human_acc: {species_code: OMA_ID}}
    target_codes = list(TARGET_SPECIES.keys())

    with gzip.open(FILE_GROUPS, 'rt') as f:
        for line in f:
            if line[0] == '#':
                continue
            parts = line.split()
            if len(parts) < 3:
                continue

            group_ids = parts[1:]

            # Find our human seed
            human_acc = None
            for gid in group_ids:
                if gid in human_oma_map:
                    human_acc = human_oma_map[gid]
                    break
            if not human_acc:
                continue

            # Collect target-species members
            members = {}
            for gid in group_ids:
                for code in target_codes:
                    if gid.startswith(code):
                        if code not in members:
                            members[code] = gid
                            orthologs_needed[gid] = (human_acc, code)
                        break

            family_structure[human_acc] = members

    log.info(f"Phase 2 done: {len(family_structure)} families, "
             f"{len(orthologs_needed)} orthologs in {time.time()-t0:.1f}s")
    return orthologs_needed, family_structure


# ── Phase 3 ──────────────────────────────────────────────────
def phase3_extract_sequences(orthologs_needed):
    """Stream oma-seqs.fa.gz, extract only needed sequences."""
    log.info(f"Phase 3: Streaming 5GB FASTA for {len(orthologs_needed)} sequences...")
    t0 = time.time()

    sequences = {}
    needed_set = set(orthologs_needed.keys())
    total_needed = len(needed_set)

    current_id = None
    is_needed = False
    current_seq = []
    lines_scanned = 0

    with gzip.open(FILE_SEQS, 'rt') as f:
        for line in f:
            lines_scanned += 1

            if line[0] == '>':
                # Save previous if needed
                if is_needed and current_id:
                    sequences[current_id] = "".join(current_seq)
                    if len(sequences) % 5000 == 0:
                        elapsed = time.time() - t0
                        log.info(f"  {len(sequences)}/{total_needed} sequences "
                                 f"({lines_scanned:,} lines, {elapsed:.0f}s)")
                    if len(sequences) >= total_needed:
                        log.info("  All sequences found — stopping early.")
                        break

                # Parse new header
                current_id = line[1:].split()[0]
                is_needed = current_id in needed_set
                current_seq = []
            elif is_needed:
                current_seq.append(line.rstrip())

        # Handle last entry
        if is_needed and current_id and current_id not in sequences:
            sequences[current_id] = "".join(current_seq)

    elapsed = time.time() - t0
    log.info(f"Phase 3 done: {len(sequences)}/{total_needed} sequences "
             f"in {elapsed:.0f}s ({lines_scanned:,} lines)")
    return sequences


# ── Phase 4 ──────────────────────────────────────────────────
def phase4_save_to_db(conn, family_structure, sequences, orthologs_needed):
    """Batch insert proteins + families."""
    log.info("Phase 4: Saving to database...")
    t0 = time.time()

    # Build protein rows
    prot_rows = []
    for oma_id, seq in sequences.items():
        _, code = orthologs_needed[oma_id]
        species_name, taxid = TARGET_SPECIES[code]
        prot_rows.append((oma_id, taxid, species_name, len(seq), seq))

    # Bulk upsert proteins
    oma_to_pid = {}
    if prot_rows:
        log.info(f"  Inserting {len(prot_rows)} proteins...")
        with conn.cursor() as cur:
            sql = """
                INSERT INTO protein (uniprot_acc, taxid, species, length, sequence)
                VALUES %s
                ON CONFLICT (uniprot_acc) DO UPDATE SET sequence = EXCLUDED.sequence
                RETURNING protein_id, uniprot_acc
            """
            psycopg2.extras.execute_values(cur, sql, prot_rows, page_size=500, fetch=True)
            for pid, acc in cur.fetchall():
                oma_to_pid[acc] = pid
        conn.commit()
        log.info(f"  Proteins committed.")

    # Build family rows
    family_rows = []
    for human_acc, members in family_structure.items():
        member_pids = {}
        for code, oma_id in members.items():
            if oma_id in oma_to_pid:
                member_pids[code] = oma_to_pid[oma_id]
        family_rows.append((
            f"Ortholog family of {human_acc}",
            human_acc,
            json.dumps(member_pids),
        ))

    if family_rows:
        log.info(f"  Inserting {len(family_rows)} families...")
        with conn.cursor() as cur:
            sql = """
                INSERT INTO lcr_family (description, source_gene, members)
                VALUES %s
            """
            psycopg2.extras.execute_values(cur, sql, family_rows, page_size=500)
        conn.commit()
        log.info(f"  Families committed.")

    log.info(f"Phase 4 done in {time.time()-t0:.1f}s")


# ── Main ─────────────────────────────────────────────────────
def main():
    t_start = time.time()
    conn = psycopg2.connect(**DB_CONFIG)

    # Load human proteins from DB
    with conn.cursor() as cur:
        cur.execute("SELECT uniprot_acc FROM protein WHERE species='Homo sapiens'")
        human_accs = {r[0] for r in cur.fetchall()}
    log.info(f"DB: {len(human_accs)} human proteins loaded.")

    if not human_accs:
        log.error("No human proteins in DB. Run Step 1 first.")
        return

    # Execute phases
    human_oma_map = phase1_map_uniprot(human_accs)
    orthologs_needed, family_structure = phase2_find_orthologs(human_oma_map)
    sequences = phase3_extract_sequences(orthologs_needed)
    phase4_save_to_db(conn, family_structure, sequences, orthologs_needed)

    conn.close()
    log.info(f"All done in {time.time()-t_start:.0f}s total.")


if __name__ == "__main__":
    main()
