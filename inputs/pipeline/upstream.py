import csv
import gzip
import hashlib
import importlib.metadata
import io
import json
import platform
import tempfile
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path

from pipeline.data import ONE_BASED, ZERO_BASED, OMA_ID, UNIPROT_ACC, Region, load_inputs, load_mappings


API_ENDPOINT = "https://lcrannotdb.lcr-lab.org/api/proteins/"
HUMAN_COLUMNS = [
    "Protein ID", "UniProtACC", "Organism", "Protein Sequence", "Protein name",
    "LCR ID", "LCR Sequence", "Start of LCR", "End of LCR", "LCR length", "Identification method",
]
REQUIRED_HUMAN = {"UniProtACC", "Organism", "Protein Sequence", "LCR Sequence",
                  "Start of LCR", "End of LCR", "Identification method"}
OMA_FILES = ("oma-uniprot.txt.gz", "oma-groups.txt.gz", "oma-seqs.fa.gz")
TARGETS = {
    "MOUSE": ("Mouse", "10090"), "RATNO": ("Rat", "10116"),
    "BOVIN": ("Cow", "9913"), "CANLF": ("Dog", "9615"),
    "PIG": ("Pig", "9823"), "MACMU": ("Macaque", "9544"),
    "PANTR": ("Chimp", "9598"), "MONDO": ("Opossum", "13616"),
    "CHICK": ("Chicken", "9031"), "XENTR": ("Xenopus", "8364"),
    "DANRE": ("Zebrafish", "7955"),
}


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def file_record(path):
    path = Path(path)
    return {"filename": path.name, "size_bytes": path.stat().st_size, "sha256": sha256(path)}


def utc_now():
    return datetime.now(timezone.utc).isoformat()


def mode(method):
    method = method.strip().lower()
    if "strict" in method:
        return "strict"
    if not method or "default" in method or "seg" in method:
        return "default"
    return "custom"


def human_reader(text):
    # The historical API returns semicolon-delimited CSV. Also accept saved comma CSV.
    text = text.lstrip("\ufeff")
    delimiter = ";" if ";" in text.partition("\n")[0] else ","
    reader = csv.DictReader(io.StringIO(text, newline=""), delimiter=delimiter)
    missing = REQUIRED_HUMAN - set(reader.fieldnames or ())
    if missing:
        raise ValueError(f"Human CSV is missing columns: {', '.join(sorted(missing))}")
    return reader


def read_human_seeds(path):
    proteins, seeds = {}, []
    for number, row in enumerate(human_reader(Path(path).read_text(encoding="utf-8-sig")), 2):
        if None in row or any(row.get(key) is None for key in REQUIRED_HUMAN):
            raise ValueError(f"Human CSV row {number}: malformed record")
        if row["Organism"].strip() != "Homo sapiens":
            continue
        acc = row["UniProtACC"].strip()
        seq = row["Protein Sequence"].strip().strip('"').strip()
        lcr = row["LCR Sequence"].strip().strip('"').strip()
        if not UNIPROT_ACC.fullmatch(acc) or not seq or not seq.isalpha() or seq != seq.upper():
            raise ValueError(f"Human CSV row {number}: invalid accession or protein sequence")
        if acc in proteins and proteins[acc] != seq:
            raise ValueError(f"Human CSV row {number}: conflicting sequences for {acc}")
        start, end = int(row["Start of LCR"]), int(row["End of LCR"])
        if lcr and 0 <= start < end <= len(seq) and seq[start:end] == lcr:
            system = ZERO_BASED
        elif lcr and 1 <= start <= end <= len(seq) and seq[start - 1:end] == lcr:
            system = ONE_BASED
        else:
            raise ValueError(f"Human CSV row {number}: coordinates do not reproduce the LCR sequence")
        proteins[acc] = seq
        seeds.append((acc, start, end, lcr, mode(row["Identification method"]), system,
                      row.get("LCR ID", "").strip()))
    if not seeds:
        raise ValueError("Human CSV has no valid Homo sapiens LCRs")
    # Keep all seed rows. Sorting provides deterministic IDs without relying on a database counter.
    return proteins, sorted(seeds)


def human_oma_mapping(path, human_accessions):
    mapping = {}
    with gzip.open(path, "rt", encoding="utf-8-sig") as handle:
        for line in handle:
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 2:
                raise ValueError(f"{path}: expected OMA ID and UniProt accession")
            oma, acc = fields[:2]
            if oma.startswith("HUMAN") and acc in human_accessions:
                if not OMA_ID.fullmatch(oma):
                    raise ValueError(f"Invalid OMA ID: {oma}")
                if oma in mapping and mapping[oma] != acc:
                    raise ValueError(f"OMA ID {oma} maps to multiple human seed accessions")
                mapping[oma] = acc
    if not mapping:
        raise ValueError("No human seed accessions matched the OMA crosswalk")
    return mapping


def select_oma_families(path, human_mapping):
    """Historical selection: first human match; first member/species; last group/human."""
    families = {}
    with gzip.open(path, "rt", encoding="utf-8-sig") as handle:
        for line in handle:
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 3:
                continue
            group, *members = fields
            human = next((human_mapping[member] for member in members if member in human_mapping), None)
            if human is None:
                continue
            selected = {}
            for member in members:
                code = next((code for code in TARGETS if member.startswith(code)), None)
                if code is not None and code not in selected:
                    if not OMA_ID.fullmatch(member):
                        raise ValueError(f"Invalid OMA group member: {member}")
                    selected[code] = member
            # Retain the original group's ID as provenance, independently of local family IDs.
            families[human] = {"oma_group_id": group, "members": selected}
    if not families:
        raise ValueError("No OMA groups matched the human seeds")
    return families


def extract_oma_sequences(path, wanted):
    sequences, current, parts = {}, None, []
    def save():
        if current not in wanted:
            return
        seq = "".join(parts)
        if current in sequences:
            raise ValueError(f"Duplicate OMA FASTA ID: {current}")
        if not seq or not seq.isalpha() or seq != seq.upper():
            raise ValueError(f"OMA FASTA {current}: expected an uppercase amino-acid sequence")
        sequences[current] = seq
    with gzip.open(path, "rt", encoding="utf-8-sig") as handle:
        for line in handle:
            if line.startswith(">"):
                save()
                current = line[1:].split()[0]
                parts = []
            elif current in wanted:
                parts.append(line.strip())
        save()
    missing = wanted - sequences.keys()
    if missing:
        raise ValueError(f"OMA FASTA is missing {len(missing)} selected proteins; examples: {sorted(missing)[:5]}")
    return sequences


def write_csv(path, headers, values):
    with Path(path).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(headers)
        writer.writerows(values)


def prepare(human_csv, oma_dir, output, uniprot_release=None, oma_release=None, build_anchors=False):
    human_csv, oma_dir, output = Path(human_csv).resolve(), Path(oma_dir).resolve(), Path(output).resolve()
    required = [human_csv, *(oma_dir / name for name in OMA_FILES)]
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise ValueError("Missing raw input files:\n" + "\n".join(missing))
    if output.exists():
        raise ValueError(f"Output already exists; choose a new directory: {output}")
    raw_stats = {path: (path.stat().st_size, path.stat().st_mtime_ns) for path in required}
    print("Reading human LCRs and OMA memberships...", flush=True)
    human_proteins, seeds = read_human_seeds(human_csv)
    mapping = human_oma_mapping(oma_dir / OMA_FILES[0], set(human_proteins))
    families = select_oma_families(oma_dir / OMA_FILES[1], mapping)
    needed = {oma for family in families.values() for oma in family["members"].values()}
    if not needed:
        raise ValueError("Selected OMA groups have no target-species members")
    print(f"Reading sequences for {len(needed):,} selected OMA proteins...", flush=True)
    orthologs = extract_oma_sequences(oma_dir / OMA_FILES[2], needed)
    human_ids = {acc: number for number, acc in enumerate(sorted(human_proteins), 1)}
    oma_ids = {oma: number for number, oma in enumerate(sorted(orthologs), len(human_ids) + 1)}
    species = {oma: code for family in families.values() for code, oma in family["members"].items()}
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".prepare-inputs-", dir=output.parent) as temporary:
        staged = Path(temporary) / "inputs"
        staged.mkdir()
        protein_rows = [[human_ids[acc], acc, "", "", "Homo sapiens", len(seq), seq]
                        for acc, seq in sorted(human_proteins.items())]
        protein_rows.extend([oma_ids[oma], "", oma, TARGETS[species[oma]][1],
                             TARGETS[species[oma]][0], len(seq), seq] for oma, seq in sorted(orthologs.items()))
        write_csv(staged / "protein.csv",
                  ["protein_id", "uniprot_acc", "oma_id", "taxid", "species", "length", "sequence"], protein_rows)
        write_csv(staged / "lcr_seq.csv",
                  ["lcr_id", "protein_id", "start_pos", "end_pos", "sequence", "seg_mode", "coordinate_system", "lcrannotationsdb_id"],
                  ([number, human_ids[acc], start, end, seq, seg, system, external]
                   for number, (acc, start, end, seq, seg, system, external) in enumerate(seeds, 1)))
        write_csv(staged / "lcr_family.csv", ["family_id", "source_gene", "members", "oma_group_id"],
                  ([number, acc, json.dumps({TARGETS[code][0]: oma_ids[oma] for code, oma in sorted(family["members"].items())}), family["oma_group_id"]]
                   for number, (acc, family) in enumerate(sorted(families.items()), 1)))
        prepared_proteins, prepared_regions, prepared_families = load_inputs(staged)
        alignment = None
        if build_anchors:
            from pipeline import anchors
            settings = anchors.AnchorSettings()
            mapped_regions, groups, alignment_counts = anchors.rebuild(
                prepared_proteins, prepared_regions, prepared_families, settings)
            # Match run_pipeline.py --rebuild-anchors: export all regions as one-based inclusive.
            mapped_regions = {lid: Region(r.lcr_id, r.protein_id, r.bounds()[0] + 1, r.bounds()[1],
                                         r.sequence, r.seg_mode, ONE_BASED)
                              for lid, r in mapped_regions.items()}
            external_ids = {number: seed[-1] for number, seed in enumerate(seeds, 1)}
            write_csv(staged / "lcr_seq.csv",
                      ["lcr_id", "protein_id", "start_pos", "end_pos", "sequence", "seg_mode", "coordinate_system", "lcrannotationsdb_id"],
                      ([r.lcr_id, r.protein_id, r.start, r.end, r.sequence, r.seg_mode, r.coordinate_system,
                        external_ids.get(lid, "")] for lid, r in sorted(mapped_regions.items())))
            write_csv(staged / "phylogeny_map.csv", ["family_id", "source_lcr_id", "lcr_id"],
                      ([fid, src, lid] for (fid, src), lids in sorted(groups.items()) for lid in sorted(lids)))
            prepared_proteins, prepared_regions, prepared_families = load_inputs(staged)
            load_mappings(staged / "phylogeny_map.csv", prepared_proteins, prepared_regions, prepared_families)
            try:
                biopython_version = importlib.metadata.version("biopython")
            except importlib.metadata.PackageNotFoundError:
                biopython_version = None
            alignment = {"settings": asdict(settings), "biopython_version": biopython_version,
                         "counts": alignment_counts, "groups": len(groups),
                         "alignment_code": file_record(anchors.__file__)}
        print("Recording raw-file checksums (the OMA FASTA may be large)...", flush=True)
        receipt = human_csv.with_suffix(human_csv.suffix + ".request.json")
        manifest = {
            "prepared_at_utc": utc_now(),
            "python_version": platform.python_version(),
            "raw_files": [file_record(path) for path in required],
            "human_request_receipt": None,
            "declared_releases": {"uniprot": uniprot_release, "oma": oma_release},
            "release_verification": "Labels are supplied by the caller; not independently verified.",
            "selection_policy": "First matching human per OMA group; first target member per species; last encountered group per human.",
            "coordinates": "Inferred by exact protein-subsequence match; stored explicitly per seed.",
            "id_policy": "New deterministic local IDs; not the original PostgreSQL IDs.",
            "counts": {"human_proteins": len(human_proteins), "human_lcrs": len(seeds),
                       "ortholog_proteins": len(orthologs), "families": len(families)},
            "prepared_files": [file_record(path) for path in sorted(staged.glob("*.csv"))],
            "preparation_code": file_record(Path(__file__)),
            "anchor_rebuild": alignment,
        }
        if receipt.is_file():
            recorded = json.loads(receipt.read_text(encoding="utf-8"))
            if not isinstance(recorded, dict) or not isinstance(recorded.get("response"), dict) or recorded["response"].get("sha256") != manifest["raw_files"][0]["sha256"]:
                raise ValueError("Saved human request receipt does not match the supplied CSV")
            manifest["human_request_receipt"] = recorded
        for path in required:
            if (path.stat().st_size, path.stat().st_mtime_ns) != raw_stats[path]:
                raise ValueError(f"Raw input changed during preparation: {path}")
        (staged / "sources.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
        staged.rename(output)
    return manifest
