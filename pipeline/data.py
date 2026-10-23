import ast
import csv
import gzip
import re
from dataclasses import dataclass
from pathlib import Path


OMA_ID = re.compile(r"[A-Z]{5}[0-9]+\Z")
UNIPROT_ACC = re.compile(
    r"(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})(?:-[0-9]+)?\Z"
)
ZERO_BASED = "0-based half-open"
ONE_BASED = "1-based inclusive"


@dataclass(frozen=True, slots=True)
class Protein:
    protein_id: int
    uniprot_acc: str
    oma_id: str
    taxid: str
    species: str
    sequence: str


@dataclass(frozen=True, slots=True)
class Region:
    lcr_id: int
    protein_id: int
    start: int
    end: int
    sequence: str
    seg_mode: str
    coordinate_system: str

    def bounds(self):
        """Return zero-based, half-open protein coordinates."""
        if self.coordinate_system == ZERO_BASED:
            return self.start, self.end
        if self.coordinate_system == ONE_BASED:
            return self.start - 1, self.end
        raise ValueError(f"LCR {self.lcr_id}: unknown coordinate system")


@dataclass(frozen=True, slots=True)
class Family:
    family_id: int
    source_gene: str
    members: dict


def rows(path, required):
    with Path(path).open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        missing = set(required) - set(reader.fieldnames or ())
        if missing:
            raise ValueError(f"{path}: missing columns {', '.join(sorted(missing))}")
        for number, row in enumerate(reader, 2):
            if None in row or any(row.get(key) is None for key in required):
                raise ValueError(f"{path}:{number}: malformed CSV row")
            yield row


def unique_add(table, key, value, path):
    if key in table:
        raise ValueError(f"{path}: duplicate ID {key}")
    table[key] = value


def load_proteins(path):
    proteins = {}
    for row in rows(path, ("protein_id", "uniprot_acc", "species", "sequence")):
        pid = int(row["protein_id"])
        acc = row["uniprot_acc"].strip()
        oma = row.get("oma_id", "").strip()
        # Historical protein.csv put OMA IDs in the UniProt column.
        if OMA_ID.fullmatch(acc):
            if oma and oma != acc:
                raise ValueError(f"Protein {pid}: conflicting OMA identifiers")
            oma, acc = acc, ""
        if acc and not UNIPROT_ACC.fullmatch(acc):
            raise ValueError(f"Protein {pid}: invalid UniProt accession {acc!r}")
        if oma and not OMA_ID.fullmatch(oma):
            raise ValueError(f"Protein {pid}: invalid OMA identifier {oma!r}")
        seq = row["sequence"].strip()
        if not seq or not seq.isalpha() or seq != seq.upper():
            raise ValueError(f"Protein {pid}: expected an uppercase amino-acid sequence")
        if row.get("length") and int(row["length"]) != len(seq):
            raise ValueError(f"Protein {pid}: length differs from sequence")
        protein = Protein(pid, acc, oma, row.get("taxid", "").strip(),
                          row["species"].strip(), seq)
        unique_add(proteins, pid, protein, path)
    return proteins


def load_regions(path, proteins):
    regions = {}
    for row in rows(path, ("lcr_id", "protein_id", "start_pos", "end_pos", "sequence", "seg_mode")):
        lid, pid = int(row["lcr_id"]), int(row["protein_id"])
        if pid not in proteins:
            raise ValueError(f"LCR {lid}: missing protein {pid}")
        protein = proteins[pid]
        system = row.get("coordinate_system", "").strip()
        start, end = int(row["start_pos"]), int(row["end_pos"])
        seq = row["sequence"].strip()
        if not system:
            # Only infer a convention when its coordinates exactly reproduce the sequence.
            if 0 <= start < end <= len(protein.sequence) and protein.sequence[start:end] == seq:
                system = ZERO_BASED
            elif 1 <= start <= end <= len(protein.sequence) and protein.sequence[start - 1:end] == seq:
                system = ONE_BASED
            else:
                raise ValueError(f"LCR {lid}: cannot infer coordinate convention")
        region = Region(lid, pid, start, end, seq, row["seg_mode"].strip(), system)
        left, right = region.bounds()
        if not (0 <= left < right <= len(protein.sequence)) or protein.sequence[left:right] != seq:
            raise ValueError(f"LCR {lid}: coordinates do not match its protein subsequence")
        unique_add(regions, lid, region, path)
    return regions


def load_families(path, proteins):
    families = {}
    accessions = {p.uniprot_acc for p in proteins.values() if p.uniprot_acc}
    for row in rows(path, ("family_id", "source_gene", "members")):
        fid = int(row["family_id"])
        # Support both original Python-dict CSV cells and compact JSON cells. Never eval.
        members = ast.literal_eval(row["members"])
        if not isinstance(members, dict):
            raise ValueError(f"Family {fid}: members must be a dictionary")
        members = {str(species): int(pid) for species, pid in members.items()}
        if row["source_gene"] not in accessions:
            raise ValueError(f"Family {fid}: missing human source protein")
        for pid in members.values():
            if pid not in proteins:
                raise ValueError(f"Family {fid}: missing member protein {pid}")
        unique_add(families, fid, Family(fid, row["source_gene"], members), path)
    return families


def load_mappings(path, proteins, regions, families):
    """Replay anchored relationships without applying later analysis/weighting columns."""
    groups = {}
    for row in rows(path, ("family_id", "source_lcr_id", "lcr_id")):
        if not row["source_lcr_id"].strip():
            continue
        fid, src, lid = (int(row[k]) for k in ("family_id", "source_lcr_id", "lcr_id"))
        if fid not in families or src not in regions or lid not in regions:
            raise ValueError(f"Mapping {fid}/{src}/{lid}: missing referenced record")
        family, source, target = families[fid], regions[src], regions[lid]
        if proteins[source.protein_id].uniprot_acc != family.source_gene:
            raise ValueError(f"Mapping {fid}/{src}/{lid}: wrong source protein")
        if proteins[source.protein_id].species != "Homo sapiens":
            raise ValueError(f"Mapping {fid}/{src}/{lid}: source must be human")
        if lid == src or target.protein_id not in family.members.values():
            raise ValueError(f"Mapping {fid}/{src}/{lid}: target is not a family ortholog")
        if target.seg_mode != "anchor":
            raise ValueError(f"Mapping {fid}/{src}/{lid}: target is not anchor-derived")
        groups.setdefault((fid, src), set()).add(lid)
    if not groups:
        raise ValueError("No anchored ortholog groups found")
    return groups


def load_inputs(directory):
    directory = Path(directory)
    proteins = load_proteins(directory / "protein.csv")
    regions = load_regions(directory / "lcr_seq.csv", proteins)
    families = load_families(directory / "lcr_family.csv", proteins)
    return proteins, regions, families


def read_crosswalk(path, proteins):
    """Read an optional local OMA-ID/UniProt two-column TSV(.gz), preserving all matches."""
    wanted_oma = {p.oma_id for p in proteins.values() if p.oma_id}
    wanted_acc = {p.uniprot_acc for p in proteins.values() if p.uniprot_acc}
    by_oma, by_uniprot = {}, {}
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8-sig") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 2:
                raise ValueError(f"{path}: expected OMA ID and UniProt accession")
            oma, acc = parts[:2]
            if not OMA_ID.fullmatch(oma) or not UNIPROT_ACC.fullmatch(acc):
                raise ValueError(f"{path}: invalid crosswalk pair {oma!r}, {acc!r}")
            if oma in wanted_oma or acc in wanted_acc:
                by_oma.setdefault(oma, set()).add(acc)
                by_uniprot.setdefault(acc, set()).add(oma)
    return by_oma, by_uniprot
