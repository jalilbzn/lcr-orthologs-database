from collections import Counter, defaultdict
from dataclasses import dataclass

from .data import ONE_BASED, Region


MAMMALS = {"Chimp", "Mouse", "Rat", "Macaque", "Cow", "Pig", "Dog", "Opossum",
           "PANTR", "MOUSE", "RATNO", "MACMU", "BOVIN", "PIG", "CANLF", "MONDO"}
SCORABLE = frozenset("ACDEFGHIKLMNPQRSTVWYBZX*")


@dataclass(frozen=True)
class AnchorSettings:
    flank_length: int = 100
    min_coverage: float = 0.5
    min_interval_length: int = 1
    max_interval_length: int = 500

    def __post_init__(self):
        if self.flank_length < 10 or not (0 < self.min_coverage <= 1):
            raise ValueError("flank_length must be >=10; min_coverage must be in (0,1]")
        if not (1 <= self.min_interval_length <= self.max_interval_length):
            raise ValueError("invalid anchor interval length limits")


def make_aligners():
    try:
        from Bio.Align import PairwiseAligner, substitution_matrices
    except ImportError as error:
        raise RuntimeError("Rebuilding anchors requires Biopython: pip install biopython") from error
    aligners = {}
    for name in ("BLOSUM62", "BLOSUM45"):
        aligner = PairwiseAligner()
        aligner.mode = "local"
        aligner.substitution_matrix = substitution_matrices.load(name)
        aligner.open_gap_score = -11
        aligner.extend_gap_score = -1
        aligners[name] = aligner
    return aligners


def sw_hit(query, target, aligner):
    """Return target [start,end) and query coverage for the first best local alignment."""
    if len(query) < 10 or len(target) < 10:
        return None
    # The legacy mapper also rejected unsupported residues through its alignment error path.
    # Preserve that policy explicitly, rather than replacing biological residues.
    if not set(query).issubset(SCORABLE) or not set(target).issubset(SCORABLE):
        return None
    # Index directly: len(alignments) can overflow for many equally scoring alignments.
    try:
        top = aligner.align(query, target)[0]
    except IndexError:
        return None
    q_blocks, t_blocks = top.aligned
    if not len(q_blocks):
        return None
    coverage = sum(int(end - start) for start, end in q_blocks) / len(query)
    return int(t_blocks[0][0]), int(t_blocks[-1][1]), round(coverage, 3)


def interval_between_anchors(source, human_sequence, target_sequence, aligner, settings):
    """Infer an ortholog interval; one-anchor fallback retains the legacy 2L+1 span."""
    start, end = source.bounds()
    left = human_sequence[max(0, start - settings.flank_length):start]
    right = human_sequence[end:end + settings.flank_length]
    left_hit = sw_hit(left, target_sequence, aligner)
    right_hit = sw_hit(right, target_sequence, aligner)
    left_ok = left_hit is not None and left_hit[2] >= settings.min_coverage
    right_ok = right_hit is not None and right_hit[2] >= settings.min_coverage
    if not left_ok and not right_ok:
        return None
    if left_ok and right_ok:
        left_bound, right_bound = left_hit[1], right_hit[0]
    elif left_ok:
        left_bound = left_hit[1]
        right_bound = min(len(target_sequence), left_bound + 2 * (end - start) + 1)
    else:
        right_bound = right_hit[0]
        left_bound = max(0, right_bound - 2 * (end - start) - 1)
    if not (settings.min_interval_length <= right_bound - left_bound <= settings.max_interval_length):
        return None
    return left_bound, right_bound


def rebuild(proteins, regions, families, settings, family_ids=None, aligners=None):
    """Recompute mappings from human seeds and cached OMA family/protein inputs."""
    aligners = make_aligners() if aligners is None else aligners
    by_acc = {p.uniprot_acc: p for p in proteins.values() if p.uniprot_acc}
    seeds = defaultdict(list)
    for region in regions.values():
        if proteins[region.protein_id].species == "Homo sapiens":
            seeds[region.protein_id].append(region)
    # Old mapped regions are inputs for replay only, never carried into a rebuild.
    output_regions = {r.lcr_id: r for values in seeds.values() for r in values}
    next_id = max(output_regions, default=0) + 1
    groups, counts = {}, Counter()
    selected = sorted(family_ids if family_ids is not None else families)
    for number, fid in enumerate(selected, 1):
        family = families[fid]
        human = by_acc[family.source_gene]
        for species, pid in sorted(family.members.items()):
            target = proteins[pid]
            name = "BLOSUM62" if species in MAMMALS or target.species in MAMMALS else "BLOSUM45"
            aligner = aligners[name]
            for source in sorted(seeds[human.protein_id], key=lambda r: (r.start, r.lcr_id)):
                counts["attempted"] += 1
                interval = interval_between_anchors(source, human.sequence, target.sequence,
                                                     aligner, settings)
                if interval is None:
                    counts["unmapped"] += 1
                    continue
                left, right = interval
                region = Region(next_id, pid, left + 1, right, target.sequence[left:right],
                                "anchor", ONE_BASED)
                output_regions[next_id] = region
                groups.setdefault((fid, source.lcr_id), set()).add(next_id)
                next_id += 1
                counts["mapped"] += 1
        if number % 100 == 0:
            print(f"Anchors: {number:,}/{len(selected):,} families; {counts['mapped']:,} mapped", flush=True)
    if not groups:
        raise ValueError("Rebuild produced no anchored ortholog groups")
    return output_regions, groups, dict(counts)
