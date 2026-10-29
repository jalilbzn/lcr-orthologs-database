import argparse
import sys
import tempfile
from pathlib import Path

from pipeline.data import ONE_BASED, Region, load_inputs, load_mappings, read_crosswalk
from pipeline.export import export_fastas, filter_groups, metadata_tables, write_metadata


ROOT = Path(__file__).resolve().parent
OUTPUT_NAMES = {"anchor_fastas_default", "anchor_fastas_strict", "metadata.xlsx"}


def parser():
    argument_parser = argparse.ArgumentParser(description=__doc__)
    argument_parser.add_argument("--inputs", type=Path, default=ROOT / "inputs" / "tables")
    argument_parser.add_argument("--output", type=Path, default=ROOT / "outputs")
    argument_parser.add_argument("--rebuild-anchors", action="store_true",
                                 help="recompute anchor mappings instead of replaying phylogeny_map.csv")
    argument_parser.add_argument("--family", type=int, action="append",
                                 help="restrict to a local family ID; repeat for several families")
    argument_parser.add_argument("--oma-uniprot", type=Path,
                                 help="optional local OMA-ID/UniProt TSV or TSV.gz crosswalk")
    argument_parser.add_argument("--min-lcr-length", type=int, default=5, metavar="AA",
                                 help="minimum source and ortholog LCR length (default: 5; 1 restores unfiltered exports)")
    argument_parser.add_argument("--overwrite", action="store_true",
                                 help="replace a previous pipeline output after successful generation")
    return argument_parser


def check_destination(inputs, output, overwrite):
    if output == inputs or output in inputs.parents or inputs in output.parents:
        raise ValueError("Input and output directories must not contain one another")
    if output.exists():
        if not output.is_dir() or output.is_symlink():
            raise ValueError("Output must be a regular directory")
        if not overwrite:
            raise ValueError(f"Output already exists: {output}. Use --overwrite to replace it.")
        unexpected = {path.name for path in output.iterdir()} - OUTPUT_NAMES
        if unexpected:
            raise ValueError(f"Output contains unrelated files; choose another directory: {sorted(unexpected)}")


def run(args, metadata_writer=write_metadata):
    if args.min_lcr_length < 1:
        raise ValueError("Minimum LCR length must be at least 1")
    if args.output.is_symlink():
        raise ValueError("Output must not be a symbolic link")
    inputs, output = args.inputs.resolve(), args.output.resolve()
    check_destination(inputs, output, args.overwrite)
    print("Reading and validating local proteins, LCRs and OMA-derived families...", flush=True)
    proteins, regions, families = load_inputs(inputs)
    selected = set(args.family) if args.family else None
    if selected is not None and selected - families.keys():
        raise ValueError(f"Unknown family IDs: {sorted(selected - families.keys())}")
    rebuild_counts = {}
    if args.rebuild_anchors:
        from pipeline.anchors import AnchorSettings, rebuild
        regions, groups, rebuild_counts = rebuild(proteins, regions, families, AnchorSettings(), selected)
        # Rebuild outputs use one consistent coordinate convention, including source seeds.
        regions = {lid: Region(r.lcr_id, r.protein_id, r.bounds()[0] + 1, r.bounds()[1],
                               r.sequence, r.seg_mode, ONE_BASED) for lid, r in regions.items()}
    else:
        groups = load_mappings(inputs / "phylogeny_map.csv", proteins, regions, families)
        if selected is not None:
            groups = {key: targets for key, targets in groups.items() if key[0] in selected}
    if not groups:
        raise ValueError("No anchored groups for the requested families")
    groups, filter_counts = filter_groups(groups, regions, args.min_lcr_length)
    print(f"Length filter >= {args.min_lcr_length} aa: "
          f"{filter_counts['short_source_groups']:,} short-source groups removed; "
          f"{filter_counts['short_ortholog_entries']:,} short ortholog entries removed from remaining groups; "
          f"{filter_counts['groups_without_orthologs']:,} groups removed after losing all orthologs. "
          f"Total removed: {filter_counts['removed_groups']:,} groups / "
          f"{filter_counts['removed_sequence_entries']:,} sequence entries.", flush=True)
    if not groups:
        raise ValueError("No anchored groups remain after length filtering")
    crosswalk = read_crosswalk(args.oma_uniprot, proteins) if args.oma_uniprot else None
    output.parent.mkdir(parents=True, exist_ok=True)
    # Stage all outputs first. Errors leave an existing successful output untouched.
    with tempfile.TemporaryDirectory(prefix=".lcr-build-", dir=output.parent) as temporary:
        staged = Path(temporary) / "outputs"
        staged.mkdir()
        counts = export_fastas(staged, groups, regions, proteins)
        tables = metadata_tables(groups, regions, proteins, crosswalk,
                                 min_lcr_length=args.min_lcr_length)
        print("Writing essential UniProt/OMA metadata...", flush=True)
        metadata_writer(staged / "metadata.xlsx", tables)
        backup = Path(temporary) / "previous"
        existed = output.exists()
        if existed:
            output.rename(backup)
        try:
            staged.rename(output)
        except BaseException:
            if existed:
                backup.rename(output)
            raise
    counts["metadata_proteins"] = len(tables["Proteins"]["rows"])
    counts.update(rebuild_counts)
    counts.update({f"length_filter_{name}": value for name, value in filter_counts.items()})
    print(f"Done: {counts['default_groups']:,} default / {counts['strict_groups']:,} strict groups; "
          f"{counts['metadata_proteins']:,} proteins in metadata.", flush=True)
    print(output, flush=True)
    return counts


def main():
    args = parser().parse_args()
    try:
        run(args)
    except (OSError, ValueError, RuntimeError) as error:
        print(f"Error: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
