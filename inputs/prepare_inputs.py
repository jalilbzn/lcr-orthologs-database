import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from inputs.pipeline.upstream import prepare


def main():
    root = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--human-lcrs", type=Path, default=root / "raw/human_lcrs.csv")
    parser.add_argument("--oma-dir", type=Path, default=root / "raw/oma")
    parser.add_argument("--output", type=Path, default=root / "inputs_rebuilt")
    parser.add_argument("--uniprot-release", help="optional declared label, not independent version verification")
    parser.add_argument("--oma-release", help="optional declared label, not independent version verification")
    parser.add_argument("--build-anchors", action="store_true",
                        help="also recompute anchor mappings and prepare all four CSV tables; requires Biopython")
    args = parser.parse_args()
    try:
        manifest = prepare(args.human_lcrs, args.oma_dir, args.output, args.uniprot_release, args.oma_release,
                           build_anchors=args.build_anchors)
    except (OSError, ValueError, RuntimeError) as error:
        print(f"Error: {error}", file=sys.stderr)
        return 1
    print(f"Prepared {args.output}: {manifest['counts']}")
    flag = "" if args.build_anchors else " --rebuild-anchors"
    print(f"Next: python run_pipeline.py --inputs <prepared-directory>{flag} --output <new-output-directory>")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
