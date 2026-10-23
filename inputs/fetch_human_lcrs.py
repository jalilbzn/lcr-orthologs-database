#!/usr/bin/env python3
import argparse
import json
import sys
from pathlib import Path
from urllib.request import Request, urlopen

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from inputs.pipeline.upstream import API_ENDPOINT, HUMAN_COLUMNS, file_record, human_reader, utc_now


def fetch(output, opener=None):
    output = Path(output)
    receipt = output.with_suffix(output.suffix + ".request.json")
    if output.exists() or receipt.exists():
        raise ValueError("Response or request receipt already exists; choose a new --output")
    payload = {"result_data": "lcrs", "columns": HUMAN_COLUMNS,
               "search_by_column": {"Organism": ["Homo sapiens"]}}
    request = Request(API_ENDPOINT, data=json.dumps(payload).encode("utf-8"), method="POST",
                      headers={"Content-Type": "application/json", "Accept": "text/csv"})
    with (opener or urlopen)(request, timeout=300) as response:
        content = response.read()
        response_headers = {key: response.headers.get(key) for key in ("Date", "ETag", "Last-Modified", "Content-Type")}
    human_reader(content.decode("utf-8-sig"))  # Reject an HTML/error response before saving it.
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("xb") as handle:
        handle.write(content)
    record = {"endpoint": API_ENDPOINT, "method": "POST", "request": payload,
              "retrieved_at_utc": utc_now(), "response_headers": response_headers,
              "response": file_record(output),
              "upstream_release": None,
              "release_note": "This is a fresh API response, not a pinned historical UniProt release."}
    with receipt.open("x", encoding="utf-8") as handle:
        handle.write(json.dumps(record, indent=2) + "\n")
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path(__file__).resolve().parent / "raw/human_lcrs.csv")
    args = parser.parse_args()
    try:
        record = fetch(args.output)
    except (OSError, ValueError) as error:
        print(f"Error: {error}", file=sys.stderr)
        return 1
    print(f"Saved {args.output}; SHA-256: {record['response']['sha256']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
