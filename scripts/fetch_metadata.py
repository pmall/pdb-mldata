import argparse
import csv
import sys
from pathlib import Path

import requests

DEFAULT_METADATA_CSV = Path("data/metadata.csv")
DEFAULT_MIN_LENGTH = 4
DEFAULT_MAX_LENGTH = 40


def build_metadata_query(min_length: int, max_length: int) -> dict[str, object]:
    """Build the RCSB query for peptide-length protein entries with at least two proteins."""
    return {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                {
                    "type": "group",
                    "logical_operator": "and",
                    "nodes": [
                        {
                            "type": "terminal",
                            "service": "text",
                            "parameters": {
                                "attribute": "entity_poly.rcsb_entity_polymer_type",
                                "operator": "exact_match",
                                "value": "Protein",
                            },
                        },
                        {
                            "type": "terminal",
                            "service": "text",
                            "parameters": {
                                "attribute": "entity_poly.rcsb_sample_sequence_length",
                                "operator": "range",
                                "value": {
                                    "from": min_length,
                                    "to": max_length,
                                    "include_lower": True,
                                    "include_upper": True,
                                },
                            },
                        },
                    ],
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_entry_info.polymer_entity_count_protein",
                        "operator": "greater_or_equal",
                        "value": 2,
                    },
                },
            ],
        },
        "request_options": {"return_all_hits": True},
        "return_type": "entry",
    }


def validate_parameters(output_csv: Path, min_length: int, max_length: int) -> None:
    """Validate CLI parameters before the network query runs."""
    if min_length < 1:
        raise ValueError("--min-length must be at least 1")
    if max_length < min_length:
        raise ValueError("--max-length must be greater than or equal to --min-length")
    if output_csv.exists() and output_csv.is_dir():
        raise ValueError("output_csv must be a file path, not a directory")


def fetch_pdb_metadata(output_csv: Path, min_length: int, max_length: int) -> None:
    """Fetch matching PDB IDs and write them to a one-column metadata CSV."""
    print(f"Fetching PDB entries (length {min_length}-{max_length}, min 2 proteins)...")
    output_csv.parent.mkdir(parents=True, exist_ok=True)

    search_url = "https://search.rcsb.org/rcsbsearch/v2/query"
    query = build_metadata_query(min_length, max_length)

    try:
        response = requests.post(search_url, json=query)
        response.raise_for_status()
    except requests.exceptions.RequestException as exc:
        raise RuntimeError(f"Error calling RCSB API: {exc}") from exc

    data = response.json()
    result_set = data.get("result_set", [])

    pdb_ids = sorted({hit["identifier"] for hit in result_set})
    print(f"Found {len(pdb_ids)} unique PDB entries.")

    with output_csv.open("w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["pdb_id"])
        for pdb_id in pdb_ids:
            writer.writerow([pdb_id])

    print(f"Saved metadata to {output_csv}.")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Fetch PDB metadata matching simple criteria."
    )

    parser.add_argument(
        "output_csv",
        nargs="?",
        type=Path,
        default=DEFAULT_METADATA_CSV,
        help="Path to the output metadata CSV (default: data/metadata.csv)",
    )

    parser.add_argument(
        "--min-length",
        type=int,
        default=DEFAULT_MIN_LENGTH,
        help="Minimum peptide query length (default: 4)",
    )
    parser.add_argument(
        "--max-length",
        type=int,
        default=DEFAULT_MAX_LENGTH,
        help="Maximum peptide query length (default: 40)",
    )

    args = parser.parse_args()
    try:
        validate_parameters(args.output_csv, args.min_length, args.max_length)
        fetch_pdb_metadata(args.output_csv, args.min_length, args.max_length)
    except (RuntimeError, ValueError) as exc:
        print(exc, file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
