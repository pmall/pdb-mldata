import argparse
import requests
import csv
import sys
from pathlib import Path

def fetch_pdb_metadata(output_csv: Path, min_length: int, max_length: int):
    """
    Fetches PDB entries matching the simple criteria:
    (Protein AND Length) AND (Protein Count >= 2)
    """
    print(f"Fetching PDB entries (length {min_length}-{max_length}, min 2 proteins)...")
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    
    search_url = "https://search.rcsb.org/rcsbsearch/v2/query"
    query = {
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
                                "value": "Protein"
                            }
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
                                    "include_upper": True
                                }
                            }
                        }
                    ]
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_entry_info.polymer_entity_count_protein",
                        "operator": "greater_or_equal",
                        "value": 2
                    }
                }
            ]
        },
        "request_options": {
            "return_all_hits": True
        },
        "return_type": "entry"
    }
    
    try:
        response = requests.post(search_url, json=query)
        response.raise_for_status()
    except requests.exceptions.RequestException as e:
        print(f"Error calling RCSB API: {e}", file=sys.stderr)
        sys.exit(1)
        
    data = response.json()
    result_set = data.get("result_set", [])
    
    # Extract the PDB IDs
    pdb_ids = sorted(list(set(hit["identifier"] for hit in result_set)))
    print(f"Found {len(pdb_ids)} unique PDB entries.")
    
    # Save to CSV
    with output_csv.open("w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["pdb_id"])
        for pdb_id in pdb_ids:
            writer.writerow([pdb_id])
            
    print(f"Saved metadata to {output_csv}.")


def main():
    parser = argparse.ArgumentParser(description="Fetch PDB metadata matching simple criteria.")
    
    parser.add_argument(
        "output_csv", 
        nargs="?", 
        default="data/metadata.csv", 
        help="Path to the output metadata CSV (default: 'data/metadata.csv')"
    )
    
    parser.add_argument("--min-length", type=int, default=4)
    parser.add_argument("--max-length", type=int, default=40)
    
    args = parser.parse_args()
    fetch_pdb_metadata(Path(args.output_csv), args.min_length, args.max_length)

if __name__ == "__main__":
    main()
