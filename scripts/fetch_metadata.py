import argparse
import requests
import csv
import sys
from pathlib import Path

def fetch_pdb_metadata(output_csv: Path, min_length: int, max_length: int):
    """
    Step 1: Uses the RCSB Search API to cast a broad net and get candidate PDB entries.
    Rules:
      1. At least two protein chains in the assembly
      2. At least one chain between min_length and max_length
    """
    print(f"Fetching PDB candidates (length {min_length}-{max_length}) using Search API...")
    
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    
    search_url = "https://search.rcsb.org/rcsbsearch/v2/query"
    query = {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_assembly_info.polymer_entity_instance_count_protein",
                        "operator": "greater_or_equal",
                        "value": 2
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
        "request_options": {
            "return_all_hits": True
        },
        "return_type": "entry"  # Direct standard PDB entries
    }
    
    try:
        response = requests.post(search_url, json=query)
        response.raise_for_status()
    except requests.exceptions.RequestException as e:
        print(f"Error calling RCSB API: {e}", file=sys.stderr)
        sys.exit(1)
        
    data = response.json()
    result_set = data.get("result_set", [])
    
    # Extract the PDB IDs directly from the entry
    pdb_ids = [hit["identifier"] for hit in result_set]
    
    print(f"Broad net caught {len(pdb_ids)} PDB entries.")
    
    # Save the candidates to the CSV
    with output_csv.open("w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["pdb_id"])
        for pdb_id in pdb_ids:
            writer.writerow([pdb_id])
            
    print(f"Saved candidate list to {output_csv}.")


def main():
    parser = argparse.ArgumentParser(description="Fetch PDB metadata matching broad rule criteria.")
    
    parser.add_argument(
        "output_csv", 
        nargs="?", 
        default="data/metadata.csv", 
        help="Path to the output metadata CSV. Defaults to 'data/metadata.csv'"
    )
    
    parser.add_argument(
        "--min-length", 
        type=int, 
        default=4, 
        help="Minimum sequence length (default: 4)"
    )
    parser.add_argument(
        "--max-length", 
        type=int, 
        default=40, 
        help="Maximum sequence length (default: 40) - room for terminal caps"
    )
    
    args = parser.parse_args()
    
    output_path = Path(args.output_csv)
    fetch_pdb_metadata(output_path, args.min_length, args.max_length)

if __name__ == "__main__":
    main()
