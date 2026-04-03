import argparse
import csv
import gzip
import io
import sys
import zipfile
from collections import Counter
from pathlib import Path

from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from tqdm import tqdm

# The 20 standard IUPAC Amino Acid 1-letter codes
STANDARD_AAS_1L = set("ARNDCEQGHILKMFPSTWYV")

def header_only_lines(f):
    """
    Yields lines from the CIF file only up to the massive atomic coordinate block.
    This cuts parsing time by ~99% because we only care about the entity metadata header.
    """
    for line in f:
        if line.startswith("_atom_site"):
            break
        yield line

def parse_assemblies(zip_path: Path, csv_path: Path, min_len: int, max_len: int):
    """
    Metadata-only pre-filter boilerplate.
    Identifies polypeptide entities meeting length and sequence criteria 
    before heavy 3D coordinate parsing and LMDB construction.
    """
    if not zip_path.exists():
        print(f"Error: {zip_path} does not exist.", file=sys.stderr)
        sys.exit(1)

    valid_peptides = []
    rejection_reasons = Counter()
    
    with zipfile.ZipFile(zip_path, "r") as zf:
        cif_files = [n for n in zf.namelist() if n.endswith(".cif.gz")]
        total_entries = len(cif_files)
        
        print(f"Found {total_entries} CIF files in the archive.")
        
        pbar = tqdm(cif_files, desc="Filtering Entities", unit="file", dynamic_ncols=True)
        for file_name in pbar:
            pdb_id = file_name.split("-")[0].upper()
            compressed_bytes = zf.read(file_name)
            
            with gzip.open(io.BytesIO(compressed_bytes), mode="rt") as f:
                try:
                    # Provide an early-stopping iterator to MMCIF2Dict so it parses ONLY the metadata header
                    mmcif_dict = MMCIF2Dict(header_only_lines(f))
                except Exception:
                    rejection_reasons["Parse Error"] += 1
                    continue
            
            # Extract entity polymer definitions
            entity_ids_raw = mmcif_dict.get("_entity_poly.entity_id", [])
            types_raw = mmcif_dict.get("_entity_poly.type", [])
            seqs_raw = mmcif_dict.get("_entity_poly.pdbx_seq_one_letter_code_can", [])
            
            # In MMCIF2Dict, if there's only one item, it might be returned as a string rather than a list
            entity_ids = [entity_ids_raw] if isinstance(entity_ids_raw, str) else entity_ids_raw
            types = [types_raw] if isinstance(types_raw, str) else types_raw
            seqs = [seqs_raw] if isinstance(seqs_raw, str) else seqs_raw
            
            if not entity_ids:
                rejection_reasons["No polymeric entities"] += 1
                continue
                
            for e_id, e_type, e_seq in zip(entity_ids, types, seqs):
                # 1. Must be a protein polypeptide
                if not e_type or "polypeptide" not in e_type.lower():
                    rejection_reasons[f"Non-polypeptide type"] += 1
                    continue
                    
                # 2. Extract clean sequence
                clean_seq = e_seq.replace("\n", "").replace(" ", "").replace("\r", "")
                seq_len = len(clean_seq)
                
                # 3. Enforce target length bounds
                if not (min_len <= seq_len <= max_len):
                    rejection_reasons["Sequence length out of bounds"] += 1
                    continue
                    
                # 4. Strict 20 Standard AAs orthography validation
                if not set(clean_seq).issubset(STANDARD_AAS_1L):
                    rejection_reasons["Contains non-standard structural elements or 'X'"] += 1
                    continue
                    
                # Valid peptide entity!
                valid_peptides.append({
                    "pdb_id": pdb_id,
                    "entity_id": e_id,
                    "sequence": clean_seq,
                    "length": seq_len
                })
                
                # TODO: 
                # Once we have validated this entity is a target peptide, we will resume parsing 
                # the heavy 3D coordinates block for its associated physical chains to compute 
                # binding interfaces and save them to the LMDB database.
            
            pbar.set_postfix({"Found": len(valid_peptides)})

    # Print parsing statistics
    print("\n" + "="*50)
    print("PRE-FILTERING COMPLETE")
    print("="*50)
    print(f"Total entries processed: {total_entries}")
    print(f"Total valid peptide entities found ({min_len}-{max_len} AAs): {len(valid_peptides)}")
    
    print("\nTop Rejection Reasons:")
    for reason, count in rejection_reasons.most_common(10):
        print(f" - {reason}: {count}")
    print("-" * 50)
    
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    with csv_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["pdb_id", "entity_id", "sequence", "length"])
        writer.writeheader()
        writer.writerows(valid_peptides)
        
    print(f"\nSaved metadata tracking CSV to {csv_path}")

def main():
    parser = argparse.ArgumentParser(
        description="Boilerplate for parsing biological assemblies metadata to identify peptide targets before 3D structural analysis."
    )
    
    parser.add_argument("zip_path", type=Path, help="Path to the assemblies ZIP archive.")
    parser.add_argument("csv_path", type=Path, help="Path to save the discovered peptide entities CSV.")
    parser.add_argument("--min-length", type=int, default=4, help="Minimum peptide length (inclusive).")
    parser.add_argument("--max-length", type=int, default=32, help="Maximum peptide length (inclusive).")
    
    args = parser.parse_args()
    parse_assemblies(args.zip_path, args.csv_path, args.min_length, args.max_length)

if __name__ == "__main__":
    main()
