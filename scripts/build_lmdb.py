import argparse
import csv
import gzip
import io
import sys
import zipfile
from collections import Counter
from pathlib import Path

from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.Data.IUPACData import protein_letters_3to1
from tqdm import tqdm

# Mapping from 3-letter to 1-letter AAs (Upper case)
AA_3TO1 = {k.upper(): v for k, v in protein_letters_3to1.items()}

# The 20 standard IUPAC Amino Acid 1-letter codes
STANDARD_AAS_1L = set("ARNDCEQGHILKMFPSTWYV")

def ensure_list(item):
    """Ensures the item is returned as a list even if it is a single string or None."""
    if item is None:
        return []
    if isinstance(item, str):
        return [item]
    return list(item)

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
            
            # Map Entity IDs to their Physical Chain IDs (Asym IDs)
            asym_ids_raw = mmcif_dict.get("_struct_asym.id", [])
            asym_entity_ids_raw = mmcif_dict.get("_struct_asym.entity_id", [])
            
            asym_ids = ensure_list(asym_ids_raw)
            asym_entity_ids = ensure_list(asym_entity_ids_raw)
            
            # Build entity_id -> [asym_id1, asym_id2, ...]
            entity_to_asym = {}
            for a_id, e_id in zip(asym_ids, asym_entity_ids):
                if e_id not in entity_to_asym:
                    entity_to_asym[e_id] = []
                entity_to_asym[e_id].append(a_id)
            
            # Map of asym_id -> observed sequence from seq_scheme
            # Extracts the standard residue sequence for each physical chain as defined in the CIF scheme.
            # Caps and solvent are not included in this scheme.
            seq_scheme_asyms = ensure_list(mmcif_dict.get("_pdbx_poly_seq_scheme.asym_id", []))
            seq_scheme_mons = ensure_list(mmcif_dict.get("_pdbx_poly_seq_scheme.mon_id", []))
            
            chain_sequences = {}
            for a_id, mon_id in zip(seq_scheme_asyms, seq_scheme_mons):
                if a_id not in chain_sequences:
                    chain_sequences[a_id] = []
                # Convert 3-letter to 1-letter, default to 'X' for non-standard
                chain_sequences[a_id].append(AA_3TO1.get(mon_id.upper(), "X"))

            # Extract entity polymer definitions
            entity_ids = ensure_list(mmcif_dict.get("_entity_poly.entity_id", []))
            types = ensure_list(mmcif_dict.get("_entity_poly.type", []))
            seqs = ensure_list(mmcif_dict.get("_entity_poly.pdbx_seq_one_letter_code_can", []))
            
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
                    
                # Identify physical chains for this entity
                target_asyms = entity_to_asym.get(e_id, [])
                if not target_asyms:
                    rejection_reasons["Entity ID not found in struct_asym"] += 1
                    continue
                
                for asym_id in target_asyms:
                    # 5. Exact Sequence Match
                    # Compare the physical chain's sequence (from seq_scheme) against the entity sequence.
                    # We tolerate structural gaps (missing coords), but the sequence identity 
                    # defined in the file must match the entity.
                    obs_seq_list = chain_sequences.get(asym_id, [])
                    obs_seq = "".join(obs_seq_list)
                    
                    if obs_seq != clean_seq:
                        rejection_reasons["Chain sequence mismatch vs Entity"] += 1
                        continue
                    
                    # Valid peptide chain!
                    valid_peptides.append({
                        "pdb_id": pdb_id,
                        "entity_id": e_id,
                        "chain": asym_id,
                        "sequence": clean_seq,
                        "length": seq_len
                    })
                    
                    # TODO: (Phase 2)
                    # Once we have validated that this specific chain matches our target, we will 
                    # resume parsing the heavy 3D coordinate block to identify its binding site 
                    # and save the refined backbone-only data to the LMDB database.
            
            pbar.set_postfix({"Found": len(valid_peptides)})

    # Print parsing statistics
    print("\n" + "="*50)
    print("PRE-FILTERING COMPLETE")
    print("="*50)
    print(f"Total entries processed: {total_entries}")
    print(f"Total valid peptide chains found ({min_len}-{max_len} AAs): {len(valid_peptides)}")
    
    print("\nTop Rejection Reasons:")
    for reason, count in rejection_reasons.most_common(10):
        print(f" - {reason}: {count}")
    print("-" * 50)
    
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    with csv_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["pdb_id", "entity_id", "chain", "sequence", "length"])
        writer.writeheader()
        writer.writerows(valid_peptides)
        
    print(f"\nSaved metadata tracking CSV to {csv_path}")

def main():
    parser = argparse.ArgumentParser(
        description="Boilerplate for parsing biological assemblies metadata to identify peptide targets before 3D structural analysis."
    )
    
    parser.add_argument(
        "zip_path", 
        nargs="?", 
        type=Path, 
        default="data/assemblies.zip",
        help="Path to the assemblies ZIP archive (default: data/assemblies.zip)."
    )
    parser.add_argument(
        "csv_path", 
        nargs="?", 
        type=Path, 
        default="data/peptide_counts.csv",
        help="Path to save the discovered peptide entities CSV (default: data/peptide_counts.csv)."
    )
    parser.add_argument("--min-length", type=int, default=4, help="Minimum peptide length (inclusive).")
    parser.add_argument("--max-length", type=int, default=32, help="Maximum peptide length (inclusive).")
    
    args = parser.parse_args()
    parse_assemblies(args.zip_path, args.csv_path, args.min_length, args.max_length)

if __name__ == "__main__":
    main()
