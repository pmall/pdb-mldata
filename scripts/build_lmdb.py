import argparse
import csv
import gzip
import io
import sys
import zipfile
from collections import Counter
from pathlib import Path

import numpy as np
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.Data.IUPACData import protein_letters_3to1
from tqdm import tqdm

# Mapping from 3-letter to 1-letter AAs (Upper case)
AA_3TO1 = {k.upper(): v for k, v in protein_letters_3to1.items()}

# The 20 standard IUPAC Amino Acid 1-letter codes
STANDARD_AAS_1L = set("ARNDCEQGHILKMFPSTWYV")

# Canonical mapping of 37 atoms commonly used in protein ML conventions.
# This ensures a fixed-size representation (L x 37 x 3) for each peptide.
ATOM37_LIST = [
    "N", "CA", "C", "CB", "O", "CG", "CG1", "CG2", "OG", "OG1", "SG", "CD",
    "CD1", "CD2", "ND1", "ND2", "OD1", "OD2", "SD", "CE", "CE1", "CE2", "CE3",
    "NE", "NE1", "NE2", "OE1", "OE2", "CH2", "NH1", "NH2", "OH", "CZ", "CZ2",
    "CZ3", "NZ", "OXT"
]
ATOM37_MAP = {name: i for i, name in enumerate(ATOM37_LIST)}

# Number of peptides to display and process in verification mode
VERIFY_COUNT = 50


def ensure_list(item):
    """Ensures the item is returned as a list even if it is a single string or None."""
    if item is None:
        return []
    if isinstance(item, str):
        return [item]
    return list(item)


def verify_peptide_extraction(peptide_info, chain_atom_map, raw_seq, start_seq_id):
    """
    Prints a detailed breakdown of the extracted peptide for verification.
    """
    pdb_id = peptide_info["pdb_id"]
    asym_id = peptide_info["chain"]
    trimmed_seq = peptide_info["sequence"]
    structure_coords = peptide_info["structure"]
    
    # Check if trimming occurred
    trimmed = start_seq_id > 1 or len(trimmed_seq) < len(raw_seq)
    
    print(f"[VERIFY] PDB: {pdb_id} | Chain: {asym_id}")
    if trimmed:
        print(f"  *** TRIMMING DETECTED ***")
        print(f"  Raw Seq: {raw_seq}")
        
    print(f"  Sequence: {trimmed_seq}")
    print(f"  Range: {start_seq_id} to {start_seq_id + len(trimmed_seq) - 1}")
    
    # Show breakdown for each residue
    for i, res_char in enumerate(trimmed_seq):
        res_atoms = chain_atom_map.get(asym_id, {}).get(str(start_seq_id + i), [])
        names_found = [a[0] for a in res_atoms if a[0] in ATOM37_MAP]
        
        is_gap = len(names_found) == 0
        status = "GAP" if is_gap else f"{len(names_found)} atoms"
        print(f"    Res {i+1} ({res_char}): {status} | {','.join(names_found[:10])}{'...' if len(names_found) > 10 else ''}")
    
    found_atoms = sum(1 for res in structure_coords for atom in res if not np.isnan(atom[0]))
    print(f"  Total Atoms Extracted: {found_atoms}")


def parse_assemblies(zip_path: Path, csv_path: Path, min_len: int, max_len: int, verify_mode: bool = False):
    """
    Parses biological assemblies to extract valid peptides and their 3D coordinates.
    Aligns the structural residues with the trimmed entity sequence.
    """
    if not zip_path.exists():
        print(f"Error: {zip_path} does not exist.", file=sys.stderr)
        sys.exit(1)

    valid_peptides = []
    rejection_reasons = Counter()
    
    with zipfile.ZipFile(zip_path, "r") as zf:
        cif_files = [n for n in zf.namelist() if n.endswith(".cif.gz")]
        total_entries = len(cif_files)
        
        if not verify_mode:
            print(f"Found {total_entries} CIF files in the archive.")
        
        pbar = tqdm(
            cif_files, 
            desc="Processing Assemblies", 
            unit="file", 
            dynamic_ncols=True, 
            disable=verify_mode
        )

        for file_name in pbar:
            pdb_id = file_name.split("-")[0].upper()
            compressed_bytes = zf.read(file_name)
            
            with gzip.open(io.BytesIO(compressed_bytes), mode="rt") as f:
                try:
                    mmcif_dict = MMCIF2Dict(f)
                except Exception:
                    rejection_reasons["Parse Error"] += 1
                    continue
            
            # 1. Map Entity IDs to their Physical Chain IDs (Asym IDs)
            asym_ids = ensure_list(mmcif_dict.get("_struct_asym.id", []))
            asym_entity_ids = ensure_list(mmcif_dict.get("_struct_asym.entity_id", []))
            
            entity_to_asym = {}
            for a_id, e_id in zip(asym_ids, asym_entity_ids):
                if e_id not in entity_to_asym:
                    entity_to_asym[e_id] = []
                entity_to_asym[e_id].append(a_id)

            # 2. Extract entity polymer definitions
            ent_poly_ids = ensure_list(mmcif_dict.get("_entity_poly.entity_id", []))
            ent_poly_types = ensure_list(mmcif_dict.get("_entity_poly.type", []))
            ent_poly_seqs = ensure_list(mmcif_dict.get("_entity_poly.pdbx_seq_one_letter_code_can", []))
            
            if not ent_poly_ids:
                rejection_reasons["No polymeric entities"] += 1
                continue
            
            # 3. Organize atomic Data Blocks
            atom_asym_id = ensure_list(mmcif_dict.get("_atom_site.label_asym_id", []))
            atom_seq_id = ensure_list(mmcif_dict.get("_atom_site.label_seq_id", []))
            atom_name = ensure_list(mmcif_dict.get("_atom_site.label_atom_id", []))
            atom_x = ensure_list(mmcif_dict.get("_atom_site.Cartn_x", []))
            atom_y = ensure_list(mmcif_dict.get("_atom_site.Cartn_y", []))
            atom_z = ensure_list(mmcif_dict.get("_atom_site.Cartn_z", []))

            # Group atoms by (asym_id, label_seq_id) for efficient lookup
            chain_atom_map = {}
            for a_id, s_id, name, x, y, z in zip(atom_asym_id, atom_seq_id, atom_name, atom_x, atom_y, atom_z):
                if a_id not in chain_atom_map:
                    chain_atom_map[a_id] = {}
                if s_id not in chain_atom_map[a_id]:
                    chain_atom_map[a_id][s_id] = []
                chain_atom_map[a_id][s_id].append((name, float(x), float(y), float(z)))

            # 4. Iterate over entities to identify peptides
            for e_id, e_type, e_seq in zip(ent_poly_ids, ent_poly_types, ent_poly_seqs):
                if not e_type or "polypeptide" not in e_type.lower():
                    continue
                    
                # Clean and Trim the theoretical sequence
                raw_seq = e_seq.replace("\n", "").replace(" ", "").replace("\r", "")
                non_std_chars = "".join(set(raw_seq) - STANDARD_AAS_1L)
                trimmed_seq = raw_seq.strip(non_std_chars)
                seq_len = len(trimmed_seq)
                
                if not (min_len <= seq_len <= max_len) or not set(trimmed_seq).issubset(STANDARD_AAS_1L):
                    continue
                
                # Calculate start/end offsets (1-indexed for CIF label_seq_id)
                prefix_len = raw_seq.find(trimmed_seq)
                start_seq_id = prefix_len + 1
                end_seq_id = start_seq_id + seq_len - 1
                
                # Identify physical chains (asym_ids) mapping to this entity
                target_asyms = entity_to_asym.get(e_id, [])
                for asym_id in target_asyms:
                    
                    # 5. Extract structure aligned with trimmed sequence
                    structure_coords = []
                    for current_id in range(start_seq_id, end_seq_id + 1):
                        residue_coords = np.full((37, 3), np.nan)
                        residue_atoms = chain_atom_map.get(asym_id, {}).get(str(current_id), [])
                        for name, x, y, z in residue_atoms:
                            if name in ATOM37_MAP:
                                residue_coords[ATOM37_MAP[name]] = [x, y, z]
                        structure_coords.append(residue_coords.tolist())
                    
                    peptide_info = {
                        "pdb_id": pdb_id,
                        "entity_id": e_id,
                        "chain": asym_id,
                        "sequence": trimmed_seq,
                        "length": seq_len,
                        "structure": structure_coords
                    }
                    valid_peptides.append(peptide_info)
                    
                    # Run verification if requested
                    if verify_mode and len(valid_peptides) <= VERIFY_COUNT:
                        verify_peptide_extraction(peptide_info, chain_atom_map, raw_seq, start_seq_id)
                    
                    # Early exit for verification mode
                    if verify_mode and len(valid_peptides) >= VERIFY_COUNT:
                        return valid_peptides
                    
            if not verify_mode:
                pbar.set_postfix({"Found": len(valid_peptides)})

    # Final Summary and Output for Normal Mode
    if not verify_mode:
        print(f"\nTotal valid peptides encountered: {len(valid_peptides)}")
        
        csv_path.parent.mkdir(parents=True, exist_ok=True)
        with csv_path.open("w", newline="") as f:
            # We only write metadata to CSV (no structure)
            fieldnames = ["pdb_id", "entity_id", "chain", "sequence", "length"]
            writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction='ignore')
            writer.writeheader()
            writer.writerows(valid_peptides)
        
        print(f"Saved metadata tracking CSV to {csv_path}")
    
    return valid_peptides


def main():
    parser = argparse.ArgumentParser(
        description="Extract aligned 3D coordinates for peptide chains and save metadata."
    )
    
    parser.add_argument(
        "zip_path", 
        nargs="?", 
        type=Path, 
        default="data/assemblies.zip",
        help="Path to the assemblies ZIP archive."
    )
    parser.add_argument(
        "csv_path", 
        nargs="?", 
        type=Path, 
        default="data/peptide_counts.csv",
        help="Path to save the discovered peptide entities CSV."
    )
    parser.add_argument("--min-length", type=int, default=4, help="Minimum peptide length.")
    parser.add_argument("--max-length", type=int, default=32, help="Maximum peptide length.")
    parser.add_argument("--verify", action="store_true", help="Print verification results and exit early.")
    
    args = parser.parse_args()
    
    # Run the parsing logic
    parse_assemblies(args.zip_path, args.csv_path, args.min_length, args.max_length, verify_mode=args.verify)


if __name__ == "__main__":
    main()
