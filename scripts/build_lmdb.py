import argparse
import gzip
import io
import sys
import zipfile
from collections import Counter
from pathlib import Path

from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.Polypeptide import protein_letters_3to1
from tqdm import tqdm

STANDARD_AAS = {
    'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 
    'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'
}

def parse_assemblies(zip_path: Path, min_len: int, max_len: int):
    """
    Parses a ZIP file containing .cif.gz biological assemblies using BioPython.
    Identifies and counts continuous peptide chains composed entirely of standard amino acids.
    """
    if not zip_path.exists():
        print(f"Error: {zip_path} does not exist.", file=sys.stderr)
        sys.exit(1)

    # Disable PDBConstructionWarnings to prevent log flooding from poorly formatted CIFs
    parser = MMCIFParser(QUIET=True)
    
    valid_peptides = []
    skipped_multi_model = 0
    rejection_reasons = Counter()
    
    with zipfile.ZipFile(zip_path, "r") as zf:
        cif_files = [n for n in zf.namelist() if n.endswith(".cif.gz")]
        total_entries = len(cif_files)
        
        print(f"Found {total_entries} CIF files in the archive.")
        
        pbar = tqdm(cif_files, desc="Parsing Assemblies", unit="file", dynamic_ncols=True)
        for file_name in pbar:
            pdb_id = file_name.split("-")[0].upper()
            compressed_bytes = zf.read(file_name)
            
            with gzip.open(io.BytesIO(compressed_bytes), mode="rt") as f:
                try:
                    structure = parser.get_structure(pdb_id, f)
                except Exception:
                    rejection_reasons["Parse Error"] += 1
                    continue
            
            # Defensive Rule: Skip complex multi-model NMR structures
            if len(structure) > 1:
                skipped_multi_model += 1
                rejection_reasons["Multiple Models"] += 1
                continue
                
            model = structure[0]
            for chain in model:
                chain_id = chain.get_id()
                is_valid_peptide = True
                seq_symbols = []
                
                residues = list(chain.get_residues())
                if not residues:
                    continue
                
                last_res_num = None
                for res in residues:
                    res_id = res.get_id()
                    res_name = res.get_resname().strip()
                    
                    # 1. Ignore solvent shells
                    if res_id[0] == "W" or res_name == "HOH":
                        continue
                        
                    # 2. Skip non-standard residues (e.g., terminal caps, synthetic modifications)
                    # Skipping terminal caps leaves the standard peptide core perfectly continuous.
                    # Skipping internal modifications intentionally causes a continuity gap below.
                    if res_name not in STANDARD_AAS:
                        rejection_reasons[f"Skipped Non-standard ({res_name})"] += 1
                        continue
                        
                    res_num = res_id[1]
                    ins_code = res_id[2]
                    
                    # 3. Reject non-standard insertion numbering
                    if ins_code != " ":
                        rejection_reasons["Insertion code attached"] += 1
                        is_valid_peptide = False
                        break
                        
                    # 4. Enforce strict continuity to ensure the sequence isn't physically broken
                    if last_res_num is not None and res_num != last_res_num + 1:
                        rejection_reasons["Sequence numbering gap"] += 1
                        is_valid_peptide = False
                        break
                    
                    last_res_num = res_num
                    
                    # 5. Extract 1-letter sequence representation
                    try:
                        symbol = protein_letters_3to1[res_name]
                        seq_symbols.append(symbol)
                    except KeyError:
                        rejection_reasons["KeyError on 3to1 map"] += 1
                        is_valid_peptide = False
                        break
                
                # 6. Finalize validity based on defined pipeline length constraints
                if is_valid_peptide:
                    chain_len = len(seq_symbols)
                    if min_len <= chain_len <= max_len:
                        sequence = "".join(seq_symbols)
                        valid_peptides.append({
                            "pdb_id": pdb_id,
                            "chain": chain_id,
                            "sequence": sequence,
                            "length": chain_len
                        })
            
            pbar.set_postfix({"Peptides": len(valid_peptides)})

    # Print parsing statistics
    print("\n" + "="*50)
    print("PARSING COMPLETE")
    print("="*50)
    print(f"Total entries processed: {total_entries}")
    print(f"Skipped multi-model entries: {skipped_multi_model}")
    print(f"Total peptide chains found ({min_len}-{max_len} AAs): {len(valid_peptides)}")
    
    print("\nTop 10 Rejection Reasons:")
    for reason, count in rejection_reasons.most_common(10):
        print(f" - {reason}: {count}")
    print("-" * 50)
    
    # Sequence preview limit to avoid UI freezing
    print("\nDiscovered Peptides:")
    for i, p in enumerate(valid_peptides):
        print(f"- {p['pdb_id']}_{p['chain']}: {p['sequence']} ({p['length']} residues)")
        if i >= 99:
            print(f"... and {len(valid_peptides) - 100} more omitted from console output.")
            break

def main():
    parser = argparse.ArgumentParser(
        description="Phase 1: Parse biological assemblies from ZIP and count peptide chains within specific lengths (Strict Amino Acids only)."
    )
    
    parser.add_argument(
        "zip_path", 
        nargs="?",
        type=Path, 
        default="data/assemblies.zip",
        help="Path to the assemblies.zip file (default: data/assemblies.zip)"
    )
    
    parser.add_argument(
        "--min-length", 
        type=int, 
        default=4, 
        help="Minimum inclusive length of a peptide chain (default: 4)"
    )
    
    parser.add_argument(
        "--max-length", 
        type=int, 
        default=32, 
        help="Maximum inclusive length of a peptide chain (default: 32)"
    )
    
    args = parser.parse_args()
    parse_assemblies(args.zip_path, args.min_length, args.max_length)

if __name__ == "__main__":
    main()
