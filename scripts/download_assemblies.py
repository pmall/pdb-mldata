import argparse
import csv
import os
import sys
import tempfile
import threading
import zipfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import requests
from requests.adapters import HTTPAdapter
from tqdm import tqdm
from urllib3.util.retry import Retry

def download_and_zip(pdb_id: str, session: requests.Session, zf: zipfile.ZipFile, zip_lock: threading.Lock, tmp_dir: str) -> str:
    """Downloads the assembly to a temp file and safely appends it to the zip. This avoids RAM bloat."""
    
    # We use the official static download server, which is extremely robust and fast.
    url = f"https://files.rcsb.org/download/{pdb_id.lower()}-assembly1.cif.gz"
    
    # We use a temporary file to avoid loading the entire file into RAM.
    # This prevents Out Of Memory (OOM) freezes.
    fd, temp_path = tempfile.mkstemp(dir=tmp_dir, suffix=".cif.gz")
    try:
        with os.fdopen(fd, "wb") as f:
            response = session.get(url, timeout=20, stream=True)
            if response.status_code == 404:
                raise ValueError("Assembly not found (404)")
            response.raise_for_status()
            
            # Stream directly to disk in chunks to strictly bound memory usage
            for chunk in response.iter_content(chunk_size=1024 * 1024):
                if chunk:
                    f.write(chunk)
                    
        file_name = f"{pdb_id.lower()}-assembly1.cif.gz"
        
        # Write to the zip file in a thread-safe manner.
        # This implicitly creates backpressure: a worker cannot start
        # a new download until its current download is safely committed to the ZIP.
        with zip_lock:
            zf.write(temp_path, arcname=file_name)
    finally:
        try:
            os.remove(temp_path)
        except OSError:
            pass
            
    return pdb_id

def download_assemblies(metadata_csv: Path, output_zip: Path, workers: int):
    """
    Reads PDB IDs from the CSV and downloads their first assemblies in parallel,
    saving the CIF files directly into a single ZIP archive.
    """
    if not metadata_csv.exists():
        print(f"Error: {metadata_csv} does not exist. Run fetch_metadata.py first.", file=sys.stderr)
        sys.exit(1)
        
    pdb_ids = []
    with metadata_csv.open("r", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            if "pdb_id" in row:
                pdb_ids.append(row["pdb_id"])
                
    total_pdbs = len(pdb_ids)
    if total_pdbs == 0:
        print("No PDB IDs found in the CSV. Exiting.")
        sys.exit(0)
        
    print(f"Found {total_pdbs} PDB entries to process.")
    
    output_zip.parent.mkdir(parents=True, exist_ok=True)
    
    # Robust static file connection pool
    session = requests.Session()
    retries = Retry(total=3, backoff_factor=0.5, status_forcelist=[500, 502, 503, 504])
    session.mount("https://", HTTPAdapter(pool_connections=workers, pool_maxsize=workers, max_retries=retries))
    
    zip_lock = threading.Lock()
    success_count = 0
    fail_count = 0
    
    mode = "a" if output_zip.exists() else "w"
    
    with zipfile.ZipFile(output_zip, mode, zipfile.ZIP_STORED) as zf:
        
        existing_files = set(zf.namelist())
        pending_pdb_ids = [
            pid for pid in pdb_ids 
            if f"{pid.lower()}-assembly1.cif.gz" not in existing_files
        ]
        
        skipped = total_pdbs - len(pending_pdb_ids)
        if skipped > 0:
            print(f"Skipping {skipped} assemblies already present in the ZIP archive.")
            success_count += skipped
            
        pending_total = len(pending_pdb_ids)
        if pending_total == 0:
            print("All entries are already downloaded!")
            return
            
        print(f"Starting parallel downloads. Saving strictly to {output_zip}...")

        # Using tqdm for a production-grade, bug-free progress tracker
        with tempfile.TemporaryDirectory() as tmp_dir:
            with tqdm(total=pending_total, desc="Downloading", unit="file", dynamic_ncols=True) as pbar:
                with ThreadPoolExecutor(max_workers=workers) as executor:
                    future_to_pdb = {
                        executor.submit(download_and_zip, pid, session, zf, zip_lock, tmp_dir): pid 
                        for pid in pending_pdb_ids
                    }
                    
                    for future in as_completed(future_to_pdb):
                        pdb_id = future_to_pdb[future]
                        try:
                            # If this succeeds, the file is already downloaded and zipped via the worker
                            future.result()
                            success_count += 1
                            
                        except Exception as e:
                            fail_count += 1
                            pbar.write(f"[WARNING] {pdb_id}: {e}")
                            
                        finally:
                            pbar.set_postfix({"Failed/Missing": fail_count})
                            pbar.update(1)

    print(f"\nFinished! Total downloaded: {success_count - skipped}, Missing/Failed: {fail_count}.")

def main():
    parser = argparse.ArgumentParser(description="Download leading biological assemblies for PDB entries.")
    
    parser.add_argument(
        "metadata_csv", 
        nargs="?", 
        default="data/metadata.csv", 
        help="Path to the input metadata CSV. Defaults to 'data/metadata.csv'"
    )
    
    parser.add_argument(
        "output_zip", 
        nargs="?", 
        default="data/assemblies.zip", 
        help="Path to the output ZIP. Defaults to 'data/assemblies.zip'"
    )
    
    parser.add_argument(
        "--workers", 
        type=int, 
        default=8, 
        help="Number of parallel download workers (default: 8)"
    )
    
    args = parser.parse_args()
    
    metadata_path = Path(args.metadata_csv)
    zip_path = Path(args.output_zip)
    
    download_assemblies(metadata_path, zip_path, args.workers)

if __name__ == "__main__":
    main()
