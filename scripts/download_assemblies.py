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

DEFAULT_METADATA_CSV = Path("data/metadata.csv")
DEFAULT_ASSEMBLIES_ZIP = Path("data/assemblies.zip")
DEFAULT_WORKERS = 8


def download_and_zip(
    pdb_id: str,
    session: requests.Session,
    zf: zipfile.ZipFile,
    zip_lock: threading.Lock,
    tmp_dir: str,
) -> str:
    """Download one assembly to a temp file and append it to the shared ZIP."""
    url = f"https://files.rcsb.org/download/{pdb_id.lower()}-assembly1.cif.gz"

    # The temporary file keeps memory bounded while many workers download in parallel.
    fd, temp_path = tempfile.mkstemp(dir=tmp_dir, suffix=".cif.gz")
    try:
        with os.fdopen(fd, "wb") as handle:
            response = session.get(url, timeout=20, stream=True)
            if response.status_code == 404:
                raise ValueError("Assembly not found (404)")
            response.raise_for_status()

            for chunk in response.iter_content(chunk_size=1024 * 1024):
                if chunk:
                    handle.write(chunk)

        file_name = f"{pdb_id.lower()}-assembly1.cif.gz"

        # ZIP writes are serialized because ZipFile is shared across workers.
        with zip_lock:
            zf.write(temp_path, arcname=file_name)
    finally:
        try:
            os.remove(temp_path)
        except OSError:
            pass

    return pdb_id


def read_pdb_ids(metadata_csv: Path) -> list[str]:
    """Read PDB IDs from the metadata CSV produced by fetch_metadata."""
    pdb_ids: list[str] = []
    with metadata_csv.open("r", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            if "pdb_id" in row:
                pdb_ids.append(row["pdb_id"])
    return pdb_ids


def validate_parameters(metadata_csv: Path, output_zip: Path, workers: int) -> None:
    """Validate CLI parameters before network downloads start."""
    if not metadata_csv.exists():
        raise ValueError(f"{metadata_csv} does not exist. Run fetch_metadata.py first.")
    if metadata_csv.is_dir():
        raise ValueError("metadata_csv must be a file path, not a directory")
    if output_zip.exists() and output_zip.is_dir():
        raise ValueError("output_zip must be a file path, not a directory")
    if workers < 1:
        raise ValueError("--workers must be at least 1")


def download_assemblies(metadata_csv: Path, output_zip: Path, workers: int) -> None:
    """Download first biological assemblies into a single ZIP archive."""
    pdb_ids = read_pdb_ids(metadata_csv)

    total_pdbs = len(pdb_ids)
    if total_pdbs == 0:
        print("No PDB IDs found in the CSV. Exiting.")
        return

    print(f"Found {total_pdbs} PDB entries to process.")

    output_zip.parent.mkdir(parents=True, exist_ok=True)

    # Retries cover transient static-file server failures without changing skip semantics.
    session = requests.Session()
    retries = Retry(total=3, backoff_factor=0.5, status_forcelist=[500, 502, 503, 504])
    session.mount(
        "https://",
        HTTPAdapter(
            pool_connections=workers, pool_maxsize=workers, max_retries=retries
        ),
    )

    zip_lock = threading.Lock()
    success_count = 0
    fail_count = 0

    mode = "a" if output_zip.exists() else "w"

    with zipfile.ZipFile(output_zip, mode, zipfile.ZIP_STORED) as zf:
        existing_files = set(zf.namelist())
        pending_pdb_ids = [
            pid
            for pid in pdb_ids
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

        with tempfile.TemporaryDirectory() as tmp_dir:
            with tqdm(
                total=pending_total, desc="Downloading", unit="file", dynamic_ncols=True
            ) as pbar:
                with ThreadPoolExecutor(max_workers=workers) as executor:
                    future_to_pdb = {
                        executor.submit(
                            download_and_zip, pid, session, zf, zip_lock, tmp_dir
                        ): pid
                        for pid in pending_pdb_ids
                    }

                    for future in as_completed(future_to_pdb):
                        pdb_id = future_to_pdb[future]
                        try:
                            future.result()
                            success_count += 1

                        except Exception as exc:
                            fail_count += 1
                            pbar.write(f"[WARNING] {pdb_id}: {exc}")

                        finally:
                            pbar.set_postfix({"Failed/Missing": fail_count})
                            pbar.update(1)

    print(
        f"\nFinished! Total downloaded: {success_count - skipped}, Missing/Failed: {fail_count}."
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Download leading biological assemblies for PDB entries."
    )

    parser.add_argument(
        "metadata_csv",
        nargs="?",
        type=Path,
        default=DEFAULT_METADATA_CSV,
        help="Path to the input metadata CSV (default: data/metadata.csv)",
    )

    parser.add_argument(
        "output_zip",
        nargs="?",
        type=Path,
        default=DEFAULT_ASSEMBLIES_ZIP,
        help="Path to the output ZIP (default: data/assemblies.zip)",
    )

    parser.add_argument(
        "--workers",
        type=int,
        default=DEFAULT_WORKERS,
        help="Number of parallel download workers (default: 8)",
    )

    args = parser.parse_args()
    try:
        validate_parameters(args.metadata_csv, args.output_zip, args.workers)
        download_assemblies(args.metadata_csv, args.output_zip, args.workers)
    except ValueError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
