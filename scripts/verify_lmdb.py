from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import cast

import lmdb
import numpy as np

from pdb_mldata.lmdb_utils import decode_lmdb_entry

DEFAULT_LMDB_PATH = Path("data/pdb_mldata.lmdb")


def json_default(obj: object) -> str:
    """Convert unknown JSON values to strings for defensive debug output."""
    if isinstance(obj, np.ndarray):
        return f"<shape {obj.shape}>"
    return str(obj)


def summarize_structures(data: object) -> object:
    """Recursively replace large structural arrays with shape summaries."""
    if isinstance(data, dict):
        new_data: dict[object, object] = {}
        for key, value in data.items():
            if key in ["structure", "b_factors", "occupancy"] and isinstance(
                value, np.ndarray
            ):
                new_data[key] = f"<shape {value.shape}>"
            else:
                new_data[key] = summarize_structures(value)
        return new_data
    if isinstance(data, list):
        return [summarize_structures(item) for item in data]
    return data


def validate_parameters(db_path: Path) -> None:
    """Validate CLI parameters before opening the LMDB environment."""
    if not db_path.exists():
        raise ValueError(f"{db_path} not found")
    if not db_path.is_dir():
        raise ValueError("db_path must be an LMDB directory path")


def verify_lmdb(db_path: Path) -> None:
    """Decode and pretty-print entries from the database, excluding full structural data."""
    env = lmdb.open(str(db_path), readonly=True)
    with env.begin() as txn:
        for key, value in txn.cursor():
            raw_data = decode_lmdb_entry(cast(bytes, value))
            clean_data = summarize_structures(raw_data)
            print(json.dumps(clean_data, indent=2))
    env.close()


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Verify LMDB database entries by dumping metadata and structural shapes."
    )
    parser.add_argument(
        "db_path",
        nargs="?",
        type=Path,
        default=DEFAULT_LMDB_PATH,
        help="Path to the LMDB database directory (default: data/pdb_mldata.lmdb)",
    )
    args = parser.parse_args()

    try:
        validate_parameters(args.db_path)
    except ValueError as exc:
        parser.error(str(exc))

    verify_lmdb(args.db_path)


if __name__ == "__main__":
    main()
