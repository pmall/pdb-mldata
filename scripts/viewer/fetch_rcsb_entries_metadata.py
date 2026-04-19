from __future__ import annotations

import argparse
import sqlite3
from dataclasses import dataclass
from pathlib import Path

from tqdm import tqdm

from pdb_mldata.viewer_rcsb import (
    DEFAULT_BATCH_SIZE,
    ENTRY_METADATA_QUERY,
    batched,
    configure_session,
    extract_entry_metadata,
    fetch_graphql,
    first_string,
)

DEFAULT_SQLITE_PATH = Path("data/pdb_mldata.sqlite")


@dataclass(frozen=True)
class EntryMetadataStats:
    entries_inserted: int
    failures: int


def validate_parameters(sqlite_path: Path, batch_size: int) -> None:
    """Validate CLI parameters before opening the viewer database."""
    if not sqlite_path.exists():
        raise ValueError(f"{sqlite_path} not found")
    if sqlite_path.is_dir():
        raise ValueError("sqlite_path must be a file path, not a directory")
    if batch_size < 1:
        raise ValueError("--batch-size must be at least 1")


def rebuild_entry_metadata_schema(connection: sqlite3.Connection) -> None:
    """Delete only entry metadata before rebuilding it from selected rows."""
    connection.executescript(
        """
        PRAGMA foreign_keys = ON;

        DROP TABLE IF EXISTS search_terms_targets;
        DROP TABLE IF EXISTS search_terms;
        DROP TABLE IF EXISTS entries_metadata;

        CREATE TABLE entries_metadata (
            pdb_id TEXT PRIMARY KEY,
            title TEXT NOT NULL,
            experimental_methods_json TEXT NOT NULL,
            resolution_combined_json TEXT NOT NULL,
            best_resolution REAL,
            deposition_date TEXT,
            initial_release_date TEXT,
            FOREIGN KEY (pdb_id) REFERENCES entries (pdb_id)
        );
        """
    )
    connection.commit()


def read_entry_ids(connection: sqlite3.Connection) -> list[str]:
    return [
        row[0]
        for row in connection.execute(
            "SELECT pdb_id FROM entries ORDER BY pdb_id"
        ).fetchall()
    ]


def insert_row(
    connection: sqlite3.Connection,
    table: str,
    columns: tuple[str, ...],
    values: tuple[object, ...],
) -> None:
    placeholders = ", ".join("?" for _ in columns)
    sql = f"INSERT INTO {table} ({', '.join(columns)}) VALUES ({placeholders})"
    connection.execute(sql, values)


def fetch_entry_metadata(
    connection: sqlite3.Connection,
    pdb_ids: list[str],
    batch_size: int,
) -> EntryMetadataStats:
    columns = (
        "pdb_id",
        "title",
        "experimental_methods_json",
        "resolution_combined_json",
        "best_resolution",
        "deposition_date",
        "initial_release_date",
    )
    session = configure_session()
    inserted = 0
    failures = 0

    progress = tqdm(
        batched(pdb_ids, batch_size),
        desc="Fetching RCSB entry metadata",
        unit="batch",
    )
    for pdb_id_batch in progress:
        try:
            payload = fetch_graphql(
                session=session,
                query=ENTRY_METADATA_QUERY,
                variables={"entryIds": pdb_id_batch},
            )
        except Exception:
            failures += len(pdb_id_batch)
            continue

        rows = payload.get("entries")
        if not isinstance(rows, list):
            failures += len(pdb_id_batch)
            continue

        returned_ids: set[str] = set()
        for row in rows:
            if not isinstance(row, dict):
                failures += 1
                continue
            pdb_id = first_string(row, "rcsb_id")
            if not pdb_id:
                failures += 1
                continue
            try:
                insert_row(
                    connection=connection,
                    table="entries_metadata",
                    columns=columns,
                    values=extract_entry_metadata(pdb_id=pdb_id, data=row),
                )
                inserted += 1
                returned_ids.add(pdb_id)
            except Exception:
                failures += 1

        failures += len(set(pdb_id_batch) - returned_ids)

    return EntryMetadataStats(entries_inserted=inserted, failures=failures)


def fetch_rcsb_entries_metadata(
    sqlite_path: Path, batch_size: int
) -> EntryMetadataStats:
    connection = sqlite3.connect(sqlite_path)
    rebuild_entry_metadata_schema(connection)
    entry_ids = read_entry_ids(connection)
    stats = fetch_entry_metadata(
        connection=connection,
        pdb_ids=entry_ids,
        batch_size=batch_size,
    )
    connection.commit()
    connection.close()
    return stats


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Fetch RCSB entry metadata for selected viewer database rows."
    )
    parser.add_argument(
        "sqlite_path",
        nargs="?",
        type=Path,
        default=DEFAULT_SQLITE_PATH,
        help="Path to the SQLite database (default: data/pdb_mldata.sqlite)",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=DEFAULT_BATCH_SIZE,
        help="Number of RCSB GraphQL IDs per request (default: 100)",
    )
    args = parser.parse_args()

    try:
        validate_parameters(sqlite_path=args.sqlite_path, batch_size=args.batch_size)
        stats = fetch_rcsb_entries_metadata(
            sqlite_path=args.sqlite_path,
            batch_size=args.batch_size,
        )
    except (sqlite3.Error, ValueError) as exc:
        parser.error(str(exc))

    print(f"Wrote RCSB entry metadata to {args.sqlite_path}")
    print(f"Entry metadata inserted: {stats.entries_inserted}")
    print(f"Entry metadata fetch failures: {stats.failures}")


if __name__ == "__main__":
    main()
