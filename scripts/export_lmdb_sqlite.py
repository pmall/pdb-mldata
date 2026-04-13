from __future__ import annotations

import argparse
import json
import sqlite3
from dataclasses import dataclass
from pathlib import Path
from typing import cast

import lmdb
from tqdm import tqdm

from pdb_mldata.lmdb_utils import ExportLmdbEntry, unpack_lmdb_entry_for_export

DEFAULT_LMDB_PATH = Path("data/pdb_mldata.lmdb")
DEFAULT_SQLITE_PATH = Path("data/pdb_mldata.sqlite")


@dataclass(frozen=True)
class ExportStats:
    entries_inserted: int
    entities_inserted: int
    pairs_inserted: int


def validate_parameters(
    db_path: Path,
    sqlite_path: Path,
    limit: int | None,
) -> None:
    """Validate CLI parameters before opening LMDB or SQLite."""
    if not db_path.exists():
        raise ValueError(f"{db_path} not found")
    if not db_path.is_dir():
        raise ValueError("db_path must be an LMDB directory path")
    if sqlite_path.exists() and sqlite_path.is_dir():
        raise ValueError("sqlite_path must be a file path, not a directory")
    if limit is not None and limit < 1:
        raise ValueError("--limit must be at least 1 when provided")


def residue_names_json(residue_names: list[str]) -> str:
    return json.dumps(residue_names, separators=(",", ":"))


def assembly_file_stem(pdb_id: str) -> str:
    return f"{pdb_id.lower()}-assembly1"


def create_schema(connection: sqlite3.Connection) -> None:
    """Create the viewer database schema for a fresh export database."""
    connection.executescript(
        """
        PRAGMA foreign_keys = ON;

        CREATE TABLE entries (
            pdb_id TEXT PRIMARY KEY,
            assembly_file_stem TEXT NOT NULL
        );

        CREATE TABLE selected_peptide_entities (
            pdb_id TEXT NOT NULL,
            entity_id TEXT NOT NULL,
            sequence TEXT NOT NULL,
            residue_names_json TEXT NOT NULL,
            PRIMARY KEY (pdb_id, entity_id),
            FOREIGN KEY (pdb_id) REFERENCES entries (pdb_id)
        );

        CREATE TABLE selected_pairs (
            pair_id INTEGER PRIMARY KEY,
            pdb_id TEXT NOT NULL,
            peptide_entity_id TEXT NOT NULL,
            peptide_chain_id TEXT NOT NULL,
            peptide_sequence TEXT NOT NULL,
            peptide_residue_names_json TEXT NOT NULL,
            receptor_entity_id TEXT NOT NULL,
            receptor_chain_id TEXT NOT NULL,
            receptor_sequence TEXT NOT NULL,
            receptor_residue_names_json TEXT NOT NULL,
            FOREIGN KEY (pdb_id) REFERENCES entries (pdb_id),
            FOREIGN KEY (pdb_id, peptide_entity_id)
                REFERENCES selected_peptide_entities (pdb_id, entity_id),
            UNIQUE (
                pdb_id,
                peptide_entity_id,
                peptide_chain_id,
                receptor_entity_id,
                receptor_chain_id
            )
        );

        CREATE INDEX idx_selected_pairs_pdb_id
            ON selected_pairs (pdb_id);
        CREATE INDEX idx_selected_pairs_peptide_chain
            ON selected_pairs (peptide_chain_id);
        CREATE INDEX idx_selected_pairs_receptor_chain
            ON selected_pairs (receptor_chain_id);
        CREATE INDEX idx_selected_pairs_peptide_entity
            ON selected_pairs (pdb_id, peptide_entity_id);
        CREATE INDEX idx_selected_pairs_receptor_entity
            ON selected_pairs (pdb_id, receptor_entity_id);
        """
    )
    connection.commit()


def delete_existing_sqlite_database(sqlite_path: Path) -> None:
    """Remove the previous export artifact before building a fresh read-only database."""
    sqlite_path.unlink(missing_ok=True)
    sqlite_path.with_name(f"{sqlite_path.name}-wal").unlink(missing_ok=True)
    sqlite_path.with_name(f"{sqlite_path.name}-shm").unlink(missing_ok=True)


def fetch_one_row(
    connection: sqlite3.Connection,
    table: str,
    columns: tuple[str, ...],
    key_columns: tuple[str, ...],
    key_values: tuple[object, ...],
) -> tuple[object, ...] | None:
    where_clause = " AND ".join(f"{column} = ?" for column in key_columns)
    sql = f"SELECT {', '.join(columns)} FROM {table} WHERE {where_clause}"
    row = connection.execute(sql, key_values).fetchone()
    if row is None:
        return None
    return tuple(row)


def insert_row(
    connection: sqlite3.Connection,
    table: str,
    columns: tuple[str, ...],
    values: tuple[object, ...],
) -> None:
    placeholders = ", ".join("?" for _ in columns)
    sql = f"INSERT INTO {table} ({', '.join(columns)}) VALUES ({placeholders})"
    connection.execute(sql, values)


def insert_entry(
    connection: sqlite3.Connection,
    entry: ExportLmdbEntry,
) -> None:
    columns = ("pdb_id", "assembly_file_stem")
    values = (entry["pdb_id"], assembly_file_stem(entry["pdb_id"]))
    insert_row(
        connection=connection,
        table="entries",
        columns=columns,
        values=values,
    )


def insert_entity(
    connection: sqlite3.Connection,
    pdb_id: str,
    entity_id: str,
    sequence: str,
    residue_names: list[str],
) -> None:
    columns = ("pdb_id", "entity_id", "sequence", "residue_names_json")
    values = (pdb_id, entity_id, sequence, residue_names_json(residue_names))
    insert_row(
        connection=connection,
        table="selected_peptide_entities",
        columns=columns,
        values=values,
    )


def fetch_pair_by_natural_key(
    connection: sqlite3.Connection,
    columns: tuple[str, ...],
    pdb_id: str,
    peptide_entity_id: str,
    peptide_chain_id: str,
    receptor_entity_id: str,
    receptor_chain_id: str,
) -> tuple[object, ...] | None:
    return fetch_one_row(
        connection=connection,
        table="selected_pairs",
        columns=columns,
        key_columns=(
            "pdb_id",
            "peptide_entity_id",
            "peptide_chain_id",
            "receptor_entity_id",
            "receptor_chain_id",
        ),
        key_values=(
            pdb_id,
            peptide_entity_id,
            peptide_chain_id,
            receptor_entity_id,
            receptor_chain_id,
        ),
    )


def insert_pair(
    connection: sqlite3.Connection,
    pdb_id: str,
    peptide_entity_id: str,
    peptide_chain_id: str,
    peptide_sequence: str,
    peptide_residue_names: list[str],
    receptor_entity_id: str,
    receptor_chain_id: str,
    receptor_sequence: str,
    receptor_residue_names: list[str],
) -> None:
    columns = (
        "pdb_id",
        "peptide_entity_id",
        "peptide_chain_id",
        "peptide_sequence",
        "peptide_residue_names_json",
        "receptor_entity_id",
        "receptor_chain_id",
        "receptor_sequence",
        "receptor_residue_names_json",
    )
    values = (
        pdb_id,
        peptide_entity_id,
        peptide_chain_id,
        peptide_sequence,
        residue_names_json(peptide_residue_names),
        receptor_entity_id,
        receptor_chain_id,
        receptor_sequence,
        residue_names_json(receptor_residue_names),
    )

    existing_by_natural_key = fetch_pair_by_natural_key(
        connection=connection,
        columns=columns,
        pdb_id=pdb_id,
        peptide_entity_id=peptide_entity_id,
        peptide_chain_id=peptide_chain_id,
        receptor_entity_id=receptor_entity_id,
        receptor_chain_id=receptor_chain_id,
    )
    if existing_by_natural_key is not None:
        raise ValueError(
            "duplicate selected pair in LMDB: "
            f"pdb_id={pdb_id}, "
            f"peptide_entity_id={peptide_entity_id}, "
            f"peptide_chain_id={peptide_chain_id}, "
            f"receptor_entity_id={receptor_entity_id}, "
            f"receptor_chain_id={receptor_chain_id}"
        )

    insert_row(
        connection=connection,
        table="selected_pairs",
        columns=columns,
        values=values,
    )


def export_entry(
    connection: sqlite3.Connection,
    entry: ExportLmdbEntry,
) -> tuple[int, int, int]:
    insert_entry(connection, entry)
    entries_inserted = 1
    entities_inserted = 0
    pairs_inserted = 0

    for entity in entry["entities"]:
        insert_entity(
            connection=connection,
            pdb_id=entry["pdb_id"],
            entity_id=entity["entity_id"],
            sequence=entity["sequence"],
            residue_names=entity["residue_names"],
        )
        entities_inserted += 1

        for pair in entity["pairs"]:
            insert_pair(
                connection=connection,
                pdb_id=entry["pdb_id"],
                peptide_entity_id=pair["peptide"]["entity_id"],
                peptide_chain_id=pair["peptide"]["chain"],
                peptide_sequence=pair["peptide"]["sequence"],
                peptide_residue_names=pair["peptide"]["residue_names"],
                receptor_entity_id=pair["receptor"]["entity_id"],
                receptor_chain_id=pair["receptor"]["chain"],
                receptor_sequence=pair["receptor"]["sequence"],
                receptor_residue_names=pair["receptor"]["residue_names"],
            )
            pairs_inserted += 1

    return (
        entries_inserted,
        entities_inserted,
        pairs_inserted,
    )


def export_lmdb_to_sqlite(
    db_path: Path,
    sqlite_path: Path,
    limit: int | None,
) -> ExportStats:
    sqlite_path.parent.mkdir(parents=True, exist_ok=True)
    delete_existing_sqlite_database(sqlite_path)
    connection = sqlite3.connect(sqlite_path)
    create_schema(connection)

    entries_inserted = 0
    entities_inserted = 0
    pairs_inserted = 0

    env = lmdb.open(str(db_path), readonly=True, lock=False)
    with env.begin() as txn:
        total_entries = txn.stat()["entries"]
        if limit is not None:
            total_entries = min(total_entries, limit)

        cursor = txn.cursor()
        progress = tqdm(cursor, total=total_entries, desc="Exporting selections")
        for entry_index, (_key, value) in enumerate(progress):
            if limit is not None and entry_index >= limit:
                break

            entry = unpack_lmdb_entry_for_export(cast(bytes, value))
            (
                entry_inserted_count,
                entity_inserted_count,
                pair_inserted_count,
            ) = export_entry(
                connection=connection,
                entry=entry,
            )
            entries_inserted += entry_inserted_count
            entities_inserted += entity_inserted_count
            pairs_inserted += pair_inserted_count
            progress.set_postfix(
                entries=entries_inserted,
                entities=entities_inserted,
                pairs=pairs_inserted,
            )

    env.close()
    connection.commit()
    connection.close()

    return ExportStats(
        entries_inserted=entries_inserted,
        entities_inserted=entities_inserted,
        pairs_inserted=pairs_inserted,
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Export parsed peptide/receptor selections from LMDB to SQLite."
    )
    parser.add_argument(
        "db_path",
        nargs="?",
        type=Path,
        default=DEFAULT_LMDB_PATH,
        help="Path to the LMDB database directory (default: data/pdb_mldata.lmdb)",
    )
    parser.add_argument(
        "sqlite_path",
        nargs="?",
        type=Path,
        default=DEFAULT_SQLITE_PATH,
        help="Path to the output SQLite database (default: data/pdb_mldata.sqlite)",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Export only the first N LMDB entries for testing",
    )
    args = parser.parse_args()

    try:
        validate_parameters(
            db_path=args.db_path,
            sqlite_path=args.sqlite_path,
            limit=args.limit,
        )
        stats = export_lmdb_to_sqlite(
            db_path=args.db_path,
            sqlite_path=args.sqlite_path,
            limit=args.limit,
        )
    except (sqlite3.Error, ValueError) as exc:
        parser.error(str(exc))

    print(f"Wrote selection export to {args.sqlite_path}")
    print(f"Entries inserted: {stats.entries_inserted}")
    print(f"Peptide entities inserted: {stats.entities_inserted}")
    print(f"Pairs inserted: {stats.pairs_inserted}")


if __name__ == "__main__":
    main()
