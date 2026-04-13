from __future__ import annotations

import argparse
import sqlite3
from pathlib import Path
from typing import TypeAlias

from tqdm import tqdm

DEFAULT_SQLITE_PATH = Path("data/pdb_mldata.sqlite")
DEFAULT_SQL_PATH = Path("data/pdb_mldata.sql")

SqlValue: TypeAlias = int | float | str | None

TABLE_COLUMNS = {
    "entries": ("pdb_id", "assembly_file_stem"),
    "selected_peptide_entities": (
        "pdb_id",
        "entity_id",
        "sequence",
        "residue_names_json",
    ),
    "selected_pairs": (
        "pair_id",
        "pdb_id",
        "peptide_entity_id",
        "peptide_chain_id",
        "peptide_sequence",
        "peptide_residue_names_json",
        "receptor_entity_id",
        "receptor_chain_id",
        "receptor_sequence",
        "receptor_residue_names_json",
    ),
    "entry_metadata": (
        "pdb_id",
        "title",
        "experimental_methods_json",
        "resolution_combined_json",
        "best_resolution",
        "deposition_date",
        "initial_release_date",
    ),
    "polymer_entity_metadata": (
        "pdb_id",
        "entity_id",
        "entity_name",
        "organism_scientific_names_json",
        "taxonomy_ids_json",
        "uniprot_accessions_json",
        "polymer_type",
        "sequence_length",
    ),
    "metadata_fetch_failures": (
        "failure_id",
        "fetched_at",
        "metadata_scope",
        "pdb_id",
        "entity_id",
        "error_message",
    ),
}
TABLE_ORDER = tuple(TABLE_COLUMNS)
POSTGRES_SCHEMA = """
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
    pair_id BIGINT PRIMARY KEY,
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

CREATE TABLE entry_metadata (
    pdb_id TEXT PRIMARY KEY,
    title TEXT NOT NULL,
    experimental_methods_json TEXT NOT NULL,
    resolution_combined_json TEXT NOT NULL,
    best_resolution DOUBLE PRECISION,
    deposition_date TEXT,
    initial_release_date TEXT,
    FOREIGN KEY (pdb_id) REFERENCES entries (pdb_id)
);

CREATE TABLE polymer_entity_metadata (
    pdb_id TEXT NOT NULL,
    entity_id TEXT NOT NULL,
    entity_name TEXT NOT NULL,
    organism_scientific_names_json TEXT NOT NULL,
    taxonomy_ids_json TEXT NOT NULL,
    uniprot_accessions_json TEXT NOT NULL,
    polymer_type TEXT NOT NULL,
    sequence_length BIGINT,
    PRIMARY KEY (pdb_id, entity_id),
    FOREIGN KEY (pdb_id) REFERENCES entries (pdb_id)
);

CREATE TABLE metadata_fetch_failures (
    failure_id BIGINT PRIMARY KEY,
    fetched_at TEXT NOT NULL,
    metadata_scope TEXT NOT NULL,
    pdb_id TEXT NOT NULL,
    entity_id TEXT,
    error_message TEXT NOT NULL,
    FOREIGN KEY (pdb_id) REFERENCES entries (pdb_id)
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
CREATE INDEX idx_polymer_entity_metadata_entity_name
    ON polymer_entity_metadata (entity_name);
CREATE INDEX idx_metadata_failures_pdb_id
    ON metadata_fetch_failures (pdb_id);
CREATE INDEX idx_metadata_failures_scope
    ON metadata_fetch_failures (metadata_scope);
"""


def validate_parameters(sqlite_path: Path, sql_path: Path) -> None:
    """Validate CLI parameters before opening the SQLite database."""
    if not sqlite_path.exists():
        raise ValueError(f"{sqlite_path} not found")
    if sqlite_path.is_dir():
        raise ValueError("sqlite_path must be a file path, not a directory")
    if sql_path.exists() and sql_path.is_dir():
        raise ValueError("sql_path must be a file path, not a directory")


def quote_identifier(identifier: str) -> str:
    return '"' + identifier.replace('"', '""') + '"'


def quote_sql_value(value: SqlValue) -> str:
    if value is None:
        return "NULL"
    if isinstance(value, int | float):
        return str(value)
    return "'" + value.replace("'", "''") + "'"


def expected_columns_sql(table: str) -> str:
    return ", ".join(quote_identifier(column) for column in TABLE_COLUMNS[table])


def validate_sqlite_schema(connection: sqlite3.Connection) -> None:
    """Fail if the SQLite database is not the viewer database we know how to dump."""
    for table, expected_columns in TABLE_COLUMNS.items():
        rows = connection.execute(f"PRAGMA table_info({quote_identifier(table)})")
        actual_columns = tuple(row[1] for row in rows.fetchall())
        if actual_columns != expected_columns:
            raise ValueError(
                f"{table} columns {actual_columns} do not match {expected_columns}"
            )


def count_dumped_rows(connection: sqlite3.Connection) -> int:
    total = 0
    for table in TABLE_ORDER:
        total += connection.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0]
    return total


def iter_table_rows(
    connection: sqlite3.Connection,
    table: str,
) -> sqlite3.Cursor:
    columns_sql = expected_columns_sql(table)
    return connection.execute(
        f"SELECT {columns_sql} FROM {quote_identifier(table)} ORDER BY rowid"
    )


def render_insert_statement(table: str, row: tuple[SqlValue, ...]) -> str:
    columns_sql = expected_columns_sql(table)
    values_sql = ", ".join(quote_sql_value(value) for value in row)
    return (
        f"INSERT INTO {quote_identifier(table)} ({columns_sql}) VALUES ({values_sql});"
    )


def dump_sqlite_to_postgres_sql(sqlite_path: Path, sql_path: Path) -> int:
    sql_path.parent.mkdir(parents=True, exist_ok=True)
    sql_path.unlink(missing_ok=True)

    connection = sqlite3.connect(sqlite_path)
    try:
        validate_sqlite_schema(connection)
        total_rows = count_dumped_rows(connection)

        statement_count = 0
        with sql_path.open("w", encoding="utf-8") as handle:
            handle.write("BEGIN;\n")
            handle.write(POSTGRES_SCHEMA.strip())
            handle.write("\n")
            statement_count += POSTGRES_SCHEMA.count("CREATE ")

            with tqdm(total=total_rows, desc="Dumping SQL", unit="row") as progress:
                for table in TABLE_ORDER:
                    for row in iter_table_rows(connection, table):
                        handle.write(render_insert_statement(table, row))
                        handle.write("\n")
                        statement_count += 1
                        progress.update(1)

            handle.write("COMMIT;\n")
            statement_count += 2
    finally:
        connection.close()

    return statement_count


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Dump the SQLite viewer database to PostgreSQL-compatible SQL."
    )
    parser.add_argument(
        "sqlite_path",
        nargs="?",
        type=Path,
        default=DEFAULT_SQLITE_PATH,
        help="Path to the SQLite database (default: data/pdb_mldata.sqlite)",
    )
    parser.add_argument(
        "sql_path",
        nargs="?",
        type=Path,
        default=DEFAULT_SQL_PATH,
        help="Path to the SQL dump file (default: data/pdb_mldata.sql)",
    )
    args = parser.parse_args()

    try:
        validate_parameters(sqlite_path=args.sqlite_path, sql_path=args.sql_path)
        statement_count = dump_sqlite_to_postgres_sql(
            sqlite_path=args.sqlite_path,
            sql_path=args.sql_path,
        )
    except (OSError, sqlite3.Error, ValueError) as exc:
        parser.error(str(exc))

    print(f"Wrote PostgreSQL SQL dump to {args.sql_path}")
    print(f"SQL statements: {statement_count}")


if __name__ == "__main__":
    main()
