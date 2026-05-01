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
    "peptides": (
        "pdb_id",
        "entity_id",
        "sequence",
        "residue_names_json",
    ),
    "chain_pairs": (
        "chain_pair_id",
        "pdb_id",
        "peptide_entity_id",
        "peptide_chain_id",
        "peptide_sequence",
        "peptide_residue_names_json",
        "peptide_binding_start",
        "peptide_binding_stop",
        "receptor_entity_id",
        "receptor_chain_id",
        "receptor_sequence",
        "receptor_residue_names_json",
        "receptor_binding_start",
        "receptor_binding_stop",
    ),
    "entries_metadata": (
        "pdb_id",
        "title",
        "experimental_methods_json",
        "resolution_combined_json",
        "best_resolution",
        "deposition_date",
        "initial_release_date",
    ),
    "peptides_metadata": (
        "pdb_id",
        "entity_id",
        "entity_name",
        "organism_scientific_names_json",
        "taxonomy_ids_json",
        "accessions_json",
        "polymer_type",
        "sequence_length",
    ),
    "targets_metadata": (
        "pdb_id",
        "entity_id",
        "entity_name",
        "organism_scientific_names_json",
        "taxonomy_ids_json",
        "accessions_json",
        "polymer_type",
        "sequence_length",
    ),
    "uniprot_metadata": (
        "accession",
        "reviewed",
        "recommended_name",
        "gene_names_json",
        "organism_scientific_name",
        "taxonomy_id",
        "function_text",
        "subcellular_locations_json",
        "keywords_json",
        "go_terms_json",
        "interpro_ids_json",
        "pfam_ids_json",
    ),
    "peptides_accessions": (
        "pdb_id",
        "entity_id",
        "accession",
    ),
    "targets_accessions": (
        "pdb_id",
        "entity_id",
        "accession",
    ),
    "search_terms": (
        "search_term_id",
        "term",
        "term_kind",
        "rank_weight",
    ),
    "search_terms_targets": (
        "search_terms_target_id",
        "search_term_id",
        "pdb_id",
        "entity_id",
    ),
}
TABLE_ORDER = tuple(TABLE_COLUMNS)
POSTGRES_SCHEMA = """
CREATE EXTENSION IF NOT EXISTS pg_trgm;

CREATE TABLE entries (
    pdb_id TEXT PRIMARY KEY,
    assembly_file_stem TEXT NOT NULL
);

CREATE TABLE peptides (
    pdb_id TEXT NOT NULL,
    entity_id TEXT NOT NULL,
    sequence TEXT NOT NULL,
    residue_names_json TEXT NOT NULL,
    PRIMARY KEY (pdb_id, entity_id),
    FOREIGN KEY (pdb_id) REFERENCES entries (pdb_id)
);

CREATE TABLE entries_metadata (
    pdb_id TEXT PRIMARY KEY,
    title TEXT NOT NULL,
    experimental_methods_json TEXT NOT NULL,
    resolution_combined_json TEXT NOT NULL,
    best_resolution DOUBLE PRECISION,
    deposition_date TEXT,
    initial_release_date TEXT,
    FOREIGN KEY (pdb_id) REFERENCES entries (pdb_id)
);

CREATE TABLE chain_pairs (
    chain_pair_id BIGINT PRIMARY KEY,
    pdb_id TEXT NOT NULL,
    peptide_entity_id TEXT NOT NULL,
    peptide_chain_id TEXT NOT NULL,
    peptide_sequence TEXT NOT NULL,
    peptide_residue_names_json TEXT NOT NULL,
    peptide_binding_start BIGINT NOT NULL,
    peptide_binding_stop BIGINT NOT NULL,
    receptor_entity_id TEXT NOT NULL,
    receptor_chain_id TEXT NOT NULL,
    receptor_sequence TEXT NOT NULL,
    receptor_residue_names_json TEXT NOT NULL,
    receptor_binding_start BIGINT NOT NULL,
    receptor_binding_stop BIGINT NOT NULL,
    FOREIGN KEY (pdb_id) REFERENCES entries (pdb_id),
    FOREIGN KEY (pdb_id, peptide_entity_id)
        REFERENCES peptides (pdb_id, entity_id),
    UNIQUE (
        pdb_id,
        peptide_entity_id,
        peptide_chain_id,
        receptor_entity_id,
        receptor_chain_id
    )
);

CREATE TABLE peptides_metadata (
    pdb_id TEXT NOT NULL,
    entity_id TEXT NOT NULL,
    entity_name TEXT NOT NULL,
    organism_scientific_names_json TEXT NOT NULL,
    taxonomy_ids_json TEXT NOT NULL,
    accessions_json TEXT NOT NULL,
    polymer_type TEXT NOT NULL,
    sequence_length BIGINT,
    PRIMARY KEY (pdb_id, entity_id),
    FOREIGN KEY (pdb_id, entity_id)
        REFERENCES peptides (pdb_id, entity_id)
);

CREATE TABLE targets_metadata (
    pdb_id TEXT NOT NULL,
    entity_id TEXT NOT NULL,
    entity_name TEXT NOT NULL,
    organism_scientific_names_json TEXT NOT NULL,
    taxonomy_ids_json TEXT NOT NULL,
    accessions_json TEXT NOT NULL,
    polymer_type TEXT NOT NULL,
    sequence_length BIGINT,
    PRIMARY KEY (pdb_id, entity_id),
    FOREIGN KEY (pdb_id) REFERENCES entries (pdb_id)
);

CREATE TABLE uniprot_metadata (
    accession TEXT PRIMARY KEY,
    reviewed BIGINT NOT NULL,
    recommended_name TEXT NOT NULL,
    gene_names_json TEXT NOT NULL,
    organism_scientific_name TEXT NOT NULL,
    taxonomy_id BIGINT,
    function_text TEXT NOT NULL,
    subcellular_locations_json TEXT NOT NULL,
    keywords_json TEXT NOT NULL,
    go_terms_json TEXT NOT NULL,
    interpro_ids_json TEXT NOT NULL,
    pfam_ids_json TEXT NOT NULL
);

CREATE TABLE peptides_accessions (
    pdb_id TEXT NOT NULL,
    entity_id TEXT NOT NULL,
    accession TEXT NOT NULL,
    FOREIGN KEY (pdb_id, entity_id)
        REFERENCES peptides (pdb_id, entity_id),
    FOREIGN KEY (pdb_id, entity_id)
        REFERENCES peptides_metadata (pdb_id, entity_id),
    FOREIGN KEY (accession) REFERENCES uniprot_metadata (accession),
    PRIMARY KEY (pdb_id, entity_id, accession)
);

CREATE TABLE targets_accessions (
    pdb_id TEXT NOT NULL,
    entity_id TEXT NOT NULL,
    accession TEXT NOT NULL,
    FOREIGN KEY (pdb_id, entity_id)
        REFERENCES targets_metadata (pdb_id, entity_id),
    FOREIGN KEY (accession) REFERENCES uniprot_metadata (accession),
    PRIMARY KEY (pdb_id, entity_id, accession)
);

CREATE TABLE search_terms (
    search_term_id BIGINT PRIMARY KEY,
    term TEXT NOT NULL,
    term_kind TEXT NOT NULL,
    rank_weight BIGINT NOT NULL,
    UNIQUE (term)
);

CREATE TABLE search_terms_targets (
    search_terms_target_id BIGINT PRIMARY KEY,
    search_term_id BIGINT NOT NULL,
    pdb_id TEXT NOT NULL,
    entity_id TEXT NOT NULL,
    FOREIGN KEY (search_term_id) REFERENCES search_terms (search_term_id),
    FOREIGN KEY (pdb_id) REFERENCES entries (pdb_id),
    FOREIGN KEY (pdb_id, entity_id)
        REFERENCES targets_metadata (pdb_id, entity_id),
    UNIQUE (search_term_id, pdb_id, entity_id)
);
"""

POSTGRES_INDEXES = """
CREATE INDEX idx_targets_metadata__entity_name
    ON targets_metadata (entity_name);
CREATE INDEX idx_peptides_metadata__entity_name
    ON peptides_metadata (entity_name);
CREATE INDEX idx_chain_pairs__peptide_chain_id
    ON chain_pairs (peptide_chain_id);
CREATE INDEX idx_chain_pairs__receptor_chain_id
    ON chain_pairs (receptor_chain_id);
CREATE INDEX idx_chain_pairs__pdb_id_peptide_entity_id
    ON chain_pairs (pdb_id, peptide_entity_id);
CREATE INDEX idx_chain_pairs__pdb_id_receptor_entity_id
    ON chain_pairs (pdb_id, receptor_entity_id);
CREATE INDEX idx_peptides_accessions__accession
    ON peptides_accessions (accession);
CREATE INDEX idx_targets_accessions__accession
    ON targets_accessions (accession);
CREATE INDEX idx_search_terms__term_kind
    ON search_terms (term_kind);
CREATE INDEX idx_search_terms_targets__search_term_id
    ON search_terms_targets (search_term_id);
CREATE INDEX idx_search_terms_targets__pdb_id_entity_id
    ON search_terms_targets (pdb_id, entity_id);
CREATE INDEX idx_search_terms__term_trgm
    ON search_terms
    USING GIN (term gin_trgm_ops);
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


def quote_copy_value(value: SqlValue) -> str:
    if value is None:
        return r"\N"
    if isinstance(value, int | float):
        return str(value)
    return (
        value.replace("\\", "\\\\")
        .replace("\t", r"\t")
        .replace("\n", r"\n")
        .replace("\r", r"\r")
    )


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


def render_copy_row(row: tuple[SqlValue, ...]) -> str:
    return "\t".join(quote_copy_value(value) for value in row)


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
                    columns_sql = expected_columns_sql(table)
                    handle.write(
                        f"COPY {quote_identifier(table)} ({columns_sql}) FROM stdin;\n"
                    )
                    statement_count += 1
                    for row in iter_table_rows(connection, table):
                        handle.write(render_copy_row(row))
                        handle.write("\n")
                        progress.update(1)
                    handle.write("\\.\n")

            handle.write(POSTGRES_INDEXES.strip())
            handle.write("\n")
            statement_count += POSTGRES_INDEXES.count("CREATE ")
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
