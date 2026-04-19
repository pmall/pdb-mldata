from __future__ import annotations

import argparse
import sqlite3
from dataclasses import dataclass
from pathlib import Path

from tqdm import tqdm

from pdb_mldata.viewer_rcsb import (
    DEFAULT_BATCH_SIZE,
    POLYMER_ENTITY_METADATA_QUERY,
    batched,
    configure_session,
    entity_graphql_id,
    extract_polymer_entity_metadata,
    fetch_graphql,
    first_string,
    parse_entity_graphql_id,
)

DEFAULT_SQLITE_PATH = Path("data/pdb_mldata.sqlite")


@dataclass(frozen=True)
class EntityMetadataStats:
    peptides_inserted: int
    targets_inserted: int
    failures: int


def validate_parameters(sqlite_path: Path, batch_size: int) -> None:
    """Validate CLI parameters before opening the viewer database."""
    if not sqlite_path.exists():
        raise ValueError(f"{sqlite_path} not found")
    if sqlite_path.is_dir():
        raise ValueError("sqlite_path must be a file path, not a directory")
    if batch_size < 1:
        raise ValueError("--batch-size must be at least 1")


def rebuild_entities_metadata_schema(connection: sqlite3.Connection) -> None:
    """Delete entity-derived metadata before rebuilding RCSB entity metadata."""
    connection.executescript(
        """
        PRAGMA foreign_keys = ON;

        DROP TABLE IF EXISTS search_terms_targets;
        DROP TABLE IF EXISTS search_terms;
        DROP TABLE IF EXISTS peptides_accessions;
        DROP TABLE IF EXISTS targets_accessions;
        DROP TABLE IF EXISTS uniprot_targets;
        DROP TABLE IF EXISTS uniprot_metadata;
        DROP TABLE IF EXISTS peptides_metadata;
        DROP TABLE IF EXISTS targets_metadata;

        CREATE TABLE peptides_metadata (
            pdb_id TEXT NOT NULL,
            entity_id TEXT NOT NULL,
            entity_name TEXT NOT NULL,
            organism_scientific_names_json TEXT NOT NULL,
            taxonomy_ids_json TEXT NOT NULL,
            accessions_json TEXT NOT NULL,
            polymer_type TEXT NOT NULL,
            sequence_length INTEGER,
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
            sequence_length INTEGER,
            PRIMARY KEY (pdb_id, entity_id),
            FOREIGN KEY (pdb_id) REFERENCES entries (pdb_id)
        );

        CREATE INDEX idx_peptides_metadata__entity_name
            ON peptides_metadata (entity_name);
        CREATE INDEX idx_targets_metadata__entity_name
            ON targets_metadata (entity_name);
        """
    )
    connection.commit()


def read_entity_ids(
    connection: sqlite3.Connection,
) -> tuple[list[tuple[str, str]], list[tuple[str, str]]]:
    peptide_ids = [
        (row[0], row[1])
        for row in connection.execute(
            """
            SELECT DISTINCT pdb_id, peptide_entity_id AS entity_id
            FROM chain_pairs
            ORDER BY pdb_id, entity_id
            """
        ).fetchall()
    ]
    target_ids = [
        (row[0], row[1])
        for row in connection.execute(
            """
            SELECT DISTINCT pdb_id, receptor_entity_id AS entity_id
            FROM chain_pairs
            ORDER BY pdb_id, entity_id
            """
        ).fetchall()
    ]
    return peptide_ids, target_ids


def insert_row(
    connection: sqlite3.Connection,
    table: str,
    columns: tuple[str, ...],
    values: tuple[object, ...],
) -> None:
    placeholders = ", ".join("?" for _ in columns)
    sql = f"INSERT INTO {table} ({', '.join(columns)}) VALUES ({placeholders})"
    connection.execute(sql, values)


def fetch_entity_metadata(
    connection: sqlite3.Connection,
    peptide_ids: list[tuple[str, str]],
    target_ids: list[tuple[str, str]],
    batch_size: int,
) -> EntityMetadataStats:
    columns = (
        "pdb_id",
        "entity_id",
        "entity_name",
        "organism_scientific_names_json",
        "taxonomy_ids_json",
        "accessions_json",
        "polymer_type",
        "sequence_length",
    )
    session = configure_session()
    peptides_inserted = 0
    targets_inserted = 0
    failures = 0
    peptide_graphql_ids = {
        entity_graphql_id(pdb_id=pdb_id, entity_id=entity_id)
        for pdb_id, entity_id in peptide_ids
    }
    target_graphql_ids = {
        entity_graphql_id(pdb_id=pdb_id, entity_id=entity_id)
        for pdb_id, entity_id in target_ids
    }
    graphql_ids = sorted(peptide_graphql_ids | target_graphql_ids)

    with tqdm(
        total=len(graphql_ids),
        desc="Fetching RCSB entity metadata",
        unit="entity",
    ) as progress:
        for entity_id_batch in batched(graphql_ids, batch_size):
            try:
                payload = fetch_graphql(
                    session=session,
                    query=POLYMER_ENTITY_METADATA_QUERY,
                    variables={"entityIds": entity_id_batch},
                )
            except Exception:
                failures += len(entity_id_batch)
                progress.update(len(entity_id_batch))
                continue

            rows = payload.get("polymer_entities")
            if not isinstance(rows, list):
                failures += len(entity_id_batch)
                progress.update(len(entity_id_batch))
                continue

            returned_ids: set[str] = set()
            for row in rows:
                if not isinstance(row, dict):
                    failures += 1
                    continue
                rcsb_id = first_string(row, "rcsb_id")
                if not rcsb_id:
                    failures += 1
                    continue
                try:
                    pdb_id, entity_id = parse_entity_graphql_id(rcsb_id)
                    values = extract_polymer_entity_metadata(
                        pdb_id=pdb_id,
                        entity_id=entity_id,
                        data=row,
                    )
                    if rcsb_id in peptide_graphql_ids:
                        insert_row(
                            connection=connection,
                            table="peptides_metadata",
                            columns=columns,
                            values=values,
                        )
                        peptides_inserted += 1
                    if rcsb_id in target_graphql_ids:
                        insert_row(
                            connection=connection,
                            table="targets_metadata",
                            columns=columns,
                            values=values,
                        )
                        targets_inserted += 1
                    returned_ids.add(rcsb_id)
                except Exception:
                    failures += 1

            failures += len(set(entity_id_batch) - returned_ids)
            progress.update(len(entity_id_batch))

    return EntityMetadataStats(
        peptides_inserted=peptides_inserted,
        targets_inserted=targets_inserted,
        failures=failures,
    )


def fetch_rcsb_entities_metadata(
    sqlite_path: Path, batch_size: int
) -> EntityMetadataStats:
    connection = sqlite3.connect(sqlite_path)
    rebuild_entities_metadata_schema(connection)
    peptide_ids, target_ids = read_entity_ids(connection)
    stats = fetch_entity_metadata(
        connection=connection,
        peptide_ids=peptide_ids,
        target_ids=target_ids,
        batch_size=batch_size,
    )
    connection.commit()
    connection.close()
    return stats


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Fetch RCSB peptide and target entity metadata for viewer rows."
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
        stats = fetch_rcsb_entities_metadata(
            sqlite_path=args.sqlite_path,
            batch_size=args.batch_size,
        )
    except (sqlite3.Error, ValueError) as exc:
        parser.error(str(exc))

    print(f"Wrote RCSB entity metadata to {args.sqlite_path}")
    print(f"Peptide metadata inserted: {stats.peptides_inserted}")
    print(f"Target metadata inserted: {stats.targets_inserted}")
    print(f"Entity metadata fetch failures: {stats.failures}")


if __name__ == "__main__":
    main()
