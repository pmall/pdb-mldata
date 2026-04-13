from __future__ import annotations

import argparse
import json
import sqlite3
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import TypeAlias

import requests
from requests.adapters import HTTPAdapter
from tqdm import tqdm
from urllib3.util.retry import Retry

DEFAULT_SQLITE_PATH = Path("data/pdb_mldata.sqlite")
DEFAULT_TIMEOUT = 20.0

JsonObject: TypeAlias = dict[str, object]


@dataclass(frozen=True)
class MetadataStats:
    entries_inserted: int
    entities_inserted: int
    failures_inserted: int


def validate_parameters(sqlite_path: Path) -> None:
    """Validate CLI parameters before opening the viewer database."""
    if not sqlite_path.exists():
        raise ValueError(f"{sqlite_path} not found")
    if sqlite_path.is_dir():
        raise ValueError("sqlite_path must be a file path, not a directory")


def json_text(value: object) -> str:
    return json.dumps(value, separators=(",", ":"), sort_keys=True)


def current_utc_timestamp() -> str:
    return datetime.now(UTC).isoformat()


def rebuild_metadata_schema(connection: sqlite3.Connection) -> None:
    """Delete only metadata tables before rebuilding metadata from selected rows."""
    connection.executescript(
        """
        PRAGMA foreign_keys = ON;

        DROP TABLE IF EXISTS metadata_fetch_failures;
        DROP TABLE IF EXISTS polymer_entity_metadata;
        DROP TABLE IF EXISTS entry_metadata;

        CREATE TABLE entry_metadata (
            pdb_id TEXT PRIMARY KEY,
            title TEXT NOT NULL,
            experimental_methods_json TEXT NOT NULL,
            resolution_combined_json TEXT NOT NULL,
            best_resolution REAL,
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
            sequence_length INTEGER,
            PRIMARY KEY (pdb_id, entity_id),
            FOREIGN KEY (pdb_id) REFERENCES entries (pdb_id)
        );

        CREATE TABLE metadata_fetch_failures (
            failure_id INTEGER PRIMARY KEY,
            fetched_at TEXT NOT NULL,
            metadata_scope TEXT NOT NULL,
            pdb_id TEXT NOT NULL,
            entity_id TEXT,
            error_message TEXT NOT NULL,
            FOREIGN KEY (pdb_id) REFERENCES entries (pdb_id)
        );

        CREATE INDEX idx_polymer_entity_metadata_entity_name
            ON polymer_entity_metadata (entity_name);
        CREATE INDEX idx_metadata_failures_pdb_id
            ON metadata_fetch_failures (pdb_id);
        CREATE INDEX idx_metadata_failures_scope
            ON metadata_fetch_failures (metadata_scope);
        """
    )
    connection.commit()


def configure_session() -> requests.Session:
    session = requests.Session()
    retries = Retry(total=3, backoff_factor=0.5, status_forcelist=[500, 502, 503, 504])
    session.mount("https://", HTTPAdapter(max_retries=retries))
    return session


def fetch_json(session: requests.Session, url: str) -> JsonObject:
    response = session.get(url, timeout=DEFAULT_TIMEOUT)
    response.raise_for_status()
    data = response.json()
    if not isinstance(data, dict):
        raise ValueError("RCSB response was not a JSON object")
    return data


def read_entry_ids(connection: sqlite3.Connection) -> list[str]:
    return [
        row[0]
        for row in connection.execute(
            "SELECT pdb_id FROM entries ORDER BY pdb_id"
        ).fetchall()
    ]


def read_polymer_entity_ids(
    connection: sqlite3.Connection,
) -> list[tuple[str, str]]:
    return [
        (row[0], row[1])
        for row in connection.execute(
            """
            SELECT DISTINCT pdb_id, entity_id
            FROM (
                SELECT pdb_id, entity_id
                FROM selected_peptide_entities
                UNION
                SELECT pdb_id, receptor_entity_id AS entity_id
                FROM selected_pairs
            )
            ORDER BY pdb_id, entity_id
            """
        ).fetchall()
    ]


def sorted_strings(values: object) -> list[str]:
    if not isinstance(values, list):
        return []
    strings = [value for value in values if isinstance(value, str)]
    return sorted(set(strings))


def sorted_ints(values: list[object]) -> list[int]:
    ints = [value for value in values if isinstance(value, int)]
    return sorted(set(ints))


def first_string(data: JsonObject, key: str) -> str:
    value = data.get(key)
    return value if isinstance(value, str) else ""


def first_int(data: JsonObject, key: str) -> int | None:
    value = data.get(key)
    return value if isinstance(value, int) else None


def best_resolution(resolutions: list[object]) -> float | None:
    numeric_resolutions = [
        value for value in resolutions if isinstance(value, int | float)
    ]
    if not numeric_resolutions:
        return None
    return float(min(numeric_resolutions))


def extract_entry_metadata(pdb_id: str, data: JsonObject) -> tuple[object, ...]:
    struct = data.get("struct")
    entry_info = data.get("rcsb_entry_info")
    accession_info = data.get("rcsb_accession_info")
    database_status = data.get("pdbx_database_status")
    exptl = data.get("exptl")

    title = first_string(struct, "title") if isinstance(struct, dict) else ""
    if not title:
        title = ""

    methods: list[str] = []
    if isinstance(exptl, list):
        for item in exptl:
            if isinstance(item, dict):
                method = item.get("method")
                if isinstance(method, str):
                    methods.append(method)
    methods = sorted(set(methods))

    resolutions: list[object] = []
    if isinstance(entry_info, dict):
        raw_resolutions = entry_info.get("resolution_combined")
        if isinstance(raw_resolutions, list):
            resolutions = raw_resolutions

    deposition_date = (
        first_string(database_status, "recvd_initial_deposition_date")
        if isinstance(database_status, dict)
        else None
    )
    initial_release_date = (
        first_string(accession_info, "initial_release_date")
        if isinstance(accession_info, dict)
        else None
    )

    return (
        pdb_id,
        title,
        json_text(methods),
        json_text(resolutions),
        best_resolution(resolutions),
        deposition_date or None,
        initial_release_date or None,
    )


def extract_uniprot_accessions(identifiers: JsonObject) -> list[str]:
    accessions = set(sorted_strings(identifiers.get("uniprot_ids")))
    reference_identifiers = identifiers.get("reference_sequence_identifiers")
    if isinstance(reference_identifiers, list):
        for item in reference_identifiers:
            if not isinstance(item, dict):
                continue
            database_name = item.get("database_name")
            accession = item.get("database_accession")
            if database_name == "UniProt" and isinstance(accession, str):
                accessions.add(accession)
    return sorted(accessions)


def extract_polymer_entity_metadata(
    pdb_id: str,
    entity_id: str,
    data: JsonObject,
) -> tuple[object, ...]:
    polymer_entity = data.get("rcsb_polymer_entity")
    entity_poly = data.get("entity_poly")
    source_organisms = data.get("rcsb_entity_source_organism")
    identifiers = data.get("rcsb_polymer_entity_container_identifiers")

    entity_name = (
        first_string(polymer_entity, "pdbx_description")
        if isinstance(polymer_entity, dict)
        else ""
    )
    polymer_type = (
        first_string(entity_poly, "rcsb_entity_polymer_type")
        if isinstance(entity_poly, dict)
        else ""
    )
    sequence_length = (
        first_int(entity_poly, "rcsb_sample_sequence_length")
        if isinstance(entity_poly, dict)
        else None
    )

    organism_names: list[str] = []
    taxonomy_values: list[object] = []
    if isinstance(source_organisms, list):
        for organism in source_organisms:
            if not isinstance(organism, dict):
                continue
            scientific_name = organism.get("scientific_name")
            taxonomy_id = organism.get("ncbi_taxonomy_id")
            if isinstance(scientific_name, str):
                organism_names.append(scientific_name)
            if taxonomy_id is not None:
                taxonomy_values.append(taxonomy_id)

    uniprot_accessions = (
        extract_uniprot_accessions(identifiers) if isinstance(identifiers, dict) else []
    )

    return (
        pdb_id,
        entity_id,
        entity_name,
        json_text(sorted(set(organism_names))),
        json_text(sorted_ints(taxonomy_values)),
        json_text(uniprot_accessions),
        polymer_type,
        sequence_length,
    )


def insert_row(
    connection: sqlite3.Connection,
    table: str,
    columns: tuple[str, ...],
    values: tuple[object, ...],
) -> None:
    placeholders = ", ".join("?" for _ in columns)
    sql = f"INSERT INTO {table} ({', '.join(columns)}) VALUES ({placeholders})"
    connection.execute(sql, values)


def insert_metadata_failure(
    connection: sqlite3.Connection,
    metadata_scope: str,
    pdb_id: str,
    entity_id: str | None,
    error_message: str,
) -> None:
    connection.execute(
        """
        INSERT INTO metadata_fetch_failures (
            fetched_at,
            metadata_scope,
            pdb_id,
            entity_id,
            error_message
        )
        VALUES (?, ?, ?, ?, ?)
        """,
        (
            current_utc_timestamp(),
            metadata_scope,
            pdb_id,
            entity_id,
            error_message,
        ),
    )


def fetch_entry_metadata(
    connection: sqlite3.Connection,
    session: requests.Session,
    pdb_ids: list[str],
) -> tuple[int, int]:
    columns = (
        "pdb_id",
        "title",
        "experimental_methods_json",
        "resolution_combined_json",
        "best_resolution",
        "deposition_date",
        "initial_release_date",
    )
    inserted = 0
    failures = 0

    for pdb_id in tqdm(pdb_ids, desc="Fetching entry metadata"):
        try:
            url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
            data = fetch_json(session=session, url=url)
            values = extract_entry_metadata(pdb_id=pdb_id, data=data)
            insert_row(
                connection=connection,
                table="entry_metadata",
                columns=columns,
                values=values,
            )
            inserted += 1
        except Exception as exc:
            failures += 1
            insert_metadata_failure(
                connection=connection,
                metadata_scope="entry",
                pdb_id=pdb_id,
                entity_id=None,
                error_message=str(exc),
            )

    return inserted, failures


def fetch_polymer_entity_metadata(
    connection: sqlite3.Connection,
    session: requests.Session,
    entity_ids: list[tuple[str, str]],
) -> tuple[int, int]:
    columns = (
        "pdb_id",
        "entity_id",
        "entity_name",
        "organism_scientific_names_json",
        "taxonomy_ids_json",
        "uniprot_accessions_json",
        "polymer_type",
        "sequence_length",
    )
    inserted = 0
    failures = 0

    for pdb_id, entity_id in tqdm(entity_ids, desc="Fetching entity metadata"):
        try:
            url = (
                "https://data.rcsb.org/rest/v1/core/polymer_entity/"
                f"{pdb_id}/{entity_id}"
            )
            data = fetch_json(session=session, url=url)
            values = extract_polymer_entity_metadata(
                pdb_id=pdb_id,
                entity_id=entity_id,
                data=data,
            )
            insert_row(
                connection=connection,
                table="polymer_entity_metadata",
                columns=columns,
                values=values,
            )
            inserted += 1
        except Exception as exc:
            failures += 1
            insert_metadata_failure(
                connection=connection,
                metadata_scope="polymer_entity",
                pdb_id=pdb_id,
                entity_id=entity_id,
                error_message=str(exc),
            )

    return inserted, failures


def fetch_sqlite_metadata(sqlite_path: Path) -> MetadataStats:
    connection = sqlite3.connect(sqlite_path)
    rebuild_metadata_schema(connection)

    entry_ids = read_entry_ids(connection)
    entity_ids = read_polymer_entity_ids(connection)

    session = configure_session()
    entry_inserted, entry_failures = fetch_entry_metadata(
        connection=connection,
        session=session,
        pdb_ids=entry_ids,
    )
    entity_inserted, entity_failures = fetch_polymer_entity_metadata(
        connection=connection,
        session=session,
        entity_ids=entity_ids,
    )
    connection.commit()
    connection.close()

    return MetadataStats(
        entries_inserted=entry_inserted,
        entities_inserted=entity_inserted,
        failures_inserted=entry_failures + entity_failures,
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Fetch RCSB core metadata for selected LMDB export rows."
    )
    parser.add_argument(
        "sqlite_path",
        nargs="?",
        type=Path,
        default=DEFAULT_SQLITE_PATH,
        help="Path to the SQLite database (default: data/pdb_mldata.sqlite)",
    )
    args = parser.parse_args()

    try:
        validate_parameters(sqlite_path=args.sqlite_path)
        stats = fetch_sqlite_metadata(sqlite_path=args.sqlite_path)
    except (sqlite3.Error, ValueError) as exc:
        parser.error(str(exc))

    print(f"Wrote metadata to {args.sqlite_path}")
    print(f"Entry metadata inserted: {stats.entries_inserted}")
    print(f"Polymer entity metadata inserted: {stats.entities_inserted}")
    print(f"Metadata fetch failures inserted: {stats.failures_inserted}")


if __name__ == "__main__":
    main()
