from __future__ import annotations

import argparse
import json
import sqlite3
from dataclasses import dataclass
from pathlib import Path
from typing import TypeAlias, TypeVar

import requests
from requests.adapters import HTTPAdapter
from tqdm import tqdm
from urllib3.util.retry import Retry

DEFAULT_SQLITE_PATH = Path("data/pdb_mldata.sqlite")
DEFAULT_TIMEOUT = 20.0
DEFAULT_BATCH_SIZE = 100
UNIPROT_SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"

JsonObject: TypeAlias = dict[str, object]
T = TypeVar("T")


@dataclass(frozen=True)
class EntityAccession:
    pdb_id: str
    entity_id: str
    accession: str


@dataclass(frozen=True)
class MetadataStats:
    accessions_inserted: int
    peptide_mappings_inserted: int
    target_mappings_inserted: int
    failures: int


def validate_parameters(sqlite_path: Path, batch_size: int) -> None:
    """Validate CLI parameters before opening the viewer database."""
    if not sqlite_path.exists():
        raise ValueError(f"{sqlite_path} not found")
    if sqlite_path.is_dir():
        raise ValueError("sqlite_path must be a file path, not a directory")
    if batch_size < 1:
        raise ValueError("--batch-size must be at least 1")


def json_text(value: object) -> str:
    return json.dumps(value, separators=(",", ":"), sort_keys=True)


def configure_session() -> requests.Session:
    session = requests.Session()
    retries = Retry(total=3, backoff_factor=0.5, status_forcelist=[500, 502, 503, 504])
    session.mount("https://", HTTPAdapter(max_retries=retries))
    return session


def batched(values: list[T], batch_size: int) -> list[list[T]]:
    return [
        values[index : index + batch_size]
        for index in range(0, len(values), batch_size)
    ]


def accession_query(accessions: list[str]) -> str:
    return "(" + " OR ".join(f"accession:{accession}" for accession in accessions) + ")"


def fetch_search_results(
    session: requests.Session,
    accessions: list[str],
) -> list[JsonObject]:
    response = session.get(
        UNIPROT_SEARCH_URL,
        params={
            "query": accession_query(accessions),
            "format": "json",
            "size": str(len(accessions)),
        },
        timeout=DEFAULT_TIMEOUT,
    )
    response.raise_for_status()
    data = response.json()
    if not isinstance(data, dict):
        raise ValueError("UniProt response was not a JSON object")
    results = data.get("results")
    if not isinstance(results, list):
        raise ValueError("UniProt search response did not contain a results list")
    return [result for result in results if isinstance(result, dict)]


def rebuild_uniprot_schema(connection: sqlite3.Connection) -> None:
    """Delete only UniProt tables before rebuilding metadata from accessions."""
    connection.executescript(
        """
        PRAGMA foreign_keys = ON;

        DROP TABLE IF EXISTS search_terms_targets;
        DROP TABLE IF EXISTS search_terms;
        DROP TABLE IF EXISTS peptides_accessions;
        DROP TABLE IF EXISTS targets_accessions;
        DROP TABLE IF EXISTS uniprot_targets;
        DROP TABLE IF EXISTS uniprot_metadata;

        CREATE TABLE uniprot_metadata (
            accession TEXT PRIMARY KEY,
            reviewed INTEGER NOT NULL,
            recommended_name TEXT NOT NULL,
            gene_names_json TEXT NOT NULL,
            organism_scientific_name TEXT NOT NULL,
            taxonomy_id INTEGER,
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

        CREATE INDEX idx_peptides_accessions__accession
            ON peptides_accessions (accession);
        CREATE INDEX idx_targets_accessions__accession
            ON targets_accessions (accession);
        """
    )
    connection.commit()


def first_string(data: JsonObject, key: str) -> str:
    value = data.get(key)
    return value if isinstance(value, str) else ""


def first_int(data: JsonObject, key: str) -> int | None:
    value = data.get(key)
    return value if isinstance(value, int) else None


def object_list(value: object) -> list[JsonObject]:
    if not isinstance(value, list):
        return []
    return [item for item in value if isinstance(item, dict)]


def nested_string(data: JsonObject, path: tuple[str, ...]) -> str:
    current: object = data
    for key in path:
        if not isinstance(current, dict):
            return ""
        current = current.get(key)
    return current if isinstance(current, str) else ""


def protein_name(data: JsonObject) -> str:
    description = data.get("proteinDescription")
    if not isinstance(description, dict):
        return ""

    recommended_name = nested_string(
        description,
        ("recommendedName", "fullName", "value"),
    )
    if recommended_name:
        return recommended_name

    submission_names = object_list(description.get("submissionNames"))
    for submission_name in submission_names:
        name = nested_string(submission_name, ("fullName", "value"))
        if name:
            return name

    return ""


def gene_names(data: JsonObject) -> list[str]:
    names: set[str] = set()
    for gene in object_list(data.get("genes")):
        gene_name = nested_string(gene, ("geneName", "value"))
        if gene_name:
            names.add(gene_name)
        for collection_name in ("synonyms", "orderedLocusNames", "orfNames"):
            for item in object_list(gene.get(collection_name)):
                value = first_string(item, "value")
                if value:
                    names.add(value)
    return sorted(names)


def comment_text(data: JsonObject, comment_type: str) -> str:
    texts: list[str] = []
    for comment in object_list(data.get("comments")):
        if first_string(comment, "commentType") != comment_type:
            continue
        for text in object_list(comment.get("texts")):
            value = first_string(text, "value")
            if value:
                texts.append(value)
    return "\n".join(dict.fromkeys(texts))


def subcellular_locations(data: JsonObject) -> list[str]:
    locations: set[str] = set()
    for comment in object_list(data.get("comments")):
        if first_string(comment, "commentType") != "SUBCELLULAR LOCATION":
            continue
        for location_group in object_list(comment.get("subcellularLocations")):
            for key in ("location", "topology", "orientation"):
                value = nested_string(location_group, (key, "value"))
                if value:
                    locations.add(value)
    return sorted(locations)


def keywords(data: JsonObject) -> list[str]:
    names = [
        first_string(keyword, "name") for keyword in object_list(data.get("keywords"))
    ]
    return sorted({name for name in names if name})


def cross_references(data: JsonObject, database: str) -> list[JsonObject]:
    refs: list[JsonObject] = []
    for reference in object_list(data.get("uniProtKBCrossReferences")):
        if first_string(reference, "database") == database:
            refs.append(reference)
    return refs


def property_value(reference: JsonObject, key: str) -> str:
    for item in object_list(reference.get("properties")):
        if first_string(item, "key") == key:
            return first_string(item, "value")
    return ""


def go_terms(data: JsonObject) -> list[JsonObject]:
    terms: list[JsonObject] = []
    for reference in cross_references(data=data, database="GO"):
        term_id = first_string(reference, "id")
        term_name = property_value(reference=reference, key="GoTerm")
        evidence = property_value(reference=reference, key="GoEvidenceType")
        terms.append({"id": term_id, "term": term_name, "evidence": evidence})
    return sorted(terms, key=lambda term: str(term["id"]))


def cross_reference_ids(data: JsonObject, database: str) -> list[str]:
    ids = [
        first_string(reference, "id")
        for reference in cross_references(data=data, database=database)
    ]
    return sorted({reference_id for reference_id in ids if reference_id})


def extract_uniprot_metadata(accession: str, data: JsonObject) -> tuple[object, ...]:
    organism = data.get("organism")
    organism_data = organism if isinstance(organism, dict) else {}
    entry_type = first_string(data, "entryType")

    return (
        accession,
        int(entry_type.lower().startswith("uniprotkb reviewed")),
        protein_name(data),
        json_text(gene_names(data)),
        first_string(organism_data, "scientificName"),
        first_int(organism_data, "taxonId"),
        comment_text(data=data, comment_type="FUNCTION"),
        json_text(subcellular_locations(data)),
        json_text(keywords(data)),
        json_text(go_terms(data)),
        json_text(cross_reference_ids(data=data, database="InterPro")),
        json_text(cross_reference_ids(data=data, database="Pfam")),
    )


def primary_accession(data: JsonObject) -> str:
    accession = first_string(data, "primaryAccession")
    if not accession:
        raise ValueError("UniProt result did not contain primaryAccession")
    return accession


def secondary_accessions(data: JsonObject) -> list[str]:
    values = data.get("secondaryAccessions")
    if not isinstance(values, list):
        return []
    return [value for value in values if isinstance(value, str)]


def result_accessions(data: JsonObject) -> set[str]:
    return {primary_accession(data), *secondary_accessions(data)}


def parse_accessions(value: str) -> list[str]:
    data = json.loads(value)
    if not isinstance(data, list):
        raise ValueError("accessions_json must contain a JSON list")
    accessions = [accession for accession in data if isinstance(accession, str)]
    if len(accessions) != len(data):
        raise ValueError("accessions_json must contain only strings")
    return sorted(set(accessions))


def read_entity_accessions(
    connection: sqlite3.Connection,
    table: str,
) -> list[EntityAccession]:
    mappings: list[EntityAccession] = []
    rows = connection.execute(
        f"""
        SELECT pdb_id, entity_id, accessions_json
        FROM {table}
        ORDER BY pdb_id, entity_id
        """
    ).fetchall()
    for pdb_id, entity_id, accessions_json in rows:
        if not isinstance(pdb_id, str) or not isinstance(entity_id, str):
            raise ValueError(f"{table} has malformed entity keys")
        if not isinstance(accessions_json, str):
            raise ValueError(f"{table}.accessions_json must be text")
        for accession in parse_accessions(accessions_json):
            mappings.append(
                EntityAccession(
                    pdb_id=pdb_id,
                    entity_id=entity_id,
                    accession=accession,
                )
            )
    return mappings


def insert_row(
    connection: sqlite3.Connection,
    table: str,
    columns: tuple[str, ...],
    values: tuple[object, ...],
) -> None:
    placeholders = ", ".join("?" for _ in columns)
    sql = f"INSERT INTO {table} ({', '.join(columns)}) VALUES ({placeholders})"
    connection.execute(sql, values)


def fetch_accession_metadata(
    connection: sqlite3.Connection,
    session: requests.Session,
    accessions: list[str],
    batch_size: int,
) -> tuple[dict[str, str], int, int]:
    columns = (
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
    )
    resolved_accessions: dict[str, str] = {}
    inserted_primary_accessions: set[str] = set()
    failures = 0

    with tqdm(
        total=len(accessions),
        desc="Fetching UniProt metadata",
        unit="accession",
    ) as progress:
        for accession_batch in batched(accessions, batch_size):
            try:
                rows = fetch_search_results(session=session, accessions=accession_batch)
            except Exception:
                failures += len(accession_batch)
                progress.update(len(accession_batch))
                continue

            batch_resolved_accessions: set[str] = set()
            for row in rows:
                try:
                    primary = primary_accession(row)
                    matching_accessions = sorted(
                        set(accession_batch) & result_accessions(row)
                    )
                    if not matching_accessions:
                        continue

                    if primary not in inserted_primary_accessions:
                        values = extract_uniprot_metadata(accession=primary, data=row)
                        insert_row(
                            connection=connection,
                            table="uniprot_metadata",
                            columns=columns,
                            values=values,
                        )
                        inserted_primary_accessions.add(primary)

                    for accession in matching_accessions:
                        resolved_accessions[accession] = primary
                        batch_resolved_accessions.add(accession)
                except Exception:
                    failures += 1

            failures += len(set(accession_batch) - batch_resolved_accessions)
            progress.update(len(accession_batch))

    return resolved_accessions, len(inserted_primary_accessions), failures


def insert_accession_mappings(
    connection: sqlite3.Connection,
    table: str,
    mappings: list[EntityAccession],
    resolved_accessions: dict[str, str],
) -> int:
    columns = ("pdb_id", "entity_id", "accession")
    rows: set[tuple[str, str, str]] = set()
    for mapping in mappings:
        primary = resolved_accessions.get(mapping.accession)
        if primary is None:
            continue
        rows.add((mapping.pdb_id, mapping.entity_id, primary))

    for row in sorted(rows):
        insert_row(
            connection=connection,
            table=table,
            columns=columns,
            values=row,
        )
    return len(rows)


def fetch_uniprot_metadata(sqlite_path: Path, batch_size: int) -> MetadataStats:
    connection = sqlite3.connect(sqlite_path)
    peptide_mappings = read_entity_accessions(
        connection=connection,
        table="peptides_metadata",
    )
    target_mappings = read_entity_accessions(
        connection=connection,
        table="targets_metadata",
    )
    accessions = sorted(
        {mapping.accession for mapping in [*peptide_mappings, *target_mappings]}
    )
    rebuild_uniprot_schema(connection)

    session = configure_session()
    resolved_accessions, metadata_inserted, failures = fetch_accession_metadata(
        connection=connection,
        session=session,
        accessions=accessions,
        batch_size=batch_size,
    )
    peptide_mappings_inserted = insert_accession_mappings(
        connection=connection,
        table="peptides_accessions",
        mappings=peptide_mappings,
        resolved_accessions=resolved_accessions,
    )
    target_mappings_inserted = insert_accession_mappings(
        connection=connection,
        table="targets_accessions",
        mappings=target_mappings,
        resolved_accessions=resolved_accessions,
    )
    connection.commit()
    connection.close()

    return MetadataStats(
        accessions_inserted=metadata_inserted,
        peptide_mappings_inserted=peptide_mappings_inserted,
        target_mappings_inserted=target_mappings_inserted,
        failures=failures,
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Fetch UniProt metadata for selected viewer database entities."
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
        help="Number of UniProt accessions per search request (default: 100)",
    )
    args = parser.parse_args()

    try:
        validate_parameters(sqlite_path=args.sqlite_path, batch_size=args.batch_size)
        stats = fetch_uniprot_metadata(
            sqlite_path=args.sqlite_path,
            batch_size=args.batch_size,
        )
    except (json.JSONDecodeError, sqlite3.Error, ValueError) as exc:
        parser.error(str(exc))

    print(f"Wrote UniProt metadata to {args.sqlite_path}")
    print(f"Accessions inserted: {stats.accessions_inserted}")
    print(f"UniProt peptide mappings inserted: {stats.peptide_mappings_inserted}")
    print(f"UniProt target mappings inserted: {stats.target_mappings_inserted}")
    print(f"UniProt metadata fetch failures: {stats.failures}")


if __name__ == "__main__":
    main()
