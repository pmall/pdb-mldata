from __future__ import annotations

import json
from typing import TypeAlias, TypeVar

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

DEFAULT_TIMEOUT = 20.0
DEFAULT_BATCH_SIZE = 100
GRAPHQL_URL = "https://data.rcsb.org/graphql"

ENTRY_METADATA_QUERY = """
query EntryMetadata($entryIds: [String!]!) {
  entries(entry_ids: $entryIds) {
    rcsb_id
    struct {
      title
    }
    rcsb_entry_info {
      resolution_combined
    }
    rcsb_accession_info {
      initial_release_date
    }
    pdbx_database_status {
      recvd_initial_deposition_date
    }
    exptl {
      method
    }
  }
}
"""

POLYMER_ENTITY_METADATA_QUERY = """
query PolymerEntityMetadata($entityIds: [String!]!) {
  polymer_entities(entity_ids: $entityIds) {
    rcsb_id
    rcsb_polymer_entity {
      pdbx_description
    }
    entity_poly {
      rcsb_entity_polymer_type
      rcsb_sample_sequence_length
    }
    rcsb_entity_source_organism {
      scientific_name
      ncbi_taxonomy_id
    }
    rcsb_polymer_entity_container_identifiers {
      uniprot_ids
      reference_sequence_identifiers {
        database_name
        database_accession
      }
    }
  }
}
"""

JsonObject: TypeAlias = dict[str, object]
T = TypeVar("T")


def json_text(value: object) -> str:
    return json.dumps(value, separators=(",", ":"), sort_keys=True)


def configure_session() -> requests.Session:
    session = requests.Session()
    retries = Retry(total=3, backoff_factor=0.5, status_forcelist=[500, 502, 503, 504])
    session.mount("https://", HTTPAdapter(max_retries=retries))
    return session


def fetch_graphql(
    session: requests.Session,
    query: str,
    variables: JsonObject,
) -> JsonObject:
    response = session.post(
        GRAPHQL_URL,
        json={"query": query, "variables": variables},
        timeout=DEFAULT_TIMEOUT,
    )
    response.raise_for_status()
    data = response.json()
    if not isinstance(data, dict):
        raise ValueError("RCSB GraphQL response was not a JSON object")
    errors = data.get("errors")
    if errors:
        raise ValueError(f"RCSB GraphQL errors: {json_text(errors)}")
    payload = data.get("data")
    if not isinstance(payload, dict):
        raise ValueError("RCSB GraphQL response did not contain a data object")
    return payload


def batched(values: list[T], batch_size: int) -> list[list[T]]:
    return [
        values[index : index + batch_size]
        for index in range(0, len(values), batch_size)
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


def extract_accessions(identifiers: JsonObject) -> list[str]:
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

    accessions = (
        extract_accessions(identifiers) if isinstance(identifiers, dict) else []
    )

    return (
        pdb_id,
        entity_id,
        entity_name,
        json_text(sorted(set(organism_names))),
        json_text(sorted_ints(taxonomy_values)),
        json_text(accessions),
        polymer_type,
        sequence_length,
    )


def entity_graphql_id(pdb_id: str, entity_id: str) -> str:
    return f"{pdb_id}_{entity_id}"


def parse_entity_graphql_id(rcsb_id: str) -> tuple[str, str]:
    parts = rcsb_id.split("_", maxsplit=1)
    if len(parts) != 2 or not parts[0] or not parts[1]:
        raise ValueError(f"malformed RCSB polymer entity id: {rcsb_id}")
    return parts[0], parts[1]
