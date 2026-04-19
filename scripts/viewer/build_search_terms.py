from __future__ import annotations

import argparse
import json
import sqlite3
from dataclasses import dataclass
from pathlib import Path
from typing import TypeAlias

from tqdm import tqdm

DEFAULT_SQLITE_PATH = Path("data/pdb_mldata.sqlite")

JsonObject: TypeAlias = dict[str, object]

TERM_WEIGHTS = {
    "pdb_id": 110,
    "accession": 100,
    "gene_name": 95,
    "protein_name": 90,
    "target_name": 80,
    "entry_title": 70,
    "organism": 55,
    "keyword": 50,
    "go_term": 45,
    "interpro": 35,
    "pfam": 35,
}


@dataclass(frozen=True)
class SearchTarget:
    pdb_id: str
    entity_id: str


@dataclass(frozen=True)
class SearchTerm:
    term: str
    term_kind: str
    rank_weight: int


SearchIndex: TypeAlias = dict[str, tuple[SearchTerm, set[SearchTarget]]]


@dataclass(frozen=True)
class SearchTermStats:
    terms_inserted: int
    targets_inserted: int


def validate_parameters(sqlite_path: Path) -> None:
    """Validate CLI parameters before opening the viewer database."""
    if not sqlite_path.exists():
        raise ValueError(f"{sqlite_path} not found")
    if sqlite_path.is_dir():
        raise ValueError("sqlite_path must be a file path, not a directory")


def rebuild_search_schema(connection: sqlite3.Connection) -> None:
    """Delete only derived search tables before rebuilding autocomplete terms."""
    connection.executescript(
        """
        PRAGMA foreign_keys = ON;

        DROP TABLE IF EXISTS search_terms_targets;
        DROP TABLE IF EXISTS search_terms;

        CREATE TABLE search_terms (
            search_term_id INTEGER PRIMARY KEY,
            term TEXT NOT NULL,
            term_kind TEXT NOT NULL,
            rank_weight INTEGER NOT NULL,
            UNIQUE (term)
        );

        CREATE TABLE search_terms_targets (
            search_terms_target_id INTEGER PRIMARY KEY,
            search_term_id INTEGER NOT NULL,
            pdb_id TEXT NOT NULL,
            entity_id TEXT NOT NULL,
            FOREIGN KEY (search_term_id) REFERENCES search_terms (search_term_id),
            FOREIGN KEY (pdb_id) REFERENCES entries (pdb_id),
            FOREIGN KEY (pdb_id, entity_id)
                REFERENCES targets_metadata (pdb_id, entity_id),
            UNIQUE (search_term_id, pdb_id, entity_id)
        );

        CREATE INDEX idx_search_terms__term_kind
            ON search_terms (term_kind);
        CREATE INDEX idx_search_terms_targets__search_term_id
            ON search_terms_targets (search_term_id);
        CREATE INDEX idx_search_terms_targets__pdb_id_entity_id
            ON search_terms_targets (pdb_id, entity_id);
        """
    )
    connection.commit()


def parse_string_list(value: str, column_name: str) -> list[str]:
    data = json.loads(value)
    if not isinstance(data, list):
        raise ValueError(f"{column_name} must contain a JSON list")
    strings = [item for item in data if isinstance(item, str)]
    if len(strings) != len(data):
        raise ValueError(f"{column_name} must contain only strings")
    return strings


def parse_go_terms(value: str) -> list[str]:
    data = json.loads(value)
    if not isinstance(data, list):
        raise ValueError("go_terms_json must contain a JSON list")

    terms: list[str] = []
    for item in data:
        if not isinstance(item, dict):
            raise ValueError("go_terms_json must contain only objects")
        term_id = item.get("id")
        term = item.get("term")
        if isinstance(term_id, str) and term_id:
            terms.append(term_id)
        if isinstance(term, str) and term:
            terms.append(term)
    return terms


def clean_term(term: str) -> str:
    return " ".join(term.split()).casefold()


def add_term(
    search_index: SearchIndex,
    term: str,
    term_kind: str,
    target: SearchTarget,
) -> None:
    clean = clean_term(term)
    if not clean:
        return
    candidate = SearchTerm(
        term=clean,
        term_kind=term_kind,
        rank_weight=TERM_WEIGHTS[term_kind],
    )
    existing = search_index.get(clean)
    if existing is None:
        search_index[clean] = (candidate, {target})
        return

    existing_term, targets = existing
    if candidate.rank_weight > existing_term.rank_weight:
        targets.add(target)
        search_index[clean] = (candidate, targets)
        return

    targets.add(target)


def read_targets(
    connection: sqlite3.Connection,
) -> dict[str, list[SearchTarget]]:
    rows = connection.execute(
        """
        SELECT DISTINCT pdb_id, entity_id
        FROM targets_accessions
        ORDER BY pdb_id, entity_id
        """
    ).fetchall()

    targets_by_pdb_id: dict[str, list[SearchTarget]] = {}
    for pdb_id, entity_id in rows:
        if not isinstance(pdb_id, str) or not isinstance(entity_id, str):
            raise ValueError("targets_accessions has malformed target keys")
        targets_by_pdb_id.setdefault(pdb_id, []).append(
            SearchTarget(pdb_id=pdb_id, entity_id=entity_id)
        )
    return targets_by_pdb_id


def add_entry_terms(
    connection: sqlite3.Connection,
    search_index: SearchIndex,
    targets_by_pdb_id: dict[str, list[SearchTarget]],
) -> None:
    rows = connection.execute(
        """
        SELECT entries.pdb_id, entries_metadata.title
        FROM entries
        LEFT JOIN entries_metadata ON entries_metadata.pdb_id = entries.pdb_id
        ORDER BY entries.pdb_id
        """
    ).fetchall()
    for pdb_id, title in tqdm(rows, desc="Building entry search terms"):
        if not isinstance(pdb_id, str):
            raise ValueError("entries.pdb_id must be text")
        targets = targets_by_pdb_id.get(pdb_id, [])
        for target in targets:
            add_term(
                search_index=search_index,
                term=pdb_id,
                term_kind="pdb_id",
                target=target,
            )
            if isinstance(title, str):
                add_term(
                    search_index=search_index,
                    term=title,
                    term_kind="entry_title",
                    target=target,
                )


def add_target_terms(
    connection: sqlite3.Connection,
    search_index: SearchIndex,
) -> None:
    rows = connection.execute(
        """
        SELECT DISTINCT
            targets_metadata.pdb_id,
            targets_metadata.entity_id,
            targets_metadata.entity_name
        FROM targets_accessions
        JOIN targets_metadata
            ON targets_metadata.pdb_id = targets_accessions.pdb_id
            AND targets_metadata.entity_id = targets_accessions.entity_id
        ORDER BY targets_metadata.pdb_id, targets_metadata.entity_id
        """
    ).fetchall()
    for pdb_id, entity_id, entity_name in tqdm(
        rows,
        desc="Building target search terms",
    ):
        if not isinstance(pdb_id, str) or not isinstance(entity_id, str):
            raise ValueError("targets_metadata has malformed target keys")
        if not isinstance(entity_name, str):
            raise ValueError("targets_metadata.entity_name must be text")

        target = SearchTarget(pdb_id=pdb_id, entity_id=entity_id)
        add_term(
            search_index=search_index,
            term=entity_name,
            term_kind="target_name",
            target=target,
        )


def add_uniprot_terms(
    connection: sqlite3.Connection,
    search_index: SearchIndex,
) -> None:
    rows = connection.execute(
        """
        SELECT
            targets_accessions.pdb_id,
            targets_accessions.entity_id,
            uniprot_metadata.accession,
            uniprot_metadata.recommended_name,
            uniprot_metadata.gene_names_json,
            uniprot_metadata.organism_scientific_name,
            uniprot_metadata.keywords_json,
            uniprot_metadata.go_terms_json,
            uniprot_metadata.interpro_ids_json,
            uniprot_metadata.pfam_ids_json
        FROM targets_accessions
        JOIN uniprot_metadata
            ON uniprot_metadata.accession = targets_accessions.accession
        ORDER BY
            targets_accessions.pdb_id,
            targets_accessions.entity_id,
            uniprot_metadata.accession
        """
    ).fetchall()

    for row in tqdm(rows, desc="Building UniProt search terms"):
        (
            pdb_id,
            entity_id,
            accession,
            recommended_name,
            gene_names_json,
            organism_scientific_name,
            keywords_json,
            go_terms_json,
            interpro_ids_json,
            pfam_ids_json,
        ) = row
        if not isinstance(pdb_id, str) or not isinstance(entity_id, str):
            raise ValueError("targets_accessions has malformed target keys")
        if not isinstance(accession, str):
            raise ValueError("uniprot_metadata.accession must be text")

        target = SearchTarget(pdb_id=pdb_id, entity_id=entity_id)
        add_term(
            search_index=search_index,
            term=accession,
            term_kind="accession",
            target=target,
        )
        if isinstance(recommended_name, str):
            add_term(
                search_index=search_index,
                term=recommended_name,
                term_kind="protein_name",
                target=target,
            )
        if isinstance(organism_scientific_name, str):
            add_term(
                search_index=search_index,
                term=organism_scientific_name,
                term_kind="organism",
                target=target,
            )

        if not isinstance(gene_names_json, str):
            raise ValueError("uniprot_metadata.gene_names_json must be text")
        for gene_name in parse_string_list(gene_names_json, "gene_names_json"):
            add_term(
                search_index=search_index,
                term=gene_name,
                term_kind="gene_name",
                target=target,
            )

        if not isinstance(keywords_json, str):
            raise ValueError("uniprot_metadata.keywords_json must be text")
        for keyword in parse_string_list(keywords_json, "keywords_json"):
            add_term(
                search_index=search_index,
                term=keyword,
                term_kind="keyword",
                target=target,
            )

        if not isinstance(go_terms_json, str):
            raise ValueError("uniprot_metadata.go_terms_json must be text")
        for go_term in parse_go_terms(go_terms_json):
            add_term(
                search_index=search_index,
                term=go_term,
                term_kind="go_term",
                target=target,
            )

        if not isinstance(interpro_ids_json, str):
            raise ValueError("uniprot_metadata.interpro_ids_json must be text")
        for interpro_id in parse_string_list(interpro_ids_json, "interpro_ids_json"):
            add_term(
                search_index=search_index,
                term=interpro_id,
                term_kind="interpro",
                target=target,
            )

        if not isinstance(pfam_ids_json, str):
            raise ValueError("uniprot_metadata.pfam_ids_json must be text")
        for pfam_id in parse_string_list(pfam_ids_json, "pfam_ids_json"):
            add_term(
                search_index=search_index,
                term=pfam_id,
                term_kind="pfam",
                target=target,
            )


def insert_search_terms(
    connection: sqlite3.Connection,
    search_index: SearchIndex,
) -> tuple[int, int]:
    ordered_terms = sorted(
        [term for term, _targets in search_index.values()],
        key=lambda item: (item.term, item.term_kind),
    )
    term_rows = [
        (term.term, term.term_kind, term.rank_weight) for term in ordered_terms
    ]
    connection.executemany(
        """
        INSERT INTO search_terms (
            term,
            term_kind,
            rank_weight
        )
        VALUES (?, ?, ?)
        """,
        term_rows,
    )

    term_ids = {
        row[1]: row[0]
        for row in connection.execute(
            """
            SELECT search_term_id, term
            FROM search_terms
            """
        ).fetchall()
    }

    target_rows: list[tuple[int, str, str]] = []
    for term in ordered_terms:
        search_term_id = term_ids[term.term]
        _search_term, targets = search_index[term.term]
        for target in sorted(
            targets,
            key=lambda item: (item.pdb_id, item.entity_id),
        ):
            target_rows.append((search_term_id, target.pdb_id, target.entity_id))

    connection.executemany(
        """
        INSERT INTO search_terms_targets (
            search_term_id,
            pdb_id,
            entity_id
        )
        VALUES (?, ?, ?)
        """,
        target_rows,
    )
    return len(term_rows), len(target_rows)


def build_search_terms(sqlite_path: Path) -> SearchTermStats:
    connection = sqlite3.connect(sqlite_path)
    targets_by_pdb_id = read_targets(connection)
    rebuild_search_schema(connection)

    search_index: SearchIndex = {}
    add_entry_terms(
        connection=connection,
        search_index=search_index,
        targets_by_pdb_id=targets_by_pdb_id,
    )
    add_target_terms(connection=connection, search_index=search_index)
    add_uniprot_terms(connection=connection, search_index=search_index)

    terms_inserted, targets_inserted = insert_search_terms(
        connection=connection,
        search_index=search_index,
    )
    connection.commit()
    connection.close()

    return SearchTermStats(
        terms_inserted=terms_inserted,
        targets_inserted=targets_inserted,
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build derived viewer autocomplete search terms."
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
        stats = build_search_terms(sqlite_path=args.sqlite_path)
    except (json.JSONDecodeError, sqlite3.Error, ValueError) as exc:
        parser.error(str(exc))

    print(f"Wrote viewer search terms to {args.sqlite_path}")
    print(f"Search terms inserted: {stats.terms_inserted}")
    print(f"Search term targets inserted: {stats.targets_inserted}")


if __name__ == "__main__":
    main()
