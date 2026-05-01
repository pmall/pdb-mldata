# Viewer Database Scripts

These scripts build and enrich the read-only viewer database. The viewer database
is a sidecar artifact derived from the LMDB and can be rebuilt from scratch.

## Run Order

```bash
uv run python scripts/viewer/export_sqlite.py
uv run python scripts/viewer/fetch_rcsb_entries_metadata.py
uv run python scripts/viewer/fetch_rcsb_entities_metadata.py
uv run python scripts/viewer/fetch_uniprot_metadata.py
uv run python scripts/viewer/build_search_terms.py
uv run python scripts/viewer/dump_postgres_sql.py
```

Use `--limit` only with `export_sqlite.py` for small local verification exports.

## Scripts

### `export_sqlite.py`

Exports peptide entities and peptide/receptor chain pairs from `data/pdb_mldata.lmdb` into
`data/pdb_mldata.sqlite`.

Rules:

- Rebuild `entries`, `peptides`, and `chain_pairs`.
- Keep peptide entities as the only first-class entity table.
- Store peptide and receptor binding spans in `chain_pairs` for viewer sequence
  highlighting.

Default parameters:

- LMDB database folder path: `data/pdb_mldata.lmdb`.
- SQLite database file path: `data/pdb_mldata.sqlite`.
- Optional processing limit for testing: process only `N` entries.

### `fetch_rcsb_entries_metadata.py`

Fetches RCSB entry metadata for rows in the SQLite viewer database.

Rules:

- Rebuild only `entries_metadata`.
- Use batched RCSB GraphQL requests instead of per-row REST calls.
- Count metadata fetch failures in console output only.

Default parameters:

- SQLite database file path: `data/pdb_mldata.sqlite`.
- RCSB GraphQL batch size: `100`.

### `fetch_rcsb_entities_metadata.py`

Fetches RCSB entity metadata for peptide and target entities in the SQLite
viewer database.

Rules:

- Rebuild `peptides_metadata` and `targets_metadata`.
- Read peptide entity IDs from `chain_pairs.peptide_entity_id`.
- Read target entity IDs from `chain_pairs.receptor_entity_id`.
- Drop derived UniProt and search tables before replacing entity metadata.
- Use batched RCSB GraphQL requests instead of per-row REST calls.
- Use one progress bar for the total number of unique RCSB polymer entities to
  fetch.
- Store RCSB-provided UniProt accession lists in `accessions_json`.
- Count metadata fetch failures in console output only.

Default parameters:

- SQLite database file path: `data/pdb_mldata.sqlite`.
- RCSB GraphQL batch size: `100`.

### `fetch_uniprot_metadata.py`

Fetches UniProt metadata for accessions listed in `peptides_metadata` and
`targets_metadata`.

Rules:

- Rebuild only the UniProt metadata tables.
- Use batched UniProt search requests instead of one request per accession.
- Use one progress bar for the number of unique UniProt accessions to fetch.
- Store one row per fetched UniProt entry in `uniprot_metadata`, keyed by the
  primary accession returned by UniProt.
- Store peptide-to-primary-accession relationships in `peptides_accessions`.
- Store target-to-primary-accession relationships in `targets_accessions`.
- Use `accession` as the accession column name everywhere.
- Do not fetch or store UniProt sequence features yet.
- Count metadata fetch failures in console output only.

Default parameters:

- SQLite database file path: `data/pdb_mldata.sqlite`.
- UniProt accession batch size: `100`.

### `build_search_terms.py`

Builds derived autocomplete tables used for term lookup and trigram search in
PostgreSQL.

Rules:

- Rebuild only the derived search tables.
- Start from `targets_accessions`; search terms are target-only.
- Store one row per normalized text term in `search_terms`.
- Normalize search terms with whitespace cleanup and case folding.
- Keep the highest-weight `term_kind` when multiple sources produce the same
  normalized term.
- Map each term to matching targets in `search_terms_targets`.
- Include PDB IDs, entry titles, target names, primary UniProt accessions,
  UniProt recommended names, gene names, organisms, keywords, GO terms,
  InterPro IDs, and Pfam IDs.
- Do not include sequences, experimental methods, resolutions, dates, chain IDs,
  or pair IDs.

Default parameter:

- SQLite database file path: `data/pdb_mldata.sqlite`.

### `dump_postgres_sql.py`

Dumps the SQLite viewer database to PostgreSQL-compatible SQL at
`data/pdb_mldata.sql`.

Rules:

- Validate the SQLite schema before dumping.
- Include RCSB metadata, UniProt metadata, and entity-to-accession mapping tables.
- Include `pg_trgm` setup and a trigram GIN index for `search_terms`.
- Bulk-load table contents with PostgreSQL `COPY ... FROM stdin`.

Default parameters:

- SQLite database file path: `data/pdb_mldata.sqlite`.
- PostgreSQL-compatible SQL file path: `data/pdb_mldata.sql`.

## Metadata Tables

`peptides` stores one row per selected peptide entity from the LMDB.

`chain_pairs` stores one row per selected peptide-chain/receptor-chain pair from
the LMDB. It includes 1-based inclusive peptide and receptor binding spans
relative to the stored trimmed chain sequences.

`peptides_metadata.accessions_json` stores the accession list reported by RCSB
for each peptide entity.

`targets_metadata.accessions_json` stores the accession list reported by RCSB for
each target entity.

`uniprot_metadata` stores one row per successfully fetched UniProt entry. Its
primary key is the primary UniProt `accession`.

`peptides_accessions` maps peptide entities to primary UniProt accessions. Its
primary key is `(pdb_id, entity_id, accession)`.

`targets_accessions` maps target entities to primary UniProt accessions. Its
primary key is `(pdb_id, entity_id, accession)`.

`search_terms` stores one row per normalized autocomplete term. It has a row id
primary key and a unique constraint on `term`.

`search_terms_targets` maps each autocomplete term to one or more target entities.
It has a row id primary key and a unique constraint on `(search_term_id, pdb_id,
entity_id)`.
