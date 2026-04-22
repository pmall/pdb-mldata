# PDB MLData

PDB MLData is a project for building a structured collection of peptide/receptor
pairs from the Protein Data Bank. It starts from public RCSB PDB entries,
focuses on curated biological assemblies, applies conservative structural
filters, and represents the accepted pairs in a compact format for
machine-learning workflows.

The project is intentionally defensive: ambiguous structural neighborhoods are
rejected instead of repaired, inferred, or silently merged. The goal is a dataset
whose inclusion rules are easy to audit and whose failures can be studied before
the filters are widened.

Live viewer for selected chain pairs:
https://pdb-viewer-ashy.vercel.app/

## What Counts As A Pair

A peptide is a protein entity whose sequence is 4 to 32 amino acids long after
trimming common terminal caps. The initial metadata query uses a wider 4 to 40
amino-acid window so capped peptides that trim down into range are not missed.

A receptor is another protein entity within 5 Angstroms of a peptide chain. The
peptide chain is kept only when it has exactly one meaningful neighboring chain
and that neighbor is a protein chain. Water and explicitly listed
crystallographic solvent or buffer molecules do not count as competing
neighbors.

Saved data follows the natural PDB hierarchy:

```text
PDB entry
  peptide entity
    peptide-chain / receptor-chain pair
```

For each saved chain, the dataset stores:

- The normalized one-letter sequence.
- Exact retained 3-letter residue names after cap trimming.
- 37 atom positions per residue in the AlphaFold atom order.
- Per-atom B-factors.
- Per-atom occupancy values.

Non-standard amino acids are retained. When Gemmi cannot convert a residue to a
one-letter amino-acid code, the normalized sequence uses `X`, while the exact
3-letter residue name remains available.

## Data Artifacts

The main dataset is an LMDB database keyed by PDB ID. Each entry contains the
accepted peptide entities for that PDB structure, and each peptide entity
contains one or more observed peptide-chain/receptor-chain pairs.

The intermediate artifacts mirror the project pipeline:

- A metadata table of candidate PDB entries.
- A ZIP archive of first biological assembly mmCIF files.
- The final peptide/receptor LMDB database.
- Optional viewer databases derived from the LMDB.

## LMDB Schema

LMDB keys are PDB IDs. Values are `msgpack` payloads with this structure:

```json
{
  "pdb_id": "",
  "entities": [
    {
      "entity_id": "",
      "sequence": "",
      "residue_names": [],
      "pairs": [
        {
          "peptide": {
            "entity_id": "",
            "chain": "",
            "sequence": "",
            "residue_names": [],
            "structure": "<bytes>",
            "b_factors": "<bytes>",
            "occupancy": "<bytes>"
          },
          "receptor": {
            "entity_id": "",
            "chain": "",
            "sequence": "",
            "residue_names": [],
            "structure": "<bytes>",
            "b_factors": "<bytes>",
            "occupancy": "<bytes>"
          }
        }
      ]
    }
  ]
}
```

Binary array conventions:

- `structure`: `float16 [N, 37, 3]`, serialized with `.tobytes()`.
- `b_factors`: `uint8 [N, 37]`, serialized with `.tobytes()`.
- `occupancy`: `uint8 [N, 37]`, occupancy multiplied by 100 before storage.
- `255`: missing-value sentinel for `uint8` quality arrays.

The shared LMDB utilities define the encoding and decoding conventions so that
the structural byte arrays can be restored as typed coordinate, B-factor, and
occupancy arrays.

## Viewer Database

The viewer database is a derived sidecar artifact for search and publication. It
projects selected LMDB metadata into relational tables, enriches entries and
entities with RCSB and UniProt metadata, and derives target-focused search terms
for autocomplete and trigram lookup.

The viewer layer is intentionally separate from the LMDB. It is a publication
and exploration surface, not the source of truth for peptide/receptor pair
selection.

## Data Policy

PDB MLData favors clear rejection over hidden repair:

- Ambiguous peptide neighborhoods are skipped.
- Multi-model entries, including NMR entries, are skipped.
- Duplicate peptide-chain/receptor-chain pair keys cause the whole entry to be
  skipped as malformed.
- Entity sequences remain ideal PDB entity sequences.
- Chain sequences and structures represent observed biological assembly chains.
- Raw PDB structures are not rewritten to force entries into the dataset.
- Required parser data must come from Gemmi parsed objects or documented APIs.

This policy is especially important for downstream machine-learning use: a
smaller, auditable dataset is preferable to a larger dataset with unclear rescue
logic.

## Status

The data-building pipeline is active and intentionally conservative. The filter
set will evolve as rejection logs and downstream reports are studied.
