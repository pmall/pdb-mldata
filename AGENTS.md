# PDB MLData Project - Agent Instructions

This file contains project-specific instructions for AI agents working on this repository. Agents must follow these rules at all times.

## 1. Code Guidelines

### 1.1 Package Manager

- Always use `uv`.
- Never use `python` or `pip` directly.
- Run Python code with `uv run python`.
- Add or update dependencies only through `uv`.

### 1.2 Script Structure

- Every executable script must define a `main()` function.
- `main()` must only parse CLI parameters, validate those parameters, and call a separate logic function.
- Business logic must live outside `main()` in modular, single-purpose functions.
- Each function must have one clear responsibility.
- Shared behavior must be moved into reusable library functions when it is used by multiple scripts.
- Encoding and decoding conventions for LMDB entries must live in one shared library location. Prefer a standard encoder and decoder instead of duplicating serialization logic across scripts.

### 1.3 CLI Parameters

- All file path parameters must be positional arguments.
- All file path parameters must have default values.
- All parameters must have default values unless they enable a truly optional feature.
- Default values for paths and settings must be consistent across the project.
- Optional functionality may omit a default when absence means "do not use this feature". Example: `--log-skip-path`.
- Validate CLI parameters in `main()` before calling the logic function.

### 1.4 Console Reporting

- Scripts should provide clear console reporting.
- Long-running operations should usually show a progress bar with progress and estimated time to completion.
- Console output should make skipped data, rejection counts, and completion status easy to understand.

### 1.5 Comment Style

- Comments must explain why the code is doing something when that reasoning is not obvious from the code itself.
- Comments should be useful to future coding agents and future maintainers.
- Avoid comments that restate the code.
- Keep comments short and precise.
- Use a single comment style across the project: complete-sentence comments immediately before the non-obvious block they explain.

### 1.6 Code Style

- Keep code modular and explicit.
- Prefer small, single-purpose functions over large procedural blocks.
- Prefer clear names over abbreviations.
- Keep data transformations explicit and local to the step that owns them.
- Follow the style already established in the repository unless it conflicts with this file.
- Ruff is the default Python formatter and import cleanup tool.
- After every Python code update, run `uv run ruff format`.
- After every Python code update, run `uv run ruff check --fix`.
- Pyright is the default Python type checker.
- After every Python code update, run `uv run pyright`.

### 1.7 Typing

- Typing must not be defensive.
- Avoid `Optional` unless `None` is a real, required state in the domain model.
- Avoid `Any`.
- Aim for exact types.
- Prefer concrete domain types, dataclasses, typed dictionaries, or explicit container types where they make the data contract clearer.

### 1.8 Code Updates

- Code updates must always be defensive.
- Do not update working parts of the code unless the user explicitly asks for that change.
- Keep edits scoped to the requested behavior.
- Do not refactor unrelated code opportunistically.
- Do not silently alter raw PDB data structures to make entries fit expectations.
- At the end of every coding update, delete dead code and remove useless imports introduced or exposed by the change.

### 1.9 Testing And Script Runs

- Avoid running full pipeline scripts for testing because they are long-running.
- When a script supports `--limit`, use `--limit` with a small number of entries for agent-side verification.
- If a meaningful verification requires a full script run, ask the developer to run it.
- Do not start long-running data processing jobs unless the user explicitly asks for a full run.

### 1.10 Git Workflow

- Never commit unless the user explicitly asks for a commit.
- Check the diff before preparing a commit message.
- Commit messages should usually start with a one-line explanation, followed by a bullet list of details.
- Commit messages must focus on meaningful information for other developers.
- Commit messages must describe meaningful changes since the last commit.
- Do not include back-and-forth session details in commit messages.
- Do not focus commit messages on low-level implementation details unless those details affect other developers.
- Never add a co-author.

## 2. Defensive Data Science Philosophy

- This project follows defensive data science.
- Ambiguity means reject the data.
- Uncertainty means reject the data.
- If an entry does not fit the established rules, either stop or log and skip it.
- Never manually merge chains or alter data to force entries to fit a mental model.
- Never make data manipulation decisions silently.
- Never apply transformations to the raw PDB data structure itself.
- Start with hard filters that keep only certain data.
- Log and study rejection causes before widening filters.
- Document every change in data strategy in `AGENTS.md` before implementing it.

## 3. Overall Project Goal

This project builds an exhaustive repository of peptide/receptor pairs from the whole Protein Data Bank.

A peptide is a protein entity with length between 4 and 32 amino acids after trimming usual N-terminal and C-terminal caps.

A receptor is another protein entity within 5 Angstroms of the peptide entity. It must be the only other meaningful chain in the peptide neighborhood. Meaningless ligand chains, such as water and the ligand chains defined in the LMDB-building script, do not count as receptor neighbors.

The final structure is hierarchical because of the structure of the PDB:

- PDB entry: entry-level record.
- Peptide entity: a PDB entry may have many peptide entities with at least one valid receptor.
- Pair: a peptide entity may have many chain-level peptide/receptor pairs.

For sequences, trim the usual caps on both C-terminus and N-terminus.

For structures, extract the 37 amino-acid atom 3D positions and store them sorted according to the AlphaFold convention.

Store residue-level quality values:

- B-factor.
- Occupancy.

Storage conventions:

- Structure coordinates are stored as `float16` bytes.
- B-factors are stored as `uint8` bytes.
- Occupancy values are multiplied by 100 and stored as `uint8` bytes.
- The value `255` is the `NaN`/`None` sentinel for `uint8` arrays.
- Entries are serialized with `msgpack.packb`.

The LMDB entry schema is:

```json
{
  "pdb_id": "",
  "entities": [
    {
      "entity_id": "",
      "sequence": "",
      "pairs": [
        {
          "peptide": {
            "entity_id": "",
            "chain": "",
            "sequence": "",
            "structure": "<bytes>",
            "b_factors": "<bytes>",
            "occupancy": "<bytes>"
          },
          "receptor": {
            "entity_id": "",
            "chain": "",
            "sequence": "",
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

Schema details:

- `pdb_id`: PDB identifier.
- `entities`: one entry per peptide entity.
- `entity_id`: PDB entity identifier.
- `sequence`: entity or raw chain sequence, depending on the object.
- `pairs`: one entry per peptide-chain/receptor-chain pair.
- `structure`: `float16 [N, 37, 3]` serialized with `.tobytes()`.
- `b_factors`: `uint8 [N, 37]` serialized with `.tobytes()`, where `255` means missing.
- `occupancy`: `uint8 [N, 37]` serialized with `.tobytes()`, where occupancy was multiplied by 100 and `255` means missing.

## 4. Current Script Pipeline

There are currently three scripts. More scripts may be added later.

### 4.1 `fetch_metadata`

Goal: query the PDB API and write all matching PDB entity IDs to a metadata CSV file.

Selection criteria:

- PDB entry has at least two protein entities.
- PDB entry has at least one peptide entity.
- For this query only, a peptide is a protein entity with length between 4 and 40 amino acids.
- The upper bound is 40 instead of 32 to include 32 amino-acid peptides with possible C-terminal or N-terminal caps.
- The value 40 may be adjusted later after studying logs from downstream scripts.

Parameters:

- Metadata CSV file path: `data/metadata.csv`.
- Minimum peptide length: `4`.
- Maximum peptide length: `40`.

### 4.2 `download_assemblies`

Goal: download the PDB entries retrieved by `fetch_metadata`.

Download rules:

- Download the first biological assembly only.
- Use the first biological assembly because it is curated data.
- Store all downloaded gz files in one ZIP archive.
- Do not write thousands of assembly files directly into a folder.

Parameters:

- Metadata CSV file path: `data/metadata.csv`.
- Assemblies ZIP file path: `data/assemblies.zip`.

### 4.3 `build_lmdb`

Goal: parse all biological assemblies downloaded by `download_assemblies` with `gemmi` and build the LMDB database.

This script owns the main data filters. These filters will be updated carefully over time after rejection logs are studied.

Current rules:

- Read assemblies from `data/assemblies.zip`.
- Skip multi-model entries, including NMR entries.
- Trim caps from all protein entities.
- Keep trimmed protein entities as peptide entities when their length is between 4 and 32 amino acids.
- Keep peptide entities only when the trimmed sequence contains standard amino acids only.
- For each chain of each peptide entity, use `scipy.spatial.KDTree` to identify all other chains in a 5 Angstrom neighborhood.
- Water and explicitly defined non-polymer ligand chains do not count as neighbors.
- If the peptide chain has exactly one meaningful neighboring chain and that chain is a protein chain, save the pair under the peptide entity.
- Skip entries that do not satisfy these rules.
- Existing LMDB database folders must be deleted before writing a new LMDB database to avoid errors caused by LMDB upsert behavior.

Work in progress:

- Protein-chain processing will be updated in future iterations.
- Non-standard amino-acid handling will be updated in future iterations.
- Any strategy change must be documented in this file before implementation.

Parameters:

- Assemblies ZIP file path: `data/assemblies.zip`.
- LMDB database folder path: `data/pdb_mldata.lmdb`.
- Minimum peptide length: `4`.
- Maximum peptide length: `32`.
- Distance threshold: `5.0`.
- Optional processing limit for testing: process only `N` entries.
