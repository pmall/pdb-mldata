# PDB MLData Project - Agent Instructions

This file contains project-specific instructions for AI agents working on this repository. Agents must follow these rules at all times.

## 1. Engineering Guidelines

### 1.1 Package And Tooling

- Always use `uv`.
- Never use `python` or `pip` directly.
- Run Python code with `uv run python`.
- Add or update dependencies only through `uv`.
- Ruff is the default Python formatter and import cleanup tool.
- Pyright is the default Python type checker.

### 1.2 Script And Module Structure

- Every executable script must define a `main()` function.
- `main()` must only parse CLI parameters, validate those parameters, and call a separate logic function.
- Business logic must live outside `main()` in modular, single-purpose functions.
- Each function must have one clear responsibility.
- Shared behavior must be moved into reusable library functions when it is used by multiple entrypoints.
- Reusable library code must live in importable package modules, not in `scripts/`.
- `scripts/` must contain executable entrypoints only.

### 1.3 CLI Parameters

- All file path parameters must be positional arguments.
- All file path parameters must have default values.
- All parameters must have default values unless they enable a truly optional feature.
- Default values for paths and settings must be consistent across the project.
- Validate CLI parameters in `main()` before calling the logic function.

### 1.4 Console Reporting

- Scripts should provide clear console reporting.
- Long-running operations should usually show a progress bar with progress and estimated time to completion.
- Console output should make skipped data, rejection counts, and completion status easy to understand.

### 1.5 Code Style

- Keep code modular and explicit.
- Prefer small, single-purpose functions over large procedural blocks.
- Prefer clear names over abbreviations.
- Keep data transformations explicit and local to the step that owns them.
- Follow the style already established in the repository unless it conflicts with this file.

### 1.6 Typing

- Typing must not be defensive.
- Avoid `Optional` unless `None` is a real, required state in the domain model.
- Avoid `Any`.
- Aim for exact types.
- Prefer concrete domain types, dataclasses, typed dictionaries, or explicit container types where they make the data contract clearer.

### 1.7 Comments

- Comments must explain why the code is doing something when that reasoning is not obvious from the code itself.
- Comments should be useful to future coding agents and future maintainers.
- Avoid comments that restate the code.
- Keep comments short and precise.
- Use a single comment style across the project: complete-sentence comments immediately before the non-obvious block they explain.

### 1.8 Code Updates

- Code updates must always be defensive.
- Do not update working parts of the code unless the user explicitly asks for that change.
- Keep edits scoped to the requested behavior.
- Do not refactor unrelated code opportunistically.
- Do not silently alter raw source data structures to make records fit expectations.

### 1.9 Verification

At the end of every coding session:

- Delete dead code and remove useless imports introduced or exposed by the change.
- Run `uv run ruff format`.
- Run `uv run ruff check --fix`.
- Run `uv run pyright`.
- Run `uv run python -m compileall -q scripts`.
- Check the diff and confirm it contains only intentional changes.

For data-processing or pipeline script verification:

- Prefer smoke tests over full pipeline runs during agent-side verification.
- When a script supports `--limit`, use `--limit` with a small number of records for agent-side verification.
- If a meaningful verification requires a full data-processing run, ask the developer to run it.
- Do not start long-running data processing jobs unless the user explicitly asks for a full run.

### 1.10 Git Workflow

- Never commit unless the user explicitly asks for a commit.
- Check the diff before preparing a commit message.
- Commit messages should usually start with a one-line explanation, followed by a bullet list of details.
- When using `git commit` from the command line, do not pass each bullet as a separate `-m` argument because that creates extra blank lines between bullets.
- Commit messages must focus on meaningful information for other developers.
- Commit messages must describe meaningful changes since the last commit.
- Do not include back-and-forth session details in commit messages.
- Do not focus commit messages on low-level implementation details unless those details affect other developers.
- Never add a co-author.

## 2. Defensive Data Science Rules

### 2.1 Core Philosophy

- This project follows defensive data science.
- Ambiguity means reject the data.
- Uncertainty means reject the data.
- If an entry does not fit the established rules, either stop or log and skip it.
- Start with hard filters that keep only certain data.
- Log and study rejection causes before widening filters.
- Document every change in data strategy in `AGENTS.md` before implementing it.

### 2.2 Data Handling

- Never manually merge chains or alter data to force entries to fit a mental model.
- Never make data manipulation decisions silently.
- Never apply transformations to the raw PDB data structure itself.
- When a parser is used for a data format, use only data exposed by that parser's parsed objects or documented parser APIs.
- Do not manually parse raw structured-data fields to rescue missing values when a parser is already responsible for that format.
- If required data is not available through the parser, stop and ask before adding fallback parsing or inference.

### 2.3 Public And Internal Data

- Ingestion from public data may reject or skip entries according to documented rules.
- Consumers of internal data must not repair, infer, deduplicate, or silently skip malformed internal data.
- Any inconsistency in internal data is a bug signal; downstream scripts must fail clearly so the source can be fixed.

## 3. Domain Overview

### 3.1 Project Goal

This project builds an exhaustive repository of peptide/receptor pairs from the whole Protein Data Bank.

A peptide is a protein entity with length between 4 and 32 amino acids after trimming usual N-terminal and C-terminal caps.

A receptor is another protein entity within 5 Angstroms of the peptide entity. It must be the only other meaningful chain in the peptide neighborhood. Meaningless ligand chains, such as water and the ligand chains defined in the LMDB-building script, do not count as receptor neighbors.

The source PDB structure is hierarchical:

- PDB entry: entry-level record.
- Peptide entity: one ideal PDB entity sequence.
- Pair: one observed peptide-chain/receptor-chain realization.

### 3.2 Documentation Locations

- Storage schemas and binary encoding conventions live in `docs/storage_schemas.md`.
- Shared LMDB encoder and decoder implementations live in `pdb_mldata/lmdb_utils.py`.
- Detailed script-specific data strategy lives in each script's top module docstring.
- Broad data-strategy changes must be reflected here; script-specific rule changes must be documented in the owning script docstring before implementation.
- Schema changes must be documented in `docs/storage_schemas.md` before implementation.

## 4. Pipeline

Core data-building scripts live under `scripts/`. Subjective dataset-selection scripts live under `scripts/curation/`. Viewer database scripts live under `scripts/viewer/`; see `scripts/viewer/README.md` for sidecar viewer database details.

| Stage | Script | Input | Output | Details |
| --- | --- | --- | --- | --- |
| Metadata fetch | `scripts/fetch_metadata.py` | RCSB API | `data/metadata.csv` | Script docstring |
| Assembly download | `scripts/download_assemblies.py` | `data/metadata.csv` | `data/assemblies.zip` | Script docstring |
| Raw LMDB build | `scripts/build_lmdb.py` | `data/assemblies.zip` | `data/pdb_mldata.lmdb` | Script docstring |
| Binding curation | `scripts/curation/filter_binding_pairs.py` | `data/pdb_mldata.lmdb` | `data/pdb_mldata_binding.lmdb` | Script docstring |
| Best-pair curation | `scripts/curation/select_best_pairs.py` | `data/pdb_mldata_binding.lmdb` | `data/pdb_mldata_best_pair.lmdb` | Script docstring |
