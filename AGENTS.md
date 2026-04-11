# PDB MLData Project - Agent Instructions

This file contains specific instructions for AI agents working on this project. Agents must strictly adhere to these rules at all times.

## 1. Project Overview

This is a biological machine learning project designed to parse PDB (Protein Data Bank) files and extract exact peptide-to-binding-site pairs into an LMDB database.

## 2. Core Philosophy & Rules (CRITICAL)

- **Defensive Data Approach**: You must adopt a strictly defensive approach when handling data.
- **Strict Compliance**: If a data entry does not perfectly fit the established rules, you must either STOP or LOG and SKIP the entry.
- **No Silenced Fixes**: **NEVER** try to manually merge chains or alter the data to force entries to fit our mental model.
- **No Silent Decisions**: Never make data manipulation decisions silently.
- **Raw Data Integrity**: Never apply transformations to the raw PDB data structure itself.
- **Package Manager**: **Always use `uv`**. NEVER use `python` or `pip` commands directly.

## 3. Script Structure & Best Practices

All Python scripts structured for this pipeline must follow these architectural guidelines:

- **Parameter Parsing**: `argparse` (or equivalent parsing) must *only* happen inside a dedicated `main()` function. 
- **Business Logic Separation**: After parsing parameters, the `main` function must invoke a separate, dedicated function containing the core logic, passing the parsed arguments.
- **CLI Arguments**: Always use **positional arguments** for input files and output files/directories.
- **Defaults**: Always provide sensible default values for optional parameters (e.g., distance thresholds, logging levels).

## 4. Data Pipeline Steps

The pipeline consists of the following consecutive scripts, strictly implementing the logical steps below:

### Script 1: PDB Metadata Query

- **Goal**: Hit the RCSB API to fetch PDB entries satisfying specific biological criteria.
- **Condition**: Retrieve entries with _at least two protein chains_, where at least one chain (the peptide) has a length between **4 and 40 amino acids (inclusive)**. 
    - **Note**: We use a higher upper bound (40) for the initial query to account for terminal caps (like ACE or NH2) that RCSB counts toward length, ensuring we don't accidentally filter out 32-mer peptides.
- **Output**: Save the matching entries' metadata into a CSV file.

### Script 2: Biological Assembly Downloader

- **Goal**: Download the required CIF files for the entries listed in the CSV.
- **Condition**: Only download the **first Biological Assembly**, as it is considered the most reliable.
- **Output**: Save the downloaded CIF files directly into a **single ZIP archive**.

### Script 3: Biological Assembly Parsing & LMDB Building

- **Goal**: Read the single ZIP file directly, parse each CIF file using `biopython`, and build the `lmdb` database.
- **Logic**:
  1. Iterate over each identified peptide (standard amino acid chain length **4-32 inclusive**).
     - **Cap Removal**: If a chain contains terminal modifications (e.g., ACE), they must be stripped/ignored to evaluate the core standard amino acid sequence length.
  2. Use `scipy.spatial.KDTree` to detect all atoms within a **5 Ångström distance** (configurable) of the peptide. _(Note: Distance calculations must use all 37 atoms of the peptide against all atoms in the assembly)_.
  3. **Rules for inclusion/skipping**:
     - If no other chains are found in the 5Å neighborhood: **SKIP**.
     - If **multiple** unique chains (other than the peptide) are found in the neighborhood: **SKIP**.
     - If the single matching chain is **not an amino acid chain** (non-polypeptide): **SKIP**.
  4. **Data Extraction**:
     - For both peptide and receptor, extract the **full structure coordinates** using the exact 37-atom AlphaFold nomenclature and residue-level quality metadata (**B-factors**, **occupancy**).
     - Coordinates and metadata must be serialized to native binary payloads (`.tobytes()`). Coordinates are `float16`, B-factors are `uint8`, and Occupancy is `uint8` (multiplied by 100).
     - Missing coordinates are represented using `NaN` (`float16`), but for `uint8` arrays (B-factors and occupancy), missing values are explicitly represented with a `255` sentinel value. The entire dataset is then serialized using `msgpack.packb`.
- **Schema**:
  The schema for each LMDB entry must exactly match the following JSON structure (with binary byte strings for arrays):
  ```json
  {
    "pdb_id": "", // PDB identifier
    "entities": [ // one entry per peptide entity
      {
        "entity_id": "", // the entity id
        "sequence": "", // the entity sequence
        "pairs": [ // one entry per peptide chain / receptor pair
          {
            "peptide": {
              "entity_id": "", // the peptide entity id
              "chain": "", // the peptide chain letter
              "sequence": "", // the raw chain sequence of this peptide
              "structure": "<bytes>", // float16 [N, 37, 3] (.tobytes())
              "b_factors": "<bytes>", // uint8 [N, 37] (.tobytes()), 255 = NaN
              "occupancy": "<bytes>"  // uint8 [N, 37] (.tobytes()), occupancy*100, 255 = NaN
            },
            "receptor": {
              "entity_id": "",
              "chain": "",
              "sequence": "",
              "structure": "<bytes>", // float16 [N, 37, 3] (.tobytes())
              "b_factors": "<bytes>", // uint8 [N, 37] (.tobytes()), 255 = NaN
              "occupancy": "<bytes>"  // uint8 [N, 37] (.tobytes()), occupancy*100, 255 = NaN
            }
          }
        ]
      }
    ]
  }
  ```

### Script 4: Filtering and Clustering

- **Goal**: Post-process and refine the constructed LMDB database.
- **Logic**:
  1. Filter the dataset to keep **only** entries where the peptide is associated with a **single receptor chain**.
  2. Analyze the sizes of these single-chain receptors to identify the "sweet spot" optimizing dataset volume versus sequence length.
  3. Filter the LMDB data based on this determined sweet spot.
  4. Cluster the receptors (clustering methodology TBD, pending user instruction).
