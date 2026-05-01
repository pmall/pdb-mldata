# Storage Schemas

This project stores LMDB values with `msgpack.packb(..., use_bin_type=True)`.
Encoding and decoding helpers live in `pdb_mldata/lmdb_utils.py`.

## Binary Array Conventions

- Structure coordinates are stored as `float16` bytes from NumPy `.tobytes()`.
- B-factors are stored as `uint8` bytes.
- Occupancy values are multiplied by 100 and stored as `uint8` bytes.
- The `uint8` value `255` is the missing-value sentinel for B-factor and occupancy arrays.
- Raw LMDB stores 37 AlphaFold atom positions per residue.
- The 37-atom order is defined by `ATOM_TYPES_37` in
  `pdb_mldata/filtering_rules.py`.

## Raw LMDB Schema

Artifacts:

- `data/pdb_mldata.lmdb`

Schema:

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
            "occupancy": "<bytes>",
            "interface_start": 1,
            "interface_end": 1,
            "contact_residues": 0,
            "contact_fraction": 0.0,
            "mean_contact_atom_b_factor": 0.0
          },
          "receptor": {
            "entity_id": "",
            "chain": "",
            "sequence": "",
            "residue_names": [],
            "structure": "<bytes>",
            "b_factors": "<bytes>",
            "occupancy": "<bytes>",
            "interface_start": 1,
            "interface_end": 1,
            "contact_residues": 0
          }
        }
      ]
    }
  ]
}
```

Schema notes:

- `pdb_id`: PDB identifier.
- `entities`: one entry per peptide entity.
- Entity-level `sequence` and `residue_names`: ideal peptide entity sequence after cap trimming.
- `pairs`: observed peptide-chain/receptor-chain pairs under that peptide entity.
- Chain-level `sequence` and `residue_names`: observed chain sequence after cap trimming.
- `structure`: `float16 [N, 37, 3]`.
- `b_factors`: `uint8 [N, 37]`, with `255` as missing.
- `occupancy`: `uint8 [N, 37]`, with `255` as missing.
- `interface_start` and `interface_end`: 1-based inclusive residue positions
  relative to the trimmed observed chain sequence.
- `contact_residues`: number of residues in that chain contacting the other
  chain under the binding-contact rule.
- Peptide-only `contact_fraction`: `contact_residues / peptide length`.
- Peptide-only `mean_contact_atom_b_factor`: mean B-factor among valid peptide
  contact atoms.

Helpers:

- `encode_lmdb_entry`
- `decode_lmdb_entry`
- `unpack_lmdb_entry_for_export`

## Best-Pair LMDB Schema

Artifact:

- `data/pdb_mldata_best_pair.lmdb`

Schema:

```json
{
  "pdb_id": "",
  "pairs": [
    {
      "entity": {
        "entity_id": "",
        "sequence": "",
        "residue_names": []
      },
      "peptide": {
        "entity_id": "",
        "chain": "",
        "sequence": "",
        "residue_names": [],
        "structure": "<bytes>",
        "b_factors": "<bytes>",
        "occupancy": "<bytes>",
        "interface_start": 1,
        "interface_end": 1,
        "contact_residues": 0,
        "contact_fraction": 0.0,
        "mean_contact_atom_b_factor": 0.0
      },
      "receptor": {
        "entity_id": "",
        "chain": "",
        "sequence": "",
        "residue_names": [],
        "structure": "<bytes>",
        "b_factors": "<bytes>",
        "occupancy": "<bytes>",
        "interface_start": 1,
        "interface_end": 1,
        "contact_residues": 0
      }
    }
  ]
}
```

Schema notes:

- `pairs`: one selected peptide-chain/receptor-chain pair per peptide entity.
- `entity`: ideal peptide entity sequence and residue names after cap trimming.
- `peptide`: selected observed peptide chain.
- `receptor`: selected observed receptor chain.
- `structure`: `float16 [N, 37, 3]`.
- `b_factors`: `uint8 [N, 37]`, with `255` as missing.
- `occupancy`: `uint8 [N, 37]`, with `255` as missing.
- Interface fields have the same meaning as in the raw LMDB schema.

Helpers:

- `encode_best_pair_lmdb_entry`
- `decode_best_pair_lmdb_entry`
- `unpack_best_pair_lmdb_entry_for_export`
