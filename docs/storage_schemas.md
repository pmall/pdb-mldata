# Storage Schemas

This project stores LMDB values with `msgpack.packb(..., use_bin_type=True)`.
Encoding and decoding helpers live in `pdb_mldata/lmdb_utils.py`.

## Binary Array Conventions

- Structure coordinates are stored as `float16` bytes from NumPy `.tobytes()`.
- B-factors are stored as `float16` bytes from NumPy `.tobytes()`.
- Occupancy values are stored as `float16` bytes from NumPy `.tobytes()`.
- Missing coordinates, B-factors, and occupancy values are stored as `NaN`.
- LMDB structure arrays store 3 backbone atom positions per residue.
- The backbone atom order is `N`, `CA`, `C`, as defined by
  `BACKBONE_ATOM_TYPES` in `pdb_mldata/structure_encoding.py`.
- Raw LMDB binding validation uses the 37 AlphaFold atom set defined by
  `ATOM_TYPES_37` in `pdb_mldata/structure_encoding.py`.

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
- `structure`: `float16 [N, 3, 3]`, in `N`, `CA`, `C` atom order.
- `b_factors`: `float16 [N, 3]`, in `N`, `CA`, `C` atom order.
- `occupancy`: `float16 [N, 3]`, in `N`, `CA`, `C` atom order.
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
- `structure`: `float16 [N, 3, 3]`, in `N`, `CA`, `C` atom order.
- `b_factors`: `float16 [N, 3]`, in `N`, `CA`, `C` atom order.
- `occupancy`: `float16 [N, 3]`, in `N`, `CA`, `C` atom order.
- Interface fields have the same meaning as in the raw LMDB schema.

Helpers:

- `encode_best_pair_lmdb_entry`
- `decode_best_pair_lmdb_entry`
- `unpack_best_pair_lmdb_entry_for_export`
