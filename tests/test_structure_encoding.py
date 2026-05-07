from __future__ import annotations

import unittest
from dataclasses import dataclass
from collections.abc import Iterator, Sequence
from typing import cast

import gemmi
import numpy as np

from pdb_mldata.filtering_rules import (
    BACKBONE_ATOM_TYPES,
    ATOM_TYPES_37,
    extract_binding_validation_structure,
    extract_structure,
)
from pdb_mldata.lmdb_utils import ChainData, decode_chain_data


@dataclass(frozen=True)
class FakePosition:
    x: float
    y: float
    z: float


@dataclass(frozen=True)
class FakeAtom:
    name: str
    pos: FakePosition
    b_iso: float | None
    occ: float | None


class FakeResidue:
    def __init__(self, atoms: list[FakeAtom]) -> None:
        self._atoms = atoms

    def __iter__(self) -> Iterator[FakeAtom]:
        return iter(self._atoms)


def fake_residues(atoms: list[FakeAtom]) -> Sequence[gemmi.Residue]:
    return cast(Sequence[gemmi.Residue], [FakeResidue(atoms)])


class StructureEncodingTests(unittest.TestCase):
    def test_extract_structure_stores_backbone_float16_arrays(self) -> None:
        residues = fake_residues(
            [
                FakeAtom("N", FakePosition(1.0, 2.0, 3.0), 11.25, 0.5),
                FakeAtom("CA", FakePosition(4.0, 5.0, 6.0), 22.5, 1.0),
                FakeAtom("O", FakePosition(7.0, 8.0, 9.0), 33.75, 0.25),
            ]
        )

        structure_bytes, b_factor_bytes, occupancy_bytes = extract_structure(residues)

        structure = np.frombuffer(structure_bytes, dtype=np.float16).reshape(-1, 3, 3)
        b_factors = np.frombuffer(b_factor_bytes, dtype=np.float16).reshape(-1, 3)
        occupancy = np.frombuffer(occupancy_bytes, dtype=np.float16).reshape(-1, 3)

        self.assertEqual(BACKBONE_ATOM_TYPES, ("N", "CA", "C"))
        np.testing.assert_array_equal(
            structure[0],
            np.asarray(
                [
                    [1.0, 2.0, 3.0],
                    [4.0, 5.0, 6.0],
                    [np.nan, np.nan, np.nan],
                ],
                dtype=np.float16,
            ),
        )
        np.testing.assert_array_equal(
            b_factors[0],
            np.asarray([11.25, 22.5, np.nan], dtype=np.float16),
        )
        np.testing.assert_array_equal(
            occupancy[0],
            np.asarray([0.5, 1.0, np.nan], dtype=np.float16),
        )

    def test_validation_structure_keeps_legacy_contact_atoms(self) -> None:
        residues = fake_residues(
            [
                FakeAtom("O", FakePosition(7.0, 8.0, 9.0), 33.75, 0.25),
            ]
        )

        structure, b_factors = extract_binding_validation_structure(residues)
        oxygen_index = ATOM_TYPES_37.index("O")

        self.assertEqual(structure.shape, (1, 37, 3))
        np.testing.assert_array_equal(
            structure[0, oxygen_index],
            np.asarray([7.0, 8.0, 9.0], dtype=np.float16),
        )
        self.assertEqual(float(b_factors[0, oxygen_index]), 34.0)

    def test_decode_chain_data_reads_new_backbone_schema(self) -> None:
        chain_data: ChainData = {
            "entity_id": "1",
            "chain": "A",
            "sequence": "ACDE",
            "residue_names": ["ALA", "CYS", "ASP", "GLU"],
            "structure": np.zeros((4, 3, 3), dtype=np.float16).tobytes(),
            "b_factors": np.ones((4, 3), dtype=np.float16).tobytes(),
            "occupancy": np.full((4, 3), 0.5, dtype=np.float16).tobytes(),
        }

        decoded = decode_chain_data(chain_data)

        self.assertEqual(cast(np.ndarray, decoded["structure"]).shape, (4, 3, 3))
        self.assertEqual(cast(np.ndarray, decoded["b_factors"]).dtype, np.float16)
        self.assertEqual(cast(np.ndarray, decoded["occupancy"]).dtype, np.float16)


if __name__ == "__main__":
    unittest.main()
