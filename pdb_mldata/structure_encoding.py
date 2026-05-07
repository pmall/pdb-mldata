from __future__ import annotations

from typing import Sequence

import gemmi
import numpy as np

ATOM_TYPES_37 = (
    "N",
    "CA",
    "C",
    "CB",
    "O",
    "CG",
    "CG1",
    "CG2",
    "OG",
    "OG1",
    "SG",
    "CD",
    "CD1",
    "CD2",
    "ND1",
    "ND2",
    "OD1",
    "OD2",
    "SD",
    "CE",
    "CE1",
    "CE2",
    "CE3",
    "NE",
    "NE1",
    "NE2",
    "OE1",
    "OE2",
    "CH2",
    "NH1",
    "NH2",
    "OH",
    "CZ",
    "CZ2",
    "CZ3",
    "NZ",
    "OXT",
)
BACKBONE_ATOM_TYPES = ("N", "CA", "C")


def extract_structure(
    residues: Sequence[gemmi.Residue],
) -> tuple[bytes, bytes, bytes]:
    """Extract persisted backbone coordinates, B-factors, and occupancy as bytes."""
    coords = []
    b_factors = []
    occupancy = []

    for residue in residues:
        atom_coords = [[np.nan, np.nan, np.nan] for _ in BACKBONE_ATOM_TYPES]
        atom_b = [np.nan] * len(BACKBONE_ATOM_TYPES)
        atom_occ = [np.nan] * len(BACKBONE_ATOM_TYPES)

        for atom in residue:
            if atom.name not in BACKBONE_ATOM_TYPES:
                continue
            atom_index = BACKBONE_ATOM_TYPES.index(atom.name)
            atom_coords[atom_index] = [atom.pos.x, atom.pos.y, atom.pos.z]

            if atom.b_iso is not None:
                atom_b[atom_index] = atom.b_iso

            if atom.occ is not None:
                atom_occ[atom_index] = atom.occ

        coords.append(atom_coords)
        b_factors.append(atom_b)
        occupancy.append(atom_occ)

    coords_arr = np.asarray(coords, dtype=np.float16)
    b_factors_arr = np.asarray(b_factors, dtype=np.float16)
    occupancy_arr = np.asarray(occupancy, dtype=np.float16)

    return coords_arr.tobytes(), b_factors_arr.tobytes(), occupancy_arr.tobytes()


def extract_binding_validation_structure(
    residues: Sequence[gemmi.Residue],
) -> tuple[np.ndarray, np.ndarray]:
    """Extract the legacy contact tensor used only for build-time validation."""
    coords = []
    b_factors = []

    for residue in residues:
        atom_coords = [[np.nan, np.nan, np.nan] for _ in ATOM_TYPES_37]
        atom_b = [255] * len(ATOM_TYPES_37)

        for atom in residue:
            if atom.name not in ATOM_TYPES_37:
                continue
            atom_index = ATOM_TYPES_37.index(atom.name)
            atom_coords[atom_index] = [atom.pos.x, atom.pos.y, atom.pos.z]

            if atom.b_iso is not None:
                atom_b[atom_index] = (
                    min(max(0, round(atom.b_iso)), 254) if atom.b_iso >= 0 else 0
                )

        coords.append(atom_coords)
        b_factors.append(atom_b)

    coords_arr = np.asarray(coords, dtype=np.float16)
    b_factors_uint8 = np.asarray(b_factors, dtype=np.uint8)
    b_factors_arr = b_factors_uint8.astype(np.float32)
    b_factors_arr[b_factors_uint8 == 255] = np.nan

    return coords_arr, b_factors_arr
