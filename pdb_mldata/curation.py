from __future__ import annotations

from dataclasses import dataclass
from typing import cast

import numpy as np
from scipy.spatial import KDTree

from pdb_mldata.lmdb_utils import ChainData

STANDARD_AMINO_ACID_CODES = frozenset("ACDEFGHIKLMNPQRSTVWY")


@dataclass(frozen=True)
class BestPairMetrics:
    valid_contact_residues: int
    valid_contact_fraction: float
    mean_valid_contact_atom_b_factor: float
    finite_peptide_residues: int
    receptor_residues: int


@dataclass(frozen=True)
class BindingInterfaceMetrics:
    peptide_atoms: int
    receptor_atoms: int
    peptide_contact_residues: int
    peptide_contact_fraction: float
    peptide_interface_start: int
    peptide_interface_end: int
    receptor_contact_residues: int
    receptor_interface_start: int
    receptor_interface_end: int
    mean_contact_atom_b_factor: float


@dataclass(frozen=True)
class BindingFilter:
    distance: float
    min_contact_residues: int
    min_contact_fraction: float
    max_contact_atom_b_factor: float


def collect_finite_atom_coordinates(
    structure: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Return finite atom coordinates and their source residue indexes."""
    finite_atom_mask = np.isfinite(structure).all(axis=2)
    residue_indexes = np.broadcast_to(
        np.arange(structure.shape[0])[:, np.newaxis],
        finite_atom_mask.shape,
    )
    return (
        structure[finite_atom_mask].astype(np.float32),
        residue_indexes[finite_atom_mask],
    )


def collect_finite_peptide_atom_data(
    structure: np.ndarray,
    b_factors: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return finite peptide atom coordinates, residue indexes, and B-factors."""
    finite_atom_mask = np.isfinite(structure).all(axis=2) & np.isfinite(b_factors)
    residue_indexes = np.broadcast_to(
        np.arange(structure.shape[0])[:, np.newaxis],
        finite_atom_mask.shape,
    )
    return (
        structure[finite_atom_mask].astype(np.float32),
        residue_indexes[finite_atom_mask],
        b_factors[finite_atom_mask].astype(np.float32),
    )


def collect_valid_contact_data(
    peptide_atom_coordinates: np.ndarray,
    peptide_atom_residue_indexes: np.ndarray,
    peptide_atom_b_factors: np.ndarray,
    receptor_atom_coordinates: np.ndarray,
    receptor_atom_residue_indexes: np.ndarray,
    distance: float,
    max_contact_atom_b_factor: float,
) -> tuple[set[int], set[int], np.ndarray]:
    """Find valid contact residues and the peptide atom B-factors that support them."""
    receptor_tree = KDTree(receptor_atom_coordinates)
    nearby_receptor_atoms_by_peptide_atom = receptor_tree.query_ball_point(
        peptide_atom_coordinates,
        r=distance,
    )

    valid_contact_residue_indexes = set()
    receptor_contact_residue_indexes = set()
    valid_contact_atom_b_factors: list[float] = []
    for peptide_atom_index, nearby_receptor_atom_indexes in enumerate(
        nearby_receptor_atoms_by_peptide_atom
    ):
        if not nearby_receptor_atom_indexes:
            continue
        peptide_atom_b_factor = peptide_atom_b_factors[peptide_atom_index]
        if peptide_atom_b_factor > max_contact_atom_b_factor:
            continue
        valid_contact_residue_indexes.add(
            int(peptide_atom_residue_indexes[peptide_atom_index])
        )
        receptor_contact_residue_indexes.update(
            int(receptor_atom_residue_indexes[receptor_atom_index])
            for receptor_atom_index in nearby_receptor_atom_indexes
        )
        valid_contact_atom_b_factors.append(float(peptide_atom_b_factor))

    return (
        valid_contact_residue_indexes,
        receptor_contact_residue_indexes,
        np.asarray(valid_contact_atom_b_factors, dtype=np.float32),
    )


def count_finite_residues(structure: np.ndarray) -> int:
    """Count residues with at least one finite atom coordinate."""
    finite_atom_mask = cast(np.ndarray, np.isfinite(structure).all(axis=2))
    return sum(1 for residue_atom_mask in finite_atom_mask if residue_atom_mask.any())


def interface_span(residue_indexes: set[int]) -> tuple[int, int]:
    if not residue_indexes:
        return 0, 0
    return min(residue_indexes) + 1, max(residue_indexes) + 1


def calculate_binding_interface_metrics(
    peptide: ChainData,
    receptor: ChainData,
    distance: float,
    max_contact_atom_b_factor: float,
) -> BindingInterfaceMetrics:
    peptide_structure = cast(np.ndarray, peptide["structure"])
    peptide_b_factors = cast(np.ndarray, peptide["b_factors"])
    receptor_structure = cast(np.ndarray, receptor["structure"])
    peptide_residues = len(peptide["sequence"])

    (
        peptide_atoms,
        peptide_residue_indexes,
        peptide_atom_b_factors,
    ) = collect_finite_peptide_atom_data(
        structure=peptide_structure,
        b_factors=peptide_b_factors,
    )
    receptor_atoms, receptor_residue_indexes = collect_finite_atom_coordinates(
        receptor_structure
    )

    if peptide_atoms.size == 0 or receptor_atoms.size == 0:
        return BindingInterfaceMetrics(
            peptide_atoms=len(peptide_atoms),
            receptor_atoms=len(receptor_atoms),
            peptide_contact_residues=0,
            peptide_contact_fraction=0.0,
            peptide_interface_start=0,
            peptide_interface_end=0,
            receptor_contact_residues=0,
            receptor_interface_start=0,
            receptor_interface_end=0,
            mean_contact_atom_b_factor=float("inf"),
        )

    (
        peptide_contact_residue_indexes,
        receptor_contact_residue_indexes,
        valid_contact_atom_b_factors,
    ) = collect_valid_contact_data(
        peptide_atom_coordinates=peptide_atoms,
        peptide_atom_residue_indexes=peptide_residue_indexes,
        peptide_atom_b_factors=peptide_atom_b_factors,
        receptor_atom_coordinates=receptor_atoms,
        receptor_atom_residue_indexes=receptor_residue_indexes,
        distance=distance,
        max_contact_atom_b_factor=max_contact_atom_b_factor,
    )
    peptide_interface_start, peptide_interface_end = interface_span(
        peptide_contact_residue_indexes
    )
    receptor_interface_start, receptor_interface_end = interface_span(
        receptor_contact_residue_indexes
    )
    mean_contact_atom_b_factor = float("inf")
    if valid_contact_atom_b_factors.size > 0:
        mean_contact_atom_b_factor = float(valid_contact_atom_b_factors.mean())

    peptide_contact_residues = len(peptide_contact_residue_indexes)
    return BindingInterfaceMetrics(
        peptide_atoms=len(peptide_atoms),
        receptor_atoms=len(receptor_atoms),
        peptide_contact_residues=peptide_contact_residues,
        peptide_contact_fraction=peptide_contact_residues / peptide_residues,
        peptide_interface_start=peptide_interface_start,
        peptide_interface_end=peptide_interface_end,
        receptor_contact_residues=len(receptor_contact_residue_indexes),
        receptor_interface_start=receptor_interface_start,
        receptor_interface_end=receptor_interface_end,
        mean_contact_atom_b_factor=mean_contact_atom_b_factor,
    )
