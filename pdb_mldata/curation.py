from __future__ import annotations

from dataclasses import dataclass
from typing import cast

import numpy as np
from scipy.spatial import KDTree

from pdb_mldata.lmdb_utils import ChainData, PairData

STANDARD_AMINO_ACID_CODES = frozenset("ACDEFGHIKLMNPQRSTVWY")


@dataclass(frozen=True)
class PeptideContentMetrics:
    peptide_residues: int
    standard_residues: int
    nonstandard_residues: int


@dataclass(frozen=True)
class PeptideContentFilter:
    min_standard_residues: int
    max_nonstandard_fraction: float
    standard_codes: frozenset[str] = STANDARD_AMINO_ACID_CODES


@dataclass(frozen=True)
class PeptideContentDecision:
    is_accepted: bool
    reason: str
    metrics: PeptideContentMetrics


@dataclass(frozen=True)
class BindingMetrics:
    peptide_residues: int
    peptide_atoms: int
    receptor_atoms: int
    valid_contact_residues: int
    valid_contact_fraction: float


@dataclass(frozen=True)
class BindingFilter:
    distance: float
    min_contact_residues: int
    min_contact_fraction: float
    max_contact_atom_b_factor: float


@dataclass(frozen=True)
class BindingDecision:
    is_accepted: bool
    reason: str
    metrics: BindingMetrics


def calculate_peptide_content_metrics(
    peptide_sequence: str,
    standard_codes: frozenset[str],
) -> PeptideContentMetrics:
    standard_residues = sum(
        1 for residue in peptide_sequence if residue in standard_codes
    )
    peptide_residues = len(peptide_sequence)
    return PeptideContentMetrics(
        peptide_residues=peptide_residues,
        standard_residues=standard_residues,
        nonstandard_residues=peptide_residues - standard_residues,
    )


def evaluate_peptide_content(
    peptide_sequence: str,
    peptide_content_filter: PeptideContentFilter,
) -> PeptideContentDecision:
    metrics = calculate_peptide_content_metrics(
        peptide_sequence=peptide_sequence,
        standard_codes=peptide_content_filter.standard_codes,
    )
    if metrics.standard_residues < peptide_content_filter.min_standard_residues:
        return PeptideContentDecision(
            is_accepted=False,
            reason="peptide_content_fewer_than_4_standard_aa",
            metrics=metrics,
        )
    if (
        metrics.nonstandard_residues
        > metrics.peptide_residues * peptide_content_filter.max_nonstandard_fraction
    ):
        return PeptideContentDecision(
            is_accepted=False,
            reason="peptide_content_too_many_nonstandard_aa",
            metrics=metrics,
        )
    return PeptideContentDecision(
        is_accepted=True,
        reason="accepted",
        metrics=metrics,
    )


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


def find_valid_contact_residue_indexes(
    peptide_atom_coordinates: np.ndarray,
    peptide_atom_residue_indexes: np.ndarray,
    peptide_atom_b_factors: np.ndarray,
    receptor_atom_coordinates: np.ndarray,
    distance: float,
    max_contact_atom_b_factor: float,
) -> set[int]:
    """Find peptide residues with at least one close, low-B-factor peptide atom."""
    receptor_tree = KDTree(receptor_atom_coordinates)
    nearby_receptor_atoms_by_peptide_atom = receptor_tree.query_ball_point(
        peptide_atom_coordinates,
        r=distance,
    )

    valid_contact_residue_indexes = set()
    for peptide_atom_index, nearby_receptor_atom_indexes in enumerate(
        nearby_receptor_atoms_by_peptide_atom
    ):
        if not nearby_receptor_atom_indexes:
            continue
        if peptide_atom_b_factors[peptide_atom_index] > max_contact_atom_b_factor:
            continue
        valid_contact_residue_indexes.add(
            int(peptide_atom_residue_indexes[peptide_atom_index])
        )

    return valid_contact_residue_indexes


def calculate_binding_metrics(
    peptide: ChainData,
    receptor: ChainData,
    distance: float,
    max_contact_atom_b_factor: float,
) -> BindingMetrics:
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
    receptor_atoms, _receptor_residue_indexes = collect_finite_atom_coordinates(
        receptor_structure
    )
    if peptide_atoms.size == 0 or receptor_atoms.size == 0:
        return BindingMetrics(
            peptide_residues=peptide_residues,
            peptide_atoms=len(peptide_atoms),
            receptor_atoms=len(receptor_atoms),
            valid_contact_residues=0,
            valid_contact_fraction=0.0,
        )

    valid_contact_residue_indexes = find_valid_contact_residue_indexes(
        peptide_atom_coordinates=peptide_atoms,
        peptide_atom_residue_indexes=peptide_residue_indexes,
        peptide_atom_b_factors=peptide_atom_b_factors,
        receptor_atom_coordinates=receptor_atoms,
        distance=distance,
        max_contact_atom_b_factor=max_contact_atom_b_factor,
    )
    valid_contact_residues = len(valid_contact_residue_indexes)
    valid_contact_fraction = valid_contact_residues / peptide_residues

    return BindingMetrics(
        peptide_residues=peptide_residues,
        peptide_atoms=len(peptide_atoms),
        receptor_atoms=len(receptor_atoms),
        valid_contact_residues=valid_contact_residues,
        valid_contact_fraction=valid_contact_fraction,
    )


def evaluate_binding_pair(
    pair: PairData,
    binding_filter: BindingFilter,
) -> BindingDecision:
    metrics = calculate_binding_metrics(
        peptide=pair["peptide"],
        receptor=pair["receptor"],
        distance=binding_filter.distance,
        max_contact_atom_b_factor=binding_filter.max_contact_atom_b_factor,
    )
    if metrics.peptide_atoms == 0 or metrics.receptor_atoms == 0:
        return BindingDecision(
            is_accepted=False,
            reason="no_usable_coordinates",
            metrics=metrics,
        )
    if metrics.valid_contact_residues < binding_filter.min_contact_residues:
        return BindingDecision(
            is_accepted=False,
            reason="too_few_contact_residues",
            metrics=metrics,
        )
    if metrics.valid_contact_fraction < binding_filter.min_contact_fraction:
        return BindingDecision(
            is_accepted=False,
            reason="too_low_contact_fraction",
            metrics=metrics,
        )
    return BindingDecision(
        is_accepted=True,
        reason="accepted",
        metrics=metrics,
    )
