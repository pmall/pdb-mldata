"""Audit why metadata entries are absent from the raw LMDB.

Goal:
- Read PDB IDs from `data/metadata.csv`.
- Read accepted raw-LMDB entry IDs from `data/pdb_mldata.lmdb`.
- Rescan only metadata IDs absent from the raw LMDB.
- Write detailed rejection reports under `data/build_rejection_audit`.

This is a diagnostic script. It does not write or modify LMDB data.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import io
import shutil
import zipfile
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Sequence

import gemmi
import lmdb
import numpy as np
from scipy.spatial import KDTree
from tqdm import tqdm

from pdb_mldata.curation import BindingFilter, calculate_binding_interface_metrics
from pdb_mldata.filtering_rules import (
    DuplicatePairError,
    NeighborRecord,
    PeptideSequenceFilter,
    analyze_neighbors,
    annotate_nonpolymer_neighbor_subtypes,
    build_assembly_atoms,
    build_subchain_entity_map,
    collect_atom_positions,
    extract_structure,
    is_amino_acid_receptor,
    is_valid_peptide_sequence,
    mark_solvent_neighbors,
    normalize_residue_names,
    pair_key,
)
from pdb_mldata.lmdb_utils import PairData, decode_chain_data

DEFAULT_METADATA_CSV = Path("data/metadata.csv")
DEFAULT_ASSEMBLIES_ZIP = Path("data/assemblies.zip")
DEFAULT_LMDB_PATH = Path("data/pdb_mldata.lmdb")
DEFAULT_OUTPUT_DIR = Path("data/build_rejection_audit")
DEFAULT_MIN_LENGTH = 4
DEFAULT_MAX_LENGTH = 32
DEFAULT_MIN_STANDARD_PEPTIDE_RESIDUES = 4
DEFAULT_MAX_NONSTANDARD_PEPTIDE_FRACTION = 0.2
DEFAULT_MIN_RECEPTOR_LENGTH = 50
DEFAULT_DISTANCE = 5.0
DEFAULT_MIN_CONTACT_RESIDUES = 4
DEFAULT_MIN_CONTACT_FRACTION = 0.5
DEFAULT_MAX_CONTACT_ATOM_B_FACTOR = 70.0


@dataclass
class AuditStats:
    metadata_entries: int = 0
    accepted_entries: int = 0
    missing_entries: int = 0
    scanned_entries: int = 0
    entry_reasons: Counter[str] = field(default_factory=Counter)
    chain_reasons: Counter[str] = field(default_factory=Counter)


@dataclass(frozen=True)
class NeighborContactStats:
    atom_contacts: int
    contacted_peptide_residues: int


def read_metadata_ids(metadata_csv: Path) -> list[str]:
    with metadata_csv.open("r", newline="") as handle:
        reader = csv.DictReader(handle)
        return sorted({row["pdb_id"].upper() for row in reader if row.get("pdb_id")})


def read_lmdb_ids(lmdb_path: Path) -> set[str]:
    accepted_ids: set[str] = set()
    env = lmdb.open(str(lmdb_path), readonly=True, lock=False, subdir=True)
    with env.begin() as txn:
        for key, _value in txn.cursor():
            accepted_ids.add(bytes(key).decode().upper())
    env.close()
    return accepted_ids


def build_zip_member_map(zip_path: Path) -> dict[str, str]:
    with zipfile.ZipFile(zip_path, "r") as zf:
        return {
            name.split("-")[0].upper(): name
            for name in zf.namelist()
            if name.endswith(".cif.gz")
        }


def parse_structure(
    zf: zipfile.ZipFile,
    filename: str,
) -> gemmi.Structure:
    with gzip.open(io.BytesIO(zf.read(filename)), "rt") as handle:
        doc = gemmi.cif.read_string(handle.read())
    return gemmi.make_structure_from_block(doc[0])


def residue_names(residues: Sequence[gemmi.Residue]) -> list[str]:
    return [residue.name for residue in residues]


def write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def reset_output_dir(output_dir: Path) -> None:
    if output_dir.exists():
        shutil.rmtree(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)


def count_neighbor_contacts(
    peptide_residues: Sequence[gemmi.Residue],
    neighbor_residues: Sequence[gemmi.Residue],
    distance: float,
) -> NeighborContactStats:
    peptide_atom_coordinates: list[list[float]] = []
    peptide_atom_residue_indexes: list[int] = []
    for residue_index, residue in enumerate(peptide_residues):
        for atom in residue:
            peptide_atom_coordinates.append([atom.pos.x, atom.pos.y, atom.pos.z])
            peptide_atom_residue_indexes.append(residue_index)

    neighbor_atoms = collect_atom_positions(neighbor_residues)
    if not peptide_atom_coordinates or neighbor_atoms.size == 0:
        return NeighborContactStats(atom_contacts=0, contacted_peptide_residues=0)

    neighbor_tree = KDTree(neighbor_atoms)
    contacts_by_peptide_atom = neighbor_tree.query_ball_point(
        np.asarray(peptide_atom_coordinates, dtype=np.float32),
        r=distance,
    )
    contacted_residue_indexes = {
        peptide_atom_residue_indexes[peptide_atom_index]
        for peptide_atom_index, contacts in enumerate(contacts_by_peptide_atom)
        if contacts
    }
    atom_contacts = sum(len(contacts) for contacts in contacts_by_peptide_atom)
    return NeighborContactStats(
        atom_contacts=atom_contacts,
        contacted_peptide_residues=len(contacted_residue_indexes),
    )


def neighbor_receptor_length(
    model: gemmi.Model,
    record: NeighborRecord,
    subchain_to_entity: dict[str, gemmi.Entity],
) -> int | str:
    entity = subchain_to_entity.get(record["subchain_id"])
    if not is_amino_acid_receptor(entity):
        return ""
    neighbor_subchain = model.get_subchain(record["subchain_id"])
    if neighbor_subchain is None:
        return ""
    span = normalize_residue_names(residue_names(list(neighbor_subchain)))
    if span is None:
        return ""
    return len(span.sequence)


def record_neighbor_rows(
    pdb_id: str,
    entity_id: str,
    peptide_chain_id: str,
    peptide_length: int,
    peptide_residues: Sequence[gemmi.Residue],
    model: gemmi.Model,
    subchain_to_entity: dict[str, gemmi.Entity],
    neighbor_records: list[NeighborRecord],
    distance: float,
    neighbor_rows: list[dict[str, object]],
) -> None:
    for record in neighbor_records:
        neighbor_subchain = model.get_subchain(record["subchain_id"])
        neighbor_residues = list(neighbor_subchain) if neighbor_subchain else []
        contact_stats = count_neighbor_contacts(
            peptide_residues=peptide_residues,
            neighbor_residues=neighbor_residues,
            distance=distance,
        )
        neighbor_rows.append(
            {
                "pdb_id": pdb_id,
                "peptide_entity_id": entity_id,
                "peptide_chain_id": peptide_chain_id,
                "peptide_length": peptide_length,
                "neighbor_subchain_id": record["subchain_id"],
                "neighbor_entity_id": record["entity_id"] or "",
                "neighbor_entity_type": record["entity_type"],
                "neighbor_polymer_type": record["polymer_type"],
                "neighbor_nonpolymer_subtype": record.get("nonpolymer_subtype", ""),
                "neighbor_is_solvent": record["is_solvent"],
                "neighbor_counts_as_competing_receptor": record[
                    "counts_as_competing_receptor"
                ],
                "neighbor_receptor_length_after_trim": neighbor_receptor_length(
                    model=model,
                    record=record,
                    subchain_to_entity=subchain_to_entity,
                ),
                "neighbor_atom_contact_count_raw": contact_stats.atom_contacts,
                "neighbor_contacted_peptide_residues_raw": (
                    contact_stats.contacted_peptide_residues
                ),
            }
        )


def chain_row(
    pdb_id: str,
    entity_id: str,
    entity_sequence: str,
    peptide_chain_id: str,
    peptide_sequence: str,
    reason: str,
    neighbor_records: list[NeighborRecord],
    receptor_subchain_id: str,
    receptor_length: int | str,
    peptide_contact_residues: int | str = "",
    peptide_contact_fraction: float | str = "",
    receptor_contact_residues: int | str = "",
    mean_contact_atom_b_factor: float | str = "",
) -> dict[str, object]:
    non_solvent_neighbors = [
        record for record in neighbor_records if not record["is_solvent"]
    ]
    competing_neighbors = [
        record for record in neighbor_records if record["counts_as_competing_receptor"]
    ]
    return {
        "pdb_id": pdb_id,
        "peptide_entity_id": entity_id,
        "entity_length": len(entity_sequence),
        "entity_sequence": entity_sequence,
        "peptide_chain_id": peptide_chain_id,
        "peptide_length": len(peptide_sequence),
        "peptide_sequence": peptide_sequence,
        "reason": reason,
        "neighbor_count": len(neighbor_records),
        "non_solvent_neighbor_count": len(non_solvent_neighbors),
        "competing_protein_neighbor_count": len(competing_neighbors),
        "receptor_subchain_id": receptor_subchain_id,
        "receptor_length": receptor_length,
        "peptide_contact_residues": peptide_contact_residues,
        "peptide_contact_fraction": peptide_contact_fraction,
        "receptor_contact_residues": receptor_contact_residues,
        "mean_contact_atom_b_factor": mean_contact_atom_b_factor,
    }


def audit_assembly(
    pdb_id: str,
    zf: zipfile.ZipFile,
    filename: str,
    peptide_filter: PeptideSequenceFilter,
    min_receptor_len: int,
    distance: float,
    binding_filter: BindingFilter,
    entry_rows: list[dict[str, object]],
    chain_rows: list[dict[str, object]],
    neighbor_rows: list[dict[str, object]],
    stats: AuditStats,
) -> None:
    try:
        structure = parse_structure(zf=zf, filename=filename)
    except Exception as exc:
        stats.entry_reasons["parse_error"] += 1
        entry_rows.append({"pdb_id": pdb_id, "reason": "parse_error", "detail": exc})
        return

    if len(structure) != 1:
        stats.entry_reasons["multi_model_entry"] += 1
        entry_rows.append(
            {"pdb_id": pdb_id, "reason": "multi_model_entry", "detail": len(structure)}
        )
        return

    model = structure[0]
    subchain_to_entity, subchain_to_entity_id = build_subchain_entity_map(structure)
    assembly_coords, atom_subchains = build_assembly_atoms(model)
    if assembly_coords.size == 0:
        stats.entry_reasons["no_assembly_atoms"] += 1
        entry_rows.append(
            {"pdb_id": pdb_id, "reason": "no_assembly_atoms", "detail": ""}
        )
        return

    tree = KDTree(assembly_coords)
    pair_keys: set[tuple[str, str, str, str, str]] = set()
    peptide_entities = 0
    peptide_chains = 0
    writable_pairs = 0
    duplicate_pair_error: DuplicatePairError | None = None

    for entity in structure.entities:
        if entity.entity_type != gemmi.EntityType.Polymer:
            continue
        if entity.polymer_type not in (
            gemmi.PolymerType.PeptideL,
            gemmi.PolymerType.PeptideD,
        ):
            continue

        entity_span = normalize_residue_names(entity.full_sequence)
        if entity_span is None:
            continue
        if not is_valid_peptide_sequence(entity_span.sequence, peptide_filter):
            continue

        peptide_entities += 1
        for subchain_id in entity.subchains:
            peptide_chains += 1
            subchain = model.get_subchain(subchain_id)
            if subchain is None:
                stats.chain_reasons["subchain_missing"] += 1
                chain_rows.append(
                    chain_row(
                        pdb_id=pdb_id,
                        entity_id=entity.name,
                        entity_sequence=entity_span.sequence,
                        peptide_chain_id=subchain_id,
                        peptide_sequence="",
                        reason="subchain_missing",
                        neighbor_records=[],
                        receptor_subchain_id="",
                        receptor_length="",
                    )
                )
                continue

            residues = list(subchain)
            peptide_span = normalize_residue_names(residue_names(residues))
            if peptide_span is None:
                stats.chain_reasons["observed_chain_normalization_failed"] += 1
                chain_rows.append(
                    chain_row(
                        pdb_id=pdb_id,
                        entity_id=entity.name,
                        entity_sequence=entity_span.sequence,
                        peptide_chain_id=subchain_id,
                        peptide_sequence="",
                        reason="observed_chain_normalization_failed",
                        neighbor_records=[],
                        receptor_subchain_id="",
                        receptor_length="",
                    )
                )
                continue
            if not is_valid_peptide_sequence(peptide_span.sequence, peptide_filter):
                stats.chain_reasons["observed_chain_peptide_filter_failed"] += 1
                chain_rows.append(
                    chain_row(
                        pdb_id=pdb_id,
                        entity_id=entity.name,
                        entity_sequence=entity_span.sequence,
                        peptide_chain_id=subchain_id,
                        peptide_sequence=peptide_span.sequence,
                        reason="observed_chain_peptide_filter_failed",
                        neighbor_records=[],
                        receptor_subchain_id="",
                        receptor_length="",
                    )
                )
                continue

            peptide_residues = list(
                residues[peptide_span.trim_start : peptide_span.trim_end + 1]
            )
            peptide_atoms = collect_atom_positions(peptide_residues)
            if peptide_atoms.size == 0:
                stats.chain_reasons["no_peptide_atoms"] += 1
                chain_rows.append(
                    chain_row(
                        pdb_id=pdb_id,
                        entity_id=entity.name,
                        entity_sequence=entity_span.sequence,
                        peptide_chain_id=subchain_id,
                        peptide_sequence=peptide_span.sequence,
                        reason="no_peptide_atoms",
                        neighbor_records=[],
                        receptor_subchain_id="",
                        receptor_length="",
                    )
                )
                continue

            neighbor_records, _competing_neighbors = analyze_neighbors(
                peptide_atoms=peptide_atoms,
                tree=tree,
                atom_subchains=atom_subchains,
                peptide_subchain_id=subchain_id,
                subchain_to_entity=subchain_to_entity,
                distance=distance,
            )
            annotate_nonpolymer_neighbor_subtypes(model, neighbor_records)
            mark_solvent_neighbors(neighbor_records)
            record_neighbor_rows(
                pdb_id=pdb_id,
                entity_id=entity.name,
                peptide_chain_id=subchain_id,
                peptide_length=len(peptide_span.sequence),
                peptide_residues=peptide_residues,
                model=model,
                subchain_to_entity=subchain_to_entity,
                neighbor_records=neighbor_records,
                distance=distance,
                neighbor_rows=neighbor_rows,
            )

            competing_neighbors = [
                record["subchain_id"]
                for record in neighbor_records
                if record["counts_as_competing_receptor"]
            ]
            non_solvent_neighbors = [
                record for record in neighbor_records if not record["is_solvent"]
            ]
            if not neighbor_records:
                reason = "no_neighbors_within_distance"
            elif not non_solvent_neighbors:
                reason = "only_solvent_or_ignored_neighbors"
            elif len(non_solvent_neighbors) > 1:
                reason = "multiple_non_solvent_neighbors"
            elif not competing_neighbors:
                reason = "no_competing_protein_neighbor"
            else:
                reason = "candidate_receptor_selected"

            if reason != "candidate_receptor_selected":
                stats.chain_reasons[reason] += 1
                chain_rows.append(
                    chain_row(
                        pdb_id=pdb_id,
                        entity_id=entity.name,
                        entity_sequence=entity_span.sequence,
                        peptide_chain_id=subchain_id,
                        peptide_sequence=peptide_span.sequence,
                        reason=reason,
                        neighbor_records=neighbor_records,
                        receptor_subchain_id="",
                        receptor_length="",
                    )
                )
                continue

            receptor_subchain_id = competing_neighbors[0]
            receptor_chain = model.get_subchain(receptor_subchain_id)
            if receptor_chain is None:
                stats.chain_reasons["receptor_chain_missing"] += 1
                chain_rows.append(
                    chain_row(
                        pdb_id=pdb_id,
                        entity_id=entity.name,
                        entity_sequence=entity_span.sequence,
                        peptide_chain_id=subchain_id,
                        peptide_sequence=peptide_span.sequence,
                        reason="receptor_chain_missing",
                        neighbor_records=neighbor_records,
                        receptor_subchain_id=receptor_subchain_id,
                        receptor_length="",
                    )
                )
                continue

            receptor_entity_id = subchain_to_entity_id.get(receptor_subchain_id, "")
            if not receptor_entity_id:
                stats.chain_reasons["receptor_entity_missing"] += 1
                chain_rows.append(
                    chain_row(
                        pdb_id=pdb_id,
                        entity_id=entity.name,
                        entity_sequence=entity_span.sequence,
                        peptide_chain_id=subchain_id,
                        peptide_sequence=peptide_span.sequence,
                        reason="receptor_entity_missing",
                        neighbor_records=neighbor_records,
                        receptor_subchain_id=receptor_subchain_id,
                        receptor_length="",
                    )
                )
                continue

            receptor_residues = list(receptor_chain)
            receptor_span = normalize_residue_names(residue_names(receptor_residues))
            if receptor_span is None:
                stats.chain_reasons["receptor_normalization_failed"] += 1
                chain_rows.append(
                    chain_row(
                        pdb_id=pdb_id,
                        entity_id=entity.name,
                        entity_sequence=entity_span.sequence,
                        peptide_chain_id=subchain_id,
                        peptide_sequence=peptide_span.sequence,
                        reason="receptor_normalization_failed",
                        neighbor_records=neighbor_records,
                        receptor_subchain_id=receptor_subchain_id,
                        receptor_length="",
                    )
                )
                continue
            receptor_length = len(receptor_span.sequence)
            if receptor_length < min_receptor_len:
                stats.chain_reasons["receptor_under_min_length"] += 1
                chain_rows.append(
                    chain_row(
                        pdb_id=pdb_id,
                        entity_id=entity.name,
                        entity_sequence=entity_span.sequence,
                        peptide_chain_id=subchain_id,
                        peptide_sequence=peptide_span.sequence,
                        reason="receptor_under_min_length",
                        neighbor_records=neighbor_records,
                        receptor_subchain_id=receptor_subchain_id,
                        receptor_length=receptor_length,
                    )
                )
                continue

            trimmed_receptor_residues = list(
                receptor_residues[receptor_span.trim_start : receptor_span.trim_end + 1]
            )
            peptide_structure, peptide_b, peptide_occ = extract_structure(
                peptide_residues
            )
            receptor_structure, receptor_b, receptor_occ = extract_structure(
                trimmed_receptor_residues
            )
            pair: PairData = {
                "peptide": {
                    "entity_id": entity.name,
                    "chain": subchain_id,
                    "sequence": peptide_span.sequence,
                    "residue_names": peptide_span.residue_names,
                    "structure": peptide_structure,
                    "b_factors": peptide_b,
                    "occupancy": peptide_occ,
                },
                "receptor": {
                    "entity_id": receptor_entity_id,
                    "chain": receptor_subchain_id,
                    "sequence": receptor_span.sequence,
                    "residue_names": receptor_span.residue_names,
                    "structure": receptor_structure,
                    "b_factors": receptor_b,
                    "occupancy": receptor_occ,
                },
            }
            interface_metrics = calculate_binding_interface_metrics(
                peptide=decode_chain_data(pair["peptide"]),
                receptor=decode_chain_data(pair["receptor"]),
                distance=binding_filter.distance,
                max_contact_atom_b_factor=binding_filter.max_contact_atom_b_factor,
            )
            if (
                interface_metrics.peptide_atoms == 0
                or interface_metrics.receptor_atoms == 0
            ):
                reason = "no_usable_coordinates"
            elif (
                interface_metrics.peptide_contact_residues
                < binding_filter.min_contact_residues
            ):
                reason = "too_few_contact_residues"
            elif (
                interface_metrics.peptide_contact_fraction
                < binding_filter.min_contact_fraction
            ):
                reason = "too_low_contact_fraction"
            else:
                reason = "pair_would_be_written"

            if reason != "pair_would_be_written":
                stats.chain_reasons[reason] += 1
                chain_rows.append(
                    chain_row(
                        pdb_id=pdb_id,
                        entity_id=entity.name,
                        entity_sequence=entity_span.sequence,
                        peptide_chain_id=subchain_id,
                        peptide_sequence=peptide_span.sequence,
                        reason=reason,
                        neighbor_records=neighbor_records,
                        receptor_subchain_id=receptor_subchain_id,
                        receptor_length=receptor_length,
                        peptide_contact_residues=(
                            interface_metrics.peptide_contact_residues
                        ),
                        peptide_contact_fraction=(
                            interface_metrics.peptide_contact_fraction
                        ),
                        receptor_contact_residues=(
                            interface_metrics.receptor_contact_residues
                        ),
                        mean_contact_atom_b_factor=(
                            interface_metrics.mean_contact_atom_b_factor
                        ),
                    )
                )
                continue

            current_pair_key = pair_key(pdb_id, pair)
            if current_pair_key in pair_keys:
                duplicate_pair_error = DuplicatePairError(*current_pair_key)
                break
            pair_keys.add(current_pair_key)
            writable_pairs += 1
            stats.chain_reasons["pair_would_be_written"] += 1
            chain_rows.append(
                chain_row(
                    pdb_id=pdb_id,
                    entity_id=entity.name,
                    entity_sequence=entity_span.sequence,
                    peptide_chain_id=subchain_id,
                    peptide_sequence=peptide_span.sequence,
                    reason="pair_would_be_written",
                    neighbor_records=neighbor_records,
                    receptor_subchain_id=receptor_subchain_id,
                    receptor_length=receptor_length,
                    peptide_contact_residues=(
                        interface_metrics.peptide_contact_residues
                    ),
                    peptide_contact_fraction=interface_metrics.peptide_contact_fraction,
                    receptor_contact_residues=(
                        interface_metrics.receptor_contact_residues
                    ),
                    mean_contact_atom_b_factor=(
                        interface_metrics.mean_contact_atom_b_factor
                    ),
                )
            )
        if duplicate_pair_error is not None:
            break

    if duplicate_pair_error is not None:
        stats.entry_reasons["duplicate_pair_key_entry_skip"] += 1
        entry_rows.append(
            {
                "pdb_id": pdb_id,
                "reason": "duplicate_pair_key_entry_skip",
                "detail": duplicate_pair_error,
            }
        )
        return
    if writable_pairs > 0:
        stats.entry_reasons["unexpected_writable_pairs"] += 1
        entry_rows.append(
            {
                "pdb_id": pdb_id,
                "reason": "unexpected_writable_pairs",
                "detail": writable_pairs,
            }
        )
        return
    if peptide_entities == 0:
        reason = "no_peptide_entities_after_entity_filter"
    elif peptide_chains == 0:
        reason = "no_observed_peptide_chains"
    else:
        reason = "no_pairs_after_chain_filtering"
    stats.entry_reasons[reason] += 1
    entry_rows.append({"pdb_id": pdb_id, "reason": reason, "detail": ""})


def write_summary(output_dir: Path, stats: AuditStats) -> None:
    lines = [
        "# Build Rejection Audit",
        "",
        f"- Metadata entries: `{stats.metadata_entries}`",
        f"- Accepted raw LMDB entries: `{stats.accepted_entries}`",
        f"- Metadata entries absent from raw LMDB: `{stats.missing_entries}`",
        f"- Scanned absent entries: `{stats.scanned_entries}`",
        "",
        "## Entry Reasons",
        "",
        "| Count | Reason |",
        "| ---: | --- |",
    ]
    for reason, count in stats.entry_reasons.most_common():
        lines.append(f"| {count} | `{reason}` |")
    lines.extend(
        ["", "## Chain Candidate Reasons", "", "| Count | Reason |", "| ---: | --- |"]
    )
    for reason, count in stats.chain_reasons.most_common():
        lines.append(f"| {count} | `{reason}` |")
    lines.append("")
    output_dir.joinpath("summary.md").write_text("\n".join(lines))


def audit_build_rejections(
    metadata_csv: Path,
    assemblies_zip: Path,
    lmdb_path: Path,
    output_dir: Path,
    peptide_filter: PeptideSequenceFilter,
    min_receptor_length: int,
    distance: float,
    binding_filter: BindingFilter,
    limit: int | None,
) -> AuditStats:
    metadata_ids = read_metadata_ids(metadata_csv)
    accepted_ids = read_lmdb_ids(lmdb_path)
    missing_ids = sorted(set(metadata_ids) - accepted_ids)
    if limit is not None:
        missing_ids = missing_ids[:limit]

    stats = AuditStats(
        metadata_entries=len(metadata_ids),
        accepted_entries=len(accepted_ids),
        missing_entries=len(set(metadata_ids) - accepted_ids),
        scanned_entries=len(missing_ids),
    )
    reset_output_dir(output_dir)
    zip_members = build_zip_member_map(assemblies_zip)

    entry_rows: list[dict[str, object]] = []
    chain_rows: list[dict[str, object]] = []
    neighbor_rows: list[dict[str, object]] = []

    with zipfile.ZipFile(assemblies_zip, "r") as zf:
        for pdb_id in tqdm(missing_ids, desc="Auditing build rejections"):
            filename = zip_members.get(pdb_id)
            if filename is None:
                stats.entry_reasons["assembly_not_in_zip"] += 1
                entry_rows.append(
                    {"pdb_id": pdb_id, "reason": "assembly_not_in_zip", "detail": ""}
                )
                continue
            audit_assembly(
                pdb_id=pdb_id,
                zf=zf,
                filename=filename,
                peptide_filter=peptide_filter,
                min_receptor_len=min_receptor_length,
                distance=distance,
                binding_filter=binding_filter,
                entry_rows=entry_rows,
                chain_rows=chain_rows,
                neighbor_rows=neighbor_rows,
                stats=stats,
            )

    write_csv(
        output_dir / "entry_rejections.csv",
        ["pdb_id", "reason", "detail"],
        entry_rows,
    )
    write_csv(
        output_dir / "chain_candidate_rejections.csv",
        [
            "pdb_id",
            "peptide_entity_id",
            "entity_length",
            "entity_sequence",
            "peptide_chain_id",
            "peptide_length",
            "peptide_sequence",
            "reason",
            "neighbor_count",
            "non_solvent_neighbor_count",
            "competing_protein_neighbor_count",
            "receptor_subchain_id",
            "receptor_length",
            "peptide_contact_residues",
            "peptide_contact_fraction",
            "receptor_contact_residues",
            "mean_contact_atom_b_factor",
        ],
        chain_rows,
    )
    write_csv(
        output_dir / "neighbor_rejections.csv",
        [
            "pdb_id",
            "peptide_entity_id",
            "peptide_chain_id",
            "peptide_length",
            "neighbor_subchain_id",
            "neighbor_entity_id",
            "neighbor_entity_type",
            "neighbor_polymer_type",
            "neighbor_nonpolymer_subtype",
            "neighbor_is_solvent",
            "neighbor_counts_as_competing_receptor",
            "neighbor_receptor_length_after_trim",
            "neighbor_atom_contact_count_raw",
            "neighbor_contacted_peptide_residues_raw",
        ],
        neighbor_rows,
    )
    write_summary(output_dir=output_dir, stats=stats)
    return stats


def validate_parameters(
    metadata_csv: Path,
    assemblies_zip: Path,
    lmdb_path: Path,
    output_dir: Path,
    min_length: int,
    max_length: int,
    min_standard_peptide_residues: int,
    max_nonstandard_peptide_fraction: float,
    min_receptor_length: int,
    distance: float,
    min_contact_residues: int,
    min_contact_fraction: float,
    max_contact_atom_b_factor: float,
    limit: int | None,
) -> None:
    if not metadata_csv.exists():
        raise ValueError(f"{metadata_csv} not found")
    if metadata_csv.is_dir():
        raise ValueError("metadata_csv must be a file path, not a directory")
    if not assemblies_zip.exists():
        raise ValueError(f"{assemblies_zip} not found")
    if assemblies_zip.is_dir():
        raise ValueError("assemblies_zip must be a file path, not a directory")
    if not lmdb_path.exists():
        raise ValueError(f"{lmdb_path} not found")
    if lmdb_path.is_file():
        raise ValueError("lmdb_path must be an LMDB directory path")
    if output_dir.exists() and output_dir.is_file():
        raise ValueError("output_dir must be a directory path, not a file")
    if min_length < 1:
        raise ValueError("--min-length must be at least 1")
    if max_length < min_length:
        raise ValueError("--max-length must be greater than or equal to --min-length")
    if min_standard_peptide_residues < 1:
        raise ValueError("--min-standard-peptide-residues must be at least 1")
    if min_standard_peptide_residues > max_length:
        raise ValueError(
            "--min-standard-peptide-residues must be less than or equal to --max-length"
        )
    if not (0 <= max_nonstandard_peptide_fraction <= 1):
        raise ValueError(
            "--max-nonstandard-peptide-fraction must be at least 0 and at most 1"
        )
    if min_receptor_length < 1:
        raise ValueError("--min-receptor-length must be at least 1")
    if max_length > min_receptor_length:
        raise ValueError(
            "--max-length must be less than or equal to --min-receptor-length"
        )
    if distance <= 0:
        raise ValueError("--distance must be greater than 0")
    if min_contact_residues < 1:
        raise ValueError("--min-contact-residues must be at least 1")
    if not (0 <= min_contact_fraction <= 1):
        raise ValueError("--min-contact-fraction must be at least 0 and at most 1")
    if max_contact_atom_b_factor < 0:
        raise ValueError("--max-contact-atom-b-factor must be at least 0")
    if limit is not None and limit < 1:
        raise ValueError("--limit must be at least 1 when provided")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Audit why metadata entries are absent from the raw LMDB."
    )
    parser.add_argument(
        "metadata_csv",
        nargs="?",
        type=Path,
        default=DEFAULT_METADATA_CSV,
        help="Path to metadata CSV (default: data/metadata.csv)",
    )
    parser.add_argument(
        "assemblies_zip",
        nargs="?",
        type=Path,
        default=DEFAULT_ASSEMBLIES_ZIP,
        help="Path to assemblies ZIP (default: data/assemblies.zip)",
    )
    parser.add_argument(
        "lmdb_path",
        nargs="?",
        type=Path,
        default=DEFAULT_LMDB_PATH,
        help="Path to raw LMDB directory (default: data/pdb_mldata.lmdb)",
    )
    parser.add_argument(
        "output_dir",
        nargs="?",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help="Directory for audit reports (default: data/build_rejection_audit)",
    )
    parser.add_argument(
        "--min-length",
        type=int,
        default=DEFAULT_MIN_LENGTH,
        help="Minimum peptide length (default: 4)",
    )
    parser.add_argument(
        "--max-length",
        type=int,
        default=DEFAULT_MAX_LENGTH,
        help="Maximum peptide length (default: 32)",
    )
    parser.add_argument(
        "--min-standard-peptide-residues",
        type=int,
        default=DEFAULT_MIN_STANDARD_PEPTIDE_RESIDUES,
        help="Minimum standard amino acids in peptide sequences (default: 4)",
    )
    parser.add_argument(
        "--max-nonstandard-peptide-fraction",
        type=float,
        default=DEFAULT_MAX_NONSTANDARD_PEPTIDE_FRACTION,
        help="Maximum non-standard peptide residue fraction (default: 0.2)",
    )
    parser.add_argument(
        "--min-receptor-length",
        type=int,
        default=DEFAULT_MIN_RECEPTOR_LENGTH,
        help="Minimum receptor length after cap trimming (default: 50)",
    )
    parser.add_argument(
        "--distance",
        type=float,
        default=DEFAULT_DISTANCE,
        help="Neighbor distance in Angstroms (default: 5.0)",
    )
    parser.add_argument(
        "--min-contact-residues",
        type=int,
        default=DEFAULT_MIN_CONTACT_RESIDUES,
        help="Minimum contacting peptide residues (default: 4)",
    )
    parser.add_argument(
        "--min-contact-fraction",
        type=float,
        default=DEFAULT_MIN_CONTACT_FRACTION,
        help="Minimum contacting peptide residue fraction (default: 0.5)",
    )
    parser.add_argument(
        "--max-contact-atom-b-factor",
        type=float,
        default=DEFAULT_MAX_CONTACT_ATOM_B_FACTOR,
        help="Maximum B-factor for a peptide atom to count as a contact (default: 70.0)",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Scan only the first N missing metadata IDs for smoke verification",
    )
    args = parser.parse_args()

    try:
        validate_parameters(
            metadata_csv=args.metadata_csv,
            assemblies_zip=args.assemblies_zip,
            lmdb_path=args.lmdb_path,
            output_dir=args.output_dir,
            min_length=args.min_length,
            max_length=args.max_length,
            min_standard_peptide_residues=args.min_standard_peptide_residues,
            max_nonstandard_peptide_fraction=args.max_nonstandard_peptide_fraction,
            min_receptor_length=args.min_receptor_length,
            distance=args.distance,
            min_contact_residues=args.min_contact_residues,
            min_contact_fraction=args.min_contact_fraction,
            max_contact_atom_b_factor=args.max_contact_atom_b_factor,
            limit=args.limit,
        )
    except ValueError as exc:
        parser.error(str(exc))

    stats = audit_build_rejections(
        metadata_csv=args.metadata_csv,
        assemblies_zip=args.assemblies_zip,
        lmdb_path=args.lmdb_path,
        output_dir=args.output_dir,
        peptide_filter=PeptideSequenceFilter(
            min_length=args.min_length,
            max_length=args.max_length,
            min_standard_residues=args.min_standard_peptide_residues,
            max_nonstandard_fraction=args.max_nonstandard_peptide_fraction,
        ),
        min_receptor_length=args.min_receptor_length,
        distance=args.distance,
        binding_filter=BindingFilter(
            distance=args.distance,
            min_contact_residues=args.min_contact_residues,
            min_contact_fraction=args.min_contact_fraction,
            max_contact_atom_b_factor=args.max_contact_atom_b_factor,
        ),
        limit=args.limit,
    )
    print(f"Metadata entries: {stats.metadata_entries}")
    print(f"Accepted raw LMDB entries: {stats.accepted_entries}")
    print(f"Metadata entries absent from raw LMDB: {stats.missing_entries}")
    print(f"Scanned absent entries: {stats.scanned_entries}")
    print(f"Wrote audit reports to {args.output_dir}")


if __name__ == "__main__":
    main()
