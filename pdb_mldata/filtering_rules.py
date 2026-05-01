"""Shared filtering rules for peptide/receptor pair admission."""

from __future__ import annotations

import gzip
import io
import zipfile
from dataclasses import dataclass
from typing import NotRequired, Sequence, TypeAlias, TypedDict

import gemmi
import numpy as np
from scipy.spatial import KDTree

from pdb_mldata.curation import BindingFilter, calculate_binding_interface_metrics
from pdb_mldata.lmdb_utils import EntityData, LmdbEntry, PairData, decode_chain_data


class NeighborRecord(TypedDict):
    subchain_id: str
    entity_id: str | None
    entity_type: str
    polymer_type: str
    is_solvent: bool
    counts_as_competing_receptor: bool
    nonpolymer_subtype: NotRequired[str | list[str] | None]


ProcessResult: TypeAlias = tuple[str, LmdbEntry] | None

# These are terminal blocking groups, not modified amino acids; internal matches
# are retained because only terminal caps should affect sequence offsets.
COMMON_TERMINAL_CAPS = frozenset(
    {
        "ACE",
        "ACY",
        "FOR",
        "NH2",
        "NME",
        "NMC",
    }
)
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
PEPTIDE_POLYMER_TYPES = (gemmi.PolymerType.PeptideL, gemmi.PolymerType.PeptideD)
STANDARD_AMINO_ACID_CODES = frozenset("ACDEFGHIKLMNPQRSTVWY")


class DuplicatePairError(Exception):
    def __init__(
        self,
        pdb_id: str,
        peptide_entity_id: str,
        peptide_chain_id: str,
        receptor_entity_id: str,
        receptor_chain_id: str,
    ) -> None:
        self.pdb_id = pdb_id
        self.peptide_entity_id = peptide_entity_id
        self.peptide_chain_id = peptide_chain_id
        self.receptor_entity_id = receptor_entity_id
        self.receptor_chain_id = receptor_chain_id
        super().__init__(
            "duplicate peptide/receptor pair "
            f"{pdb_id} "
            f"{peptide_entity_id}:{peptide_chain_id}->"
            f"{receptor_entity_id}:{receptor_chain_id}"
        )


@dataclass(frozen=True)
class NormalizedResidueSpan:
    sequence: str
    residue_names: list[str]
    trim_start: int
    trim_end: int


@dataclass(frozen=True)
class PeptideSequenceFilter:
    min_length: int
    max_length: int
    min_standard_residues: int
    max_nonstandard_fraction: float
    standard_codes: frozenset[str] = STANDARD_AMINO_ACID_CODES


# Non-polymer subtypes that are crystallographic solvent/buffer artifacts with no
# biological role in the binding interaction (solvents, buffers, PEG, preservatives,
# reducing agents, detergents). Treated as solvent: excluded from the neighbor count.
NONPOLYMER_SOLVENT_SUBTYPES = frozenset(
    {
        # Organic solvents & cryoprotectants
        "GOL",
        "EDO",
        "DMS",
        "MPD",
        "IPA",
        "EOH",
        "MOH",
        "ACN",
        "DCE",
        "DIO",
        "BU3",
        "TBU",
        "OCT",
        "BNZ",
        "ERY",
        "CCN",
        "ETA",
        "HEX",
        "OXY",
        # Buffer inorganic anions
        "SO4",
        "PO4",
        "NO3",
        "CO3",
        "ACT",
        "ACY",
        "FMT",
        "NH4",
        "NH3",
        "CO2",
        "AZI",
        "CIT",
        "TFA",
        # Organic buffer molecules
        "MES",
        "EPE",
        "TRS",
        "IMD",
        "URE",
        "B3P",
        "BTB",
        # PEG fragments & polyethers
        "PG4",
        "PEG",
        "PGV",
        "PGE",
        "1PE",
        "2PE",
        "3PE",
        "PE4",
        "PE5",
        "PE8",
        "P6G",
        "C8E",
        "PG0",
        "PGO",
        "15P",
        "16P",
        "1PG",
        "PEV",
        "PEE",
        "PEF",
        "PGR",
        "PGT",
        "XPE",
        "P4G",
        "P3A",
        "JEF",
        # Reducing agents
        "DTT",
        "BME",
        # Detergents (crystallization)
        "LMT",
        # Preservatives & crystallization additives
        "IPH",
        "CRS",
        "SPM",
        "SPD",
        "MLI",
        "TLA",
        "ACE",
        "NH2",
        "MPO",
        "TME",
    }
)


def enum_name(value: object) -> str:
    text = str(value)
    if "." in text:
        return text.split(".")[-1]
    return text


def normalize_residue_names(
    residue_names: Sequence[str],
) -> NormalizedResidueSpan | None:
    if not residue_names:
        return None

    start = 0
    end = len(residue_names) - 1
    while start <= end and residue_names[start].upper() in COMMON_TERMINAL_CAPS:
        start += 1
    while end >= start and residue_names[end].upper() in COMMON_TERMINAL_CAPS:
        end -= 1
    if start > end:
        return None

    trimmed_names = [name.upper() for name in residue_names[start : end + 1]]
    sequence = gemmi.one_letter_code(trimmed_names)
    return NormalizedResidueSpan(
        sequence=sequence,
        residue_names=trimmed_names,
        trim_start=start,
        trim_end=end,
    )


def is_valid_peptide_sequence(
    sequence: str,
    peptide_filter: PeptideSequenceFilter,
) -> bool:
    if not (peptide_filter.min_length <= len(sequence) <= peptide_filter.max_length):
        return False

    standard_residues = sum(
        1 for residue in sequence if residue in peptide_filter.standard_codes
    )
    if standard_residues < peptide_filter.min_standard_residues:
        return False

    nonstandard_residues = len(sequence) - standard_residues
    return not (
        nonstandard_residues > len(sequence) * peptide_filter.max_nonstandard_fraction
    )


def extract_structure(
    residues: Sequence[gemmi.Residue],
) -> tuple[bytes, bytes, bytes]:
    """Extract 37-atom coordinates plus residue-level B-factor and occupancy as bytes."""
    coords = []
    b_factors = []
    occupancy = []

    for residue in residues:
        atom_coords = [[np.nan, np.nan, np.nan] for _ in ATOM_TYPES_37]
        atom_b = [255] * len(ATOM_TYPES_37)
        atom_occ = [255] * len(ATOM_TYPES_37)

        for atom in residue:
            if atom.name not in ATOM_TYPES_37:
                continue
            atom_index = ATOM_TYPES_37.index(atom.name)
            atom_coords[atom_index] = [atom.pos.x, atom.pos.y, atom.pos.z]

            if atom.b_iso is not None:
                b_val = min(max(0, round(atom.b_iso)), 254) if atom.b_iso >= 0 else 0
                atom_b[atom_index] = b_val

            if atom.occ is not None:
                occ_val = min(max(0, round(atom.occ * 100)), 100)
                atom_occ[atom_index] = occ_val

        coords.append(atom_coords)
        b_factors.append(atom_b)
        occupancy.append(atom_occ)

    coords_arr = np.asarray(coords, dtype=np.float16)
    b_factors_arr = np.asarray(b_factors, dtype=np.uint8)
    occupancy_arr = np.asarray(occupancy, dtype=np.uint8)

    return coords_arr.tobytes(), b_factors_arr.tobytes(), occupancy_arr.tobytes()


def collect_atom_positions(residues: Sequence[gemmi.Residue]) -> np.ndarray:
    positions = []
    for residue in residues:
        for atom in residue:
            positions.append([atom.pos.x, atom.pos.y, atom.pos.z])
    if not positions:
        return np.empty((0, 3), dtype=np.float32)
    return np.asarray(positions, dtype=np.float32)


def build_subchain_entity_map(
    structure: gemmi.Structure,
) -> tuple[dict[str, gemmi.Entity], dict[str, str]]:
    subchain_to_entity = {}
    subchain_to_entity_id = {}

    for entity in structure.entities:
        for subchain_id in entity.subchains:
            subchain_to_entity[subchain_id] = entity
            subchain_to_entity_id[subchain_id] = entity.name

    return subchain_to_entity, subchain_to_entity_id


def pair_key(pdb_id: str, pair: PairData) -> tuple[str, str, str, str, str]:
    return (
        pdb_id,
        pair["peptide"]["entity_id"],
        pair["peptide"]["chain"],
        pair["receptor"]["entity_id"],
        pair["receptor"]["chain"],
    )


def build_assembly_atoms(model: gemmi.Model) -> tuple[np.ndarray, list[str]]:
    all_coords = []
    atom_subchains = []
    seen_subchains = set()

    for chain in model:
        for span in chain.subchains():
            subchain_id = span.subchain_id()
            if subchain_id in seen_subchains:
                continue
            seen_subchains.add(subchain_id)

            for residue in span:
                for atom in residue:
                    all_coords.append([atom.pos.x, atom.pos.y, atom.pos.z])
                    atom_subchains.append(subchain_id)

    if not all_coords:
        return np.empty((0, 3), dtype=np.float32), []

    return np.asarray(all_coords, dtype=np.float32), atom_subchains


def is_amino_acid_receptor(entity: gemmi.Entity | None) -> bool:
    return bool(
        entity is not None
        and entity.entity_type == gemmi.EntityType.Polymer
        and entity.polymer_type in PEPTIDE_POLYMER_TYPES
    )


def classify_neighbor_subchain(
    subchain_id: str,
    peptide_subchain_id: str,
    subchain_to_entity: dict[str, gemmi.Entity],
) -> NeighborRecord:
    entity = subchain_to_entity.get(subchain_id)
    entity_id = entity.name if entity is not None else None
    entity_type = enum_name(entity.entity_type) if entity is not None else "Unknown"
    polymer_type = enum_name(entity.polymer_type) if entity is not None else "Unknown"
    counts_as_competing_receptor = (
        subchain_id != peptide_subchain_id and is_amino_acid_receptor(entity)
    )
    return {
        "subchain_id": subchain_id,
        "entity_id": entity_id,
        "entity_type": entity_type,
        "polymer_type": polymer_type,
        "is_solvent": entity is not None
        and entity.entity_type == gemmi.EntityType.Water,
        "counts_as_competing_receptor": counts_as_competing_receptor,
    }


def analyze_neighbors(
    peptide_atoms: np.ndarray,
    tree: KDTree,
    atom_subchains: list[str],
    peptide_subchain_id: str,
    subchain_to_entity: dict[str, gemmi.Entity],
    distance: float,
) -> tuple[list[NeighborRecord], list[str]]:
    if peptide_atoms.size == 0:
        return [], []

    neighbors = tree.query_ball_point(peptide_atoms, r=distance)
    nearby_subchains = set()

    for hits in neighbors:
        for atom_index in hits:
            subchain_id = atom_subchains[atom_index]
            if subchain_id != peptide_subchain_id:
                nearby_subchains.add(subchain_id)

    neighbor_records = [
        classify_neighbor_subchain(subchain_id, peptide_subchain_id, subchain_to_entity)
        for subchain_id in sorted(nearby_subchains)
    ]
    competing = [
        record["subchain_id"]
        for record in neighbor_records
        if record["counts_as_competing_receptor"]
    ]
    return neighbor_records, competing


def annotate_nonpolymer_neighbor_subtypes(
    model: gemmi.Model,
    neighbor_records: list[NeighborRecord],
) -> None:
    for record in neighbor_records:
        if record["entity_type"] != "NonPolymer":
            continue
        neighbor_subchain = model.get_subchain(record["subchain_id"])
        neighbor_residue_names = (
            sorted({residue.name for residue in neighbor_subchain})
            if neighbor_subchain
            else []
        )
        record["nonpolymer_subtype"] = (
            neighbor_residue_names[0]
            if len(neighbor_residue_names) == 1
            else (neighbor_residue_names or None)
        )


def mark_solvent_neighbors(neighbor_records: list[NeighborRecord]) -> None:
    for record in neighbor_records:
        if (
            record["entity_type"] == "NonPolymer"
            and record.get("nonpolymer_subtype") in NONPOLYMER_SOLVENT_SUBTYPES
        ):
            record["is_solvent"] = True


def choose_receptor_from_neighbors(
    neighbor_records: list[NeighborRecord],
    competing_neighbors: list[str],
) -> str | None:
    if not neighbor_records:
        return None

    non_solvent_neighbors = [
        record for record in neighbor_records if not record["is_solvent"]
    ]
    if not non_solvent_neighbors:
        return None

    if len(non_solvent_neighbors) > 1:
        return None

    if len(competing_neighbors) == 0:
        return None

    return competing_neighbors[0]


def process_assembly(
    zf: zipfile.ZipFile,
    filename: str,
    peptide_filter: PeptideSequenceFilter,
    min_receptor_len: int,
    distance: float,
    binding_filter: BindingFilter,
) -> ProcessResult:
    pdb_id = filename.split("-")[0].upper()

    try:
        with gzip.open(io.BytesIO(zf.read(filename)), "rt") as handle:
            doc = gemmi.cif.read_string(handle.read())
            structure = gemmi.make_structure_from_block(doc[0])
    except Exception:
        return None

    if len(structure) != 1:
        return None

    model = structure[0]
    subchain_to_entity, subchain_to_entity_id = build_subchain_entity_map(structure)
    assembly_coords, atom_subchains = build_assembly_atoms(model)
    if assembly_coords.size == 0:
        return None

    tree = KDTree(assembly_coords)
    entry: LmdbEntry = {"pdb_id": pdb_id, "entities": []}
    found_pairs = False
    pair_keys: set[tuple[str, str, str, str, str]] = set()

    for entity in structure.entities:
        if entity.entity_type != gemmi.EntityType.Polymer:
            continue
        if entity.polymer_type not in PEPTIDE_POLYMER_TYPES:
            continue

        entity_span = normalize_residue_names(entity.full_sequence)
        if entity_span is None:
            continue
        if not is_valid_peptide_sequence(entity_span.sequence, peptide_filter):
            continue

        entity_entry: EntityData = {
            "entity_id": entity.name,
            "sequence": entity_span.sequence,
            "residue_names": entity_span.residue_names,
            "pairs": [],
        }

        for subchain_id in entity.subchains:
            subchain = model.get_subchain(subchain_id)
            if subchain is None:
                continue

            residues = list(subchain)
            peptide_span = normalize_residue_names(
                [residue.name for residue in residues]
            )
            if peptide_span is None:
                continue
            if not is_valid_peptide_sequence(peptide_span.sequence, peptide_filter):
                continue

            peptide_residues = list(
                residues[peptide_span.trim_start : peptide_span.trim_end + 1]
            )
            peptide_atoms = collect_atom_positions(peptide_residues)
            neighbor_records, competing_neighbors = analyze_neighbors(
                peptide_atoms=peptide_atoms,
                tree=tree,
                atom_subchains=atom_subchains,
                peptide_subchain_id=subchain_id,
                subchain_to_entity=subchain_to_entity,
                distance=distance,
            )
            annotate_nonpolymer_neighbor_subtypes(model, neighbor_records)
            mark_solvent_neighbors(neighbor_records)
            receptor_subchain_id = choose_receptor_from_neighbors(
                neighbor_records=neighbor_records,
                competing_neighbors=competing_neighbors,
            )
            if receptor_subchain_id is None:
                continue
            receptor_chain = model.get_subchain(receptor_subchain_id)
            if receptor_chain is None:
                continue

            receptor_entity_id = subchain_to_entity_id.get(receptor_subchain_id, "")
            if not receptor_entity_id:
                continue

            receptor_residues = list(receptor_chain)
            receptor_span = normalize_residue_names(
                [residue.name for residue in receptor_residues]
            )
            if receptor_span is None:
                continue
            if len(receptor_span.sequence) < min_receptor_len:
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
                continue
            if (
                interface_metrics.peptide_contact_residues
                < binding_filter.min_contact_residues
            ):
                continue
            if (
                interface_metrics.peptide_contact_fraction
                < binding_filter.min_contact_fraction
            ):
                continue

            pair["peptide"]["interface_start"] = (
                interface_metrics.peptide_interface_start
            )
            pair["peptide"]["interface_end"] = interface_metrics.peptide_interface_end
            pair["peptide"]["contact_residues"] = (
                interface_metrics.peptide_contact_residues
            )
            pair["peptide"]["contact_fraction"] = (
                interface_metrics.peptide_contact_fraction
            )
            pair["peptide"]["mean_contact_atom_b_factor"] = (
                interface_metrics.mean_contact_atom_b_factor
            )
            pair["receptor"]["interface_start"] = (
                interface_metrics.receptor_interface_start
            )
            pair["receptor"]["interface_end"] = interface_metrics.receptor_interface_end
            pair["receptor"]["contact_residues"] = (
                interface_metrics.receptor_contact_residues
            )

            current_pair_key = pair_key(pdb_id, pair)
            if current_pair_key in pair_keys:
                raise DuplicatePairError(*current_pair_key)
            pair_keys.add(current_pair_key)
            entity_entry["pairs"].append(pair)
            found_pairs = True

        if entity_entry["pairs"]:
            entry["entities"].append(entity_entry)

    if not found_pairs:
        return None
    return pdb_id, entry
