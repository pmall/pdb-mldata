"""
PDB Biological Assembly Parser & LMDB Builder.

The identification logic uses mmCIF label subchains end-to-end. We do not
merge author chains. Common terminal caps are trimmed from entity and chain
sequences, while internal non-standard residues are retained.
"""

from __future__ import annotations

import argparse
import gzip
import io
import shutil
import zipfile
from dataclasses import dataclass
from pathlib import Path
from typing import NotRequired, Sequence, TypeAlias, TypedDict

import gemmi
import lmdb
import numpy as np
from scipy.spatial import KDTree
from tqdm import tqdm

from pdb_mldata.lmdb_utils import EntityData, LmdbEntry, PairData, encode_lmdb_entry

DEFAULT_ASSEMBLIES_ZIP = Path("data/assemblies.zip")
DEFAULT_LMDB_PATH = Path("data/pdb_mldata.lmdb")
DEFAULT_MIN_LENGTH = 4
DEFAULT_MAX_LENGTH = 32
DEFAULT_DISTANCE = 5.0


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


@dataclass(frozen=True)
class NormalizedResidueSpan:
    sequence: str
    residue_names: list[str]
    trim_start: int
    trim_end: int


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
    min_len: int,
    max_len: int,
    distance: float,
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

    for entity in structure.entities:
        if entity.entity_type != gemmi.EntityType.Polymer:
            continue
        if entity.polymer_type not in PEPTIDE_POLYMER_TYPES:
            continue

        entity_span = normalize_residue_names(entity.full_sequence)
        if entity_span is None:
            continue
        if not (min_len <= len(entity_span.sequence) <= max_len):
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
            for _rec in neighbor_records:
                if _rec["entity_type"] == "NonPolymer":
                    neighbor_subchain = model.get_subchain(_rec["subchain_id"])
                    neighbor_residue_names = (
                        sorted({residue.name for residue in neighbor_subchain})
                        if neighbor_subchain
                        else []
                    )
                    _rec["nonpolymer_subtype"] = (
                        neighbor_residue_names[0]
                        if len(neighbor_residue_names) == 1
                        else (neighbor_residue_names or None)
                    )
            # Flag non-polymer neighbors whose subtype is a known crystallographic
            # solvent artifact — they play no role in binding and should not count
            # as competing neighbors.
            for _rec in neighbor_records:
                if (
                    _rec["entity_type"] == "NonPolymer"
                    and _rec.get("nonpolymer_subtype") in NONPOLYMER_SOLVENT_SUBTYPES
                ):
                    _rec["is_solvent"] = True
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
            entity_entry["pairs"].append(pair)
            found_pairs = True

        if entity_entry["pairs"]:
            entry["entities"].append(entity_entry)

    if not found_pairs:
        return None
    return pdb_id, entry


def build_lmdb(
    zip_path: Path,
    output_path: Path,
    min_len: int,
    max_len: int,
    distance: float,
    limit: int | None,
) -> None:
    if output_path.exists():
        shutil.rmtree(output_path)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    env = lmdb.open(str(output_path), map_size=10**11)

    total_pairs = 0
    total_entities = 0

    with zipfile.ZipFile(zip_path, "r") as zf:
        files = [name for name in zf.namelist() if name.endswith(".cif.gz")]
        if limit is not None:
            files = files[:limit]

        progress = tqdm(files, desc="Building LMDB")

        for filename in progress:
            result = process_assembly(zf, filename, min_len, max_len, distance)

            if result is None:
                progress.set_postfix(ents=total_entities, pairs=total_pairs)
                continue

            pdb_id, data = result
            total_entities += len(data["entities"])
            for entity in data["entities"]:
                total_pairs += len(entity["pairs"])
            progress.set_postfix(ents=total_entities, pairs=total_pairs)

            with env.begin(write=True) as txn:
                txn.put(pdb_id.encode(), encode_lmdb_entry(data))

    env.close()


def validate_parameters(
    zip_path: Path,
    output_path: Path,
    min_length: int,
    max_length: int,
    distance: float,
    limit: int | None,
) -> None:
    """Validate CLI parameters before destructive LMDB creation starts."""
    if not zip_path.exists():
        raise ValueError(f"{zip_path} not found")
    if zip_path.is_dir():
        raise ValueError("zip_path must be a file path, not a directory")
    if output_path.exists() and output_path.is_file():
        raise ValueError("output_path must be an LMDB directory path, not a file")
    if min_length < 1:
        raise ValueError("--min-length must be at least 1")
    if max_length < min_length:
        raise ValueError("--max-length must be greater than or equal to --min-length")
    if distance <= 0:
        raise ValueError("--distance must be greater than 0")
    if limit is not None and limit < 1:
        raise ValueError("--limit must be at least 1 when provided")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build the peptide/receptor LMDB from assembly ZIP data."
    )
    parser.add_argument(
        "zip_path",
        nargs="?",
        type=Path,
        default=DEFAULT_ASSEMBLIES_ZIP,
        help="Path to the input assemblies ZIP (default: data/assemblies.zip)",
    )
    parser.add_argument(
        "output_path",
        nargs="?",
        type=Path,
        default=DEFAULT_LMDB_PATH,
        help="Path to the output LMDB directory (default: data/pdb_mldata.lmdb)",
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
        "--distance",
        type=float,
        default=DEFAULT_DISTANCE,
        help="Neighbor distance in Angstroms (default: 5.0)",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Process only the first N entries for testing",
    )
    args = parser.parse_args()

    try:
        validate_parameters(
            zip_path=args.zip_path,
            output_path=args.output_path,
            min_length=args.min_length,
            max_length=args.max_length,
            distance=args.distance,
            limit=args.limit,
        )
    except ValueError as exc:
        parser.error(str(exc))

    build_lmdb(
        zip_path=args.zip_path,
        output_path=args.output_path,
        min_len=args.min_length,
        max_len=args.max_length,
        distance=args.distance,
        limit=args.limit,
    )


if __name__ == "__main__":
    main()
