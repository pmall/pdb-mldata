"""
PDB Biological Assembly Parser & LMDB Builder.

The identification logic uses mmCIF label subchains end-to-end. We do not
merge author chains. Terminal non-standard residues are always trimmed from
both peptide and receptor chains, while internal receptor modifications are kept.
"""

import argparse
from collections import Counter
import gzip
import io
import json
import msgpack
import sys
import shutil
import zipfile
from pathlib import Path
from typing import Counter as CounterType, Dict, List, Optional, Sequence, Tuple

import gemmi
import lmdb
import numpy as np
from scipy.spatial import KDTree
from tqdm import tqdm

THREE_TO_ONE = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
}
STANDARD_AA_NAMES = set(THREE_TO_ONE)
ATOM_TYPES_37 = (
    "N", "CA", "C", "CB", "O", "CG", "CG1", "CG2", "OG", "OG1", "SG", "CD",
    "CD1", "CD2", "ND1", "ND2", "OD1", "OD2", "SD", "CE", "CE1", "CE2", "CE3",
    "NE", "NE1", "NE2", "OE1", "OE2", "CH2", "NH1", "NH2", "OH", "CZ", "CZ2",
    "CZ3", "NZ", "OXT"
)
PEPTIDE_POLYMER_TYPES = (gemmi.PolymerType.PeptideL, gemmi.PolymerType.PeptideD)

# Non-polymer subtypes that are crystallographic solvent/buffer artifacts with no
# biological role in the binding interaction (solvents, buffers, PEG, preservatives,
# reducing agents, detergents). Treated as solvent: excluded from the neighbor count.
NONPOLYMER_SOLVENT_SUBTYPES = frozenset({
    # Organic solvents & cryoprotectants
    "GOL", "EDO", "DMS", "MPD", "IPA", "EOH", "MOH", "ACN", "DCE", "DIO",
    "BU3", "TBU", "OCT", "BNZ", "ERY", "CCN", "ETA", "HEX", "OXY",
    # Buffer inorganic anions
    "SO4", "PO4", "NO3", "CO3", "ACT", "ACY", "FMT", "NH4", "NH3",
    "CO2", "AZI", "CIT", "TFA",
    # Organic buffer molecules
    "MES", "EPE", "TRS", "IMD", "URE", "B3P", "BTB",
    # PEG fragments & polyethers
    "PG4", "PEG", "PGV", "PGE", "1PE", "2PE", "3PE", "PE4", "PE5", "PE8",
    "P6G", "C8E", "PG0", "PGO", "15P", "16P", "1PG", "PEV", "PEE", "PEF",
    "PGR", "PGT", "XPE", "P4G", "P3A", "JEF",
    # Reducing agents
    "DTT", "BME",
    # Detergents (crystallization)
    "LMT",
    # Preservatives & crystallization additives
    "IPH", "CRS", "SPM", "SPD", "MLI", "TLA", "ACE", "NH2", "MPO", "TME",
})



def make_diagnostic_record(
    pdb_id: str,
    entity_id: Optional[str],
    peptide_subchain_id: Optional[str],
    status: str,
    skip_reason_code: Optional[str],
    skip_reason_detail: Optional[str],
    peptide: Optional[Dict] = None,
    neighbors: Optional[List[Dict]] = None,
    entry_peptide_screen: Optional[Dict] = None,
) -> Dict:
    return {
        "pdb_id": pdb_id,
        "entity_id": entity_id,
        "peptide_subchain_id": peptide_subchain_id,
        "status": status,
        "skip_reason_code": skip_reason_code,
        "skip_reason_detail": skip_reason_detail,
        "peptide": peptide,
        "neighbors": neighbors or [],
        "entry_peptide_screen": entry_peptide_screen,
    }


def make_peptide_details(
    raw_residue_count: int,
    trimmed_start: Optional[int] = None,
    trimmed_end: Optional[int] = None,
    trimmed_sequence: Optional[str] = None,
    trimmed_length: Optional[int] = None,
) -> Dict:
    return {
        "raw_residue_count": raw_residue_count,
        "trimmed_start": trimmed_start,
        "trimmed_end": trimmed_end,
        "trimmed_sequence": trimmed_sequence,
        "trimmed_length": trimmed_length,
    }


def enum_name(value) -> str:
    text = str(value)
    if "." in text:
        return text.split(".")[-1]
    return text


def is_standard_amino_acid(residue: gemmi.Residue) -> bool:
    return residue.name in STANDARD_AA_NAMES


def trim_terminal_caps(
    residues: Sequence[gemmi.Residue],
) -> Tuple[Optional[Tuple[List[gemmi.Residue], int, int, str]], Optional[str]]:
    """Return the contiguous standard-AA core after trimming terminal caps."""
    if not residues:
        return None, "empty_peptide_subchain"

    flags = [is_standard_amino_acid(residue) for residue in residues]
    if not any(flags):
        return None, "no_standard_residues_after_trim"

    start = next(i for i, flag in enumerate(flags) if flag)
    end = len(flags) - 1 - next(i for i, flag in enumerate(reversed(flags)) if flag)
    core = list(residues[start : end + 1])

    if not core:
        return None, "no_standard_residues_after_trim"
    if any(not is_standard_amino_acid(residue) for residue in core):
        return None, "internal_nonstandard_residue_in_peptide"

    sequence = "".join(THREE_TO_ONE[residue.name] for residue in core)
    return (core, start, end, sequence), None


def trim_receptor_terminal_caps(
    residues: Sequence[gemmi.Residue],
) -> Optional[Tuple[List[gemmi.Residue], str]]:
    """Always trim terminal non-standard caps, but keep internal modified residues."""
    if not residues:
        return None

    flags = [is_standard_amino_acid(residue) for residue in residues]
    if not any(flags):
        return None

    start = next(i for i, flag in enumerate(flags) if flag)
    end = len(flags) - 1 - next(i for i, flag in enumerate(reversed(flags)) if flag)
    core = list(residues[start : end + 1])
    if not core:
        return None

    sequence = "".join(THREE_TO_ONE.get(residue.name, "X") for residue in core)
    return core, sequence


def extract_structure(
    residues: Sequence[gemmi.Residue],
) -> Tuple[bytes, bytes, bytes]:
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
) -> Tuple[Dict[str, gemmi.Entity], Dict[str, str]]:
    subchain_to_entity = {}
    subchain_to_entity_id = {}

    for entity in structure.entities:
        for subchain_id in entity.subchains:
            subchain_to_entity[subchain_id] = entity
            subchain_to_entity_id[subchain_id] = entity.name

    return subchain_to_entity, subchain_to_entity_id


def build_assembly_atoms(model: gemmi.Model) -> Tuple[np.ndarray, List[str]]:
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


def is_amino_acid_receptor(entity: Optional[gemmi.Entity]) -> bool:
    return bool(
        entity is not None
        and entity.entity_type == gemmi.EntityType.Polymer
        and entity.polymer_type in PEPTIDE_POLYMER_TYPES
    )


def classify_neighbor_subchain(
    subchain_id: str,
    peptide_subchain_id: str,
    subchain_to_entity: Dict[str, gemmi.Entity],
) -> Dict:
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
        "is_solvent": entity is not None and entity.entity_type == gemmi.EntityType.Water,
        "counts_as_competing_receptor": counts_as_competing_receptor,
    }


def analyze_neighbors(
    peptide_atoms: np.ndarray,
    tree: KDTree,
    atom_subchains: List[str],
    peptide_subchain_id: str,
    subchain_to_entity: Dict[str, gemmi.Entity],
    distance: float
) -> Tuple[List[Dict], List[str]]:
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


def rejection_neighbors(neighbor_records: List[Dict]) -> List[Dict]:
    return [
        {
            "subchain_id": record["subchain_id"],
            "entity_id": record["entity_id"],
            "entity_type": record["entity_type"],
            "polymer_type": record["polymer_type"],
            "nonpolymer_subtype": record.get("nonpolymer_subtype"),
            "counts_as_competing_receptor": record["counts_as_competing_receptor"],
        }
        for record in neighbor_records
    ]


def choose_receptor_from_neighbors(
    neighbor_records: List[Dict],
    competing_neighbors: List[str],
    distance: float,
) -> Tuple[Optional[str], Optional[str], Optional[str]]:
    if not neighbor_records:
        return None, "no_neighbors_in_shell", f"No other subchains found within {distance}A of peptide"

    non_solvent_neighbors = [record for record in neighbor_records if not record["is_solvent"]]
    if not non_solvent_neighbors:
        return None, "only_solvent_neighbors", f"Only solvent subchains found within {distance}A of peptide"

    if len(non_solvent_neighbors) > 1:
        neighbor_types = ", ".join(
            f"{record['subchain_id']}:{record['entity_type']}/{record['polymer_type']}"
            for record in non_solvent_neighbors
        )
        return (
            None,
            "too_many_neighbors",
            f"{len(non_solvent_neighbors)} non-solvent nearby subchains in {distance}A shell: {neighbor_types}",
        )

    only_neighbor = non_solvent_neighbors[0]
    if len(competing_neighbors) == 0:
        return (
            None,
            "single_neighbor_wrong_type",
            "Single nearby subchain is not a competing peptide/polymer neighbor: "
            f"{only_neighbor['subchain_id']}:{only_neighbor['entity_type']}/{only_neighbor['polymer_type']}",
        )

    return competing_neighbors[0], None, None


def process_assembly(
    zf: zipfile.ZipFile,
    filename: str,
    min_len: int,
    max_len: int,
    distance: float,
) -> Tuple[Optional[Tuple[str, Dict]], List[Dict]]:
    pdb_id = filename.split("-")[0].upper()
    diagnostics = []
    found_valid_peptide_candidate = False
    peptide_screen_counts: CounterType[str] = Counter()
    peptide_polymer_subchains_seen = 0

    try:
        with gzip.open(io.BytesIO(zf.read(filename)), "rt") as handle:
            doc = gemmi.cif.read_string(handle.read())
            structure = gemmi.make_structure_from_block(doc[0])
    except Exception as exc:
        diagnostics.append(
            make_diagnostic_record(
                pdb_id=pdb_id,
                entity_id=None,
                peptide_subchain_id=None,
                status="skipped",
                skip_reason_code="read_error",
                skip_reason_detail=str(exc),
            )
        )
        return None, diagnostics

    if len(structure) != 1:
        diagnostics.append(
            make_diagnostic_record(
                pdb_id=pdb_id,
                entity_id=None,
                peptide_subchain_id=None,
                status="skipped",
                skip_reason_code="multiple_models",
                skip_reason_detail=f"Structure has {len(structure)} models",
            )
        )
        return None, diagnostics

    model = structure[0]
    subchain_to_entity, subchain_to_entity_id = build_subchain_entity_map(structure)
    assembly_coords, atom_subchains = build_assembly_atoms(model)
    if assembly_coords.size == 0:
        diagnostics.append(
            make_diagnostic_record(
                pdb_id=pdb_id,
                entity_id=None,
                peptide_subchain_id=None,
                status="skipped",
                skip_reason_code="no_assembly_atoms",
                skip_reason_detail="No atoms found across assembly subchains",
            )
        )
        return None, diagnostics

    tree = KDTree(assembly_coords)
    entry = {"pdb_id": pdb_id, "entities": []}
    found_pairs = False

    for entity in structure.entities:
        if entity.entity_type != gemmi.EntityType.Polymer:
            continue
        if entity.polymer_type not in PEPTIDE_POLYMER_TYPES:
            continue

        entity_entry = {"entity_id": entity.name, "sequence": "", "pairs": []}
        entity_sequence = None

        for subchain_id in entity.subchains:
            peptide_polymer_subchains_seen += 1
            subchain = model.get_subchain(subchain_id)
            if subchain is None:
                peptide_screen_counts["missing_subchain"] += 1
                continue

            residues = list(subchain)
            peptide_core, peptide_failure = trim_terminal_caps(residues)
            if peptide_core is None:
                peptide_screen_counts[peptide_failure] += 1
                continue

            peptide_residues, trim_start, trim_end, peptide_sequence = peptide_core
            peptide_details = make_peptide_details(
                raw_residue_count=len(residues),
                trimmed_start=trim_start,
                trimmed_end=trim_end,
                trimmed_sequence=peptide_sequence,
                trimmed_length=len(peptide_sequence),
            )
            if not (min_len <= len(peptide_sequence) <= max_len):
                if len(peptide_sequence) < min_len:
                    peptide_screen_counts["length_too_short"] += 1
                else:
                    peptide_screen_counts["length_too_long"] += 1
                continue

            found_valid_peptide_candidate = True
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
                    _sc = model.get_subchain(_rec["subchain_id"])
                    _names = sorted({r.name for r in _sc}) if _sc else []
                    _rec["nonpolymer_subtype"] = _names[0] if len(_names) == 1 else (_names or None)
            # Flag non-polymer neighbors whose subtype is a known crystallographic
            # solvent artifact — they play no role in binding and should not count
            # as competing neighbors.
            for _rec in neighbor_records:
                if (
                    _rec["entity_type"] == "NonPolymer"
                    and _rec.get("nonpolymer_subtype") in NONPOLYMER_SOLVENT_SUBTYPES
                ):
                    _rec["is_solvent"] = True
            receptor_subchain_id, neighbor_reason_code, neighbor_reason_detail = choose_receptor_from_neighbors(
                neighbor_records=neighbor_records,
                competing_neighbors=competing_neighbors,
                distance=distance,
            )
            if receptor_subchain_id is None:
                diagnostics.append(
                    make_diagnostic_record(
                        pdb_id=pdb_id,
                        entity_id=entity.name,
                        peptide_subchain_id=subchain_id,
                        status="skipped",
                        skip_reason_code=neighbor_reason_code,
                        skip_reason_detail=neighbor_reason_detail,
                        peptide=peptide_details,
                        neighbors=rejection_neighbors(neighbor_records),
                    )
                )
                continue
            receptor_chain = model.get_subchain(receptor_subchain_id)
            if receptor_chain is None:
                diagnostics.append(
                    make_diagnostic_record(
                        pdb_id=pdb_id,
                        entity_id=entity.name,
                        peptide_subchain_id=subchain_id,
                        status="skipped",
                        skip_reason_code="missing_receptor_subchain",
                        skip_reason_detail=f"Receptor subchain {receptor_subchain_id} not found in model",
                        peptide=peptide_details,
                        neighbors=rejection_neighbors(neighbor_records),
                    )
                )
                continue

            receptor_core = trim_receptor_terminal_caps(list(receptor_chain))
            if receptor_core is None:
                diagnostics.append(
                    make_diagnostic_record(
                        pdb_id=pdb_id,
                        entity_id=entity.name,
                        peptide_subchain_id=subchain_id,
                        status="skipped",
                        skip_reason_code="receptor_no_standard_residue_after_trim",
                        skip_reason_detail=(
                            f"Receptor subchain {receptor_subchain_id} has no standard residues after trimming"
                        ),
                        peptide=peptide_details,
                        neighbors=rejection_neighbors(neighbor_records),
                    )
                )
                continue

            receptor_residues, receptor_sequence = receptor_core
            receptor_entity_id = subchain_to_entity_id.get(receptor_subchain_id, "")
            if not receptor_entity_id:
                diagnostics.append(
                    make_diagnostic_record(
                        pdb_id=pdb_id,
                        entity_id=entity.name,
                        peptide_subchain_id=subchain_id,
                        status="skipped",
                        skip_reason_code="missing_receptor_entity_id",
                        skip_reason_detail=f"No entity id found for receptor subchain {receptor_subchain_id}",
                        peptide=peptide_details,
                        neighbors=rejection_neighbors(neighbor_records),
                    )
                )
                continue

            peptide_structure, peptide_b, peptide_occ = extract_structure(peptide_residues)
            receptor_structure, receptor_b, receptor_occ = extract_structure(receptor_residues)

            pair = {
                "peptide": {
                    "entity_id": entity.name,
                    "chain": subchain_id,
                    "sequence": peptide_sequence,
                    "structure": peptide_structure,
                    "b_factors": peptide_b,
                    "occupancy": peptide_occ,
                },
                "receptor": {
                    "entity_id": receptor_entity_id,
                    "chain": receptor_subchain_id,
                    "sequence": receptor_sequence,
                    "structure": receptor_structure,
                    "b_factors": receptor_b,
                    "occupancy": receptor_occ,
                },
            }
            entity_entry["pairs"].append(pair)
            found_pairs = True
            diagnostics.append(
                make_diagnostic_record(
                    pdb_id=pdb_id,
                    entity_id=entity.name,
                    peptide_subchain_id=subchain_id,
                    status="passed",
                    skip_reason_code=None,
                    skip_reason_detail=None,
                    peptide=peptide_details,
                    neighbors=[],
                )
            )

            if entity_sequence is None:
                entity_sequence = peptide_sequence
            elif entity_sequence != peptide_sequence:
                entity_sequence = ""

        if entity_entry["pairs"]:
            entity_entry["sequence"] = entity_sequence or ""
            entry["entities"].append(entity_entry)

    if not found_valid_peptide_candidate:
        screen_summary = {
            "polymer_peptide_subchains_seen": peptide_polymer_subchains_seen,
            "missing_subchain": peptide_screen_counts["missing_subchain"],
            "empty_peptide_subchain": peptide_screen_counts["empty_peptide_subchain"],
            "no_standard_residues_after_trim": peptide_screen_counts["no_standard_residues_after_trim"],
            "internal_nonstandard_residue_in_peptide": peptide_screen_counts["internal_nonstandard_residue_in_peptide"],
            "length_too_short": peptide_screen_counts["length_too_short"],
            "length_too_long": peptide_screen_counts["length_too_long"],
        }
        diagnostics.append(
            make_diagnostic_record(
                pdb_id=pdb_id,
                entity_id=None,
                peptide_subchain_id=None,
                status="skipped",
                skip_reason_code="no_valid_peptide_in_entry",
                skip_reason_detail=(
                    f"No peptide subchain passed terminal-cap trimming and length filter [{min_len}, {max_len}]"
                ),
                entry_peptide_screen=screen_summary,
            )
        )

    if not found_pairs:
        return None, diagnostics
    return (pdb_id, entry), diagnostics


def build_lmdb(
    zip_path: Path,
    output_path: Path,
    min_len: int,
    max_len: int,
    distance: float,
    limit: Optional[int],
    skip_log_path: Optional[Path],
) -> None:
    if not zip_path.exists():
        print(f"Error: {zip_path} not found")
        sys.exit(1)

    if output_path.exists():
        print(f"Clearing existing database at {output_path}...")
        shutil.rmtree(output_path)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    if skip_log_path is not None:
        skip_log_path.parent.mkdir(parents=True, exist_ok=True)
    env = lmdb.open(str(output_path), map_size=10**11)

    total_pairs = 0
    total_entities = 0
    assemblies_read = 0
    entries_with_valid_peptide = 0
    entries_without_valid_peptide = 0
    candidate_counter = 0
    accepted_candidate_counter = 0
    assembly_skip_reason_counts: CounterType[str] = Counter()
    candidate_skip_reason_counts: CounterType[str] = Counter()

    with zipfile.ZipFile(zip_path, "r") as zf, (
        skip_log_path.open("w") if skip_log_path is not None else open("/dev/null", "w")
    ) as skip_log_handle:
        files = [name for name in zf.namelist() if name.endswith(".cif.gz")]
        if limit is not None:
            files = files[:limit]

        print(f"Processing {len(files)} entries...")
        progress = tqdm(files, desc="Building LMDB")

        for filename in progress:
            assemblies_read += 1
            result, diagnostics = process_assembly(zf, filename, min_len, max_len, distance)

            for diagnostic in diagnostics:
                if skip_log_path is not None:
                    skip_log_handle.write(json.dumps(diagnostic) + "\n")
                if diagnostic["skip_reason_code"] == "no_valid_peptide_in_entry":
                    entries_without_valid_peptide += 1
                    assembly_skip_reason_counts[diagnostic["skip_reason_code"]] += 1
                    continue

                if diagnostic["entity_id"] is not None or diagnostic["peptide_subchain_id"] is not None:
                    candidate_counter += 1
                    if diagnostic["status"] == "passed":
                        accepted_candidate_counter += 1
                    elif diagnostic["skip_reason_code"] is not None:
                        candidate_skip_reason_counts[diagnostic["skip_reason_code"]] += 1
                elif diagnostic["status"] == "skipped" and diagnostic["skip_reason_code"] is not None:
                    assembly_skip_reason_counts[diagnostic["skip_reason_code"]] += 1

            if any(
                diagnostic["entity_id"] is not None or diagnostic["peptide_subchain_id"] is not None
                for diagnostic in diagnostics
            ):
                entries_with_valid_peptide += 1

            if result is None:
                progress.set_postfix(ents=total_entities, pairs=total_pairs)
                continue

            pdb_id, data = result
            total_entities += len(data["entities"])
            for entity in data["entities"]:
                total_pairs += len(entity["pairs"])
            progress.set_postfix(ents=total_entities, pairs=total_pairs)

            with env.begin(write=True) as txn:
                txn.put(pdb_id.encode(), msgpack.packb(data, use_bin_type=True))

    print(f"\nDone: {total_pairs} pairs from {total_entities} entities")
    print(f"Assemblies read: {assemblies_read}")
    print(f"Entries with at least one valid peptide candidate: {entries_with_valid_peptide}")
    print(f"Entries with no valid peptide candidate: {entries_without_valid_peptide}")
    print(f"Valid peptide candidates evaluated: {candidate_counter}")
    print(f"Accepted peptide/receptor pairs: {accepted_candidate_counter}")
    if assembly_skip_reason_counts:
        print("Assembly-level skip reasons:")
        for reason_code, count in assembly_skip_reason_counts.most_common():
            print(f"  {reason_code}: {count}")
    if candidate_skip_reason_counts:
        print("Candidate-level rejection reasons:")
        for reason_code, count in candidate_skip_reason_counts.most_common():
            print(f"  {reason_code}: {count}")
    env.close()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("zip_path", nargs="?", type=Path, default="data/assemblies.zip")
    parser.add_argument("output_path", nargs="?", type=Path, default="data/pdb_mldata.lmdb")
    parser.add_argument("--min-length", type=int, default=4)
    parser.add_argument("--max-length", type=int, default=32)
    parser.add_argument("--distance", type=float, default=5.0)
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument("--skip-log-path", type=Path, default=None)
    args = parser.parse_args()

    build_lmdb(
        zip_path=args.zip_path,
        output_path=args.output_path,
        min_len=args.min_length,
        max_len=args.max_length,
        distance=args.distance,
        limit=args.limit,
        skip_log_path=args.skip_log_path,
    )


if __name__ == "__main__":
    main()
