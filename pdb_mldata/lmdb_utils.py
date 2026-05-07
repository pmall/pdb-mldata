from __future__ import annotations

from typing import NotRequired, TypeAlias, TypedDict, cast

import numpy as np
import msgpack

BytePayload: TypeAlias = bytes | np.ndarray


class ChainData(TypedDict):
    entity_id: str
    chain: str
    sequence: str
    residue_names: list[str]
    structure: BytePayload
    b_factors: BytePayload
    occupancy: BytePayload
    interface_start: NotRequired[int]
    interface_end: NotRequired[int]
    contact_residues: NotRequired[int]
    contact_fraction: NotRequired[float]
    mean_contact_atom_b_factor: NotRequired[float]


class PairData(TypedDict):
    peptide: ChainData
    receptor: ChainData


class EntityData(TypedDict):
    entity_id: str
    sequence: str
    residue_names: list[str]
    pairs: list[PairData]


class LmdbEntry(TypedDict):
    pdb_id: str
    entities: list[EntityData]


class EntitySummaryData(TypedDict):
    entity_id: str
    sequence: str
    residue_names: list[str]


class BestPairData(TypedDict):
    entity: EntitySummaryData
    peptide: ChainData
    receptor: ChainData


class BestPairLmdbEntry(TypedDict):
    pdb_id: str
    pairs: list[BestPairData]


class ExportChainData(TypedDict):
    entity_id: str
    chain: str
    sequence: str
    residue_names: list[str]
    interface_start: NotRequired[int]
    interface_end: NotRequired[int]
    contact_residues: NotRequired[int]
    contact_fraction: NotRequired[float]
    mean_contact_atom_b_factor: NotRequired[float]


class ExportPairData(TypedDict):
    peptide: ExportChainData
    receptor: ExportChainData


class ExportEntityData(TypedDict):
    entity_id: str
    sequence: str
    residue_names: list[str]
    pairs: list[ExportPairData]


class ExportLmdbEntry(TypedDict):
    pdb_id: str
    entities: list[ExportEntityData]


class ExportBestPairData(TypedDict):
    entity: EntitySummaryData
    peptide: ExportChainData
    receptor: ExportChainData


class ExportBestPairLmdbEntry(TypedDict):
    pdb_id: str
    pairs: list[ExportBestPairData]


def encode_lmdb_entry(entry: LmdbEntry) -> bytes:
    """Serialize one LMDB entry with the project-wide msgpack convention."""
    return cast(bytes, msgpack.packb(entry, use_bin_type=True))


def encode_best_pair_lmdb_entry(entry: BestPairLmdbEntry) -> bytes:
    """Serialize one best-pair LMDB entry with the project-wide msgpack convention."""
    return cast(bytes, msgpack.packb(entry, use_bin_type=True))


def select_export_chain_data(chain_data: ChainData) -> ExportChainData:
    """Select chain fields needed by downstream database exports."""
    export_data: ExportChainData = {
        "entity_id": chain_data["entity_id"],
        "chain": chain_data["chain"],
        "sequence": chain_data["sequence"],
        "residue_names": chain_data["residue_names"],
    }
    if "interface_start" in chain_data:
        export_data["interface_start"] = chain_data["interface_start"]
    if "interface_end" in chain_data:
        export_data["interface_end"] = chain_data["interface_end"]
    if "contact_residues" in chain_data:
        export_data["contact_residues"] = chain_data["contact_residues"]
    if "contact_fraction" in chain_data:
        export_data["contact_fraction"] = chain_data["contact_fraction"]
    if "mean_contact_atom_b_factor" in chain_data:
        export_data["mean_contact_atom_b_factor"] = chain_data[
            "mean_contact_atom_b_factor"
        ]
    return export_data


def unpack_lmdb_entry_for_export(entry_bytes: bytes) -> ExportLmdbEntry:
    """Unpack one LMDB payload without decoding heavy structural byte arrays."""
    data = cast(LmdbEntry, msgpack.unpackb(entry_bytes, raw=False))
    export_entities: list[ExportEntityData] = []

    for entity in data.get("entities", []):
        export_pairs: list[ExportPairData] = []
        for pair in entity.get("pairs", []):
            export_pairs.append(
                {
                    "peptide": select_export_chain_data(pair["peptide"]),
                    "receptor": select_export_chain_data(pair["receptor"]),
                }
            )

        export_entities.append(
            {
                "entity_id": entity["entity_id"],
                "sequence": entity["sequence"],
                "residue_names": entity["residue_names"],
                "pairs": export_pairs,
            }
        )

    return {"pdb_id": data["pdb_id"], "entities": export_entities}


def unpack_best_pair_lmdb_entry_for_export(
    entry_bytes: bytes,
) -> ExportBestPairLmdbEntry:
    """Unpack one best-pair payload without decoding heavy structural byte arrays."""
    data = cast(BestPairLmdbEntry, msgpack.unpackb(entry_bytes, raw=False))
    export_pairs: list[ExportBestPairData] = []

    for pair in data.get("pairs", []):
        export_pairs.append(
            {
                "entity": pair["entity"],
                "peptide": select_export_chain_data(pair["peptide"]),
                "receptor": select_export_chain_data(pair["receptor"]),
            }
        )

    return {"pdb_id": data["pdb_id"], "pairs": export_pairs}


def decode_chain_data(chain_data: ChainData) -> ChainData:
    """Decode one chain payload back to backbone coordinate arrays."""
    decoded = chain_data.copy()

    structure = np.frombuffer(
        cast(bytes, chain_data["structure"]), dtype=np.float16
    ).copy()
    decoded["structure"] = structure.reshape(-1, 3, 3)

    decoded["b_factors"] = np.frombuffer(
        cast(bytes, chain_data["b_factors"]), dtype=np.float16
    ).reshape(-1, 3)

    decoded["occupancy"] = np.frombuffer(
        cast(bytes, chain_data["occupancy"]), dtype=np.float16
    ).reshape(-1, 3)

    return decoded


def decode_lmdb_entry(entry_bytes: bytes) -> LmdbEntry:
    """Unpack one LMDB payload and restore compressed structural arrays."""
    data = cast(LmdbEntry, msgpack.unpackb(entry_bytes, raw=False))

    for entity in data.get("entities", []):
        for pair in entity.get("pairs", []):
            if "peptide" in pair:
                pair["peptide"] = decode_chain_data(pair["peptide"])
            if "receptor" in pair:
                pair["receptor"] = decode_chain_data(pair["receptor"])

    return data


def decode_best_pair_lmdb_entry(entry_bytes: bytes) -> BestPairLmdbEntry:
    """Unpack one best-pair LMDB payload and restore compressed structural arrays."""
    data = cast(BestPairLmdbEntry, msgpack.unpackb(entry_bytes, raw=False))

    for pair in data.get("pairs", []):
        pair["peptide"] = decode_chain_data(pair["peptide"])
        pair["receptor"] = decode_chain_data(pair["receptor"])

    return data
