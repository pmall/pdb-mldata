from __future__ import annotations

from typing import TypeAlias, TypedDict, cast

import numpy as np
import msgpack

BytePayload: TypeAlias = bytes | np.ndarray


class ChainData(TypedDict):
    entity_id: str
    chain: str
    sequence: str
    structure: BytePayload
    b_factors: BytePayload
    occupancy: BytePayload


class PairData(TypedDict):
    peptide: ChainData
    receptor: ChainData


class EntityData(TypedDict):
    entity_id: str
    sequence: str
    pairs: list[PairData]


class LmdbEntry(TypedDict):
    pdb_id: str
    entities: list[EntityData]


def encode_lmdb_entry(entry: LmdbEntry) -> bytes:
    """Serialize one LMDB entry with the project-wide msgpack convention."""
    return cast(bytes, msgpack.packb(entry, use_bin_type=True))


def decode_chain_data(chain_data: ChainData) -> ChainData:
    """Decode one chain payload back to arrays and restore sentinel NaN values."""
    decoded = chain_data.copy()

    structure = np.frombuffer(
        cast(bytes, chain_data["structure"]), dtype=np.float16
    ).copy()
    decoded["structure"] = structure.reshape(-1, 37, 3)

    b_factors_uint8 = np.frombuffer(
        cast(bytes, chain_data["b_factors"]), dtype=np.uint8
    ).reshape(-1, 37)
    b_factors = b_factors_uint8.astype(np.float32)
    b_factors[b_factors_uint8 == 255] = np.nan
    decoded["b_factors"] = b_factors

    occupancy_uint8 = np.frombuffer(
        cast(bytes, chain_data["occupancy"]), dtype=np.uint8
    ).reshape(-1, 37)
    occupancy = occupancy_uint8.astype(np.float32)
    is_missing = occupancy_uint8 == 255
    occupancy = occupancy / 100.0
    occupancy[is_missing] = np.nan
    decoded["occupancy"] = occupancy

    return decoded


def decode_lmdb_entry(entry_bytes: bytes) -> LmdbEntry:
    """Unpack one LMDB payload and restore compressed structural arrays."""
    data = cast(LmdbEntry, msgpack.unpackb(entry_bytes))

    for entity in data.get("entities", []):
        for pair in entity.get("pairs", []):
            if "peptide" in pair:
                pair["peptide"] = decode_chain_data(pair["peptide"])
            if "receptor" in pair:
                pair["receptor"] = decode_chain_data(pair["receptor"])

    return data
