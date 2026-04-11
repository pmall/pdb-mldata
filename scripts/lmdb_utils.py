import msgpack
import numpy as np

def decode_chain_data(chain_data: dict) -> dict:
    """
    Decodes the numpy bytes for a single chain back to usable float arrays.
    Restores NaN values efficiently from the 255 sentinels and rescales occupancy.
    """
    decoded = chain_data.copy()
    
    # 1. Structure is directly float16, NaN already supported natively
    structure = np.frombuffer(chain_data["structure"], dtype=np.float16).copy()
    decoded["structure"] = structure.reshape(-1, 37, 3)

    # 2. B-factors were compressed to uint8, 255 meant NaN
    b_factors_uint8 = np.frombuffer(chain_data["b_factors"], dtype=np.uint8).reshape(-1, 37)
    b_factors = b_factors_uint8.astype(np.float32)
    b_factors[b_factors_uint8 == 255] = np.nan
    decoded["b_factors"] = b_factors

    # 3. Occupancy was * 100 and compressed to uint8, 255 meant NaN
    occupancy_uint8 = np.frombuffer(chain_data["occupancy"], dtype=np.uint8).reshape(-1, 37)
    occupancy = occupancy_uint8.astype(np.float32)
    # Perform division, but restore NaNs using a boolean mask
    is_missing = (occupancy_uint8 == 255)
    occupancy = occupancy / 100.0
    occupancy[is_missing] = np.nan
    decoded["occupancy"] = occupancy

    return decoded


def decode_lmdb_entry(entry_bytes: bytes) -> dict:
    """
    Unpacks a msgpack-compressed LMDB payload and fully restores all
    compressed arrays, sequence offsets, and structural NaN values.
    """
    data = msgpack.unpackb(entry_bytes)
    
    # Recursively decode structures inside the dataset
    for entity in data.get("entities", []):
        for pair in entity.get("pairs", []):
            if "peptide" in pair:
                pair["peptide"] = decode_chain_data(pair["peptide"])
            if "receptor" in pair:
                pair["receptor"] = decode_chain_data(pair["receptor"])
                
    return data
