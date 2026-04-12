from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import cast

import lmdb
import numpy as np
from tqdm import tqdm

from pdb_mldata.lmdb_utils import ChainData, EntityData, decode_lmdb_entry

DEFAULT_LMDB_PATH = Path("data/pdb_mldata.lmdb")
VALID_SEQUENCE_CODES = frozenset("ACDEFGHIKLMNPQRSTUVWYOX")


@dataclass(frozen=True)
class SanityStats:
    entries: int
    entities: int
    pairs: int
    chains: int
    failures: int


def validate_parameters(db_path: Path, max_failures: int) -> None:
    """Validate CLI parameters before scanning the LMDB."""
    if not db_path.exists():
        raise ValueError(f"{db_path} not found")
    if not db_path.is_dir():
        raise ValueError("db_path must be an LMDB directory path")
    if max_failures < 1:
        raise ValueError("--max-failures must be at least 1")


def add_failure(
    failures: list[str],
    max_failures: int,
    context: str,
    message: str,
) -> None:
    if len(failures) < max_failures:
        failures.append(f"{context}: {message}")


def check_sequence_fields(
    sequence: str,
    residue_names: list[str],
    context: str,
    failures: list[str],
    max_failures: int,
) -> None:
    if not sequence:
        add_failure(failures, max_failures, context, "empty sequence")
    if len(sequence) != len(residue_names):
        add_failure(
            failures,
            max_failures,
            context,
            f"sequence length {len(sequence)} != residue_names length {len(residue_names)}",
        )

    invalid_codes = sorted(set(sequence) - VALID_SEQUENCE_CODES)
    if invalid_codes:
        add_failure(
            failures,
            max_failures,
            context,
            f"invalid sequence codes: {','.join(invalid_codes)}",
        )


def check_chain_arrays(
    chain: ChainData,
    context: str,
    failures: list[str],
    max_failures: int,
) -> None:
    sequence_length = len(chain["sequence"])
    structure = cast(np.ndarray, chain["structure"])
    b_factors = cast(np.ndarray, chain["b_factors"])
    occupancy = cast(np.ndarray, chain["occupancy"])

    expected_structure_shape = (sequence_length, 37, 3)
    expected_quality_shape = (sequence_length, 37)
    if structure.shape != expected_structure_shape:
        add_failure(
            failures,
            max_failures,
            context,
            f"structure shape {structure.shape} != {expected_structure_shape}",
        )
    if b_factors.shape != expected_quality_shape:
        add_failure(
            failures,
            max_failures,
            context,
            f"b_factors shape {b_factors.shape} != {expected_quality_shape}",
        )
    if occupancy.shape != expected_quality_shape:
        add_failure(
            failures,
            max_failures,
            context,
            f"occupancy shape {occupancy.shape} != {expected_quality_shape}",
        )


def check_chain(
    chain: ChainData,
    context: str,
    failures: list[str],
    max_failures: int,
) -> None:
    check_sequence_fields(
        sequence=chain["sequence"],
        residue_names=chain["residue_names"],
        context=context,
        failures=failures,
        max_failures=max_failures,
    )
    check_chain_arrays(
        chain=chain,
        context=context,
        failures=failures,
        max_failures=max_failures,
    )


def check_entity(
    entity: EntityData,
    pdb_id: str,
    failures: list[str],
    max_failures: int,
) -> int:
    entity_context = f"{pdb_id} entity {entity['entity_id']}"
    check_sequence_fields(
        sequence=entity["sequence"],
        residue_names=entity["residue_names"],
        context=entity_context,
        failures=failures,
        max_failures=max_failures,
    )

    chain_count = 0
    for pair_index, pair in enumerate(entity["pairs"]):
        peptide_context = (
            f"{entity_context} pair {pair_index} peptide {pair['peptide']['chain']}"
        )
        receptor_context = (
            f"{entity_context} pair {pair_index} receptor {pair['receptor']['chain']}"
        )
        check_chain(
            chain=pair["peptide"],
            context=peptide_context,
            failures=failures,
            max_failures=max_failures,
        )
        check_chain(
            chain=pair["receptor"],
            context=receptor_context,
            failures=failures,
            max_failures=max_failures,
        )
        chain_count += 2

    return chain_count


def check_lmdb_sanity(db_path: Path, max_failures: int) -> SanityStats:
    failures: list[str] = []
    entry_count = 0
    entity_count = 0
    pair_count = 0
    chain_count = 0

    env = lmdb.open(str(db_path), readonly=True, lock=False)
    with env.begin() as txn:
        total_entries = txn.stat()["entries"]
        cursor = txn.cursor()
        for key, value in tqdm(cursor, total=total_entries, desc="Checking LMDB"):
            entry = decode_lmdb_entry(cast(bytes, value))
            pdb_id = entry.get("pdb_id", cast(bytes, key).decode())
            entry_count += 1

            if not entry["entities"]:
                add_failure(
                    failures,
                    max_failures,
                    pdb_id,
                    "entry has no peptide entities",
                )

            entity_count += len(entry["entities"])
            for entity in entry["entities"]:
                pair_count += len(entity["pairs"])
                chain_count += check_entity(
                    entity=entity,
                    pdb_id=pdb_id,
                    failures=failures,
                    max_failures=max_failures,
                )

    env.close()

    print(f"Entries: {entry_count}")
    print(f"Peptide entities: {entity_count}")
    print(f"Pairs: {pair_count}")
    print(f"Chains checked: {chain_count}")
    print(f"Failures: {len(failures)}")

    if failures:
        print("")
        print("Failures:")
        for failure in failures:
            print(f"- {failure}")
        if len(failures) == max_failures:
            print(f"- Stopped reporting after --max-failures={max_failures}")

    return SanityStats(
        entries=entry_count,
        entities=entity_count,
        pairs=pair_count,
        chains=chain_count,
        failures=len(failures),
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Check strict LMDB sequence, residue-name, and structure invariants."
    )
    parser.add_argument(
        "db_path",
        nargs="?",
        type=Path,
        default=DEFAULT_LMDB_PATH,
        help="Path to the LMDB database directory (default: data/pdb_mldata.lmdb)",
    )
    parser.add_argument(
        "--max-failures",
        type=int,
        default=50,
        help="Maximum number of failures to print before suppressing details",
    )
    args = parser.parse_args()

    try:
        validate_parameters(db_path=args.db_path, max_failures=args.max_failures)
    except ValueError as exc:
        parser.error(str(exc))

    stats = check_lmdb_sanity(
        db_path=args.db_path,
        max_failures=args.max_failures,
    )
    if stats.failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
