from __future__ import annotations

import argparse
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import cast

import lmdb
import msgpack
from tqdm import tqdm

from pdb_mldata.curation import BestPairMetrics, calculate_best_pair_metrics
from pdb_mldata.lmdb_utils import (
    BestPairData,
    BestPairLmdbEntry,
    EntityData,
    EntitySummaryData,
    LmdbEntry,
    PairData,
    decode_chain_data,
    encode_best_pair_lmdb_entry,
)

DEFAULT_INPUT_LMDB_PATH = Path("data/pdb_mldata_binding.lmdb")
DEFAULT_OUTPUT_LMDB_PATH = Path("data/pdb_mldata_best_pair.lmdb")
DEFAULT_DISTANCE = 5.0
DEFAULT_MAX_CONTACT_ATOM_B_FACTOR = 70.0


@dataclass(frozen=True)
class RankedPair:
    pair: PairData
    metrics: BestPairMetrics


@dataclass
class BestPairStats:
    entries_read: int = 0
    entries_written: int = 0
    entities_read: int = 0
    pairs_read: int = 0
    pairs_written: int = 0
    surplus_pairs_removed: int = 0
    single_pair_entities: int = 0
    ranked_entities: int = 0
    stable_id_ties: int = 0


def unpack_lmdb_entry(entry_bytes: bytes) -> LmdbEntry:
    return cast(LmdbEntry, msgpack.unpackb(entry_bytes, raw=False))


def decode_pair_for_metrics(pair: PairData) -> PairData:
    return {
        "peptide": decode_chain_data(pair["peptide"]),
        "receptor": decode_chain_data(pair["receptor"]),
    }


def calculate_ranked_pair(
    pair: PairData,
    distance: float,
    max_contact_atom_b_factor: float,
) -> RankedPair:
    decoded_pair = decode_pair_for_metrics(pair)
    return RankedPair(
        pair=pair,
        metrics=calculate_best_pair_metrics(
            peptide=decoded_pair["peptide"],
            receptor=decoded_pair["receptor"],
            distance=distance,
            max_contact_atom_b_factor=max_contact_atom_b_factor,
        ),
    )


def quality_rank_key(ranked_pair: RankedPair) -> tuple[int, float, float, int, int]:
    metrics = ranked_pair.metrics
    return (
        -metrics.valid_contact_residues,
        -metrics.valid_contact_fraction,
        metrics.mean_valid_contact_atom_b_factor,
        -metrics.finite_peptide_residues,
        metrics.receptor_residues,
    )


def stable_rank_key(ranked_pair: RankedPair) -> tuple[str, str, str]:
    pair = ranked_pair.pair
    return (
        pair["peptide"]["chain"],
        pair["receptor"]["entity_id"],
        pair["receptor"]["chain"],
    )


def rank_key(
    ranked_pair: RankedPair,
) -> tuple[int, float, float, int, int, str, str, str]:
    return (*quality_rank_key(ranked_pair), *stable_rank_key(ranked_pair))


def select_best_pair(
    entity: EntityData,
    distance: float,
    max_contact_atom_b_factor: float,
    stats: BestPairStats,
) -> PairData:
    if not entity["pairs"]:
        raise ValueError(f"entity {entity['entity_id']} has no pairs")

    if len(entity["pairs"]) == 1:
        stats.single_pair_entities += 1
    else:
        stats.ranked_entities += 1
        stats.surplus_pairs_removed += len(entity["pairs"]) - 1

    ranked_pairs = [
        calculate_ranked_pair(
            pair=pair,
            distance=distance,
            max_contact_atom_b_factor=max_contact_atom_b_factor,
        )
        for pair in entity["pairs"]
    ]
    best_ranked_pair = min(ranked_pairs, key=rank_key)
    best_quality_key = quality_rank_key(best_ranked_pair)
    quality_tie_count = sum(
        1
        for ranked_pair in ranked_pairs
        if quality_rank_key(ranked_pair) == best_quality_key
    )
    if quality_tie_count > 1:
        stats.stable_id_ties += 1

    return best_ranked_pair.pair


def build_best_pair_data(entity: EntityData, pair: PairData) -> BestPairData:
    entity_summary: EntitySummaryData = {
        "entity_id": entity["entity_id"],
        "sequence": entity["sequence"],
        "residue_names": entity["residue_names"],
    }
    return {
        "entity": entity_summary,
        "peptide": pair["peptide"],
        "receptor": pair["receptor"],
    }


def select_best_pairs_for_entry(
    entry: LmdbEntry,
    distance: float,
    max_contact_atom_b_factor: float,
    stats: BestPairStats,
) -> BestPairLmdbEntry:
    best_pairs: list[BestPairData] = []

    stats.entries_read += 1
    stats.entities_read += len(entry["entities"])
    for entity in entry["entities"]:
        stats.pairs_read += len(entity["pairs"])
        pair = select_best_pair(
            entity=entity,
            distance=distance,
            max_contact_atom_b_factor=max_contact_atom_b_factor,
            stats=stats,
        )
        best_pairs.append(build_best_pair_data(entity=entity, pair=pair))

    return {"pdb_id": entry["pdb_id"], "pairs": best_pairs}


def build_best_pair_lmdb(
    input_path: Path,
    output_path: Path,
    distance: float,
    max_contact_atom_b_factor: float,
    limit: int | None,
) -> BestPairStats:
    if output_path.exists():
        shutil.rmtree(output_path)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    input_env = lmdb.open(str(input_path), readonly=True, lock=False)
    output_env = lmdb.open(str(output_path), map_size=10**11)
    stats = BestPairStats()

    with input_env.begin() as input_txn:
        total_entries = input_txn.stat()["entries"]
        if limit is not None:
            total_entries = min(total_entries, limit)

        cursor = input_txn.cursor()
        progress = tqdm(cursor, total=total_entries, desc="Selecting best pairs")
        for entry_index, (key, value) in enumerate(progress):
            if limit is not None and entry_index >= limit:
                break

            entry = unpack_lmdb_entry(cast(bytes, value))
            best_pair_entry = select_best_pairs_for_entry(
                entry=entry,
                distance=distance,
                max_contact_atom_b_factor=max_contact_atom_b_factor,
                stats=stats,
            )
            stats.entries_written += 1
            stats.pairs_written += len(best_pair_entry["pairs"])
            with output_env.begin(write=True) as output_txn:
                output_txn.put(
                    cast(bytes, key),
                    encode_best_pair_lmdb_entry(best_pair_entry),
                )

            progress.set_postfix(
                pairs=stats.pairs_written,
                removed=stats.surplus_pairs_removed,
                ties=stats.stable_id_ties,
            )

    input_env.close()
    output_env.close()
    return stats


def validate_parameters(
    input_path: Path,
    output_path: Path,
    distance: float,
    max_contact_atom_b_factor: float,
    limit: int | None,
) -> None:
    """Validate CLI parameters before destructive LMDB creation starts."""
    if not input_path.exists():
        raise ValueError(f"{input_path} not found")
    if not input_path.is_dir():
        raise ValueError("input_path must be an LMDB directory path")
    if input_path.resolve() == output_path.resolve():
        raise ValueError("input_path and output_path must be different LMDB paths")
    if output_path.exists() and output_path.is_file():
        raise ValueError("output_path must be an LMDB directory path, not a file")
    if distance <= 0:
        raise ValueError("--distance must be greater than 0")
    if max_contact_atom_b_factor <= 0:
        raise ValueError("--max-contact-atom-b-factor must be greater than 0")
    if limit is not None and limit < 1:
        raise ValueError("--limit must be at least 1 when provided")


def print_stats(stats: BestPairStats, output_path: Path) -> None:
    print(f"Wrote best-pair LMDB to {output_path}")
    print(f"PDB entries: {stats.entries_read} -> {stats.entries_written}")
    print(f"Peptide entities read: {stats.entities_read}")
    print(f"Peptide/receptor chain pairs: {stats.pairs_read} -> {stats.pairs_written}")
    print(f"Surplus pairs removed: {stats.surplus_pairs_removed}")
    print(f"Entities already containing one pair: {stats.single_pair_entities}")
    print(f"Entities ranked for best pair: {stats.ranked_entities}")
    print(f"Ranking ties resolved by stable IDs: {stats.stable_id_ties}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Select one deterministic best peptide/receptor pair per peptide entity."
    )
    parser.add_argument(
        "input_path",
        nargs="?",
        type=Path,
        default=DEFAULT_INPUT_LMDB_PATH,
        help=(
            "Path to the binding-filtered LMDB directory "
            "(default: data/pdb_mldata_binding.lmdb)"
        ),
    )
    parser.add_argument(
        "output_path",
        nargs="?",
        type=Path,
        default=DEFAULT_OUTPUT_LMDB_PATH,
        help=(
            "Path to the output best-pair LMDB directory "
            "(default: data/pdb_mldata_best_pair.lmdb)"
        ),
    )
    parser.add_argument(
        "--distance",
        type=float,
        default=DEFAULT_DISTANCE,
        help="Contact distance in Angstroms (default: 5.0)",
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
        help="Process only the first N LMDB entries for testing",
    )
    args = parser.parse_args()

    try:
        validate_parameters(
            input_path=args.input_path,
            output_path=args.output_path,
            distance=args.distance,
            max_contact_atom_b_factor=args.max_contact_atom_b_factor,
            limit=args.limit,
        )
    except ValueError as exc:
        parser.error(str(exc))

    stats = build_best_pair_lmdb(
        input_path=args.input_path,
        output_path=args.output_path,
        distance=args.distance,
        max_contact_atom_b_factor=args.max_contact_atom_b_factor,
        limit=args.limit,
    )
    print_stats(stats=stats, output_path=args.output_path)


if __name__ == "__main__":
    main()
