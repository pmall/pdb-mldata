"""Select one deterministic best peptide/receptor pair per peptide entity.

Goal:
- Read raw peptide/receptor pairs from `data/pdb_mldata.lmdb`.
- Rank all saved chain-level pairs under each peptide entity.
- Write one selected pair per peptide entity to `data/pdb_mldata_best_pair.lmdb`.

Ranking rules:
- Fail clearly if an input peptide entity has no pairs.
- Prefer the pair with the most stored peptide contact residues.
- Then prefer the highest stored peptide contact fraction.
- Then prefer the lowest stored mean B-factor among contact peptide atoms.
- Then prefer the most finite peptide residues.
- Then prefer the shorter receptor sequence.
- Use stable IDs as the final deterministic tie-breaker: peptide chain ID,
  receptor entity ID, then receptor chain ID.

Output behavior:
- Delete an existing output LMDB folder before writing to avoid LMDB upsert
  issues.
- The output schema changes because each retained record now represents one
  selected peptide-chain/receptor-chain pair for one peptide entity.
- See `docs/storage_schemas.md` for the best-pair LMDB schema.

Default parameters:
- Input LMDB: `data/pdb_mldata.lmdb`.
- Output LMDB: `data/pdb_mldata_best_pair.lmdb`.
- Optional `--limit` for smoke verification.
"""

from __future__ import annotations

import argparse
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import cast

import lmdb
import msgpack
import numpy as np
from tqdm import tqdm

from pdb_mldata.curation import BestPairMetrics, count_finite_residues
from pdb_mldata.lmdb_utils import (
    BestPairData,
    BestPairLmdbEntry,
    ChainData,
    EntityData,
    EntitySummaryData,
    LmdbEntry,
    PairData,
    decode_chain_data,
    encode_best_pair_lmdb_entry,
)

DEFAULT_INPUT_LMDB_PATH = Path("data/pdb_mldata.lmdb")
DEFAULT_OUTPUT_LMDB_PATH = Path("data/pdb_mldata_best_pair.lmdb")


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


def require_int_field(chain: ChainData, chain_name: str, field_name: str) -> int:
    value = cast(dict[str, object], chain).get(field_name)
    if type(value) is not int:
        raise ValueError(f"{chain_name}.{field_name} must be stored as an integer")
    return value


def require_float_field(chain: ChainData, chain_name: str, field_name: str) -> float:
    value = cast(dict[str, object], chain).get(field_name)
    if isinstance(value, bool) or not isinstance(value, int | float):
        raise ValueError(f"{chain_name}.{field_name} must be stored as a number")
    return float(value)


def require_peptide_interface_fields(peptide: ChainData) -> tuple[int, float, float]:
    for field in (
        "contact_residues",
        "contact_fraction",
        "mean_contact_atom_b_factor",
    ):
        if field not in peptide:
            raise ValueError(f"peptide.{field} is missing from stored pair")
    return (
        require_int_field(peptide, "peptide", "contact_residues"),
        require_float_field(peptide, "peptide", "contact_fraction"),
        require_float_field(peptide, "peptide", "mean_contact_atom_b_factor"),
    )


def calculate_ranked_pair(pair: PairData) -> RankedPair:
    decoded_pair = decode_pair_for_metrics(pair)
    finite_peptide_residues = count_finite_residues(
        cast(np.ndarray, decoded_pair["peptide"]["structure"])
    )
    contact_residues, contact_fraction, mean_contact_atom_b_factor = (
        require_peptide_interface_fields(pair["peptide"])
    )
    return RankedPair(
        pair=pair,
        metrics=BestPairMetrics(
            valid_contact_residues=contact_residues,
            valid_contact_fraction=contact_fraction,
            mean_valid_contact_atom_b_factor=mean_contact_atom_b_factor,
            finite_peptide_residues=finite_peptide_residues,
            receptor_residues=len(pair["receptor"]["sequence"]),
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
    stats: BestPairStats,
) -> PairData:
    if not entity["pairs"]:
        raise ValueError(f"entity {entity['entity_id']} has no pairs")

    if len(entity["pairs"]) == 1:
        stats.single_pair_entities += 1
    else:
        stats.ranked_entities += 1
        stats.surplus_pairs_removed += len(entity["pairs"]) - 1

    ranked_pairs = [calculate_ranked_pair(pair=pair) for pair in entity["pairs"]]
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
    stats: BestPairStats,
) -> BestPairLmdbEntry:
    best_pairs: list[BestPairData] = []

    stats.entries_read += 1
    stats.entities_read += len(entry["entities"])
    for entity in entry["entities"]:
        stats.pairs_read += len(entity["pairs"])
        pair = select_best_pair(
            entity=entity,
            stats=stats,
        )
        best_pairs.append(build_best_pair_data(entity=entity, pair=pair))

    return {"pdb_id": entry["pdb_id"], "pairs": best_pairs}


def build_best_pair_lmdb(
    input_path: Path,
    output_path: Path,
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
            "Path to the raw peptide/receptor LMDB directory "
            "(default: data/pdb_mldata.lmdb)"
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
            limit=args.limit,
        )
    except ValueError as exc:
        parser.error(str(exc))

    stats = build_best_pair_lmdb(
        input_path=args.input_path,
        output_path=args.output_path,
        limit=args.limit,
    )
    print_stats(stats=stats, output_path=args.output_path)


if __name__ == "__main__":
    main()
