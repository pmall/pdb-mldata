"""Filter raw peptide/receptor pairs into the binding-curated LMDB.

Goal:
- Read chain-level pairs from `data/pdb_mldata.lmdb`.
- Reject pairs that do not look like useful peptide/receptor binding samples.
- Write accepted pairs to `data/pdb_mldata_binding.lmdb` with the same schema.

Filtering rules:
- Apply peptide-chain content filtering before distance/contact filtering.
- Treat only `ACDEFGHIKLMNPQRSTVWY` as standard one-letter amino-acid codes.
- Reject pairs whose peptide chain has fewer than 4 standard amino-acid residues.
- Reject pairs whose peptide chain has more than 20 percent non-standard
  one-letter codes.
- Count `X`, `U`, and every other code outside the standard alphabet as
  non-standard.
- Use all finite atom coordinates from the stored 37-atom arrays for contact
  evaluation.
- A peptide residue counts as contacting the receptor when any stored peptide
  atom is within 5 Angstroms of any stored receptor atom and that same peptide
  atom has B-factor less than or equal to 70.
- Keep a pair only when at least 4 peptide residues and at least half of peptide
  residues have qualifying contact atoms.
- Do not use occupancy in this curation rule.
- Reject pairs with no usable peptide or receptor coordinates.
- Drop peptide entities and entries with no accepted pairs.
- Keep all accepted pairs; this script does not select the best chain pair.

Output behavior:
- Delete an existing output LMDB folder before writing to avoid LMDB upsert
  issues.
- See `docs/storage_schemas.md` for the raw/binding LMDB schema.

Default parameters:
- Input LMDB: `data/pdb_mldata.lmdb`.
- Output LMDB: `data/pdb_mldata_binding.lmdb`.
- Distance threshold: 5.0 Angstroms.
- Minimum contacting peptide residues: 4.
- Minimum contacting peptide residue fraction: 0.5.
- Maximum contact peptide-atom B-factor: 70.0.
- Minimum standard peptide-chain residues: 4.
- Maximum non-standard peptide-chain residue fraction: 0.2.
- Optional `--limit` for smoke verification.
"""

from __future__ import annotations

import argparse
import shutil
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import cast

import lmdb
import msgpack
from tqdm import tqdm

from pdb_mldata.curation import (
    BindingFilter,
    PeptideContentFilter,
    evaluate_binding_pair,
    evaluate_peptide_content,
)
from pdb_mldata.lmdb_utils import (
    EntityData,
    LmdbEntry,
    PairData,
    decode_chain_data,
    encode_lmdb_entry,
)

DEFAULT_INPUT_LMDB_PATH = Path("data/pdb_mldata.lmdb")
DEFAULT_OUTPUT_LMDB_PATH = Path("data/pdb_mldata_binding.lmdb")
DEFAULT_DISTANCE = 5.0
DEFAULT_MIN_CONTACT_RESIDUES = 4
DEFAULT_MIN_CONTACT_FRACTION = 0.5
DEFAULT_MAX_CONTACT_ATOM_B_FACTOR = 70.0
DEFAULT_MIN_STANDARD_PEPTIDE_RESIDUES = 4
DEFAULT_MAX_NONSTANDARD_PEPTIDE_FRACTION = 0.2


@dataclass
class BindingFilterStats:
    entries_read: int = 0
    entities_read: int = 0
    pairs_read: int = 0
    entries_written: int = 0
    entities_written: int = 0
    pairs_written: int = 0
    entities_dropped: int = 0
    entries_dropped: int = 0
    rejected_pairs: int = 0
    rejection_reasons: Counter[str] = field(default_factory=Counter)


def decode_pair_for_metrics(pair: PairData) -> PairData:
    return {
        "peptide": decode_chain_data(pair["peptide"]),
        "receptor": decode_chain_data(pair["receptor"]),
    }


def filter_entity_pairs(
    entity: EntityData,
    peptide_content_filter: PeptideContentFilter,
    binding_filter: BindingFilter,
    stats: BindingFilterStats,
) -> EntityData | None:
    filtered_pairs: list[PairData] = []

    for pair in entity["pairs"]:
        stats.pairs_read += 1
        peptide_content_decision = evaluate_peptide_content(
            peptide_sequence=pair["peptide"]["sequence"],
            peptide_content_filter=peptide_content_filter,
        )
        if not peptide_content_decision.is_accepted:
            stats.rejected_pairs += 1
            stats.rejection_reasons[peptide_content_decision.reason] += 1
            continue

        decoded_pair = decode_pair_for_metrics(pair)
        decision = evaluate_binding_pair(
            pair=decoded_pair,
            binding_filter=binding_filter,
        )
        if decision.is_accepted:
            filtered_pairs.append(pair)
            continue

        stats.rejected_pairs += 1
        stats.rejection_reasons[decision.reason] += 1

    if not filtered_pairs:
        stats.entities_dropped += 1
        return None

    return {
        "entity_id": entity["entity_id"],
        "sequence": entity["sequence"],
        "residue_names": entity["residue_names"],
        "pairs": filtered_pairs,
    }


def filter_lmdb_entry(
    entry: LmdbEntry,
    peptide_content_filter: PeptideContentFilter,
    binding_filter: BindingFilter,
    stats: BindingFilterStats,
) -> LmdbEntry | None:
    filtered_entities: list[EntityData] = []

    stats.entries_read += 1
    stats.entities_read += len(entry["entities"])
    for entity in entry["entities"]:
        filtered_entity = filter_entity_pairs(
            entity=entity,
            peptide_content_filter=peptide_content_filter,
            binding_filter=binding_filter,
            stats=stats,
        )
        if filtered_entity is not None:
            filtered_entities.append(filtered_entity)

    if not filtered_entities:
        stats.entries_dropped += 1
        return None

    return {
        "pdb_id": entry["pdb_id"],
        "entities": filtered_entities,
    }


def unpack_lmdb_entry(entry_bytes: bytes) -> LmdbEntry:
    return cast(LmdbEntry, msgpack.unpackb(entry_bytes, raw=False))


def build_binding_filtered_lmdb(
    input_path: Path,
    output_path: Path,
    peptide_content_filter: PeptideContentFilter,
    binding_filter: BindingFilter,
    limit: int | None,
) -> BindingFilterStats:
    if output_path.exists():
        shutil.rmtree(output_path)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    input_env = lmdb.open(str(input_path), readonly=True, lock=False)
    output_env = lmdb.open(str(output_path), map_size=10**11)
    stats = BindingFilterStats()

    with input_env.begin() as input_txn:
        total_entries = input_txn.stat()["entries"]
        if limit is not None:
            total_entries = min(total_entries, limit)

        cursor = input_txn.cursor()
        progress = tqdm(cursor, total=total_entries, desc="Filtering binding pairs")
        for entry_index, (key, value) in enumerate(progress):
            if limit is not None and entry_index >= limit:
                break

            entry = unpack_lmdb_entry(cast(bytes, value))
            filtered_entry = filter_lmdb_entry(
                entry=entry,
                peptide_content_filter=peptide_content_filter,
                binding_filter=binding_filter,
                stats=stats,
            )
            if filtered_entry is not None:
                stats.entries_written += 1
                stats.entities_written += len(filtered_entry["entities"])
                stats.pairs_written += sum(
                    len(entity["pairs"]) for entity in filtered_entry["entities"]
                )
                with output_env.begin(write=True) as output_txn:
                    output_txn.put(cast(bytes, key), encode_lmdb_entry(filtered_entry))

            progress.set_postfix(
                ents=stats.entities_written,
                pairs=stats.pairs_written,
                rej_pairs=stats.rejected_pairs,
            )

    input_env.close()
    output_env.close()
    return stats


def validate_parameters(
    input_path: Path,
    output_path: Path,
    distance: float,
    min_contact_residues: int,
    min_contact_fraction: float,
    max_contact_atom_b_factor: float,
    min_standard_peptide_residues: int,
    max_nonstandard_peptide_fraction: float,
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
    if min_contact_residues < 1:
        raise ValueError("--min-contact-residues must be at least 1")
    if not (0 < min_contact_fraction <= 1):
        raise ValueError("--min-contact-fraction must be greater than 0 and at most 1")
    if max_contact_atom_b_factor <= 0:
        raise ValueError("--max-contact-atom-b-factor must be greater than 0")
    if min_standard_peptide_residues < 1:
        raise ValueError("--min-standard-peptide-residues must be at least 1")
    if not (0 <= max_nonstandard_peptide_fraction <= 1):
        raise ValueError(
            "--max-nonstandard-peptide-fraction must be at least 0 and at most 1"
        )
    if limit is not None and limit < 1:
        raise ValueError("--limit must be at least 1 when provided")


def print_stats(stats: BindingFilterStats, output_path: Path) -> None:
    print(f"Wrote binding-filtered LMDB to {output_path}")
    print(f"PDB entries: {stats.entries_read} -> {stats.entries_written}")
    print(f"Peptide entities: {stats.entities_read} -> {stats.entities_written}")
    print(f"Peptide/receptor chain pairs: {stats.pairs_read} -> {stats.pairs_written}")
    print(f"Rejected peptide/receptor chain pairs: {stats.rejected_pairs}")
    for reason, count in stats.rejection_reasons.most_common():
        print(f"Rejected {reason}: {count}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Filter LMDB peptide/receptor pairs by subjective binding-contact rules."
    )
    parser.add_argument(
        "input_path",
        nargs="?",
        type=Path,
        default=DEFAULT_INPUT_LMDB_PATH,
        help="Path to the input LMDB directory (default: data/pdb_mldata.lmdb)",
    )
    parser.add_argument(
        "output_path",
        nargs="?",
        type=Path,
        default=DEFAULT_OUTPUT_LMDB_PATH,
        help=(
            "Path to the output binding-filtered LMDB directory "
            "(default: data/pdb_mldata_binding.lmdb)"
        ),
    )
    parser.add_argument(
        "--distance",
        type=float,
        default=DEFAULT_DISTANCE,
        help="Contact distance in Angstroms (default: 5.0)",
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
        "--min-standard-peptide-residues",
        type=int,
        default=DEFAULT_MIN_STANDARD_PEPTIDE_RESIDUES,
        help="Minimum standard amino acids in the peptide chain (default: 4)",
    )
    parser.add_argument(
        "--max-nonstandard-peptide-fraction",
        type=float,
        default=DEFAULT_MAX_NONSTANDARD_PEPTIDE_FRACTION,
        help="Maximum non-standard peptide-chain residue fraction (default: 0.2)",
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
            min_contact_residues=args.min_contact_residues,
            min_contact_fraction=args.min_contact_fraction,
            max_contact_atom_b_factor=args.max_contact_atom_b_factor,
            min_standard_peptide_residues=args.min_standard_peptide_residues,
            max_nonstandard_peptide_fraction=args.max_nonstandard_peptide_fraction,
            limit=args.limit,
        )
    except ValueError as exc:
        parser.error(str(exc))

    stats = build_binding_filtered_lmdb(
        input_path=args.input_path,
        output_path=args.output_path,
        peptide_content_filter=PeptideContentFilter(
            min_standard_residues=args.min_standard_peptide_residues,
            max_nonstandard_fraction=args.max_nonstandard_peptide_fraction,
        ),
        binding_filter=BindingFilter(
            distance=args.distance,
            min_contact_residues=args.min_contact_residues,
            min_contact_fraction=args.min_contact_fraction,
            max_contact_atom_b_factor=args.max_contact_atom_b_factor,
        ),
        limit=args.limit,
    )
    print_stats(stats=stats, output_path=args.output_path)


if __name__ == "__main__":
    main()
