from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import cast

import lmdb
from tqdm import tqdm

from pdb_mldata.lmdb_utils import decode_lmdb_entry

DEFAULT_LMDB_PATH = Path("data/pdb_mldata.lmdb")
DEFAULT_REPORT_PATH = Path("data/lmdb_report.md")


@dataclass(frozen=True)
class LmdbCounts:
    pdb_entries: int
    peptide_entities: int
    pairs: int
    unique_peptide_entity_sequences: int
    unique_peptide_entity_residue_name_sequences: int
    unique_pair_peptide_chain_sequences: int


def validate_parameters(db_path: Path, report_path: Path) -> None:
    """Validate CLI parameters before opening the LMDB environment."""
    if not db_path.exists():
        raise ValueError(f"{db_path} not found")
    if not db_path.is_dir():
        raise ValueError("db_path must be an LMDB directory path")
    if report_path.exists() and report_path.is_dir():
        raise ValueError("report_path must be a file path, not a directory")


def collect_lmdb_counts(db_path: Path) -> LmdbCounts:
    pdb_entries = 0
    peptide_entities = 0
    pairs = 0
    peptide_entity_sequences: set[str] = set()
    peptide_entity_residue_name_sequences: set[tuple[str, ...]] = set()
    pair_peptide_chain_sequences: set[str] = set()

    env = lmdb.open(str(db_path), readonly=True, lock=False)
    with env.begin() as txn:
        total_entries = txn.stat()["entries"]
        cursor = txn.cursor()
        for _key, value in tqdm(cursor, total=total_entries, desc="Reporting LMDB"):
            entry = decode_lmdb_entry(cast(bytes, value))
            pdb_entries += 1
            peptide_entities += len(entry["entities"])
            for entity in entry["entities"]:
                pairs += len(entity["pairs"])
                peptide_entity_sequences.add(entity["sequence"])
                peptide_entity_residue_name_sequences.add(
                    tuple(entity["residue_names"])
                )
                for pair in entity["pairs"]:
                    pair_peptide_chain_sequences.add(pair["peptide"]["sequence"])
    env.close()

    return LmdbCounts(
        pdb_entries=pdb_entries,
        peptide_entities=peptide_entities,
        pairs=pairs,
        unique_peptide_entity_sequences=len(peptide_entity_sequences),
        unique_peptide_entity_residue_name_sequences=len(
            peptide_entity_residue_name_sequences
        ),
        unique_pair_peptide_chain_sequences=len(pair_peptide_chain_sequences),
    )


def render_markdown_report(counts: LmdbCounts, db_path: Path) -> str:
    return "\n".join(
        [
            "# PDB MLData LMDB Report",
            "",
            f"Source LMDB: `{db_path}`",
            "",
            "## Dataset Size",
            "",
            "The LMDB hierarchy is `PDB entry -> peptide entity -> peptide/receptor pair`.",
            "",
            "| Level | Meaning | Count |",
            "| --- | --- | ---: |",
            f"| PDB entries | LMDB records keyed by PDB ID. | {counts.pdb_entries} |",
            f"| Peptide entities | PDB entry/entity records saved after filtering. | {counts.peptide_entities} |",
            f"| Pairs | Observed peptide-chain/receptor-chain pairs saved under those entities. | {counts.pairs} |",
            "",
            "## Peptide Entity Sequence Diversity",
            "",
            "These counts describe peptide entity sequences after cap trimming and Gemmi one-letter conversion.",
            "",
            "| Metric | Meaning | Count |",
            "| --- | --- | ---: |",
            f"| Unique peptide entity sequences | Unique one-letter peptide entity sequences. | {counts.unique_peptide_entity_sequences} |",
            f"| Unique exact peptide entity residue-name sequences | Unique peptide entity sequences represented by exact residue names, including non-standard residue names. | {counts.unique_peptide_entity_residue_name_sequences} |",
            "",
            "## Peptide Chains In Saved Pairs",
            "",
            "A peptide entity can have more than one saved peptide-chain/receptor-chain pair. This section counts the peptide chain side of those pairs.",
            "",
            "| Metric | Meaning | Count |",
            "| --- | --- | ---: |",
            f"| Unique peptide-chain sequences in saved pairs | Unique one-letter sequences for peptide chains that appear in saved pairs. | {counts.unique_pair_peptide_chain_sequences} |",
            "",
        ]
    )


def write_lmdb_report(db_path: Path, report_path: Path) -> None:
    counts = collect_lmdb_counts(db_path)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(render_markdown_report(counts, db_path), encoding="utf-8")
    print(f"Wrote LMDB report to {report_path}")
    print(f"PDB entries: {counts.pdb_entries}")
    print(f"Peptide entities: {counts.peptide_entities}")
    print(f"Pairs: {counts.pairs}")
    print(f"Unique peptide entity sequences: {counts.unique_peptide_entity_sequences}")
    print(
        "Unique peptide entity residue-name sequences: "
        f"{counts.unique_peptide_entity_residue_name_sequences}"
    )
    print(
        "Unique peptide-chain sequences in saved pairs: "
        f"{counts.unique_pair_peptide_chain_sequences}"
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Report basic counts from the peptide/receptor LMDB."
    )
    parser.add_argument(
        "db_path",
        nargs="?",
        type=Path,
        default=DEFAULT_LMDB_PATH,
        help="Path to the LMDB database directory (default: data/pdb_mldata.lmdb)",
    )
    parser.add_argument(
        "report_path",
        nargs="?",
        type=Path,
        default=DEFAULT_REPORT_PATH,
        help="Path to the Markdown report file (default: data/lmdb_report.md)",
    )
    args = parser.parse_args()

    try:
        validate_parameters(db_path=args.db_path, report_path=args.report_path)
    except ValueError as exc:
        parser.error(str(exc))

    write_lmdb_report(db_path=args.db_path, report_path=args.report_path)


if __name__ == "__main__":
    main()
