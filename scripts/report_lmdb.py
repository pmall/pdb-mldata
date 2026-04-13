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
STANDARD_AMINO_ACID_RESIDUES = frozenset(
    {
        "ALA",
        "ARG",
        "ASN",
        "ASP",
        "CYS",
        "GLN",
        "GLU",
        "GLY",
        "HIS",
        "ILE",
        "LEU",
        "LYS",
        "MET",
        "PHE",
        "PRO",
        "SER",
        "THR",
        "TRP",
        "TYR",
        "VAL",
    }
)


@dataclass(frozen=True)
class LmdbCounts:
    pdb_entries: int
    peptide_entities: int
    pairs: int
    unique_peptide_entity_sequences: int
    unique_standard_peptide_entity_residue_name_sequences: int
    unique_peptide_entity_residue_name_sequences: int
    unique_pair_peptide_chain_sequences: int
    unique_standard_pair_peptide_chain_residue_name_sequences: int
    unique_pair_peptide_chain_residue_name_sequences: int


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
    standard_peptide_entity_residue_name_sequences: set[tuple[str, ...]] = set()
    peptide_entity_residue_name_sequences: set[tuple[str, ...]] = set()
    pair_peptide_chain_sequences: set[str] = set()
    standard_pair_peptide_chain_residue_name_sequences: set[tuple[str, ...]] = set()
    pair_peptide_chain_residue_name_sequences: set[tuple[str, ...]] = set()

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
                entity_residue_name_sequence = tuple(entity["residue_names"])
                if all(
                    residue_name in STANDARD_AMINO_ACID_RESIDUES
                    for residue_name in entity_residue_name_sequence
                ):
                    standard_peptide_entity_residue_name_sequences.add(
                        entity_residue_name_sequence
                    )
                peptide_entity_residue_name_sequences.add(entity_residue_name_sequence)
                for pair in entity["pairs"]:
                    pair_peptide_chain_sequences.add(pair["peptide"]["sequence"])
                    pair_residue_name_sequence = tuple(pair["peptide"]["residue_names"])
                    if all(
                        residue_name in STANDARD_AMINO_ACID_RESIDUES
                        for residue_name in pair_residue_name_sequence
                    ):
                        standard_pair_peptide_chain_residue_name_sequences.add(
                            pair_residue_name_sequence
                        )
                    pair_peptide_chain_residue_name_sequences.add(
                        pair_residue_name_sequence
                    )
    env.close()

    return LmdbCounts(
        pdb_entries=pdb_entries,
        peptide_entities=peptide_entities,
        pairs=pairs,
        unique_peptide_entity_sequences=len(peptide_entity_sequences),
        unique_standard_peptide_entity_residue_name_sequences=len(
            standard_peptide_entity_residue_name_sequences
        ),
        unique_peptide_entity_residue_name_sequences=len(
            peptide_entity_residue_name_sequences
        ),
        unique_pair_peptide_chain_sequences=len(pair_peptide_chain_sequences),
        unique_standard_pair_peptide_chain_residue_name_sequences=len(
            standard_pair_peptide_chain_residue_name_sequences
        ),
        unique_pair_peptide_chain_residue_name_sequences=len(
            pair_peptide_chain_residue_name_sequences
        ),
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
            "These counts describe peptide entity sequences after cap trimming.",
            "",
            "| Metric | Meaning | Count |",
            "| --- | --- | ---: |",
            f"| Unique Gemmi one-letter sequences | Unique peptide entity sequences after Gemmi one-letter translation. | {counts.unique_peptide_entity_sequences} |",
            f"| Unique standard 3-letter residue-name sequences | Unique exact peptide entity residue-name sequences containing only the 20 standard amino acids. | {counts.unique_standard_peptide_entity_residue_name_sequences} |",
            f"| Unique 3-letter residue-name sequences | Unique exact peptide entity residue-name sequences, including non-standard residue names. | {counts.unique_peptide_entity_residue_name_sequences} |",
            "",
            "## Peptide-Chain Sequence Diversity",
            "",
            "These counts describe observed peptide chains in saved peptide/receptor pairs after cap trimming.",
            "",
            "| Metric | Meaning | Count |",
            "| --- | --- | ---: |",
            f"| Unique Gemmi one-letter sequences | Unique observed peptide-chain sequences after Gemmi one-letter translation. | {counts.unique_pair_peptide_chain_sequences} |",
            f"| Unique standard 3-letter residue-name sequences | Unique exact observed peptide-chain residue-name sequences containing only the 20 standard amino acids. | {counts.unique_standard_pair_peptide_chain_residue_name_sequences} |",
            f"| Unique 3-letter residue-name sequences | Unique exact observed peptide-chain residue-name sequences, including non-standard residue names. | {counts.unique_pair_peptide_chain_residue_name_sequences} |",
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
    print(
        f"Unique entity one-letter sequences: {counts.unique_peptide_entity_sequences}"
    )
    print(
        "Unique standard entity residue-name sequences: "
        f"{counts.unique_standard_peptide_entity_residue_name_sequences}"
    )
    print(
        "Unique entity residue-name sequences: "
        f"{counts.unique_peptide_entity_residue_name_sequences}"
    )
    print(
        "Unique peptide-chain one-letter sequences in saved pairs: "
        f"{counts.unique_pair_peptide_chain_sequences}"
    )
    print(
        "Unique standard peptide-chain residue-name sequences in saved pairs: "
        f"{counts.unique_standard_pair_peptide_chain_residue_name_sequences}"
    )
    print(
        "Unique peptide-chain residue-name sequences in saved pairs: "
        f"{counts.unique_pair_peptide_chain_residue_name_sequences}"
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
