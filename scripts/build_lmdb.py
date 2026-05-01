"""Build the raw peptide/receptor pair LMDB from downloaded assemblies.

Goal:
- Read biological assemblies from `data/assemblies.zip`.
- Parse each assembly with Gemmi.
- Write peptide entity records and chain-level peptide/receptor pairs to
  `data/pdb_mldata.lmdb`.

Selection and parsing rules:
- Use mmCIF label subchains end-to-end; do not merge author chains.
- Skip multi-model entries, including NMR entries.
- Trim common terminal caps from protein entity and observed chain sequences.
- Keep trimmed protein entities as peptide entities when their normalized
  entity sequence passes the peptide sequence filter.
- Keep observed peptide chains only when their normalized sequence passes the
  peptide sequence filter after cap trimming.
- Keep internal non-standard amino acids.
- Treat only `ACDEFGHIKLMNPQRSTVWY` as standard one-letter amino-acid codes
  for peptide sequence filtering.
- Reject peptide entity and observed chain sequences with fewer than 4 standard
  amino-acid residues.
- Reject peptide entity and observed chain sequences with more than 20 percent
  non-standard one-letter codes.
- Keep a pair only when at least 4 peptide residues and at least half of peptide
  residues have qualifying contact atoms to the receptor.
- A contact atom must be within 5 Angstroms of the receptor and have B-factor
  less than or equal to 70.
- Normalize sequences only through Gemmi parser APIs.
- Do not manually parse mmCIF parent fields or infer parent amino acids.
- Normalize residues Gemmi cannot convert to a one-letter amino-acid code to `X`.
- Store exact retained 3-letter residue names beside normalized sequences.
- For each peptide chain, identify other chains within 5 Angstroms using
  `scipy.spatial.KDTree`.
- Water and explicitly defined non-polymer ligand chains do not count as
  receptor neighbors.
- Save a pair only when the peptide chain has exactly one meaningful
  neighboring chain and that neighbor is a protein chain.
- Keep a protein neighbor as a receptor only when its normalized sequence has
  at least the configured minimum receptor length after cap trimming.
- If an entry produces duplicate peptide-chain/receptor-chain pair keys, skip
  the whole entry as malformed.
- Do not collapse duplicate pairs to force an entry into the LMDB.

Structure rules:
- Entity sequences represent ideal PDB entity sequences.
- Chain sequences and structures represent observed assembly chains.
- Structures are trimmed with their observed chain sequences.
- Structure arrays store 37 AlphaFold atom positions; extra atoms on
  non-standard amino acids are discarded.

Output behavior:
- Delete an existing output LMDB folder before writing to avoid LMDB upsert
  issues.
- See `docs/storage_schemas.md` for the raw LMDB schema.

Default parameters:
- Assemblies ZIP: `data/assemblies.zip`.
- Output LMDB: `data/pdb_mldata.lmdb`.
- Minimum peptide length: 4.
- Maximum peptide length: 32.
- Minimum standard peptide residues: 4.
- Maximum non-standard peptide residue fraction: 0.2.
- Minimum receptor length: 50.
- Distance threshold: 5.0 Angstroms.
- Minimum contacting peptide residues: 4.
- Minimum contacting peptide residue fraction: 0.5.
- Maximum contact peptide-atom B-factor: 70.0.
- Optional `--limit` for smoke verification.
"""

from __future__ import annotations

import argparse
import shutil
import zipfile
from pathlib import Path

import lmdb
from tqdm import tqdm

from pdb_mldata.curation import BindingFilter
from pdb_mldata.filtering_rules import (
    DuplicatePairError,
    PeptideSequenceFilter,
    process_assembly,
)
from pdb_mldata.lmdb_utils import encode_lmdb_entry

DEFAULT_ASSEMBLIES_ZIP = Path("data/assemblies.zip")
DEFAULT_LMDB_PATH = Path("data/pdb_mldata.lmdb")
DEFAULT_MIN_LENGTH = 4
DEFAULT_MAX_LENGTH = 32
DEFAULT_MIN_STANDARD_PEPTIDE_RESIDUES = 4
DEFAULT_MAX_NONSTANDARD_PEPTIDE_FRACTION = 0.2
DEFAULT_MIN_RECEPTOR_LENGTH = 50
DEFAULT_DISTANCE = 5.0
DEFAULT_MIN_CONTACT_RESIDUES = 4
DEFAULT_MIN_CONTACT_FRACTION = 0.5
DEFAULT_MAX_CONTACT_ATOM_B_FACTOR = 70.0


def build_lmdb(
    zip_path: Path,
    output_path: Path,
    peptide_filter: PeptideSequenceFilter,
    min_receptor_len: int,
    distance: float,
    binding_filter: BindingFilter,
    limit: int | None,
) -> None:
    if output_path.exists():
        shutil.rmtree(output_path)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    env = lmdb.open(str(output_path), map_size=10**11)

    total_pairs = 0
    total_entities = 0
    duplicate_pair_skips = 0

    with zipfile.ZipFile(zip_path, "r") as zf:
        files = [name for name in zf.namelist() if name.endswith(".cif.gz")]
        if limit is not None:
            files = files[:limit]

        progress = tqdm(files, desc="Building LMDB")

        for filename in progress:
            try:
                result = process_assembly(
                    zf,
                    filename,
                    peptide_filter,
                    min_receptor_len,
                    distance,
                    binding_filter,
                )
            except DuplicatePairError as exc:
                duplicate_pair_skips += 1
                progress.write(f"[WARNING] Skipping {exc.pdb_id}: {exc}")
                progress.set_postfix(
                    ents=total_entities,
                    pairs=total_pairs,
                    dup_pairs=duplicate_pair_skips,
                )
                continue

            if result is None:
                progress.set_postfix(
                    ents=total_entities,
                    pairs=total_pairs,
                    dup_pairs=duplicate_pair_skips,
                )
                continue

            pdb_id, data = result
            total_entities += len(data["entities"])
            for entity in data["entities"]:
                total_pairs += len(entity["pairs"])
            progress.set_postfix(
                ents=total_entities,
                pairs=total_pairs,
                dup_pairs=duplicate_pair_skips,
            )

            with env.begin(write=True) as txn:
                txn.put(pdb_id.encode(), encode_lmdb_entry(data))

    env.close()
    print(f"Duplicate-pair entries skipped: {duplicate_pair_skips}")


def validate_parameters(
    zip_path: Path,
    output_path: Path,
    min_length: int,
    max_length: int,
    min_standard_peptide_residues: int,
    max_nonstandard_peptide_fraction: float,
    min_receptor_length: int,
    distance: float,
    min_contact_residues: int,
    min_contact_fraction: float,
    max_contact_atom_b_factor: float,
    limit: int | None,
) -> None:
    """Validate CLI parameters before destructive LMDB creation starts."""
    if not zip_path.exists():
        raise ValueError(f"{zip_path} not found")
    if zip_path.is_dir():
        raise ValueError("zip_path must be a file path, not a directory")
    if output_path.exists() and output_path.is_file():
        raise ValueError("output_path must be an LMDB directory path, not a file")
    if min_length < 1:
        raise ValueError("--min-length must be at least 1")
    if max_length < min_length:
        raise ValueError("--max-length must be greater than or equal to --min-length")
    if min_standard_peptide_residues < 1:
        raise ValueError("--min-standard-peptide-residues must be at least 1")
    if min_standard_peptide_residues > max_length:
        raise ValueError(
            "--min-standard-peptide-residues must be less than or equal to --max-length"
        )
    if not (0 <= max_nonstandard_peptide_fraction <= 1):
        raise ValueError(
            "--max-nonstandard-peptide-fraction must be at least 0 and at most 1"
        )
    if min_receptor_length < 1:
        raise ValueError("--min-receptor-length must be at least 1")
    if max_length > min_receptor_length:
        raise ValueError(
            "--max-length must be less than or equal to --min-receptor-length"
        )
    if distance <= 0:
        raise ValueError("--distance must be greater than 0")
    if min_contact_residues < 1:
        raise ValueError("--min-contact-residues must be at least 1")
    if not (0 <= min_contact_fraction <= 1):
        raise ValueError("--min-contact-fraction must be at least 0 and at most 1")
    if max_contact_atom_b_factor < 0:
        raise ValueError("--max-contact-atom-b-factor must be at least 0")
    if limit is not None and limit < 1:
        raise ValueError("--limit must be at least 1 when provided")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build the peptide/receptor LMDB from assembly ZIP data."
    )
    parser.add_argument(
        "zip_path",
        nargs="?",
        type=Path,
        default=DEFAULT_ASSEMBLIES_ZIP,
        help="Path to the input assemblies ZIP (default: data/assemblies.zip)",
    )
    parser.add_argument(
        "output_path",
        nargs="?",
        type=Path,
        default=DEFAULT_LMDB_PATH,
        help="Path to the output LMDB directory (default: data/pdb_mldata.lmdb)",
    )
    parser.add_argument(
        "--min-length",
        type=int,
        default=DEFAULT_MIN_LENGTH,
        help="Minimum peptide length (default: 4)",
    )
    parser.add_argument(
        "--max-length",
        type=int,
        default=DEFAULT_MAX_LENGTH,
        help="Maximum peptide length (default: 32)",
    )
    parser.add_argument(
        "--min-standard-peptide-residues",
        type=int,
        default=DEFAULT_MIN_STANDARD_PEPTIDE_RESIDUES,
        help="Minimum standard amino acids in peptide sequences (default: 4)",
    )
    parser.add_argument(
        "--max-nonstandard-peptide-fraction",
        type=float,
        default=DEFAULT_MAX_NONSTANDARD_PEPTIDE_FRACTION,
        help="Maximum non-standard peptide residue fraction (default: 0.2)",
    )
    parser.add_argument(
        "--min-receptor-length",
        type=int,
        default=DEFAULT_MIN_RECEPTOR_LENGTH,
        help="Minimum receptor length after cap trimming (default: 50)",
    )
    parser.add_argument(
        "--distance",
        type=float,
        default=DEFAULT_DISTANCE,
        help="Neighbor distance in Angstroms (default: 5.0)",
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
        "--limit",
        type=int,
        default=None,
        help="Process only the first N entries for testing",
    )
    args = parser.parse_args()

    try:
        validate_parameters(
            zip_path=args.zip_path,
            output_path=args.output_path,
            min_length=args.min_length,
            max_length=args.max_length,
            min_standard_peptide_residues=args.min_standard_peptide_residues,
            max_nonstandard_peptide_fraction=args.max_nonstandard_peptide_fraction,
            min_receptor_length=args.min_receptor_length,
            distance=args.distance,
            min_contact_residues=args.min_contact_residues,
            min_contact_fraction=args.min_contact_fraction,
            max_contact_atom_b_factor=args.max_contact_atom_b_factor,
            limit=args.limit,
        )
    except ValueError as exc:
        parser.error(str(exc))

    build_lmdb(
        zip_path=args.zip_path,
        output_path=args.output_path,
        peptide_filter=PeptideSequenceFilter(
            min_length=args.min_length,
            max_length=args.max_length,
            min_standard_residues=args.min_standard_peptide_residues,
            max_nonstandard_fraction=args.max_nonstandard_peptide_fraction,
        ),
        min_receptor_len=args.min_receptor_length,
        distance=args.distance,
        binding_filter=BindingFilter(
            distance=args.distance,
            min_contact_residues=args.min_contact_residues,
            min_contact_fraction=args.min_contact_fraction,
            max_contact_atom_b_factor=args.max_contact_atom_b_factor,
        ),
        limit=args.limit,
    )


if __name__ == "__main__":
    main()
