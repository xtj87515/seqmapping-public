#!/usr/bin/env python3
"""
Remove UniParc sequences (UniRef90_UPI*) and append unique SwissProt sequences matched by UniProt accession.
Also outputs a list of appended UniProt IDs with sequence lengths.

SwissProt headers like:
    >sp|Q6GZX3|002L_FRG3G ...
become:
    >Q6GZX3

Uniref90 headers like:
  >UniRef90_A0A5A9P0L4 ...
become:
  >A0A5A9P0L4 ...

NOTE:
  - We keep the full UniRef90 identifier token after "UniRef90_" up to whitespace.
    This avoids collapsing distinct IDs such as Q8WZ42, Q8WZ42-2, Q8WZ42.1, Q8WZ42_foo
    into a single header (which caused duplicates like multiple ">Q8WZ42" records in preprocessed_uniref90.fasta).
  - For SwissProt appending, we still track canonical UniProt accessions (split at - . _)
    so we do not append SwissProt records that are already represented by UniRef entries.
"""

from __future__ import annotations

import argparse
import logging
import re
import time
from dataclasses import dataclass
from pathlib import Path

from seqmapping.utils.logging import get_logger
from seqmapping.utils.paths import (
    INTERMEDIATE_DIR,
    UNIREF90_FASTA,
    UNIPROT_SPROT_FASTA,
    UNIREF90_PROCESSED_FASTA,
    ADDED_SPROT_IDS_TSV,
    ensure_directories,
)

# Regex patterns
# UniRef90: capture the FULL token after UniRef90_ up to first whitespace (Option A)
UNIREF90_ID_PATTERN = re.compile(r"^>UniRef90_([^\s]+)")
SPROT_ACC_PATTERN = re.compile(r"^>sp\|([A-Z0-9]+)\|")


@dataclass(frozen=True)
class RemoveStats:
    removed: int
    kept: int
    # Full UniRef IDs written (e.g., Q8WZ42-2)
    seen_uniref_ids: set[str]
    # Canonical UniProt accessions observed (e.g., Q8WZ42)
    seen_uniprot_accessions: set[str]


@dataclass(frozen=True)
class AppendStats:
    added: int


def canonical_uniprot_acc(token: str) -> str:
    """
    Convert a UniRef token to a canonical UniProt-like accession for "presence" checks.
    Examples:
      Q8WZ42-2 -> Q8WZ42
      Q8WZ42.1 -> Q8WZ42
      Q8WZ42_foo -> Q8WZ42
    """
    return re.split(r"[-._]", token, maxsplit=1)[0]


def remove_uniparc_seqs(input_fasta: Path, output_fasta: Path, logger: logging.Logger) -> RemoveStats:
    # Remove UniParc (UniRef90_UPI*) sequences and write kept ones; collect IDs.
    removed = 0
    kept = 0
    seen_uniref_ids: set[str] = set()
    seen_uniprot_accessions: set[str] = set()

    skip = False
    write_current = False  # whether to write sequence lines for current record

    logger.info("Processing UniRef90 FASTA: %s", input_fasta)
    output_fasta.parent.mkdir(parents=True, exist_ok=True)

    with input_fasta.open("rt") as fin, output_fasta.open("wt") as fout:
        for line in fin:
            if line.startswith(">"):
                # reset record flags
                skip = False
                write_current = False

                if line.startswith(">UniRef90_UPI"):
                    skip = True
                    removed += 1
                    continue

                # keep
                kept += 1

                # Normalize header: >UniRef90_<TOKEN> -> ><TOKEN>  (full token, no truncation)
                m = UNIREF90_ID_PATTERN.match(line)
                if m:
                    token = m.group(1)

                    # If UniRef90 FASTA contains exact duplicate token headers, skip duplicates
                    if token in seen_uniref_ids:
                        skip = True
                        continue

                    seen_uniref_ids.add(token)
                    seen_uniprot_accessions.add(canonical_uniprot_acc(token))
                    fout.write(f">{token}\n")
                    write_current = True
                else:
                    # Fallback: strip UniRef90_ prefix but keep rest of line as-is
                    fout.write(re.sub(r"^>UniRef90_", ">", line))
                    write_current = True

            elif not skip and write_current:
                fout.write(line)

    logger.info("Removed UniParc sequences: %d", removed)
    logger.info("Kept non-UniParc sequences: %d", kept)
    logger.info("Unique UniRef IDs written: %d", len(seen_uniref_ids))
    logger.info("Unique canonical UniProt accessions observed: %d", len(seen_uniprot_accessions))
    return RemoveStats(
        removed=removed,
        kept=kept,
        seen_uniref_ids=seen_uniref_ids,
        seen_uniprot_accessions=seen_uniprot_accessions,
    )


def _flush_sprot_record(fout, fidout, header_line: str, seq_lines: list[str]) -> int:
    # Writes one SwissProt record to output FASTA and TSV.
    # header_line must already be simplified (e.g., '>Q6GZX3\n').
    # Returns sequence length.

    fout.write(header_line)
    for l in seq_lines:
        fout.write(l)

    seq_len = sum(len(l.strip()) for l in seq_lines)
    acc = header_line[1:].strip()
    fidout.write(f"{acc}\t{seq_len}\n")
    return seq_len


def append_unique_sprot(
    sprot_fasta: Path,
    output_fasta: Path,
    added_ids_tsv: Path,
    seen_uniprot_accessions: set[str],
    logger: logging.Logger,
) -> AppendStats:

    # Append only SwissProt sequences whose accession is not already present in UniRef90.
    # Simplify headers to just '>ACC'.
    # Record added accessions and sequence lengths in TSV file.

    logger.info("Appending unique SwissProt sequences from: %s", sprot_fasta)

    if not sprot_fasta.exists():
        logger.warning("SwissProt FASTA not found, skipping: %s", sprot_fasta)
        added_ids_tsv.parent.mkdir(parents=True, exist_ok=True)
        added_ids_tsv.write_text("sseqid\tsseq_length\n")
        return AppendStats(added=0)

    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    added_ids_tsv.parent.mkdir(parents=True, exist_ok=True)

    added = 0
    current_acc: str | None = None
    current_seq_lines: list[str] = []
    keep_current = False

    with sprot_fasta.open("rt") as fin, output_fasta.open("at") as fout, added_ids_tsv.open("wt") as fidout:
        fidout.write("sseqid\tsseq_length\n")

        for line in fin:
            if line.startswith(">"):
                # flush previous record if we were keeping it
                if keep_current and current_acc is not None:
                    _flush_sprot_record(
                        fout=fout,
                        fidout=fidout,
                        header_line=f">{current_acc}\n",
                        seq_lines=current_seq_lines,
                    )
                    added += 1

                # reset for new record
                current_seq_lines = []
                current_acc = None
                keep_current = False

                m = SPROT_ACC_PATTERN.match(line)
                if not m:
                    continue  # skip malformed header

                acc = m.group(1)

                # Skip if UniRef90 already contains this canonical accession
                if acc in seen_uniprot_accessions:
                    continue

                current_acc = acc
                keep_current = True

            else:
                if keep_current:
                    current_seq_lines.append(line)

        # flush last record
        if keep_current and current_acc is not None:
            _flush_sprot_record(
                fout=fout,
                fidout=fidout,
                header_line=f">{current_acc}\n",
                seq_lines=current_seq_lines,
            )
            added += 1

    logger.info("Appended SwissProt unique sequences: %d", added)
    logger.info("Added IDs + lengths written to: %s", added_ids_tsv)
    return AppendStats(added=added)


def build_argparser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Remove UniParc sequences from UniRef90 FASTA and append missing SwissProt sequences."
    )

    # - raw inputs in data/raw
    # - outputs in data/intermediate
    p.add_argument(
        "--uniref90",
        type=Path,
        default=UNIREF90_FASTA,
        help="Input UniRef90 FASTA (default: data/raw/uniref90/uniref90.fasta)",
    )
    p.add_argument(
        "--sprot",
        type=Path,
        default=UNIPROT_SPROT_FASTA,
        help="Input UniProt Swiss-Prot FASTA (default: data/raw/uniprot/uniprot_sprot.fasta)",
    )
    p.add_argument(
        "--out-fasta",
        type=Path,
        default=UNIREF90_PROCESSED_FASTA,
        help="Output processed FASTA (default: data/intermediate/uniref90_processed.fasta)",
    )
    p.add_argument(
        "--added-ids",
        type=Path,
        default=ADDED_SPROT_IDS_TSV,
        help="Output TSV of added IDs + lengths (default: data/intermediate/added_sprot_ids.tsv)",
    )
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_argparser().parse_args(argv)

    ensure_directories()

    # Make sure log file name matches the src filename (e.g., preprocess_uniref90.log)
    module_stem = Path(__file__).stem  # preprocess_uniref90
    logger = get_logger(__name__, module_stem)

    start = time.time()
    logger.info("=== preprocess_uniref90 started ===")
    logger.info("UniRef90 input: %s", args.uniref90)
    logger.info("SwissProt input: %s", args.sprot)
    logger.info("Processed output FASTA: %s", args.out_fasta)
    logger.info("Added IDs TSV: %s", args.added_ids)

    if not args.uniref90.exists():
        logger.error("UniRef90 input FASTA not found: %s", args.uniref90)
        return 2

    try:
        stats = remove_uniparc_seqs(args.uniref90, args.out_fasta, logger)
        sprot_stats = append_unique_sprot(
            args.sprot,
            args.out_fasta,
            args.added_ids,
            stats.seen_uniprot_accessions,
            logger,
        )

        logger.info(
            "Removed: %d | Kept: %d | Appended SwissProt: %d",
            stats.removed,
            stats.kept,
            sprot_stats.added,
        )
        logger.info("Total duration: %.2f seconds", time.time() - start)
        logger.info("=== preprocess_uniref90 finished ===")
        return 0
    except Exception:
        logger.exception("Unhandled error")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
