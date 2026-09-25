#!/usr/bin/env python3

"""
Pipeline summary (two-pass DIAMOND search, restart-safe):
1. Initialize workflow and directories
    - Ensure all required directories exist (raw, intermediate, hits, logs).
    - Initialize centralized logging and record run metadata.
2. Load PDB sequence features table
    - Input: data/intermediate/pdb_features.tsv
    - Read columns:
        * qseqid
        * qseq_length
        * qseq_n_source_organisms
        * qseq_tax_name
        * qseq_description
3. Build Pass 2 seed FASTA (filter-based selection)
    - Evaluate all PDB unique sequences against the following rules:
    - A sequence is selected for Pass 2 if ANY condition is met:
        * qseq_length <= MIN_QUERY_LENGTH (default: 30)
        * qseq_n_source_organisms > MAX_ORGANISMS (default: 1) OR missing
        * qseq_tax_name OR qseq_description contains:
            - "Human immunodeficiency" OR "HIV" (case-insensitive)
        * qseq_tax_name OR qseq_description contains:
            - "zinc" OR "thioredoxin" (case-insensitive)
    - Write selected sequences to:
        * data/intermediate/pdb_filtered_seqs_pass2.fasta
    - Create restart marker:
        * data/intermediate/filter_pass2.done
4. Create DIAMOND database (if needed)
    - Input FASTA:
        * data/intermediate/uniref90_processed.fasta
    - Output database:
        * data/intermediate/uniref90_processed.dmnd
    - Run `diamond makedb` only if database is missing or forced.
    - Validate database file and create restart marker:
        * data/intermediate/diamond_makedb.done
5. Run DIAMOND blastp — Pass 1 (all sequences)
    - Input FASTA:
        * data/intermediate/pdb_unique_sequences.fasta
    - Output TSV:
        * data/hits/diamond_hits_pass1.tsv
    - DIAMOND parameters:
        * Fully specified via config.py
        * Includes explicit defaults (e.g. BLOSUM62, masking enabled,
          default e-value, single best HSP per target)
    - Write output atomically (tmp → rename).
    - Validate output and create restart marker:
        * data/intermediate/diamond_pass1.done
6. Append no-hit queries to Pass 2 FASTA (one-time, idempotent)
    - Identify all qseqid values absent from Pass 1 results.
    - Append corresponding FASTA records to:
        * data/intermediate/pdb_filtered_seqs_pass2.fasta
    - Append is performed atomically (rewrite + replace).
    - Write markers:
        * data/intermediate/pdb_nohit.done 
        * data/intermediate/append_nohits.done   (restart control)
7. Run DIAMOND blastp — Pass 2 (filtered + no-hit sequences)
    - Input FASTA:
        * data/intermediate/pdb_filtered_seqs_pass2.fasta
    - Output TSV:
        * data/hits/diamond_hits_pass2.tsv
    - DIAMOND parameters:
        * Very-sensitive mode
        * High e-value threshold
        * PAM30 scoring matrix
        * Masking disabled
        * All HSPs reported
    - Write output atomically (tmp → rename).
    - Validate output and create restart marker:
        * data/intermediate/diamond_pass2.done
8. Workflow management and restart behavior
    - Each major step is guarded by a .done marker.
    - On restart (default behavior):
        * Completed steps are skipped automatically.
        * Outputs are validated before skipping.
    - CLI controls:
        * --resume / --no-resume
        * --force (rerun all steps)
        * --force-steps <step1 step2 ...>
        * --rebuild-pass2-fasta
    - All external command stdout/stderr is captured in step-specific log files.
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
import time
from pathlib import Path
from typing import Iterable
import shutil 

import pandas as pd

from seqmapping.utils.config import (
    DIAMOND_BIN,
    DIAMOND_OUTFMT_FIELDS,
    DIAMOND_PASS1,
    DIAMOND_PASS2,
    DIAMOND_THREADS,
    MAX_ORGANISMS,
    MIN_QUERY_LENGTH,
)
from seqmapping.utils.logging import get_logger, start_run
from seqmapping.utils.paths import (
    # Inputs
    PDB_FEATURES_TSV,
    PDB_UNIQUE_FASTA,
    UNIREF90_PROCESSED_FASTA,
    # Outputs
    DIAMOND_HITS_PASS1_TSV,
    DIAMOND_HITS_PASS2_TSV,
    PDB_FILTERED_SEQS_PASS2_FASTA,
    PDB_NOHIT_MARKER,
    UNIREF90_DIAMOND_DB_PREFIX,
    # Done markers
    FILTER_PASS2_DONE,
    DIAMOND_MAKEDB_DONE,
    DIAMOND_PASS1_DONE,
    APPEND_NOHITS_DONE,
    DIAMOND_PASS2_DONE,
    # Logging / dirs
    LOG_DIR,
    ensure_directories,
)

# Logs
MAKEDB_LOG: Path = LOG_DIR / "run_diamond_makedb.log"
SEARCH1_LOG: Path = LOG_DIR / "run_diamond_pass1.log"
SEARCH2_LOG: Path = LOG_DIR / "run_diamond_pass2.log"


# -----------------------------------------------------------------------------
# Atomic IO helpers + basic validation
# -----------------------------------------------------------------------------
def atomic_write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(text, encoding="utf-8")
    os.replace(tmp, path)


def atomic_copy_replace(src_tmp: Path, dst_final: Path) -> None:
    dst_final.parent.mkdir(parents=True, exist_ok=True)
    os.replace(src_tmp, dst_final)


def is_nonempty_file(path: Path) -> bool:
    try:
        return path.exists() and path.is_file() and path.stat().st_size > 0
    except OSError:
        return False


def looks_like_fasta(path: Path) -> bool:
    # Cheap sanity check: file non-empty and has at least one '>' header.
    if not is_nonempty_file(path):
        return False
    try:
        with path.open("r", encoding="utf-8") as f:
            for _ in range(2000):  # scan first N lines
                line = f.readline()
                if not line:
                    break
                if line.startswith(">"):
                    return True
        return False
    except OSError:
        return False


def looks_like_tsv(path: Path, min_lines: int = 1) -> bool:
    # Cheap sanity check: non-empty and has at least `min_lines` lines.
    if not is_nonempty_file(path):
        return False
    try:
        n = 0
        with path.open("r", encoding="utf-8") as f:
            for line in f:
                if line.strip():
                    n += 1
                if n >= min_lines:
                    return True
        return False
    except OSError:
        return False


def dmnd_path(db_prefix: Path) -> Path:
    return Path(f"{db_prefix}.dmnd")

def require_diamond(logger) -> None:
    """
    Fail fast if DIAMOND is missing/not runnable.
    Log the resolved path and version for reproducibility.
    """
    if DIAMOND_BIN == "diamond":
        resolved = shutil.which("diamond")
        if resolved is None:
            raise RuntimeError(
                "DIAMOND not found on PATH. Install with:\n"
                "  sudo apt-get update && sudo apt-get install -y diamond-aligner\n"
                "or set DIAMOND_BIN in config.py to an absolute path."
            )
        diamond_cmd = "diamond"
    else:
        resolved = DIAMOND_BIN
        diamond_cmd = DIAMOND_BIN
        if not Path(resolved).exists():
            raise RuntimeError(f"DIAMOND_BIN does not exist: {resolved}")

    r = subprocess.run([diamond_cmd, "--version"], capture_output=True, text=True, check=True)
    logger.info(f"Using DIAMOND: {resolved}")
    logger.info(f"DIAMOND version: {r.stdout.strip()}")


# -----------------------------------------------------------------------------
# Step bookkeeping
# -----------------------------------------------------------------------------
STEP_FILTER = "filter_pass2_seed"
STEP_MAKEDB = "makedb"
STEP_PASS1 = "pass1"
STEP_APPEND = "append_nohits"
STEP_PASS2 = "pass2"

ALL_STEPS = [STEP_FILTER, STEP_MAKEDB, STEP_PASS1, STEP_APPEND, STEP_PASS2]

STEP_DONE_FILES = {
    STEP_FILTER: FILTER_PASS2_DONE,
    STEP_MAKEDB: DIAMOND_MAKEDB_DONE,
    STEP_PASS1: DIAMOND_PASS1_DONE,
    STEP_APPEND: APPEND_NOHITS_DONE,
    STEP_PASS2: DIAMOND_PASS2_DONE,
}


def mark_done(path: Path, msg: str) -> None:
    atomic_write_text(path, msg + "\n")


def clear_done(path: Path) -> None:
    if path.exists():
        path.unlink()


def should_run_step(step: str, args, logger) -> bool:
    """
    Resume logic:
      - If --force: run everything
      - If step listed in --force-steps: run it
      - Else if --resume and .done exists: skip
      - Else run
    """
    if args.force:
        return True
    if step in args.force_steps:
        return True
    done_file = STEP_DONE_FILES[step]
    if args.resume and done_file.exists():
        logger.info(f"Skipping {step}: done marker exists ({done_file})")
        return False
    return True


# -----------------------------------------------------------------------------
# DIAMOND command build/run
# -----------------------------------------------------------------------------
def build_diamond_blastp_cmd(
    *,
    db_prefix: Path,
    query_fasta: Path,
    output_tsv: Path,
    pass_cfg: dict,
) -> list[str]:
    """
    Build a 'diamond blastp' command entirely from config.py.

    pass_cfg schema:
      - sensitivity_flag: str (e.g. "--mid-sensitive")
      - max_target_seqs: int
      - evalue: float|int
      - masking: int (0/1)
      - matrix: str
      - max_hsps: int
    """
    cmd: list[str] = [
        DIAMOND_BIN,
        "blastp",
        "-d",
        str(db_prefix),
        "-q",
        str(query_fasta),
        "-o",
        str(output_tsv),
        "--threads",
        str(DIAMOND_THREADS),
        "--outfmt",
        "6",
        *DIAMOND_OUTFMT_FIELDS,
    ]

    sens = pass_cfg.get("sensitivity_flag")
    if sens:
        cmd.append(str(sens))

    cmd.extend(["--max-target-seqs", str(pass_cfg["max_target_seqs"])])
    cmd.extend(["-e", str(pass_cfg["evalue"])])
    cmd.extend(["--masking", str(pass_cfg["masking"])])
    cmd.extend(["--matrix", str(pass_cfg["matrix"])])
    cmd.extend(["--max-hsps", str(pass_cfg["max_hsps"])])

    return cmd


def run_command_capture_stdout_stderr(
    logger,
    cmd: list[str],
    log_file: Path,
    step_name: str,
) -> None:
    logger.info(f"Running {step_name} ...")
    start = time.time()

    log_file.parent.mkdir(parents=True, exist_ok=True)
    try:
        with log_file.open("a", encoding="utf-8") as lf:
            subprocess.run(cmd, stdout=lf, stderr=subprocess.STDOUT, check=True)
    except subprocess.CalledProcessError as e:
        logger.error(f"ERROR in {step_name}: exit {e.returncode}. See {log_file}")
        raise RuntimeError(f"{step_name} failed (exit={e.returncode})") from e

    logger.info(f"{step_name} completed in {time.time() - start:.2f} sec.")


def run_diamond_blastp_atomic(
    logger,
    *,
    db_prefix: Path,
    query_fasta: Path,
    output_tsv: Path,
    pass_cfg: dict,
    log_file: Path,
    step_name: str,
) -> None:
    # Run DIAMOND blastp writing to a temp file first, then atomically renaming to final.
    # Prevents half-written TSVs from being treated as complete after crashes.
    
    output_tsv.parent.mkdir(parents=True, exist_ok=True)
    tmp = output_tsv.with_suffix(output_tsv.suffix + ".tmp")

    # Ensure no stale tmp from a prior crash
    if tmp.exists():
        tmp.unlink()

    cmd = build_diamond_blastp_cmd(
        db_prefix=db_prefix,
        query_fasta=query_fasta,
        output_tsv=tmp,  # write temp
        pass_cfg=pass_cfg,
    )
    run_command_capture_stdout_stderr(logger, cmd, log_file, step_name)

    # Validate tmp looks OK
    if not looks_like_tsv(tmp, min_lines=1):
        raise RuntimeError(f"{step_name} produced invalid TSV: {tmp}")

    atomic_copy_replace(tmp, output_tsv)
    logger.info(f"{step_name} output committed atomically: {output_tsv}")


# -----------------------------------------------------------------------------
# FASTA building (atomic)
# -----------------------------------------------------------------------------
def write_selected_fasta_atomic(
    *,
    src_fasta: Path,
    dst_fasta: Path,
    allowed_ids: set[str],
) -> int:
  
    # Write selected sequences atomically: write to tmp then os.replace(). Returns number of sequences written.
    dst_fasta.parent.mkdir(parents=True, exist_ok=True)
    tmp = dst_fasta.with_suffix(dst_fasta.suffix + ".tmp")
    if tmp.exists():
        tmp.unlink()

    written = 0
    with src_fasta.open("r", encoding="utf-8") as fin, tmp.open("w", encoding="utf-8") as fout:
        write_seq = False
        seq_lines: list[str] = []

        for line in fin:
            if line.startswith(">"):
                if write_seq and seq_lines:
                    fout.writelines(seq_lines)
                    written += 1
                seq_id = line[1:].strip().split()[0]
                write_seq = seq_id in allowed_ids
                seq_lines = [line] if write_seq else []
            else:
                if write_seq:
                    seq_lines.append(line)

        if write_seq and seq_lines:
            fout.writelines(seq_lines)
            written += 1

    if written == 0:
        tmp.unlink(missing_ok=True)
        raise RuntimeError(f"Filtered FASTA would be empty: {dst_fasta}")

    if not looks_like_fasta(tmp):
        tmp.unlink(missing_ok=True)
        raise RuntimeError(f"Filtered FASTA failed validation: {tmp}")

    atomic_copy_replace(tmp, dst_fasta)
    return written


def append_nohits_fasta_atomic(
    *,
    src_fasta: Path,
    dst_fasta: Path,
    nohit_ids: set[str],
) -> int:
    """
    Append no-hit sequences atomically by creating a new file:
      new = old + appended
    then replace old.
    Returns number appended.
    """
    if not dst_fasta.exists():
        raise RuntimeError(f"Pass2 FASTA missing before append: {dst_fasta}")
    if not looks_like_fasta(dst_fasta):
        raise RuntimeError(f"Pass2 FASTA invalid before append: {dst_fasta}")

    tmp = dst_fasta.with_suffix(dst_fasta.suffix + ".tmp")
    if tmp.exists():
        tmp.unlink()

    appended = 0
    # Copy old content, then append new records
    with dst_fasta.open("r", encoding="utf-8") as fin_old, tmp.open("w", encoding="utf-8") as fout:
        for line in fin_old:
            fout.write(line)

        write_seq = False
        seq_lines: list[str] = []
        with src_fasta.open("r", encoding="utf-8") as fin_src:
            for line in fin_src:
                if line.startswith(">"):
                    if write_seq and seq_lines:
                        fout.writelines(seq_lines)
                        appended += 1
                    seq_id = line[1:].strip().split()[0]
                    write_seq = seq_id in nohit_ids
                    seq_lines = [line] if write_seq else []
                else:
                    if write_seq:
                        seq_lines.append(line)

            if write_seq and seq_lines:
                fout.writelines(seq_lines)
                appended += 1

    if not looks_like_fasta(tmp):
        tmp.unlink(missing_ok=True)
        raise RuntimeError(f"Pass2 FASTA failed validation after append: {tmp}")

    atomic_copy_replace(tmp, dst_fasta)
    return appended


# -----------------------------------------------------------------------------
# Steps
# -----------------------------------------------------------------------------
def step_filter_pass2_seed(logger, args) -> None:
    logger.info(f"Loading PDB features table: {PDB_FEATURES_TSV}")
    df = pd.read_csv(PDB_FEATURES_TSV, sep="\t", dtype=str)

    df["qseq_length"] = pd.to_numeric(df["qseq_length"], errors="coerce")
    df["qseq_n_source_organisms"] = pd.to_numeric(df["qseq_n_source_organisms"], errors="coerce")

    tax_desc = (
        df["qseq_tax_name"].fillna("").str.lower()
        + " "
        + df["qseq_description"].fillna("").str.lower()
    )

    fail_length = df["qseq_length"] <= MIN_QUERY_LENGTH
    fail_n_source = df["qseq_n_source_organisms"].isna() | (df["qseq_n_source_organisms"] > MAX_ORGANISMS)
    fail_hiv = tax_desc.str.contains("human immunodeficiency|hiv", na=False)
    fail_zinc = tax_desc.str.contains("zinc|thioredoxin", na=False)

    fail_any = fail_length | fail_n_source | fail_hiv | fail_zinc
    fail_ids = set(df.loc[fail_any, "qseqid"].dropna().tolist())

    logger.info(f"Initial filter step: {len(fail_ids)} sequences flagged for Pass 2.")
    n_written = write_selected_fasta_atomic(
        src_fasta=PDB_UNIQUE_FASTA,
        dst_fasta=PDB_FILTERED_SEQS_PASS2_FASTA,
        allowed_ids=fail_ids,
    )
    logger.info(f"Wrote {n_written} sequences: {PDB_FILTERED_SEQS_PASS2_FASTA}")

    mark_done(FILTER_PASS2_DONE, f"OK: wrote {n_written} filtered sequences")
    logger.info(f"Marked done: {FILTER_PASS2_DONE}")


def step_makedb(logger, args) -> None:
    dmnd = dmnd_path(UNIREF90_DIAMOND_DB_PREFIX)

    # Build DB atomically is DIAMOND’s job; we just validate afterward.
    run_command_capture_stdout_stderr(
        logger=logger,
        cmd=[
            DIAMOND_BIN,
            "makedb",
            "--in",
            str(UNIREF90_PROCESSED_FASTA),
            "-d",
            str(UNIREF90_DIAMOND_DB_PREFIX),
            "--threads",
            str(DIAMOND_THREADS),
        ],
        log_file=MAKEDB_LOG,
        step_name="Diamond makedb",
    )

    if not is_nonempty_file(dmnd):
        raise RuntimeError(f"DIAMOND DB missing/empty after makedb: {dmnd}")

    mark_done(DIAMOND_MAKEDB_DONE, f"OK: db present {dmnd.name}")
    logger.info(f"Marked done: {DIAMOND_MAKEDB_DONE}")


def step_pass1(logger, args) -> None:
    run_diamond_blastp_atomic(
        logger,
        db_prefix=UNIREF90_DIAMOND_DB_PREFIX,
        query_fasta=PDB_UNIQUE_FASTA,
        output_tsv=DIAMOND_HITS_PASS1_TSV,
        pass_cfg=DIAMOND_PASS1,
        log_file=SEARCH1_LOG,
        step_name="Diamond blastp Pass 1",
    )
    mark_done(DIAMOND_PASS1_DONE, f"OK: wrote {DIAMOND_HITS_PASS1_TSV.name}")
    logger.info(f"Marked done: {DIAMOND_PASS1_DONE}")


def step_append_nohits(logger, args) -> None:
    # Uses pass1 hits to compute which qseqid never appeared as qseqid in results
    if not looks_like_tsv(DIAMOND_HITS_PASS1_TSV, min_lines=1):
        raise RuntimeError(f"Pass 1 TSV invalid/missing: {DIAMOND_HITS_PASS1_TSV}")

    # The filtered FASTA must exist before appending
    if not looks_like_fasta(PDB_FILTERED_SEQS_PASS2_FASTA):
        raise RuntimeError(f"Pass 2 FASTA seed missing/invalid: {PDB_FILTERED_SEQS_PASS2_FASTA}")

    # Optional: support legacy PDB_NOHIT_MARKER for backward compatibility
    if PDB_NOHIT_MARKER.exists() and not APPEND_NOHITS_DONE.exists():
        logger.info(f"Legacy no-hit marker exists ({PDB_NOHIT_MARKER}); creating .done marker too.")

    logger.info("Computing hit set from Pass 1 results ...")
    hit_set: set[str] = set()
    with DIAMOND_HITS_PASS1_TSV.open("r", encoding="utf-8") as f:
        for line in f:
            if line.strip():
                hit_set.add(line.split("\t", 1)[0])

    # Compute no-hit IDs by scanning FASTA headers
    nohit_ids: set[str] = set()
    with PDB_UNIQUE_FASTA.open("r", encoding="utf-8") as fin:
        for line in fin:
            if line.startswith(">"):
                seq_id = line[1:].strip().split()[0]
                if seq_id not in hit_set:
                    nohit_ids.add(seq_id)

    logger.info(f"No-hit queries: {len(nohit_ids)}")

    appended = append_nohits_fasta_atomic(
        src_fasta=PDB_UNIQUE_FASTA,
        dst_fasta=PDB_FILTERED_SEQS_PASS2_FASTA,
        nohit_ids=nohit_ids,
    )

    # Write both markers: project legacy + new done
    atomic_write_text(PDB_NOHIT_MARKER, f"Appended {appended} no-hit queries\n")
    mark_done(APPEND_NOHITS_DONE, f"OK: appended {appended} no-hit sequences")
    logger.info(f"Appended {appended} sequences; marked done: {APPEND_NOHITS_DONE}")


def step_pass2(logger, args) -> None:
    if not looks_like_fasta(PDB_FILTERED_SEQS_PASS2_FASTA):
        raise RuntimeError(f"Pass 2 FASTA missing/invalid: {PDB_FILTERED_SEQS_PASS2_FASTA}")

    run_diamond_blastp_atomic(
        logger,
        db_prefix=UNIREF90_DIAMOND_DB_PREFIX,
        query_fasta=PDB_FILTERED_SEQS_PASS2_FASTA,
        output_tsv=DIAMOND_HITS_PASS2_TSV,
        pass_cfg=DIAMOND_PASS2,
        log_file=SEARCH2_LOG,
        step_name="Diamond blastp Pass 2",
    )
    mark_done(DIAMOND_PASS2_DONE, f"OK: wrote {DIAMOND_HITS_PASS2_TSV.name}")
    logger.info(f"Marked done: {DIAMOND_PASS2_DONE}")


# -----------------------------------------------------------------------------
# CLI + orchestration
# -----------------------------------------------------------------------------
def parse_args(argv: list[str]) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Two-pass DIAMOND search with restart-safe step markers and atomic outputs."
    )

    p.add_argument(
        "--resume",
        action="store_true",
        default=True,
        help="Resume based on .done markers (default: True).",
    )
    p.add_argument(
        "--no-resume",
        dest="resume",
        action="store_false",
        help="Do not use .done markers to skip steps.",
    )

    p.add_argument(
        "--force",
        action="store_true",
        help="Rerun ALL steps (clears all .done markers first).",
    )

    p.add_argument(
        "--force-steps",
        nargs="*",
        default=[],
        choices=ALL_STEPS,
        help=f"Rerun specific steps only. Choices: {', '.join(ALL_STEPS)}",
    )

    p.add_argument(
        "--rebuild-pass2-fasta",
        action="store_true",
        help="Rebuild pass2 fasta seed and re-append nohits (clears related markers).",
    )

    return p.parse_args(argv)


def clear_markers_for_rebuild_pass2(logger) -> None:
    """
    If user wants to rebuild the pass2 FASTA, we must clear:
      - FILTER_PASS2_DONE (seed creation)
      - APPEND_NOHITS_DONE (append step)
      - PDB_NOHIT_MARKER (legacy marker)
      - DIAMOND_PASS2_DONE (downstream depends on fasta content)
    """
    for f in (FILTER_PASS2_DONE, APPEND_NOHITS_DONE, DIAMOND_PASS2_DONE, PDB_NOHIT_MARKER):
        if f.exists():
            logger.info(f"Clearing marker: {f}")
            f.unlink()


def clear_all_done_markers(logger) -> None:
    for f in STEP_DONE_FILES.values():
        if f.exists():
            logger.info(f"Clearing done marker: {f}")
            f.unlink()
    if PDB_NOHIT_MARKER.exists():
        logger.info(f"Clearing legacy marker: {PDB_NOHIT_MARKER}")
        PDB_NOHIT_MARKER.unlink()


def main() -> None:
    args = parse_args(sys.argv[1:])
    ensure_directories()

    logger = get_logger(__name__, "run_diamond_search")
    start_run(logger, "run_diamond_search", argv=sys.argv)
    
    require_diamond(logger)

    if args.force:
        clear_all_done_markers(logger)

    if args.rebuild_pass2_fasta:
        clear_markers_for_rebuild_pass2(logger)
        # Ensure we rerun these steps even in resume mode
        if STEP_FILTER not in args.force_steps:
            args.force_steps.append(STEP_FILTER)
        if STEP_APPEND not in args.force_steps:
            args.force_steps.append(STEP_APPEND)

    start_time = time.time()

    # Step: filter pass2 seed
    if should_run_step(STEP_FILTER, args, logger):
        step_filter_pass2_seed(logger, args)
    else:
        # If skipping, still validate the expected output exists
        if not looks_like_fasta(PDB_FILTERED_SEQS_PASS2_FASTA):
            raise RuntimeError(f"Filter step skipped but output missing/invalid: {PDB_FILTERED_SEQS_PASS2_FASTA}")

    # Step: makedb
    if should_run_step(STEP_MAKEDB, args, logger):
        step_makedb(logger, args)
    else:
        dmnd = dmnd_path(UNIREF90_DIAMOND_DB_PREFIX)
        if not is_nonempty_file(dmnd):
            raise RuntimeError(f"makedb step skipped but db missing/empty: {dmnd}")

    # Step: pass1
    if should_run_step(STEP_PASS1, args, logger):
        step_pass1(logger, args)
    else:
        if not looks_like_tsv(DIAMOND_HITS_PASS1_TSV, min_lines=1):
            raise RuntimeError(f"Pass1 skipped but TSV missing/invalid: {DIAMOND_HITS_PASS1_TSV}")

    # Step: append nohits
    if should_run_step(STEP_APPEND, args, logger):
        step_append_nohits(logger, args)
    else:
        # Ensure the FASTA exists and looks valid
        if not looks_like_fasta(PDB_FILTERED_SEQS_PASS2_FASTA):
            raise RuntimeError(f"Append step skipped but pass2 FASTA invalid: {PDB_FILTERED_SEQS_PASS2_FASTA}")

    # Step: pass2
    if should_run_step(STEP_PASS2, args, logger):
        step_pass2(logger, args)
    else:
        if not looks_like_tsv(DIAMOND_HITS_PASS2_TSV, min_lines=1):
            raise RuntimeError(f"Pass2 skipped but TSV missing/invalid: {DIAMOND_HITS_PASS2_TSV}")

    logger.info(f"Workflow complete in {time.time() - start_time:.2f} seconds.")


if __name__ == "__main__":
    main()
