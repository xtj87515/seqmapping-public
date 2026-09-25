#!/usr/bin/env python3

"""
Workflow
--------
1. Read `pairs_for_emboss.tsv` (pairs previously flagged as suspicious in the
   SeqMapping vs DW comparison).

2. For each (qseqid, sseqid) pair:
   - Retrieve query sequence from PDB_UNIQUE_FASTA
   - Retrieve subject sequence from UNIREF90_PROCESSED_FASTA
   - Submit pairwise alignment to EBI EMBOSS-water REST API
   - Wait for job completion and download the alignment report
   - Extract EMBOSS alignment length

3. Compare alignment lengths:
   - SeqMapping (length_seqmap)
   - DW/SIFTS (length_dw)
   - EMBOSS (length_emboss)

   Compute:
       length_ratio_seqmap_dw     = |length_seqmap - length_dw| / length_dw
       length_ratio_seqmap_emboss = |length_seqmap - length_emboss| / length_emboss

4. Outputs (written incrementally, resume-safe):
   - sifts_seqmapping_emboss_comparison.tsv  : numeric comparison table
   - sifts_seqmapping_emboss_comparison.md   : full EMBOSS alignments
   - sifts_seqmapping_emboss_comparison.txt  : report

Note:
If the script stops (API outage, network, killed job), rerunning it will NOT
restart from scratch. Previously completed pairs are detected via `length_emboss`
and skipped; only unfinished pairs are resubmitted.

"""

from __future__ import annotations

import os
import sys
import time
import logging
from pathlib import Path
from typing import Dict, Optional, List

import pandas as pd
import requests

from seqmapping.utils.paths import (
    ensure_benchmark_directories,
    PDB_UNIQUE_FASTA,
    UNIREF90_PROCESSED_FASTA,
    SIFTS_SEQMAPPING_EMBOSS_TSV,
    SIFTS_SEQMAPPING_EMBOSS_MD,
    SIFTS_SEQMAPPING_EMBOSS_REPORT,
    BENCHMARKS_DATA_DIR,
    BENCHMARKS_LOG_DIR,
    BENCHMARKS_REPORT_DIR,
)

from seqmapping.utils.logging import get_benchmark_logger, start_run


# ====================== CONFIGURATION ======================
BASE_URL = "https://www.ebi.ac.uk/Tools/services/rest/emboss_water"
USER_EMAIL = "tongji.xing@rcsb.org"
JOB_TITLE = "SIFTS_SeqMapping_EMBOSS_Water"

LIST_EMBOSS = BENCHMARKS_DATA_DIR / "pairs_for_emboss.tsv"

# Parameters for EMBOSS water
WATER_PARAMS = {
    "gapopen": "10.0",
    "gapext": "1.0",       # correct param name for EBI REST
    "matrix": "EBLOSUM62",
    "outformat": "out",
}

# ====================== DEBUG / CONTROL FLAGS ======================
DEBUG_MODE = False
DEBUG_LIMIT = 15
LOG_EVERY = 30

POLL_INTERVAL_SEC = 5
POLL_TIMEOUT_SEC = 300

# ---- throttling (optional; tune as needed)
SUBMIT_SLEEP_SEC = 0.25        # sleep after each submission attempt
POLL_EXTRA_SLEEP_SEC = 0.10    # extra sleep after each status check
RESULT_FETCH_SLEEP_SEC = 0.10  # sleep before fetching large result


# ====================== UTIL ======================
def suppress_print() -> None:
    sys.stdout = open(os.devnull, "w")


def restore_print() -> None:
    sys.stdout = sys.__stdout__


def throttle_sleep(seconds: float) -> None:
    if seconds and seconds > 0:
        time.sleep(seconds)


def safe_ratio_abs(num: float, den: float) -> Optional[float]:
    """
    abs(num) / den, returning None if den <= 0 or inputs are missing.
    """
    try:
        if den is None:
            return None
        den_f = float(den)
        if den_f <= 0:
            return None
        return abs(float(num)) / den_f
    except Exception:
        return None


def load_fasta_dict(fasta_file: str | Path) -> Dict[str, str]:
    """
    Load FASTA file into a dict: {sequence_id: sequence_string}
    Supports plain or gzipped FASTA (*.gz).
    """
    fasta_file = str(fasta_file)
    start_time = time.time()
    logging.getLogger(__name__).info(f"Loading FASTA file: {fasta_file}")

    sequences: Dict[str, str] = {}
    current_id: Optional[str] = None
    current_seq: List[str] = []
    seq_count = 0

    if fasta_file.endswith(".gz"):
        import gzip
        open_func = gzip.open
        mode = "rt"
    else:
        open_func = open
        mode = "r"

    try:
        with open_func(fasta_file, mode) as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line:
                    continue

                if line.startswith(">"):
                    # flush previous
                    if current_id is not None:
                        seq_str = "".join(current_seq)
                        sequences[current_id] = seq_str
                        seq_count += 1
                        current_seq = []

                    header_parts = line[1:].split()
                    if header_parts:
                        current_id = header_parts[0]
                    else:
                        current_id = f"seq_{seq_count + 1}"
                        logging.getLogger(__name__).warning(
                            f"Empty header at line {line_num}, assigning ID: {current_id}"
                        )

                    if seq_count > 0 and seq_count % 5_000_000 == 0:
                        elapsed = time.time() - start_time
                        logging.getLogger(__name__).info(
                            f"Loaded {seq_count:,} sequences in {elapsed:.1f}s..."
                        )
                else:
                    current_seq.append(line.upper())

            # last record
            if current_id is not None and current_seq:
                seq_str = "".join(current_seq)
                sequences[current_id] = seq_str
                seq_count += 1

        elapsed = time.time() - start_time
        logging.getLogger(__name__).info(
            f"Loaded {seq_count:,} sequences in {elapsed:.2f} seconds"
        )
        return sequences

    except Exception as e:
        logging.getLogger(__name__).error(f"Failed to load FASTA file {fasta_file}: {e}")
        raise


def submit_job(pair_id: str, seq1: str, seq2: str) -> Optional[str]:
    url = f"{BASE_URL}/run"
    payload = {
        "email": USER_EMAIL,
        "title": JOB_TITLE,
        "asequence": f">{pair_id}_A\n{seq1}",
        "bsequence": f">{pair_id}_B\n{seq2}",
        **WATER_PARAMS,
    }
    try:
        r = requests.post(url, data=payload, timeout=60)
        r.raise_for_status()
        return r.text.strip()
    except Exception:
        return None


def check_status(job_id: str) -> str:
    url = f"{BASE_URL}/status/{job_id}"
    r = requests.get(url, timeout=30)
    r.raise_for_status()
    return r.text.strip()


def get_result(job_id: str, result_type: str = "out") -> str:
    url = f"{BASE_URL}/result/{job_id}/{result_type}"
    r = requests.get(url, timeout=60)
    r.raise_for_status()
    return r.text


def wait_for_completion(job_id: str, interval: int = POLL_INTERVAL_SEC, timeout: int = POLL_TIMEOUT_SEC) -> None:
    start = time.time()
    while True:
        status = check_status(job_id)
        throttle_sleep(POLL_EXTRA_SLEEP_SEC)

        if status == "FINISHED":
            return
        if status in {"ERROR", "FAILURE"}:
            raise RuntimeError(f"Job {job_id} failed with status {status}")
        if time.time() - start > timeout:
            raise TimeoutError(f"Job {job_id} timed out after {timeout}s")

        time.sleep(interval)


def parse_alignment_length(report_text: str) -> Optional[int]:
    for line in report_text.splitlines():
        if line.strip().startswith("# Length:"):
            try:
                return int(line.strip().split()[2])
            except Exception:
                return None
    return None


def replace_headers(report_text: str, qid: str, sid: str) -> str:
    report_text = report_text.replace(f"{qid}_{sid}_A", qid)
    report_text = report_text.replace(f"{qid}_{sid}_B", sid)
    return report_text


# ====================== MAIN ======================
def main() -> None:
    ensure_benchmark_directories()
    BENCHMARKS_DATA_DIR.mkdir(parents=True, exist_ok=True)
    BENCHMARKS_LOG_DIR.mkdir(parents=True, exist_ok=True)
    BENCHMARKS_REPORT_DIR.mkdir(parents=True, exist_ok=True)

    logger = get_benchmark_logger(__file__)
    start_run(logger, run_name="compare_seqmapping_with_emboss_water", argv=None)

    # Keep console quiet
    restore_print()
    print("[INFO] Starting EMBOSS Water batch alignment...")
    suppress_print()

    t0 = time.time()

    logger.info("Inputs:")
    logger.info("  LIST_EMBOSS                    = %s", str(LIST_EMBOSS))
    logger.info("  PDB_UNIQUE_FASTA               = %s", str(PDB_UNIQUE_FASTA))
    logger.info("  UNIREF90_PROCESSED_FASTA       = %s", str(UNIREF90_PROCESSED_FASTA))
    logger.info("Outputs:")
    logger.info("  SIFTS_SEQMAPPING_EMBOSS_TSV     = %s", str(SIFTS_SEQMAPPING_EMBOSS_TSV))
    logger.info("  SIFTS_SEQMAPPING_EMBOSS_MD      = %s", str(SIFTS_SEQMAPPING_EMBOSS_MD))
    logger.info("  SIFTS_SEQMAPPING_EMBOSS_REPORT  = %s", str(SIFTS_SEQMAPPING_EMBOSS_REPORT))

    logger.info("Loading FASTA files into memory (dict lookup)...")
    # NOTE: may require large RAM for UniRef90.
    qseq_db = load_fasta_dict(PDB_UNIQUE_FASTA)
    sseq_db = load_fasta_dict(UNIREF90_PROCESSED_FASTA)

    logger.info("Reading input list for EMBOSS water...")
    in_df = pd.read_csv(LIST_EMBOSS, sep="\t", dtype=str)
    if "qseqid" not in in_df.columns or "sseqid" not in in_df.columns:
        raise RuntimeError(f"Input file missing required columns qseqid/sseqid: {LIST_EMBOSS}")

    # Rename the input ratio column for output consistency
    if "length_diff_ratio" in in_df.columns and "length_ratio_seqmap_dw" not in in_df.columns:
        in_df = in_df.rename(columns={"length_diff_ratio": "length_ratio_seqmap_dw"})

    if DEBUG_MODE:
        in_df = in_df.head(DEBUG_LIMIT).copy()
        logger.info("[DEBUG] Limiting to first %d pairs", DEBUG_LIMIT)

    in_df["pair_key"] = in_df["qseqid"].astype(str) + "\t" + in_df["sseqid"].astype(str)

    # Load or initialize output TSV (resume-safe by key)
    if SIFTS_SEQMAPPING_EMBOSS_TSV.exists():
        prev_df = pd.read_csv(SIFTS_SEQMAPPING_EMBOSS_TSV, sep="\t", dtype=str)
        # Rebuild key for older versions
        if "pair_key" not in prev_df.columns:
            if "qseqid" in prev_df.columns and "sseqid" in prev_df.columns:
                prev_df["pair_key"] = prev_df["qseqid"].astype(str) + "\t" + prev_df["sseqid"].astype(str)
            else:
                raise RuntimeError("Existing output TSV lacks qseqid/sseqid; cannot resume safely.")

        # We keep internal resume columns in-memory only
        for col in ["length_emboss", "job_id_emboss", "emboss_status", "emboss_error"]:
            if col not in prev_df.columns:
                prev_df[col] = pd.NA

        out_df = in_df.merge(
            prev_df[["pair_key", "length_emboss", "job_id_emboss", "emboss_status", "emboss_error"]],
            on="pair_key",
            how="left",
        )
        logger.info("Resuming from existing TSV (%d rows).", len(out_df))
    else:
        out_df = in_df.copy()
        out_df["length_emboss"] = pd.NA
        out_df["job_id_emboss"] = pd.NA
        out_df["emboss_status"] = pd.NA
        out_df["emboss_error"] = pd.NA
        SIFTS_SEQMAPPING_EMBOSS_TSV.parent.mkdir(parents=True, exist_ok=True)
        # write initial TSV (final columns will be cleaned before every write)
        _write_clean_tsv(out_df, SIFTS_SEQMAPPING_EMBOSS_TSV)
        logger.info("Initialized output TSV (%d rows).", len(out_df))

    total = len(out_df)
    completed = int(pd.to_numeric(out_df["length_emboss"], errors="coerce").notna().sum())
    logger.info("Already completed: %d/%d pairs", completed, total)

    if not SIFTS_SEQMAPPING_EMBOSS_MD.exists():
        SIFTS_SEQMAPPING_EMBOSS_MD.parent.mkdir(parents=True, exist_ok=True)
        SIFTS_SEQMAPPING_EMBOSS_MD.write_text("# EMBOSS Water Alignment Results\n\n", encoding="utf-8")

    # Process each pair
    for idx, row in out_df.iterrows():
        qid = str(row["qseqid"])
        sid = str(row["sseqid"])
        pair_id = f"{qid}_{sid}"

        # skip already processed
        if str(row.get("length_emboss", "")).strip() not in {"", "NA", "NaN", "nan", "None"}:
            continue

        qseq = qseq_db.get(qid)
        sseq = sseq_db.get(sid)
        if not qseq or not sseq:
            out_df.at[idx, "emboss_status"] = "missing_sequence"
            out_df.at[idx, "emboss_error"] = f"Missing sequence(s): qseq={'Y' if qseq else 'N'} sseq={'Y' if sseq else 'N'}"
            out_df.at[idx, "length_emboss"] = pd.NA
            _write_clean_tsv(out_df, SIFTS_SEQMAPPING_EMBOSS_TSV)

            with open(SIFTS_SEQMAPPING_EMBOSS_MD, "a", encoding="utf-8") as f:
                f.write(f"## {pair_id}\nMissing sequence(s) in FASTA DB(s).\n\n")
            continue

        # ---- throttling before submission
        throttle_sleep(SUBMIT_SLEEP_SEC)

        job_id = submit_job(pair_id, qseq, sseq)
        if not job_id:
            out_df.at[idx, "emboss_status"] = "submit_failed"
            out_df.at[idx, "emboss_error"] = "Job submission failed"
            out_df.at[idx, "length_emboss"] = pd.NA
            _write_clean_tsv(out_df, SIFTS_SEQMAPPING_EMBOSS_TSV)

            with open(SIFTS_SEQMAPPING_EMBOSS_MD, "a", encoding="utf-8") as f:
                f.write(f"## {pair_id}\nJob submission failed.\n\n")
            continue

        out_df.at[idx, "job_id_emboss"] = job_id
        out_df.at[idx, "emboss_status"] = "submitted"
        out_df.at[idx, "emboss_error"] = pd.NA
        _write_clean_tsv(out_df, SIFTS_SEQMAPPING_EMBOSS_TSV)

        try:
            wait_for_completion(job_id)
            throttle_sleep(RESULT_FETCH_SLEEP_SEC)

            report_text = get_result(job_id, "out")
            report_text = replace_headers(report_text, qid, sid)

            align_length = parse_alignment_length(report_text)
            out_df.at[idx, "length_emboss"] = align_length
            out_df.at[idx, "emboss_status"] = "finished"
            out_df.at[idx, "emboss_error"] = pd.NA

            # compute derived ratio: abs(length_seqmap - length_emboss) / length_emboss
            length_seqmap = row.get("length_seqmap")
            out_df.at[idx, "length_ratio_seqmap_emboss"] = safe_ratio_abs(
                (float(length_seqmap) - float(align_length)) if (length_seqmap not in [None, "", "NA", "nan", "NaN"] and align_length is not None) else 0.0,
                float(align_length) if align_length is not None else 0.0,
            )

            with open(SIFTS_SEQMAPPING_EMBOSS_MD, "a", encoding="utf-8") as f:
                f.write(f"## Alignment: {qid} vs {sid}\n\n```text\n{report_text}\n```\n\n")

            _write_clean_tsv(out_df, SIFTS_SEQMAPPING_EMBOSS_TSV)
            completed += 1

            if completed % LOG_EVERY == 0:
                logger.info("Processed %d/%d pairs. Last=%s length_emboss=%s", completed, total, pair_id, str(align_length))

        except Exception as e:
            out_df.at[idx, "emboss_status"] = "failed"
            out_df.at[idx, "emboss_error"] = str(e)
            out_df.at[idx, "length_emboss"] = pd.NA
            _write_clean_tsv(out_df, SIFTS_SEQMAPPING_EMBOSS_TSV)

            with open(SIFTS_SEQMAPPING_EMBOSS_MD, "a", encoding="utf-8") as f:
                f.write(f"## {pair_id}\nFailed: {e}\n\n")
            continue

    # Final summary report (TXT)
    finished = int((out_df["emboss_status"] == "finished").sum())
    failed = int((out_df["emboss_status"] == "failed").sum())
    missing_seq = int((out_df["emboss_status"] == "missing_sequence").sum())
    submit_failed = int((out_df["emboss_status"] == "submit_failed").sum())
    remaining = int(total - pd.to_numeric(out_df["length_emboss"], errors="coerce").notna().sum())

    txt_lines = []
    txt_lines.append("SIFTS SeqMapping vs EMBOSS-water Comparison (Batch Run)")
    txt_lines.append("=" * 56)
    txt_lines.append("")
    txt_lines.append(f"Input pairs: {LIST_EMBOSS}")
    txt_lines.append(f"Q FASTA:     {PDB_UNIQUE_FASTA}")
    txt_lines.append(f"S FASTA:     {UNIREF90_PROCESSED_FASTA}")
    txt_lines.append("")
    txt_lines.append(f"Total pairs:        {total:,}")
    txt_lines.append(f"Finished:           {finished:,}")
    txt_lines.append(f"Failed:             {failed:,}")
    txt_lines.append(f"Missing sequence:   {missing_seq:,}")
    txt_lines.append(f"Submit failed:      {submit_failed:,}")
    txt_lines.append(f"Remaining (NA len): {remaining:,}")
    txt_lines.append("")
    txt_lines.append("Outputs")
    txt_lines.append("-" * 7)
    txt_lines.append(f"TSV: {SIFTS_SEQMAPPING_EMBOSS_TSV}")
    txt_lines.append(f"MD:  {SIFTS_SEQMAPPING_EMBOSS_MD}")
    txt_lines.append(f"TXT: {SIFTS_SEQMAPPING_EMBOSS_REPORT}")
    txt_lines.append("")

    SIFTS_SEQMAPPING_EMBOSS_REPORT.parent.mkdir(parents=True, exist_ok=True)
    SIFTS_SEQMAPPING_EMBOSS_REPORT.write_text("\n".join(txt_lines), encoding="utf-8")

    elapsed = time.time() - t0
    logger.info("Completed EMBOSS batch (%0.2fs). Finished=%d/%d", elapsed, finished, total)

    restore_print()
    print(f"[INFO] Finished EMBOSS batch. Completed={finished:,}/{total:,}")
    print(f"[INFO] TSV: {SIFTS_SEQMAPPING_EMBOSS_TSV}")
    print(f"[INFO] MD:  {SIFTS_SEQMAPPING_EMBOSS_MD}")
    print(f"[INFO] TXT: {SIFTS_SEQMAPPING_EMBOSS_REPORT}")
    suppress_print()


def _write_clean_tsv(out_df: pd.DataFrame, out_path: Path) -> None:
    """
    Write TSV without internal/resume columns:
      - pair_key, job_id_emboss, emboss_status, emboss_error
    Keep:
      - length_emboss (and place it after length_seqmap and length_dw)

    Also ensures:
      - length_diff_ratio is renamed to length_ratio_seqmap_dw
      - length_ratio_seqmap_emboss exists (may be NA until EMBOSS completes)
    """
    df = out_df.copy()

    # ensure rename (in case older TSV is loaded)
    if "length_diff_ratio" in df.columns and "length_ratio_seqmap_dw" not in df.columns:
        df = df.rename(columns={"length_diff_ratio": "length_ratio_seqmap_dw"})

    # ensure derived column exists
    if "length_ratio_seqmap_emboss" not in df.columns:
        df["length_ratio_seqmap_emboss"] = pd.NA

    # compute length_ratio_seqmap_emboss where possible
    if "length_emboss" in df.columns and "length_seqmap" in df.columns:
        seqmap = pd.to_numeric(df["length_seqmap"], errors="coerce")
        emboss = pd.to_numeric(df["length_emboss"], errors="coerce")
        mask = emboss.notna() & (emboss > 0) & seqmap.notna()
        df.loc[mask, "length_ratio_seqmap_emboss"] = (
            (seqmap.loc[mask] - emboss.loc[mask]).abs() / emboss.loc[mask]
        ).astype(float)

    # drop internal columns (keep length_emboss!)
    drop_cols = [c for c in ["pair_key", "job_id_emboss", "emboss_status", "emboss_error"] if c in df.columns]
    df = df.drop(columns=drop_cols, errors="ignore")

    # reorder: put length_emboss right after length_seqmap and length_dw if present
    if all(c in df.columns for c in ["length_seqmap", "length_dw", "length_emboss"]):
        cols = list(df.columns)
        cols.remove("length_emboss")
        insert_at = cols.index("length_dw") + 1
        cols.insert(insert_at, "length_emboss")
        df = df[cols]

    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, sep="\t", index=False, encoding="utf-8")

if __name__ == "__main__":
    main()
