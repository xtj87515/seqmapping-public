#!/usr/bin/env python3
from __future__ import annotations

import time
from pathlib import Path
import pandas as pd

from seqmapping.utils.paths import (
    ensure_benchmark_directories,
    ALIGN_COMPARISON_TSV,  
    ALIGN_COMPARISON_REPORT,
    BENCHMARKS_DATA_DIR,
    BENCHMARKS_LOG_DIR,
    BENCHMARKS_REPORT_DIR,
)

from seqmapping.utils.logging import get_benchmark_logger, start_run


# ============================ ADDITIONAL OUTPUT NAMES 
LIST_MISSING_DW = BENCHMARKS_DATA_DIR / "pairs_missing_dw.tsv"
LIST_EMBOSS = BENCHMARKS_DATA_DIR / "pairs_for_emboss.tsv"
LIST_START_ONLY = BENCHMARKS_DATA_DIR / "pairs_large_start_diff.tsv"
LIST_REPETITIVE = BENCHMARKS_DATA_DIR / "repetitive_subjects.tsv"
LIST_RECUR_SSEQID = BENCHMARKS_DATA_DIR / "top_problem_sseqids.tsv"
LIST_RECUR_QSEQID = BENCHMARKS_DATA_DIR / "top_problem_qseqids.tsv"


# thresholds (tune later if needed)
LENGTH_DIFF_CUTOFF = 0.05
QSTART_DIFF_CUTOFF = 0.10
SSTART_DIFF_CUTOFF = 0.10
REPEAT_MIN_COUNT = 25


# ======================================================
def safe_ratio(num, den):
    out = pd.Series([pd.NA] * len(num), index=num.index, dtype="Float64")
    mask = den > 0
    out.loc[mask] = (num.loc[mask].abs() / den.loc[mask]).astype("Float64")
    return out


def fmt_int(n: int) -> str:
    return f"{n:,}"


def main():
    ensure_benchmark_directories()
    BENCHMARKS_REPORT_DIR.mkdir(parents=True, exist_ok=True)

    logger = get_benchmark_logger(__file__)
    start_run(logger, run_name="alignment_comparison_report", argv=None)

    t0 = time.time()

    logger.info("Loading ALIGN_COMPARISON_TSV table")
    df = pd.read_csv(ALIGN_COMPARISON_TSV, sep="\t", dtype={"qseqid": str, "sseqid": str})

    # -------------------- numeric conversion
    num_cols = [
        "qseq_length", "sseq_length",
        "length_seqmap", "length_dw",
        "qstart", "entity_beg_seq_id",
        "sstart", "ref_beg_seq_id",
    ]
    for c in num_cols:
        df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0).astype(int)

    # -------------------- compute ratios
    df["length_diff_ratio"] = safe_ratio(df["length_seqmap"] - df["length_dw"], df["length_dw"])
    df["qstart_diff_ratio"] = safe_ratio(df["qstart"] - df["entity_beg_seq_id"], df["qseq_length"])
    df["sstart_diff_ratio"] = safe_ratio(df["sstart"] - df["ref_beg_seq_id"], df["sseq_length"])

    # flags
    df["has_length_mismatch"] = df["any_mismatch"].str.contains("length_seqmap", na=False)
    df["missing_dw"] = (df["any_mismatch"] == "missing data") | (df["length_dw"] == 0)

    # -------------------- dataset summary
    n_pairs = len(df)
    n_queries = df["qseqid"].nunique()
    n_all_match = int((df["any_mismatch"] == "all match").sum())
    n_any_mismatch = int(n_pairs - n_all_match)

    logger.info(f"pairs={n_pairs:,} queries={n_queries:,}")

    # -------------------- inspection buckets
    missing_dw = df[df["missing_dw"]]
    missing_dw.to_csv(LIST_MISSING_DW, sep="\t", index=False, encoding="utf-8")

    # EMBOSS candidates (major length disagreements)
    emboss_candidates = df[
        (df["has_length_mismatch"]) &
        (df["length_diff_ratio"] > LENGTH_DIFF_CUTOFF) &
        (~df["missing_dw"])
    ].copy()
    emboss_candidates.sort_values("length_diff_ratio", ascending=False)\
        .to_csv(LIST_EMBOSS, sep="\t", index=False, encoding="utf-8")

    # start-only mismatches (no length mismatch) with large start ratios
    start_only = df[
        (~df["has_length_mismatch"]) &
        (
            (df["qstart_diff_ratio"] >= QSTART_DIFF_CUTOFF) |
            (df["sstart_diff_ratio"] >= SSTART_DIFF_CUTOFF)
        )
    ].copy()
    start_only.to_csv(LIST_START_ONLY, sep="\t", index=False, encoding="utf-8")

    # repetitive subjects
    high_sstart = df[df["sstart_diff_ratio"] >= SSTART_DIFF_CUTOFF].copy()
    sseq_counts = high_sstart["sseqid"].value_counts()
    repetitive_ids = sseq_counts[sseq_counts >= REPEAT_MIN_COUNT].index
    repetitive_df = high_sstart[high_sstart["sseqid"].isin(repetitive_ids)].copy()
    repetitive_df.to_csv(LIST_REPETITIVE, sep="\t", index=False, encoding="utf-8")

    # recurrence statistics
    mismatch_df = df[df["any_mismatch"] != "all match"].copy()
    mismatch_df["sseqid"].value_counts().head(50)\
        .rename_axis("sseqid").reset_index(name="count")\
        .to_csv(LIST_RECUR_SSEQID, sep="\t", index=False, encoding="utf-8")

    mismatch_df["qseqid"].value_counts().head(50)\
        .rename_axis("qseqid").reset_index(name="count")\
        .to_csv(LIST_RECUR_QSEQID, sep="\t", index=False, encoding="utf-8")

    # -------------------- mismatch breakdown
    mismatch_counts = df["any_mismatch"].value_counts()

    # -------------------- report (TXT, not Markdown)
    logger.info("Writing TXT report")

    lines = []
    lines.append("SeqMapping vs DW Alignment Comparison Report")
    lines.append("=" * 46)
    lines.append("")
    lines.append(f"Input: {ALIGN_COMPARISON_TSV}")
    lines.append("")
    lines.append(f"Queries (unique qseqid): {fmt_int(n_queries)}")
    lines.append(f"Pairs (rows):            {fmt_int(n_pairs)}")
    lines.append(f"All match:               {fmt_int(n_all_match)}")
    lines.append(f"Any mismatch:            {fmt_int(n_any_mismatch)}")
    lines.append("")

    lines.append("Definitions")
    lines.append("-" * 11)
    lines.append("length_diff_ratio = abs(length_seqmap - length_dw) / length_dw     (length_dw > 0)")
    lines.append("qstart_diff_ratio = abs(qstart - entity_beg_seq_id) / qseq_length  (qseq_length > 0)")
    lines.append("sstart_diff_ratio = abs(sstart - ref_beg_seq_id) / sseq_length     (sseq_length > 0)")
    lines.append("")

    lines.append("any_mismatch breakdown")
    lines.append("-" * 22)
    # Print as numbered list 
    for i, (label, count) in enumerate(mismatch_counts.items(), start=1):
        lines.append(f"{i:>2}  {label:<30}  {fmt_int(int(count))}")
    lines.append(f"{'':>2}  {'Total:':<30}  {fmt_int(n_pairs)}")
    lines.append("")

    lines.append("Focus analyses")
    lines.append("-" * 14)
    lines.append(f"Length mismatch > {LENGTH_DIFF_CUTOFF:.2f}: {fmt_int(len(emboss_candidates))} pairs  (see: {LIST_EMBOSS.name})")
    lines.append(f"Missing DW mapping:       {fmt_int(len(missing_dw))} pairs  (see: {LIST_MISSING_DW.name})")
    lines.append(f"Large start discrepancies:{fmt_int(len(start_only))} pairs  (see: {LIST_START_ONLY.name})")
    lines.append(f"Repetitive subject suspects (min {REPEAT_MIN_COUNT} occurrences among high sstart): {fmt_int(len(repetitive_df))} pairs  (see: {LIST_REPETITIVE.name})")
    lines.append("")

    lines.append("Generated lists")
    lines.append("-" * 15)
    lines.append(f"- {LIST_EMBOSS}  (candidates for EMBOSS-water / 3rd party)")
    lines.append(f"- {LIST_MISSING_DW}  (DW missing mappings; likely SIFTS/DW gaps or accession issues)")
    lines.append(f"- {LIST_START_ONLY}  (start-only discrepancies without length mismatch)")
    lines.append(f"- {LIST_REPETITIVE}  (repetitive subjects likely causing alternate placements)")
    lines.append(f"- {LIST_RECUR_SSEQID}  (top recurring problematic sseqids among mismatches)")
    lines.append(f"- {LIST_RECUR_QSEQID}  (top recurring problematic qseqids among mismatches)")
    lines.append("")

    lines.append("Next steps")
    lines.append("-" * 20)
    lines.append(f"1) Run EMBOSS water on pairs in: {LIST_EMBOSS.name}")
    lines.append("   - Compare EMBOSS alignment length (and optionally starts) to SeqMapping vs DW.")
    lines.append("2) Review top recurring sseqids (often repeats / low-complexity / paralogs):")
    lines.append(f"   - See: {LIST_RECUR_SSEQID.name}")
    lines.append("3) Inspect missing DW pairs (may indicate missing SIFTS segments or accession mismaps):")
    lines.append(f"   - See: {LIST_MISSING_DW.name}")
    lines.append("4) For large start discrepancies, check whether differences are explainable by:")
    lines.append("   - terminal gap conventions")
    lines.append("   - repeat regions allowing multiple equally-good alignments")
    lines.append("   - chimeras / engineered constructs")
    lines.append("")

    ALIGN_COMPARISON_REPORT.write_text("\n".join(lines), encoding="utf-8")

    elapsed = time.time() - t0
    logger.info(f"Report generation complete ({elapsed:.2f}s)")
    logger.info(f"Wrote report: {ALIGN_COMPARISON_REPORT}")


if __name__ == "__main__":
    main()
