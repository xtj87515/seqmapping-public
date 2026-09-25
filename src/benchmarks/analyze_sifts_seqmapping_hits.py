#!/usr/bin/env python3

# Simple text benchmark report for SIFTS vs seqmapping evaluation.


from __future__ import annotations

import time
import pandas as pd

from seqmapping.utils.paths import (
    ensure_benchmark_directories,
    BENCHMARKS_LOG_DIR,
    SEQMAPPING_BENCHMARK_REPORT,
    SIFTS_SEQMAPPING_HITS_TSV,
)

from seqmapping.utils.logging import get_benchmark_logger, start_run


STAT_LABEL_ORDER = [
    "no_diamond_no_sifts",
    "no_diamond_yes_sifts",
    "yes_diamond_no_sifts",
    "yes_diamond_yes_sifts_FullMatches",
    "yes_diamond_yes_sifts_NoMatches",
    "yes_diamond_yes_sifts_PartialMatches",
]


# --------------------- logging ---------------------

def setup_logging():
    ensure_benchmark_directories()
    BENCHMARKS_LOG_DIR.mkdir(parents=True, exist_ok=True)

    logger = get_benchmark_logger(__file__)
    start_run(logger, run_name="seqmapping benchmark report", argv=None)

    logger.info("Generating seqmapping benchmark report.")
    logger.info("Input:  %s", str(SIFTS_SEQMAPPING_HITS_TSV))
    logger.info("Output: %s", str(SEQMAPPING_BENCHMARK_REPORT))
    return logger


def section(title: str) -> str:
    return f"\n{'='*80}\n{title}\n{'='*80}\n"


def df_to_text(df: pd.DataFrame) -> str:
    return df.to_string(index=False) + "\n"


# --------------------- helpers ---------------------

def read_data(logger) -> pd.DataFrame:
    t0 = time.time()
    df = pd.read_csv(SIFTS_SEQMAPPING_HITS_TSV, sep="\t", dtype=str, keep_default_na=False)
    logger.info("Loaded TSV rows=%d cols=%d (%.2fs)", len(df), len(df.columns), time.time() - t0)

    numeric_cols = [
        "n_sifts_hits",
        "n_siftsHits_in_top1filtered",
        "n_siftsHits_in_top5filtered",
        "n_sseqid_in_DiamondRaw",
    ]
    for c in numeric_cols:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0).astype(int)

    # Guarantee these exist as strings for filtering/grouping
    for c in ["stat_label", "qseq_source_type", "qseq_domain_name"]:
        if c not in df.columns:
            df[c] = ""

    return df


def counts_with_percent_ordered(df: pd.DataFrame, column: str, order: list[str] | None = None) -> pd.DataFrame:

    # Count values of `column`, compute percent, and optionally impose a fixed order.
    # If `order` is provided, missing categories are included with count=0.

    total = len(df)
    vc = df[column].value_counts(dropna=False)

    if order is None:
        out = vc.rename_axis(column).reset_index(name="count")
        out["percent"] = (out["count"] / total * 100).round(2)
        return out.sort_values("count", ascending=False)

    # fixed ordering
    rows = []
    for key in order:
        cnt = int(vc.get(key, 0))
        rows.append({column: key, "count": cnt, "percent": round((cnt / total * 100), 2)})

    return pd.DataFrame(rows)


# --------------------- analyses ---------------------

def overall_stats(df: pd.DataFrame) -> str:
    out = counts_with_percent_ordered(df, "stat_label", STAT_LABEL_ORDER)
    return df_to_text(out)


def nonsynthetic_stats(df: pd.DataFrame) -> str:
    df2 = df[df["qseq_source_type"] != "synthetic"]
    out = counts_with_percent_ordered(df2, "stat_label", STAT_LABEL_ORDER)
    return df_to_text(out)


def topk_recovery(df: pd.DataFrame) -> str:
    df = df[(df["qseq_source_type"] != "synthetic") &
            (df["stat_label"] == "yes_diamond_yes_sifts_FullMatches")]

    total = len(df)
    if total == 0:
        return "No FullMatches found.\n"

    top1 = (df["n_siftsHits_in_top1filtered"] > 0).sum()
    top5 = (df["n_siftsHits_in_top5filtered"] > 0).sum()
    rescued = ((df["n_siftsHits_in_top1filtered"] == 0) &
               (df["n_siftsHits_in_top5filtered"] > 0)).sum()

    return f"""
FullMatches (non-synthetic): {total:,}

Top1 recovered expected ID: {top1:,} ({top1/total*100:.2f}%)
Top5 recovered expected ID: {top5:,} ({top5/total*100:.2f}%)
Recovered by Top5 only (ranking issue): {rescued:,} ({rescued/total*100:.2f}%)
"""


def coverage_analysis(df: pd.DataFrame) -> str:
    df = df[df["stat_label"] == "yes_diamond_yes_sifts_FullMatches"].copy()
    if len(df) == 0:
        return "No FullMatches available for coverage analysis.\n"

    df["expected_size"] = df["n_sifts_hits"].clip(lower=1)
    df["coverage"] = df["n_siftsHits_in_top5filtered"] / df["expected_size"]

    bins = pd.cut(df["coverage"], bins=[0, 0.25, 0.5, 0.75, 1.0], include_lowest=True)
    cov = bins.value_counts().reset_index()
    cov.columns = ["coverage_range", "count"]
    cov["percent"] = (cov["count"] / len(df) * 100).round(2)

    return df_to_text(cov.sort_values("coverage_range"))


def domain_breakdown(df: pd.DataFrame) -> str:
    df = df[df["qseq_source_type"] != "synthetic"]
    out = df["qseq_domain_name"].value_counts(dropna=False).rename_axis("qseq_domain_name").reset_index(name="count")
    out["percent"] = (out["count"] / len(df) * 100).round(2)
    return df_to_text(out)


def failure_by_domain(df: pd.DataFrame) -> str:
    fails = df[df["stat_label"] == "yes_diamond_yes_sifts_NoMatches"]
    if len(fails) == 0:
        return "No NoMatch cases detected.\n"
    out = fails["qseq_domain_name"].value_counts(dropna=False).rename_axis("qseq_domain_name").reset_index(name="count")
    out["percent"] = (out["count"] / len(fails) * 100).round(2)
    return df_to_text(out)


# --------------------- report ---------------------

def build_report(df: pd.DataFrame) -> str:
    text: list[str] = []
    text.append("SIFTS SEQUENCE MAPPING BENCHMARK REPORT\n")
    text.append(f"Input dataset: {SIFTS_SEQMAPPING_HITS_TSV}\n")
    text.append(f"Total queries analyzed: {len(df):,}\n")

    text.append(section("1. Overall Match Classification (fixed order)"))
    text.append(overall_stats(df))

    text.append(section("2. Non-Synthetic Only (fixed order)"))
    text.append(nonsynthetic_stats(df))

    text.append(section("3. Top-K Filtered Hit Recovery (non-synthetic FullMatches)"))
    text.append(topk_recovery(df))

    text.append(section("4. Coverage of Expected IDs (Top5) among FullMatches"))
    text.append(coverage_analysis(df))

    text.append(section("5. Domain Distribution (non-synthetic)"))
    text.append(domain_breakdown(df))

    text.append(section("6. NoMatches by Domain"))
    text.append(failure_by_domain(df))

    text.append(section("Interpretation Guide"))
    text.append(
        "FullMatches     : sequence search agrees with SIFTS mapping\n"
        "PartialMatches  : some expected proteins found\n"
        "NoMatches       : sequence search contradicts SIFTS mapping\n"
        "Top5 rescue     : ranking issue, not search sensitivity\n"
    )

    return "".join(text)


# --------------------- main ---------------------

def main():
    logger = setup_logging()
    overall_t0 = time.time()

    df = read_data(logger)
    report = build_report(df)

    SEQMAPPING_BENCHMARK_REPORT.parent.mkdir(parents=True, exist_ok=True)
    SEQMAPPING_BENCHMARK_REPORT.write_text(report)

    logger.info("Wrote report: %s", str(SEQMAPPING_BENCHMARK_REPORT))
    logger.info("Done (%.2fs)", time.time() - overall_t0)

    print(f"Benchmark report written to: {SEQMAPPING_BENCHMARK_REPORT}")


if __name__ == "__main__":
    main()
