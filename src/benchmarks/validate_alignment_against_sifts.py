#!/usr/bin/env python3
"""
Workflow (DW vs SeqMapping alignment comparison)
- Load SIFTS_SEQMAPPING_HITS_TSV and keep only rows where:
    1) stat_label == "yes_diamond_yes_sifts_FullMatches"
    2) qseq_source_type != "synthetic"
    3) all accessions in `sifts_hits` are fully contained in `siftsHits_in_top5filtered`
       (both columns can contain 1 or multiple IDs).
- From those rows, build an allowlist of (qseqid, sseqid) pairs from `sifts_hits`.
- Stream-parse ALIGNED_FILTERED_HITS_TSV (records span multiple lines due to alignment_pretty).
- For each aligned (qseqid, sseqid) record:
    - Skip unless (qseqid, sseqid) is in the allowlist.
    - Query DW Mongo (SIFTS provenance) to get DW alignment length and starts.
    - Compare DW vs SeqMapping (length, qstart, sstart) and record mismatch status.
- Merge in qseq_length (PDB_FEATURES_TSV) and sseq_length (UNIREF90_FEATURES_TSV).
- Write the final comparison TSV (ALIGN_COMPARISON).
- Log allowlisted pairs missing from ALIGNED_FILTERED_HITS_TSV.
"""

from __future__ import annotations

import sys
import time
import logging
import re
from pathlib import Path
from typing import Iterator, Dict, Any, Optional, List, Tuple, Set

import pandas as pd
from pymongo import MongoClient

from seqmapping.utils.config import (
    DW_MONGO_URI,
    DW_MONGO_DB,
    DW_MONGO_COLLECTION,
)

from seqmapping.utils.paths import (
    ALIGNED_FILTERED_HITS_TSV,
    SIFTS_SEQMAPPING_HITS_TSV,
    PDB_FEATURES_TSV,
    UNIREF90_FEATURES_TSV,
    ALIGN_COMPARISON_TSV,
    ensure_benchmark_directories,
)

from seqmapping.utils.logging import get_benchmark_logger, start_run


# ---------------------- Tunables ----------------------
LOG_INTERVAL = 1000
UNIREF_CHUNKSIZE = 2_000_000  # for large UniRef features file

TARGET_STAT_LABEL = "yes_diamond_yes_sifts_FullMatches"


# -----------------------------------------------------------------------------
# DW query
# -----------------------------------------------------------------------------
def get_dw_alignment_lengths(
    collection,
    qseqid: str,
    sseqid: str,
) -> Optional[Tuple[int, Optional[int], Optional[int]]]:
    doc = collection.find_one(
        {
            "rcsb_id": qseqid,
            "rcsb_polymer_entity_align.reference_database_accession": sseqid,
            "rcsb_polymer_entity_align.provenance_source": "SIFTS",
        },
        {"rcsb_polymer_entity_align.$": 1, "_id": 0},
    )

    if not doc or "rcsb_polymer_entity_align" not in doc:
        return None

    align_obj = doc["rcsb_polymer_entity_align"][0]
    aligned_regions = align_obj.get("aligned_regions", [])
    if not aligned_regions:
        return None

    total_length = int(sum(int(r.get("length", 0) or 0) for r in aligned_regions))
    entity_beg_seq_id = aligned_regions[0].get("entity_beg_seq_id")
    ref_beg_seq_id = aligned_regions[0].get("ref_beg_seq_id")

    try:
        entity_beg_seq_id = int(entity_beg_seq_id) if entity_beg_seq_id is not None else None
    except Exception:
        entity_beg_seq_id = None

    try:
        ref_beg_seq_id = int(ref_beg_seq_id) if ref_beg_seq_id is not None else None
    except Exception:
        ref_beg_seq_id = None

    return total_length, entity_beg_seq_id, ref_beg_seq_id


# -----------------------------------------------------------------------------
# Parsing ALIGNED_FILTERED_HITS_TSV with embedded newlines in alignment_pretty
# -----------------------------------------------------------------------------
ALIGNED_HEADERS = [
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "bitscore",
    "qcovhsp",
    "scovhsp",
    "alignment_pretty",
]

EXPECTED_TAB_COLS = 14  # including alignment_pretty as the last column on the "start line"


def _is_record_start_line(line: str) -> bool:
    if not line:
        return False
    if line[0].isspace():
        return False
    return line.count("\t") >= (EXPECTED_TAB_COLS - 1)


def iter_aligned_records(path: str | Path, logger: logging.Logger) -> Iterator[Dict[str, Any]]:
    path = Path(path)
    with path.open("r", encoding="utf-8", errors="replace") as f:
        header = f.readline()
        if not header:
            return

        header_cols = header.rstrip("\n").split("\t")
        if header_cols[:13] != ALIGNED_HEADERS[:13]:
            logger.warning(
                "Unexpected header columns in ALIGNED_FILTERED_HITS_TSV. "
                "Proceeding with parser assuming smith_waterman_align.py format."
            )

        current_start: Optional[str] = None
        current_cont: List[str] = []

        for raw in f:
            line = raw.rstrip("\n")

            if _is_record_start_line(line):
                if current_start is not None:
                    yield _parse_record_block(current_start, current_cont, logger)
                current_start = line
                current_cont = []
            else:
                if current_start is None:
                    continue
                current_cont.append(line)

        if current_start is not None:
            yield _parse_record_block(current_start, current_cont, logger)


def _parse_record_block(start_line: str, cont_lines: List[str], logger: logging.Logger) -> Dict[str, Any]:
    parts = start_line.split("\t")
    if len(parts) < 13:
        logger.warning(f"Malformed record start line (too few cols): {start_line[:200]!r}")
        parts = (parts + [""] * 13)[:13] + [""]

    if len(parts) == 13:
        parts.append("")
    elif len(parts) > 14:
        parts = parts[:13] + ["\t".join(parts[13:])]

    base = dict(zip(ALIGNED_HEADERS[:13], parts[:13]))
    pretty_first = parts[13] if len(parts) > 13 else ""
    pretty_rest = "\n".join(cont_lines) if cont_lines else ""
    pretty = pretty_first
    if pretty_rest:
        pretty = (pretty + "\n" + pretty_rest) if pretty else pretty_rest

    base["alignment_pretty"] = pretty
    return base


# -----------------------------------------------------------------------------
# Allowed pairs: FullMatches, non-synthetic, and SIFTS hits within top5 filtered
# -----------------------------------------------------------------------------
_SPLIT_RE = re.compile(r"[,\s;|]+")


def _parse_id_list(cell: Any) -> List[str]:
    if cell is None:
        return []
    s = str(cell).strip()
    if s == "" or s.lower() == "nan":
        return []

    if (s.startswith("[") and s.endswith("]")) or (s.startswith("(") and s.endswith(")")):
        s = s[1:-1].strip()

    s = s.replace('"', "").replace("'", "")
    return [p.strip() for p in _SPLIT_RE.split(s) if p.strip()]


def _norm_text(x: Any) -> str:
    return str(x).strip().lower() if x is not None else ""


def load_allowed_pairs(
    labeled_hits_tsv: str | Path,
    logger: logging.Logger,
    stat_label: str = TARGET_STAT_LABEL,
    sifts_col: str = "sifts_hits",
    top5_col: str = "siftsHits_in_top5filtered",
    source_type_col: str = "qseq_source_type",
) -> Dict[str, Set[str]]:
    """
    Keep rows where:
      - stat_label == TARGET_STAT_LABEL (case-insensitive, stripped)
      - qseq_source_type != "synthetic" (case-insensitive, stripped)
      - set(sifts_hits) ⊆ set(siftsHits_in_top5filtered)

    Allowlist pairs are then (qseqid, each sifts_hit).
    """
    labeled_hits_tsv = Path(labeled_hits_tsv)
    logger.info(f"Loading labeled hits table from {labeled_hits_tsv}")

    usecols = ["qseqid", "stat_label", sifts_col, top5_col, source_type_col]
    df = pd.read_csv(labeled_hits_tsv, sep="\t", usecols=usecols, dtype=str)

    before = len(df)

    want = (stat_label or "").strip().lower()
    df["stat_label_norm"] = df["stat_label"].astype(str).str.strip().str.lower()
    df["source_type_norm"] = df[source_type_col].astype(str).str.strip().str.lower()

    # stat_label filter
    df = df[df["stat_label_norm"] == want].copy()

    # non-synthetic filter
    df = df[df["source_type_norm"] != "synthetic"].copy()

    # subset filter: sifts_hits within top5
    keep_mask = []
    for _, row in df.iterrows():
        sifts_ids = set(_parse_id_list(row.get(sifts_col)))
        top5_ids = set(_parse_id_list(row.get(top5_col)))
        keep_mask.append(bool(sifts_ids) and sifts_ids.issubset(top5_ids))

    df = df.loc[keep_mask].copy()

    after = len(df)
    logger.info(
        f"Filtered rows: start={before:,}, after stat_label+non-synthetic+subset={after:,}"
    )

    allow: Dict[str, Set[str]] = {}
    n_pairs = 0

    for _, row in df.iterrows():
        q = str(row.get("qseqid") or "").strip()
        if not q:
            continue
        sifts_ids = _parse_id_list(row.get(sifts_col))
        if not sifts_ids:
            continue
        sset = allow.setdefault(q, set())
        before_sz = len(sset)
        sset.update(sifts_ids)
        n_pairs += (len(sset) - before_sz)

    logger.info(f"Built allowlist: {len(allow):,} qseqids, {n_pairs:,} pairs")

    shown = 0
    for q, sset in allow.items():
        for s in sorted(sset):
            logger.info(f"ALLOWLIST example: {q}\t{s}")
            shown += 1
            if shown >= 10:
                break
        if shown >= 10:
            break

    return allow


# -----------------------------------------------------------------------------
# Feature loaders
# -----------------------------------------------------------------------------
def load_pdb_lengths(pdb_features_tsv: str | Path, logger: logging.Logger) -> pd.DataFrame:
    logger.info(f"Loading PDB features from {pdb_features_tsv}")
    df = pd.read_csv(pdb_features_tsv, sep="\t", usecols=["qseqid", "qseq_length"])
    df["qseq_length"] = pd.to_numeric(df["qseq_length"], errors="coerce")
    return df


def load_uniref_lengths(uniref_features_tsv: str | Path, logger: logging.Logger) -> pd.DataFrame:
    logger.info(f"Loading UniRef90 features (sseqid, sseq_length) from {uniref_features_tsv}")
    usecols = ["sseqid", "sseq_length"]
    chunks: List[pd.DataFrame] = []
    total_rows = 0
    t0 = time.time()

    for chunk in pd.read_csv(uniref_features_tsv, sep="\t", usecols=usecols, chunksize=UNIREF_CHUNKSIZE):
        chunks.append(chunk)
        total_rows += len(chunk)
        if total_rows and (total_rows % (10 * UNIREF_CHUNKSIZE) == 0):
            elapsed = time.time() - t0
            logger.info(f"  Loaded {total_rows:,} UniRef rows so far ({elapsed/60:.1f} min)")

    df = pd.concat(chunks, ignore_index=True) if chunks else pd.DataFrame(columns=usecols)
    df["sseq_length"] = pd.to_numeric(df["sseq_length"], errors="coerce")
    logger.info(f"Finished loading UniRef features: {len(df):,} entries total")
    return df


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------
def main() -> None:
    ensure_benchmark_directories()

    log_stem = Path(__file__).stem
    logger = get_benchmark_logger(__file__)
    start_run(logger, run_name=log_stem, argv=sys.argv)

    t0 = time.time()
    logger.info("=== Starting DW vs alignment comparison ===")
    logger.info(f"Input aligned hits: {ALIGNED_FILTERED_HITS_TSV}")
    logger.info(f"Labeled hits table (filter source): {SIFTS_SEQMAPPING_HITS_TSV}")
    logger.info(f"Target stat_label: {TARGET_STAT_LABEL!r}")
    logger.info(f"Output: {ALIGN_COMPARISON}")
    logger.info(f"Mongo: db={DW_MONGO_DB}, collection={DW_MONGO_COLLECTION}")

    # Mongo connection
    try:
        client = MongoClient(DW_MONGO_URI)
        db = client[DW_MONGO_DB]
        collection = db[DW_MONGO_COLLECTION]
    except Exception as e:
        logger.error(f"Failed to connect to MongoDB: {e}")
        logger.exception(e)
        raise SystemExit(1)

    # Build allowlist from SIFTS_SEQMAPPING_HITS_TSV under requested constraints
    try:
        allowed = load_allowed_pairs(
            SIFTS_SEQMAPPING_HITS_TSV,
            logger,
            stat_label=TARGET_STAT_LABEL,
            sifts_col="sifts_hits",
            top5_col="siftsHits_in_top5filtered",
            source_type_col="qseq_source_type",
        )
    except Exception as e:
        logger.error(f"Failed to load/build allowlist from labeled hits table: {e}")
        logger.exception(e)
        raise SystemExit(1)

    expected_pairs: Set[Tuple[str, str]] = set()
    for q, sset in allowed.items():
        for s in sset:
            expected_pairs.add((q, s))
    seen_pairs: Set[Tuple[str, str]] = set()

    # Load feature lengths
    try:
        pdb_features = load_pdb_lengths(PDB_FEATURES_TSV, logger)
    except Exception as e:
        logger.error(f"Failed to load PDB features: {e}")
        logger.exception(e)
        raise SystemExit(1)

    try:
        uniref_features = load_uniref_lengths(UNIREF90_FEATURES_TSV, logger)
    except Exception as e:
        logger.error(f"Failed to load UniRef features: {e}")
        logger.exception(e)
        raise SystemExit(1)

    # Process aligned pairs
    results: List[List[Any]] = []
    processed_total = 0
    processed_kept = 0
    skipped_not_allowed = 0

    try:
        for rec in iter_aligned_records(ALIGNED_FILTERED_HITS_TSV, logger):
            processed_total += 1

            qseqid = str(rec["qseqid"])
            sseqid = str(rec["sseqid"])

            sset = allowed.get(qseqid)
            if not sset or sseqid not in sset:
                skipped_not_allowed += 1
                continue

            seen_pairs.add((qseqid, sseqid))
            processed_kept += 1

            length_seqmap = _to_int(rec.get("length"))
            qstart = _to_int(rec.get("qstart"))
            sstart = _to_int(rec.get("sstart"))

            dw_data = get_dw_alignment_lengths(collection, qseqid, sseqid)
            if not dw_data:
                results.append([
                    qseqid, sseqid,
                    length_seqmap, None,
                    qstart, None,
                    sstart, None,
                    "missing data",
                ])
            else:
                length_dw, entity_beg_seq_id, ref_beg_seq_id = dw_data

                mismatches = []
                if length_seqmap != (length_dw if length_dw is not None else length_seqmap):
                    mismatches.append("length_seqmap")
                if qstart != (entity_beg_seq_id if entity_beg_seq_id is not None else qstart):
                    mismatches.append("qstart")
                if sstart != (ref_beg_seq_id if ref_beg_seq_id is not None else sstart):
                    mismatches.append("sstart")

                mismatch_status = "all match" if not mismatches else ",".join(mismatches)

                results.append([
                    qseqid, sseqid,
                    length_seqmap, length_dw,
                    qstart, entity_beg_seq_id,
                    sstart, ref_beg_seq_id,
                    mismatch_status,
                ])

            if LOG_INTERVAL and (processed_kept % LOG_INTERVAL == 0):
                elapsed = time.time() - t0
                logger.info(
                    f"Kept/processed {processed_kept:,} allowlisted pairs "
                    f"(total parsed {processed_total:,}, skipped {skipped_not_allowed:,}), "
                    f"elapsed {elapsed/60:.1f} min"
                )

    except Exception as e:
        logger.error(f"Failed while parsing/processing aligned hits: {e}")
        logger.exception(e)
        raise SystemExit(1)

    logger.info(f"Total parsed records from aligned TSV: {processed_total:,}")
    logger.info(f"Total kept (allowlisted) pairs compared: {processed_kept:,}")
    logger.info(f"Total skipped (not in allowlist): {skipped_not_allowed:,}")

    missing_in_aligned = expected_pairs - seen_pairs
    logger.info(f"Expected allowlisted pairs: {len(expected_pairs):,}")
    logger.info(f"Found allowlisted pairs in aligned TSV: {len(seen_pairs):,}")
    logger.info(f"Missing allowlisted pairs in aligned TSV: {len(missing_in_aligned):,}")
    for q, s in sorted(missing_in_aligned)[:25]:
        logger.warning(f"MISSING in ALIGNED_FILTERED_HITS_TSV: {q}\t{s}")

    out_df = pd.DataFrame(
        results,
        columns=[
            "qseqid", "sseqid",
            "length_seqmap", "length_dw",
            "qstart", "entity_beg_seq_id",
            "sstart", "ref_beg_seq_id",
            "any_mismatch",
        ],
    )

    out_df = out_df.merge(pdb_features, on="qseqid", how="left")
    out_df = out_df.merge(uniref_features, on="sseqid", how="left")

    out_df = out_df[
        [
            "qseqid", "sseqid",
            "qseq_length", "sseq_length",
            "length_seqmap", "length_dw",
            "qstart", "entity_beg_seq_id",
            "sstart", "ref_beg_seq_id",
            "any_mismatch",
        ]
    ]

    int_cols = [
        "qseq_length", "sseq_length",
        "length_seqmap", "length_dw",
        "qstart", "entity_beg_seq_id",
        "sstart", "ref_beg_seq_id",
    ]
    for col in int_cols:
        out_df[col] = pd.to_numeric(out_df[col], errors="coerce").fillna(0).astype(int)

    out_path = Path(ALIGN_COMPARISON)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(out_path, sep="\t", index=False)

    elapsed_total = (time.time() - t0) / 60
    logger.info(f"Comparison complete. Results written to {out_path}")
    logger.info(f"Total time: {elapsed_total:.1f} min")


def _to_int(x: Any) -> int:
    try:
        if x is None:
            return 0
        if isinstance(x, (int,)):
            return int(x)
        s = str(x).strip()
        if s == "" or s.lower() == "nan":
            return 0
        return int(float(s))
    except Exception:
        return 0


if __name__ == "__main__":
    main()
