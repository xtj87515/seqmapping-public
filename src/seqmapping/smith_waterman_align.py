#!/usr/bin/env python3

# We previously used pyfaidx 0.8.1.4 with Python 3.8, which performed well. However, the newer pyfaidx 0.9.0.3
# that comes with Python 3.12 introduced extra overhead and significantly slowed down the alignment.
# A custom FASTA parser is instead implemented that executes much faster.

from __future__ import annotations

import sys
import time
import logging
from pathlib import Path
from typing import List, Tuple, Optional, Dict

import pandas as pd
import parasail

from seqmapping.utils.config import (
    ALIGNMENT_LINE_WIDTH,
    ALIGNMENT_LOG_INTERVAL,
    SW_MATRIX_NAME,
    SW_GAP_OPEN,
    SW_GAP_EXTEND,
    MAX_FILTER_RANK,
)
from seqmapping.utils.logging import get_logger, start_run
from seqmapping.utils.paths import (
    PDB_UNIQUE_FASTA,
    UNIREF90_PROCESSED_FASTA,
    FILTERED_HITS_TSV,
    ALIGNED_FILTERED_HITS_TSV,
)


# -----------------------------------------------------------------------------
# Custom FASTA parser
# -----------------------------------------------------------------------------
def load_fasta_dict(fasta_file: str | Path) -> Dict[str, str]:
    # Load FASTA file into a dict: {sequence_id: sequence_string}
    # Supports plain or gzipped FASTA (*.gz).

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

                    if seq_count > 0 and seq_count % 5000000 == 0:
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


# -----------------------------------------------------------------------------
# Alignment utilities
# -----------------------------------------------------------------------------
def get_scoring_matrix(name: str):
    # Map config string -> parasail matrix object.

    key = name.strip().lower()
    if key == "blosum62":
        return parasail.blosum62
    if key == "blosum80":
        return parasail.blosum80
    if key == "blosum45":
        return parasail.blosum45
    if key == "pam30":
        return parasail.pam30
    raise ValueError(f"Unsupported SW_MATRIX_NAME={name!r}")


SCORING_MATRIX = get_scoring_matrix(SW_MATRIX_NAME)


def create_pretty_alignment(
    aligned_q: str,
    aligned_s: str,
    qseqid: str,
    sseqid: str,
    qstart: int,
    sstart: int,
    line_length: int = ALIGNMENT_LINE_WIDTH,
) -> str:

    # Create a human-readable alignment block with IDs and 1-based positions.
    # Returned string begins with a newline, so it can live in a TSV last column.

    if not aligned_q or not aligned_s:
        return ""

    match_line = "".join("|" if a == b else " " for a, b in zip(aligned_q, aligned_s))

    lines: List[str] = []
    total_length = len(aligned_q)

    for i in range(0, total_length, line_length):
        chunk_end = min(i + line_length, total_length)

        q_chunk = aligned_q[i:chunk_end]
        m_chunk = match_line[i:chunk_end]
        s_chunk = aligned_s[i:chunk_end]

        # gaps before this chunk
        q_gaps_before = aligned_q[:i].count("-")
        s_gaps_before = aligned_s[:i].count("-")

        q_pos_start = qstart + i - q_gaps_before
        s_pos_start = sstart + i - s_gaps_before

        q_pos_end = q_pos_start + len(q_chunk.replace("-", "")) - 1
        s_pos_end = s_pos_start + len(s_chunk.replace("-", "")) - 1

        q_line = f"{qseqid:15s} {q_pos_start:4d} {q_chunk} {q_pos_end:4d}"
        m_line = f"{'':20s} {m_chunk}"
        s_line = f"{sseqid:15s} {s_pos_start:4d} {s_chunk} {s_pos_end:4d}"

        lines.extend([q_line, m_line, s_line, ""])

    return "\n" + "\n".join(lines).strip()


def align_pair(
    qseqid: str, sseqid: str, query_seqs: Dict[str, str], subject_seqs: Dict[str, str]
) -> Optional[dict]:

    q_seq = query_seqs.get(qseqid)
    if q_seq is None:
        logging.getLogger(__name__).debug(f"Query sequence {qseqid} not found.")
        return None

    s_seq = subject_seqs.get(sseqid)
    if s_seq is None:
        logging.getLogger(__name__).debug(f"Subject sequence {sseqid} not found.")
        return None

    try:
        func = getattr(parasail, "sw_trace_striped_stats_16")
        result = func(q_seq, s_seq, SW_GAP_OPEN, SW_GAP_EXTEND, SCORING_MATRIX)
    except AttributeError:
        func = getattr(parasail, "sw_trace_striped_16")
        result = func(q_seq, s_seq, SW_GAP_OPEN, SW_GAP_EXTEND, SCORING_MATRIX)

    tb = getattr(result, "traceback", None)
    aligned_q = getattr(tb, "query", "") if tb else ""
    aligned_s = getattr(tb, "ref", "") if tb else ""

    stats = getattr(result, "stats", None)
    if stats:
        matches = int(stats.matches)
        mismatches = int(stats.mismatches)
        gapopens = int(stats.openings)
        aln_length = int(stats.length)
    else:
        aln_length = len(aligned_q)
        matches = sum(a == b for a, b in zip(aligned_q, aligned_s))
        mismatches = sum(a != b for a, b in zip(aligned_q, aligned_s))
        # NOTE: this is not true "gap opens" without stats, but keeps a fallback signal
        gapopens = aligned_q.count("-") + aligned_s.count("-")

    pident = (matches / aln_length * 100.0) if aln_length else 0.0
    bitscore = float(result.score)
    qcovhsp = (aln_length / len(q_seq) * 100.0) if q_seq else 0.0
    scovhsp = (aln_length / len(s_seq) * 100.0) if s_seq else 0.0

    q_aligned = len(aligned_q.replace("-", "")) if aligned_q else 0
    s_aligned = len(aligned_s.replace("-", "")) if aligned_s else 0

    end_q = getattr(result, "end_query", None)  # 0-based inclusive end
    end_s = getattr(result, "end_ref", None)

    if end_q is not None and q_aligned > 0:
        qstart = end_q - q_aligned + 2  # -> 1-based inclusive
        qend = end_q + 1
    else:
        qstart, qend = 1, len(q_seq)

    if end_s is not None and s_aligned > 0:
        sstart = end_s - s_aligned + 2
        send = end_s + 1
    else:
        sstart, send = 1, len(s_seq)

    qstart = max(1, qstart)
    qend = min(len(q_seq), qend)
    sstart = max(1, sstart)
    send = min(len(s_seq), send)

    pretty = create_pretty_alignment(aligned_q, aligned_s, qseqid, sseqid, qstart, sstart)

    return {
        "qseqid": qseqid,
        "sseqid": sseqid,
        "pident": round(pident, 4),
        "length": aln_length,
        "mismatch": mismatches,
        "gapopen": gapopens,
        "qstart": qstart,
        "qend": qend,
        "sstart": sstart,
        "send": send,
        "bitscore": round(bitscore, 1),
        "qcovhsp": round(qcovhsp, 2),
        "scovhsp": round(scovhsp, 2),
        "alignment_pretty": pretty,
    }


def format_output_row(row: dict) -> str:
    cols = [
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
    ]

    out = []
    for col in cols:
        v = row.get(col, "")
        if isinstance(v, float):
            if col == "pident":
                out.append(f"{v:.4f}")
            elif col in ("qcovhsp", "scovhsp"):
                out.append(f"{v:.2f}")
            elif col == "bitscore":
                out.append(f"{v:.1f}")
            else:
                out.append(str(v))
        else:
            out.append(str(v))

    tabular = "\t".join(out)
    pretty = row.get("alignment_pretty", "")
    return f"{tabular}\t{pretty}"


def process_pairs(
    pairs: List[Tuple[str, str]],
    query_seqs: Dict[str, str],
    subject_seqs: Dict[str, str],
    output_tsv: str | Path,
    logger: logging.Logger,
) -> None:
    total = len(pairs)
    start_time = time.time()
    logger.info(f"Starting alignment of {total:,} pairs (single-threaded)...")

    results: List[dict] = []
    missing_queries = set()
    missing_subjects = set()

    for i, (q, s) in enumerate(pairs, 1):
        res = align_pair(q, s, query_seqs, subject_seqs)
        if res:
            results.append(res)
        else:
            if q not in query_seqs:
                missing_queries.add(q)
            if s not in subject_seqs:
                missing_subjects.add(s)

        if ALIGNMENT_LOG_INTERVAL and (i % ALIGNMENT_LOG_INTERVAL == 0):
            elapsed = time.time() - start_time
            rate = i / elapsed if elapsed > 0 else 0.0
            logger.info(
                f"Processed {i:,}/{total:,} pairs ({i/total*100:.2f}%), "
                f"{rate:.1f} pairs/sec, elapsed {elapsed/60:.1f} min"
            )

    if missing_queries:
        logger.warning(f"Missing {len(missing_queries)} query sequences. First few: {list(missing_queries)[:5]}")
    if missing_subjects:
        logger.warning(f"Missing {len(missing_subjects)} subject sequences. First few: {list(missing_subjects)[:5]}")

    output_tsv = Path(output_tsv)
    output_tsv.parent.mkdir(parents=True, exist_ok=True)

    headers = [
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

    with output_tsv.open("w", newline="") as f:
        f.write("\t".join(headers) + "\n")
        for row in results:
            f.write(format_output_row(row) + "\n")

    logger.info(f"Saved {len(results):,} alignment results to {output_tsv}")


def main() -> None:
    log_stem = Path(__file__).stem
    logger = get_logger(__name__, log_stem=log_stem)
    start_run(logger, run_name=log_stem, argv=sys.argv)

    logger.info("--- Starting SW pairwise alignment ---")
    logger.info("Configuration:")
    logger.info(f"  MAX_FILTER_RANK: {MAX_FILTER_RANK}")
    logger.info(f"  ALIGNMENT_LINE_WIDTH: {ALIGNMENT_LINE_WIDTH}")
    logger.info(f"  ALIGNMENT_LOG_INTERVAL: {ALIGNMENT_LOG_INTERVAL}")
    logger.info(f"  SW_MATRIX_NAME: {SW_MATRIX_NAME}")
    logger.info(f"  SW_GAP_OPEN: {SW_GAP_OPEN}")
    logger.info(f"  SW_GAP_EXTEND: {SW_GAP_EXTEND}")
    logger.info(f"  FILTERED_HITS_TSV: {FILTERED_HITS_TSV}")
    logger.info(f"  PDB_UNIQUE_FASTA: {PDB_UNIQUE_FASTA}")
    logger.info(f"  UNIREF90_PROCESSED_FASTA: {UNIREF90_PROCESSED_FASTA}")
    logger.info(f"  ALIGNED_FILTERED_HITS_TSV: {ALIGNED_FILTERED_HITS_TSV}")

    try:
        cols_to_read = ["qseqid", "sseqid", "filter_rank"]
        df_hits = pd.read_csv(FILTERED_HITS_TSV, sep="\t", usecols=cols_to_read)

        df_filtered = df_hits[df_hits["filter_rank"] <= MAX_FILTER_RANK].copy()
        pairs = list(zip(df_filtered["qseqid"], df_filtered["sseqid"]))

        all_query_ids = set(df_filtered["qseqid"].unique())
        all_subject_ids = set(df_filtered["sseqid"].unique())

        logger.info(f"Total hits in file: {len(df_hits):,}")
        logger.info(f"Hits with filter_rank <= {MAX_FILTER_RANK}: {len(df_filtered):,}")
        logger.info(f"Unique query sequences: {len(all_query_ids):,}")
        logger.info(f"Unique subject sequences: {len(all_subject_ids):,}")
        logger.info(f"Total pairs to align: {len(pairs):,}")

        rank_counts = df_filtered["filter_rank"].value_counts().sort_index()
        for rank, count in rank_counts.items():
            logger.info(f"  filter_rank {rank}: {count:,} hits")

    except Exception as e:
        logger.error(f"Failed to read hits file: {e}")
        logger.exception(e)
        raise SystemExit(1)

    try:
        logger.info("Loading query sequences (PDB unique)...")
        query_seqs = load_fasta_dict(PDB_UNIQUE_FASTA)

        logger.info("Loading subject sequences (UniRef90 processed)...")
        subject_seqs = load_fasta_dict(UNIREF90_PROCESSED_FASTA)

    except Exception as e:
        logger.error(f"Failed to load FASTA sequences: {e}")
        logger.exception(e)
        raise SystemExit(1)

    all_query_ids = set(df_filtered["qseqid"].unique())
    all_subject_ids = set(df_filtered["sseqid"].unique())

    missing_queries = all_query_ids - set(query_seqs.keys())
    missing_subjects = all_subject_ids - set(subject_seqs.keys())

    if missing_queries:
        logger.warning(
            f"Missing {len(missing_queries)} query IDs in FASTA ({len(missing_queries)/max(1,len(all_query_ids))*100:.1f}%)"
        )
        logger.info(f"First 5 missing queries: {list(missing_queries)[:5]}")

    if missing_subjects:
        logger.warning(
            f"Missing {len(missing_subjects)} subject IDs in FASTA ({len(missing_subjects)/max(1,len(all_subject_ids))*100:.1f}%)"
        )
        logger.info(f"First 5 missing subjects: {list(missing_subjects)[:5]}")

    valid_pairs: List[Tuple[str, str]] = []
    removed = 0
    for q, s in pairs:
        if q in query_seqs and s in subject_seqs:
            valid_pairs.append((q, s))
        else:
            removed += 1

    if removed:
        logger.warning(f"Removed {removed:,} pairs due to missing sequences")
        logger.info(f"Valid pairs for alignment: {len(valid_pairs):,}")

    process_pairs(valid_pairs, query_seqs, subject_seqs, ALIGNED_FILTERED_HITS_TSV, logger)

    logger.info("--- Alignment complete ---")


if __name__ == "__main__":
    main()
