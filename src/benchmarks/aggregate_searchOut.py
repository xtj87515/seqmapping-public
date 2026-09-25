#!/usr/bin/env python3

"""
Memory-safe aggregation of two Diamond search output files, deduplicating sseqids per qseqid.
Modified to clean sseqids:
- Extracts the central ID part from pipe-separated formats (e.g., 'Q0KIY5' from 'sp|Q0KIY5|MYG_KOGBR').
- Removes the 'UniRef90_' prefix from other IDs.
Supports test mode (process first N entries) and header control.
"""

import csv
import time
import os
import sys
from datetime import datetime

from seqmapping.utils.paths import (
    DIAMOND_HITS_PASS1_TSV,
    DIAMOND_HITS_PASS2_TSV,
    AGGREGATED_HITS_TSV,
    ensure_benchmark_directories,
)

from seqmapping.utils.logging import get_benchmark_logger, start_run


# Log name 
logger = get_benchmark_logger(__file__)

HAS_HEADER = False           # Set to True if the first row is a header in DIAMOND_HITS_PASS1_TSV and DIAMOND_HITS_PASS2_TSV
TEST_MODE = False            # Set to True to run a short test
TEST_LIMIT = 100             # Number of rows to process per file in test mode
PROGRESS_INTERVAL = 10_000_000

ensure_benchmark_directories()
start_run(logger, argv=sys.argv)

logger.info("Starting aggregation process for search_out files.")
logger.info(f"Test mode: {TEST_MODE} (limit {TEST_LIMIT} rows per file)")


# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------
def clean_sseqid(s):
    """
    Cleans the sseqid string.
    1. If it contains '|' (e.g., 'sp|Q0KIY5|MYG_KOGBR'), extract the middle part ('Q0KIY5').
    2. Otherwise, remove the 'UniRef90_' prefix.
    """
    if '|' in s:
        parts = s.split('|')
        if len(parts) >= 2:
            return parts[1]
        return s
    return s.replace("UniRef90_", "")


def process_files(output_file):
    start_time = time.time()
    qseqid_to_sseqids = {}
    total_lines = 0

    try:
        # ---------- PASS 1 ----------
        logger.info(f"Processing file 1/2: {DIAMOND_HITS_PASS1_TSV}")
        with open(DIAMOND_HITS_PASS1_TSV, 'r', encoding='utf-8') as fin:
            reader = csv.reader(fin, delimiter='\t')

            # Skip header if flagged
            if HAS_HEADER:
                next(reader, None)

            for row in reader:
                if len(row) < 2:
                    continue

                total_lines += 1
                qseqid = row[0]
                sseqid = row[1]

                if qseqid not in qseqid_to_sseqids:
                    qseqid_to_sseqids[qseqid] = set()
                qseqid_to_sseqids[qseqid].add(sseqid)

                if total_lines % PROGRESS_INTERVAL == 0:
                    logger.info(f"Processed {total_lines:,} lines so far...")

                # Stop early if in test mode
                if TEST_MODE and total_lines >= TEST_LIMIT:
                    logger.info(f"Test mode limit ({TEST_LIMIT}) reached. Stopping early.")
                    break

        # ---------- PASS 2 ----------
        if not (TEST_MODE and total_lines >= TEST_LIMIT):
            logger.info(f"Processing file 2/2: {DIAMOND_HITS_PASS2_TSV}")
            with open(DIAMOND_HITS_PASS2_TSV, 'r', encoding='utf-8') as fin:
                reader = csv.reader(fin, delimiter='\t')

                # Skip header if flagged
                if HAS_HEADER:
                    next(reader, None)

                for row in reader:
                    if len(row) < 2:
                        continue

                    total_lines += 1
                    qseqid = row[0]
                    sseqid = row[1]

                    if qseqid not in qseqid_to_sseqids:
                        qseqid_to_sseqids[qseqid] = set()
                    qseqid_to_sseqids[qseqid].add(sseqid)

                    if total_lines % PROGRESS_INTERVAL == 0:
                        logger.info(f"Processed {total_lines:,} lines so far...")

                    # Stop early if in test mode
                    if TEST_MODE and total_lines >= TEST_LIMIT:
                        logger.info(f"Test mode limit ({TEST_LIMIT}) reached. Stopping early.")
                        break

        logger.info(f"Finished reading all input files. Total lines processed: {total_lines:,}")
        logger.info(f"Unique qseqids found: {len(qseqid_to_sseqids):,}")

        # Write aggregated output
        with open(output_file, 'w', newline='', encoding='utf-8') as fout:
            writer = csv.writer(fout, delimiter='\t')
            writer.writerow(['qseqid', 'sseqids'])

            for qseqid, sseqids in qseqid_to_sseqids.items():
                cleaned_sseqids = [clean_sseqid(s) for s in sorted(sseqids)]
                writer.writerow([qseqid, ','.join(cleaned_sseqids)])

        elapsed = time.time() - start_time
        logger.info(f"Aggregation complete in {elapsed:.2f} seconds ({elapsed/60:.2f} minutes).")
        logger.info(f"Output written to {output_file}")

    except Exception as e:
        logger.exception(f"ERROR during aggregation: {e}")


if __name__ == "__main__":
    # Suppress console output
    try:
        sys.stdout = open(os.devnull, 'w')
        sys.stderr = open(os.devnull, 'w')
    except Exception:
        pass

    process_files(AGGREGATED_HITS_TSV)
