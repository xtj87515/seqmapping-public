#!/usr/bin/env python3
from __future__ import annotations

"""
Key Features:
- Loads UniRef90 and PDB feature metadata from MongoDB.
- Fetches all unique qseqids via distinct "qseqid".
- Chunks qseqids for multiprocessing.
- Robust multiprocessing-safe logging:
  - All worker logs (including tracebacks) are routed to the parent via a Queue.
  - Parent writes to file AND prints to console.
- Worker exceptions are no longer swallowed silently; they are logged with full traceback.

Default scoring rules:
filter_score=pident+0.08×qcovhsp+taxid_match+name_match+parent_match+sprot_bonus+refprot_bonus
Where:
taxid_match = 2.0 if qseq_taxid == sseq_taxid
name_match = 1.5 if qseq_tax_name and sseq_tax_name match (species/genus)
parent_match = 1.0 if qseq_parent_taxid == sseq_parent_taxid
sprot_bonus = 2.0 if sseq_is_uniprot_sprot == 1
refprot_bonus = 0.5 if sseq_is_ref_proteome == 1
"""

import os
import sys
import time
import logging
import logging.handlers
import multiprocessing as mp
from multiprocessing import Pool, cpu_count

import pandas as pd
import pymongo

from seqmapping.utils.config import (
    # Mongo (seqmapping DB)
    MONGO_URI,
    MONGO_DB_NAME,
    # Filter collection names (as constants)
    FILTER_PDB_COLLECTION,
    FILTER_UNIREF_COLLECTION,
    FILTER_SEARCH_COLLECTION,
    # Parallel defaults
    FILTER_QSEQID_CHUNK_SIZE,
    FILTER_NUM_PROCESSES,
    # Scoring / selection
    FILTER_QCOV_ALPHA,
    TAXONOMY_MATCH_BONUS,
    ORGANISM_NAME_MATCH_BONUS,
    PARENT_TAXID_MATCH_BONUS,
    SWISSPROT_REFERENCE_BONUS,
    REFPROT_BONUS,
    MAX_HITS_PER_QUERY,
    # Logging
    FILTER_LOG_EVERY,
)
from seqmapping.utils.logging import start_run
from seqmapping.utils.paths import FILTERED_HITS_TSV, LOG_DIR, ensure_directories


# -----------------------------------------------------------------------------
# Multiprocessing-safe logging (QueueHandler + QueueListener)
# -----------------------------------------------------------------------------
def _log_formatter() -> logging.Formatter:
    return logging.Formatter(
        "%(asctime)s %(levelname)s %(processName)s %(name)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


def setup_parent_logging(log_file_path: str) -> tuple[logging.handlers.QueueListener, "mp.Queue"]:
    """
    Parent process:
      - Creates a multiprocessing queue
      - Starts a QueueListener that writes to BOTH file and console
    """
    # IMPORTANT: using mp.Manager().Queue for broad compatibility
    manager = mp.Manager()
    log_queue = manager.Queue(-1)

    fmt = _log_formatter()

    # File handler (listener writes to it)
    os.makedirs(os.path.dirname(log_file_path), exist_ok=True)
    fh = logging.FileHandler(log_file_path, encoding="utf-8")
    fh.setLevel(logging.INFO)
    fh.setFormatter(fmt)

    # Suppress console printout
    listener = logging.handlers.QueueListener(log_queue, fh, respect_handler_level=True)
    listener.start()

    return listener, log_queue


def setup_worker_logging(log_queue) -> None:
    """
    Worker process initializer:
      - Routes ALL logs to the parent's queue via QueueHandler
      - Captures unhandled exceptions in the worker
    """
    qh = logging.handlers.QueueHandler(log_queue)
    root = logging.getLogger()
    root.handlers = []  # avoid duplicates
    root.addHandler(qh)
    root.setLevel(logging.INFO)

    def _excepthook(exc_type, exc, tb):
        logging.getLogger("worker").error("Unhandled exception", exc_info=(exc_type, exc, tb))

    sys.excepthook = _excepthook


# Use one named logger
logger = logging.getLogger("filter_hits")


# -----------------------------------------------------------------------------
# Output schema
# -----------------------------------------------------------------------------
def get_output_columns() -> list[str]:
    return [
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
        "evalue",
        "bitscore",
        "qcovhsp",
        "scovhsp",
        "qseq_length",
        "qseq_taxid",
        "qseq_tax_name",
        "qseq_domain_name",
        "qseq_n_source_organisms",
        "qseq_parent_taxid",
        "qseq_grandparent_taxid",
        "qseq_source_type",
        "qseq_description",
        "sseq_taxid",
        "sseq_tax_name",
        "sseq_length",
        "sseq_is_uniprot_sprot",
        "sseq_is_ref_proteome",
        "sseq_parent_taxid",
        "sseq_grandparent_taxid",
        "filter_score",
        "filter_rank",
    ]


# -----------------------------------------------------------------------------
# Utilities
# -----------------------------------------------------------------------------
def normalize_sseqid(sseqid: object) -> object:
    # Convert values like "sp|Q0KIY5|MYG_KOGBR" → "Q0KIY5"
    if isinstance(sseqid, str) and "|" in sseqid:
        parts = sseqid.split("|")
        for part in parts:
            if part and part[0] in ("P", "Q", "O") and len(part) in (6, 10):
                return part
    return sseqid


def is_same_organism_vectorized(q_names: pd.Series, s_names: pd.Series) -> pd.Series:
    q_norm = q_names.fillna("").astype("string").str.lower().str.replace(r"[^a-z0-9\s;]", "", regex=True)
    s_norm = s_names.fillna("").astype("string").str.lower().str.replace(r"[^a-z0-9\s]", "", regex=True)
    q_gs = q_norm.str.split(";").str[0].str.split().str[:2].str.join(" ")
    s_gs = s_norm.str.split().str[:2].str.join(" ")
    return (q_gs == s_gs) & (q_gs != "")


def calculate_score(df: pd.DataFrame) -> pd.Series:
    score = pd.Series(0.0, index=df.index)

    score += df.get("pident", 0).fillna(0).astype(float)
    score += float(FILTER_QCOV_ALPHA) * df.get("qcovhsp", 0).fillna(0).astype(float)

    if "qseq_taxid" in df.columns and "sseq_taxid" in df.columns:
        score += (df["qseq_taxid"].astype("string") == df["sseq_taxid"].astype("string")).fillna(False) * float(
            TAXONOMY_MATCH_BONUS
        )

    if "qseq_tax_name" in df.columns and "sseq_tax_name" in df.columns:
        score += is_same_organism_vectorized(df["qseq_tax_name"], df["sseq_tax_name"]).fillna(False) * float(
            ORGANISM_NAME_MATCH_BONUS
        )

    if "qseq_parent_taxid" in df.columns and "sseq_parent_taxid" in df.columns:
        score += (
            (df["qseq_parent_taxid"].astype("string") == df["sseq_parent_taxid"].astype("string")).fillna(False)
            * float(PARENT_TAXID_MATCH_BONUS)
        )

    if "sseq_is_uniprot_sprot" in df.columns:
        score += (df["sseq_is_uniprot_sprot"].fillna(0).astype("Int64") == 1) * float(SWISSPROT_REFERENCE_BONUS)

    if "sseq_is_ref_proteome" in df.columns:
        score += (df["sseq_is_ref_proteome"].fillna(0).astype("Int64") == 1) * float(REFPROT_BONUS)

    return score.astype(float)


# -----------------------------------------------------------------------------
# Worker (one Mongo client per task)
# -----------------------------------------------------------------------------
def process_chunk(qseqid_chunk: list[str]) -> pd.DataFrame | None:
    client = None
    try:
        client = pymongo.MongoClient(
            MONGO_URI,
            retryWrites=True,
            retryReads=True,
            serverSelectionTimeoutMS=30_000,
            connectTimeoutMS=30_000,
            socketTimeoutMS=120_000,
        )
        db = client[MONGO_DB_NAME]

        search_docs = list(db[FILTER_SEARCH_COLLECTION].find({"qseqid": {"$in": qseqid_chunk}}, {"_id": 0}))
        if not search_docs:
            return None

        search_df = pd.DataFrame(search_docs)
        if search_df.empty:
            return None

        if "sseqid" in search_df.columns:
            search_df["sseqid"] = search_df["sseqid"].apply(normalize_sseqid)

        pdb_docs = list(db[FILTER_PDB_COLLECTION].find({"qseqid": {"$in": qseqid_chunk}}, {"_id": 0}))
        pdb_df = pd.DataFrame(pdb_docs)

        sseqids = search_df["sseqid"].dropna().unique().tolist() if "sseqid" in search_df.columns else []
        uniref_docs = list(db[FILTER_UNIREF_COLLECTION].find({"sseqid": {"$in": sseqids}}, {"_id": 0}))
        uniref_df = pd.DataFrame(uniref_docs)

        merged_df = search_df.merge(pdb_df, on="qseqid", how="left").merge(uniref_df, on="sseqid", how="left")
        if merged_df.empty:
            return None

        merged_df["filter_score"] = calculate_score(merged_df)

        n = int(MAX_HITS_PER_QUERY)
        top_df = (
            merged_df.sort_values(["qseqid", "filter_score"], ascending=[True, False], kind="mergesort")
            .groupby("qseqid", as_index=False, sort=False)
            .head(n)
        )

        top_df["filter_rank"] = top_df.groupby("qseqid")["filter_score"].rank(method="dense", ascending=False).astype(
            "Int64"
        )

        return top_df

    except Exception:
        logger.exception(
            "Worker failed. chunk_size=%d first_qseqid=%s",
            len(qseqid_chunk),
            qseqid_chunk[0] if qseqid_chunk else "NA",
        )
        return None

    finally:
        if client is not None:
            try:
                client.close()
            except Exception:
                logger.exception("Failed to close Mongo client in worker.")


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------
def main() -> None:
    ensure_directories()

    try:
        LOG_DIR.mkdir(parents=True, exist_ok=True)  # if LOG_DIR is a pathlib.Path
    except AttributeError:
        # If LOG_DIR is a string, ensure directory exists
        os.makedirs(str(LOG_DIR), exist_ok=True)

    log_file = str(LOG_DIR / "filter_hits.log") if hasattr(LOG_DIR, "__truediv__") else os.path.join(str(LOG_DIR), "filter_hits.log")

    # Start listener first
    listener, log_queue = setup_parent_logging(log_file)

    # IMPORTANT: Parent must ALSO emit logs into the queue
    qh = logging.handlers.QueueHandler(log_queue)
    root = logging.getLogger()
    root.handlers = []  # avoid duplicate prints
    root.addHandler(qh)
    root.setLevel(logging.INFO)

    # Make sure our named logger propagates to root
    logger.setLevel(logging.INFO)
    logger.propagate = True

    # Log unhandled exceptions in the parent too
    def _parent_excepthook(exc_type, exc, tb):
        logger.error("Unhandled exception in parent", exc_info=(exc_type, exc, tb))

    sys.excepthook = _parent_excepthook

    try:
        start_run(logger, "filter_hits_parallel", argv=sys.argv)

        out_path = FILTERED_HITS_TSV
        out_path.parent.mkdir(parents=True, exist_ok=True)

        # Always overwrite output file
        if out_path.exists():
            out_path.unlink()
            logger.info("Removed existing output: %s", out_path)

        # Fetch all unique qseqids
        client = pymongo.MongoClient(
            MONGO_URI,
            retryWrites=True,
            retryReads=True,
            serverSelectionTimeoutMS=30_000,
            connectTimeoutMS=30_000,
            socketTimeoutMS=120_000,
        )
        try:
            db = client[MONGO_DB_NAME]

            logger.info("Fetching all unique qseqids via distinct() ...")
            all_qseqids = db[FILTER_SEARCH_COLLECTION].distinct("qseqid")
            total = len(all_qseqids)
            logger.info("Found %s unique qseqids.", f"{total:,}")

        finally:
            client.close()

        # Write header
        cols = get_output_columns()
        with out_path.open("w", encoding="utf-8") as out_f:
            out_f.write("\t".join(cols) + "\n")

        # Chunk qseqids
        q_chunk = int(FILTER_QSEQID_CHUNK_SIZE)
        chunks = [all_qseqids[i : i + q_chunk] for i in range(0, total, q_chunk)]
        total_chunks = len(chunks)

        # Determine processes
        if int(FILTER_NUM_PROCESSES) <= 0:
            nprocs = cpu_count()
        else:
            nprocs = int(FILTER_NUM_PROCESSES)

        logger.info("Starting parallel processing with %d workers on %d chunks ...", nprocs, total_chunks)
        t0 = time.time()
        processed = 0

        with Pool(
            processes=nprocs,
            initializer=setup_worker_logging,
            initargs=(log_queue,),
            maxtasksperchild=200,
        ) as pool:
            for result in pool.imap_unordered(process_chunk, chunks):
                processed += 1
                if result is not None and not result.empty:
                    result = result.reindex(columns=cols)
                    result.to_csv(out_path, sep="\t", index=False, header=False, mode="a", na_rep="NA")

                log_every = int(FILTER_LOG_EVERY) if "FILTER_LOG_EVERY" in globals() else 10
                if processed % log_every == 0 or processed == total_chunks:
                    elapsed = time.time() - t0
                    logger.info(
                        "Processed %d/%d chunks (%.1f%%) - Elapsed: %.1f min",
                        processed,
                        total_chunks,
                        100.0 * processed / max(total_chunks, 1),
                        elapsed / 60.0,
                    )

        logger.info("Parallel processing complete.")
        logger.info("Output written to: %s", out_path)
        logger.info("Logs written to: %s", log_file)
        logger.info("filter_hits finished successfully.")

    finally:
        # Ensure listener stops and flushes
        try:
            listener.stop()
        except Exception:
            pass


if __name__ == "__main__":
    main()
