#!/usr/bin/env python3

# If the program crashed mid-load and you don’t want to start over, run without dropping:
# python -m seqmapping.write_hits_to_mongodb --collections uniref90_features --no-drop &

# If run with drop, checkpoints need to be reset:
# python -m seqmapping.write_hits_to_mongodb --collections uniref90_features --drop --reset-checkpoints &


"""
Main features:
- Loads seqmapping TSV outputs into MongoDB: pdb_features, uniref90_features,
  and DIAMOND search_out (pass1 + pass2) using a single FILE_CONFIG spec.
- Supports partial runs / restart: per-(collection,file) JSON checkpoints store
  rows_done + chunk_num so a crash can resume with --no-drop (no full reload).
- Reads TSVs in pandas chunks; for headerless DIAMOND TSVs, supplies explicit
  column names (DIAMOND_OUTFMT_FIELDS) and appends a "source" tag (pass1/pass2).
- Normalizes DIAMOND sseqid by stripping "UniRef90_" prefix before insertion.
- Optional checkpoint management: --no-checkpoints disables resume tracking;
  --reset-checkpoints deletes checkpoint files for selected collections.
- Creates configured indexes after each collection load (compound indexes supported).
- CLI controls: choose collections, drop vs no-drop, chunksize, batch size,
  retry budget, and console progress display (tqdm).

"""



from __future__ import annotations

import argparse
import json
import math
import sys
import time
from pathlib import Path
from typing import Any, Iterable, Optional

import pandas as pd
import pymongo
from pymongo.collection import Collection
from pymongo.errors import AutoReconnect, BulkWriteError, NetworkTimeout
from tqdm import tqdm

from seqmapping.utils.config import (
    LOG_CHUNK_INTERVAL,
    MONGO_BATCH_ORDERED,
    MONGO_CHUNKSIZE,
    MONGO_COLLECTIONS_TO_LOAD_DEFAULT,
    MONGO_CONNECT_TIMEOUT_MS,
    MONGO_DB_NAME,
    MONGO_INSERT_BATCH_SIZE,          # NEW
    MONGO_MAX_RETRIES,                # NEW
    MONGO_RETRY_BACKOFF_MAX_SEC,      # NEW
    MONGO_RETRY_BACKOFF_SEC,          # NEW
    MONGO_ENABLE_CHECKPOINTS,         # NEW
    MONGO_SERVER_SELECTION_TIMEOUT_MS,
    MONGO_SOCKET_TIMEOUT_MS,           
    MONGO_URI,
    MONGO_WAITQUEUE_TIMEOUT_MS,
    SUPPRESS_CONSOLE_PROGRESS,
    DIAMOND_OUTFMT_FIELDS,
)
from seqmapping.utils.logging import get_logger, start_run
from seqmapping.utils.paths import (
    DIAMOND_HITS_PASS1_TSV,
    DIAMOND_HITS_PASS2_TSV,
    MONGO_CHECKPOINT_DIR,             # NEW
    PDB_FEATURES_TSV,
    UNIREF90_FEATURES_TSV,
    ensure_directories,
)

# ----------------------------
# Collection/file configuration
# ----------------------------
# NOTE:
# - For DIAMOND hits TSVs, there is no header row; we supply DIAMOND_OUTFMT_FIELDS.
# - We add a "source" field ("pass1"/"pass2") to the search_out collection.

FILE_CONFIG: list[dict[str, Any]] = [
    {
        "collection_name": "pdb_features",
        "file_paths": [PDB_FEATURES_TSV],
        "has_header": True,
        "columns": [
            "qseqid",
            "qseq_length",
            "qseq_taxid",
            "qseq_tax_name",
            "qseq_domain_name",
            "qseq_n_source_organisms",
            "qseq_parent_taxid",
            "qseq_grandparent_taxid",
            "qseq_source_type",
            "qseq_description",
        ],
        "dtype": {
            "qseqid": "string",
            "qseq_length": "Int64",
            "qseq_taxid": "string",
            "qseq_tax_name": "string",
            "qseq_domain_name": "string",
            "qseq_n_source_organisms": "Int64",
            "qseq_parent_taxid": "string",
            "qseq_grandparent_taxid": "string",
            "qseq_source_type": "string",
            "qseq_description": "string",
        },
        "indexes": [
            [("qseqid", pymongo.ASCENDING)],
        ],
    },
    {
        "collection_name": "uniref90_features",
        "file_paths": [UNIREF90_FEATURES_TSV],
        "has_header": True,
        "columns": [
            "sseqid",
            "sseq_taxid",
            "sseq_tax_name",
            "sseq_length",
            "sseq_is_uniprot_sprot",
            "sseq_is_ref_proteome",
            "sseq_parent_taxid",
            "sseq_grandparent_taxid",
        ],
        "dtype": {
            "sseqid": "string",
            "sseq_taxid": "string",
            "sseq_tax_name": "string",
            "sseq_length": "Int64",
            "sseq_is_uniprot_sprot": "Int64",
            "sseq_is_ref_proteome": "Int64",
            "sseq_parent_taxid": "string",
            "sseq_grandparent_taxid": "string",
        },
        "indexes": [
            [("sseqid", pymongo.ASCENDING)],
        ],
    },
    {
        "collection_name": "search_out",
        "file_paths": [DIAMOND_HITS_PASS1_TSV, DIAMOND_HITS_PASS2_TSV],
        "has_header": False,
        "columns": DIAMOND_OUTFMT_FIELDS + ["source"],
        "dtype": {
            "qseqid": "string",
            "sseqid": "string",
            "pident": "float64",
            "length": "Int64",
            "mismatch": "Int64",
            "gapopen": "Int64",
            "qstart": "Int64",
            "qend": "Int64",
            "sstart": "Int64",
            "send": "Int64",
            "evalue": "float64",
            "bitscore": "float64",
            "qcovhsp": "float64",
            "scovhsp": "float64",
            "source": "string",
        },
        "indexes": [
            [("qseqid", pymongo.ASCENDING), ("sseqid", pymongo.ASCENDING)],
            [("qseqid", pymongo.ASCENDING)],
            [("sseqid", pymongo.ASCENDING)],
        ],
    },
]


# ----------------------------
# CLI
# ----------------------------
def parse_args(argv: list[str]) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Load seqmapping TSV outputs into MongoDB (robust + resumable).")

    p.add_argument(
        "--collections",
        nargs="*",
        default=None,
        help=(
            "Collections to load (subset of: pdb_features, uniref90_features, search_out). "
            "If omitted, uses MONGO_COLLECTIONS_TO_LOAD_DEFAULT from config.py. "
            "If set to empty list explicitly (i.e. `--collections` with no args), loads all."
        ),
    )
    p.add_argument(
        "--drop",
        action="store_true",
        default=False,
        help="Drop existing collections before loading (default: False).",
    )
    p.add_argument(
        "--no-drop",
        dest="drop",
        action="store_false",
        help="Do not drop collections before loading. Useful for resuming from checkpoints.",
    )
    p.add_argument("--mongo-uri", default=MONGO_URI)
    p.add_argument("--db-name", default=MONGO_DB_NAME)
    p.add_argument("--chunksize", type=int, default=MONGO_CHUNKSIZE)
    p.add_argument(
        "--insert-batch-size",
        type=int,
        default=MONGO_INSERT_BATCH_SIZE,
        help="Number of documents per insert_many() call (smaller is more robust).",
    )
    p.add_argument(
        "--max-retries",
        type=int,
        default=MONGO_MAX_RETRIES,
        help="Max retries for transient Mongo network errors (AutoReconnect/NetworkTimeout).",
    )
    p.add_argument(
        "--show-progress",
        action="store_true",
        default=not SUPPRESS_CONSOLE_PROGRESS,
        help="Show tqdm progress on console.",
    )
    p.add_argument(
        "--no-checkpoints",
        dest="checkpoints",
        action="store_false",
        default=MONGO_ENABLE_CHECKPOINTS,
        help="Disable checkpointing/resume.",
    )
    p.add_argument(
        "--reset-checkpoints",
        action="store_true",
        default=False,
        help="Delete checkpoint files for selected collections before loading.",
    )

    return p.parse_args(argv)


# ----------------------------
# Mongo connection
# ----------------------------
def connect_mongo(uri: str, logger) -> pymongo.MongoClient:
    logger.info(f"Connecting to MongoDB: {uri}")
    client = pymongo.MongoClient(
        uri,
        connectTimeoutMS=MONGO_CONNECT_TIMEOUT_MS,
        socketTimeoutMS=MONGO_SOCKET_TIMEOUT_MS,  
        serverSelectionTimeoutMS=MONGO_SERVER_SELECTION_TIMEOUT_MS,
        waitQueueTimeoutMS=MONGO_WAITQUEUE_TIMEOUT_MS,
        retryWrites=True,
        retryReads=True,
    )
    client.admin.command("ping")
    logger.info("MongoDB connection OK (ping succeeded).")
    return client


# ----------------------------
# Checkpoints
# ----------------------------
def checkpoint_path(collection_name: str, file_path: Path) -> Path:
    safe = file_path.name.replace("/", "_").replace(".", "_")
    return MONGO_CHECKPOINT_DIR / f"{collection_name}__{safe}.json"


def load_checkpoint(collection_name: str, file_path: Path, *, enabled: bool) -> Optional[dict[str, int]]:
    if not enabled:
        return None
    cp = checkpoint_path(collection_name, file_path)
    if not cp.exists():
        return None
    try:
        return json.loads(cp.read_text(encoding="utf-8"))
    except Exception:
        return None


def save_checkpoint(collection_name: str, file_path: Path, rows_done: int, chunk_num: int, *, enabled: bool) -> None:
    if not enabled:
        return
    MONGO_CHECKPOINT_DIR.mkdir(parents=True, exist_ok=True)
    cp = checkpoint_path(collection_name, file_path)
    cp.write_text(json.dumps({"rows_done": int(rows_done), "chunk_num": int(chunk_num)}, indent=2), encoding="utf-8")


def maybe_reset_checkpoints(selected_cfgs: list[dict[str, Any]], *, enabled: bool, logger) -> None:
    if not enabled:
        return
    if not MONGO_CHECKPOINT_DIR.exists():
        return
    for cfg in selected_cfgs:
        cname = cfg["collection_name"]
        for fp in cfg["file_paths"]:
            cp = checkpoint_path(cname, Path(fp))
            if cp.exists():
                cp.unlink()
                logger.info(f"Deleted checkpoint: {cp}")


# ----------------------------
# Data helpers
# ----------------------------
def normalize_sseqid(chunk: pd.DataFrame) -> pd.DataFrame:
    # strip UniRef90_ prefix if present.
    if "sseqid" in chunk.columns:
        chunk["sseqid"] = chunk["sseqid"].astype("string").str.replace("UniRef90_", "", regex=False)
    return chunk


def to_records(chunk: pd.DataFrame) -> list[dict[str, Any]]:
    chunk = chunk.where(pd.notnull(chunk), None)
    return chunk.to_dict("records")


def iter_tsv_chunks_resumable(
    file_path: Path,
    *,
    chunksize: int,
    columns: list[str],
    dtype: dict[str, Any],
    has_header: bool,
    resume_rows: int,
) -> Iterable[pd.DataFrame]:
    
    # Resume by skipping `resume_rows` data rows efficiently by advancing the file handle.
    # - If has_header=True: skip header line + resume_rows data lines.
    # - If has_header=False: skip resume_rows lines.
    
    if resume_rows < 0:
        resume_rows = 0

    if has_header:
        f = file_path.open("r", encoding="utf-8")
        _ = f.readline()  # consume header line
    
        for _ in range(resume_rows):
            if not f.readline():
                break
    
        # header already consumed -> force header=None and supply names
        return pd.read_csv(
            f,
            sep="\t",
            header=None,
            names=columns,
            dtype=dtype,
            chunksize=chunksize,
            low_memory=False,
        )


    # headerless (DIAMOND)
    read_cols = [c for c in columns if c != "source"]
    read_dtype = {k: v for k, v in dtype.items() if k != "source"}

    f = file_path.open("r", encoding="utf-8")
    for _ in range(resume_rows):
        if not f.readline():
            break

    return pd.read_csv(
        f,
        sep="\t",
        header=None,
        names=read_cols,
        dtype=read_dtype,
        chunksize=chunksize,
        low_memory=False,
    )


# ----------------------------
# Robust insert
# ----------------------------
def insert_many_with_retry(
    logger,
    collection: Collection,
    batch: list[dict[str, Any]],
    *,
    context: str,
    max_retries: int,
) -> int:
    
    # Returns number of docs considered inserted for accounting.
    # On BulkWriteError, counts successful inserts (len(batch) - writeErrors).
    # On transient errors, retries with exponential backoff.
    
    attempt = 0
    backoff = float(MONGO_RETRY_BACKOFF_SEC)

    while True:
        try:
            collection.insert_many(batch, ordered=MONGO_BATCH_ORDERED)
            return len(batch)
        except BulkWriteError as e:
            write_errors = e.details.get("writeErrors", []) if e.details else []
            n_failed = len(write_errors)
            n_ok = max(0, len(batch) - n_failed)
            logger.warning(f"{context}: BulkWriteError: {n_failed} write errors; continuing.")
            return n_ok
        except (AutoReconnect, NetworkTimeout) as e:
            attempt += 1
            if attempt > max_retries:
                logger.error(f"{context}: exceeded retries ({max_retries}) on {type(e).__name__}: {e}")
                raise
            sleep_s = min(backoff, float(MONGO_RETRY_BACKOFF_MAX_SEC))
            logger.warning(
                f"{context}: transient Mongo error {type(e).__name__}: {e}. "
                f"Retry {attempt}/{max_retries} in {sleep_s:.1f}s"
            )
            time.sleep(sleep_s)
            backoff *= 2.0


# ----------------------------
# Indexes
# ----------------------------
def create_indexes(logger, collection: Collection, index_specs: list[list[tuple[str, int]]]) -> None:
    existing = collection.index_information()
    for spec in index_specs:
        name = "_".join([f"{f}_{d}" for f, d in spec])
        if name in existing:
            logger.info(f"Index exists on {collection.name}: {name} (skipping)")
            continue
        try:
            collection.create_index(spec, name=name, background=True)
            logger.info(f"Created index on {collection.name}: {name}")
        except Exception as e:
            logger.warning(f"Failed to create index on {collection.name} ({name}): {e}")


# ----------------------------
# Load routine
# ----------------------------
def load_one_collection(
    logger,
    *,
    db,
    collection_name: str,
    file_paths: list[Path],
    chunksize: int,
    insert_batch_size: int,
    max_retries: int,
    columns: list[str],
    dtype: dict[str, Any],
    has_header: bool,
    index_specs: list[list[tuple[str, int]]],
    show_progress: bool,
    drop: bool,
    checkpoints: bool,
) -> None:
    collection: Collection = db[collection_name]

    if drop:
        logger.info(f"Dropping collection for clean reload: {collection_name}")
        collection.drop()

    total_inserted = 0
    t0 = time.time()

    for fp in file_paths:
        if not fp.exists():
            raise FileNotFoundError(f"Missing input file for {collection_name}: {fp}")

        # Determine source tag for DIAMOND outputs
        source_tag = None
        if collection_name == "search_out":
            if fp.name.endswith("pass1.tsv"):
                source_tag = "pass1"
            elif fp.name.endswith("pass2.tsv"):
                source_tag = "pass2"
            else:
                source_tag = fp.stem

        # Resume position
        cp = load_checkpoint(collection_name, fp, enabled=checkpoints)
        resume_rows = int(cp["rows_done"]) if cp and "rows_done" in cp else 0
        if resume_rows > 0:
            logger.info(f"Resuming {collection_name} from checkpoint: rows_done={resume_rows:,} file={fp.name}")

        logger.info(f"Loading file -> {collection_name}: {fp} (source={source_tag})")

        chunk_iter = iter_tsv_chunks_resumable(
            fp,
            chunksize=chunksize,
            columns=columns,
            dtype=dtype,
            has_header=has_header,
            resume_rows=resume_rows,
        )

        pbar = tqdm(
            desc=f"Loading {collection_name}",
            unit="rows",
            disable=not show_progress,
            file=sys.stdout,
        )

        chunk_num = 0
        last_logged = 0
        rows_done_this_file = 0

        for chunk in chunk_iter:
            chunk_num += 1

            if collection_name == "search_out":
                chunk = normalize_sseqid(chunk)
                chunk["source"] = source_tag

            records = to_records(chunk)
            if not records:
                continue

            # Insert in smaller batches
            inserted_this_chunk = 0
            n_batches = math.ceil(len(records) / insert_batch_size)
            for b in range(n_batches):
                start = b * insert_batch_size
                batch = records[start : start + insert_batch_size]
                ctx = f"{collection_name} file={fp.name} chunk={chunk_num} batch={b+1}/{n_batches}"
                inserted_this_chunk += insert_many_with_retry(
                    logger,
                    collection,
                    batch,
                    context=ctx,
                    max_retries=max_retries,
                )

            total_inserted += inserted_this_chunk
            rows_done_this_file += len(records)
            pbar.update(len(records))

            # Save checkpoint after each chunk
            save_checkpoint(
                collection_name,
                fp,
                rows_done=resume_rows + rows_done_this_file,
                chunk_num=chunk_num,
                enabled=checkpoints,
            )

            if (chunk_num == 1) or ((chunk_num - 1) % LOG_CHUNK_INTERVAL == 0):
                logger.info(
                    f"{collection_name}: processed chunk {chunk_num} "
                    f"(rows={len(records)}). rows_done={resume_rows + rows_done_this_file:,} "
                    f"total_inserted={total_inserted:,}"
                )
                last_logged = chunk_num

        pbar.close()

        if chunk_num == 0:
            logger.info(f"{collection_name}: file empty, no chunks processed: {fp}")
        elif chunk_num != last_logged:
            logger.info(
                f"{collection_name}: processed final chunk {chunk_num}. "
                f"rows_done={resume_rows + rows_done_this_file:,} total_inserted={total_inserted:,}"
            )

    logger.info(f"Creating indexes for {collection_name} ...")
    create_indexes(logger, collection, index_specs)

    final_count = collection.count_documents({})
    logger.info(
        f"Finished {collection_name}: final_count={final_count:,}, "
        f"inserted_this_run={total_inserted:,}, elapsed={time.time() - t0:.2f}s"
    )


def main() -> None:
    args = parse_args(sys.argv[1:])
    ensure_directories()
    if args.checkpoints:
        MONGO_CHECKPOINT_DIR.mkdir(parents=True, exist_ok=True)

    logger = get_logger(__name__, "write_hits_to_mongodb")
    start_run(logger, "write_hits_to_mongodb", argv=sys.argv)

    # Determine collections to load
    if args.collections is None:
        collections_to_load = list(MONGO_COLLECTIONS_TO_LOAD_DEFAULT)
    else:
        collections_to_load = args.collections  # empty list means load all

    if collections_to_load == []:
        selected = FILE_CONFIG
        logger.info("No collections specified -> loading ALL configured collections.")
    else:
        selected = [c for c in FILE_CONFIG if c["collection_name"] in set(collections_to_load)]
        if not selected:
            raise ValueError(
                f"No valid collections matched: {collections_to_load}. "
                f"Valid: {[c['collection_name'] for c in FILE_CONFIG]}"
            )
        logger.info(f"Loading selected collections: {[c['collection_name'] for c in selected]}")

    if args.reset_checkpoints:
        logger.info("Resetting checkpoints for selected collections ...")
        maybe_reset_checkpoints(selected, enabled=args.checkpoints, logger=logger)

    client = None
    try:
        client = connect_mongo(args.mongo_uri, logger)
        db = client[args.db_name]

        for cfg in selected:
            load_one_collection(
                logger,
                db=db,
                collection_name=cfg["collection_name"],
                file_paths=[Path(p) for p in cfg["file_paths"]],
                chunksize=args.chunksize,
                insert_batch_size=args.insert_batch_size,
                max_retries=args.max_retries,
                columns=cfg["columns"],
                dtype=cfg["dtype"],
                has_header=cfg["has_header"],
                index_specs=cfg["indexes"],
                show_progress=args.show_progress,
                drop=args.drop,
                checkpoints=args.checkpoints,
            )

        logger.info("MongoDB load completed successfully.")

    except Exception as e:
        logger.error(f"MongoDB load FAILED: {e}", exc_info=True)
        raise
    finally:
        if client is not None:
            client.close()
            logger.info("MongoDB connection closed.")


if __name__ == "__main__":
    main()
