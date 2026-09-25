#!/usr/bin/env python3
"""
Workflow:
- Download SIFTS: If missing, download and unzip the SIFTS chain-to-UniProt mapping file.
- Extract chains from RCSB MongoDB: For each polymer entity (rcsb_id), retrieve its author chain IDs.
- Load SIFTS mapping: Read the SIFTS table containing: PDB chain → UniProt accession
- Join the two sources
   - Match MongoDB chains with SIFTS chains using: (PDB code, chain ID)
   - Convert: polymer entity → UniProt accession
- Output
   - Remove missing and duplicate pairs.
   - Write final mapping:
         qseqid   sseqid
   - Save as SIFTS_MAPPING_TSV.
"""


from __future__ import annotations

import sys
import time
import datetime
import gzip
import shutil
from pathlib import Path

import pandas as pd
import pymongo
import requests

from seqmapping.utils.config import (
    DW_MONGO_URI,
    DW_MONGO_DB,
    DW_MONGO_COLLECTION,
    validate_dw_mongo_config,
    SIFTS_URL,
)

from seqmapping.utils.paths import (
    PDB_CHAIN_UNIPROT_TSV,
    SIFTS_MAPPING_TSV,
    ensure_benchmark_directories,
)

from seqmapping.utils.logging import (
    get_benchmark_logger,
    start_run,
)


logger = get_benchmark_logger(__file__)


def download_sifts_file() -> None:
    """Download and extract SIFTS if missing."""
    if PDB_CHAIN_UNIPROT_TSV.exists():
        logger.info(f"SIFTS file already exists: {PDB_CHAIN_UNIPROT_TSV}")
        return

    gz_path = PDB_CHAIN_UNIPROT_TSV.with_suffix(".tsv.gz")

    try:
        logger.info(f"Downloading SIFTS file from {SIFTS_URL}")

        with requests.get(SIFTS_URL, stream=True, timeout=600) as r:
            r.raise_for_status()
            with open(gz_path, "wb") as f:
                for chunk in r.iter_content(chunk_size=8192):
                    if chunk:
                        f.write(chunk)

        logger.info("Download complete, extracting...")

        with gzip.open(gz_path, "rb") as f_in, open(PDB_CHAIN_UNIPROT_TSV, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)

        gz_path.unlink()
        logger.info("Extraction complete.")

    except Exception:
        logger.exception("Failed to download or extract SIFTS file")
        sys.exit(1)



def main() -> None:
    start_time = time.time()

    ensure_benchmark_directories()
    validate_dw_mongo_config()

    start_run(logger, argv=sys.argv)

    try:
        # Step 1
        download_sifts_file()

        # Step 2: MongoDB
        logger.info(f"Connecting to MongoDB at {DW_MONGO_URI}")
        client = pymongo.MongoClient(DW_MONGO_URI)
        db = client[DW_MONGO_DB]
        collection = db[DW_MONGO_COLLECTION]

        pipeline = [
            {"$project": {"_id": 0, "rcsb_polymer_entity_container_identifiers.auth_asym_ids": 1, "rcsb_id": 1}}
        ]

        rows = []
        for doc in collection.aggregate(pipeline):
            rcsb_id = doc.get("rcsb_id", "")
            auth_asym_ids = doc.get("rcsb_polymer_entity_container_identifiers", {}).get("auth_asym_ids", [])

            for author_id in auth_asym_ids:
                pdb_code = rcsb_id.split("_")[0].upper()
                rows.append((rcsb_id, pdb_code, author_id))

        df_auth = pd.DataFrame(rows, columns=["rcsb_id", "PDB", "author_id"])
        logger.info(f"Extracted {len(df_auth):,} rows from MongoDB")

        # Step 3
        logger.info("Loading SIFTS mapping")
        sifts_df = pd.read_csv(PDB_CHAIN_UNIPROT_TSV, sep="\t", comment="#", dtype=str)
        sifts_df["PDB"] = sifts_df["PDB"].str.upper()
        logger.info(f"SIFTS rows: {len(sifts_df):,}")

        # Step 4
        logger.info("Joining Mongo results with SIFTS")
        merged = (
            sifts_df
            .merge(df_auth, left_on=["PDB", "CHAIN"], right_on=["PDB", "author_id"], how="left")
            .rename(columns={"rcsb_id": "qseqid", "SP_PRIMARY": "sseqid"})
            [["qseqid", "sseqid"]]
            .dropna(subset=["qseqid", "sseqid"])
            .drop_duplicates(subset=["qseqid", "sseqid"])
        )

        logger.info(f"Final merged rows: {len(merged):,}")

        merged.to_csv(SIFTS_MAPPING_TSV, sep="\t", index=False)
        logger.info(f"Wrote output to {SIFTS_MAPPING_TSV}")

        elapsed = time.time() - start_time
        logger.info(f"Completed successfully in {elapsed:.2f} seconds")

    except Exception:
        logger.exception("Fatal error")
        sys.exit(1)


if __name__ == "__main__":
    main()
