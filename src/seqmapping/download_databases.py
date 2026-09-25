#!/usr/bin/env python3
"""
Download and prepare external databases required by the pipeline.

Design goals: robust downloads (streaming), safe re-runs (skip if file exists).
"""

from __future__ import annotations

import gzip
import tarfile
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import requests

# Deterministic output locations via seqmapping.utils.paths.
from seqmapping.utils.logging import get_logger
from seqmapping.utils.paths import (
    RAW_PDB_DIR,
    RAW_UNIPROT_DIR,
    RAW_UNIREF90_DIR,
    RAW_TAXONOMY_DIR,
    INTERMEDIATE_DIR,
    ensure_directories,
)

CHUNK_SIZE = 8192


# -----------------------------------------------------------------------------
# Data model
# -----------------------------------------------------------------------------
@dataclass(frozen=True)
class DatabaseSpec:
    name: str
    url: str
    filename: str
    target_dir: Path
    gz: bool = False         # URL serves .gz that should be decompressed
    taxdump: bool = False    # special: NCBI taxdump.tar.gz extract + generate parent map


# -----------------------------------------------------------------------------
# Networking 
# -----------------------------------------------------------------------------
def make_session() -> requests.Session:
    # Create a requests session
    s = requests.Session()
    s.headers.update({"User-Agent": "seqmapping/0.1 (download_databases.py)"})
    return s


def get_remote_version(url: str, session: requests.Session) -> str:
    # Try to fetch version / timestamp info
    try:
        # UniProt: relnotes contains release line
        if "uniprot.org" in url and "current_release" in url:
            relnotes = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/relnotes.txt"
            r = session.get(relnotes, timeout=30)
            r.raise_for_status()
            line0 = r.text.splitlines()[0].strip() if r.text else ""
            return line0 or "Version info not available"

        # Generic: use HEAD Last-Modified
        r = session.head(url, timeout=30, allow_redirects=True)
        r.raise_for_status()
        return r.headers.get("Last-Modified", "Version info not available")
    except Exception:
        return "Version info not available"


def download_file(url: str, output_path: Path, session: requests.Session, logger) -> None:
    """Download a file with streaming I/O."""
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if output_path.exists() and output_path.stat().st_size > 0:
        logger.info(f"{output_path.name} already exists - skipping download")
        return

    logger.info(f"Remote version: {get_remote_version(url, session)}")

    with session.get(url, stream=True, timeout=120) as response:
        response.raise_for_status()
        with output_path.open("wb") as fh:
            for chunk in response.iter_content(chunk_size=CHUNK_SIZE):
                if chunk:
                    fh.write(chunk)

    logger.info(f"Downloaded {output_path.name}")


# -----------------------------------------------------------------------------
# Extraction 
# -----------------------------------------------------------------------------
def extract_gzip(gz_path: Path, output_path: Path, logger) -> None:
    """Extract .gz file using streaming I/O."""
    if output_path.exists() and output_path.stat().st_size > 0:
        logger.info(f"{output_path.name} already exists - skipping decompression")
        try:
            gz_path.unlink(missing_ok=True)
        except Exception:
            pass
        return

    output_path.parent.mkdir(parents=True, exist_ok=True)

    with gzip.open(gz_path, "rb") as f_in, output_path.open("wb") as f_out:
        for chunk in iter(lambda: f_in.read(CHUNK_SIZE), b""):
            f_out.write(chunk)

    try:
        gz_path.unlink(missing_ok=True)
    except Exception:
        pass

    logger.info(f"Decompressed {gz_path.name} -> {output_path.name}")


def extract_taxdump(tar_path: Path, output_dir: Path, logger) -> None:
    # Extract NCBI taxdump tarball.
    output_dir.mkdir(parents=True, exist_ok=True)

    tar_name = tar_path.name
    with tarfile.open(tar_path, "r:gz") as tar:
        tar.extractall(path=output_dir)

    try:
        tar_path.unlink(missing_ok=True)
    except Exception:
        pass

    logger.info(f"Extracted {tar_name} to {output_dir}")


def generate_ncbi_parent_taxids(taxdump_dir: Path, output_tsv: Path, logger) -> None:
    # Generate TSV with columns: taxid, parent_taxid, grandparent_taxid
    nodes_file = taxdump_dir / "nodes.dmp"
    if not nodes_file.exists():
        raise FileNotFoundError(f"{nodes_file} not found")

    taxid_to_parent: dict[str, str] = {}

    # nodes.dmp is ASCII-ish but be defensive
    with nodes_file.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            parts = line.split("\t|\t")
            if len(parts) >= 2:
                taxid_to_parent[parts[0].strip()] = parts[1].strip()

    output_tsv.parent.mkdir(parents=True, exist_ok=True)
    with output_tsv.open("w", encoding="utf-8") as out:
        out.write("taxid\tparent_taxid\tgrandparent_taxid\n")
        for taxid, parent in taxid_to_parent.items():
            grandparent = taxid_to_parent.get(parent, "")
            out.write(f"{taxid}\t{parent}\t{grandparent}\n")

    logger.info(f"Generated NCBI parent/grandparent taxids -> {output_tsv.name}")


# -----------------------------------------------------------------------------
# Pipeline
# -----------------------------------------------------------------------------
def iter_databases() -> Iterable[DatabaseSpec]:
    return [
        DatabaseSpec(
            name="PDB protein sequences",
            url="http://bl-east.rcsb.org/4-coastal/pdb_protein_sequence_all.fasta-A.gz",
            filename="pdb_protein_sequences.fasta",
            target_dir=RAW_PDB_DIR,
            gz=True,
        ),
        DatabaseSpec(
            name="UniProt Swiss-Prot FASTA",
            url="https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz",
            filename="uniprot_sprot.fasta",
            target_dir=RAW_UNIPROT_DIR,
            gz=True,
        ),
        DatabaseSpec(
            name="UniProt Swiss-Prot XML",
            url="https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.xml.gz",
            filename="uniprot_sprot.xml",
            target_dir=RAW_UNIPROT_DIR,
            gz=True,
        ),
        DatabaseSpec(
            name="NCBI taxdump",
            url="https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz",
            filename="taxdump.tar.gz",
            target_dir=RAW_TAXONOMY_DIR,
            taxdump=True,
        ),
        DatabaseSpec(
            name="UniProt reference proteomes",
            url=(
                "https://rest.uniprot.org/proteomes/stream?"
                "compressed=true&fields=upid,organism,organism_id,protein_count"
                "&format=tsv&query=(*)+AND+(proteome_type:1)"
            ),
            filename="uniprot_ref_proteomes.tsv",
            target_dir=RAW_UNIPROT_DIR,
            gz=True,
        ),
        DatabaseSpec(
            name="UniRef90 FASTA",
            url="https://ftp.uniprot.org/pub/databases/uniprot/current_release/uniref/uniref90/uniref90.fasta.gz",
            filename="uniref90.fasta",
            target_dir=RAW_UNIREF90_DIR,
            gz=True,
        ),
        DatabaseSpec(
            name="UniRef90 XML",
            url="https://ftp.uniprot.org/pub/databases/uniprot/current_release/uniref/uniref90/uniref90.xml.gz",
            filename="uniref90.xml",
            target_dir=RAW_UNIREF90_DIR,
            gz=True,
        ),
    ]


def process_database(db: DatabaseSpec, session: requests.Session, logger) -> None:
    # Download/extract one database entry.
    db.target_dir.mkdir(parents=True, exist_ok=True)

    output_path = db.target_dir / db.filename
    download_path = output_path.with_suffix(output_path.suffix + ".gz") if db.gz else output_path

    start = time.time()

    download_file(db.url, download_path, session, logger)

    if db.gz:
        extract_gzip(download_path, output_path, logger)

    if db.taxdump:
        # taxdump is a tar.gz already; extract into RAW_TAXONOMY_DIR
        extract_taxdump(download_path, RAW_TAXONOMY_DIR, logger)
        generate_ncbi_parent_taxids(
            RAW_TAXONOMY_DIR,
            INTERMEDIATE_DIR / "ncbi_parent_taxids.tsv",
            logger,
        )

    elapsed = time.time() - start
    logger.info(f"Finished {db.name} in {elapsed:.2f} seconds")


def main() -> int:
    # Returns: 0 on success, 1 if any database failed.
    ensure_directories()

    logger = get_logger(__name__, "download_databases")
    logger.info("=== Database download pipeline started ===")

    session = make_session()

    any_failed = False
    for db in iter_databases():
        try:
            process_database(db, session, logger)
        except Exception as exc:
            any_failed = True
            logger.error(f"Failed processing {db.name}: {exc}", exc_info=True)

    logger.info("=== All database operations completed ===")
    return 1 if any_failed else 0


if __name__ == "__main__":
    # Never print tracebacks to console; if something fatal happens, exit nonzero silently.
    try:
        raise SystemExit(main())
    except SystemExit:
        raise
    except Exception:
        # If logger isn't available for some reason, still do not print.
        try:
            logger = get_logger(__name__, "download_databases")
            logger.error("Fatal error in __main__", exc_info=True)
        except Exception:
            pass
        raise SystemExit(1)
