#!/usr/bin/env python3
"""
PDB feature extraction pipeline:
  
1) Read raw PDB FASTA (data/raw/pdb/pdb_protein_sequences.fasta)
2) Write intermediate FASTAs:
   - data/intermediate/pdb_unique_sequences.fasta
   - data/intermediate/pdb_duplicate_sequences.fasta
3) Query MongoDB for organism + description metadata
4) Right-join representative sequence lengths with MongoDB data
5) Write intermediate features TSV:
   - data/intermediate/pdb_features.tsv

"""

from __future__ import annotations

import csv
import sys
import time
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, Iterable, Tuple

import pymongo
from Bio import SeqIO

from seqmapping.utils.config import (
    MIN_QUERY_LENGTH,
    MAX_ORGANISMS,
    DW_MONGO_URI,
    DW_MONGO_DB,
    DW_MONGO_COLLECTION,
    validate_dw_mongo_config,
)
from seqmapping.utils.logging import get_logger, start_run
from seqmapping.utils.paths import (
    PDB_FASTA,
    PDB_UNIQUE_FASTA,
    PDB_DUPLICATE_FASTA,
    PDB_FEATURES_TSV,
    NCBI_PARENT_TAXIDS_TSV,
    ensure_directories,
)

MULTI_VALUE_SEPARATOR = ";"


# -----------------------------------------------------------------------------
# Step 1: FASTA processing
# -----------------------------------------------------------------------------
def process_fasta(input_file: Path, unique_output: Path, duplicate_output: Path, logger) -> Dict[str, int]:
    # Extract representative unique sequences and duplicate groups.
    # Returns: lengths_map: {representative_qseqid: sequence_length}
    
    t0 = time.time()
    logger.info(f"Processing FASTA: {input_file}")

    if not input_file.exists():
        raise FileNotFoundError(f"Input FASTA file not found: {input_file}")

    unique_output.parent.mkdir(parents=True, exist_ok=True)
    duplicate_output.parent.mkdir(parents=True, exist_ok=True)

    sequence_groups: dict[str, list[str]] = defaultdict(list)
    total = 0
    kept = 0

    with input_file.open("r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            total += 1
            seq = str(record.seq)
            if len(seq) < MIN_QUERY_LENGTH:
                continue
            sequence_groups[seq].append(record.id)
            kept += 1

    logger.info(
        f"Total records: {total}; kept (len >= {MIN_QUERY_LENGTH}): {kept}; "
        f"unique sequences: {len(sequence_groups)}"
    )

    lengths_map: dict[str, int] = {}
    unique_count = 0
    duplicate_group_count = 0

    with unique_output.open("w") as unique_out, duplicate_output.open("w") as dup_out:
        for seq_content, ids_list in sequence_groups.items():
            length = len(seq_content)

            if len(ids_list) == 1:
                qseqid = ids_list[0]
                unique_out.write(f">{qseqid}\n{seq_content}\n")
                lengths_map[qseqid] = length
                unique_count += 1
                continue

            representative_id = ids_list[0]
            unique_out.write(f">{representative_id}\n{seq_content}\n")
            lengths_map[representative_id] = length
            unique_count += 1

            duplicate_header = "|".join(ids_list)
            dup_out.write(f">{duplicate_header}\n{seq_content}\n")
            duplicate_group_count += 1

    logger.info(f"Unique representatives written: {unique_count}")
    logger.info(f"Duplicate groups written: {duplicate_group_count}")
    logger.info(f"FASTA processing completed in {time.time() - t0:.2f} seconds")

    return lengths_map


# -----------------------------------------------------------------------------
# Step 2: taxonomy
# -----------------------------------------------------------------------------
def load_parent_taxids(file_path: Path, logger) -> Dict[str, Dict[str, str]]:
    # Load NCBI parent & grandparent taxids into dict: taxid -> {parent, grandparent}.
    logger.info(f"Loading parent & grandparent taxids from {file_path}")

    if not file_path.exists():
        raise FileNotFoundError(
            f"{file_path} not found. Run download_databases first to generate it."
        )

    parent_map: dict[str, dict[str, str]] = {}
    with file_path.open("r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            taxid = row["taxid"].strip()
            parent_map[taxid] = {
                "parent": row["parent_taxid"].strip(),
                "grandparent": row["grandparent_taxid"].strip(),
            }

    logger.info(f"Loaded {len(parent_map)} taxid mappings")
    return parent_map


# -----------------------------------------------------------------------------
# Step 3: Mongo processing
# -----------------------------------------------------------------------------
def iter_mongo_docs(logger) -> Iterable[Dict[str, Any]]:
    # Yield MongoDB docs needed for this pipeline.
    validate_dw_mongo_config()

    logger.info(
        f"Connecting to MongoDB: db={DW_MONGO_DB} collection={DW_MONGO_COLLECTION}"
    )

    client = pymongo.MongoClient(DW_MONGO_URI)
    collection = client[DW_MONGO_DB][DW_MONGO_COLLECTION]

    pipeline = [
        {"$match": {"rcsb_entity_source_organism.provenance_source": "Primary Data"}},
        {
            "$project": {
                "_id": 0,
                "rcsb_id": 1,
                "rcsb_entity_source_organism.ncbi_taxonomy_id": 1,
                "rcsb_entity_source_organism.ncbi_scientific_name": 1,
                "rcsb_entity_source_organism.ncbi_parent_scientific_name": 1,
                "rcsb_entity_source_organism.source_type": 1,
                "rcsb_polymer_entity.pdbx_description": 1,
            }
        },
    ]

    yield from collection.aggregate(pipeline, allowDiskUse=True)


def process_data(
    mongo_data: Iterable[Dict[str, Any]],
    parent_taxid_map: Dict[str, Dict[str, str]],
    logger,
) -> Dict[str, Dict[str, Any]]:
    # Process MongoDB documents into per-rcsb_id sets.
    data: dict[str, dict[str, Any]] = defaultdict(
        lambda: {
            "tax_ids": set(),
            "sci_names": set(),
            "parent_sci_names": set(),
            "qseq_parent_taxids": set(),
            "qseq_grandparent_taxids": set(),
            "source_type": set(),
            "descriptions": set(),
        }
    )

    doc_count = 0
    for doc in mongo_data:
        doc_count += 1
        rcsb_id = str(doc.get("rcsb_id", "")).strip()
        if not rcsb_id:
            continue

        organisms = doc.get("rcsb_entity_source_organism", []) or []
        if MAX_ORGANISMS > 0 and len(organisms) > MAX_ORGANISMS:
            organisms = organisms[:MAX_ORGANISMS]

        for org in organisms:
            tax_id = str(org.get("ncbi_taxonomy_id", "")).strip()
            sci_name = str(org.get("ncbi_scientific_name", "")).strip()
            parent_sci_name = str(org.get("ncbi_parent_scientific_name", "")).strip()
            source_type = str(org.get("source_type", "")).strip()

            if tax_id and tax_id.lower() not in {"na", "null", "none"}:
                data[rcsb_id]["tax_ids"].add(tax_id)

                parent_info = parent_taxid_map.get(tax_id)
                if parent_info:
                    p = parent_info.get("parent", "").strip()
                    gp = parent_info.get("grandparent", "").strip()
                    if p and p.lower() not in {"na", "null", "none"}:
                        data[rcsb_id]["qseq_parent_taxids"].add(p)
                    if gp and gp.lower() not in {"na", "null", "none"}:
                        data[rcsb_id]["qseq_grandparent_taxids"].add(gp)

            if sci_name and sci_name.lower() not in {"na", "null", "none"}:
                data[rcsb_id]["sci_names"].add(sci_name)

            if parent_sci_name and parent_sci_name.lower() not in {"na", "null", "none"}:
                data[rcsb_id]["parent_sci_names"].add(parent_sci_name)

            if source_type and source_type.lower() not in {"na", "null", "none"}:
                data[rcsb_id]["source_type"].add(source_type)

        desc = str(doc.get("rcsb_polymer_entity", {}).get("pdbx_description", "")).strip()
        if desc and desc.lower() not in {"na", "null", "none"}:
            data[rcsb_id]["descriptions"].add(desc)

    logger.info(f"Processed {doc_count} MongoDB docs")
    return data


# -----------------------------------------------------------------------------
# Step 4: write TSV (right join)
# -----------------------------------------------------------------------------
def write_joined_output(
    output_tsv: Path,
    data: Dict[str, Dict[str, Any]],
    lengths_map: Dict[str, int],
    logger,
) -> None:
    # Right-join lengths_map with Mongo data and write TSV.
    logger.info(f"Writing joined output TSV: {output_tsv}")
    output_tsv.parent.mkdir(parents=True, exist_ok=True)

    fieldnames = [
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
    ]

    def join_set(values: Any) -> str:
        if not values:
            return "NA"
        return MULTI_VALUE_SEPARATOR.join(sorted(values))

    with output_tsv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()

        for rcsb_id, qseq_length in lengths_map.items():
            org_data = data.get(rcsb_id)

            if not org_data:
                writer.writerow(
                    {
                        "qseqid": rcsb_id,
                        "qseq_length": qseq_length,
                        "qseq_taxid": "NA",
                        "qseq_tax_name": "NA",
                        "qseq_domain_name": "NA",
                        "qseq_n_source_organisms": "0",
                        "qseq_parent_taxid": "NA",
                        "qseq_grandparent_taxid": "NA",
                        "qseq_source_type": "NA",
                        "qseq_description": "NA",
                    }
                )
                continue

            tax_ids = org_data.get("tax_ids", set())
            writer.writerow(
                {
                    "qseqid": rcsb_id,
                    "qseq_length": qseq_length,
                    "qseq_taxid": join_set(tax_ids),
                    "qseq_tax_name": join_set(org_data.get("sci_names", set())),
                    "qseq_domain_name": join_set(org_data.get("parent_sci_names", set())),
                    "qseq_n_source_organisms": str(len(tax_ids)) if tax_ids else "0",
                    "qseq_parent_taxid": join_set(org_data.get("qseq_parent_taxids", set())),
                    "qseq_grandparent_taxid": join_set(org_data.get("qseq_grandparent_taxids", set())),
                    "qseq_source_type": join_set(org_data.get("source_type", set())),
                    "qseq_description": join_set(org_data.get("descriptions", set())),
                }
            )

    logger.info("Joined TSV written successfully")



def main() -> int:
    ensure_directories()

    logger = get_logger(__name__, "extract_pdb_features")
    start_run(logger, "extract_pdb_features", argv=sys.argv)

    t0 = time.time()
    logger.info("=== Pipeline started ===")
    logger.info(f"INPUT: {PDB_FASTA}")
    logger.info(f"INPUT: {NCBI_PARENT_TAXIDS_TSV}")
    logger.info(f"OUTPUT (intermediate): {PDB_UNIQUE_FASTA}")
    logger.info(f"OUTPUT (intermediate): {PDB_DUPLICATE_FASTA}")
    logger.info(f"OUTPUT (intermediate features): {PDB_FEATURES_TSV}")

    try:
        lengths_map = process_fasta(PDB_FASTA, PDB_UNIQUE_FASTA, PDB_DUPLICATE_FASTA, logger)
        parent_taxid_map = load_parent_taxids(NCBI_PARENT_TAXIDS_TSV, logger)
        mongo_data = iter_mongo_docs(logger)
        processed_data = process_data(mongo_data, parent_taxid_map, logger)
        write_joined_output(PDB_FEATURES_TSV, processed_data, lengths_map, logger)
        logger.info(f"=== Pipeline complete in {time.time() - t0:.2f} seconds ===")
        return 0
    except Exception as exc:
        logger.error(f"Pipeline failed: {exc}", exc_info=True)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
