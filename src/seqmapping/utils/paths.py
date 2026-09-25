"""
Centralized filesystem paths for the seqmapping project.
"""

from __future__ import annotations

from pathlib import Path

# -----------------------------------------------------------------------------
# Project root
# File: project-root/src/seqmapping/utils/paths.py
# parents[0] = utils
# parents[1] = seqmapping
# parents[2] = src
# parents[3] = project-root
# -----------------------------------------------------------------------------
PROJECT_ROOT = Path(__file__).resolve().parents[3]
SRC_DIR = PROJECT_ROOT / "src"

# -----------------------------------------------------------------------------
# Top-level directories
# Folder rules:
# - data/raw/          : downloaded, immutable data
# - data/intermediate/ : all derived non-hit artifacts (features, taxonomy, etc.)
# - data/hits/         : DIAMOND hits, filtered hits, aligned hits
# -----------------------------------------------------------------------------
DATA_DIR = PROJECT_ROOT / "data"
RAW_DIR = DATA_DIR / "raw"
INTERMEDIATE_DIR = DATA_DIR / "intermediate"
HITS_DIR = DATA_DIR / "hits"

LOG_DIR = PROJECT_ROOT / "logs"
WORKFLOWS_DIR = PROJECT_ROOT / "workflows"


# -----------------------------------------------------------------------------
# Raw data subdirectories
# -----------------------------------------------------------------------------
RAW_PDB_DIR = RAW_DIR / "pdb"
RAW_UNIPROT_DIR = RAW_DIR / "uniprot"
RAW_UNIREF90_DIR = RAW_DIR / "uniref90"
RAW_TAXONOMY_DIR = RAW_DIR / "taxonomy"

# -----------------------------------------------------------------------------
# Raw input files
# -----------------------------------------------------------------------------
PDB_FASTA = RAW_PDB_DIR / "pdb_protein_sequences.fasta"

UNIPROT_SPROT_FASTA = RAW_UNIPROT_DIR / "uniprot_sprot.fasta"
UNIPROT_SPROT_XML = RAW_UNIPROT_DIR / "uniprot_sprot.xml"
UNIPROT_REF_PROTEOMES = RAW_UNIPROT_DIR / "uniprot_ref_proteomes.tsv"

UNIREF90_FASTA = RAW_UNIREF90_DIR / "uniref90.fasta"
UNIREF90_XML = RAW_UNIREF90_DIR / "uniref90.xml"

# -----------------------------------------------------------------------------
# Intermediate files
# -----------------------------------------------------------------------------
# Taxonomy-derived
NCBI_PARENT_TAXIDS_TSV = INTERMEDIATE_DIR / "ncbi_parent_taxids.tsv"

# Feature tables
PDB_FEATURES_TSV = INTERMEDIATE_DIR / "pdb_features.tsv"
UNIREF90_FEATURES_TSV = INTERMEDIATE_DIR / "uniref90_features.tsv"

# FASTA preprocessing
PDB_UNIQUE_FASTA = INTERMEDIATE_DIR / "pdb_unique_sequences.fasta"
PDB_DUPLICATE_FASTA = INTERMEDIATE_DIR / "pdb_duplicate_sequences.fasta"

UNIREF90_PROCESSED_FASTA = INTERMEDIATE_DIR / "uniref90_processed.fasta"
ADDED_SPROT_IDS_TSV = INTERMEDIATE_DIR / "added_sprot_ids.tsv"

# DIAMOND search workflow intermediates
PDB_FILTERED_SEQS_PASS2_FASTA = INTERMEDIATE_DIR / "pdb_filtered_seqs_pass2.fasta"
PDB_NOHIT_MARKER = INTERMEDIATE_DIR / "pdb_nohit.done"

# DIAMOND DB prefix (DIAMOND will create <prefix>.dmnd)
UNIREF90_DIAMOND_DB_PREFIX = INTERMEDIATE_DIR / "uniref90_processed"

# --- Step markers (.done) for restart/resume ---
FILTER_PASS2_DONE = INTERMEDIATE_DIR / "filter_pass2.done"
DIAMOND_MAKEDB_DONE = INTERMEDIATE_DIR / "diamond_makedb.done"
DIAMOND_PASS1_DONE = INTERMEDIATE_DIR / "diamond_pass1.done"
APPEND_NOHITS_DONE = INTERMEDIATE_DIR / "append_nohits.done"
DIAMOND_PASS2_DONE = INTERMEDIATE_DIR / "diamond_pass2.done"

# --- MongoDB loader checkpoints ---
MONGO_CHECKPOINT_DIR = INTERMEDIATE_DIR / "mongo_checkpoints"

# -----------------------------------------------------------------------------
# Hit profiles
# -----------------------------------------------------------------------------
DIAMOND_HITS_PASS1_TSV = HITS_DIR / "diamond_hits_pass1.tsv"
DIAMOND_HITS_PASS2_TSV = HITS_DIR / "diamond_hits_pass2.tsv"

FILTERED_HITS_TSV = HITS_DIR / "filtered_hits.tsv"
ALIGNED_FILTERED_HITS_TSV = HITS_DIR / "filtered_hits_aligned.tsv"


# -----------------------------------------------------------------------------
# Benchmarks
# - src/benchmarks/  : benchmark code 
# - benchmarks/*     : benchmark artifacts (data/logs)
# -----------------------------------------------------------------------------
BENCHMARKS_CODE_DIR = SRC_DIR / "benchmarks"          # benchmarks src
BENCHMARKS_DIR = PROJECT_ROOT / "benchmarks"          # artifacts root

BENCHMARKS_DATA_DIR = BENCHMARKS_DIR / "data"
BENCHMARKS_LOG_DIR = BENCHMARKS_DIR / "logs"
BENCHMARKS_REPORT_DIR = BENCHMARKS_DIR / "report"

UNIREF90_CLUSTERS_TSV = BENCHMARKS_DATA_DIR / "uniref90_clusters.tsv"
PDB_CHAIN_UNIPROT_TSV = BENCHMARKS_DATA_DIR / "pdb_chain_uniprot.tsv"
SIFTS_MAPPING_TSV = BENCHMARKS_DATA_DIR / "sifts_mapping.tsv"
AGGREGATED_HITS_TSV = BENCHMARKS_DATA_DIR / "aggregated_hits.tsv"
SIFTS_SEQMAPPING_HITS_TSV = BENCHMARKS_DATA_DIR / "sifts_seqmapping_hits_withLabels.tsv"
SEQMAPPING_BENCHMARK_REPORT = BENCHMARKS_REPORT_DIR / "seqmapping_benchmark_report.txt"
ALIGN_COMPARISON_TSV = BENCHMARKS_DATA_DIR / "sifts_seqmapping_align_comparison.tsv"
ALIGN_COMPARISON_REPORT = BENCHMARKS_REPORT_DIR / "alignment_comparison_report.txt"

SIFTS_SEQMAPPING_EMBOSS_TSV = BENCHMARKS_DATA_DIR / "sifts_seqmapping_emboss_comparison.tsv"
SIFTS_SEQMAPPING_EMBOSS_MD = BENCHMARKS_DATA_DIR / "sifts_seqmapping_emboss_comparison.md"
SIFTS_SEQMAPPING_EMBOSS_REPORT = BENCHMARKS_REPORT_DIR / "sifts_seqmapping_emboss_comparison.txt"


# -----------------------------------------------------------------------------
# Directory initialization
# -----------------------------------------------------------------------------
def ensure_directories() -> None:
    # Create required directory structure. Safe to call multiple times.
    for path in (
        # main pipeline dirs
        RAW_DIR,
        INTERMEDIATE_DIR,
        MONGO_CHECKPOINT_DIR,
        HITS_DIR,
        LOG_DIR,
        RAW_PDB_DIR,
        RAW_UNIPROT_DIR,
        RAW_UNIREF90_DIR,
        RAW_TAXONOMY_DIR,
        # benchmark dirs
        BENCHMARKS_DIR,
        BENCHMARKS_DATA_DIR,
        BENCHMARKS_LOG_DIR,
        BENCHMARKS_REPORT_DIR,
    ):
        path.mkdir(parents=True, exist_ok=True)
        
        
def ensure_benchmark_directories() -> None:
    BENCHMARKS_DIR.mkdir(parents=True, exist_ok=True)
    BENCHMARKS_DATA_DIR.mkdir(parents=True, exist_ok=True)
    BENCHMARKS_LOG_DIR.mkdir(parents=True, exist_ok=True)
    BENCHMARKS_REPORT_DIR.mkdir(parents=True, exist_ok=True)

