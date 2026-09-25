# SeqMapping: Protein Sequence Mapping Pipeline

A comprehensive bioinformatics pipeline for mapping PDB protein sequences to UniRef90 clusters using DIAMOND and Smith-Waterman alignment, with robust MongoDB integration and SIFTS-based (Structure Integration with Function, Taxonomy and Sequences resource) benchmarking.

The project is intended as a modern alternative to SIFTS-style sequence cross-referencing, providing transparent scoring, reproducibility, and scalable processing for tens of millions of alignments.

------------------------------------------------------------------------

## Key features

-   restart-safe workflow
-   scalable to millions of sequences
-   taxonomy-aware scoring
-   reproducible mapping
-   MongoDB storage
-   full benchmarking against SIFTS

------------------------------------------------------------------------

## High-level pipeline

```         
download databases
        ↓
preprocess UniRef90
        ↓
extract features (PDB + UniRef90 + SWISS-prot)
        ↓
DIAMOND homology search (2-pass)
        ↓
MongoDB storage (features + Diamond hits output)
        ↓
taxonomy-aware filtering
        ↓
Smith-Waterman alignment validation
        ↓
SIFTS comparison benchmark
```

------------------------------------------------------------------------

## Installation

### Prerequisites

-   Python 3.11+

-   MongoDB 4.4+

-   DIAMOND aligner 2.0+

-   650GB+ disk space (for full UniRef90 database)

### External dependency: DIAMOND

Install:

``` bash
sudo apt-get update
sudo apt-get install -y diamond-aligner
```

Verify:

``` bash
diamond --version
which diamond
```

------------------------------------------------------------------------

## MongoDB configuration

### Data-warehouse Mongo (RCSB metadata)

```         
export DW_MONGO_URI="mongodb://USER:PASSWORD@HOST:27017/

# Optional overrides
export DW_MONGO_DB="dw"
export DW_MONGO_COLLECTION="core_polymer_entity"
```

Used for organism and description metadata during PDB feature extraction.

### SeqMapping Mongo (results database)

Configured in:

```         
seqmapping/utils/config.py
```

Stores:

-   pdb_features
-   uniref90_features
-   search_out (DIAMOND hits)

------------------------------------------------------------------------

## Quickstart

After cloning the repository, and installing Python and DIAMOND:

``` bash
# 1. Create environment
python -m venv .venv
source .venv/bin/activate
pip install -U pip

# 2. Install project
pip install -e .

# 3. Configure RCSB metadata database
export DW_MONGO_URI="mongodb://USER:PASSWORD@HOST:27017/"

# 4. Download required biological databases
python -m seqmapping.download_databases &

# 5. Prepare UniRef database
python -m seqmapping.preprocess_uniref90 &

# 6. Extract sequence features
python -m seqmapping.extract_pdb_features &
python -m seqmapping.extract_uniref90_features &

# 7. Run homology search
python -m seqmapping.run_diamond_search &

# 8. Load results into MongoDB
python -m seqmapping.write_hits_to_mongodb &

# 9. Rank hits
python -m seqmapping.filter_hits &

# 10. Validate with Smith-Waterman alignment
python -m seqmapping.smith_waterman_align &

# 11. Benchmark vs SIFTS
python -m benchmarks.map_uniref_rcsb_id &   # Map PDB chains to UniProt IDs using SIFTS
python -m benchmarks.aggregate_searchOut &  # Aggregate DIAMOND search results per query
python -m benchmarks.combine_sifts_seqmapping_hits & # Combine SIFTS and seqmapping results
python -m benchmarks.analyze_sifts_seqmapping_hits & # Generate benchmark report
python -m benchmarks.validate_alignment_against_sifts & # Benchmark refined alignments
python -m benchmarks.analyze_align_comparison & # Generate benchmark report for refined alignments
python -m benchmarks.compare_seqmapping_with_emboss_water & # Benchmark pairs previously flagged as suspicious in the SeqMapping vs SIFTS comparison
```

At completion, the benchmark report will be produced in the `benchmarks/report` directory.

------------------------------------------------------------------------

## Repository structure

```         
src/seqmapping
├── download_databases.py
├── extract_pdb_features.py
├── extract_uniref90_features.py
├── preprocess_uniref90.py
├── run_diamond_search.py
├── filter_hits.py
├── smith_waterman_align.py
├── write_hits_to_mongodb.py
└── utils/
    ├── config.py
    ├── logging.py
    └── paths.py
```

Benchmarking tools live in:

```         
src/benchmarks/
```

------------------------------------------------------------------------

## Full runtime directory tree

```         
.
├── data
│   ├── raw
│   │   ├── pdb/
│   │   ├── uniprot/
│   │   ├── uniref90/
│   │   └── taxonomy/
│   ├── intermediate
│   │   ├── pdb_features.tsv
│   │   ├── uniref90_features.tsv
│   │   ├── pdb_unique_sequences.fasta
│   │   ├── pdb_duplicate_sequences.fasta
│   │   ├── uniref90_processed.fasta
│   │   └── *.done (restart markers)
│   └── hits
│       ├── diamond_hits_pass1.tsv
│       ├── diamond_hits_pass2.tsv
│       ├── filtered_hits.tsv
│       └── filtered_hits_aligned.tsv
├── logs
│   ├── seqmapping.download_databases.log
│   ├── seqmapping.extract_pdb_features.log
│   ├── seqmapping.extract_uniref90_features.log
│   ├── seqmapping.preprocess_uniref90.log
│   ├── run_diamond_pass1.log
│   ├── run_diamond_pass2.log
│   ├── seqmapping.filter_hits.log
│   ├── seqmapping.smith_waterman_align.log
│   └── seqmapping.write_hits_to_mongodb.log
└── benchmarks
    ├── data
    |   ├── aggregated_hits.tsv
    |   ├── pdb_chain_uniprot.tsv
    |   ├── sifts_mapping.tsv
    |   ├── sifts_seqmapping_hits_withLabels.tsv
    |   └── uniref90_clusters.tsv
    ├── logs
        ├── aggregate_searchOut.log
        ├── analyze_sifts_seqmapping_hits.log
        ├── combine_sifts_seqmapping_hits.log
        └── map_uniref_rcsb_id.log
    └── report
        └── seqmapping_benchmark_report.txt
```

------------------------------------------------------------------------

## Benchmarking vs SIFTS

The benchmark compares SeqMapping results with SIFTS expectations.

Each query is classified into one of six categories:

| Category | Meaning |
|--------------------|----------------------------------------------------|
| **no_diamond_no_sifts** | Neither DIAMOND nor SIFTS found a mapping (likely unmappable sequence). |
| **no_diamond_yes_sifts** | SIFTS expects a mapping but homology search found none (possible sequence differences or search sensitivity limits). |
| **yes_diamond_no_sifts** | Homology detected a mapping but SIFTS does not provide one (potentially missing annotation in SIFTS). |
| **yes_diamond_yes_sifts_FullMatches** | Homology search recovers all expected SIFTS targets (strong agreement). |
| **yes_diamond_yes_sifts_NoMatches** | Both methods map the sequence but to different proteins (true disagreement). |
| **yes_diamond_yes_sifts_PartialMatches** | Only some expected SIFTS mappings are recovered (partial agreement). |

------------------------------------------------------------------------

## Logging

All modules write to file-only logs:

```         
logs/<module>.log
```

This enables reproducible and restartable runs.

------------------------------------------------------------------------

## Typical runtime

| Step               | Time             |
|--------------------|------------------|
| feature extraction | hours            |
| DIAMOND search     | hours            |
| MongoDB load       | hours            |
| filtering          | minutes to hours |
| SW alignment       | minutes          |
| benchmarking       | minutes          |

(Highly dependent on CPU count and storage throughput.)

------------------------------------------------------------------------

## Citation

If you use this pipeline in a publication, please cite:

> SeqMapping: large-scale mapping of PDB protein sequences to UniProt using homology search, taxonomy-aware filtering, and alignment validation??
