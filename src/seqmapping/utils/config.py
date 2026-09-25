import os

# -----------------------------------------------------------------------------
# MongoDB import
# -----------------------------------------------------------------------------
DW_MONGO_URI = os.environ.get("DW_MONGO_URI", "").strip()
DW_MONGO_DB = os.environ.get("DW_MONGO_DB", "dw").strip()
DW_MONGO_COLLECTION = os.environ.get("DW_MONGO_COLLECTION", "core_polymer_entity").strip()


def validate_dw_mongo_config() -> None:
    """
    Ensure DW MongoDB configuration is available.
    """
    if not DW_MONGO_URI:
        raise RuntimeError(
            "Missing DW MongoDB configuration.\n\n"
            "Set the following environment variable:\n"
            "  DW_MONGO_URI=mongodb://USER:PASSWORD@HOST:27017/?options\n\n"
            "Optional:\n"
            "  DW_MONGO_DB (default: dw)\n"
            "  DW_MONGO_COLLECTION (default: core_polymer_entity)"
        )


# -----------------------------------------------------------------------------
# General pipeline 
# -----------------------------------------------------------------------------
DEFAULT_N_THREADS = 10
RANDOM_SEED = 42


# -----------------------------------------------------------------------------
# Feature extraction 
# -----------------------------------------------------------------------------
WRITE_UNIREF90_CLUSTERS = True


# -----------------------------------------------------------------------------
# DIAMOND search 
# -----------------------------------------------------------------------------
# Path to DIAMOND binary.
# Default: "diamond" (resolved from PATH, e.g. system install via apt).
# Can be set to an absolute path for locked-down environments.
DIAMOND_BIN = "diamond"

DIAMOND_THREADS = 16

DIAMOND_OUTFMT_FIELDS = [
    "qseqid", "sseqid", "pident", "length", "mismatch",
    "gapopen", "qstart", "qend", "sstart", "send", "evalue",
    "bitscore", "qcovhsp", "scovhsp",
]

# Pass 1 
DIAMOND_PASS1 = {
    "sensitivity_flag": "--mid-sensitive",
    "max_target_seqs": 5000,

    # Formerly implicit defaults, now explicit:
    "evalue": 0.001,
    "masking": 1,            # default behavior: tantan repeat masking enabled
    "matrix": "BLOSUM62",    # default matrix
    "max_hsps": 1,           # default policy: top HSP only (explicitly set to 1)
}

# Pass 2
DIAMOND_PASS2 = {
    "sensitivity_flag": "--very-sensitive",
    "max_target_seqs": 30000,

    "evalue": 100,
    "masking": 0,
    "matrix": "PAM30",
    "max_hsps": 0,           # report all alternative HSPs
}

MAX_ORGANISMS = 1

# -----------------------------------------------------------------------------
# MongoDB export
# -----------------------------------------------------------------------------
MONGO_URI = "mongodb://128.6.159.216:27017/"
MONGO_DB_NAME = "seqmappingV6"


# Default collections to load if user doesn't specify via CLI.
# Set to [] to load ALL configured collections.
MONGO_COLLECTIONS_TO_LOAD_DEFAULT = ["uniref90_features", "pdb_features", "search_out"]

# Insert tuning
MONGO_CHUNKSIZE = 500000
MONGO_BATCH_ORDERED = False

# Connection timeouts (ms)
MONGO_CONNECT_TIMEOUT_MS = 3600000
MONGO_SOCKET_TIMEOUT_MS = 3600000
MONGO_SERVER_SELECTION_TIMEOUT_MS = 3600000
MONGO_WAITQUEUE_TIMEOUT_MS = 3600000

# Console progress
SUPPRESS_CONSOLE_PROGRESS = True
LOG_CHUNK_INTERVAL = 10


# MongoDB load robustness
"""
Controls how MongoDB bulk loads are made robust and restartable by inserting data in smaller batches, 
retrying transient connection failures with exponential backoff, and recording progress checkpoints 
so long-running imports can safely resume instead of restarting after interruptions.
"""
MONGO_INSERT_BATCH_SIZE = 10000        # docs per insert_many call
MONGO_MAX_RETRIES = 8                  # retries on AutoReconnect/NetworkTimeout
MONGO_RETRY_BACKOFF_SEC = 2.0          # base backoff
MONGO_RETRY_BACKOFF_MAX_SEC = 60.0     # cap
MONGO_ENABLE_CHECKPOINTS = True


# -----------------------------------------------------------------------------
# Filtering thresholds
# -----------------------------------------------------------------------------
MIN_QUERY_LENGTH = 30

# How many hits to keep per qseqid after scoring
MAX_HITS_PER_QUERY = 5

FILTER_QCOV_ALPHA = 0.08
TAXONOMY_MATCH_BONUS = 2.0
ORGANISM_NAME_MATCH_BONUS = 1.5
PARENT_TAXID_MATCH_BONUS = 1.0
SWISSPROT_REFERENCE_BONUS = 2.0   
REFPROT_BONUS = 0.5               

# Modify FILTER_QSEQID_CHUNK_SIZE and FILTER_NUM_PROCESSES to less, to avoid OOM issues
FILTER_QSEQID_CHUNK_SIZE = 180
FILTER_NUM_PROCESSES = 14   # 0 => use cpu_count()
FILTER_LOG_EVERY = 20

# Mongo collection names used by filter stage
FILTER_PDB_COLLECTION = "pdb_features"
FILTER_UNIREF_COLLECTION = "uniref90_features"
FILTER_SEARCH_COLLECTION = "search_out"

# -----------------------------------------------------------------------------
# Smith Waterman alignment pretty-print output 
# -----------------------------------------------------------------------------
ALIGNMENT_LINE_WIDTH = 100     # default 
ALIGNMENT_LOG_INTERVAL = 100000 
SW_MATRIX_NAME = "blosum62"    # keep as string; script maps to parasail matrix
SW_GAP_OPEN = 10
SW_GAP_EXTEND = 1
MAX_FILTER_RANK = 5  # Only process hits with filter_rank <= this value (1..5)


# -----------------------------------------------------------------------------
# Benchmarks
# -----------------------------------------------------------------------------
SIFTS_URL = "https://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/tsv/pdb_chain_uniprot.tsv.gz"
