"""
Centralized logging utilities for the seqmapping pipeline.

Design goals:
- File-only logging (no console output)
- Consistent formatting across all modules
"""

from __future__ import annotations

import logging
import os
import socket
from datetime import datetime
from pathlib import Path
from typing import Optional

from seqmapping.utils.paths import LOG_DIR

# -----------------------------------------------------------------------------
# Formatting
# -----------------------------------------------------------------------------
LOG_FORMAT = "%(asctime)s - %(levelname)s - %(message)s"
DATE_FORMAT = "%Y-%m-%d %H:%M:%S"


# -----------------------------------------------------------------------------
# Core setup
# -----------------------------------------------------------------------------
def setup_logging(
    log_file: str | Path,
    level: int = logging.INFO,
    logger_name: str = "seqmapping",
) -> logging.Logger:
   
    log_path = Path(log_file)
    log_path.parent.mkdir(parents=True, exist_ok=True)

    logger = logging.getLogger(logger_name)
    logger.setLevel(level)

    # Prevent bubbling to root logger (which might print to console)
    logger.propagate = False

    # Clear handlers so repeated calls don't duplicate output
    logger.handlers.clear()

    formatter = logging.Formatter(fmt=LOG_FORMAT, datefmt=DATE_FORMAT)

    file_handler = logging.FileHandler(
        log_path,
        mode="a",            # append across runs
        encoding="utf-8",    # safe for unicode, regardless of terminal locale
        delay=True,
    )
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    # Prevent logging internals from printing handler errors to stderr
    logging.raiseExceptions = False

    return logger


def get_logger(
    module_name: str,
    log_stem: Optional[str] = None,
    level: int = logging.INFO,
) -> logging.Logger:
    """
    Return a logger for a module, writing to logs/<log_stem>.log.

    Examples:
      logger = get_logger(__name__)
        -> logs/seqmapping.download_databases.log

      logger = get_logger(__name__, "download_databases")
        -> logs/download_databases.log
    """
    if log_stem is None:
        log_stem = module_name

    log_file = LOG_DIR / f"{log_stem}.log"
    return setup_logging(log_file=log_file, level=level, logger_name=module_name)
  
# -----------------------------------------------------------------------------
# Benchmark logging
# -----------------------------------------------------------------------------
from seqmapping.utils.paths import BENCHMARKS_LOG_DIR

def get_benchmark_logger(
    script_file: str | Path,
    level: int = logging.INFO,
) -> logging.Logger:
    """
    Create a logger whose filename matches the script name,
    writing into benchmarks/logs/.

    Example:
        map_uniref_rcsb_id.py -> benchmarks/logs/map_uniref_rcsb_id.log
    """
    script_path = Path(script_file)
    log_stem = script_path.stem
    log_file = BENCHMARKS_LOG_DIR / f"{log_stem}.log"

    return setup_logging(log_file=log_file, level=level, logger_name=log_stem)



# -----------------------------------------------------------------------------
# Run separators
# -----------------------------------------------------------------------------
def start_run(
    logger: logging.Logger,
    run_name: Optional[str] = None,
    *,
    argv: Optional[list[str]] = None,
) -> None:
    # Insert a blank line + a "RUN START" header so multiple runs are seperated
    ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    host = socket.gethostname()
    pid = os.getpid()

    label = run_name or logger.name
    cmd = ""
    if argv:
        cmd = " | cmd=" + " ".join(argv)

    # Blank line to separate between re-runs
    logger.info("")
    logger.info(f"=== RUN START: {label} | time={ts} | host={host} | pid={pid}{cmd} ===")
