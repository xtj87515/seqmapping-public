#!/usr/bin/env python3
import csv
import time
from collections import defaultdict

"""
Workflow (Build SIFTS + SeqMapping hit summary table)
- Load PDB_FEATURES_TSV:
    - Create a per-qseqid feature dictionary (length, taxonomy, source type, description, etc.).

- Load SIFTS_MAPPING_TSV:
    - Build sifts_map: qseqid -> list of UniProt accessions (SIFTS hits).

- Compute wanted UniProt IDs from SIFTS:
    - Collect the union of all SIFTS hit accessions across qseqids (wanted_ids).

- Load UniRef90 cluster mappings (filtered to wanted_ids only):
    - Identify which wanted_ids are UniRef90 representatives (uniref_direct).
    - Map wanted member IDs -> representative IDs (uniref_cluster).

- Load Diamond aggregated hits (AGGREGATED_HITS_TSV):
    - Build diamond_map: qseqid -> set of raw Diamond hit accessions.
    - Build diamond_counts: qseqid -> number of raw Diamond hits.

- Load curated “added Swiss-Prot IDs” (ADDED_SPROT_IDS_TSV):
    - Build added_sprot set for special-case inclusion / prioritization.

- Load filtered hit ranks (FILTERED_HITS_TSV) and aggregate:
    - For each qseqid, collect ranked hits with filter_rank <= 5.
    - Derive:
        - fh_top1[qseqid] = top1 filtered hit
        - fh_top5_list[qseqid] = list of top5 filtered hits
        - fh_top5_str[qseqid] = comma-joined top5 string (for output)

- For each qseqid in PDB_FEATURES_TSV, compute final SIFTS-vs-Diamond status and overlaps:
    - If qseqid has no SIFTS hits:
        - Label as no_diamond_no_sifts or yes_diamond_no_sifts.
    - Else:
        - Classify each SIFTS hit into:
            - direct UniRef90 rep (in_ur90),
            - UniRef90 cluster member (in_cluster -> store representative),
            - added Swiss-Prot (in_sprot).
        - Build combined_hits:
            - prioritize added Swiss-Prot if present; otherwise use UniRef90 reps + clusters.
        - Intersect combined_hits with Diamond raw hits:
            - produces sifts_hitsORclustersORsprotIDs_inRaw (and its count).
        - Assign stat_label:
            - no_diamond_yes_sifts / yes_diamond_yes_sifts_FullMatches /
              yes_diamond_yes_sifts_NoMatches / yes_diamond_yes_sifts_PartialMatches.

    - Integrate filtered-hit overlap:
        - filtered_hits column = top5 filtered hits (fh_top5_str)
        - siftsHits_in_top1filtered = top1 filtered hit if it is in the “SIFTS-derived Diamond hits”
        - siftsHits_in_top5filtered = intersection of top5 filtered hits with “SIFTS-derived Diamond hits”

- Write one consolidated output table (SIFTS_SEQMAPPING_HITS_TSV):
    - Single final TSV with SIFTS summary columns + selected PDB feature columns + filtered-hit overlap columns.
    - No intermediate TSVs are written.
"""


from seqmapping.utils.paths import (
    ensure_benchmark_directories,
    BENCHMARKS_DATA_DIR,
    BENCHMARKS_LOG_DIR,
    PDB_FEATURES_TSV,
    UNIREF90_CLUSTERS_TSV,
    ADDED_SPROT_IDS_TSV,
    AGGREGATED_HITS_TSV,
    SIFTS_MAPPING_TSV,
    FILTERED_HITS_TSV,
    SIFTS_SEQMAPPING_HITS_TSV,
)

from seqmapping.utils.logging import get_benchmark_logger, start_run


# ====================== config & logging
csv.field_size_limit(1_000_000_000)  # 1GB

def setup_logging():
    ensure_benchmark_directories()
    BENCHMARKS_DATA_DIR.mkdir(parents=True, exist_ok=True)
    BENCHMARKS_LOG_DIR.mkdir(parents=True, exist_ok=True)

    logger = get_benchmark_logger(__file__)
    start_run(logger, run_name="combine_sifts_searchOut + combine_sifts_filteredHits", argv=None)

    logger.info("Starting integrated SIFTS pipeline (single final output; no intermediates).")
    logger.info("Inputs:")
    logger.info("  PDB_FEATURES_TSV          = %s", str(PDB_FEATURES_TSV))
    logger.info("  UNIREF90_CLUSTERS_TSV     = %s", str(UNIREF90_CLUSTERS_TSV))
    logger.info("  AGGREGATED_HITS_TSV     = %s", str(AGGREGATED_HITS_TSV))
    logger.info("  SIFTS_MAPPING_TSV         = %s", str(SIFTS_MAPPING_TSV))
    logger.info("  ADDED_SPROT_IDS_TSV       = %s", str(ADDED_SPROT_IDS_TSV))
    logger.info("  FILTERED_HITS_TSV         = %s", str(FILTERED_HITS_TSV))
    logger.info("Output:")
    logger.info("  SIFTS_SEQMAPPING_HITS_TSV               = %s", str(SIFTS_SEQMAPPING_HITS_TSV))
    return logger


def step_start(logger, step_name: str):
    logger.info(f"{step_name}: start")
    return time.time()


def step_done(logger, step_name: str, t0: float, extra: str = ""):
    elapsed = time.time() - t0
    if extra:
        logger.info(f"{step_name}: done ({extra}, {elapsed:.2f}s)\n")
    else:
        logger.info(f"{step_name}: done ({elapsed:.2f}s)\n")


# ====================== loaders
def load_pdb_features(logger):
    step = "load_pdb_features"
    t0 = step_start(logger, step)

    pdb_features = {}
    with open(PDB_FEATURES_TSV, "r") as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader)

        for i, row in enumerate(reader, start=1):
            if not row:
                continue
            qseqid = row[0]
            pdb_features[qseqid] = {col: val for col, val in zip(header, row)}

            if i % 1_000_000 == 0:
                logger.info(f"{step}: read {i:,} rows (qseqids={len(pdb_features):,})")

    step_done(logger, step, t0, extra=f"{len(pdb_features):,} qseqids")
    return pdb_features, header


def build_sifts_mapping(logger):
    step = "build_sifts_mapping"
    t0 = step_start(logger, step)

    sifts_map = defaultdict(list)
    total_pairs = 0

    with open(SIFTS_MAPPING_TSV, "r") as f:
        reader = csv.reader(f, delimiter="\t")
        for i, row in enumerate(reader, start=1):
            if len(row) >= 2:
                qseqid, hit = row[0], row[1]
                sifts_map[qseqid].append(hit)
                total_pairs += 1

            if i % 5_000_000 == 0:
                logger.info(
                    f"{step}: {i:,} lines (qseqids={len(sifts_map):,}, pairs={total_pairs:,})"
                )

    step_done(logger, step, t0, extra=f"{len(sifts_map):,} qseqids, {total_pairs:,} pairs")
    return sifts_map


def compute_wanted_ids_from_sifts(logger, sifts_map):
    step = "compute_wanted_ids"
    t0 = step_start(logger, step)

    wanted_ids = set()
    total_hits = 0

    for i, hits in enumerate(sifts_map.values(), start=1):
        total_hits += len(hits)
        wanted_ids.update(hits)

        if i % 1_000_000 == 0:
            logger.info(
                f"{step}: scanned {i:,} qseqids (unique={len(wanted_ids):,}, total_hits={total_hits:,})"
            )

    step_done(logger, step, t0, extra=f"unique={len(wanted_ids):,}, total_hits={total_hits:,}")
    return wanted_ids


def build_uniref_mappings_filtered(logger, wanted_ids):
    """
    UniRef90 file format:
      col1: sseqid (representative)
      col2: member_count
      col3: comma-separated member_ids (may be empty)

    Build mappings ONLY for IDs in wanted_ids.

    Returns:
      uniref_direct: set of wanted IDs that appear as representatives (col1)
      uniref_cluster: dict mapping (wanted member_id) -> representative_id
    """
    step = "build_uniref_filtered"
    t0 = step_start(logger, f"{step} (wanted_ids={len(wanted_ids):,})")

    uniref_direct = set()
    uniref_cluster = {}  # member_id -> rep_id (single value)

    lines_seen = 0
    nonzero_mc = 0
    members_parsed = 0
    members_matched = 0

    with open(UNIREF90_CLUSTERS_TSV, "r") as f:
        header = next(f, None)
        if header:
            logger.info(f"{step}: header: {header.strip()}")

        for line_num, line in enumerate(f, start=2):
            lines_seen += 1
            line = line.rstrip("\n")
            if not line:
                continue

            # split into at most 3 fields
            parts = line.split("\t", 2)
            if len(parts) < 2:
                continue

            rep_id = parts[0]
            if not rep_id:
                continue

            if rep_id in wanted_ids:
                uniref_direct.add(rep_id)
                uniref_cluster.setdefault(rep_id, rep_id)

            mc_str = parts[1]
            if mc_str != "0":
                nonzero_mc += 1
                if len(parts) == 3 and parts[2]:
                    for member_id in parts[2].split(","):
                        if not member_id:
                            continue
                        members_parsed += 1
                        if member_id in wanted_ids:
                            members_matched += 1
                            uniref_cluster.setdefault(member_id, rep_id)

            if line_num % 2_000_000 == 0:
                logger.info(
                    f"{step}: {line_num:,} lines "
                    f"(direct={len(uniref_direct):,}, mapped={len(uniref_cluster):,}, "
                    f"nonzero_mc={nonzero_mc:,}, members_parsed={members_parsed:,}, members_matched={members_matched:,})"
                )

    step_done(
        logger,
        step,
        t0,
        extra=f"direct={len(uniref_direct):,}, mapped={len(uniref_cluster):,}, scanned_lines={lines_seen:,}",
    )
    return uniref_direct, uniref_cluster


def build_diamond_mapping(logger):
    step = "build_diamond_mapping"
    t0 = step_start(logger, step)

    diamond_map = defaultdict(set)
    diamond_counts = defaultdict(int)

    with open(AGGREGATED_HITS_TSV, "r") as f:
        for line_num, line in enumerate(f, start=1):
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                qseqid = parts[0]
                sseqids = parts[1].split(",") if parts[1] else []
                diamond_map[qseqid].update(sseqids)
                diamond_counts[qseqid] = len(sseqids)

            if line_num % 1_000_000 == 0:
                logger.info(f"{step}: {line_num:,} lines (qseqids={len(diamond_map):,})")

    step_done(logger, step, t0, extra=f"qseqids={len(diamond_map):,}")
    return diamond_map, diamond_counts


def load_added_sprot_ids(logger):
    step = "load_added_sprot_ids"
    t0 = step_start(logger, step)

    added_sprot = set()
    with open(ADDED_SPROT_IDS_TSV, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for i, row in enumerate(reader, start=1):
            sseqid = row.get("sseqid")
            if sseqid:
                added_sprot.add(sseqid)

            if i % 1_000_000 == 0:
                logger.info(f"{step}: read {i:,} rows (unique={len(added_sprot):,})")

    step_done(logger, step, t0, extra=f"added_sprot={len(added_sprot):,}")
    return added_sprot


def aggregate_filtered_hits(logger):
    """
    Combine_sifts_filteredHits.py logic:
    - sort by filter_rank
    - for each qseqid: keep top1 sseqid, top5 list, and a comma-string of top5
    """
    step = "aggregate_filtered_hits"
    t0 = step_start(logger, step)

    top1 = {}
    top5_list = {}
    top5_str = {}

    # streaming-friendly top5 collector: store up to 5 (rank, sseqid) per qseqid
    buckets = defaultdict(list)

    with open(FILTERED_HITS_TSV, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if "qseqid" not in reader.fieldnames or "sseqid" not in reader.fieldnames or "filter_rank" not in reader.fieldnames:
            raise RuntimeError(
                "filtered_hits.tsv must have columns: qseqid, sseqid, filter_rank"
            )

        for i, row in enumerate(reader, start=1):
            q = row.get("qseqid")
            s = row.get("sseqid")
            r = row.get("filter_rank")

            if not q or not s or r is None:
                continue

            try:
                rank = int(r)
            except Exception:
                continue

            # Keep only ranks that could possibly matter for top5
            if rank > 5:
                continue

            buckets[q].append((rank, s))

            if i % 2_000_000 == 0:
                logger.info(f"{step}: read {i:,} lines (qseqids_seen={len(buckets):,})")

    # finalize: sort within each qseqid and produce top1/top5
    for q, items in buckets.items():
        items.sort(key=lambda x: x[0])
        sseqids = [s for _, s in items][:5]
        if not sseqids:
            continue
        top1[q] = sseqids[0]
        top5_list[q] = sseqids
        top5_str[q] = ",".join(sseqids)

    step_done(logger, step, t0, extra=f"qseqids={len(top1):,}")
    return top1, top5_list, top5_str


# ====================== processing
def process_sifts_data(
    logger,
    pdb_features,
    feature_header,
    sifts_map,
    uniref_direct,
    uniref_cluster,
    diamond_map,
    diamond_counts,
    added_sprot,
    fh_top1,
    fh_top5_list,
    fh_top5_str,
):
    step = "process_sifts_data"
    t0 = step_start(logger, step)

    # Column sets (kept in original style from combine_sifts_filteredHits.py)
    sifts_columns = [
        "qseqid",
        "sifts_hits",
        "n_sifts_hits",
        "sifts_hits_inUr90",
        "n_sifts_hits_inUr90",
        "sifts_hits_inUr90Cluster",
        "n_sifts_hits_inUr90Cluster",
        "sifts_hits_in_addedSprot",
        "n_sifts_hits_in_addedSprot",
        "sifts_hitsORclustersORsprotIDs_inRaw",
        "n_sifts_hitsORclustersORsprotIDs_inRaw",
        "n_sseqid_in_DiamondRaw",
        "n_sifts_hits_inUr90OrClusterOrSprot",
        "stat_label",
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
    filtered_columns = [
        "filtered_hits",
        "siftsHits_in_top1filtered",
        "n_siftsHits_in_top1filtered",
        "siftsHits_in_top5filtered",
        "n_siftsHits_in_top5filtered",
    ]
    output_columns = sifts_columns + filtered_columns

    # stream output (avoid holding all rows in memory)
    SIFTS_SEQMAPPING_HITS_TSV.parent.mkdir(parents=True, exist_ok=True)
    with open(SIFTS_SEQMAPPING_HITS_TSV, "w") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(output_columns)

        for idx, qseqid in enumerate(pdb_features.keys(), start=1):
            sifts_hits = sifts_map.get(qseqid, [])
            features = pdb_features.get(qseqid, {})

            # ------------------ original searchOut logic ------------------
            if not sifts_hits:
                stat_label = (
                    "no_diamond_no_sifts"
                    if qseqid not in diamond_counts
                    else "yes_diamond_no_sifts"
                )

                sifts_hits_str = "NA"
                n_sifts_hits = "NA"

                ur90_str, n_in_ur90 = "NA", "NA"
                cluster_str, n_in_cluster = "NA", "NA"
                sprot_str, n_in_sprot = "NA", "NA"

                diamond_str, n_in_diamond = "NA", "NA"

                n_sseqid_in_diamondraw = diamond_counts.get(qseqid, "NA")
                n_sifts_hits_in_ur90_or_cluster_or_sprot = "NA"
            else:
                n_sifts_hits = len(sifts_hits)
                sifts_hits_str = ",".join(sifts_hits)

                in_ur90, in_cluster, in_sprot = [], [], []

                # ----- CLASSIFY EACH SIFTS HIT -----
                for hit in sifts_hits:
                    if hit in uniref_direct:
                        in_ur90.append(hit)
                        continue
                    elif hit in uniref_cluster:
                        rep_id = uniref_cluster[hit]
                        in_cluster.append(rep_id)
                        if hit in added_sprot:
                            in_sprot.append(hit)
                        continue
                    elif hit in added_sprot:
                        in_sprot.append(hit)

                # ----- STRINGIFY -----
                ur90_str = ",".join(in_ur90) if in_ur90 else "NA"
                n_in_ur90 = len(in_ur90) if in_ur90 else "NA"

                cluster_str = ",".join(in_cluster) if in_cluster else "NA"
                n_in_cluster = len(in_cluster) if in_cluster else "NA"

                sprot_str = ",".join(in_sprot) if in_sprot else "NA"
                n_in_sprot = len(in_sprot) if in_sprot else "NA"

                # ----- PRIORITY COMBINATION -----
                if in_sprot:
                    combined_hits = set(in_sprot)
                else:
                    combined_hits = set(in_ur90 + in_cluster)

                n_sifts_hits_in_ur90_or_cluster_or_sprot = len(combined_hits) if combined_hits else "NA"

                # ----- FILTER BY DIAMOND RESULTS -----
                diamond_hits = diamond_map.get(qseqid, set())
                filtered_hits_vs_diamond = [hit for hit in combined_hits if hit in diamond_hits]

                if filtered_hits_vs_diamond:
                    diamond_str = ",".join(sorted(filtered_hits_vs_diamond))
                    n_in_diamond = len(filtered_hits_vs_diamond)
                else:
                    diamond_str = "NA"
                    n_in_diamond = "NA"

                # ----- LABEL -----
                if qseqid not in diamond_counts:
                    stat_label = "no_diamond_yes_sifts"
                else:
                    if n_in_diamond == len(combined_hits) and n_in_diamond != "NA":
                        stat_label = "yes_diamond_yes_sifts_FullMatches"
                    elif n_in_diamond == "NA" or n_in_diamond == 0:
                        stat_label = "yes_diamond_yes_sifts_NoMatches"
                    else:
                        stat_label = "yes_diamond_yes_sifts_PartialMatches"

                n_sseqid_in_diamondraw = diamond_counts.get(qseqid, "NA")

            # ------------------ filteredHits overlap logic (integrated) ------------------
            # filtered_hits column = comma-joined top5 filtered hits for that qseqid
            filtered_hits_str = fh_top5_str.get(qseqid, "NA")
            top1_filtered_hit = fh_top1.get(qseqid)
            top5_filtered_hits_list = fh_top5_list.get(qseqid)

            # matching logic uses sifts_hitsORclustersORsprotIDs_inRaw
            sifts_hits_or_raw = str(diamond_str) if diamond_str is not None else "NA"
            if not sifts_hits_or_raw or sifts_hits_or_raw.lower() == "nan":
                sifts_hits_or_raw = "NA"

            if sifts_hits_or_raw == "NA" or top1_filtered_hit is None:
                siftsHits_in_top1filtered = "NA"
                n_siftsHits_in_top1filtered = 0
            else:
                sifts_set = set(sifts_hits_or_raw.split(","))
                if top1_filtered_hit in sifts_set:
                    siftsHits_in_top1filtered = top1_filtered_hit
                    n_siftsHits_in_top1filtered = 1
                else:
                    siftsHits_in_top1filtered = "NA"
                    n_siftsHits_in_top1filtered = 0

            if sifts_hits_or_raw == "NA" or not isinstance(top5_filtered_hits_list, list) or not top5_filtered_hits_list:
                siftsHits_in_top5filtered = "NA"
                n_siftsHits_in_top5filtered = 0
            else:
                sifts_set = set(sifts_hits_or_raw.split(","))
                matches = sorted(sifts_set.intersection(top5_filtered_hits_list))
                if matches:
                    siftsHits_in_top5filtered = ",".join(matches)
                    n_siftsHits_in_top5filtered = len(matches)
                else:
                    siftsHits_in_top5filtered = "NA"
                    n_siftsHits_in_top5filtered = 0

            # ------------------ output row (curated column subset) ------------------
            row_dict = {
                "qseqid": qseqid,
                "sifts_hits": sifts_hits_str,
                "n_sifts_hits": n_sifts_hits,
                "sifts_hits_inUr90": ur90_str,
                "n_sifts_hits_inUr90": n_in_ur90,
                "sifts_hits_inUr90Cluster": cluster_str,
                "n_sifts_hits_inUr90Cluster": n_in_cluster,
                "sifts_hits_in_addedSprot": sprot_str,
                "n_sifts_hits_in_addedSprot": n_in_sprot,
                "sifts_hitsORclustersORsprotIDs_inRaw": diamond_str,
                "n_sifts_hitsORclustersORsprotIDs_inRaw": n_in_diamond,
                "n_sseqid_in_DiamondRaw": n_sseqid_in_diamondraw,
                "n_sifts_hits_inUr90OrClusterOrSprot": n_sifts_hits_in_ur90_or_cluster_or_sprot,
                "stat_label": stat_label,
                # selected PDB feature columns (pulled from pdb_features dict)
                "qseq_length": features.get("qseq_length", "NA"),
                "qseq_taxid": features.get("qseq_taxid", "NA"),
                "qseq_tax_name": features.get("qseq_tax_name", "NA"),
                "qseq_domain_name": features.get("qseq_domain_name", "NA"),
                "qseq_n_source_organisms": features.get("qseq_n_source_organisms", "NA"),
                "qseq_parent_taxid": features.get("qseq_parent_taxid", "NA"),
                "qseq_grandparent_taxid": features.get("qseq_grandparent_taxid", "NA"),
                "qseq_source_type": features.get("qseq_source_type", "NA"),
                "qseq_description": features.get("qseq_description", "NA"),
                # filtered hit columns
                "filtered_hits": filtered_hits_str,
                "siftsHits_in_top1filtered": siftsHits_in_top1filtered,
                "n_siftsHits_in_top1filtered": n_siftsHits_in_top1filtered,
                "siftsHits_in_top5filtered": siftsHits_in_top5filtered,
                "n_siftsHits_in_top5filtered": n_siftsHits_in_top5filtered,
            }

            writer.writerow([row_dict.get(col, "NA") for col in output_columns])

            if idx % 10_000 == 0:
                logger.info(f"{step}: processed {idx:,} qseqids")

    step_done(logger, step, t0, extra=f"rows_written={len(pdb_features):,}")


# ====================== MAIN ======================
def main():
    logger = setup_logging()
    overall_t0 = time.time()

    try:
        pdb_features, feature_header = load_pdb_features(logger)

        sifts_map = build_sifts_mapping(logger)
        wanted_ids = compute_wanted_ids_from_sifts(logger, sifts_map)

        uniref_direct, uniref_cluster = build_uniref_mappings_filtered(logger, wanted_ids)

        diamond_map, diamond_counts = build_diamond_mapping(logger)
        added_sprot = load_added_sprot_ids(logger)

        fh_top1, fh_top5_list, fh_top5_str = aggregate_filtered_hits(logger)

        process_sifts_data(
            logger,
            pdb_features,
            feature_header,
            sifts_map,
            uniref_direct,
            uniref_cluster,
            diamond_map,
            diamond_counts,
            added_sprot,
            fh_top1,
            fh_top5_list,
            fh_top5_str,
        )

        elapsed = time.time() - overall_t0
        logger.info(f"pipeline finished successfully (total={elapsed:.2f}s)")

    except Exception as e:
        logger.error(f"Fatal error: {str(e)}", exc_info=True)
        raise


if __name__ == "__main__":
    main()
