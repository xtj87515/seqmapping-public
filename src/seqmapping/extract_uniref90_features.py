#!/usr/bin/env python3
"""
Full pipeline:
 - Parse UniRef90 XML -> data/intermediate/uniref90_features.tsv (streaming)
 - Optionally write benchmarks/data/uniref90_clusters.tsv
 - Read data/intermediate/added_sprot_ids.tsv and append those UniProt accessions' features
   by parsing data/raw/uniprot/uniprot_sprot.xml (streaming).

NOTE:
 - The UniRef90 clusters TSV will contain NO UniParc(UPI) IDs anywhere in the member_ids list.
 - member_count written to clusters TSV is the filtered member count (after removing UPI IDs), so it matches the member_ids list.
"""

from __future__ import annotations

import csv
import re
import time
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Dict, Set, Tuple

from seqmapping.utils.logging import get_logger, start_run
from seqmapping.utils.paths import (
    ADDED_SPROT_IDS_TSV,
    BENCHMARKS_DIR,
    NCBI_PARENT_TAXIDS_TSV,
    UNIREF90_FEATURES_TSV,
    UNIREF90_XML,
    UNIPROT_REF_PROTEOMES,
    UNIPROT_SPROT_FASTA,
    UNIPROT_SPROT_XML,
    UNIREF90_CLUSTERS_TSV,
    ensure_directories,
)
from seqmapping.utils.config import WRITE_UNIREF90_CLUSTERS


# -----------------------------------------------------------------------------
# Utility loaders
# -----------------------------------------------------------------------------
def load_parent_maps(path: Path, logger) -> Tuple[Dict[str, str], Dict[str, str]]:
    # Returns (parent_map, grandparent_map) where keys are taxid strings.
    parent: Dict[str, str] = {}
    grandparent: Dict[str, str] = {}

    if not path.exists():
        logger.warning("Parent taxids file not found: %s", path)
        return parent, grandparent

    with path.open("r", encoding="utf-8", errors="replace") as fh:
        _header = fh.readline()  # consume header row
        for line in fh:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            taxid = parts[0].strip()
            parent_tax = parts[1].strip()
            parent[taxid] = parent_tax
            if len(parts) >= 3:
                grandparent[taxid] = parts[2].strip()

    logger.info(
        "Loaded parent map entries: %d, grandparent map entries: %d",
        len(parent),
        len(grandparent),
    )
    return parent, grandparent


def load_ref_taxids(path: Path, logger) -> Set[str]:
    # Produce a set of taxids that are reference proteomes.
    s: Set[str] = set()
    if not path.exists():
        logger.warning("Ref proteomes file not found: %s", path)
        return s

    with path.open("r", encoding="utf-8", errors="replace") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        header_lower = [h.lower() for h in header]

        idx = None
        for name in ("organism id", "organism_id", "organismid", "organism-id"):
            if name in header_lower:
                idx = header_lower.index(name)
                break

        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if not parts or (len(parts) == 1 and not parts[0].strip()):
                continue

            if idx is not None and idx < len(parts):
                v = parts[idx].strip()
                if v:
                    s.add(v)
            else:
                # heuristic: first numeric-looking cell
                for p in parts:
                    p = p.strip()
                    if p.isdigit():
                        s.add(p)
                        break

    logger.info("Loaded %d reference-proteome taxids (heuristic)", len(s))
    return s


def extract_swissprot_ids_from_fasta(fasta_path: Path, logger) -> Set[str]:
    # Grab UniProt accessions from FASTA headers like >sp|P12345|...
    ids: Set[str] = set()
    if not fasta_path.exists():
        logger.warning("Swiss-Prot FASTA not found: %s", fasta_path)
        return ids

    pat = re.compile(r"^>sp\|([A-Z0-9]+)\|")
    with fasta_path.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                m = pat.match(line)
                if m:
                    ids.add(m.group(1))

    logger.info("Extracted %d SwissProt accessions from FASTA", len(ids))
    return ids


# -----------------------------------------------------------------------------
# Helpers: identify and filter UniParc/UPI IDs in UniRef member lists
# -----------------------------------------------------------------------------
def _is_uniparc_upi(value: str | None) -> bool:
    # UniParc accessions are UPI########

    if not value:
        return False
    v = value.strip()
    return v.startswith("UPI")


def _extract_member_identifier(db_ref, ns: str | None) -> str:
    # Prefer UniProtKB accession property if present; else fall back to dbReference @id.

    acc_prop = (
        db_ref.find(f"{ns}property[@type='UniProtKB accession']") if ns else db_ref.find("property[@type='UniProtKB accession']")
    )
    if acc_prop is not None:
        return (acc_prop.attrib.get("value", "") or "").strip()

    return (db_ref.attrib.get("id", "") or "").strip()


def _collect_filtered_member_ids(elem, ns: str | None) -> Tuple[list[str], int]:
    # Return (filtered_member_ids, filtered_upi_count). Ensures NO IDs starting with 'UPI'.

    filtered: list[str] = []
    filtered_upi = 0

    members = elem.findall(f"{ns}member") if ns else elem.findall("member")
    for m in members:
        db_ref = m.find(f"{ns}dbReference") if ns else m.find("dbReference")
        if db_ref is None:
            continue

        mid = _extract_member_identifier(db_ref, ns)
        if _is_uniparc_upi(mid):
            filtered_upi += 1
            continue

        # Sometimes UniRef member dbReference could contain multiple identifiers or empty; keep only non-empty non-UPI values.
        if mid:
            filtered.append(mid)

    return filtered, filtered_upi


# -----------------------------------------------------------------------------
# Core: parse UniRef90 and write features/clusters
# -----------------------------------------------------------------------------
def parse_uniref90_and_write(
    xml_file: Path,
    features_out: Path,
    clusters_out: Path,
    swissprot_ids: Set[str],
    ref_taxids: Set[str],
    parent_map: Dict[str, str],
    grandparent_map: Dict[str, str],
    write_clusters: bool = True,
    logger=None,
) -> None:
    """
    Stream through UniRef90 XML, extract representative member info and (optionally) clusters,
    write features TSV and clusters TSV incrementally (memory efficient).

    - Skips UniParc-based UniRef entries whose entry id is 'UniRef90_UPI...'.
    - If write_clusters is enabled, skips member_ids with ALL UPI IDs removed.
    """

    if not xml_file.exists():
        logger.error("UniRef90 XML not found: %s", xml_file)
        raise FileNotFoundError(str(xml_file))

    features_out.parent.mkdir(parents=True, exist_ok=True)
    with features_out.open("w", newline="", encoding="utf-8") as fh_features:
        writer_features = csv.writer(fh_features, delimiter="\t")
        writer_features.writerow(
            [
                "sseqid",
                "sseq_taxid",
                "sseq_tax_name",
                "sseq_length",
                "sseq_is_uniprot_sprot",
                "sseq_is_ref_proteome",
                "sseq_parent_taxid",
                "sseq_grandparent_taxid",
            ]
        )

        fh_clusters = None
        writer_clusters = None
        if write_clusters:
            clusters_out.parent.mkdir(parents=True, exist_ok=True)
            fh_clusters = clusters_out.open("w", newline="", encoding="utf-8")
            writer_clusters = csv.writer(fh_clusters, delimiter="\t")
            writer_clusters.writerow(["sseqid", "member_count", "member_ids"])

        ns = None
        total = 0
        kept = 0
        upi_entries_skipped = 0
        members_upi_filtered = 0

        start = time.time()

        context = ET.iterparse(str(xml_file), events=("end",))
        for _event, elem in context:
            tag = elem.tag

            if ns is None and "}" in tag:
                ns_uri = tag.split("}")[0].strip("{")
                ns = "{" + ns_uri + "}"

            # process only entry elements
            if (ns and tag == f"{ns}entry") or (ns is None and tag.endswith("entry")):
                total += 1

                entry_id = elem.attrib.get("id", "")
                if not entry_id:
                    elem.clear()
                    continue

                # Filter out UniParc-based UniRef entries: UniRef90_UPI...
                if entry_id.startswith("UniRef90_UPI"):
                    upi_entries_skipped += 1
                    elem.clear()
                    continue

                sseqid = entry_id.replace("UniRef90_", "", 1)

                taxid = ""
                organism = ""
                sseq_length = ""

                rep = elem.find(f"{ns}representativeMember") if ns else elem.find("representativeMember")
                if rep is not None:
                    dbref = rep.find(f"{ns}dbReference") if ns else rep.find("dbReference")
                    if dbref is not None:
                        props = dbref.findall(f"{ns}property") if ns else dbref.findall("property")
                        for prop in props:
                            t = prop.attrib.get("type", "")
                            if t == "NCBI taxonomy":
                                taxid = prop.attrib.get("value", "")
                            elif t == "source organism":
                                organism = prop.attrib.get("value", "")
                            elif t == "length":
                                sseq_length = prop.attrib.get("value", "")

                is_sprot = 1 if sseqid in swissprot_ids else 0
                is_ref = 1 if (taxid and taxid in ref_taxids) else 0
                parent_tax = parent_map.get(taxid, "") if taxid else ""
                grandparent_tax = grandparent_map.get(taxid, "") if taxid else ""

                writer_features.writerow(
                    [sseqid, taxid, organism, sseq_length, is_sprot, is_ref, parent_tax, grandparent_tax]
                )
                kept += 1

                if write_clusters and writer_clusters is not None:
                    member_ids, upi_filtered = _collect_filtered_member_ids(elem, ns)
                    members_upi_filtered += upi_filtered

                    # NOTE: member_count matches filtered member_ids list (no UPI anywhere).
                    writer_clusters.writerow([sseqid, str(len(member_ids)), ",".join(member_ids)])

                elem.clear()

                if total % 1_000_000 == 0:
                    elapsed = time.time() - start
                    logger.info(
                        "Parsed %d entries (kept %d, skipped UPI entries %d, filtered UPI members %d) in %.1fs",
                        total,
                        kept,
                        upi_entries_skipped,
                        members_upi_filtered,
                        elapsed,
                    )

        if fh_clusters is not None:
            fh_clusters.close()

        elapsed = time.time() - start
        logger.info(
            "Finished parsing UniRef90: total=%d, kept=%d, skipped_UPI_entries=%d, filtered_UPI_members=%d, time=%.1fs",
            total,
            kept,
            upi_entries_skipped,
            members_upi_filtered,
            elapsed,
        )


# -----------------------------------------------------------------------------
# Core: append Swiss-Prot features for added IDs
# -----------------------------------------------------------------------------
def read_added_ids(added_file: Path, logger) -> Dict[str, int | None]:
    # Read added_sprot_ids.tsv expected columns "sseqid\t sseq_length"
    # Returns dict accession -> int(length) or None if missing.
    d: Dict[str, int | None] = {}
    if not added_file.exists():
        logger.warning("Added IDs file not found: %s", added_file)
        return d

    with added_file.open("r", encoding="utf-8", errors="replace") as fh:
        _header = fh.readline()  # skip header
        for line in fh:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) >= 2:
                acc = parts[0].strip()
                try:
                    length = int(parts[1])
                except Exception:
                    length = None
                d[acc] = length
            elif len(parts) == 1:
                d[parts[0].strip()] = None

    logger.info("Loaded %d accessions from %s", len(d), added_file)
    return d


def append_swissprot_features_from_added(
    added_dict: Dict[str, int | None],
    sprot_xml: Path,
    features_out: Path,
    ref_taxids: Set[str],
    parent_map: Dict[str, str],
    grandparent_map: Dict[str, str],
    logger=None,
) -> int:
    # Stream parse Swiss-Prot XML and append features for accessions in added_dict.
    # Write appended rows into features_out (append mode).
    if not added_dict:
        logger.info("No added Swiss-Prot accessions to append.")
        return 0

    if not sprot_xml.exists():
        logger.error("Swiss-Prot XML not found: %s", sprot_xml)
        return 0

    # Ensure features file exists and has header
    if not features_out.exists():
        features_out.parent.mkdir(parents=True, exist_ok=True)
        with features_out.open("w", newline="", encoding="utf-8") as fh:
            writer = csv.writer(fh, delimiter="\t")
            writer.writerow(
                [
                    "sseqid",
                    "sseq_taxid",
                    "sseq_tax_name",
                    "sseq_length",
                    "sseq_is_uniprot_sprot",
                    "sseq_is_ref_proteome",
                    "sseq_parent_taxid",
                    "sseq_grandparent_taxid",
                ]
            )

    found: Set[str] = set()
    appended = 0
    total = 0
    start = time.time()

    with features_out.open("a", newline="", encoding="utf-8") as fh_out:
        writer = csv.writer(fh_out, delimiter="\t")

        ns = None
        context = ET.iterparse(str(sprot_xml), events=("end",))
        for _event, elem in context:
            tag = elem.tag

            if ns is None and "}" in tag:
                ns_uri = tag.split("}")[0].strip("{")
                ns = "{" + ns_uri + "}"

            if (ns and tag == f"{ns}entry") or (ns is None and tag.endswith("entry")):
                total += 1

                acc_elems = elem.findall(f"{ns}accession") if ns else elem.findall("accession")
                accs = [ae.text.strip() for ae in acc_elems if ae is not None and ae.text]
                matched = [a for a in accs if a in added_dict]
                if not matched:
                    elem.clear()
                    continue

                taxid = ""
                org_name = ""

                org = elem.find(f"{ns}organism") if ns else elem.find("organism")
                if org is not None:
                    dbrefs = org.findall(f"{ns}dbReference") if ns else org.findall("dbReference")
                    for dbref in dbrefs:
                        if dbref.attrib.get("type", "").lower() in ("ncbi taxonomy", "ncbi_taxonomy"):
                            taxid = dbref.attrib.get("id", "")
                            break

                    names = org.findall(f"{ns}name") if ns else org.findall("name")
                    for name_el in names:
                        if name_el.attrib.get("type", "") == "scientific":
                            org_name = name_el.text or ""
                            break
                    if not org_name and names:
                        org_name = names[0].text or ""

                seq_elem = elem.find(f"{ns}sequence") if ns else elem.find("sequence")
                seq_len_attr = None
                if seq_elem is not None:
                    seq_len_attr = seq_elem.attrib.get("length")

                for acc in matched:
                    added_len = added_dict.get(acc)
                    if added_len is not None:
                        seq_len = added_len
                    else:
                        seq_len = int(seq_len_attr) if (seq_len_attr and seq_len_attr.isdigit()) else ""

                    is_sprot = 1
                    is_ref = 1 if (taxid and taxid in ref_taxids) else 0
                    parent_tax = parent_map.get(taxid, "") if taxid else ""
                    grandparent_tax = grandparent_map.get(taxid, "") if taxid else ""

                    writer.writerow([acc, taxid, org_name, seq_len, is_sprot, is_ref, parent_tax, grandparent_tax])
                    appended += 1
                    found.add(acc)

                elem.clear()

                if total % 500_000 == 0:
                    elapsed = time.time() - start
                    logger.info("Scanned %d SwissProt entries, appended %d, elapsed %.1fs", total, appended, elapsed)

    # Add placeholder rows for any IDs not found in Swiss-Prot XML
    not_found = set(added_dict.keys()) - found
    if not_found:
        with features_out.open("a", newline="", encoding="utf-8") as fh_out:
            writer2 = csv.writer(fh_out, delimiter="\t")
            for acc in not_found:
                seq_len = added_dict.get(acc, "")
                writer2.writerow([acc, "", "", seq_len, 1, 0, "", ""])
                appended += 1
        logger.info(
            "%d accessions from %s not found in Swiss-Prot XML; added placeholder rows.",
            len(not_found),
            ADDED_SPROT_IDS_TSV,
        )

    logger.info(
        "Appended total %d Swiss-Prot feature rows (found %d, placeholders %d)",
        appended,
        len(found),
        len(not_found),
    )
    return appended


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------
def main() -> int:
    ensure_directories()

    module_stem = Path(__file__).stem
    logger = get_logger(__name__, module_stem)
    start_run(logger, run_name=module_stem, argv=None)

    t0 = time.time()
    logger.info("Pipeline started")

    # Inputs (from seqmapping.utils.paths / config)
    uniref90_xml = UNIREF90_XML
    sprot_xml = UNIPROT_SPROT_XML
    added_ids = ADDED_SPROT_IDS_TSV
    parent_taxids = NCBI_PARENT_TAXIDS_TSV
    ref_proteomes = UNIPROT_REF_PROTEOMES
    sprot_fasta = UNIPROT_SPROT_FASTA

    # Outputs
    out_features = UNIREF90_FEATURES_TSV
    out_clusters = UNIREF90_CLUSTERS_TSV

    logger.info("UniRef90 XML: %s", uniref90_xml)
    logger.info("Swiss-Prot XML: %s", sprot_xml)
    logger.info("Added IDs TSV: %s", added_ids)
    logger.info("Parent taxids TSV: %s", parent_taxids)
    logger.info("Reference proteomes TSV: %s", ref_proteomes)
    logger.info("Swiss-Prot FASTA: %s", sprot_fasta)
    logger.info("Output features TSV: %s", out_features)
    if WRITE_UNIREF90_CLUSTERS:
        logger.info("Output clusters TSV: %s", out_clusters)

    parent_map, grandparent_map = load_parent_maps(parent_taxids, logger)
    ref_taxids = load_ref_taxids(ref_proteomes, logger)
    swissprot_ids = extract_swissprot_ids_from_fasta(sprot_fasta, logger)

    # UniRef90 parse (always rebuild, or gate this if you want)
    if out_features.exists():
        logger.info("%s exists, rebuilding (config-driven run).", out_features)

    logger.info(
        "Parsing UniRef90 XML and writing features%s.",
        " (and clusters)" if WRITE_UNIREF90_CLUSTERS else "",
    )

    parse_uniref90_and_write(
        xml_file=uniref90_xml,
        features_out=out_features,
        clusters_out=out_clusters,
        swissprot_ids=swissprot_ids,
        ref_taxids=ref_taxids,
        parent_map=parent_map,
        grandparent_map=grandparent_map,
        write_clusters=WRITE_UNIREF90_CLUSTERS,
        logger=logger,
    )

    # Append Swiss-Prot additions
    added = read_added_ids(added_ids, logger)
    if added:
        appended_count = append_swissprot_features_from_added(
            added_dict=added,
            sprot_xml=sprot_xml,
            features_out=out_features,
            ref_taxids=ref_taxids,
            parent_map=parent_map,
            grandparent_map=grandparent_map,
            logger=logger,
        )
        logger.info("Appended %d Swiss-Prot feature rows", appended_count)
    else:
        logger.info("No added Swiss-Prot IDs to append")

    elapsed = time.time() - t0
    logger.info("Pipeline complete in %.1f seconds", elapsed)
    return 0



if __name__ == "__main__":
    raise SystemExit(main())
