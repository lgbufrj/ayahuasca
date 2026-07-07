"""
map_gene_ids.py
─────────────────────────────────────────────────────────────────────────────
Maps phased genome gene IDs → non-phased genome gene IDs using BLASTN
of CDS sequences. Since BLAST is run against the phased genome but
transcriptome abundances refer to the non-phased genome, this bridge
is needed to merge the two data sources.

Usage:
    python scripts/genomics/map_gene_ids.py [--overwrite] [--workers N]

The script iterates over all non-reference organisms (tucunaca, caupuri)
and all reference-species × protein combinations, producing one mapping
CSV per organism:
    proteins/<prot>/analysis/genomic/blast/<ooi>/<ref_species>/gene_id_map.csv
─────────────────────────────────────────────────────────────────────────────
"""

from __future__ import annotations

import csv
import logging
import os
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from Bio import SeqIO

from data import GENOME_PATH, PROTEINS_PATH, organisms, proteins
from logger import build_logger

log = build_logger("map_gene_ids")

EXPECT_CUTOFF = 1e-10
IDENTITY_CUTOFF = 95.0
OVERWRITE = "--overwrite" in sys.argv


def get_protein_organisms(proteins):
    result = []
    for prot_name, info in proteins.items():
        for species, org_info in info["organisms"].items():
            uniprot = org_info["uniprot_id"]
            fasta_name = f"{prot_name}_{uniprot}"
            result.append((prot_name, species, fasta_name))
    return result


def extract_hit_ids(blast_csv: str) -> list[str]:
    """Read the BLAST CSV and return unique hit_ids (phased gene IDs)."""
    ids: list[str] = []
    if not os.path.exists(blast_csv):
        return ids
    with open(blast_csv) as f:
        reader = csv.DictReader(f, delimiter=";")
        for row in reader:
            ids.append(row["hit_id"])
    # Remove duplicates while preserving order
    seen: set[str] = set()
    unique: list[str] = []
    for i in ids:
        if i not in seen:
            seen.add(i)
            unique.append(i)
    return unique


def extract_cds_sequences(
    cds_fasta: str, hit_ids: list[str]
) -> list[tuple[str, str]]:
    """Extract CDS sequences for the given hit_ids from the FASTA file.

    Returns list of (id, sequence_string) tuples.
    """
    seq_by_id: dict[str, str] = {}
    for rec in SeqIO.parse(cds_fasta, "fasta"):
        rec_id = rec.id.split("|")[0]
        seq_by_id[rec_id] = str(rec.seq)

    result: list[tuple[str, str]] = []
    for hid in hit_ids:
        seq = seq_by_id.get(hid)
        if seq:
            result.append((hid, seq))
        else:
            log.warning("CDS not found for phased gene %s", hid)
    return result


def run_blastn_and_parse(
    query_fasta: str, db_path: str, max_targets: int = 1
) -> list[dict]:
    """Run blastn with the given query FASTA against the database.

    Returns list of dicts with keys: qseqid, sseqid, pident, length,
    mismatch, evalue, bitscore.
    """
    cmd = [
        "blastn",
        "-query", query_fasta,
        "-db", db_path,
        "-outfmt",
        "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore",
        "-max_target_seqs", str(max_targets),
        "-max_hsps", "1",
        "-evalue", str(EXPECT_CUTOFF),
    ]
    try:
        result = subprocess.run(
            cmd, capture_output=True, text=True, check=True
        )
    except subprocess.CalledProcessError as e:
        log.error("blastn failed: %s", e.stderr.strip())
        return []

    rows: list[dict] = []
    for line in result.stdout.strip().splitlines():
        if not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) < 7:
            continue
        rows.append({
            "qseqid": parts[0],
            "sseqid": parts[1],
            "pident": float(parts[2]),
            "length": int(parts[3]),
            "mismatch": int(parts[4]),
            "evalue": float(parts[9]),
            "bitscore": float(parts[10]),
        })
    return rows


def build_mapping(
    organism: str, hit_ids: list[str]
) -> dict[str, str]:
    """Map phased → non-phased gene IDs by BLASTN of CDS sequences.

    Returns dict {phased_id: non_phased_id}.
    """
    phased_cds = f"{GENOME_PATH}/{organism}/phased/cds.fasta"
    non_phased_db = f"{GENOME_PATH}/{organism}/non_phased/blast/cds_non_phased_db"

    if not os.path.exists(phased_cds):
        log.warning("Phased CDS not found for %s: %s", organism, phased_cds)
        return {}

    if not any(
        os.path.exists(f"{non_phased_db}.{ext}")
        for ext in ["nsq", "nin", "nhr"]
    ):
        log.warning(
            "Non-phased BLAST DB not found for %s: %s", organism, non_phased_db
        )
        return {}

    sequences = extract_cds_sequences(phased_cds, hit_ids)
    if not sequences:
        log.warning("No CDS sequences found for the given hit_ids in %s", organism)
        return {}

    # Write query sequences to a temp FASTA
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".fasta", delete=False
    ) as tmp:
        tmp_path = tmp.name
        for sid, seq in sequences:
            tmp.write(f">{sid}\n{seq}\n")

    try:
        blast_results = run_blastn_and_parse(tmp_path, non_phased_db)
    finally:
        os.unlink(tmp_path)

    mapping: dict[str, str] = {}
    for row in blast_results:
        phased = row["qseqid"]
        non_phased = row["sseqid"].split("|")[0] if "|" in row["sseqid"] else row["sseqid"]
        # Only accept high-identity matches
        if row["pident"] >= IDENTITY_CUTOFF:
            # Keep the best hit (first due to -max_target_seqs 1)
            if phased not in mapping:
                mapping[phased] = non_phased

    return mapping


def process_organism(organism: str) -> None:
    """Process one organism: find all BLAST CSVs, build mappings."""
    prots = get_protein_organisms(proteins)

    for prot_name, ref_species, fasta_name in prots:
        blast_csv = (
            f"{PROTEINS_PATH}/{prot_name}/analysis/genomic/blast"
            f"/{organism}/{ref_species}/{fasta_name}_blast_filtered.csv"
        )

        if not os.path.exists(blast_csv):
            continue

        hit_ids = extract_hit_ids(blast_csv)
        if not hit_ids:
            log.info(
                "No hit IDs in %s for %s / %s", blast_csv, prot_name, organism
            )
            continue

        log.info(
            "Mapping %d phased hits for %s / %s / %s",
            len(hit_ids), organism, ref_species, prot_name,
        )

        mapping = build_mapping(organism, hit_ids)

        if not mapping:
            log.warning("No mapping found for %s / %s / %s", organism, ref_species, prot_name)
            continue

        # Write mapping CSV
        out_dir = (
            f"{PROTEINS_PATH}/{prot_name}/analysis/genomic/blast"
            f"/{organism}/{ref_species}"
        )
        os.makedirs(out_dir, exist_ok=True)
        out_path = f"{out_dir}/gene_id_map.csv"

        if os.path.exists(out_path) and not OVERWRITE:
            log.info("Mapping already exists, skipping: %s", out_path)
            continue

        with open(out_path, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["phased_id", "non_phased_id"])
            for phased_id, non_phased_id in mapping.items():
                writer.writerow([phased_id, non_phased_id])

        mapped_count = len(mapping)
        log.info(
            "Mapped %d / %d IDs → %s",
            mapped_count, len(hit_ids), out_path,
        )


def main() -> int:
    oois = [
        org
        for org, meta in organisms.items()
        if not meta.get("ref", False)
    ]

    if not oois:
        log.warning("No non-reference organisms found.")
        return 0

    for organism in oois:
        log.info("Processing organism: %s", organism)
        process_organism(organism)

    log.info("Done mapping gene IDs.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
