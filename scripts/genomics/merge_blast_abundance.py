"""
merge_blast_abundance.py
─────────────────────────────────────────────────────────────────────────────
Merges BLAST results with transcriptome abundance data.

With --genome-type phased (default), maps phased gene IDs to non-phased
IDs via gene_id_map.csv before looking up abundances.

With --genome-type non_phased, the BLAST hit IDs are already non-phased
and are used directly for abundance lookup.

Output: one CSV per organism under paper/tables/ with columns:
  enzyme, hit_id, non_phased_id, expect, bit_score, coverage,
  identity, similarity, gaps, td, sp, chr, hit_description,
  transcript_coverage, FPKM, TPM

Usage:
    python scripts/genomics/merge_blast_abundance.py [--overwrite] [-g phased|non_phased]
─────────────────────────────────────────────────────────────────────────────
"""

from __future__ import annotations

import csv
import logging
import os
import sys

from data import (
    GENOME_PATH,
    PAPER_PATH,
    PROTEINS_PATH,
    TRANSCRIPTOME_PATH,
    organisms,
    proteins,
)
from logger import build_logger

log = build_logger("merge_blast_abundance")

OVERWRITE = "--overwrite" in sys.argv

GENOME_TYPE = "phased"
for i, arg in enumerate(sys.argv):
    if arg in ("--genome-type", "-g") and i + 1 < len(sys.argv):
        GENOME_TYPE = sys.argv[i + 1]


TRANSCRIPTOME_FILENAMES: dict[str, str] = {
    "tucunaca": "tucunaca.tab",
    "caupuri": "caupuri.tab",
}


def get_protein_organisms(proteins):
    result = []
    for prot_name, info in proteins.items():
        for species, org_info in info["organisms"].items():
            uniprot = org_info["uniprot_id"]
            fasta_name = f"{prot_name}_{uniprot}"
            result.append((prot_name, species, fasta_name))
    return result


def load_mapping(map_csv: str) -> dict[str, str]:
    """Load phased → non-phased gene ID mapping from CSV.

    Returns dict {phased_id: non_phased_id}.
    """
    mapping: dict[str, str] = {}
    if not os.path.exists(map_csv):
        return mapping
    with open(map_csv) as f:
        reader = csv.DictReader(f)
        for row in reader:
            mapping[row["phased_id"]] = row["non_phased_id"]
    log.info("Loaded %d mapping entries from %s", len(mapping), map_csv)
    return mapping


def load_abundances(
    organism: str,
) -> dict[str, dict[str, float]]:
    """Load transcriptome abundances from the StringTie2 .tab file.

    Returns dict {gene_id: {coverage, fpkm, tpm}}.
    """
    filename = TRANSCRIPTOME_FILENAMES.get(organism)
    if not filename:
        log.warning("No transcriptome filename configured for %s", organism)
        return {}

    ab_path = f"{TRANSCRIPTOME_PATH}/{organism}/{filename}"
    if not os.path.exists(ab_path):
        log.warning("Transcriptome file not found: %s", ab_path)
        return {}

    abundances: dict[str, dict[str, float]] = {}
    with open(ab_path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            gene_id = row["Gene ID"].strip()
            try:
                abundances[gene_id] = {
                    "coverage": float(row.get("Coverage", 0)),
                    "fpkm": float(row.get("FPKM", 0)),
                    "tpm": float(row.get("TPM", 0)),
                }
            except (ValueError, KeyError) as e:
                log.debug("Skipping row for %s: %s", gene_id, e)

    log.info("Loaded %d abundance entries for %s", len(abundances), organism)
    return abundances


def process_organism(organism: str) -> None:
    """Process one organism: merge all BLAST CSVs with abundances."""
    prots = get_protein_organisms(proteins)

    # Load the abundance data once per organism
    abundances = load_abundances(organism)
    if not abundances:
        log.warning("No abundance data for %s, skipping.", organism)
        return

    merged_rows: list[dict] = []

    for prot_name, ref_species, fasta_name in prots:
        blast_csv = (
            f"{PROTEINS_PATH}/{prot_name}/analysis/genomic/blast"
            f"/{organism}/{ref_species}/{fasta_name}_blast_filtered.csv"
        )

        if not os.path.exists(blast_csv):
            log.debug("BLAST CSV not found: %s", blast_csv)
            continue

        if GENOME_TYPE == "phased":
            map_csv = (
                f"{PROTEINS_PATH}/{prot_name}/analysis/genomic/blast"
                f"/{organism}/{ref_species}/gene_id_map.csv"
            )
            mapping = load_mapping(map_csv)
        else:
            mapping = {}

        # Read BLAST CSV
        with open(blast_csv) as f:
            reader = csv.DictReader(f, delimiter=";")
            for row in reader:
                hit_id = row["hit_id"]

                if GENOME_TYPE == "phased":
                    non_phased_id = mapping.get(hit_id, "")
                    lookup_id = non_phased_id.rsplit(".", 1)[0] if non_phased_id else ""
                else:
                    non_phased_id = hit_id
                    lookup_id = hit_id.rsplit(".", 1)[0]

                ab = abundances.get(lookup_id, {}) if lookup_id else {}

                merged_rows.append({
                    "enzyme": row.get("enzyme", ""),
                    "hit_id": hit_id,
                    "non_phased_id": non_phased_id,
                    "expect": row.get("expect", ""),
                    "hit_len": row.get("hit_len", ""),
                    "bit_score": row.get("bit_score", ""),
                    "coverage": row.get("coverage", ""),
                    "identity": row.get("identity", ""),
                    "similarity": row.get("similarity", ""),
                    "gaps": row.get("gaps", ""),
                    "td": row.get("td", ""),
                    "sp": row.get("sp", ""),
                    "chr": row.get("chr", ""),
                    "hit_description": row.get("hit_description", ""),
                    "transcript_coverage": ab.get("coverage", ""),
                    "FPKM": ab.get("fpkm", ""),
                    "TPM": ab.get("tpm", ""),
                })

    if not merged_rows:
        log.warning("No merged rows for %s", organism)
        return

    # Write merged CSV to paper/tables/
    tables_dir = f"{PAPER_PATH}/tables"
    os.makedirs(tables_dir, exist_ok=True)
    out_path = f"{tables_dir}/{organism}_blast_abundance.csv"

    if os.path.exists(out_path) and not OVERWRITE:
        log.info("Output already exists, skipping: %s", out_path)
        log.info("  Use --overwrite to force regeneration.")
        return

    fieldnames = [
        "enzyme", "hit_id", "non_phased_id", "expect", "hit_len",
        "bit_score", "coverage", "identity", "similarity", "gaps",
        "td", "sp", "chr", "hit_description",
        "transcript_coverage", "FPKM", "TPM",
    ]

    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(merged_rows)

    log.info(
        "Wrote %d rows → %s", len(merged_rows), out_path,
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

    log.info("Done merging BLAST results with abundances.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
