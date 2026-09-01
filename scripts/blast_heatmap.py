import os
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import seaborn as sns
from Bio.Blast import NCBIXML

import data
from logger import build_logger

warnings.filterwarnings('ignore')

log = build_logger("blast_heatmap")

BASE = Path(data.BASE_PATH)
PROJECT = Path(data.PROJECT_PATH)
PROTEINS = Path(data.PROTEINS_PATH)
PAPER = Path(data.PAPER_PATH)
FIGURES = PAPER / "figures"
FIGURES.mkdir(parents=True, exist_ok=True)

sns.set_style("whitegrid")
sns.set_context("paper", font_scale=1.2)

TARGET_ORGS = ["tucunaca", "caupuri"]
ORG_LABELS = {meta["species"].split()[0]: org for org, meta in data.organisms.items()}
ORG_DISPLAY = {
    "tucunaca": "B. caapi (Tucunaca)",
    "caupuri": "B. caapi (Caupuri)",
}

REF_SPECIES_LABEL = {
    "arabidopsis": "A. thaliana",
    "rice": "O. sativa",
    "catharanthus": "C. roseus",
    "tabaco": "N. tabacum",
}


def collect_blast_data():
    rows = []
    for prot_key, prot_data in data.proteins.items():
        ptn_abbrev = prot_data["abbreviation"]
        for ref_org_key, ref_org_data in prot_data["organisms"].items():
            uniprot_id = ref_org_data["uniprot_id"]
            fasta_name = f"{prot_key}_{uniprot_id}"
            for organism in TARGET_ORGS:
                xml_path = (
                    PROTEINS / prot_key / "analysis" / "genomic" / "blast"
                    / organism / ref_org_key
                    / f"{fasta_name}.xml"
                )
                if not xml_path.exists():
                    continue
                try:
                    hits = []
                    with open(xml_path) as f:
                        for blast_record in NCBIXML.parse(f):
                            query_len = blast_record.query_length
                            for alignment in blast_record.alignments:
                                for hsp in alignment.hsps:
                                    hits.append({
                                        "bit_score": hsp.bits,
                                    })

                    best = max(hits, key=lambda h: h["bit_score"]) if hits else None
                    rows.append({
                        "enzyme": ptn_abbrev,
                        "ref_org": ref_org_key,
                        "ref_label": REF_SPECIES_LABEL.get(ref_org_key, ref_org_key),
                        "target_organism": organism,
                        "best_bit_score": best["bit_score"] if best else np.nan,
                    })
                except Exception as e:
                    log.warning("Error parsing %s: %s", xml_path, e)
    return pd.DataFrame(rows)


ENZYME_LABELS = {
    "TDC": "TDC",
    "STS": "STR",
    "SARD4": "SARD4",
    "T5H": "T5H",
    "ASMT": "ASMT",
    "POX": "Peroxidase",
}


def main():
    log.info("Collecting BLAST XML data from run_blast.py output...")
    df = collect_blast_data()
    log.info("Collected %d rows", len(df))
    if df.empty:
        log.warning("No data collected. Run run_blast.py first.")
        return

    df = df[df["ref_org"] != "tabaco"].copy()

    df["enzyme_label"] = df["enzyme"].map(ENZYME_LABELS)
    df["row_label"] = df["enzyme_label"] + " (" + df["ref_label"] + ")"
    df["col_label"] = df["target_organism"].map(ORG_DISPLAY)

    pivot = df.pivot_table(
        index="row_label", columns="col_label",
        values="best_bit_score", aggfunc="max"
    )

    fig, ax = plt.subplots(figsize=(6, max(5, len(pivot) * 0.45)))

    colors = [(0.05, 0.05, 0.6), (0.95, 0.95, 0.95), (0.6, 0.02, 0.02)]
    cmap = LinearSegmentedColormap.from_list("blue_red", colors)

    sns.heatmap(pivot, annot=True, fmt=".0f", cmap=cmap,
                xticklabels=True, yticklabels=True,
                linewidths=0.5, center=None,
                cbar_kws={"label": "Best bit score"},
                ax=ax, annot_kws={"fontsize": 7})

    ax.set_title("BLAST bit score", fontsize=13, fontweight="bold")
    ax.set_ylabel("Query enzyme (reference species)", fontsize=11)
    ax.set_xlabel("Target genome", fontsize=11)
    ax.tick_params(axis="x", rotation=45)
    ax.tick_params(axis="y", rotation=0)

    plt.tight_layout()
    fig.savefig(FIGURES / "blast_heatmap.png", dpi=300, bbox_inches="tight")
    log.info("Saved %s", FIGURES / "blast_heatmap.png")
    plt.close(fig)


if __name__ == "__main__":
    main()
