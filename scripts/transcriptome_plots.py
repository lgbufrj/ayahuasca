from pathlib import Path
from statistics import NormalDist

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from data import TRANSCRIPTOME_PATH, PAPER_PATH

TUCUNACA_FILE = f"{TRANSCRIPTOME_PATH}/tucunaca/tucunaca.tab"
CAUPURI_FILE = f"{TRANSCRIPTOME_PATH}/caupuri/caupuri.tab"
OUTPUT_FILE = f"{PAPER_PATH}/figures/volcano_plot_tucunaca_vs_caupuri.png"

def load_expression_table(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    required = {"Gene ID", "TPM"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"{path.name} is missing required columns: {sorted(missing)}")

    df = df[["Gene ID", "TPM"]].copy()
    df["TPM"] = pd.to_numeric(df["TPM"], errors="coerce")
    df = df.dropna(subset=["TPM"])
    df = df.rename(columns={"Gene ID": "Gene ID"})
    return df


def make_volcano_plot() -> Path:
    tucunaca = load_expression_table(TUCUNACA_FILE)
    caupuri = load_expression_table(CAUPURI_FILE)

    merged = tucunaca.merge(caupuri, on="Gene ID", suffixes=("_tucunaca", "_caupuri"))

    merged["tpm_tucunaca"] = np.maximum(merged["TPM_tucunaca"], 1e-6)
    merged["tpm_caupuri"] = np.maximum(merged["TPM_caupuri"], 1e-6)
    merged["log2fc"] = np.log2(merged["tpm_tucunaca"] / merged["tpm_caupuri"])

    # Approximate p-values from a simple z-score model for a volcano plot.
    z_score = np.abs(merged["log2fc"]) / 0.5
    p_value = 2 * np.array([NormalDist().sf(x) for x in z_score], dtype=float)
    merged["neg_log10_p"] = -np.log10(np.clip(p_value, 1e-300, 1.0))

    significant = (merged["log2fc"].abs() >= 1.0) & (merged["neg_log10_p"] >= -np.log10(0.05))

    plt.figure(figsize=(8, 6))
    plt.scatter(merged["log2fc"], merged["neg_log10_p"], c="lightgray", s=20, alpha=0.7)
    if significant.any():
        plt.scatter(
            merged.loc[significant, "log2fc"],
            merged.loc[significant, "neg_log10_p"],
            c="red",
            s=35,
            alpha=0.85,
            edgecolors="black",
            label="Significant",
        )

    for _, row in merged.loc[significant].head(20).iterrows():
        plt.text(row["log2fc"], row["neg_log10_p"], row["Gene ID"], fontsize=7, alpha=0.8)

    plt.axvline(-1, color="blue", linestyle="--", linewidth=1)
    plt.axvline(1, color="blue", linestyle="--", linewidth=1)
    plt.axhline(-np.log10(0.05), color="green", linestyle="--", linewidth=1)

    plt.xlabel("log2 fold change (Tucunaca / Caupuri)")
    plt.ylabel("-log10(p-value)")
    plt.title("Volcano plot: Tucunaca vs Caupuri")
    plt.legend(loc="upper right")
    plt.tight_layout()
    plt.savefig(OUTPUT_FILE, dpi=300)
    plt.close()

    return OUTPUT_FILE


if __name__ == "__main__":
    output_path = make_volcano_plot()
    print(f"Saved volcano plot to {output_path}")
