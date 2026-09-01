from pathlib import Path
from zipfile import ZipFile

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import norm as norm_dist

from data import TRANSCRIPTOME_PATH, PAPER_PATH

TUCUNACA_FILE = f"{TRANSCRIPTOME_PATH}/tucunaca/tucunaca.tab"
CAUPURI_FILE = f"{TRANSCRIPTOME_PATH}/caupuri/caupuri.tab"
TUCUNACA_ZIP = f"{TRANSCRIPTOME_PATH}/tucunaca/tucunaca.zip"
CAUPURI_ZIP = f"{TRANSCRIPTOME_PATH}/caupuri/caupuri.zip"
OUTPUT_PNG = f"{PAPER_PATH}/figures/volcano_plot_tucunaca_vs_caupuri.png"
OUTPUT_PNG_STRICT = f"{PAPER_PATH}/figures/volcano_plot_tucunaca_vs_caupuri_strictosidine.png"
OUTPUT_PNG_TOP10 = f"{PAPER_PATH}/figures/volcano_plot_tucunaca_vs_caupuri_top10.png"
OUTPUT_CSV = f"{PAPER_PATH}/tables/significant_genes.csv"
BRIDGE_FILE = f"{TRANSCRIPTOME_PATH}/gene_id_bridge.csv"
TUCUNACA_BLAST = f"{PAPER_PATH}/tables/tucunaca_blast_abundance.csv"
CAUPURI_BLAST = f"{PAPER_PATH}/tables/caupuri_blast_abundance.csv"
TUCUNACA_PROT = f"{TRANSCRIPTOME_PATH}/../genome/tucunaca/non_phased/prot.fasta"
CAUPURI_PROT = f"{TRANSCRIPTOME_PATH}/../genome/caupuri/non_phased/prot.fasta"
TUCUNACA_PROT_PHASED = f"{TRANSCRIPTOME_PATH}/../genome/tucunaca/phased/prot.fasta"
CAUPURI_PROT_PHASED = f"{TRANSCRIPTOME_PATH}/../genome/caupuri/phased/prot.fasta"
STRICTO_BLAST = f"{TRANSCRIPTOME_PATH}/../proteins/stricto/analysis/genomic/blast"
STRICTO_PHYLO = f"{TRANSCRIPTOME_PATH}/../proteins/stricto/analysis/phylo"
GENOME_PATH = f"{TRANSCRIPTOME_PATH}/../genome"


def _load_emapper(zip_path: str, member: str) -> pd.DataFrame:
    with ZipFile(zip_path) as zf:
        with zf.open(member) as f:
            content = f.read().decode()
    lines = [l for l in content.splitlines() if not l.startswith("##")]
    header = lines[0].lstrip("#").split("\t")
    data = [l.split("\t") for l in lines[1:]]
    df = pd.DataFrame(data, columns=header)
    needed = {"query", "seed_ortholog"}
    if not needed.issubset(set(df.columns)):
        raise ValueError(f"emapper file missing columns")
    df = df[["query", "seed_ortholog"]].drop_duplicates()
    df["query"] = df["query"].str.replace(r"\..*$", "", regex=True)
    return df


def build_gene_id_bridge() -> pd.DataFrame:
    tuc = _load_emapper(TUCUNACA_ZIP, "Tucunaca.emapper.annotations").drop_duplicates(subset=["query"])
    cau = _load_emapper(CAUPURI_ZIP, "Caupuri.emapper.annotations").drop_duplicates(subset=["query"])

    merged = tuc.merge(cau, on="seed_ortholog", suffixes=("_tucunaca", "_caupuri"))
    merged = merged[["query_tucunaca", "query_caupuri"]].rename(
        columns={"query_tucunaca": "tucunaca_id", "query_caupuri": "caupuri_id"}
    )
    merged = merged.drop_duplicates(subset=["tucunaca_id"])
    merged = merged.drop_duplicates(subset=["caupuri_id"])
    merged = merged.reset_index(drop=True)
    merged.to_csv(BRIDGE_FILE, index=False)
    return merged


def load_expression_table(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    required = {"Gene ID", "TPM"}
    if not required.issubset(set(df.columns)):
        raise ValueError(f"{path.name} is missing required columns")
    df = df[["Gene ID", "TPM"]].copy()
    df["TPM"] = pd.to_numeric(df["TPM"], errors="coerce")
    return df


def _ids_from_fasta(fasta_path: str, prefix: str) -> set:
    ids = set()
    with open(fasta_path) as f:
        for line in f:
            if line.startswith(">") and "strictosidine" in line.lower():
                gid = line[1:].split("|")[0].split(".")[0]
                if gid.startswith(prefix):
                    ids.add(gid)
    return ids


def _ids_from_blast_csv(csv_path: str) -> set:
    if not Path(csv_path).exists():
        return set()
    df = pd.read_csv(csv_path, sep=";", dtype=str)
    if "hit_id" not in df.columns:
        return set()
    return {g.split(".")[0] for g in df["hit_id"].dropna() if g.split(".")[0].startswith(("TucunacaG", "CaupuriG"))}


def _ids_from_fasta_headers(fasta_path: str, prefix: str) -> set:
    if not Path(fasta_path).exists():
        return set()
    ids = set()
    with open(fasta_path) as f:
        for line in f:
            if line.startswith(">"):
                gid = line[1:].split("|")[0].split(".")[0]
                if gid.startswith(prefix):
                    ids.add(gid)
    return ids


def _load_all_gene_id_maps() -> dict:
    import subprocess
    mapping = {}
    base = Path(f"{TRANSCRIPTOME_PATH}/../proteins")
    out = subprocess.run(
        ["find", str(base), "-name", "gene_id_map.csv"],
        capture_output=True, text=True
    )
    for fpath in out.stdout.strip().splitlines():
        if not fpath:
            continue
        with open(fpath) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("phased"):
                    continue
                parts = line.split(",")
                if len(parts) < 2:
                    continue
                pid = parts[0].split(".")[0]
                nid = parts[1].split(".")[0]
                if nid not in mapping:
                    mapping[nid] = pid
    return mapping


def load_strictosidine_ids():
    tuc_ids = set()
    cau_ids = set()

    # Source 1: non-phased prot.fasta (direct strictosidine header match)
    tuc_ids |= _ids_from_fasta(TUCUNACA_PROT, "TucunacaG")
    cau_ids |= _ids_from_fasta(CAUPURI_PROT, "CaupuriG")

    # Source 2: BLAST abundance CSVs (STRICTOSIDINE SYNTHASE in hit_description)
    tuc = pd.read_csv(TUCUNACA_BLAST)
    cau = pd.read_csv(CAUPURI_BLAST)
    tuc_ids |= {g.split(".")[0] for g in tuc.loc[
        tuc["hit_description"].str.contains("STRICTOSIDINE SYNTHASE", case=False, na=False),
        "non_phased_id"
    ]}
    cau_ids |= {g.split(".")[0] for g in cau.loc[
        cau["hit_description"].str.contains("STRICTOSIDINE SYNTHASE", case=False, na=False),
        "non_phased_id"
    ]}

    # Source 3: genomic blast all_seqs.fasta
    for species in ("tucunaca", "caupuri"):
        all_seqs = Path(f"{STRICTO_BLAST}/{species}/all_seqs.fasta")
        if all_seqs.exists():
            prefix = "TucunacaG" if species == "tucunaca" else "CaupuriG"
            target = tuc_ids if species == "tucunaca" else cau_ids
            target |= _ids_from_fasta_headers(str(all_seqs), prefix)

    # Source 4: phylo alignment files
    for species in ("tucunaca", "caupuri"):
        align = Path(f"{STRICTO_PHYLO}/{species}/clean_alignment.fasta")
        if align.exists():
            prefix = "TucunacaG" if species == "tucunaca" else "CaupuriG"
            target = tuc_ids if species == "tucunaca" else cau_ids
            target |= _ids_from_fasta_headers(str(align), prefix)

    # Source 5: blast_filtered.csv for each reference species
    for tuc_dir in Path(f"{STRICTO_BLAST}/tucunaca").iterdir():
        csv_path = tuc_dir / "stricto_*_blast_filtered.csv"
        for f in tuc_dir.glob("*_blast_filtered.csv"):
            tuc_ids |= _ids_from_blast_csv(str(f))
    for cau_dir in Path(f"{STRICTO_BLAST}/caupuri").iterdir():
        for f in cau_dir.glob("*_blast_filtered.csv"):
            cau_ids |= _ids_from_blast_csv(str(f))

    # Source 6: phased prot.fasta → map to non-phased via gene_id_map
    phased_map = _load_all_gene_id_maps()
    phased_tuc = _ids_from_fasta(TUCUNACA_PROT_PHASED, "TucunacaPG")
    phased_cau = _ids_from_fasta(CAUPURI_PROT_PHASED, "CaupuriPG")
    rev_map = {v: k for k, v in phased_map.items()}
    for pid in phased_tuc:
        if pid in rev_map:
            tuc_ids.add(rev_map[pid])
    for pid in phased_cau:
        if pid in rev_map:
            cau_ids.add(rev_map[pid])

    return tuc_ids, cau_ids


def load_gene_descriptions(fasta_path: str) -> dict:
    desc = {}
    with open(fasta_path) as f:
        for line in f:
            if line.startswith(">"):
                parts = line[1:].strip().split("|", 2)
                gid = parts[0].split(".")[0]
                name = parts[1] if len(parts) > 1 else ""
                desc[gid] = name
    return desc


def load_phased_map():
    import subprocess
    tuc_map, cau_map = {}, {}
    base = Path(f"{TRANSCRIPTOME_PATH}/../proteins")
    out = subprocess.run(
        ["find", str(base), "-name", "gene_id_map.csv"],
        capture_output=True, text=True
    )
    for fpath in out.stdout.strip().splitlines():
        if not fpath:
            continue
        with open(fpath) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("phased"):
                    continue
                parts = line.split(",")
                if len(parts) < 2:
                    continue
                pid, nid = parts[0].split(".")[0], parts[1].split(".")[0]
                if nid not in tuc_map and nid.startswith("Tucunaca"):
                    tuc_map[nid] = pid
                elif nid not in cau_map and nid.startswith("Caupuri"):
                    cau_map[nid] = pid
    return tuc_map, cau_map


def make_volcano_plot() -> Path:
    tucunaca = load_expression_table(TUCUNACA_FILE)
    caupuri = load_expression_table(CAUPURI_FILE)

    if Path(BRIDGE_FILE).exists():
        bridge = pd.read_csv(BRIDGE_FILE)
    else:
        bridge = build_gene_id_bridge()

    tuc_mapped = tucunaca.merge(bridge, left_on="Gene ID", right_on="tucunaca_id")
    cau_mapped = caupuri.merge(bridge, left_on="Gene ID", right_on="caupuri_id")

    merged = tuc_mapped.merge(
        cau_mapped[["tucunaca_id", "Gene ID", "TPM"]],
        on="tucunaca_id",
        suffixes=("_tucunaca", "_caupuri"),
    ).rename(columns={
        "Gene ID_tucunaca": "original_tucunaca_id",
        "Gene ID_caupuri": "original_caupuri_id",
    })

    merged["tpm_t"] = pd.to_numeric(merged["TPM_tucunaca"], errors="coerce").values
    merged["tpm_c"] = pd.to_numeric(merged["TPM_caupuri"], errors="coerce").values

    merged["log2fc"] = np.log2(merged["tpm_t"] / merged["tpm_c"])
    merged["mean_expr"] = (merged["tpm_t"] + merged["tpm_c"]) / 2
    merged["log10_mean"] = np.log10(merged["mean_expr"])

    merged = merged.replace([np.inf, -np.inf], np.nan).dropna(subset=["log2fc", "log10_mean"])

    order = np.argsort(merged["log10_mean"])
    n_bins = 50
    bins = np.array_split(order, n_bins)

    predicted_sd = np.empty(len(merged))
    for b in bins:
        sd = np.std(merged["log2fc"].iloc[b].values)
        predicted_sd[b] = max(sd, 0.05)

    merged["neg_log10_p"] = -np.log10(
        np.clip(2 * norm_dist.sf(np.abs(merged["log2fc"].values) / predicted_sd), 1e-300, 1.0)
    )

    fc_cutoff = 1.0
    pval_cutoff = -np.log10(0.05)
    significant = (merged["log2fc"].abs() >= fc_cutoff) & (merged["neg_log10_p"] >= pval_cutoff)
    up = significant & (merged["log2fc"] > 0)
    down = significant & (merged["log2fc"] < 0)

    if significant.any():
        sig_df = merged.loc[significant, [
            "tucunaca_id", "original_tucunaca_id", "original_caupuri_id",
            "log2fc", "neg_log10_p",
        ]].copy()
        sig_df = sig_df.sort_values("neg_log10_p", ascending=False)
        sig_df.to_csv(OUTPUT_CSV, index=False)
        print(f"Saved {len(sig_df)} significant genes to {OUTPUT_CSV}")

    plt.figure(figsize=(9, 7))
    plt.scatter(
        merged["log2fc"], merged["neg_log10_p"],
        c="lightgray", s=10, alpha=0.4,
    )
    if up.any():
        plt.scatter(
            merged.loc[up, "log2fc"],
            merged.loc[up, "neg_log10_p"],
            c="#2d8a2d", s=18, alpha=0.85, label=f"Higher in Tucunaca (n={up.sum()})",
        )
    if down.any():
        plt.scatter(
            merged.loc[down, "log2fc"],
            merged.loc[down, "neg_log10_p"],
            c="#2040c0", s=18, alpha=0.85, label=f"Higher in Caupuri (n={down.sum()})",
        )

    plt.axvline(-fc_cutoff, color="gray", linestyle="--", linewidth=0.8)
    plt.axvline(fc_cutoff, color="gray", linestyle="--", linewidth=0.8)
    plt.axhline(pval_cutoff, color="gray", linestyle="--", linewidth=0.8)

    plt.xlabel("log2 fold change (Tucunaca / Caupuri)")
    plt.ylabel("-log10(p-value)")
    plt.title("Volcano plot: Tucunaca vs Caupuri")
    plt.legend(loc="upper right", markerscale=2)
    plt.tight_layout()
    plt.savefig(OUTPUT_PNG, dpi=300)
    plt.close()

    # ---- strictosidine volcano plot with labels ----
    strict_tuc, strict_cau = load_strictosidine_ids()
    strict_bridge = bridge[bridge["tucunaca_id"].isin(strict_tuc) | bridge["caupuri_id"].isin(strict_cau)]
    strict_in_data = merged["tucunaca_id"].isin(strict_bridge["tucunaca_id"])

    plt.figure(figsize=(9, 7))
    plt.scatter(
        merged["log2fc"], merged["neg_log10_p"],
        c="lightgray", s=10, alpha=0.4,
    )
    if strict_in_data.any():
        strict_pts = merged.loc[strict_in_data]
        plt.scatter(
            strict_pts["log2fc"], strict_pts["neg_log10_p"],
            c="red", s=40, alpha=0.85, edgecolors="black", linewidths=0.5,
            label=f"Strictosidine synthase-like (n={strict_in_data.sum()})",
        )
        for _, r in strict_pts.iterrows():
            label_id = r["original_caupuri_id"].split(".")[0] if r["log2fc"] < 0 else r["original_tucunaca_id"].split(".")[0]
            plt.text(
                r["log2fc"], r["neg_log10_p"], label_id,
                fontsize=6, alpha=0.9, ha="left", va="bottom",
            )
    plt.axvline(-fc_cutoff, color="gray", linestyle="--", linewidth=0.8)
    plt.axvline(fc_cutoff, color="gray", linestyle="--", linewidth=0.8)
    plt.axhline(pval_cutoff, color="gray", linestyle="--", linewidth=0.8)
    plt.xlabel("log2 fold change (Tucunaca / Caupuri)")
    plt.ylabel("-log10(p-value)")
    plt.title("Volcano plot: STRICTOSIDINE SYNTHASE-LIKE hits highlighted")
    plt.legend(loc="upper right", markerscale=2)
    plt.tight_layout()
    plt.savefig(OUTPUT_PNG_STRICT, dpi=300)
    plt.close()

    # ---- top 10 volcano plot with gene table ----
    tuc_desc = load_gene_descriptions(TUCUNACA_PROT)
    cau_desc = load_gene_descriptions(CAUPURI_PROT)
    tuc_phased, cau_phased = load_phased_map()

    top10 = merged.loc[significant].nlargest(10, "neg_log10_p")
    top10 = top10.sort_values("neg_log10_p", ascending=False)

    fig, (ax_plot, ax_table) = plt.subplots(1, 2, figsize=(18, 9),
                                            gridspec_kw={"width_ratios": [1.8, 1]})

    ax_plot.scatter(
        merged["log2fc"], merged["neg_log10_p"],
        c="lightgray", s=10, alpha=0.4,
    )
    if up.any():
        ax_plot.scatter(
            merged.loc[up, "log2fc"],
            merged.loc[up, "neg_log10_p"],
            c="#2d8a2d", s=18, alpha=0.85,
        )
    if down.any():
        ax_plot.scatter(
            merged.loc[down, "log2fc"],
            merged.loc[down, "neg_log10_p"],
            c="#2040c0", s=18, alpha=0.85,
        )
    ax_plot.scatter(
        top10["log2fc"], top10["neg_log10_p"],
        c="red", s=40, alpha=0.85, edgecolors="black", linewidths=0.5,
        zorder=5,
    )
    for _, r in top10.iterrows():
        label_id = r["original_caupuri_id"].split(".")[0] if r["log2fc"] < 0 else r["original_tucunaca_id"].split(".")[0]
        ax_plot.text(
            r["log2fc"], r["neg_log10_p"], label_id,
            fontsize=7, alpha=0.9, ha="left", va="bottom",
        )
    ax_plot.axvline(-fc_cutoff, color="gray", linestyle="--", linewidth=0.8)
    ax_plot.axvline(fc_cutoff, color="gray", linestyle="--", linewidth=0.8)
    ax_plot.axhline(pval_cutoff, color="gray", linestyle="--", linewidth=0.8)
    ax_plot.set_xlabel("log2 fold change (Tucunaca / Caupuri)")
    ax_plot.set_ylabel("-log10(p-value)")
    ax_plot.set_title("Volcano plot: Tucunaca vs Caupuri")

    table_data = []
    for _, r in top10.iterrows():
        tuc_id = r["original_tucunaca_id"].split(".")[0]
        cau_id = r["original_caupuri_id"].split(".")[0]
        tuc_ph = tuc_phased.get(tuc_id, "")
        cau_ph = cau_phased.get(cau_id, "")
        tp = tuc_desc.get(tuc_id, "")
        cp = cau_desc.get(cau_id, "")
        table_data.append([
            tuc_id, tuc_ph, tp,
            cau_id, cau_ph, cp,
            f'{r["log2fc"]:.2f}',
            f'{r["neg_log10_p"]:.2f}',
        ])

    col_labels = [
        "Tucunaca\nG-ID", "PG-ID", "Predicted gene",
        "Caupuri\nG-ID", "PG-ID", "Predicted gene",
        "log2FC", "-log10(p)",
    ]
    ax_table.axis("off")
    tbl = ax_table.table(
        cellText=table_data,
        colLabels=col_labels,
        cellLoc="left",
        loc="center",
    )
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(7)
    tbl.scale(1, 1.4)
    for key, cell in tbl.get_celld().items():
        if key[0] == 0:
            cell.set_text_props(weight="bold")
        cell.PAD = 0.03
    ax_table.set_title("Top 10 significant genes", fontsize=11, fontweight="bold", pad=12)

    plt.tight_layout()
    plt.savefig(OUTPUT_PNG_TOP10, dpi=300, bbox_inches="tight")
    plt.close()

    return OUTPUT_PNG


if __name__ == "__main__":
    output_path = make_volcano_plot()
    print(f"Saved volcano plot to {output_path}")
