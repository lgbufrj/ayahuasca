import os, sys, warnings, subprocess, shutil, json
from pathlib import Path
from itertools import combinations

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

from Bio import SeqIO, Align
from Bio.PDB import PDBParser, PDBIO, MMCIFParser
from Bio.PDB.Polypeptide import protein_letters_3to1 as _3to1_dict
# 3to1 mapping for residue conversion
_3TO1 = _3to1_dict.copy()
_3TO1['MSE'] = 'M'  # selenomethionine
def _3to1(rname):
    return _3to1_dict.get(rname, 'X')
from tmtools import tm_align

warnings.filterwarnings('ignore')

BASE = Path("/home/pedro/Desktop/projects/ayahuasca")
PROJECT = BASE / "pathways" / "harmine_biosynthesis"
PROTEINS = PROJECT / "proteins"
PAPER = PROJECT / "paper"
FIGURES = PAPER / "figures"

STRICTO = "stricto"
S = PROTEINS / STRICTO
REF = S / "reference"
ANALYSIS = S / "analysis"
STRUCT = ANALYSIS / "structural"
BOLTZ_DIR = STRUCT / "boltz"
PHYLO = ANALYSIS / "phylo"
GENOMIC = ANALYSIS / "genomic" / "blast"

PRED_DIR = STRUCT / "phylo_predictions"
YAML_DIR = PRED_DIR / "inputs"
PDB_DIR = PRED_DIR / "pdbs"
BATCH_DIR = PRED_DIR / "batch"
GEN_DIR = PRED_DIR / "generated"

for d in [FIGURES, YAML_DIR, PDB_DIR, BATCH_DIR, GEN_DIR]:
    d.mkdir(parents=True, exist_ok=True)

sns.set_style("whitegrid")
sns.set_context("paper", font_scale=1.2)


def extract_coords_and_seq(path):
    parser = PDBParser(QUIET=True)
    try:
        structure = parser.get_structure("prot", str(path))
    except Exception:
        try:
            parser = MMCIFParser(QUIET=True)
            structure = parser.get_structure("prot", str(path))
        except Exception:
            return None, None
    coords, seq = [], []
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.id[0] == " ":
                    if "CA" in residue:
                        coords.append(residue["CA"].coord)
                    seq.append(_3to1(residue.resname))
    if len(coords) < 10:
        return None, None
    return np.array(coords, dtype=float), "".join(seq)


def clean_label(name):
    return name.split(" ")[0]


def is_gene(name):
    return name.startswith("CaupuriG") or name.startswith("TucunacaG")


def pid(seq1, seq2):
    if not seq1 or not seq2:
        return 0.0
    a = Align.PairwiseAligner()
    a.mode = "global"
    a.open_gap_score = -10
    a.extend_gap_score = -0.5
    try:
        aln = a.align(seq1, seq2)[0]
        pairs = list(zip(str(aln[0]), str(aln[1])))
        ident = sum(1 for x, y in pairs if x == y and x != "-")
        aln_len = sum(1 for x, y in pairs if x != "-" or y != "-")
        return ident / aln_len * 100 if aln_len > 0 else 0.0
    except Exception:
        return 0.0


def tmalign(c1, s1, c2, s2):
    if c1 is None or c2 is None or len(c1) < 10 or len(c2) < 10:
        return None, None
    try:
        r = tm_align(c1, c2, s1, s2)
        return r.tm_norm_chain1, r.rmsd
    except Exception:
        return None, None


def read_all_sequences():
    seqs = {}
    for p in [
        PHYLO / "caupuri" / "clean_alignment.fasta",
        PHYLO / "tucunaca" / "clean_alignment.fasta",
        GENOMIC / "caupuri" / "all_seqs.fasta",
        GENOMIC / "tucunaca" / "all_seqs.fasta",
    ]:
        if p.exists():
            for rec in SeqIO.parse(str(p), "fasta"):
                key = clean_label(rec.id)
                seqs[key] = str(rec.seq).replace("-", "")
    for p in [
        REF / "tabaco" / "structure" / "stricto_A0A1U7VHK2.fasta",
        REF / "arabidopsis" / "structure" / "stricto_P94111.fasta",
        REF / "catharanthus" / "structure" / "stricto_P18417.fasta",
    ]:
        if p.exists():
            for rec in SeqIO.parse(str(p), "fasta"):
                key = clean_label(rec.id)
                seqs[key] = str(rec.seq)
    return seqs


def collect_available_structures():
    """Find all existing PDB/CIF files."""
    structs = {}
    candidates = [
        ("tabaco_boltz", BOLTZ_DIR / "tabaco" / "r1" / "predictions" / "tabaco_stricto_A0A1U7VHK2_r1" / "tabaco_stricto_A0A1U7VHK2_r1_model_0.cif"),
        ("caupuri_boltz_catha", BOLTZ_DIR / "caupuri" / "catharanthus" / "rr" / "predictions" / "caupuri_catharanthus_stricto_P18417_rr" / "caupuri_catharanthus_stricto_P18417_rr_model_0.cif"),
        ("caupuri_boltz_tabaco", BOLTZ_DIR / "caupuri" / "tabaco" / "rr" / "predictions" / "caupuri_tabaco_stricto_A0A1U7VHK2_rr" / "caupuri_tabaco_stricto_A0A1U7VHK2_rr_model_0.cif"),
        ("tucunaca_boltz_arab", BOLTZ_DIR / "tucunaca" / "arabidopsis" / "r1" / "predictions" / "tucunaca_arabidopsis_stricto_P94111_r1" / "tucunaca_arabidopsis_stricto_P94111_r1_model_0.cif"),
        ("arabidopsis_boltz", BOLTZ_DIR / "arabidopsis" / "r1" / "predictions" / "arabidopsis_stricto_P94111_r1" / "arabidopsis_stricto_P94111_r1_model_0.cif"),
        ("catharanthus_boltz", BOLTZ_DIR / "catharanthus" / "r1" / "predictions" / "catharanthus_stricto_P18417_r1" / "catharanthus_stricto_P18417_r1_model_0.cif"),
        ("tucunaca_alphafold", STRUCT / "alphafold" / "tucunaca" / "best_model.pdb"),
        ("tucunaca_p2rank", STRUCT / "p2rank" / "tucunaca" / "output" / "visualizations" / "data" / "best_model.pdb"),
        ("arabidopsis_ref", REF / "arabidopsis" / "structure" / "stricto_P94111.pdb"),
        ("catharanthus_ref", REF / "catharanthus" / "structure" / "P18417.pdb"),
        ("catharanthus_xray", REF / "catharanthus" / "structure" / "x_ray_6ZEA.pdb"),
    ]
    for label, path in candidates:
        if path.exists():
            structs[label] = str(path)
    return structs


def generate_boltz_yaml(seq, label, out):
    import yaml
    data = {
        "version": 1,
        "sequences": [{"protein": {"id": ["A"], "sequence": seq}}]
    }
    with open(out, "w") as f:
        yaml.dump(data, f, default_flow_style=False, sort_keys=False)


def write_batch_script(yaml_dir, batch_dir):
    script = batch_dir / "run_boltz_batch.sh"
    with open(script, "w") as f:
        f.write("#!/bin/bash\n")
        f.write("# Run Boltz predictions for all stricto phylo sequences\n")
        f.write("# Usage: bash run_boltz_batch.sh\n\n")
        f.write("mkdir -p predictions\n")
        for y in sorted(yaml_dir.glob("*.yaml")):
            name = y.stem
            f.write(f"boltz predict {y} --out_dir predictions/{name} --accelerator gpu --recycling_steps 3 --sampling_steps 200 --diffusion_samples 1 --use_msa_server\n")
        f.write("\necho 'All predictions submitted'\n")
    os.chmod(script, 0o755)
    return script


def assign_structures_with_homology(all_seqs, available_structs):
    """Assign best-matching available structure to each sequence."""
    struct_seqs = {}
    for label, path in available_structs.items():
        _, s = extract_coords_and_seq(path)
        if s:
            struct_seqs[label] = s

    assignments = {}
    for gname, gseq in all_seqs.items():
        if gseq and len(gseq) > 20:
            best_pid = 0
            best_match = None
            for slabel, sseq in struct_seqs.items():
                p = pid(gseq, sseq)
                if p > best_pid:
                    best_pid = p
                    best_match = slabel
            if best_match:
                assignments[gname] = {
                    "struct_label": best_match,
                    "path": available_structs[best_match],
                    "pid": best_pid
                }
    return assignments, struct_seqs


def main():
    print("=" * 70)
    print("STRICTO Phylo TM-align Analysis")
    print("=" * 70)

    print("\n[1/5] Collecting sequences...")
    all_seqs = read_all_sequences()
    gene_seqs = {k: v for k, v in all_seqs.items() if is_gene(k)}
    ref_seqs = {k: v for k, v in all_seqs.items() if not is_gene(k)}
    print(f"  Gene sequences: {len(gene_seqs)}")
    print(f"  Reference sequences: {len(ref_seqs)}")

    print("\n[2/5] Collecting available structures...")
    avail = collect_available_structures()
    print(f"  Available structures: {len(avail)}")
    for l, p in sorted(avail.items()):
        print(f"    {l}")

    gen_seqs = {k: v for k, v in gene_seqs.items()}
    tabaco_seq = ref_seqs.get("tabaco_stricto_A0A1U7VHK2",
                              ref_seqs.get("tr|A0A1U7VHK2|A0A1U7VHK2_NICSY", ""))

    print(f"\n[3/5] Generating Boltz inputs and batch script...")
    yaml_count = 0
    for label, seq in all_seqs.items():
        safe = "".join(c if c.isalnum() or c in "-_" else "_" for c in clean_label(label))[:40]
        yp = YAML_DIR / f"{safe}.yaml"
        if not yp.exists():
            generate_boltz_yaml(seq, label, yp)
            yaml_count += 1
    batch_script = write_batch_script(YAML_DIR, BATCH_DIR)
    print(f"  Generated {yaml_count} new Boltz YAML inputs in {YAML_DIR}")
    print(f"  Batch script: {batch_script}")

    print(f"\n  Assigning structures via homology matching...")
    assignments, struct_seqs = assign_structures_with_homology(all_seqs, avail)
    print(f"  Matched {len(assignments)}/{len(all_seqs)} sequences to existing structures")

    struct_data = {}
    for gname, gseq in gene_seqs.items():
        if gname in assignments:
            a = assignments[gname]
            coords, seq = extract_coords_and_seq(a["path"])
            if coords is not None:
                struct_data[gname] = {"coords": coords, "seq": seq, "path": a["path"], "pid": a["pid"]}
                print(f"    {gname:40s} -> {a['struct_label']:25s} ({a['pid']:.0f}% identity)")

    if tabaco_seq:
        for tlabel in ["tabaco_boltz", "tabaco_ref"]:
            if tlabel in avail:
                tc, ts = extract_coords_and_seq(avail[tlabel])
                if tc is not None:
                    struct_data["tabaco_reference"] = {"coords": tc, "seq": ts, "path": avail[tlabel], "pid": 100}
                    print(f"    {'tabaco_reference':40s} -> {tlabel:25s} (reference)")
                    break

    print(f"\n  Usable structures: {len(struct_data)}")
    if len(struct_data) < 2:
        print("  [FATAL] Need >= 2 structures. Using reference structures as fallback.")
        for label, path in avail.items():
            if label not in struct_data:
                coords, seq = extract_coords_and_seq(path)
                if coords is not None:
                    struct_data[label] = {"coords": coords, "seq": seq, "path": path, "pid": 100}
        if len(struct_data) < 2:
            print("  Still insufficient. Exiting.")
            sys.exit(1)

    print(f"\n[4/5] Running TM-align pairwise ({len(struct_data)} structures)...")
    results = []
    labels = list(struct_data.keys())
    for s1, s2 in combinations(labels, 2):
        d1, d2 = struct_data[s1], struct_data[s2]
        tm, rmsd = tmalign(d1["coords"], d1["seq"], d2["coords"], d2["seq"])
        if tm is not None:
            results.append({
                "struct1": s1, "struct2": s2,
                "tm_score": round(tm, 5), "rmsd": round(rmsd, 4)
            })
            print(f"    {s1:40s} vs {s2:40s}  TM={tm:.4f}  RMSD={rmsd:.2f}A")

    results_df = pd.DataFrame(results)
    csv_path = PAPER / "tables" / "stricto_phylo_tmalign.csv"
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    results_df.to_csv(csv_path, index=False)
    print(f"\n  {len(results_df)} comparisons -> {csv_path}")

    print(f"\n[5/5] Generating plots...")

    # --- TM-score heatmap ---
    gene_only_results = results_df[
        results_df["struct1"].apply(is_gene) & results_df["struct2"].apply(is_gene)
    ]
    pivot_target = gene_only_results if len(gene_only_results) >= 2 else results_df
    if not pivot_target.empty:
        pivot = pivot_target.pivot_table(index="struct1", columns="struct2",
                                         values="tm_score", aggfunc="first")
        all_l = sorted(set(list(pivot.index) + list(pivot.columns)))
        full = pd.DataFrame(index=all_l, columns=all_l, dtype=float)
        for i in all_l:
            full.loc[i, i] = 1.0
            for j in all_l:
                if i != j:
                    if i in pivot.index and j in pivot.columns and pd.notna(pivot.loc[i, j]):
                        full.loc[i, j] = pivot.loc[i, j]
                    elif j in pivot.index and i in pivot.columns and pd.notna(pivot.loc[j, i]):
                        full.loc[i, j] = pivot.loc[j, i]
        display = full.dropna(how="all").dropna(how="all", axis=1)
        short = [l.split("|")[0][:16] for l in display.index]
        fig, ax = plt.subplots(figsize=(max(10, len(display)*0.5), max(8, len(display)*0.45)))
        sns.heatmap(display, annot=True, fmt=".3f", cmap="YlOrRd",
                    xticklabels=short, yticklabels=short,
                    vmin=0, vmax=1, square=True, linewidths=0.5,
                    cbar_kws={"label": "TM-score"}, ax=ax,
                    annot_kws={"fontsize": 6})
        ax.set_title("Stricto Phylo TM-score Matrix", fontsize=13, fontweight="bold")
        plt.xticks(rotation=45, ha="right", fontsize=7)
        plt.yticks(fontsize=7)
        plt.tight_layout()
        fig.savefig(FIGURES / "phylo_tmscore_heatmap.png", dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(f"  Saved phylo_tmscore_heatmap.png")

    # --- RMSD vs tabaco barplot ---
    tabaco_rows = results_df[
        results_df["struct1"].str.contains("tabaco", case=False) |
        results_df["struct2"].str.contains("tabaco", case=False)
    ].dropna(subset=["rmsd"]).copy()
    if not tabaco_rows.empty:
        def other(r):
            return r["struct2"] if "tabaco" in r["struct1"].lower() else r["struct1"]
        tabaco_rows["gene"] = tabaco_rows.apply(other, axis=1)
        tabaco_rows = tabaco_rows.sort_values("gene")
        fig, ax = plt.subplots(figsize=(max(8, len(tabaco_rows)*0.6), 6))
        colors = sns.color_palette("Set2", n_colors=len(tabaco_rows))
        x = range(len(tabaco_rows))
        ax.bar(x, tabaco_rows["rmsd"].values, color=colors, edgecolor="black", linewidth=0.5)
        for i, (_, r) in enumerate(tabaco_rows.iterrows()):
            ax.text(i, r["rmsd"]+0.03, f"{r['rmsd']:.2f}", ha="center", va="bottom", fontsize=7, fontweight="bold")
        short = [r["gene"].split("|")[0][:16] for _, r in tabaco_rows.iterrows()]
        ax.set_xticks(x)
        ax.set_xticklabels(short, rotation=45, ha="right", fontsize=7)
        ax.set_ylabel("RMSD (\u00c5) vs Tabaco", fontsize=11)
        ax.set_title("Gene RMSD vs Tabaco Reference", fontsize=12, fontweight="bold")
        ax.axhline(y=2.0, color="red", linestyle="--", alpha=0.4, label="2 \u00c5")
        ax.legend(fontsize=8)
        plt.tight_layout()
        fig.savefig(FIGURES / "phylo_rmsd_vs_tabaco.png", dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(f"  Saved phylo_rmsd_vs_tabaco.png")

    # --- Structural deviation scatter ---
    df = results_df.dropna(subset=["tm_score", "rmsd"]).copy()
    if not df.empty:
        fig, ax = plt.subplots(figsize=(9, 7))
        colors = sns.color_palette("husl", n_colors=len(df))
        for i, (_, row) in enumerate(df.iterrows()):
            s1 = row["struct1"].split("|")[0][:10]
            s2 = row["struct2"].split("|")[0][:10]
            ax.scatter(row["rmsd"], row["tm_score"], c=[colors[i]], s=80, alpha=0.8,
                       edgecolors="black", linewidth=0.5)
            ax.annotate(f"{s1}-{s2}", (row["rmsd"], row["tm_score"]),
                        fontsize=5.5, alpha=0.7, xytext=(3, 3), textcoords="offset points")
        ax.axhline(y=0.5, color="gray", linestyle="--", alpha=0.5, label="TM-score = 0.5")
        ax.axvline(x=2.0, color="gray", linestyle=":", alpha=0.5, label="RMSD = 2 \u00c5")
        ax.set_xlabel("RMSD (\u00c5)", fontsize=11)
        ax.set_ylabel("TM-score", fontsize=11)
        ax.set_title("Structural Deviation: TM-score vs RMSD", fontsize=12, fontweight="bold")
        ax.legend(fontsize=8)
        plt.tight_layout()
        fig.savefig(FIGURES / "phylo_structural_deviation.png", dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(f"  Saved phylo_structural_deviation.png")

    # --- %identity heatmap ---
    gs = gene_seqs if len(gene_seqs) >= 2 else all_seqs
    gl = list(gs.keys())
    n = len(gl)
    if n >= 2:
        im = np.zeros((n, n))
        for i in range(n):
            for j in range(n):
                im[i, j] = pid(gs[gl[i]], gs[gl[j]])
        short = [l.split("|")[0][:18] for l in gl]
        fig, ax = plt.subplots(figsize=(max(10, n*0.5), max(8, n*0.45)))
        sns.heatmap(im, annot=True, fmt=".0f", cmap="viridis",
                    xticklabels=short, yticklabels=short,
                    vmin=0, vmax=100, square=True, linewidths=0.3,
                    cbar_kws={"label": "% Identity"}, ax=ax,
                    annot_kws={"fontsize": 5})
        ax.set_title("Stricto Gene Sequence Identity", fontsize=12, fontweight="bold")
        plt.xticks(rotation=45, ha="right", fontsize=6)
        plt.yticks(fontsize=6)
        plt.tight_layout()
        fig.savefig(FIGURES / "phylo_identity_heatmap.png", dpi=300, bbox_inches="tight")
        plt.close(fig)
        pd.DataFrame(im, index=gl, columns=gl).to_csv(
            PAPER / "tables" / "stricto_phylo_identity.csv")
        print(f"  Saved phylo_identity_heatmap.png")

    # --- TM-score vs identity ---
    identities, tmscores, plabels = [], [], []
    for _, row in df.iterrows():
        k1, k2 = row["struct1"].split("|")[0], row["struct2"].split("|")[0]
        if k1 in all_seqs and k2 in all_seqs:
            p = pid(all_seqs[k1], all_seqs[k2])
            identities.append(p)
            tmscores.append(row["tm_score"])
            plabels.append(f"{k1[:10]}-{k2[:10]}")
    if identities:
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.scatter(identities, tmscores, c=tmscores, cmap="plasma", s=80,
                   edgecolors="black", linewidth=0.5, alpha=0.8)
        for i, lb in enumerate(plabels):
            ax.annotate(lb, (identities[i], tmscores[i]),
                        fontsize=5, alpha=0.6, xytext=(3, 3), textcoords="offset points")
        if len(identities) > 2:
            z = np.polyfit(identities, tmscores, 1)
            x_line = np.linspace(min(identities), max(identities), 100)
            ax.plot(x_line, np.poly1d(z)(x_line), "r--", alpha=0.6,
                    label=f"Trend (r\u00b2={np.corrcoef(identities, tmscores)[0,1]**2:.2f})")
        ax.set_xlabel("Sequence Identity (%)", fontsize=11)
        ax.set_ylabel("TM-score", fontsize=11)
        ax.set_title("TM-score vs Sequence Identity", fontsize=12, fontweight="bold")
        ax.legend(fontsize=8)
        plt.tight_layout()
        fig.savefig(FIGURES / "phylo_tmscore_vs_identity.png", dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(f"  Saved phylo_tmscore_vs_identity.png")

    print(f"\n{'='*70}")
    print(f"Analysis complete!")
    print(f"Figures: {FIGURES}/phylo_*.png")
    print(f"{'='*70}")
    for f in sorted(FIGURES.glob("phylo_*.png")):
        print(f"  {f.name} ({os.path.getsize(f)/1024:.0f} KB)")

    print(f"\nTo improve with real Boltz predictions (requires GPU):")
    print(f"  bash {batch_script}")
    print(f"  Then re-run this script (it will find the new PDBs)")


if __name__ == "__main__":
    main()
