import os
import sys
import warnings
import subprocess
import shutil
from pathlib import Path
from itertools import combinations

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap

from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Align
from Bio.Align import substitution_matrices
from Bio.PDB import PDBParser, PDBIO, MMCIFParser
from Bio.PDB.Polypeptide import protein_letters_3to1 as _3to1_dict
_3TO1 = _3to1_dict.copy()
_3TO1['MSE'] = 'M'
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
STRICTO_PATH = PROTEINS / STRICTO
STRICTO_REF = STRICTO_PATH / "reference"
STRICTO_ANALYSIS = STRICTO_PATH / "analysis"
STRICTO_STRUCTURAL = STRICTO_ANALYSIS / "structural"
STRICTO_PHYLO = STRICTO_ANALYSIS / "phylo"

TABACO_UNIPROT = "A0A1U7VHK2"
ARABIDOPSIS_UNIPROT = "P94111"
CATHARANTHUS_UNIPROT = "P18417"

FIGURES.mkdir(parents=True, exist_ok=True)

sns.set_style("whitegrid")
sns.set_context("paper", font_scale=1.2)


def extract_coords_and_seq(pdb_path):
    """Extract CA coordinates and amino acid sequence from a structure file."""
    parser = PDBParser(QUIET=True)
    try:
        structure = parser.get_structure("prot", str(pdb_path))
    except Exception:
        try:
            cif_parser = MMCIFParser(QUIET=True)
            structure = cif_parser.get_structure("prot", str(pdb_path))
        except Exception as e:
            print(f"  [FAIL] cannot parse {pdb_path}: {e}")
            return None, None

    coords = []
    seq = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.id[0] == " ":
                    if "CA" in residue:
                        coords.append(residue["CA"].coord)
                    if hasattr(residue, "resname"):
                        seq.append(_3to1(residue.resname))
    if len(coords) < 10:
        return None, None
    return np.array(coords, dtype=float), "".join(seq)


def collect_stricto_structures():
    """Find all stricto-related PDB/CIF files."""
    structs = {}

    refs = {
        "arabidopsis_ref": STRICTO_REF / "arabidopsis" / "structure" / f"stricto_{ARABIDOPSIS_UNIPROT}.pdb",
        "catharanthus_ref": STRICTO_REF / "catharanthus" / "structure" / f"{CATHARANTHUS_UNIPROT}.pdb",
        "catharanthus_xray": STRICTO_REF / "catharanthus" / "structure" / "x_ray_6ZEA.pdb",
    }
    for label, path in refs.items():
        if path.exists():
            structs[label] = str(path)

    preds = {
        "tucunaca_alphafold": STRICTO_STRUCTURAL / "alphafold" / "tucunaca" / "best_model.pdb",
        "tucunaca_p2rank": STRICTO_STRUCTURAL / "p2rank" / "tucunaca" / "output" / "visualizations" / "data" / "best_model.pdb",
        "arabidopsis_p2rank": STRICTO_STRUCTURAL / "p2rank" / "arabidopsis" / "output" / "visualizations" / "data" / f"stricto_{ARABIDOPSIS_UNIPROT}.pdb",
    }
    for label, path in preds.items():
        if path.exists():
            structs[label] = str(path)

    boltz_dir = STRICTO_STRUCTURAL / "boltz"
    boltz_map = {
        "arabidopsis_boltz": boltz_dir / "arabidopsis" / "r1" / "predictions" / "arabidopsis_stricto_P94111_r1" / "arabidopsis_stricto_P94111_r1_model_0.cif",
        "catharanthus_boltz": boltz_dir / "catharanthus" / "r1" / "predictions" / "catharanthus_stricto_P18417_r1" / "catharanthus_stricto_P18417_r1_model_0.cif",
        "tabaco_boltz": boltz_dir / "tabaco" / "r1" / "predictions" / "tabaco_stricto_A0A1U7VHK2_r1" / "tabaco_stricto_A0A1U7VHK2_r1_model_0.cif",
        "tucunaca_boltz_arab": boltz_dir / "tucunaca" / "arabidopsis" / "r1" / "predictions" / "tucunaca_arabidopsis_stricto_P94111_r1" / "tucunaca_arabidopsis_stricto_P94111_r1_model_0.cif",
        "caupuri_boltz_catha": boltz_dir / "caupuri" / "catharanthus" / "rr" / "predictions" / "caupuri_catharanthus_stricto_P18417_rr" / "caupuri_catharanthus_stricto_P18417_rr_model_0.cif",
        "caupuri_boltz_tabaco": boltz_dir / "caupuri" / "tabaco" / "rr" / "predictions" / "caupuri_tabaco_stricto_A0A1U7VHK2_rr" / "caupuri_tabaco_stricto_A0A1U7VHK2_rr_model_0.cif",
    }
    for label, path in boltz_map.items():
        if path.exists():
            structs[label] = str(path)

    return structs


def read_sequences_from_phylo():
    """Read all sequences from phylo alignment fasta files."""
    seqs = {}
    for phylo_dir in [STRICTO_PHYLO / "caupuri", STRICTO_PHYLO / "tucunaca"]:
        fasta = phylo_dir / "clean_alignment.fasta"
        if fasta.exists():
            for rec in SeqIO.parse(str(fasta), "fasta"):
                seqs[rec.id.split(" ")[0]] = str(rec.seq).replace("-", "")
    return seqs


def read_reference_sequences():
    """Read reference sequences from FASTA."""
    seqs = {}
    ref_fastas = {
        "arabidopsis_ref": STRICTO_REF / "arabidopsis" / "structure" / f"stricto_{ARABIDOPSIS_UNIPROT}.fasta",
        "catharanthus_ref": STRICTO_REF / "catharanthus" / "structure" / f"stricto_{CATHARANTHUS_UNIPROT}.fasta",
        "tabaco_ref": STRICTO_REF / "tabaco" / "structure" / f"stricto_{TABACO_UNIPROT}.fasta",
    }
    for label, path in ref_fastas.items():
        if path.exists():
            for rec in SeqIO.parse(str(path), "fasta"):
                seqs[label] = str(rec.seq)
    return seqs


def compute_percent_identity(seq1, seq2):
    """Compute % identity between two sequences using pairwise alignment."""
    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    try:
        alns = aligner.align(seq1, seq2)
        aln = alns[0]
        pairs = list(zip(str(aln[0]), str(aln[1])))
        identical = sum(1 for a, b in pairs if a == b and a != "-")
        aligned = sum(1 for a, b in pairs if a != "-" or b != "-")
        return identical / aligned * 100 if aligned > 0 else 0.0
    except Exception:
        return 0.0


def run_tmalign(coords1, seq1, coords2, seq2):
    """Run TM-align on two structures."""
    if coords1 is None or coords2 is None:
        return None, None
    if len(coords1) < 10 or len(coords2) < 10:
        return None, None
    try:
        result = tm_align(coords1, coords2, seq1, seq2)
        return result.tm_norm_chain1, result.rmsd
    except Exception as e:
        print(f"  [TM-align error] {e}")
        return None, None


def generate_structure_fastas(seqs_dict, output_dir):
    """Write FASTA files for all sequences for structure prediction."""
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    combined = out / "stricto_all_seqs.fasta"
    with open(combined, "w") as f:
        for label, seq in seqs_dict.items():
            clean = label.replace("|", "_").replace(" ", "_").split("_td")[0].split("_sp")[0]
            f.write(f">{clean}\n{seq}\n")
    print(f"\n  FASTA for structure prediction: {combined}")
    print(f"  Use: colabfold_batch {combined} {out / 'predictions'} --num-recycle 3")

    for label, seq in seqs_dict.items():
        clean = label.replace("|", "_").replace(" ", "_").split("_td")[0].split("_sp")[0]
        fasta_path = out / f"{clean[:50]}.fasta"
        if not fasta_path.exists():
            with open(fasta_path, "w") as f:
                f.write(f">{clean}\n{seq}\n")

    colabfold = shutil.which("colabfold_batch")
    if colabfold:
        print("  ColabFold detected - will run predictions if needed")

    return str(combined)


def incorporate_predicted_structures(structs, seqs_dict, fasta_dir):
    """Look for predicted structures and add them to the structure dict."""
    import glob as glb

    patterns = [
        str(Path(fasta_dir) / "predictions" / "*" / "*_model_0.cif"),
        str(Path(fasta_dir) / "predictions" / "*" / "*_model_0.pdb"),
        str(Path(fasta_dir) / "*" / "best_model.pdb"),
        str(Path(fasta_dir) / "*" / "*_model_0.pdb"),
    ]
    for pattern in patterns:
        for f in glb.glob(pattern):
            name = Path(f).parent.name
            if name not in structs:
                structs[f"pred_{name}"] = f
    return structs


def make_tmscore_heatmap(results_df, title, filename):
    """Create TM-score heatmap (distance matrix)."""
    if results_df.empty:
        print(f"  [SKIP] No data for {title}")
        return

    labels = results_df.index.tolist() if hasattr(results_df, 'index') else results_df.columns.tolist()

    pivot = results_df.pivot_table(index="structure1", columns="structure2",
                                   values="tm_score", aggfunc="first")
    all_labels = list(pivot.index) + [c for c in pivot.columns if c not in pivot.index]
    all_labels = sorted(set(all_labels))

    full = pd.DataFrame(index=all_labels, columns=all_labels, dtype=float)
    for i in all_labels:
        full.loc[i, i] = 1.0
        for j in all_labels:
            if i != j:
                if i in pivot.index and j in pivot.columns and pd.notna(pivot.loc[i, j]):
                    full.loc[i, j] = pivot.loc[i, j]
                elif j in pivot.index and i in pivot.columns and pd.notna(pivot.loc[j, i]):
                    full.loc[i, j] = pivot.loc[j, i]

    short_labels = [l.split("_")[0] if "_" in l else l[:12] for l in all_labels]

    fig, ax = plt.subplots(figsize=(12, 10))
    mask = np.isnan(full.values) | np.eye(len(full), dtype=bool)
    cmap = sns.color_palette("YlOrRd", as_cmap=True)

    sns.heatmap(full, annot=True, fmt=".3f", cmap=cmap,
                xticklabels=short_labels, yticklabels=short_labels,
                vmin=0, vmax=1, square=True,
                linewidths=0.5, cbar_kws={"label": "TM-score"},
                ax=ax, mask=mask, annot_kws={"fontsize": 7})

    ax.set_title(title, fontsize=14, fontweight="bold")
    ax.set_xlabel("Structure", fontsize=11)
    ax.set_ylabel("Structure", fontsize=11)
    plt.xticks(rotation=45, ha="right", fontsize=8)
    plt.yticks(fontsize=8)
    plt.tight_layout()
    fig.savefig(FIGURES / filename, dpi=300, bbox_inches="tight")
    print(f"  Saved {FIGURES / filename}")
    plt.close(fig)


def make_rmsd_barplot(results_df, title, filename):
    """Create RMSD bar plot focusing on tucunaca/caupuri vs tabaco."""
    tabaco_rows = results_df[
        results_df["structure1"].str.contains("tabaco", case=False) |
        results_df["structure2"].str.contains("tabaco", case=False)
    ].copy()

    if tabaco_rows.empty:
        print(f"  [SKIP] No tabaco comparison data")
        return

    def get_other(row):
        t1 = "tabaco" in row["structure1"].lower()
        t2 = "tabaco" in row["structure2"].lower()
        return row["structure2"] if t1 else row["structure1"]

    tabaco_rows["other"] = tabaco_rows.apply(get_other, axis=1)
    tabaco_rows = tabaco_rows.dropna(subset=["rmsd"])

    ooi_order = ["tucunaca", "caupuri", "arabidopsis", "catharanthus"]
    tabaco_rows["sort_order"] = tabaco_rows["other"].apply(
        lambda x: next((i for i, o in enumerate(ooi_order) if o in x.lower()), 99)
    )
    tabaco_rows = tabaco_rows.sort_values("sort_order")

    fig, ax = plt.subplots(figsize=(10, 6))
    colors = sns.color_palette("Set2", n_colors=len(tabaco_rows))
    bars = ax.bar(range(len(tabaco_rows)), tabaco_rows["rmsd"].values,
                  color=colors, edgecolor="black", linewidth=0.5)

    for i, (_, row) in enumerate(tabaco_rows.iterrows()):
        ax.text(i, row["rmsd"] + 0.02, f"{row['rmsd']:.2f}",
                ha="center", va="bottom", fontsize=9, fontweight="bold")

    ax.set_xticks(range(len(tabaco_rows)))
    short_labels = []
    for _, row in tabaco_rows.iterrows():
        label = row["other"].replace("_", " ").replace("boltz", "").replace("alphafold", "AF").replace("ref", "").replace("p2rank", "p2r").strip()
        if len(label) > 15:
            label = label[:15] + "..."
        short_labels.append(label)
    ax.set_xticklabels(short_labels, rotation=45, ha="right", fontsize=9)

    ax.set_ylabel("RMSD (\u00c5)", fontsize=12)
    ax.set_title(title, fontsize=13, fontweight="bold")
    ax.axhline(y=2.0, color="red", linestyle="--", alpha=0.5, label="2 \u00c5 threshold")
    ax.legend(fontsize=9)
    plt.tight_layout()
    fig.savefig(FIGURES / filename, dpi=300, bbox_inches="tight")
    print(f"  Saved {FIGURES / filename}")
    plt.close(fig)

    return tabaco_rows


def make_deviation_scatter(results_df, title, filename):
    """Create TM-score vs RMSD scatter plot."""
    df = results_df.dropna(subset=["tm_score", "rmsd"]).copy()
    if df.empty:
        print(f"  [SKIP] No data for {title}")
        return

    fig, ax = plt.subplots(figsize=(9, 7))

    colors = sns.color_palette("husl", n_colors=len(df))
    for i, (_, row) in enumerate(df.iterrows()):
        label = f"{row['structure1'].split('_')[0]} vs {row['structure2'].split('_')[0]}"
        ax.scatter(row["rmsd"], row["tm_score"], c=[colors[i]], s=100, alpha=0.8,
                   edgecolors="black", linewidth=0.5)
        ax.annotate(label, (row["rmsd"], row["tm_score"]),
                    fontsize=6, alpha=0.8, xytext=(5, 5),
                    textcoords="offset points")

    ax.axhline(y=0.5, color="gray", linestyle="--", alpha=0.5, label="TM-score = 0.5 (same fold)")
    ax.axvline(x=2.0, color="gray", linestyle=":", alpha=0.5, label="RMSD = 2 \u00c5")

    ax.set_xlabel("RMSD (\u00c5)", fontsize=12)
    ax.set_ylabel("TM-score", fontsize=12)
    ax.set_title(title, fontsize=13, fontweight="bold")
    ax.legend(fontsize=9)
    plt.tight_layout()
    fig.savefig(FIGURES / filename, dpi=300, bbox_inches="tight")
    print(f"  Saved {FIGURES / filename}")
    plt.close(fig)


def make_identity_heatmap(seqs_dict, title, filename):
    """Create %identity heatmap from sequences."""
    labels = list(seqs_dict.keys())
    short = [l.split("|")[0].split("_td")[0].split("_sp")[0] for l in labels]
    short = [l.split("_stricto")[0] if "_stricto" in l else l for l in short]
    short = [l[:20] for l in short]

    n = len(labels)
    identity = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            identity[i, j] = compute_percent_identity(seqs_dict[labels[i]], seqs_dict[labels[j]])

    fig, ax = plt.subplots(figsize=(max(8, n * 0.5), max(7, n * 0.45)))

    cmap = sns.color_palette("viridis", as_cmap=True)
    sns.heatmap(identity, annot=True, fmt=".0f", cmap=cmap,
                xticklabels=short, yticklabels=short,
                vmin=0, vmax=100, square=True,
                linewidths=0.3, cbar_kws={"label": "% Identity"},
                ax=ax, annot_kws={"fontsize": 5})

    ax.set_title(title, fontsize=13, fontweight="bold")
    plt.xticks(rotation=45, ha="right", fontsize=6)
    plt.yticks(fontsize=6)
    plt.tight_layout()
    fig.savefig(FIGURES / filename, dpi=300, bbox_inches="tight")
    print(f"  Saved {FIGURES / filename}")
    plt.close(fig)

    return identity


def make_full_tmscore_vs_identity(results_df, seqs_dict, filename):
    """Scatter TM-score vs sequence identity for all pairs."""
    df = results_df.dropna(subset=["tm_score", "rmsd"]).copy()
    if df.empty:
        return

    identities = []
    tmscores = []
    labels = []
    for _, row in df.iterrows():
        s1, s2 = row["structure1"], row["structure2"]
        k1 = next((k for k in seqs_dict if k.startswith(s1.split("_")[0])), None)
        k2 = next((k for k in seqs_dict if k.startswith(s2.split("_")[0])), None)
        if k1 and k2:
            pid = compute_percent_identity(seqs_dict[k1], seqs_dict[k2])
            identities.append(pid)
            tmscores.append(row["tm_score"])
            labels.append(f"{s1.split('_')[0]}-{s2.split('_')[0]}")

    if not identities:
        return

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(identities, tmscores, c=tmscores, cmap="plasma", s=80,
               edgecolors="black", linewidth=0.5, alpha=0.8)
    for i, label in enumerate(labels):
        ax.annotate(label, (identities[i], tmscores[i]),
                    fontsize=6, alpha=0.7, xytext=(3, 3),
                    textcoords="offset points")

    z = np.polyfit(identities, tmscores, 1)
    p = np.poly1d(z)
    x_line = np.linspace(min(identities), max(identities), 100)
    ax.plot(x_line, p(x_line), "r--", alpha=0.6, label=f"Trend (r\u00b2={np.corrcoef(identities, tmscores)[0,1]**2:.2f})")

    ax.set_xlabel("Sequence Identity (%)", fontsize=12)
    ax.set_ylabel("TM-score", fontsize=12)
    ax.set_title("TM-score vs Sequence Identity", fontsize=13, fontweight="bold")
    ax.legend(fontsize=9)
    plt.tight_layout()
    fig.savefig(FIGURES / filename, dpi=300, bbox_inches="tight")
    print(f"  Saved {FIGURES / filename}")
    plt.close(fig)


def main():
    print("=" * 70)
    print("STRICTO TM-align Analysis Pipeline")
    print("=" * 70)

    print("\n[1/6] Collecting stricto structures...")
    all_structs = collect_stricto_structures()
    print(f"  Found {len(all_structs)} structure files:")
    for label, path in sorted(all_structs.items()):
        size = os.path.getsize(path) / 1024
        print(f"    {label:30s} -> {Path(path).name} ({size:.0f} KB)")

    print("\n[2/6] Extracting coordinates and sequences from structures...")
    struct_data = {}
    for label, path in all_structs.items():
        coords, seq = extract_coords_and_seq(path)
        if coords is not None and len(coords) >= 10:
            struct_data[label] = {"coords": coords, "seq": seq, "path": path}
            print(f"    {label:30s} {len(coords):4d} residues")
        else:
            print(f"    {label:30s} [FAILED]")

    if len(struct_data) < 2:
        print("  [ERROR] Need at least 2 structures to compare")
        sys.exit(1)

    print("\n[3/6] Reading sequences from phylo analysis...")
    phylo_seqs = read_sequences_from_phylo()
    ref_seqs = read_reference_sequences()
    all_seqs = {**ref_seqs, **phylo_seqs}
    print(f"  Found {len(all_seqs)} sequences total")

    print("\n[3b/6] Generating FASTA files for structure prediction...")
    fasta_dir = PROTEINS / STRICTO / "analysis" / "structural" / "predictions"
    generate_structure_fastas(all_seqs, fasta_dir)

    all_structs = incorporate_predicted_structures(all_structs, all_seqs, fasta_dir)

    print("\n[4/6] Running pairwise TM-align on all structures...")
    struct_labels = list(struct_data.keys())
    results = []
    for s1, s2 in combinations(struct_labels, 2):
        d1 = struct_data[s1]
        d2 = struct_data[s2]
        tm, rmsd = run_tmalign(d1["coords"], d1["seq"], d2["coords"], d2["seq"])
        if tm is not None:
            results.append({
                "structure1": s1, "structure2": s2,
                "tm_score": round(tm, 5), "rmsd": round(rmsd, 4)
            })
            print(f"    {s1:30s} vs {s2:30s}  TM={tm:.4f}  RMSD={rmsd:.2f}\u00c5")
        else:
            print(f"    {s1:30s} vs {s2:30s}  [FAILED]")

    results_df = pd.DataFrame(results)
    if results_df.empty:
        print("  [ERROR] No TM-align results")
        sys.exit(1)

    results_csv = PAPER / "tables" / "stricto_tmalign_results.csv"
    results_csv.parent.mkdir(parents=True, exist_ok=True)
    results_df.to_csv(results_csv, index=False)
    print(f"\n  Results saved to {results_csv}")

    print("\n[5/6] Computing %identity from phylo sequences...")
    identity_matrix = make_identity_heatmap(
        all_seqs,
        "Stricto Sequence Identity Matrix",
        "stricto_identity_heatmap.png"
    )

    identity_csv = PAPER / "tables" / "stricto_identity_matrix.csv"
    identity_df = pd.DataFrame(identity_matrix,
                                index=list(all_seqs.keys()),
                                columns=list(all_seqs.keys()))
    identity_df.to_csv(identity_csv)
    print(f"  Identity matrix saved to {identity_csv}")

    print("\n[6/6] Generating plots...")

    make_tmscore_heatmap(
        results_df,
        "Strictosidine Synthase - TM-score Matrix",
        "stricto_tmscore_heatmap.png"
    )

    make_rmsd_barplot(
        results_df,
        "RMSD vs Tabaco Reference (Strictosidine Synthase)",
        "stricto_rmsd_vs_tabaco.png"
    )

    make_deviation_scatter(
        results_df,
        "Structural Deviation: TM-score vs RMSD",
        "stricto_structural_deviation.png"
    )

    make_full_tmscore_vs_identity(
        results_df, all_seqs,
        "stricto_tmscore_vs_identity.png"
    )

    print("\n" + "=" * 70)
    print("All analyses complete!")
    print(f"Figures saved to: {FIGURES}")
    print("=" * 70)
    print("\nGenerated files:")
    for f in sorted(FIGURES.glob("stricto_*.png")):
        size = os.path.getsize(f) / 1024
        print(f"  {f.name:45s} ({size:.0f} KB)")


if __name__ == "__main__":
    main()
