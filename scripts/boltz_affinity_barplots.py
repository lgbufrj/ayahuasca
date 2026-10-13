
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from mpl_toolkits.mplot3d import Axes3D


from data import PAPER_PATH

# ============================================
# CONFIGURATION
# ============================================
# Output files and figures paths
CSV = Path(
    f"{PAPER_PATH}/tables/boltz_affinity_results.csv"
)
OUTPUT_DIR = Path(
    f"{PAPER_PATH}/figures/supplementary/singles"
)

# All available variables (pick 1 to 3 of these):
VARIABLES = [
    "binding_affinity",
    "kd",
    "dG",
    "IC50",
    "delta_affinity",
    "delta_kd",
    "delta_dG",
    "delta_IC50",
    "confidence",
    "delta_confidence",
    "plddt",
    "delta_plddt",
    "iplddt",
    "delta_iplddt",
    "pde",
    "delta_pde",
    "ipde",
    "delta_ipde",
    "ligand_distance"
]

# PLOTS: each entry is a tuple of 1-3 variables
# 1 variable -> barplot, 2 variables -> scatterplot, 3 -> 3D scatter
PLOTS = [
    ("ipde", "iplddt"),
    ("dG",),
    ("confidence","dG"),
    ("dG", "ligand_distance"),
]

# All available organisms:
ORGANISMS = ["ecoli", "pseudomonas", "arabidopsis", "catharanthus", "rice", "tabaco", "tucunaca", "yurema"]

# SELECTED ORGANISMS (one or more)
SELECTED_ORGANISMS = [
    "ecoli",
    "pseudomonas"
]

# Colors used for r1 and rr when only one organism is selected
RXN_COLORS = {"r1": "#4C72B0", "rr": "#DD8452", "r2": "#55A868", "r3": "#C44E52", "r4": "#8172B3", "r5": "#937860"}

# Up to 6 different colors used for the organism + reaction key
# when two or more organisms are selected
COMBO_COLORS = [
    "#4C72B0",
    "#DD8452",
    "#55A868",
    "#C44E52",
    "#8172B3",
    "#937860",
    "#777777",
    "#D9C426",
    "#8172B3",
    "#69E72F",
    "#777777",
]

# ============================================
# DATA LOADING
# ============================================
df_all = pd.read_csv(CSV, skipinitialspace=True)
df_all.columns = [c.strip() for c in df_all.columns]
for col in df_all.select_dtypes(include="object").columns:
    df_all[col] = df_all[col].str.strip()
df_all = df_all.sort_values(
    ["organism", "id", "rxn_id"], kind="stable"
)

# ============================================
# HELPERS
# ============================================
def select_best_models(df):
    """Select best model per (protein, organism, id, rxn_id).
    Best = lowest dG (more negative is better), tie-breaker lowest ligand_distance.
    """
    if df.empty:
        return df
    df_sel = df.copy()
    # Ensure numeric
    if "dG" in df_sel.columns:
        df_sel["dG"] = pd.to_numeric(df_sel["dG"], errors="coerce")
    if "ligand_distance" in df_sel.columns:
        df_sel["ligand_distance"] = pd.to_numeric(df_sel["ligand_distance"], errors="coerce")

    # Drop rows missing key criteria if present
    has_criteria = df_sel[["dG", "ligand_distance"]].notna().all(axis=1)
    df_work = df_sel[has_criteria].copy() if has_criteria.any() else df_sel.copy()

    # Sort by dG ascending (most negative first), then ligand_distance ascending
    df_work = df_work.sort_values(["dG", "ligand_distance"], ascending=[True, True])

    # Keep first (best) per group
    best = df_work.groupby(["protein", "organism", "id", "rxn_id"], as_index=False).first()

    # Preserve original column order if possible
    cols = [c for c in df_sel.columns if c in best.columns]
    return best[cols] if cols else best


def build_combos(df, orgs):
    # Builds the (organism, reaction) pairs present in the data
    combos = []
    for org in orgs:
        for rxn in ("r1", "rr", "r2", "r3", "r4", "r5"):
            if ((df["organism"] == org) & (df["rxn_id"] == rxn)).any():
                combos.append((org, rxn))
    return combos


def row_colors(df, orgs):
    # One color per row, based on the (organism, reaction) combo
    combos = build_combos(df, orgs)
    combo_color = {
        combo: COMBO_COLORS[i % len(COMBO_COLORS)] for i, combo in enumerate(combos)
    }
    return [combo_color[(row["organism"], row["rxn_id"])] for _, row in df.iterrows()]


def legend_handles(df, orgs):
    # Legend key: r1/rr colors for a single organism,
    # organism + reaction colors (up to 6) for two or more organisms
    if len(orgs) == 1:
        return [
            plt.Rectangle((0, 0), 1, 1, facecolor=c, edgecolor="black", label=rxn)
            for rxn, c in RXN_COLORS.items()
        ]
    combos = build_combos(df, orgs)
    return [
        plt.Rectangle(
            (0, 0),
            1,
            1,
            facecolor=COMBO_COLORS[i % len(COMBO_COLORS)],
            edgecolor="black",
            label=f"{org} {rxn}",
        )
        for i, (org, rxn) in enumerate(combos)
    ]


# ============================================
# PLOTTING
# ============================================
def bar_plot(df, protein, orgs, variables):
    # Barplot used when exactly 1 variable is selected
    fig, ax = plt.subplots(figsize=(12, 5))

    labels = [
        f"{row['id']}\n{row['rxn_id']}"
        for _, row in df[["id", "rxn_id"]].iterrows()
    ]
    bar_colors = [RXN_COLORS[row["rxn_id"]] for _, row in df.iterrows()]

    ax.bar(range(len(df)), df[variables[0]], color=bar_colors, edgecolor="black", linewidth=0.3)

    # Per-reaction max/min y-lines
    for rxn in sorted(df["rxn_id"].unique()):
        rxn_vals = df.loc[df["rxn_id"] == rxn, variables[0]]
        max_val = rxn_vals.max()
        min_val = rxn_vals.min()
        ax.axhline(max_val, color=RXN_COLORS[rxn], linestyle="--", linewidth=0.8, alpha=0.7)
        ax.text(len(df) - 0.3, max_val, f" max {rxn}", fontsize=5, va="center", ha="left",
                color=RXN_COLORS[rxn], fontweight="bold")
        ax.axhline(min_val, color=RXN_COLORS[rxn], linestyle=":", linewidth=0.8, alpha=0.7)
        ax.text(len(df) - 0.3, min_val, f" min {rxn}", fontsize=5, va="center", ha="left",
                color=RXN_COLORS[rxn], fontweight="bold")

    # Protein group visual clustering (alternating shading + dividers)
    id_ranges = []
    prev_id = None
    start = 0
    for i, (_, row) in enumerate(df.iterrows()):
        if row["id"] != prev_id:
            if prev_id is not None:
                id_ranges.append((prev_id, start, i - 1))
            start = i
            prev_id = row["id"]
    if prev_id is not None:
        id_ranges.append((prev_id, start, len(df) - 1))

    for idx, (uid, s, e) in enumerate(id_ranges):
        if idx % 2 == 1:
            ax.axvspan(s - 0.5, e + 0.5, ymin=-0.15, ymax=1.0, alpha=0.07, color="gray", clip_on=False)
        if idx > 0:
            ax.axvline(s - 0.5, color="gray", linestyle="-", linewidth=0.8, alpha=0.5, ymin=-0.15, ymax=1.0, clip_on=False)

    # Zero reference line for delta variables
    if variables[0].startswith("delta"):
        ax.axhline(0.0, color="grey", linestyle="--", linewidth=0.8)

    ax.set_xticks(range(len(df)))
    ax.set_xticklabels(labels, fontsize=6, rotation=90)
    ax.set_ylabel(variables[0], fontsize=9)
    ax.margins(x=0.01)

    # Plot legend/key (reaction colors only)
    fig.legend(
        handles=[
            plt.Rectangle((0, 0), 1, 1, facecolor=RXN_COLORS[rxn], edgecolor="black", label=rxn)
            for rxn in sorted(df["rxn_id"].unique())
        ],
        loc="upper right",
        ncol=1,
        fontsize=8,
        frameon=False,
    )
    title = f"{protein} | {variables[0]} - organisms: {', '.join(orgs)}"
    fig.suptitle(title, fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.93))

    out = OUTPUT_DIR / f"boltz_barplot_{protein}_{'_'.join(variables)}.svg"
    fig.set_size_inches(len(df)/3, 6)
    fig.savefig(out, dpi=300)
    plt.close(fig)
    print(f"saved {out}")


def scatter_plot(df, protein, orgs, variables):
    # Scatterplot used when exactly 2 variables are selected
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(df[variables[0]], df[variables[1]], c=row_colors(df, orgs), s=40)

    # Label most extreme (x, y) pairs per reaction
    for rxn in sorted(df["rxn_id"].unique()):
        rxn_data = df[df["rxn_id"] == rxn]
        if len(rxn_data) < 2:
            continue
        x_min, x_max = rxn_data[variables[0]].min(), rxn_data[variables[0]].max()
        y_min, y_max = rxn_data[variables[1]].min(), rxn_data[variables[1]].max()
        x_range = x_max - x_min
        y_range = y_max - y_min
        if x_range > 0 and y_range > 0:
            combined = (rxn_data[variables[0]] - x_min) / x_range + (rxn_data[variables[1]] - y_min) / y_range
        elif x_range > 0:
            combined = (rxn_data[variables[0]] - x_min) / x_range
        elif y_range > 0:
            combined = (rxn_data[variables[1]] - y_min) / y_range
        else:
            continue
        for label, idx in [("max", combined.idxmax()), ("min", combined.idxmin())]:
            row = df.loc[idx]
            ax.annotate(
                f"{row['id']}",
                xy=(row[variables[0]], row[variables[1]]),
                xytext=(5, 5 if label == "max" else -10), textcoords="offset points",
                fontsize=5, color=RXN_COLORS[rxn], fontweight="bold",
            )

    ax.set_xlabel(variables[0], fontsize=10)
    ax.set_ylabel(variables[1], fontsize=10)
    # Plot legend/key
    fig.legend(
        handles=legend_handles(df, orgs),
        loc="upper right",
        ncol=1,
        fontsize=8,
        frameon=False,
    )
    fig.suptitle(
        f"{protein} | {variables[0]} vs {variables[1]} - organisms: {', '.join(orgs)}",
        fontsize=12,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.93))

    out = OUTPUT_DIR / f"boltz_scatter_{protein}_{'_'.join(variables)}.svg"
    fig.savefig(out, dpi=300)
    plt.close(fig)
    print(f"saved {out}")


def plot_3d(df, protein, orgs, variables):
    # 3D scatterplot used when exactly 3 variables are selected
    fig = plt.figure(figsize=(9, 7))
    ax = fig.add_subplot(111, projection="3d")
    ax.scatter(
        df[variables[0]],
        df[variables[1]],
        df[variables[2]],
        c=row_colors(df, orgs),
        s=40,
    )

    ax.set_xlabel(variables[0], fontsize=10)
    ax.set_ylabel(variables[1], fontsize=10)
    ax.set_zlabel(variables[2], fontsize=10)
    # Plot legend/key
    handles = legend_handles(df, orgs)
    if handles:
        fig.legend(
            handles=handles,
            loc="upper right",
            ncol=1,
            fontsize=8,
            frameon=False,
        )
    fig.suptitle(
        f"{protein} | {variables[0]} vs {variables[1]} vs {variables[2]} "
        f"- organisms: {', '.join(orgs)}",
        fontsize=12,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.93))

    out = OUTPUT_DIR / f"boltz_3d_{protein}_{'_'.join(variables)}.svg"
    fig.savefig(out, dpi=300, )
    plt.close(fig)
    print(f"saved {out}")


# ============================================
# SUBPLOT FUNCTIONS
# ============================================
def bar_subplot(ax, df, variable):
    bar_colors = [RXN_COLORS[row["rxn_id"]] for _, row in df.iterrows()]
    labels = [f"{row['id']}\n{row['rxn_id']}" for _, row in df[["id", "rxn_id"]].iterrows()]

    ax.bar(range(len(df)), df[variable], color=bar_colors, edgecolor="black", linewidth=0.3)

    for rxn in sorted(df["rxn_id"].unique()):
        rxn_vals = df.loc[df["rxn_id"] == rxn, variable]
        max_val = rxn_vals.max()
        min_val = rxn_vals.min()
        ax.axhline(max_val, color=RXN_COLORS[rxn], linestyle="--", linewidth=0.8, alpha=0.7)
        ax.text(len(df) - 0.3, max_val, f" max {rxn}", fontsize=5, va="center", ha="left",
                color=RXN_COLORS[rxn], fontweight="bold")
        ax.axhline(min_val, color=RXN_COLORS[rxn], linestyle=":", linewidth=0.8, alpha=0.7)
        ax.text(len(df) - 0.3, min_val, f" min {rxn}", fontsize=5, va="center", ha="left",
                color=RXN_COLORS[rxn], fontweight="bold")

    id_ranges = []
    prev_id = None
    start = 0
    for i, (_, row) in enumerate(df.iterrows()):
        if row["id"] != prev_id:
            if prev_id is not None:
                id_ranges.append((prev_id, start, i - 1))
            start = i
            prev_id = row["id"]
    if prev_id is not None:
        id_ranges.append((prev_id, start, len(df) - 1))

    for idx, (uid, s, e) in enumerate(id_ranges):
        if idx % 2 == 1:
            ax.axvspan(s - 0.5, e + 0.5, ymin=-0.15, ymax=1.0, alpha=0.07, color="gray", clip_on=False)
        if idx > 0:
            ax.axvline(s - 0.5, color="gray", linestyle="-", linewidth=0.8, alpha=0.5, ymin=-0.15, ymax=1.0, clip_on=False)

    if variable.startswith("delta"):
        ax.axhline(0.0, color="grey", linestyle="--", linewidth=0.8)

    ax.set_xticks(range(len(df)))
    ax.set_xticklabels(labels, fontsize=5, rotation=90)
    ax.set_ylabel(variable, fontsize=8)
    ax.set_title(variable, fontsize=9, fontweight="bold")
    ax.margins(x=0.01)


def scatter_subplot(ax, df, variables):
    colors = [RXN_COLORS[row["rxn_id"]] for _, row in df.iterrows()]
    ax.scatter(df[variables[0]], df[variables[1]], c=colors, s=30, edgecolors="black", linewidths=0.3)

    # Label most extreme (x, y) pairs per reaction
    for rxn in sorted(df["rxn_id"].unique()):
        rxn_data = df[df["rxn_id"] == rxn]
        if len(rxn_data) < 2:
            continue
        x_min, x_max = rxn_data[variables[0]].min(), rxn_data[variables[0]].max()
        y_min, y_max = rxn_data[variables[1]].min(), rxn_data[variables[1]].max()
        x_range = x_max - x_min
        y_range = y_max - y_min
        if x_range > 0 and y_range > 0:
            combined = (rxn_data[variables[0]] - x_min) / x_range + (rxn_data[variables[1]] - y_min) / y_range
        elif x_range > 0:
            combined = (rxn_data[variables[0]] - x_min) / x_range
        elif y_range > 0:
            combined = (rxn_data[variables[1]] - y_min) / y_range
        else:
            continue
        for label, idx in [("max", combined.idxmax()), ("min", combined.idxmin())]:
            row = df.loc[idx]
            ax.annotate(
                f"{row['id']} ({row['model']})",
                xy=(row[variables[0]], row[variables[1]]),
                xytext=(5, 5 if label == "max" else -10), textcoords="offset points",
                fontsize=5, color=RXN_COLORS[rxn], fontweight="bold",
            )

    ax.set_xlabel(variables[0], fontsize=8)
    ax.set_ylabel(variables[1], fontsize=8)
    ax.set_title(f"{variables[0]} vs {variables[1]}", fontsize=9, fontweight="bold")


# ============================================
# MAIN EXECUTION
# ============================================
# Configuration validation
if not SELECTED_ORGANISMS:
    raise ValueError("SELECTED_ORGANISMS must contain at least one organism.")

unknown_orgs = [o for o in SELECTED_ORGANISMS if o not in ORGANISMS]
if unknown_orgs:
    raise ValueError(f"Unknown organisms in SELECTED_ORGANISMS: {unknown_orgs}")

for entry in PLOTS:
    unknown_vars = [v for v in entry if v not in VARIABLES]
    if unknown_vars:
        raise ValueError(f"Unknown variables in PLOTS entry {entry}: {unknown_vars}")

if __name__ == "__main__":
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    proteins = sorted(df_all["protein"].dropna().unique())
    for protein in proteins:
        df = df_all[
            (df_all["protein"] == protein)
            & df_all["organism"].isin(SELECTED_ORGANISMS)
            & df_all["rxn_id"].isin(["r1", "rr", "r2", "r3", "r4", "r5"])
        ].copy()
        # Select best model per (protein, organism, id, rxn_id)
        df = select_best_models(df)

        valid_plots = []
        for plot_vars in PLOTS:
            df_plot = df.dropna(subset=list(plot_vars))
            if not df_plot.empty:
                valid_plots.append((plot_vars, df_plot))

        if not valid_plots:
            print(f"No data for {protein} with the selected configuration.")
            continue

        n = len(valid_plots)
        ncols = min(n, 2)
        nrows = (n + ncols - 1) // ncols
        max_bars = max(len(df_plot) for _, df_plot in valid_plots)
        fig_width = max(max_bars / 4, 6 * ncols)

        fig, axes = plt.subplots(nrows, ncols, figsize=(fig_width, 5 * nrows))
        if n == 1:
            axes = [axes]
        else:
            axes = list(axes.flatten())

        for i, (plot_vars, df_plot) in enumerate(valid_plots):
            ax = axes[i]
            if len(plot_vars) == 1:
                bar_subplot(ax, df_plot, plot_vars[0])
            elif len(plot_vars) == 2:
                scatter_subplot(ax, df_plot, plot_vars)

        for j in range(i + 1, len(axes)):
            axes[j].set_visible(False)

        handles = [
            plt.Rectangle((0, 0), 1, 1, facecolor=RXN_COLORS[rxn], edgecolor="black", label=rxn)
            for rxn in sorted(df["rxn_id"].unique())
        ]
        fig.legend(handles=handles, loc="upper right", ncol=1, fontsize=8, frameon=False)

        fig.suptitle(f"{protein} - organisms: {', '.join(SELECTED_ORGANISMS)}", fontsize=12)
        fig.tight_layout(rect=(0, 0, 1, 0.93))

        out = OUTPUT_DIR / f"boltz_multiplot_{protein}.svg"
        fig.savefig(out, dpi=300)
        plt.close(fig)
        print(f"saved {out}")