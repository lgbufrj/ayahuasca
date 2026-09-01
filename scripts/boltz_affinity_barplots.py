
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from mpl_toolkits.mplot3d import Axes3D

# ============================================
# CONFIGURATION
# ============================================
# Output files and figures paths
CSV = Path(
    "/home/pedro/Desktop/projects/ayahuasca/pathways/harmine_biosynthesis"
    "/paper/tables/boltz_affinity_results.csv"
)
OUTPUT_DIR = Path(
    "/home/pedro/Desktop/projects/ayahuasca/pathways/harmine_biosynthesis"
    "/paper/figures/supplementary"
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
]

# SELECTED VARIABLES (1 variable = barplot, 2 = scatterplot, 3 = 3D scatter)
SELECTED_VARIABLES = [
    "dG",
    "delta_dG", #tem que ter um ddG entre o rr e r1 (dg_rxn)
    "plddt"
    
]

# All available organisms:
ORGANISMS = ["arabidopsis", "catharanthus", "tabaco"]

# SELECTED ORGANISMS (one or more)
SELECTED_ORGANISMS = [
    "tabaco",
    "catharanthus",
    "arabidopsis"
]

# Colors used for r1 and rr when only one organism is selected
RXN_COLORS = {"r1": "#4C72B0", "rr": "#DD8452"}

# Up to 6 different colors used for the organism + reaction key
# when two or more organisms are selected
COMBO_COLORS = [
    "#4C72B0",
    "#DD8452",
    "#55A868",
    "#C44E52",
    "#8172B3",
    "#937860",
]

# ============================================
# DATA LOADING
# ============================================
df_all = pd.read_csv(CSV, skipinitialspace=True)
df_all.columns = [c.strip() for c in df_all.columns]
for col in df_all.select_dtypes(include="object").columns:
    df_all[col] = df_all[col].str.strip()
df_all = df_all.sort_values(
    ["reference_organism", "id", "rxn_id"], kind="stable"
)

# ============================================
# HELPERS
# ============================================
def build_combos(df, orgs):
    # Builds the (organism, reaction) pairs present in the data
    combos = []
    for org in orgs:
        for rxn in ("r1", "rr"):
            if ((df["reference_organism"] == org) & (df["rxn_id"] == rxn)).any():
                combos.append((org, rxn))
    return combos


def row_colors(df, orgs):
    # One color per row, based on the (organism, reaction) combo
    combos = build_combos(df, orgs)
    combo_color = {
        combo: COMBO_COLORS[i % len(COMBO_COLORS)] for i, combo in enumerate(combos)
    }
    return [combo_color[(row["reference_organism"], row["rxn_id"])] for _, row in df.iterrows()]


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
def bar_plot(df, orgs, variables):
    # Barplot used when exactly 1 variable is selected
    fig, ax = plt.subplots(figsize=(12, 5))

    labels = [
        f"{row['id']}\n{row['rxn_id']}"
        for _, row in df[["id", "rxn_id"]].iterrows()
    ]
    colors = RXN_COLORS if len(orgs) == 1 else None
    if colors is None:
        combos = build_combos(df, orgs)
        combo_color = {
            combo: COMBO_COLORS[i % len(COMBO_COLORS)] for i, combo in enumerate(combos)
        }
        colors = {
            (row["reference_organism"], row["rxn_id"]): combo_color[
                (row["reference_organism"], row["rxn_id"])
            ]
            for _, row in df.iterrows()
        }
    bar_colors = [
        RXN_COLORS[row["rxn_id"]] if len(orgs) == 1 else colors[
            (row["reference_organism"], row["rxn_id"])
        ]
        for _, row in df.iterrows()
    ]

    ax.bar(range(len(df)), df[variables[0]], color=bar_colors, edgecolor="black", linewidth=0.3)

    # Zero reference line for delta variables
    if variables[0].startswith("delta"):
        ax.axhline(0.0, color="grey", linestyle="--", linewidth=0.8)

    ax.set_xticks(range(len(df)))
    ax.set_xticklabels(labels, fontsize=6, rotation=90)
    ax.set_ylabel(variables[0], fontsize=9)
    ax.margins(x=0.01)

    # Plot legend/key
    fig.legend(
        handles=legend_handles(df, orgs),
        loc="upper right",
        ncol=1,
        fontsize=8,
        frameon=False,
    )
    title = f"{variables[0]} - organisms: {', '.join(orgs)}"
    fig.suptitle(title, fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.93))

    out = OUTPUT_DIR / f"boltz_barplot_{'_'.join(variables)}.png"
    fig.savefig(out, dpi=300)
    plt.close(fig)
    print(f"saved {out}")


def scatter_plot(df, orgs, variables):
    # Scatterplot used when exactly 2 variables are selected
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(df[variables[0]], df[variables[1]], c=row_colors(df, orgs), s=40)

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
        f"{variables[0]} vs {variables[1]} - organisms: {', '.join(orgs)}",
        fontsize=12,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.93))

    out = OUTPUT_DIR / f"boltz_scatter_{'_'.join(variables)}.png"
    fig.savefig(out, dpi=300)
    plt.close(fig)
    print(f"saved {out}")


def plot_3d(df, orgs, variables):
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
        f"{variables[0]} vs {variables[1]} vs {variables[2]} "
        f"- organisms: {', '.join(orgs)}",
        fontsize=12,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.93))

    out = OUTPUT_DIR / f"boltz_3d_{'_'.join(variables)}.png"
    fig.savefig(out, dpi=300)
    plt.close(fig)
    print(f"saved {out}")


# ============================================
# MAIN EXECUTION (no edits needed below)
# ============================================
# Configuration validation
if not 1 <= len(SELECTED_VARIABLES) <= 3:
    raise ValueError("SELECTED_VARIABLES must contain between 1 and 3 variables.")
if not SELECTED_ORGANISMS:
    raise ValueError("SELECTED_ORGANISMS must contain at least one organism.")

unknown_vars = [v for v in SELECTED_VARIABLES if v not in VARIABLES]
if unknown_vars:
    raise ValueError(f"Unknown variables in SELECTED_VARIABLES: {unknown_vars}")

unknown_orgs = [o for o in SELECTED_ORGANISMS if o not in ORGANISMS]
if unknown_orgs:
    raise ValueError(f"Unknown organisms in SELECTED_ORGANISMS: {unknown_orgs}")

if __name__ == "__main__":
    # Select only the rows of the configured organisms and reactions r1/rr
    df = df_all[
        df_all["reference_organism"].isin(SELECTED_ORGANISMS)
        & df_all["rxn_id"].isin(["r1", "rr"])
    ].copy()

    # Drop rows without data for the selected variables
    df = df.dropna(subset=SELECTED_VARIABLES)

    if df.empty:
        print("No data available for the selected configuration.")
    else:
        OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

        # 1 variable -> barplot, 2 variables -> scatterplot, 3 -> 3D scatter
        if len(SELECTED_VARIABLES) == 1:
            bar_plot(df, SELECTED_ORGANISMS, SELECTED_VARIABLES)
        elif len(SELECTED_VARIABLES) == 2:
            scatter_plot(df, SELECTED_ORGANISMS, SELECTED_VARIABLES)
        else:
            plot_3d(df, SELECTED_ORGANISMS, SELECTED_VARIABLES)