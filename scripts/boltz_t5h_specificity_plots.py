from __future__ import annotations

import json
import math
from pathlib import Path
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

T5H_DIR = Path(
    "/home/pedro/Desktop/projects/ayahuasca/pathways/harmine_biosynthesis"
    "/proteins/t5h/analysis/structural/boltz"
)
OUT_DIR = Path(
    "/home/pedro/Desktop/projects/ayahuasca/pathways/harmine_biosynthesis"
    "/paper/figures"
)
OUT_DIR.mkdir(parents=True, exist_ok=True)

TEMPERATURE_C = 20
R = 1.9872042


def pred_to_dG_kcal(pred_value: float) -> float:
    kd = math.pow(10, pred_value) / 1_000_000
    T = 273.15 + TEMPERATURE_C
    return R * T * math.log(kd) / 1000


def pred_to_pIC50(pred_value: float) -> float:
    return 6.0 - pred_value


def load_affinity(path: Path) -> float | None:
    try:
        data = json.loads(path.read_text())
        return data.get("affinity_pred_value")
    except Exception:
        return None


def parse_t5h_results(results_path: Path) -> dict:
    rel = results_path.relative_to(T5H_DIR)
    parts = rel.parts

    if len(parts) == 3:
        predictor, target, reaction = parts
    elif len(parts) == 2:
        predictor, reaction = parts
        target = predictor  # self-prediction
    else:
        raise ValueError(f"Unexpected path: {results_path}")

    return {
        "predictor": predictor,
        "target": target,
        "reaction": reaction,
        "source": f"{predictor}_{target}",
        "label": f"{predictor} \u2192 {target}" if predictor != target else f"{predictor} self",
    }


affinity_files = sorted(T5H_DIR.rglob("affinity_*.json"))
print(f"Found {len(affinity_files)} affinity JSON files\n")

records = defaultdict(dict)
for aff_path in affinity_files:
    # path: .../{predictor}/{target}/{reaction}/predictions/{name}/affinity_*.json
    results_dir = aff_path.parent.parent.parent  # .../{predictor}/{target}/{reaction}/
    parsed = parse_t5h_results(results_dir)
    rxn = parsed["reaction"]
    src = parsed["source"]
    val = load_affinity(aff_path)
    if val is None:
        continue
    records[src][rxn] = {
        "value": val,
        "pIC50": pred_to_pIC50(val),
        "dG": pred_to_dG_kcal(val),
        "label": parsed["label"],
        "predictor": parsed["predictor"],
        "target": parsed["target"],
    }

REACTIONS = ["r1", "r2", "r3"]
COLORS = {"r1": "#E74C3C", "r2": "#3498DB", "r3": "#2ECC71"}
MARKERS = {"r1": "o", "r2": "s", "r3": "^"}

print(f"\n{'Label':35s} rxn  r(n)_aff  rr_aff   ratio  ratio_pIC")
print("-" * 80)

all_points = []
for rxn in REACTIONS:
    for src, rec in sorted(records.items()):
        if rxn not in rec or "rr" not in rec:
            continue
        r_val = rec[rxn]
        rr_val = rec["rr"]
        ratio = r_val["value"] / rr_val["value"]
        ratio_pIC = r_val["pIC50"] / rr_val["pIC50"]
        print(f"  {rec[rxn]['label']:35s} {rxn:3s} {r_val['value']:8.4f} {rr_val['value']:7.4f} {ratio:7.4f} {ratio_pIC:9.4f}")
        all_points.append({
            "label": rec[rxn]["label"],
            "r_val": r_val["value"],
            "ratio": ratio,
            "rxn": rxn,
        })

fig, ax = plt.subplots(figsize=(9, 7))

for pt in all_points:
    rxn = pt["rxn"]
    ax.scatter(
        pt["r_val"], pt["ratio"],
        s=180, marker=MARKERS[rxn], zorder=5,
        color=COLORS[rxn], edgecolors="black",
        linewidth=0.6, alpha=0.8,
    )

ax.set_xlabel("affinity", fontsize=12)
ax.set_ylabel("r(n)/rr", fontsize=12)
ax.set_title("T5H Boltz specificity by reaction", fontsize=13)
ax.axhline(1.0, color="grey", linestyle="--", alpha=0.5, linewidth=0.8)
ax.margins(0.15)

from matplotlib.patches import FancyBboxPatch

formula_text = (
    "Formula:\n"
    "x = affinity(r(n))\n"
    "y = affinity(r(n)) / affinity(rr)\n"
    "affinity = Boltz predicted pKd\n"
    "(lower = stronger binding)"
)
ax.text(
    1.02, 0.98, formula_text,
    transform=ax.transAxes,
    fontsize=8,
    verticalalignment="top",
    family="monospace",
    bbox=dict(boxstyle="round,pad=0.5", facecolor="ivory", alpha=0.9,
              edgecolor="grey"),
)

legend_elements = [
    Line2D([0], [0], marker=MARKERS["r1"], color="w", markerfacecolor=COLORS["r1"],
           markeredgecolor="black", markersize=8, label="r1 = TTH"),
    Line2D([0], [0], marker=MARKERS["r2"], color="w", markerfacecolor=COLORS["r2"],
           markeredgecolor="black", markersize=8, label="r2 = designed binder"),
    Line2D([0], [0], marker=MARKERS["r3"], color="w", markerfacecolor=COLORS["r3"],
           markeredgecolor="black", markersize=8, label="r3 = designed cofactor+binder"),
    Line2D([0], [0], marker="", color="w", label="rr = reference"),
]
legend = ax.legend(
    handles=legend_elements, loc="upper left",
    bbox_to_anchor=(1.02, 0.55),  # below the formula box
    fontsize=7.5, framealpha=0.9, edgecolor="grey",
    title=None,
)

fig.tight_layout()
out_path = OUT_DIR / "t5h_specificity_combined.png"
fig.savefig(out_path, dpi=300)
plt.close(fig)
print(f"\nPlot saved: {out_path}")
