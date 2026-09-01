import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

BOLTZ_DIR = Path(
    "/home/pedro/Desktop/projects/ayahuasca/pathways/harmine_biosynthesis"
    "/proteins/stricto/analysis/structural/boltz"
)
OUTPUT = Path(
    "/home/pedro/Desktop/projects/ayahuasca/pathways/harmine_biosynthesis"
    "/paper/figures/boltz_specificity_plot.png"
)

TEMPERATURE_C = 20
R = 1.9872042

OOI_PREFIXES = {"caupuri", "tucunaca"}

# ---------------------------------------------------------------------------
# Thermodynamic conversion
# ---------------------------------------------------------------------------

def pred_to_dG_kcal(pred_value: float) -> float:
    kd = math.pow(10, pred_value) / 1_000_000
    T = 273.15 + TEMPERATURE_C
    return R * T * math.log(kd) / 1000

def pred_to_pIC50(pred_value: float) -> float:
    return 6.0 - pred_value

# ---------------------------------------------------------------------------
# Collect all affinity JSONs
# ---------------------------------------------------------------------------

def find_affinity_files(root: Path) -> list[tuple[Path, str]]:
    files = []
    for p in root.rglob("affinity_*.json"):
        name = p.stem.replace("affinity_", "")
        files.append((p, name))
    return files

def load_affinity(path: Path) -> float | None:
    try:
        data = json.loads(path.read_text())
        return data.get("affinity_pred_value")
    except Exception:
        return None

def parse_name(name: str) -> dict:
    parts = name.split("_")
    rxn = parts[-1]
    assert rxn in ("r1", "rr"), f"Unknown reaction: {rxn}"
    core_parts = parts[:-1]
    if core_parts[0] in OOI_PREFIXES:
        return {"ref": "_".join(core_parts[1:]), "reaction": rxn, "source": core_parts[0]}
    return {"ref": "_".join(core_parts), "reaction": rxn, "source": "self"}

affinity_files = find_affinity_files(BOLTZ_DIR)
print(f"Found {len(affinity_files)} affinity JSON files\n")

records = {}
for path, name in affinity_files:
    val = load_affinity(path)
    if val is None:
        continue
    info = parse_name(name)
    ref = info["ref"]
    rxn = info["reaction"]
    src = info["source"]
    if ref not in records:
        records[ref] = {"r1": {}, "rr": {}}
    records[ref][rxn][src] = val

# ---------------------------------------------------------------------------
# Impute missing rr values
# Use the r1/rr ratio from known pairs to estimate rr for r1-only entries
# ---------------------------------------------------------------------------

known_pairs = []
for ref in records:
    rec = records[ref]
    if rec["r1"] and rec["rr"]:
        r1_val = rec["r1"].get("self") or list(rec["r1"].values())[0]
        rr_val = rec["rr"].get("self") or list(rec["rr"].values())[0]
        known_pairs.append((r1_val, rr_val))

avg_r1_div_rr = 1.0
if known_pairs:
    r1_pIC50s = [pred_to_pIC50(r1) for r1, rr in known_pairs]
    rr_pIC50s = [pred_to_pIC50(rr) for r1, rr in known_pairs]
    ratios = [r1p / rrp for r1p, rrp in zip(r1_pIC50s, rr_pIC50s)]
    avg_r1_div_rr = sum(ratios) / len(ratios)
    print(f"Known r1/rr ratios (pIC50): {[f'{r:.4f}' for r in ratios]}")
    print(f"Average r1/rr ratio: {avg_r1_div_rr:.4f}")

# Build final point list
points = []
imputed = []

for ref in sorted(records):
    rec = records[ref]
    for src, r1_val in rec["r1"].items():
        if rec["rr"]:
            rr_val = list(rec["rr"].values())[0]
            imputed_flag = False
        else:
            r1_pIC50_tmp = pred_to_pIC50(r1_val)
            rr_pIC50_est = r1_pIC50_tmp / avg_r1_div_rr
            rr_val = 6.0 - rr_pIC50_est
            imputed_flag = True

        r1_pIC50 = pred_to_pIC50(r1_val)
        rr_pIC50 = pred_to_pIC50(rr_val)
        ratio = r1_pIC50 / rr_pIC50
        r1_dG = pred_to_dG_kcal(r1_val)
        rr_dG = pred_to_dG_kcal(rr_val)

        points.append({
            "label": ref, "source": src,
            "r1_pIC50": r1_pIC50, "rr_pIC50": rr_pIC50, "ratio": ratio,
            "r1_dG": r1_dG, "rr_dG": rr_dG,
            "imputed": imputed_flag,
        })
        if imputed_flag:
            imputed.append(ref)

# ---------------------------------------------------------------------------
# Print results table
# ---------------------------------------------------------------------------

print(f"\n{'Label':55s} {'src':12s} {'r1_pIC50':>10s} {'rr_pIC50':>10s} {'ratio':>8s} {'r1_ΔG':>10s} {'rr_ΔG':>10s} {'note':>10s}")
print("-" * 125)

for pt in points:
    note = "imputed" if pt["imputed"] else ""
    print(f"{pt['label']:55s} {pt['source']:12s} {pt['r1_pIC50']:10.4f} {pt['rr_pIC50']:10.4f} {pt['ratio']:8.4f} {pt['r1_dG']:10.4f} {pt['rr_dG']:10.4f} {note:>10s}")

# ---------------------------------------------------------------------------
# Plot
# x = pIC50_r1, y = pIC50_r1 / pIC50_rr
# ---------------------------------------------------------------------------

fig, ax = plt.subplots(figsize=(9, 7))

for pt in points:
    x = pt["r1_pIC50"]
    y = pt["ratio"]
    if pt["imputed"]:
        marker = "s"
        color = "lightcoral"
        ec = "firebrick"
        s = 140
    else:
        marker = "o"
        color = "steelblue"
        ec = "darkblue"
        s = 160
    ax.scatter(x, y, s=s, marker=marker, zorder=5,
               color=color, edgecolors=ec, linewidth=0.8, alpha=0.85)
    label = pt["label"]
    if pt["source"] != "self":
        label += f" ({pt['source']})"
    ax.annotate(
        label, (x, y),
        textcoords="offset points", xytext=(10, 10),
        fontsize=7.5, fontweight="bold",
        arrowprops=dict(arrowstyle="->", color="grey", lw=0.5),
    )

ax.set_xlabel("pIC\u2085\u2080 r1 reaction  (\u2212log\u2081\u2080IC\u2085\u2080, M)", fontsize=12)
ax.set_ylabel("pIC\u2085\u2080(r1) / pIC\u2085\u2080(rr)  (specificity ratio)", fontsize=12)
ax.set_title("StriCTO Boltz specificity: r1 vs rr binding affinity", fontsize=13)
ax.axhline(1.0, color="grey", linestyle="--", alpha=0.5, linewidth=0.8)
ax.margins(0.2)

from matplotlib.lines import Line2D
legend_elements = [
    Line2D([0], [0], marker="o", color="w", markerfacecolor="steelblue",
           markeredgecolor="darkblue", markersize=10, label="measured r1 + rr"),
    Line2D([0], [0], marker="s", color="w", markerfacecolor="lightcoral",
           markeredgecolor="firebrick", markersize=10, label="rr imputed"),
]
ax.legend(handles=legend_elements, loc="best", fontsize=9)

fig.tight_layout()
OUTPUT.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(OUTPUT, dpi=300)
print(f"\nPlot saved to {OUTPUT}")
