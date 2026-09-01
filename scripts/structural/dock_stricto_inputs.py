"""
dock_stricto_inputs.py
───────────────────────────────────────────────────────────────────────────
Dock all 32 tucunaca Boltz-input YAMLs using Vina (CPU-only).
Uses existing Boltz CIF structures for reference proteins and the
AlphaFold model for tucunaca's own sequences.
───────────────────────────────────────────────────────────────────────────
"""

import json, math, re, subprocess, sys, pickle, copy
from pathlib import Path
from collections import defaultdict

import numpy as np
import yaml
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

BASE = Path("/home/pedro/Desktop/projects/ayahuasca/pathways/harmine_biosynthesis")
INPUTS_DIR = BASE / "proteins" / "stricto" / "analysis" / "structural" / "boltz" / "tucunaca" / "inputs"
BOLTZ_DIR = BASE / "proteins" / "stricto" / "analysis" / "structural" / "boltz"
AF_MODEL = BASE / "proteins" / "stricto" / "analysis" / "structural" / "alphafold" / "tucunaca" / "best_model.pdb"
OUT_DIR = BASE / "proteins" / "stricto" / "analysis" / "structural" / "vina" / "tucunaca_dock"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def find_cif(protein_key: str, reaction: str) -> Path | None:
    """Find a Boltz CIF for a given protein key (checks r1 first, then rr)."""
    orgs = ["arabidopsis", "catharanthus", "tabaco", "caupuri", "tucunaca"]
    for rxn in [reaction, "r1"]:
        for org in orgs:
            base = BOLTZ_DIR / org / rxn / "predictions"
            if not base.exists():
                continue
            for subdir in base.iterdir():
                if not subdir.is_dir():
                    continue
                cifs = list(subdir.glob("*.cif"))
                for cif in cifs:
                    if protein_key in cif.stem:
                        return cif
    return None


def parse_cif_atoms(cif_path: Path) -> list[dict]:
    """Parse ATOM and HETATM records from a Boltz CIF file."""
    with open(cif_path) as f:
        lines = f.readlines()

    header = []
    header_start = header_end = None
    for i, line in enumerate(lines):
        if line.startswith("_atom_site."):
            if header_start is None:
                header_start = i
            col = line.strip().split(".")[1].strip()
            header.append(col)
        elif header_start is not None and header_end is None:
            header_end = i
            break

    if header_end is None:
        return []

    data_start = header_end
    data_end = None
    for i in range(data_start, len(lines)):
        stripped = lines[i].strip()
        if stripped == "#" or stripped.startswith("loop_"):
            data_end = i
            break

    atoms = []
    for line in lines[data_start:data_end]:
        stripped = line.strip()
        if not stripped:
            continue
        parts = stripped.split()
        if len(parts) >= len(header):
            rec = dict(zip(header, parts))
            atoms.append(rec)
    return atoms


def cif_to_pdb(cif_path: Path, pdb_path: Path, chain: str = "A") -> Path:
    """Extract protein chain A from a Boltz CIF file, save as PDB."""
    pdb_path.parent.mkdir(parents=True, exist_ok=True)
    atoms = parse_cif_atoms(cif_path)

    atom_i = 0
    res_i = 0
    last_seq = None
    lines_out = []
    for a in atoms:
        if a.get("group_PDB") != "ATOM":
            continue
        if a.get("label_asym_id") != chain:
            continue
        label_comp_id = a.get("label_comp_id", "UNK")
        if label_comp_id in ("LIG1", "LIG2"):
            continue
        atom_i += 1
        seq = a.get("label_seq_id", "0")
        if seq != last_seq:
            res_i += 1
            last_seq = seq

        elem = a.get("type_symbol", "X")
        atom_name = a.get("label_atom_id", "X")

        aname = f" {atom_name:<3s}" if len(atom_name) < 4 else atom_name[:4]

        x = float(a.get("Cartn_x", 0))
        y = float(a.get("Cartn_y", 0))
        z = float(a.get("Cartn_z", 0))
        b = float(a.get("B_iso_or_equiv", 0))

        pdb_line = (
            f"ATOM  {atom_i:5d} {aname:<4s} {label_comp_id:>3s}"
            f" A{res_i:4d}    {x:8.3f}{y:8.3f}{z:8.3f}"
            f"  1.00{b:6.2f}          {elem:>2s}\n"
        )
        lines_out.append(pdb_line)

    pdb_path.write_text("".join(lines_out) + "END\n")
    print(f"    Wrote {len(lines_out)} atoms to {pdb_path.name}")
    return pdb_path


def extract_ligand_coords(cif_path: Path, ligand_id: str = "LIG2") -> tuple[float, float, float] | None:
    """Extract center of a ligand's coordinates from Boltz CIF."""
    atoms = parse_cif_atoms(cif_path)
    coords = []
    for a in atoms:
        if a.get("label_comp_id") == ligand_id:
            x = float(a.get("Cartn_x", 0))
            y = float(a.get("Cartn_y", 0))
            z = float(a.get("Cartn_z", 0))
            coords.append((x, y, z))
    if not coords:
        return None
    arr = np.array(coords)
    return tuple(arr.mean(axis=0).tolist())


def protein_center(pdb_path: Path) -> tuple[float, float, float]:
    """Compute center of mass from PDB coordinates."""
    xs, ys, zs = [], [], []
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("ATOM"):
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                xs.append(x)
                ys.append(y)
                zs.append(z)
    if xs:
        return (sum(xs)/len(xs), sum(ys)/len(ys), sum(zs)/len(zs))
    return (0.0, 0.0, 0.0)


def smi_to_sdf(smi: str, name: str, out_path: Path) -> Path:
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smi}")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    AllChem.MMFFOptimizeMolecule(mol)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    Chem.MolToMolFile(mol, str(out_path))
    return out_path


def pdb_to_pdbqt(pdb_path: Path, pdbqt_path: Path) -> None:
    """Convert PDB to PDBQT using Open Babel (works for Boltz CIF-extracted PDBs)."""
    pdbqt_path.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        "obabel", str(pdb_path),
        "-O", str(pdbqt_path),
        "-xr",
        "--partialcharge", "gasteiger",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"    [WARN] obabel failed: {result.stderr.strip()[:200]}")
    else:
        print(f"    Receptor PDBQT: {pdbqt_path}")


def lig_sdf_to_pdbqt(sdf_path: Path, pdbqt_path: Path) -> None:
    """Convert ligand SDF to PDBQT using mk_prepare_ligand.py."""
    pdbqt_path.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        "mk_prepare_ligand.py",
        "-i", str(sdf_path),
        "-o", str(pdbqt_path),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"    [WARN] mk_prepare_ligand: {result.stderr.strip()[:200]}")


def run_vina(receptor_pdbqt: Path, ligand_pdbqt: Path,
             center: tuple, size: tuple = (20, 20, 20),
             output: Path | None = None,
             exhaustiveness: int = 8, num_modes: int = 10) -> float:
    """Run Vina docking, return best affinity (kcal/mol)."""
    output = output or OUT_DIR / "temp_out.pdbqt"
    output.parent.mkdir(parents=True, exist_ok=True)

    cx, cy, cz = center
    sx, sy, sz = size
    config = (
        f"receptor = {receptor_pdbqt}\n"
        f"ligand   = {ligand_pdbqt}\n"
        f"center_x = {cx:.3f}\n"
        f"center_y = {cy:.3f}\n"
        f"center_z = {cz:.3f}\n"
        f"size_x   = {sx}\n"
        f"size_y   = {sy}\n"
        f"size_z   = {sz}\n"
        f"exhaustiveness = {exhaustiveness}\n"
        f"num_modes      = {num_modes}\n"
    )
    config_path = output.with_suffix(".cfg")
    config_path.write_text(config)

    log_path = output.with_suffix(".log")
    cmd = ["vina", "--config", str(config_path), "--out", str(output), "--log", str(log_path)]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        return 999.0

    if log_path.exists():
        log_text = log_path.read_text()
    else:
        log_text = result.stderr

    best_aff = 999.0
    for line in log_text.splitlines():
        m = re.match(r"\s*\d+\s+(-?\d+\.?\d*)", line)
        if m:
            val = float(m.group(1))
            if val < best_aff:
                best_aff = val
    return best_aff


def pred_to_pIC50(pred_value: float) -> float:
    return 6.0 - pred_value / 1.36  # approximate conversion


# ── pipeline ────────────────────────────────────────────────────────────────

print("=" * 80)
print("StriCTO docking pipeline (CPU / Vina)")
print(f"Inputs: {INPUTS_DIR}")
print("=" * 80)

yaml_files = sorted(INPUTS_DIR.glob("*.yaml"))
print(f"\nFound {len(yaml_files)} YAML input files\n")

results = []

for yaml_path in yaml_files:
    name = yaml_path.stem
    print(f"\n── {name} ──")

    with open(yaml_path) as f:
        data = yaml.safe_load(f)

    # Extract protein sequence and reaction
    protein_seq = ""
    for ent in data.get("sequences", []):
        if "protein" in ent:
            protein_seq = ent["protein"]["sequence"]
    rxn = "rr" if "_rr" in name else "r1"

    # Protein key (name without _r1/_rr)
    protein_key = name.rsplit("_", 1)[0]

    # ── find / prepare receptor ──
    cif_path = find_cif(protein_key, rxn)
    pdb_path = OUT_DIR / "receptors" / f"{protein_key}.pdb"
    use_af = False

    if cif_path and cif_path.exists():
        print(f"  CIF: {cif_path}")
        if not pdb_path.exists():
            cif_to_pdb(cif_path, pdb_path)
        center = extract_ligand_coords(cif_path, "LIG2")
        if center is None:
            center = extract_ligand_coords(cif_path, "LIG1")
        grid_size = (20, 20, 20)
        print(f"  Grid center (from LIG): ({center[0]:.1f}, {center[1]:.1f}, {center[2]:.1f})")
    elif AF_MODEL.exists() and ("TucunacaG" in protein_key or "tucunaca" in protein_key.lower()):
        use_af = True
        if not pdb_path.exists():
            pdb_path.parent.mkdir(parents=True, exist_ok=True)
            pdb_path.write_bytes(AF_MODEL.read_bytes())
            print(f"  Copied AlphaFold model ({AF_MODEL.stat().st_size} bytes)")
        else:
            print(f"  AlphaFold model (cached)")
        center = protein_center(pdb_path)
        grid_size = (40, 40, 40)
        print(f"  Grid center (protein COM): ({center[0]:.1f}, {center[1]:.1f}, {center[2]:.1f})")
    else:
        print(f"  [SKIP] No structure for {protein_key}")
        continue

    # ── receptor PDBQT ──
    pdbqt_path = pdb_path.with_suffix(".pdbqt")
    if not pdbqt_path.exists():
        pdb_to_pdbqt(pdb_path, pdbqt_path)
    if not pdbqt_path.exists():
        print(f"  [SKIP] Failed to prepare receptor PDBQT")
        continue

    # ── ligand ──
    lig_smi = None
    for ent in data.get("sequences", []):
        if "ligand" in ent and ent["ligand"].get("id") == "C":
            lig_smi = ent["ligand"]["smiles"]
            break
    if not lig_smi:
        print(f"  [SKIP] No ligand C SMILES")
        continue

    lig_name = f"{name}_C"
    lig_sdf = OUT_DIR / "ligands" / f"{lig_name}.sdf"
    lig_pdbqt = OUT_DIR / "ligands" / f"{lig_name}.pdbqt"

    if not lig_sdf.exists():
        try:
            smi_to_sdf(lig_smi, lig_name, lig_sdf)
        except Exception as e:
            print(f"  [SKIP] Ligand SDF: {e}")
            continue

    if not lig_pdbqt.exists():
        lig_sdf_to_pdbqt(lig_sdf, lig_pdbqt)

    if not lig_pdbqt.exists():
        print(f"  [SKIP] Failed to prepare ligand PDBQT")
        continue

    # ── dock ──
    out_pdbqt = OUT_DIR / "docked" / f"{name}.pdbqt"
    try:
        best_aff = run_vina(pdbqt_path, lig_pdbqt, center=center,
                            size=grid_size, output=out_pdbqt)
    except Exception as e:
        print(f"  [SKIP] Docking: {e}")
        continue

    pic50 = pred_to_pIC50(best_aff) if best_aff < 50 else 0.0
    print(f"  \u0394G = {best_aff:.2f} kcal/mol  pIC\u2085\u2080 = {pic50:.3f}")

    results.append({
        "name": name,
        "protein_key": protein_key,
        "reaction": rxn,
        "affinity": best_aff,
        "pIC50": pic50,
        "grid_center": center,
    })

# ── report & plot ────────────────────────────────────────────────────────────
print("\n" + "=" * 80)
print(f"Completed: {len(results)} / {len(yaml_files)} docking runs")
print("=" * 80)

by_protein = defaultdict(lambda: {"r1": None, "rr": None})
for r in results:
    by_protein[r["protein_key"]][r["reaction"]] = r

print(f"\n{'Protein':55s} {'r1_aff':>8s} {'rr_aff':>8s} {'r1_pIC50':>9s} {'rr_pIC50':>9s} {'ratio':>8s}")
print("-" * 100)

plot_points = []
known_ratios = []

for pk, rec in sorted(by_protein.items()):
    r1_r = rec.get("r1")
    rr_r = rec.get("rr")
    v = lambda x: f"{x['affinity']:.2f}" if x else "N/A"
    p = lambda x: f"{x['pIC50']:.3f}" if x else "N/A"
    r1_str, rr_str = v(r1_r), v(rr_r)
    r1_p_str, rr_p_str = p(r1_r), p(rr_r)

    if r1_r and rr_r:
        rat = r1_r["pIC50"] / rr_r["pIC50"] if rr_r["pIC50"] else 0
        ratio_s = f"{rat:.4f}"
        plot_points.append((pk, r1_r["pIC50"], rr_r["pIC50"], rat, False))
        known_ratios.append(rat)
    else:
        ratio_s = "N/A"

    print(f"{pk:55s} {r1_str:>8s} {rr_str:>8s} {r1_p_str:>9s} {rr_p_str:>9s} {ratio_s:>8s}")

# Impute missing rr
if known_ratios:
    avg_ratio = sum(known_ratios) / len(known_ratios)
    print(f"\nAverage r1/rr pIC50 ratio from {len(known_ratios)} pairs: {avg_ratio:.4f}")
    n_imputed = 0
    for pk, rec in sorted(by_protein.items()):
        r1_r = rec.get("r1")
        rr_r = rec.get("rr")
        if r1_r and not rr_r:
            rr_pIC50_est = r1_r["pIC50"] / avg_ratio
            plot_points.append((pk, r1_r["pIC50"], rr_pIC50_est, avg_ratio, True))
            n_imputed += 1
            print(f"  Imputed rr for {pk}: pIC50_rr = {rr_pIC50_est:.3f}")
    print(f"Imputed {n_imputed} missing rr values")

# Plot
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

fig, ax = plt.subplots(figsize=(12, 9))

for pk, r1_p, rr_p, ratio, is_imp in plot_points:
    x = r1_p
    y = ratio
    if is_imp:
        marker = "s"
        c = "lightcoral"
        ec = "firebrick"
        label_short = pk.replace("_stricto_", "\n")
    else:
        marker = "o"
        c = "steelblue"
        ec = "darkblue"
        label_short = pk.replace("_stricto_", "\n")
    ax.scatter(x, y, s=180, marker=marker, zorder=5,
               color=c, edgecolors=ec, linewidth=0.8, alpha=0.85)
    ax.annotate(pk, (x, y),
                textcoords="offset points", xytext=(10, 10),
                fontsize=6, fontweight="bold",
                arrowprops=dict(arrowstyle="->", color="grey", lw=0.4))

ax.set_xlabel("pIC\u2085\u2080 (r1 reaction)", fontsize=12)
ax.set_ylabel("pIC\u2085\u2080(r1) / pIC\u2085\u2080(rr)  (specificity ratio)", fontsize=12)
ax.set_title("StriCTO Vina docking specificity (tucunaca inputs)", fontsize=13)
ax.axhline(1.0, color="grey", linestyle="--", alpha=0.5, linewidth=0.8)
ax.margins(0.2)

legend_elements = [
    Line2D([0],[0], marker="o", color="w", markerfacecolor="steelblue",
           markeredgecolor="darkblue", markersize=10, label="measured r1+rr"),
    Line2D([0],[0], marker="s", color="w", markerfacecolor="lightcoral",
           markeredgecolor="firebrick", markersize=10, label="rr imputed"),
]
ax.legend(handles=legend_elements, loc="best", fontsize=9)
fig.tight_layout()

plot_path = BASE / "paper" / "figures" / "vina_specificity_plot_tucunaca_inputs.png"
plot_path.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(plot_path, dpi=300)
print(f"\nPlot saved to {plot_path}")

# Save results
results_path = OUT_DIR / "docking_results.json"
results_path.write_text(json.dumps(results, indent=2, default=str))
print(f"Results saved to {results_path}")
