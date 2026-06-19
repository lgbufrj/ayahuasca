import csv
import json
import math
from pathlib import Path

from data import PROTEINS_PATH, proteins, PAPER_PATH


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

TEMPERATURE_C = 20          # °C  — change as needed
R = 1.9872042               # cal/(mol·K)
OUTPUT_CSV = f"{PAPER_PATH}/tables/boltz_affinity_results.csv"

CSV_COLUMNS = [
    "protein", "organism", "reference_organism", "cofactor", "substrate",
    "binding_affinity", "kd", "dG", "IC50",
    "delta_affinity", "delta_kd", "delta_dG", "delta_IC50",
    "confidence", "delta_confidence",
    "plddt", "delta_plddt",
    "iplddt", "delta_iplddt",
    "pde", "delta_pde",
    "ipde", "delta_ipde",
    "temperature",
]


# ---------------------------------------------------------------------------
# Thermodynamic helpers
# ---------------------------------------------------------------------------

def pred_to_kd(pred_value: float) -> float:
    """affinity_pred_value (log10 IC50 in µM) → Kd in mol/L.
    Boltz uses IC50 in µM so 10^pred_value µM = 10^pred_value * 1e-6 M.
    We treat Kd ≈ IC50 (competitive inhibition, [S] << Km approximation).
    """
    return math.pow(10, pred_value) / 1_000_000


def kd_to_dG(kd: float, temperature_c: float) -> float:
    """Kd (mol/L) → ΔG in kJ/mol."""
    T = 273.15 + temperature_c
    return R * T * math.log(kd) / 1000  # kcal/mol


def pred_to_IC50_uM(pred_value: float) -> float:
    """affinity_pred_value → IC50 in µM (direct from model definition)."""
    return math.pow(10, pred_value)


# ---------------------------------------------------------------------------
# File helpers
# ---------------------------------------------------------------------------

def load_json(path: Path) -> dict | None:
    if not path.exists():
        print(f"  [WARN] Missing: {path}")
        return None
    with open(path) as f:
        return json.load(f)


def read_affinity(path: Path) -> dict | None:
    data = load_json(path)
    if data is None:
        return None
    pred = data.get("affinity_pred_value")
    if pred is None:
        print(f"  [WARN] affinity_pred_value missing in {path}")
        return None
    kd = pred_to_kd(pred)
    return {
        "binding_affinity": pred,
        "kd":               kd,
        "dG":               kd_to_dG(kd, TEMPERATURE_C),
        "IC50":             pred_to_IC50_uM(pred),
    }


def read_confidence(path: Path) -> dict | None:
    data = load_json(path)
    if data is None:
        return None
    return {
        "confidence": data.get("confidence_score"),
        "plddt":      data.get("complex_plddt"),
        "iplddt":     data.get("complex_iplddt"),
        "pde":        data.get("complex_pde"),
        "ipde":       data.get("complex_ipde"),
    }


# ---------------------------------------------------------------------------
# Cofactor / substrate resolution
# ---------------------------------------------------------------------------

PROTON_ALIASES = {"H+", "H⁺", "H(+)", "proton"}

def resolve_cofactor_substrate(rxn: dict, cofactor_list: list[str]) -> tuple[str, str]:
    """
    From rxn['substrates'], find which one is the cofactor (present in the
    protein's cofactor list) and which is the substrate, ignoring protons.
    Returns (cofactor_name, substrate_name) — either may be "" if not found.
    """
    substrates = [
        s for s in rxn.get("substrates", [])
        if s not in PROTON_ALIASES
    ]
    cofactor_set = set(cofactor_list)
    cofactor = next((s for s in substrates if s in cofactor_set), "")
    substrate = next((s for s in substrates if s not in cofactor_set), "")
    return cofactor, substrate


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def build_file_paths(
    boltz_analysis_dir: str,
    protein_name: str,
    uniprot_id: str,
    rxn_id: str,
    ref_organism: str,
    ooi: str | None,
) -> tuple[Path, Path]:
    if ooi:
        results_dir = Path(f"{boltz_analysis_dir}/{ooi}/{ref_organism}/{rxn_id}")
        affinity_file    = results_dir / f"predictions/{ooi}_{ref_organism}_{protein_name}_{uniprot_id}_{rxn_id}/affinity_{ooi}_{ref_organism}_{protein_name}_{uniprot_id}_{rxn_id}.json"
        confidence_file  = results_dir / f"predictions/{ooi}_{ref_organism}_{protein_name}_{uniprot_id}_{rxn_id}/confidence_{ooi}_{ref_organism}_{protein_name}_{uniprot_id}_{rxn_id}_model_0.json"
    else:
        results_dir = Path(f"{boltz_analysis_dir}/{ref_organism}/{rxn_id}")
        affinity_file    = results_dir / f"predictions/{ref_organism}_{protein_name}_{uniprot_id}_{rxn_id}/affinity_{ref_organism}_{protein_name}_{uniprot_id}_{rxn_id}.json"
        confidence_file  = results_dir / f"predictions/{ref_organism}_{protein_name}_{uniprot_id}_{rxn_id}/confidence_{ref_organism}_{protein_name}_{uniprot_id}_{rxn_id}_model_0.json"
    return affinity_file, confidence_file


def delta(ooi_val, ref_val):
    if ooi_val is None or ref_val is None:
        return None
    return ooi_val - ref_val


def main():
    rows = []

    for protein_name, protein_data in proteins.items():
        boltz_analysis_dir = f"{PROTEINS_PATH}/{protein_name}/analysis/structural/boltz"

        cofactor_list  = protein_data.get("cofactors", [])
        oois           = protein_data.get("oois", [])

        for rxn in protein_data["reactions"]:
            rxn_id = rxn["id"]
            cofactor, substrate = resolve_cofactor_substrate(rxn, cofactor_list)

            for ref_organism, organism_data in protein_data["organisms"].items():
                uniprot_id = organism_data["uniprot_id"]

                # --- reference organism row ---
                ref_aff_file, ref_conf_file = build_file_paths(
                    boltz_analysis_dir, protein_name, uniprot_id,
                    rxn_id, ref_organism, ooi=None,
                )
                ref_aff  = read_affinity(ref_aff_file)
                ref_conf = read_confidence(ref_conf_file)

                # Build reference row (deltas are all 0 for the ref itself)
                ref_row = {
                    "protein":           protein_name,
                    "organism":          ref_organism,
                    "reference_organism": ref_organism,
                    "cofactor":          cofactor,
                    "substrate":         substrate,
                    # affinity metrics
                    "binding_affinity":  ref_aff["binding_affinity"] if ref_aff else None,
                    "kd":                ref_aff["kd"]               if ref_aff else None,
                    "dG":                ref_aff["dG"]               if ref_aff else None,
                    "IC50":              ref_aff["IC50"]             if ref_aff else None,
                    "delta_affinity":    0.0 if ref_aff  else None,
                    "delta_kd":          0.0 if ref_aff  else None,
                    "delta_dG":          0.0 if ref_aff  else None,
                    "delta_IC50":        0.0 if ref_aff  else None,
                    # confidence metrics
                    "confidence":        ref_conf["confidence"] if ref_conf else None,
                    "delta_confidence":  0.0 if ref_conf else None,
                    "plddt":             ref_conf["plddt"]      if ref_conf else None,
                    "delta_plddt":       0.0 if ref_conf else None,
                    "iplddt":            ref_conf["iplddt"]     if ref_conf else None,
                    "delta_iplddt":      0.0 if ref_conf else None,
                    "pde":               ref_conf["pde"]        if ref_conf else None,
                    "delta_pde":         0.0 if ref_conf else None,
                    "ipde":              ref_conf["ipde"]       if ref_conf else None,
                    "delta_ipde":        0.0 if ref_conf else None,
                    "temperature":       TEMPERATURE_C,
                }
                rows.append(ref_row)

                # --- OOI rows ---
                for ooi in oois:
                    ooi_aff_file, ooi_conf_file = build_file_paths(
                        boltz_analysis_dir, protein_name, uniprot_id,
                        rxn_id, ref_organism, ooi=ooi,
                    )
                    ooi_aff  = read_affinity(ooi_aff_file)
                    ooi_conf = read_confidence(ooi_conf_file)

                    ooi_row = {
                        "protein":            protein_name,
                        "organism":           ooi,
                        "reference_organism": ref_organism,
                        "cofactor":           cofactor,
                        "substrate":          substrate,
                        # affinity metrics
                        "binding_affinity":   ooi_aff["binding_affinity"] if ooi_aff else None,
                        "kd":                 ooi_aff["kd"]               if ooi_aff else None,
                        "dG":                 ooi_aff["dG"]               if ooi_aff else None,
                        "IC50":               ooi_aff["IC50"]             if ooi_aff else None,
                        "delta_affinity":     delta(ooi_aff["binding_affinity"] if ooi_aff else None,
                                                    ref_aff["binding_affinity"] if ref_aff else None),
                        "delta_kd":           delta(ooi_aff["kd"] if ooi_aff else None,
                                                    ref_aff["kd"] if ref_aff else None),
                        "delta_dG":           delta(ooi_aff["dG"] if ooi_aff else None,
                                                    ref_aff["dG"] if ref_aff else None),
                        "delta_IC50":         delta(ooi_aff["IC50"] if ooi_aff else None,
                                                    ref_aff["IC50"] if ref_aff else None),
                        # confidence metrics
                        "confidence":         ooi_conf["confidence"] if ooi_conf else None,
                        "delta_confidence":   delta(ooi_conf["confidence"] if ooi_conf else None,
                                                    ref_conf["confidence"] if ref_conf else None),
                        "plddt":              ooi_conf["plddt"]  if ooi_conf else None,
                        "delta_plddt":        delta(ooi_conf["plddt"] if ooi_conf else None,
                                                    ref_conf["plddt"] if ref_conf else None),
                        "iplddt":             ooi_conf["iplddt"] if ooi_conf else None,
                        "delta_iplddt":       delta(ooi_conf["iplddt"] if ooi_conf else None,
                                                    ref_conf["iplddt"] if ref_conf else None),
                        "pde":                ooi_conf["pde"]    if ooi_conf else None,
                        "delta_pde":          delta(ooi_conf["pde"] if ooi_conf else None,
                                                    ref_conf["pde"] if ref_conf else None),
                        "ipde":               ooi_conf["ipde"]   if ooi_conf else None,
                        "delta_ipde":         delta(ooi_conf["ipde"] if ooi_conf else None,
                                                    ref_conf["ipde"] if ref_conf else None),
                        "temperature":        TEMPERATURE_C,
                    }
                    rows.append(ooi_row)

    # Write CSV
    with open(OUTPUT_CSV, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=CSV_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)

    print(f"\nWrote {len(rows)} rows to {OUTPUT_CSV}")


if __name__ == "__main__":
    main()