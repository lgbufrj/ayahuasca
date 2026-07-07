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

TOTAL_MODELS = 25           # number of structure models Boltz produces per prediction
MODELS_N = 5                # how many of the most-confident models to keep per structure

CSV_COLUMNS = [
    "protein", "organism", "reference_organism", "cofactor", "substrate",
    "model", "model_rank",
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
    conf = data.get("confidence_score")
    if conf is None:
        print(f"  [WARN] confidence_score missing in {path}")
        return None
    return {
        "confidence": conf,
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
# Path builders
# ---------------------------------------------------------------------------

def build_pred_dir(
    boltz_analysis_dir: str,
    protein_name: str,
    uniprot_id: str,
    rxn_id: str,
    ref_organism: str,
    ooi: str | None,
) -> tuple[Path, str]:
    """Returns (predictions_dir, base_name) for either the reference or an OOI."""
    if ooi:
        results_dir = Path(f"{boltz_analysis_dir}/{ooi}/{ref_organism}/{rxn_id}")
        base_name = f"{ooi}_{ref_organism}_{protein_name}_{uniprot_id}_{rxn_id}"
    else:
        results_dir = Path(f"{boltz_analysis_dir}/{ref_organism}/{rxn_id}")
        base_name = f"{ref_organism}_{protein_name}_{uniprot_id}_{rxn_id}"
    return results_dir / f"predictions/{base_name}", base_name


def build_affinity_path(pred_dir: Path, base_name: str) -> Path:
    return pred_dir / f"affinity_{base_name}.json"


def build_confidence_path(pred_dir: Path, base_name: str, model_idx: int) -> Path:
    return pred_dir / f"confidence_{base_name}_model_{model_idx}.json"


# ---------------------------------------------------------------------------
# Top-N confidence model selection
# ---------------------------------------------------------------------------

def get_top_confidence_models(
    pred_dir: Path,
    base_name: str,
    n: int = MODELS_N,
    total_models: int = TOTAL_MODELS,
) -> list[tuple[int, dict]]:
    """
    Reads confidence_*_model_{i}.json for i in [0, total_models), and returns
    the top `n` (model_idx, confidence_dict) pairs sorted by confidence_score
    descending. Models with missing/unreadable confidence files are skipped.
    """
    candidates = []
    for model_idx in range(total_models):
        conf_path = build_confidence_path(pred_dir, base_name, model_idx)
        conf = read_confidence(conf_path)
        if conf is not None:
            candidates.append((model_idx, conf))

    candidates.sort(key=lambda item: item[1]["confidence"], reverse=True)
    return candidates[:n]


def delta(ooi_val, ref_val):
    if ooi_val is None or ref_val is None:
        return None
    return ooi_val - ref_val


# ---------------------------------------------------------------------------
# Row construction
# ---------------------------------------------------------------------------

def make_base_row(protein_name, organism, ref_organism, cofactor, substrate,
                   model_idx, model_rank, aff, conf, ref_aff=None, ref_conf=None):
    """
    Builds one CSV row. If ref_aff/ref_conf are None, this row IS the
    reference (so all deltas are 0 relative to itself). Otherwise deltas
    are computed against the supplied reference values (matched by rank).
    """
    return {
        "protein":            protein_name,
        "organism":           organism,
        "reference_organism": ref_organism,
        "cofactor":           cofactor,
        "substrate":          substrate,
        "model":              model_idx,
        "model_rank":         model_rank,
        "binding_affinity":   aff["binding_affinity"] if aff else None,
        "kd":                 aff["kd"]               if aff else None,
        "dG":                 aff["dG"]               if aff else None,
        "IC50":               aff["IC50"]             if aff else None,
        "delta_affinity":     delta(aff["binding_affinity"] if aff else None,
                                     ref_aff["binding_affinity"] if ref_aff else None) if ref_aff is not None else (0.0 if aff else None),
        "delta_kd":           delta(aff["kd"] if aff else None,
                                     ref_aff["kd"] if ref_aff else None) if ref_aff is not None else (0.0 if aff else None),
        "delta_dG":           delta(aff["dG"] if aff else None,
                                     ref_aff["dG"] if ref_aff else None) if ref_aff is not None else (0.0 if aff else None),
        "delta_IC50":         delta(aff["IC50"] if aff else None,
                                     ref_aff["IC50"] if ref_aff else None) if ref_aff is not None else (0.0 if aff else None),
        "confidence":         conf["confidence"] if conf else None,
        "delta_confidence":   delta(conf["confidence"] if conf else None,
                                     ref_conf["confidence"] if ref_conf else None) if ref_conf is not None else (0.0 if conf else None),
        "plddt":              conf["plddt"]  if conf else None,
        "delta_plddt":        delta(conf["plddt"] if conf else None,
                                     ref_conf["plddt"] if ref_conf else None) if ref_conf is not None else (0.0 if conf else None),
        "iplddt":             conf["iplddt"] if conf else None,
        "delta_iplddt":       delta(conf["iplddt"] if conf else None,
                                     ref_conf["iplddt"] if ref_conf else None) if ref_conf is not None else (0.0 if conf else None),
        "pde":                conf["pde"]    if conf else None,
        "delta_pde":          delta(conf["pde"] if conf else None,
                                     ref_conf["pde"] if ref_conf else None) if ref_conf is not None else (0.0 if conf else None),
        "ipde":               conf["ipde"]   if conf else None,
        "delta_ipde":         delta(conf["ipde"] if conf else None,
                                     ref_conf["ipde"] if ref_conf else None) if ref_conf is not None else (0.0 if conf else None),
        "temperature":        TEMPERATURE_C,
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

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

                # --- reference organism: affinity + top-N confidence models ---
                ref_pred_dir, ref_base_name = build_pred_dir(
                    boltz_analysis_dir, protein_name, uniprot_id,
                    rxn_id, ref_organism, ooi=None,
                )
                ref_aff = read_affinity(build_affinity_path(ref_pred_dir, ref_base_name))
                ref_top_models = get_top_confidence_models(ref_pred_dir, ref_base_name)

                if not ref_top_models:
                    print(f"  [WARN] No confidence files found for {ref_base_name}")

                for rank, (model_idx, conf) in enumerate(ref_top_models, start=1):
                    rows.append(make_base_row(
                        protein_name, ref_organism, ref_organism, cofactor, substrate,
                        model_idx, rank, ref_aff, conf,
                        ref_aff=None, ref_conf=None,  # reference row: deltas are 0
                    ))

                # --- OOI rows: paired by confidence rank against the reference ---
                for ooi in oois:
                    ooi_pred_dir, ooi_base_name = build_pred_dir(
                        boltz_analysis_dir, protein_name, uniprot_id,
                        rxn_id, ref_organism, ooi=ooi,
                    )
                    ooi_aff = read_affinity(build_affinity_path(ooi_pred_dir, ooi_base_name))
                    ooi_top_models = get_top_confidence_models(ooi_pred_dir, ooi_base_name)

                    if not ooi_top_models:
                        print(f"  [WARN] No confidence files found for {ooi_base_name}")

                    n_rows = max(len(ooi_top_models), 1)
                    for rank in range(1, n_rows + 1):
                        ooi_model_idx, ooi_conf = (
                            ooi_top_models[rank - 1] if rank - 1 < len(ooi_top_models)
                            else (None, None)
                        )
                        # pair with the reference's same-rank model, if it exists
                        ref_conf = (
                            ref_top_models[rank - 1][1] if rank - 1 < len(ref_top_models)
                            else None
                        )

                        rows.append(make_base_row(
                            protein_name, ooi, ref_organism, cofactor, substrate,
                            ooi_model_idx, rank, ooi_aff, ooi_conf,
                            ref_aff=ref_aff, ref_conf=ref_conf,
                        ))

    # Write CSV
    with open(OUTPUT_CSV, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=CSV_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)

    print(f"\nWrote {len(rows)} rows to {OUTPUT_CSV}")


if __name__ == "__main__":
    main()