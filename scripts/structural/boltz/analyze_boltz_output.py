import csv
import math
import os
from pathlib import Path
import math
import pandas as pd

from data import PROTEINS_PATH, proteins, PAPER_PATH


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

TEMPERATURE_C = 20            # °C  — change as needed
R = 1.9872042                 # cal/(mol·K)
OUTPUT_CSV = f"{PAPER_PATH}/tables/boltz_affinity_results.csv"

MODELS_N = 25                   # top-N most-confident models kept per (organism, id, reaction)
REFERENCE_REACTION_ID = "rr"   # reaction id used as the baseline for delta calculations

# Quality filters applied to each model BEFORE building the final dataframe.
# A model is kept only if it passes ALL of them; missing values drop the row.
CONFIDENCE_MIN = 0.0                 # keep models with confidence >= 0.6
IPLDDT_MIN = 0.25                     # keep models with iplddt >= 0.7
LIGAND_DISTANCE_MAX = 15.0           # keep models with ligand_distance <= 15

CSV_COLUMNS = [
    "protein", "organism", "id", "rxn_id",
    "cofactor", "substrate",
    "model", "model_rank",
    "binding_affinity", "kd", "dG", "IC50",
    "delta_affinity", "delta_kd", "delta_dG", "delta_IC50",
    "confidence", "delta_confidence",
    "plddt", "delta_plddt",
    "iplddt", "delta_iplddt",
    "pde", "delta_pde",
    "ipde", "delta_ipde",
    "ligand_distance", "delta_ligand_distance",
    "temperature",
]


# ---------------------------------------------------------------------------
# Thermodynamic helpers
# ---------------------------------------------------------------------------

def pred_to_kd(pred_value: float) -> float:
    """affinity_pred_value (log10 IC50 in µM) → Kd in mol/L."""
    if pred_value is None or pd.isna(pred_value):
        return None
    return math.pow(10, pred_value) / 1_000_000


def kd_to_dG(kd: float, temperature_c: float) -> float:
    """Kd (mol/L) → ΔG in kcal/mol."""
    if kd is None or pd.isna(kd):
        return None
    T = 273.15 + temperature_c
    return R * T * math.log(kd) / 1000


def pred_to_IC50_uM(pred_value: float) -> float:
    """affinity_pred_value → IC50 in µM (direct from model definition)."""
    if pred_value is None or pd.isna(pred_value):
        return None
    return math.pow(10, pred_value)


# ---------------------------------------------------------------------------
# Cofactor / substrate resolution (unchanged logic)
# ---------------------------------------------------------------------------

PROTON_ALIASES = {"H+", "H⁺", "H(+)", "proton"}

def resolve_cofactor_substrate(rxn: dict, cofactor_list: list[str]) -> tuple[str, str]:
    """
    From rxn['substrates'], find which one is the cofactor (present in the
    protein's cofactor list) and which is the substrate, ignoring protons.
    """
    substrates = [
        s for s in rxn.get("substrates", [])
        if s not in PROTON_ALIASES
    ]
    cofactor_set = set(cofactor_list)
    cofactor = next((s for s in substrates if s in cofactor_set), "")
    substrate = next((s for s in substrates if s not in cofactor_set), "")
    return cofactor, substrate


def delta(query_val, ref_val):
    if query_val is None or ref_val is None:
        return None
    if pd.isna(query_val) or pd.isna(ref_val):
        return None
    return query_val - ref_val

def ratio(query_val, ref_val):
    if query_val is None or ref_val is None:
        return None
    if pd.isna(query_val) or pd.isna(ref_val):
        return None
    return math.log(math.fabs(query_val)) / math.log(math.fabs(ref_val)) if ref_val != 0 else None


# ---------------------------------------------------------------------------
# Loading + deriving affinity columns from the summary CSV
# ---------------------------------------------------------------------------

def load_summary(boltz_summary_path: Path) -> pd.DataFrame:
    if not boltz_summary_path.exists():
        print(f"  [WARN] Missing: {boltz_summary_path}")
        return pd.DataFrame()

    # Check whether the first line is a usable header (some summaries carry a
    # stray format token before it, e.g. "bi0inf0").
    header_only = pd.read_csv(boltz_summary_path, skipinitialspace=True, nrows=0)
    looks_ok = len(header_only.columns) > 1

    if looks_ok:
        df = pd.read_csv(boltz_summary_path, skipinitialspace=True)
    else:
        # Skip everything up to but not including the line that contains the
        # affinity column header, and use that line as the CSV header.
        raw = boltz_summary_path.read_text(encoding="utf-8", errors="replace")
        lines = raw.splitlines()
        header_idx = next(
            (i for i, line in enumerate(lines) if "affinity_affinity_pred_value" in line),
            None,
        )
        if header_idx is None:
            print(f"  [WARN] No valid header: {boltz_summary_path}")
            return pd.DataFrame()
        df = pd.read_csv(boltz_summary_path, skipinitialspace=True, skiprows=header_idx)

    df.columns = [str(c).strip() for c in df.columns]
    for col in df.columns:
        if df[col].dtype == object:
            df[col] = df[col].str.strip()
    return df


def compute_affinity_columns(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["binding_affinity"] = df["affinity_affinity_pred_value"]
    df["kd"] = df["binding_affinity"].apply(pred_to_kd)
    df["dG"] = df["kd"].apply(lambda k: kd_to_dG(k, TEMPERATURE_C))
    df["IC50"] = df["binding_affinity"].apply(pred_to_IC50_uM)
    return df


def top_confidence_models(df_group: pd.DataFrame, n: int) -> pd.DataFrame:
    """Top-n rows of a (organism, id, reaction_id) group, by confidence_confidence_score desc."""
    sub = df_group.dropna(subset=["confidence_confidence_score"])
    sub = sub.sort_values("confidence_confidence_score", ascending=False)
    return sub.head(n).reset_index(drop=True)


def apply_quality_filters(df: pd.DataFrame) -> pd.DataFrame:
    """
    Keep only models passing the quality thresholds (rows with missing values
    in any of the filtered columns are dropped):
      - confidence_confidence_score >= CONFIDENCE_MIN
      - confidence_complex_iplddt   >= IPLDDT_MIN
      - dist_ligand_1_ligand_2      <= LIGAND_DISTANCE_MAX
    Applying this before grouping/top-N selection reduces memory usage
    downstream by discarding low-quality models early.
    """
    sub = df.copy()
    n_before = len(sub)

    min_filters = [
        ("confidence_confidence_score", CONFIDENCE_MIN),
        ("confidence_complex_iplddt", IPLDDT_MIN),
    ]
    for col, threshold in min_filters:
        if col not in sub.columns:
            print(f"  [WARN] Filter column {col!r} not present — skipping")
            continue
        sub = sub[sub[col] >= threshold]

    max_filters = [
        ("dist_ligand_1_ligand_2", LIGAND_DISTANCE_MAX),
    ]
    for col, threshold in max_filters:
        if col not in sub.columns:
            print(f"  [WARN] Filter column {col!r} not present — skipping")
            continue
        sub = sub[sub[col] <= threshold]

    removed = n_before - len(sub)
    if removed:
        print(f"  [INFO] Dropped {removed}/{n_before} models by quality filters")
    return sub


# ---------------------------------------------------------------------------
# Row construction
# ---------------------------------------------------------------------------

DELTA_KEYS = [
    "delta_affinity", "delta_kd", "delta_dG", "delta_IC50",
    "delta_confidence", "delta_plddt", "delta_iplddt",
    "delta_pde", "delta_ipde", "delta_ligand_distance",
]


def make_row(protein_name, organism, oid, rxn_id, cofactor, substrate,
             model_row, model_rank, is_reference, ref_row=None):
    aff = {
        "binding_affinity": model_row["binding_affinity"],
        "kd": model_row["kd"],
        "dG": model_row["dG"],
        "IC50": model_row["IC50"],
    }
    conf = {
        "confidence": model_row["confidence_confidence_score"],
        "plddt": model_row.get("confidence_complex_plddt"),
        "iplddt": model_row.get("confidence_complex_iplddt"),
        "pde": model_row.get("confidence_complex_pde"),
        "ipde": model_row.get("confidence_complex_ipde"),
    }
    ligand_distance = model_row.get("dist_ligand_1_ligand_2")

    if is_reference:
        # This row IS the reference reaction ("rr") -> deltas relative to itself are 0
        deltas = {k: 0.0 for k in DELTA_KEYS}
    elif ref_row is not None:
        ref_aff = {
            "binding_affinity": ref_row["binding_affinity"],
            "kd": ref_row["kd"],
            "dG": ref_row["dG"],
            "IC50": ref_row["IC50"],
        }
        ref_conf = {
            "confidence": ref_row["confidence_confidence_score"],
            "plddt": ref_row.get("confidence_complex_plddt"),
            "iplddt": ref_row.get("confidence_complex_iplddt"),
            "pde": ref_row.get("confidence_complex_pde"),
            "ipde": ref_row.get("confidence_complex_ipde"),
        }
        ref_ligand_distance = ref_row.get("dist_ligand_1_ligand_2")

        deltas = {
            "delta_affinity":        delta(aff["binding_affinity"], ref_aff["binding_affinity"]),
            "delta_kd":              delta(aff["kd"], ref_aff["kd"]),
            "delta_dG":              ratio(aff["dG"], ref_aff["dG"]),
            "delta_IC50":            delta(aff["IC50"], ref_aff["IC50"]),
            "delta_confidence":      delta(conf["confidence"], ref_conf["confidence"]),
            "delta_plddt":           delta(conf["plddt"], ref_conf["plddt"]),
            "delta_iplddt":          delta(conf["iplddt"], ref_conf["iplddt"]),
            "delta_pde":             delta(conf["pde"], ref_conf["pde"]),
            "delta_ipde":            delta(conf["ipde"], ref_conf["ipde"]),
            "delta_ligand_distance": delta(ligand_distance, ref_ligand_distance),
        }
    else:
        # no matching reference model available (missing/short rr group)
        deltas = {k: None for k in DELTA_KEYS}

    return {
        "protein": protein_name,
        "organism": organism,
        "id": oid,
        "rxn_id": rxn_id,
        "cofactor": cofactor,
        "substrate": substrate,
        "model": int(model_row["model_id"]),
        "model_rank": model_rank,
        "binding_affinity": aff["binding_affinity"],
        "kd": aff["kd"],
        "dG": aff["dG"],
        "IC50": aff["IC50"],
        "confidence": conf["confidence"],
        "plddt": conf["plddt"],
        "iplddt": conf["iplddt"],
        "pde": conf["pde"],
        "ipde": conf["ipde"],
        "ligand_distance": ligand_distance,
        "temperature": TEMPERATURE_C,
        **deltas,
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    rows = []

    for protein_name, protein_data in proteins.items():
        boltz_summary_path = Path(
            f"{PROTEINS_PATH}/{protein_name}/analysis/structural/boltz/boltz_summary.csv"
        )
        df = load_summary(boltz_summary_path)
        if df.empty:
            continue
        df = compute_affinity_columns(df)
        df = apply_quality_filters(df)
        if df.empty:
            continue

        cofactor_list = protein_data.get("cofactors", [])
        reactions_by_id = {rxn["id"]: rxn for rxn in protein_data["reactions"]}

        for organism, org_df in df.groupby("organism"):
            for oid, id_df in org_df.groupby("protein_or_gene_id"):

                # top-N confidence models per reaction ("r1", "r2", "r3", "rr", ...)
                top_by_rxn = {
                    rxn_id: top_confidence_models(rxn_df, MODELS_N)
                    for rxn_id, rxn_df in id_df.groupby("reaction_id")
                }

                ref_top = top_by_rxn.get(REFERENCE_REACTION_ID)
                if ref_top is None or ref_top.empty:
                    print(f"  [WARN] No reference reaction '{REFERENCE_REACTION_ID}' "
                          f"for {protein_name}/{organism}/{oid}")

                for rxn_id, top_models in top_by_rxn.items():
                    if top_models.empty:
                        print(f"  [WARN] No confidence data for "
                              f"{protein_name}/{organism}/{oid}/{rxn_id}")
                        continue

                    rxn_info = reactions_by_id.get(rxn_id)
                    if rxn_info is not None:
                        cofactor, substrate = resolve_cofactor_substrate(rxn_info, cofactor_list)
                    else:
                        cofactor, substrate = "", ""

                    is_reference = (rxn_id == REFERENCE_REACTION_ID)

                    for rank, (_, model_row) in enumerate(top_models.iterrows(), start=1):
                        ref_row = None
                        if not is_reference and ref_top is not None and rank - 1 < len(ref_top):
                            ref_row = ref_top.iloc[rank - 1]

                        rows.append(make_row(
                            protein_name, organism, oid, rxn_id, cofactor, substrate,
                            model_row, rank, is_reference, ref_row=ref_row,
                        ))

    os.makedirs(Path(OUTPUT_CSV).parent, exist_ok=True)
    with open(OUTPUT_CSV, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=CSV_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)

    print(f"\nWrote {len(rows)} rows to {OUTPUT_CSV}")


if __name__ == "__main__":
    main()