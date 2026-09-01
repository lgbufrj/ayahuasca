#!/usr/bin/env python3
"""
generate_boltz_csv.py

Walks a directory of Boltz output folders named like:

    boltz_results_{organism}_{enzyme_name}_{protein_or_gene_id}_{reaction_id}

each containing a `predictions/` folder with (at least):

    affinity_{organism}_{enzyme_name}_{protein_or_gene_id}_{reaction_id}.json
    confidence_{...}_model_{N}.json
    {...}_model_{N}.cif

and produces a single CSV with one row per (organism, enzyme, protein/gene id,
reaction, model) combination, including:

  - identifiers parsed from the folder name
  - all fields from the affinity_*.json file
  - all fields from the confidence_*_model_*.json file (flattened)
  - ligand info parsed directly out of the CIF's HETATM records:
    per-ligand atom count, centroid (x, y, z), mean B-factor/pLDDT,
    plus pairwise centroid-to-centroid distances between ligands
  - basic polymer info (number of polymer residues) parsed from the CIF

Usage:
    python generate_boltz_csv.py --outputs-dir ./outputs --output-csv boltz_summary.csv
"""

import argparse
import csv
import json
import os
import re
from collections import defaultdict
from itertools import combinations
from math import sqrt

FOLDER_PREFIX = "boltz_results_"

MODEL_RE_JSON = re.compile(r"_model_(\d+)\.json$")
MODEL_RE_CIF = re.compile(r"_model_(\d+)\.cif$")


# --------------------------------------------------------------------------
# Folder name parsing
# --------------------------------------------------------------------------
def parse_folder_name(folder_name):
    """
    boltz_results_{organism}_{enzyme_name}_{protein_or_gene_id}_{reaction_id}

    organism and enzyme_name are assumed to be single underscore-free tokens;
    protein_or_gene_id may itself contain underscores/dots, so we take
    everything between enzyme_name and the final token (reaction_id).
    """
    if not folder_name.startswith(FOLDER_PREFIX):
        return None

    remainder = folder_name[len(FOLDER_PREFIX):]
    parts = remainder.split("_")

    if len(parts) < 4:
        # Not enough parts to safely parse - skip gracefully
        return {
            "organism": parts[0] if len(parts) > 0 else "",
            "enzyme_name": parts[1] if len(parts) > 1 else "",
            "protein_or_gene_id": "_".join(parts[2:-1]) if len(parts) > 3 else "",
            "reaction_id": parts[-1] if len(parts) > 0 else "",
        }

    organism = parts[0]
    enzyme_name = parts[1]
    reaction_id = parts[-1]
    protein_or_gene_id = "_".join(parts[2:-1])

    return {
        "organism": organism,
        "enzyme_name": enzyme_name,
        "protein_or_gene_id": protein_or_gene_id,
        "reaction_id": reaction_id,
    }


# --------------------------------------------------------------------------
# File discovery
# --------------------------------------------------------------------------
def find_predictions_dir(boltz_result_dir):
    """Find the `predictions` folder inside a boltz_results_* directory."""
    direct = os.path.join(boltz_result_dir, "predictions")
    if os.path.isdir(direct):
        return direct
    # fall back: search recursively in case of nesting
    for root, dirs, _ in os.walk(boltz_result_dir):
        if os.path.basename(root) == "predictions":
            return root
    return None


def find_files(predictions_dir):
    """Return (affinity_json_path_or_None, {model_idx: confidence_json_path}, {model_idx: cif_path})."""
    affinity_path = None
    confidence_by_model = {}
    cif_by_model = {}

    for root, _dirs, files in os.walk(predictions_dir):
        for fname in files:
            fpath = os.path.join(root, fname)
            if fname.startswith("affinity_") and fname.endswith(".json"):
                affinity_path = fpath
            elif fname.startswith("confidence_") and fname.endswith(".json"):
                m = MODEL_RE_JSON.search(fname)
                if m:
                    confidence_by_model[int(m.group(1))] = fpath
            elif fname.endswith(".cif"):
                m = MODEL_RE_CIF.search(fname)
                if m:
                    cif_by_model[int(m.group(1))] = fpath

    return affinity_path, confidence_by_model, cif_by_model


# --------------------------------------------------------------------------
# JSON loading / flattening
# --------------------------------------------------------------------------
def load_json(path):
    if not path:
        return {}
    with open(path) as f:
        return json.load(f)


def flatten_dict(d, prefix, out):
    """Flatten nested dicts (e.g. chains_ptm, pair_chains_iptm) into out with prefixed keys."""
    for k, v in d.items():
        key = f"{prefix}_{k}"
        if isinstance(v, dict):
            flatten_dict(v, key, out)
        else:
            out[key] = v


# --------------------------------------------------------------------------
# CIF parsing
# --------------------------------------------------------------------------
def parse_cif(cif_path):
    """
    Parse an mmCIF/ModelCIF file produced by Boltz and return:
        {
          "polymer_residue_count": int,
          "ligands": {
              (auth_asym_id, comp_id): {
                  "num_atoms": int,
                  "centroid": (x, y, z),
                  "mean_bfactor": float,
              },
              ...
          }
        }

    We rely on the fixed atom_site column order used by Boltz/ModelCIF output:
        group_PDB id type_symbol label_atom_id label_alt_id label_comp_id
        label_seq_id auth_seq_id pdbx_PDB_ins_code label_asym_id
        Cartn_x Cartn_y Cartn_z occupancy label_entity_id
        auth_asym_id auth_comp_id B_iso_or_equiv pdbx_PDB_model_num
    """
    ligand_atoms = defaultdict(list)  # (auth_asym_id, comp_id) -> list of (x, y, z, bfactor)
    polymer_residue_ids = defaultdict(set)  # auth_asym_id -> set of auth_seq_id (for polymer chains)

    with open(cif_path) as f:
        for line in f:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            tokens = line.split()
            if len(tokens) < 19:
                continue

            group_pdb = tokens[0]
            comp_id = tokens[5]
            auth_seq_id = tokens[7]
            try:
                x, y, z = float(tokens[10]), float(tokens[11]), float(tokens[12])
                bfactor = float(tokens[17])
            except ValueError:
                continue
            auth_asym_id = tokens[15]

            if group_pdb == "HETATM":
                ligand_atoms[(auth_asym_id, comp_id)].append((x, y, z, bfactor))
            else:
                polymer_residue_ids[auth_asym_id].add(auth_seq_id)

    ligands = {}
    for key, atoms in ligand_atoms.items():
        n = len(atoms)
        cx = sum(a[0] for a in atoms) / n
        cy = sum(a[1] for a in atoms) / n
        cz = sum(a[2] for a in atoms) / n
        mean_b = sum(a[3] for a in atoms) / n
        ligands[key] = {
            "num_atoms": n,
            "centroid": (cx, cy, cz),
            "mean_bfactor": mean_b,
        }

    polymer_residue_count = sum(len(v) for v in polymer_residue_ids.values())

    return {
        "polymer_residue_count": polymer_residue_count,
        "ligands": ligands,
    }


def euclidean(p1, p2):
    return sqrt(sum((a - b) ** 2 for a, b in zip(p1, p2)))


# --------------------------------------------------------------------------
# Row building
# --------------------------------------------------------------------------
def build_rows_for_folder(boltz_result_dir):
    folder_name = os.path.basename(os.path.normpath(boltz_result_dir))
    ids = parse_folder_name(folder_name)
    if ids is None:
        return []

    predictions_dir = find_predictions_dir(boltz_result_dir)
    if predictions_dir is None:
        print(f"  [!] No predictions/ folder found in {boltz_result_dir}, skipping.")
        return []

    affinity_path, confidence_by_model, cif_by_model = find_files(predictions_dir)
    affinity_data = load_json(affinity_path)

    model_indices = sorted(set(confidence_by_model) | set(cif_by_model))
    if not model_indices:
        print(f"  [!] No confidence/cif files found in {predictions_dir}, skipping.")
        return []

    rows = []
    for model_idx in model_indices:
        row = {
            "organism": ids["organism"],
            "enzyme_name": ids["enzyme_name"],
            "protein_or_gene_id": ids["protein_or_gene_id"],
            "reaction_id": ids["reaction_id"],
            "model_id": model_idx,
            "source_folder": folder_name,
        }

        # --- affinity ---
        for k, v in affinity_data.items():
            row[f"affinity_{k}"] = v

        # --- confidence ---
        confidence_path = confidence_by_model.get(model_idx)
        confidence_data = load_json(confidence_path)
        flat_conf = {}
        for k, v in confidence_data.items():
            if isinstance(v, dict):
                flatten_dict(v, k, flat_conf)
            else:
                flat_conf[k] = v
        for k, v in flat_conf.items():
            row[f"confidence_{k}"] = v

        # --- CIF: ligand coords + polymer length ---
        cif_path = cif_by_model.get(model_idx)
        if cif_path:
            cif_data = parse_cif(cif_path)
            row["polymer_residue_count"] = cif_data["polymer_residue_count"]

            # Sort ligands by auth_asym_id for a stable ligand_1, ligand_2, ... ordering
            ligand_items = sorted(cif_data["ligands"].items(), key=lambda kv: kv[0][0])
            centroids = []
            for i, ((auth_asym_id, comp_id), info) in enumerate(ligand_items, start=1):
                prefix = f"ligand_{i}"
                row[f"{prefix}_chain"] = auth_asym_id
                row[f"{prefix}_comp_id"] = comp_id
                row[f"{prefix}_num_atoms"] = info["num_atoms"]
                row[f"{prefix}_centroid_x"] = round(info["centroid"][0], 4)
                row[f"{prefix}_centroid_y"] = round(info["centroid"][1], 4)
                row[f"{prefix}_centroid_z"] = round(info["centroid"][2], 4)
                row[f"{prefix}_mean_bfactor_plddt"] = round(info["mean_bfactor"], 3)
                centroids.append((f"ligand_{i}", info["centroid"]))

            # pairwise ligand-ligand centroid distances
            for (name_a, c_a), (name_b, c_b) in combinations(centroids, 2):
                row[f"dist_{name_a}_{name_b}"] = round(euclidean(c_a, c_b), 4)
        else:
            print(f"  [!] No CIF found for model {model_idx} in {predictions_dir}")

        rows.append(row)

    return rows


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--outputs-dir", default="./outputs", help="Directory containing boltz_results_* folders")
    parser.add_argument("--output-csv", default="boltz_summary.csv", help="Path to write the summary CSV")
    args = parser.parse_args()

    outputs_dir = args.outputs_dir
    if not os.path.isdir(outputs_dir):
        raise SystemExit(f"Outputs directory not found: {outputs_dir}")

    result_dirs = sorted(
        os.path.join(outputs_dir, d)
        for d in os.listdir(outputs_dir)
        if d.startswith(FOLDER_PREFIX) and os.path.isdir(os.path.join(outputs_dir, d))
    )

    if not result_dirs:
        raise SystemExit(f"No '{FOLDER_PREFIX}*' folders found under {outputs_dir}")

    all_rows = []
    for result_dir in result_dirs:
        print(f"Processing {os.path.basename(result_dir)} ...")
        all_rows.extend(build_rows_for_folder(result_dir))

    if not all_rows:
        raise SystemExit("No rows were generated - check folder/file naming against the script's assumptions.")

    # Union of all columns across rows (rows may have different numbers of ligands)
    fieldnames = []
    seen = set()
    # Force a sensible, stable column order for the "core" fields first
    core_order = [
        "organism", "enzyme_name", "protein_or_gene_id", "reaction_id",
        "model_id", "source_folder", "polymer_residue_count",
    ]
    for name in core_order:
        if name not in seen:
            fieldnames.append(name)
            seen.add(name)
    for row in all_rows:
        for key in row.keys():
            if key not in seen:
                fieldnames.append(key)
                seen.add(key)

    with open(args.output_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in all_rows:
            writer.writerow(row)

    print(f"\nWrote {len(all_rows)} row(s), {len(fieldnames)} column(s) to {args.output_csv}")


if __name__ == "__main__":
    main()