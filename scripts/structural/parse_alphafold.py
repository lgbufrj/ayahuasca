"""
parse_alphafold.py
──────────────────────────────────────────────────────────────────────────────
Parses AlphaFold/ColabFold zip outputs and prepares a single best structure.

Pipeline
────────
1. unzip *.zip into a folder (if not already extracted)
2. inspect confidence JSON files
3. select best model
4. copy best .cif to parent folder
5. convert .cif → .pdb

Expected AlphaFold/ColabFold contents
────────────────────────────────────
Typical files:
    model_0.cif
    model_1.cif
    ...
    model_0_confidences.json
    model_1_confidences.json

Selection strategy
──────────────────
- Uses ranking_score if present
- Otherwise uses mean pLDDT
- Falls back to model_0

Outputs
───────
example/
    prediction.zip
    prediction/
        model_0.cif
        ...
    best_model.cif
    best_model.pdb

Requirements
────────────
pip install biopython

Optional:
    obabel installed in PATH
──────────────────────────────────────────────────────────────────────────────
"""

from __future__ import annotations

import json
import shutil
import subprocess
import zipfile
import re

from pathlib import Path
from statistics import mean

from Bio.PDB import MMCIFParser
from Bio.PDB import PDBIO


# ──────────────────────────────────────────────────────────────────────────────
# Extraction
# ──────────────────────────────────────────────────────────────────────────────

def unzip_prediction(zip_file: Path) -> Path:
    """
    Extract zip into a same-name folder if not already extracted.

    Returns extraction directory.
    """
    extract_dir = zip_file.with_suffix("")

    if extract_dir.exists():
        print(f"[skip] already extracted: {extract_dir}")
        return extract_dir

    print(f"[extract] {zip_file.name}")

    with zipfile.ZipFile(zip_file, "r") as zf:
        zf.extractall(extract_dir)

    return extract_dir


# ──────────────────────────────────────────────────────────────────────────────
# Confidence parsing
# ──────────────────────────────────────────────────────────────────────────────

def confidence_score(conf_json: Path) -> float:
    """
    Return confidence score from AlphaFold confidence JSON.

    Priority:
        1. ranking_score
        2. mean pLDDT
    """
    with open(conf_json) as f:
        data = json.load(f)

    if "ranking_score" in data:
        return float(data["ranking_score"])

    if "plddt" in data:
        return float(mean(data["plddt"]))

    return -1.0


def find_best_model(prediction_dir: Path) -> Path:

    confidence_files = sorted(
        prediction_dir.glob("*_confidences_*.json")
    )

    # fallback
    if not confidence_files:

        cifs = sorted(
            prediction_dir.glob("*_model_0.cif")
        )

        if cifs:
            print("[warn] no confidence files found, using model_0")
            return cifs[0]

        raise FileNotFoundError(
            f"No confidence files or model_0.cif in {prediction_dir}"
        )

    best_score = float("-inf")
    best_model = None

    for conf in confidence_files:

        score = confidence_score(conf)

        # extract trailing model number
        match = re.search(r'_(\d+)\.json$', conf.name)

        if not match:
            print(f"[warn] could not parse model number: {conf.name}")
            continue

        model_idx = match.group(1)

        cif_candidates = list(
            prediction_dir.glob(f"*_model_{model_idx}.cif")
        )

        if not cif_candidates:
            print(f"[warn] no cif found for model {model_idx}")
            continue

        cif_file = cif_candidates[0]

        print(f"[score] {cif_file.name}: {score:.3f}")

        if score > best_score:
            best_score = score
            best_model = cif_file

    if best_model is None:
        raise RuntimeError(
            f"Could not determine best model in {prediction_dir}"
        )

    print(f"[best] {best_model.name} ({best_score:.3f})")

    return best_model

# ──────────────────────────────────────────────────────────────────────────────
# CIF → PDB conversion
# ──────────────────────────────────────────────────────────────────────────────

def cif_to_pdb(cif_file: Path, pdb_file: Path) -> None:
    """
    Convert mmCIF → PDB using Biopython.
    """
    print(f"[convert] {cif_file.name} → {pdb_file.name}")

    parser = MMCIFParser(QUIET=True)

    structure = parser.get_structure(
        cif_file.stem,
        str(cif_file)
    )

    io = PDBIO()
    io.set_structure(structure)
    io.save(str(pdb_file))

# ──────────────────────────────────────────────────────────────────────────────
# Cleanup
# ──────────────────────────────────────────────────────────────────────────────

def cleanup_extracted(prediction_dir: Path):
    if prediction_dir.exists():
        print(f"[cleanup] removing {prediction_dir}")
        shutil.rmtree(prediction_dir)

# ──────────────────────────────────────────────────────────────────────────────
# Main processing
# ──────────────────────────────────────────────────────────────────────────────

def process_prediction_zip(zip_file: str | Path, cleanup: bool = True) -> tuple[Path, Path]:
    """
    Full pipeline for one AlphaFold zip.

    Returns:
        (best_cif, best_pdb)
    """
    zip_file = Path(zip_file)

    if not zip_file.exists():
        raise FileNotFoundError(zip_file)

    # 1. unzip
    prediction_dir = unzip_prediction(zip_file)

    # 2. select best model
    best_cif_source = find_best_model(prediction_dir)

    # 3. copy best cif to parent folder
    best_cif = zip_file.parent / "best_model.cif"

    shutil.copy2(best_cif_source, best_cif)

    print(f"[copy] best CIF → {best_cif}")

    # 4. convert to pdb
    best_pdb = zip_file.parent / "best_model.pdb"

    cif_to_pdb(best_cif, best_pdb)

    print(f"[done] best PDB → {best_pdb}")
    
    # 5. Cleanup
    if cleanup:
        cleanup_extracted(prediction_dir)
        print(f"[done] Cleaned up.")

    return best_cif, best_pdb


# ──────────────────────────────────────────────────────────────────────────────
# CLI
# ──────────────────────────────────────────────────────────────────────────────

if __name__ == "__main__":

    # import argparse

    # parser = argparse.ArgumentParser(
    #     description="Parse AlphaFold zip outputs"
    # )

    # parser.add_argument(
    #     "zip_files",
    #     nargs="+",
    #     help="AlphaFold zip files"
    # )

    # args = parser.parse_args()

    # for zip_file in args.zip_files:

    zip_files = ["/home/pedro/Desktop/projects/ayahuasca/pathways/harmine_biosynthesis/proteins/asmt/analysis/structural/alphafold/tucunaca/asmt_tucunacag00000034930_1.zip"]

    for zip_file in zip_files:

        print("\n" + "=" * 80)
        print(f"Processing: {zip_file}")
        print("=" * 80)

        try:
            process_prediction_zip(zip_file)

        except Exception as exc:
            print(f"[error] {zip_file}: {exc}")