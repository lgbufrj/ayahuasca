"""
merge_structural_analysis.py
─────────────────────────────────────────────────────────────────────────────
Merges per-protein docking results with structural alignment data into a
single structural_analysis.csv.

Sources
    docking_results.csv       <PROTEINS_PATH>/<protein>/analysis/structural/vina/
    structural_alignments.csv <PROTEINS_PATH>/<protein>/analysis/structural/alignment/

Output
    structural_analysis.csv   <PROTEINS_PATH>/<protein>/analysis/structural/

Join key
    docking_results.organism  ↔  structural_alignments.ooi

Delta affinity
    delta_affinity = reference_affinity − organism_affinity
    Positive → organism binds more strongly than the reference.
    The reference organism's affinity is looked up from docking_results.csv
    itself, matching on (protein, reference_organism, cofactor, substrate, mode).

The output keeps all docking modes by default; use --best-only to keep only
mode 1 (best affinity) per substrate × organism pair.

Usage
    python merge_structural_analysis.py                   # all proteins
    python merge_structural_analysis.py --protein asmt    # single protein
    python merge_structural_analysis.py --best-only       # mode 1 only
─────────────────────────────────────────────────────────────────────────────
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

import data
from logger import build_logger

log = build_logger("merge_structural")

# ──────────────────────────────────────────────────────────────────────────────
# Paths
# ──────────────────────────────────────────────────────────────────────────────

def _docking_path(protein: str) -> Path:
    return (
        Path(data.PROTEINS_PATH)
        / protein / "analysis" / "structural" / "vina"
        / "docking_results.csv"
    )

def _alignment_path(protein: str) -> Path:
    return (
        Path(data.PROTEINS_PATH)
        / protein / "analysis" / "structural" / "alignment"
        / "structural_alignments.csv"
    )

def _output_path(protein: str) -> Path:
    return (
        Path(data.PROTEINS_PATH)
        / protein / "analysis" / "structural"
        / "structural_analysis.csv"
    )


# ──────────────────────────────────────────────────────────────────────────────
# Delta affinity
# ──────────────────────────────────────────────────────────────────────────────

def _add_delta_affinity(merged: pd.DataFrame) -> pd.DataFrame:
    """
    Add a delta_affinity column: reference_affinity − organism_affinity.

    The reference affinity is looked up from the same DataFrame by matching
    on (protein, reference_organism, cofactor, substrate, mode), where
    reference_organism acts as the organism key for the reference row.
    """
    # Build a lookup table: one affinity value per (protein, organism, cofactor, substrate, mode)
    ref_lookup = (
        merged[["protein", "organism", "cofactor", "substrate", "mode", "affinity"]]
        .drop_duplicates()
        .rename(columns={
            "organism":  "reference_organism",
            "affinity":  "reference_affinity",
        })
    )

    merged = merged.merge(
        ref_lookup,
        on  = ["protein", "reference_organism", "cofactor", "substrate", "mode"],
        how = "left",
    )

    n_missing = merged["reference_affinity"].isna().sum()
    if n_missing:
        log.warning(
            "%d rows could not be matched to a reference affinity "
            "(reference organism may not have been docked)",
            n_missing,
        )

    merged["delta_affinity"] = merged["reference_affinity"] - merged["affinity"]
    merged.drop(columns=["reference_affinity"], inplace=True)

    return merged


# ──────────────────────────────────────────────────────────────────────────────
# Core merge
# ──────────────────────────────────────────────────────────────────────────────

def merge_protein(protein: str, *, best_only: bool = False) -> Path | None:
    """
    Merge docking results with structural alignment for *protein*.
    Returns the output path on success, None if any source file is missing.
    """
    dock_path  = _docking_path(protein)
    align_path = _alignment_path(protein)
    out_path   = _output_path(protein)

    # ── Load ─────────────────────────────────────────────────────────────────
    missing = [p for p in (dock_path, align_path) if not p.exists()]
    if missing:
        for p in missing:
            log.warning("missing  %s", p)
        return None

    try:
        docking   = pd.read_csv(dock_path)
        alignment = pd.read_csv(align_path)
    except pd.errors.EmptyDataError:
        log.warning("Empty data error")
        return None
    except Exception as exc:
        log.warning("Error reading CSV files: %s", exc)
        return None

    log.debug("docking   rows: %d", len(docking))
    log.debug("alignment rows: %d", len(alignment))

    # ── Optionally filter to best mode ────────────────────────────────────────
    if best_only:
        docking = docking[docking["mode"] == 1].copy()
        docking.drop(columns=["rmsd_lb", "rmsd_ub"], inplace=True)
        log.debug("best-only filter → %d rows", len(docking))

    # ── Merge docking + alignment ─────────────────────────────────────────────
    # alignment.ooi is the target organism; docking.organism is the same field.
    merged = docking.merge(
        alignment.rename(columns={"ooi": "organism"}),
        on       = ["protein", "organism"],
        how      = "left",
        validate = "many_to_one",
    )

    n_unmatched = merged["tm_score"].isna().sum()
    if n_unmatched:
        unmatched = merged.loc[merged["tm_score"].isna(), "organism"].unique()
        log.warning(
            "%d rows unmatched after join (organisms with no alignment data): %s",
            n_unmatched, ", ".join(unmatched),
        )

    # ── Delta affinity ────────────────────────────────────────────────────────
    merged = _add_delta_affinity(merged)

    # ── Column order ──────────────────────────────────────────────────────────
    if best_only:
        col_order = [
            "protein", "organism", "reference_organism",
            "cofactor", "substrate",
            "mode", "affinity", "delta_affinity",
            "rmsd", "tm_score",
        ]
    else:
        col_order = [
            "protein", "organism", "reference_organism",
            "cofactor", "substrate",
            "mode", "affinity", "rmsd_lb", "rmsd_ub",
            "delta_affinity",
            "rmsd", "tm_score",
        ]

    extras = [c for c in merged.columns if c not in col_order]
    merged = merged[col_order + extras]

    # ── Write ─────────────────────────────────────────────────────────────────
    out_path.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(out_path, index=False)

    log.info("wrote %d rows → %s", len(merged), out_path)
    return out_path


# ──────────────────────────────────────────────────────────────────────────────
# CLI
# ──────────────────────────────────────────────────────────────────────────────

def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Merge docking results with structural alignment data.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument(
        "--protein", "-p",
        help="Protein key to process (default: all proteins in data.proteins)",
    )
    p.add_argument(
        "--best-only", "-b",
        action="store_true",
        help="Keep only mode 1 (best affinity) per substrate × organism pair",
    )
    p.add_argument("--verbose", "-v", action="store_true")
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    import logging
    args = parse_args(argv)

    if args.verbose:
        log.setLevel(logging.DEBUG)

    targets = [args.protein] if args.protein else list(data.proteins.keys())

    failed = []
    for prot in targets:
        if prot not in data.proteins:
            log.error("Unknown protein key: %s", prot)
            return 1
        log.info("Merging %s…", prot)
        result = merge_protein(prot, best_only=args.best_only)
        if result is None:
            failed.append(prot)

    if failed:
        log.error("Failed (missing source files): %s", ", ".join(failed))
        return 1

    log.info("Done.")
    return 0


if __name__ == "__main__":
    sys.exit(main())