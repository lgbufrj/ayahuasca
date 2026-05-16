"""
tables.py
─────────────────────────────────────────────────────────────────────────────
Generates paper tables for:
  • Thermodynamics  — per-compound ΔG and per-reaction ΔG
  • BLAST           — best hits for every protein × organism combination

All outputs land in  <PAPER_PATH>/tables/
─────────────────────────────────────────────────────────────────────────────
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd

import data
from logger import build_logger

log = build_logger("tables")

TABLES_PATH = Path(data.PAPER_PATH) / "tables"
TABLES_PATH.mkdir(parents=True, exist_ok=True)

# ──────────────────────────────────────────────────────────────────────────────
# Helpers
# ──────────────────────────────────────────────────────────────────────────────

def _save(df: pd.DataFrame, filename: str) -> None:
    """Write *df* to TABLES_PATH/<filename> and log the result."""
    dest = TABLES_PATH / filename
    df.to_csv(dest, index=False)
    log.info("Saved %d rows → %s", len(df), dest)


def _target_organisms() -> list[str]:
    """Organism keys that are NOT reference genomes."""
    return [org for org, meta in data.organisms.items() if not meta["ref"]]


# ──────────────────────────────────────────────────────────────────────────────
# Thermodynamics
# ──────────────────────────────────────────────────────────────────────────────

def _load_delta_g(compound_key: str) -> float:
    """
    Return the Gibbs free energy (kcal/mol) for *compound_key*, or 0.0 if the
    thermodynamics file is missing.
    """
    cpd        = data.compounds[compound_key]
    pubchem_id = cpd["pubchem_id"]
    path       = (
        Path(data.COMPOUNDS_PATH)
        / compound_key / "thermo"
        / f"{compound_key}_{pubchem_id}.json"
    )

    try:
        with path.open() as fh:
            return json.load(fh).get("g_minus_eel_kcalpmol", 0.0)
    except FileNotFoundError:
        log.warning("Thermo data missing for %s (CID %s)", compound_key, pubchem_id)
        return 0.0


def build_compound_thermo_table() -> pd.DataFrame:
    """One row per compound: name, PubChem ID, ΔG."""
    rows: list[dict[str, Any]] = []

    for key, cpd in data.compounds.items():
        rows.append({
            "Compound":                   cpd["name"],
            "PubChem ID":                 cpd["pubchem_id"],
            "Gibbs Free Energy (kcal/mol)": _load_delta_g(key),
        })

    return pd.DataFrame(rows)


def build_reaction_thermo_table(enzymes: dict) -> pd.DataFrame:
    """One row per reaction: enzyme, equation, calculated ΔG, reference flag."""
    rows: list[dict[str, Any]] = []

    for enzyme_data in enzymes.values():
        for rxn in enzyme_data["reactions"]:
            g_substrates = sum(_load_delta_g(s) for s in rxn["substrates"])
            g_products   = sum(_load_delta_g(p) for p in rxn["products"])
            delta_g      = round(g_products - g_substrates, 2)

            rows.append({
                "Enzyme":                  enzyme_data["abbreviation"],
                "Reaction":                f"{' + '.join(rxn['substrates'])} → {' + '.join(rxn['products'])}",
                "Calculated ΔG (kcal/mol)": delta_g,
                "Reference Reaction":       "Ref" if rxn["ref"] else "",
            })

    return pd.DataFrame(rows)


# ──────────────────────────────────────────────────────────────────────────────
# BLAST
# ──────────────────────────────────────────────────────────────────────────────

def build_blast_table() -> pd.DataFrame:
    """
    Best BLAST hit (first row of the filtered CSV) for every
    (protein × reference-organism × target-organism) combination.
    """
    rows:    list[dict[str, Any]] = []
    targets: list[str]            = _target_organisms()

    for prot_key, prot_data in data.proteins.items():
        ptn_abbrev = prot_data["abbreviation"]

        for ref_org_key, ref_org_data in prot_data["organisms"].items():
            ref_species = data.organisms[ref_org_key]["species"]
            uniprot_id  = ref_org_data["uniprot_id"]

            for organism in targets:
                blast_path = (
                    Path(data.PROTEINS_PATH)
                    / prot_key / "analysis" / "genomic" / "blast"
                    / organism / ref_org_key
                    / f"{prot_key}_{uniprot_id}_blast_filtered.csv"
                )

                try:
                    hit = pd.read_csv(blast_path, sep=";").iloc[0]
                    log.debug(
                        "hit   %s / %s (%s) → %s",
                        ptn_abbrev, ref_species, ref_org_key, organism,
                    )
                    rows.append({
                        "Enzyme":             ptn_abbrev,
                        "Reference Organism": ref_species,
                        "UniProt ID":         uniprot_id,
                        "Organism":           organism,
                        "Hit ID":             hit["hit_id"],
                        "Length":             hit["hit_len"],
                        "E-value":            hit["expect"],
                        "Bit Score":          hit["bit_score"],
                    })

                except FileNotFoundError:
                    log.warning(
                        "miss  %s / %s (%s) → %s — filtered CSV not found",
                        ptn_abbrev, ref_species, ref_org_key, organism,
                    )
                except (IndexError, KeyError) as exc:
                    log.error(
                        "parse error for %s / %s → %s: %s",
                        ptn_abbrev, organism, blast_path.name, exc,
                    )

    return pd.DataFrame(rows)


# ──────────────────────────────────────────────────────────────────────────────
# Entry point
# ──────────────────────────────────────────────────────────────────────────────

def main() -> None:
    log.info("Building thermodynamics tables…")
    _save(build_compound_thermo_table(),          "thermodynamics_compounds_table.csv")
    _save(build_reaction_thermo_table(data.proteins), "thermodynamics_table.csv")

    log.info("Building BLAST table…")
    _save(build_blast_table(), "blast_table.csv")

    log.info("All tables written to %s", TABLES_PATH)


if __name__ == "__main__":
    main()