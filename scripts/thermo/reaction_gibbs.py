"""
reaction_gibbs.py
─────────────────────────────────────────────────────────────────────────────
Reaction Gibbs free energies from the per-compound ORCA results collected in
<PAPER_PATH>/tables/<CURRENT_PROJECT>_gibbs.csv.

    ΔG_rxn = Σ G(products) − Σ G(substrates)          (kcal/mol)
    ΔΔG    = ΔG_rxn − ΔG_ref                          (ref: the enzyme's
                                                       reaction with ref: true)

G_kcal already holds the 1 M standard-state shift and the water solvent term.
The proton row (h+_pH<ref>) is shifted to the reaction pH:

    G_H+(pH) = G_H+(pH_ref) − RT·ln(10)·(pH − pH_ref)

Reaction pH, first match wins: rxn["pH"] → rxn["compartment"] →
enzyme["compartment"] (pH from data.compartments) → pH_ref.

Run from scripts/:  python -m thermo.reaction_gibbs [--ph 7.2]
─────────────────────────────────────────────────────────────────────────────
"""

from __future__ import annotations

import argparse
import math
import re
from collections import Counter
from functools import lru_cache
from pathlib import Path
from typing import Any

import pandas as pd

import data

GIBBS_CSV_PATH = Path(data.PAPER_PATH) / "tables" / f"{data.CURRENT_PROJECT}_gibbs.csv"

KCALPMOL_TO_KJPMOL = 4.184
R_KCALPMOLK        = 1.987204259e-3

PROTON = "h+"


def load_compound_gibbs(path: Path = GIBBS_CSV_PATH) -> pd.DataFrame:
    """Per-compound Gibbs table, indexed by compound key (as used in the YAML)."""
    return pd.read_csv(path, dtype={"cid": "string"}).set_index("compound")


def proton_reference_ph(gibbs: pd.DataFrame) -> float:
    """pH the proton row was computed at, parsed from its species name (h+_pH7)."""
    match = re.search(r"_pH([\d.]+)$", gibbs.at[PROTON, "species"])
    if match is None:
        raise ValueError(f"Cannot read reference pH from {gibbs.at[PROTON, 'species']!r}")
    return float(match.group(1))


def species_gibbs(species: str, gibbs: pd.DataFrame, ph: float) -> float:
    """G (kcal/mol) of *species*; the proton is shifted to *ph*."""
    g = gibbs.at[species, "G_kcal"]
    if species == PROTON:
        temperature = gibbs.at[PROTON, "T_K"]
        g -= R_KCALPMOLK * temperature * math.log(10) * (ph - proton_reference_ph(gibbs))
    return g


def reaction_ph(rxn: dict, enzyme_data: dict, default: float) -> float:
    """pH of *rxn*: explicit pH, then reaction/enzyme compartment, then *default*."""
    if rxn.get("pH") is not None:
        return float(rxn["pH"])
    compartment = rxn.get("compartment") or enzyme_data.get("compartment")
    if compartment:
        return data.compartments[compartment]["pH"]
    return default


@lru_cache(maxsize=None)
def species_composition(species: str) -> tuple[Counter, int] | None:
    """(element counts, charge) of the geometry ORCA optimised; None if not found."""
    if species == PROTON:
        return Counter({"H": 1}), 1

    inp = next((Path(data.COMPOUNDS_PATH) / species / "quantum").glob(f"{species}_*_optfreq.inp"), None)
    if inp is None:
        return None
    match = re.search(r"^\s*\*\s*xyz\s+(-?\d+)\s+\d+\s*\n(.*?)\n\s*\*", inp.read_text(), re.S | re.M)
    if match is None:
        return None
    atoms = Counter(line.split()[0] for line in match.group(2).splitlines() if line.strip())
    return atoms, int(match.group(1))


def reaction_imbalance(rxn: dict) -> str:
    """Atoms and charge of products − substrates, e.g. "H+1 q+1"; "" if balanced, "?" if unknown."""
    delta  = Counter()
    charge = 0
    for side, sign in ((rxn["substrates"], -1), (rxn["products"], 1)):
        for s in side:
            comp = species_composition(s)
            if comp is None:
                return "?"
            atoms, q = comp
            for element, n in atoms.items():
                delta[element] += sign * n
            charge += sign * q

    parts = [f"{el}{n:+d}" for el, n in sorted(delta.items()) if n]
    if charge:
        parts.append(f"q{charge:+d}")
    return " ".join(parts)


def _flag(value: Any) -> bool:
    """Normalise YAML truthy flags (True / "true" / "yes" / 1)."""
    if isinstance(value, str):
        return value.strip().lower() in {"true", "yes", "1"}
    return bool(value)


def reaction_delta_g(rxn: dict, gibbs: pd.DataFrame, ph: float) -> float:
    """ΔG (kcal/mol) of *rxn* at *ph*; NaN if any species has no Gibbs energy."""
    species = rxn["substrates"] + rxn["products"]
    if any(s not in gibbs.index or pd.isna(gibbs.at[s, "G_kcal"]) for s in species):
        return float("nan")

    g_substrates = sum(species_gibbs(s, gibbs, ph) for s in rxn["substrates"])
    g_products   = sum(species_gibbs(p, gibbs, ph) for p in rxn["products"])
    return g_products - g_substrates


def _flagged_species(rxn: dict, gibbs: pd.DataFrame) -> list[str]:
    """Species of *rxn* that are missing or whose ORCA status is not 'ok'."""
    flagged = []
    for s in rxn["substrates"] + rxn["products"]:
        if s not in gibbs.index:
            flagged.append(f"{s} (missing)")
        elif gibbs.at[s, "status"] not in {"ok", "analytic"}:
            flagged.append(f"{s} ({gibbs.at[s, 'status']})")
    return flagged


def build_reaction_gibbs_table(
    enzymes: dict, gibbs: pd.DataFrame | None = None, ph: float | None = None
) -> pd.DataFrame:
    """
    One row per reaction: enzyme, equation, pH, ΔG and ΔΔG (kcal/mol and kJ/mol).

    ΔΔG is measured against the enzyme's reference reaction (``rxn["ref"]``);
    it is NaN when the enzyme has no reference reaction. *ph* overrides the
    per-reaction pH for every reaction.
    """
    if gibbs is None:
        gibbs = load_compound_gibbs()
    default_ph = proton_reference_ph(gibbs)

    rows: list[dict[str, Any]] = []

    for enzyme_data in enzymes.values():
        reactions = enzyme_data["reactions"]
        phs       = [
            ph if ph is not None else reaction_ph(rxn, enzyme_data, default_ph)
            for rxn in reactions
        ]
        delta_gs  = [reaction_delta_g(rxn, gibbs, p) for rxn, p in zip(reactions, phs)]
        ref_g     = next(
            (g for rxn, g in zip(reactions, delta_gs) if _flag(rxn.get("ref"))),
            float("nan"),
        )

        for rxn, rxn_ph, delta_g in zip(reactions, phs, delta_gs):
            ddg = delta_g - ref_g

            rows.append({
                "Enzyme":                   enzyme_data["abbreviation"],
                "Rxn ID":                   rxn["id"],
                "Reaction":                 f"{' + '.join(rxn['substrates'])} → {' + '.join(rxn['products'])}",
                "pH":                       rxn_ph,
                "ΔG (kcal/mol)":            round(delta_g, 2),
                "ΔG (kJ/mol)":              round(delta_g * KCALPMOL_TO_KJPMOL, 2),
                "ΔΔG (kcal/mol)":           round(ddg, 2),
                "ΔΔG (kJ/mol)":             round(ddg * KCALPMOL_TO_KJPMOL, 2),
                "Imbalance":                reaction_imbalance(rxn),
                "Flagged species":          "; ".join(_flagged_species(rxn, gibbs)),
            })

    return pd.DataFrame(rows)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Reaction ΔG from per-compound ORCA Gibbs energies.")
    parser.add_argument("--ph", type=float, help="Use this pH for every reaction instead of the compartment pH.")
    args = parser.parse_args()

    pd.set_option("display.width", 200)
    pd.set_option("display.max_columns", None)
    table = build_reaction_gibbs_table(data.proteins, ph=args.ph)
    print(table.to_string(index=False))

    n_nan = table["ΔG (kcal/mol)"].isna().sum()
    if n_nan:
        print(f"\n{n_nan} reaction(s) without ΔG — see 'Flagged species'.")
    n_unbalanced = (table["Imbalance"] != "").sum()
    if n_unbalanced:
        print(f"{n_unbalanced} reaction(s) not balanced in atoms/charge — see 'Imbalance'.")
