"""
parse_docking_results.py
─────────────────────────────────────────────────────────────────────────────
Parses AutoDock Vina .log files for every protein × organism × substrate
combination and writes a single docking_results.csv per protein to:

    <PROTEINS_PATH>/<protein>/analysis/structural/vina/docking_results.csv

Columns
    protein       protein key       (e.g. asmt)
    organism      organism key      (e.g. arabidopsis)
    cofactor      cofactor key or   (e.g. sam)
    substrate     substrate key     (e.g. nas)
    mode          binding mode rank (1 = best)
    affinity      kcal/mol
    rmsd_lb       RMSD lower bound  (vs best mode)
    rmsd_ub       RMSD upper bound  (vs best mode)

Usage
    python parse_docking_results.py                  # all proteins
    python parse_docking_results.py --protein asmt   # single protein
─────────────────────────────────────────────────────────────────────────────
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterator

import data
from logger import build_logger

log = build_logger("parse_docking")

# ──────────────────────────────────────────────────────────────────────────────
# Output schema
# ──────────────────────────────────────────────────────────────────────────────

FIELDNAMES = [
    "protein", "organism", "cofactor", "substrate",
    "mode", "affinity", "rmsd_lb", "rmsd_ub",
]

# ──────────────────────────────────────────────────────────────────────────────
# Data model
# ──────────────────────────────────────────────────────────────────────────────

@dataclass
class DockingMode:
    mode:     int
    affinity: float   # kcal/mol
    rmsd_lb:  float
    rmsd_ub:  float


@dataclass
class DockingResult:
    protein:   str
    organism:  str
    cofactor:  str        # empty string if no cofactor
    substrate: str
    modes:     list[DockingMode] = field(default_factory=list)

    def as_rows(self) -> list[dict]:
        return [
            {
                "protein":   self.protein,
                "organism":  self.organism,
                "cofactor":  self.cofactor,
                "substrate": self.substrate,
                "mode":      m.mode,
                "affinity":  m.affinity,
                "rmsd_lb":   m.rmsd_lb,
                "rmsd_ub":   m.rmsd_ub,
            }
            for m in self.modes
        ]


# ──────────────────────────────────────────────────────────────────────────────
# Log parsing
# ──────────────────────────────────────────────────────────────────────────────

# Matches lines like:  "   1         -6.8      0.000      0.000"
_MODE_RE = re.compile(
    r"^\s*(\d+)\s+(-?\d+\.?\d*)\s+(\d+\.?\d*)\s+(\d+\.?\d*)\s*$"
)


def parse_log(log_path: Path) -> list[DockingMode]:
    """Extract all mode rows from a Vina .log file."""
    modes: list[DockingMode] = []

    try:
        text = log_path.read_text()
    except FileNotFoundError:
        return modes

    for line in text.splitlines():
        m = _MODE_RE.match(line)
        if m:
            modes.append(DockingMode(
                mode     = int(m.group(1)),
                affinity = float(m.group(2)),
                rmsd_lb  = float(m.group(3)),
                rmsd_ub  = float(m.group(4)),
            ))

    return modes


# ──────────────────────────────────────────────────────────────────────────────
# Job discovery  (mirrors vina.py logic)
# ──────────────────────────────────────────────────────────────────────────────

def _target_organisms() -> list[str]:
    return [org for org, meta in data.organisms.items()]


def _iter_results(protein_key: str) -> Iterator[DockingResult]:
    """
    Yield one DockingResult per (organism × substrate) combination,
    regardless of whether the log file is present.
    """
    prot_data  = data.proteins[protein_key]
    cofactors  = prot_data.get("cofactors", [])
    cof_key    = cofactors[0] if cofactors else ""

    # Collect unique non-cofactor substrates across all reactions
    substrates: list[str] = []
    seen: set[str] = set()
    for rxn in prot_data["reactions"]:
        for sub in rxn["substrates"]:
            if sub not in cofactors and sub not in seen:
                seen.add(sub)
                substrates.append(sub)

    prot_base = Path(data.PROTEINS_PATH) / protein_key

    for organism in _target_organisms():
        vina_base = prot_base / "analysis" / "structural" / "vina" / organism

        for substrate in substrates:
            # Reconstruct the run directory name (mirrors vina.py)
            if cof_key:
                stem    = f"{protein_key}_{cof_key}_{substrate}"
            else:
                stem    = f"{protein_key}_{substrate}"

            log_path = vina_base / stem / f"{stem}.log"

            modes = parse_log(log_path)

            if modes:
                log.debug("parsed %2d modes  %s / %s → %s",
                          len(modes), protein_key, organism, substrate)
            else:
                log.warning("no results  %s / %s → %s  (%s)",
                            protein_key, organism, substrate, log_path)

            yield DockingResult(
                protein   = protein_key,
                organism  = organism,
                cofactor  = cof_key,
                substrate = substrate,
                modes     = modes,
            )


# ──────────────────────────────────────────────────────────────────────────────
# Writing
# ──────────────────────────────────────────────────────────────────────────────

def write_results(protein_key: str) -> Path:
    """Collect all results for *protein_key* and write the CSV."""
    out_dir  = Path(data.PROTEINS_PATH) / protein_key / "analysis" / "structural" / "vina"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "docking_results.csv"

    rows: list[dict] = []
    for result in _iter_results(protein_key):
        rows.extend(result.as_rows())

    with out_path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=FIELDNAMES)
        writer.writeheader()
        writer.writerows(rows)

    log.info("wrote %d rows → %s", len(rows), out_path)
    return out_path


# ──────────────────────────────────────────────────────────────────────────────
# CLI
# ──────────────────────────────────────────────────────────────────────────────

def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Parse Vina docking logs into a single CSV per protein.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument(
        "--protein", "-p",
        help="Protein key to process (default: all proteins in data.proteins)",
    )
    p.add_argument("--verbose", "-v", action="store_true")
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    import logging
    args = parse_args(argv)

    if args.verbose:
        log.setLevel(logging.DEBUG)

    targets = [args.protein] if args.protein else list(data.proteins.keys())

    for prot in targets:
        if prot not in data.proteins:
            log.error("Unknown protein key: %s", prot)
            return 1
        log.info("Processing %s…", prot)
        write_results(prot)

    log.info("Done.")
    return 0


if __name__ == "__main__":
    sys.exit(main())