"""
vina.py
─────────────────────────────────────────────────────────────────────────────
Prepares receptor/ligand with Meeko and runs AutoDock Vina docking.

Pipeline (with cofactor)
  1. Receptor    : .pdb  → mk_prepare_receptor → .pdbqt       (skip if exists)
  2. Cofactor    : .sdf  → mk_prepare_ligand   → .pdbqt       (skip if exists)
  3. Cofactor dock: receptor + cofactor         → <prot>_<cof>.pdbqt  (skip if exists)
  4. Complex     : receptor.pdbqt + best cofactor pose → complex .pdbqt
  5. Per reaction: .sdf → mk_prepare_ligand → .pdbqt          (skip if exists)
                   complex + substrate        → <prot>_<cof>_<sub>.pdbqt

Pipeline (without cofactor)
  1. Receptor    : .pdb  → mk_prepare_receptor → .pdbqt       (skip if exists)
  2. Per reaction: .sdf → mk_prepare_ligand → .pdbqt          (skip if exists)
                   receptor + substrate       → <prot>_<sub>.pdbqt

Usage (module)
    from vina import ProteinDockingRunner, run_protein

Usage (CLI)
    python vina.py --protein asmt --organism arabidopsis
    python vina.py --protein asmt --organism arabidopsis --dry-run
─────────────────────────────────────────────────────────────────────────────
"""

from __future__ import annotations

import argparse
import json
import subprocess, os, sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional
import logging
import data
from logger import build_logger

log = build_logger("vina")

# ──────────────────────────────────────────────────────────────────────────────
# Constants
# ──────────────────────────────────────────────────────────────────────────────

VINA_PATH           = "vina"
MK_PREPARE_RECEPTOR = "mk_prepare_receptor.py"
MK_PREPARE_LIGAND   = "mk_prepare_ligand.py"

DEFAULT_EXHAUSTIVENESS = 8
DEFAULT_NUM_MODES      = 10
DEFAULT_ENERGY_RANGE   = 3

# ──────────────────────────────────────────────────────────────────────────────
# Data models
# ──────────────────────────────────────────────────────────────────────────────

@dataclass
class ReceptorPaths:
    """Paths for the prepared protein receptor."""
    pdb:   Path
    pdbqt: Path


@dataclass
class LigandPaths:
    """Paths for a prepared ligand (cofactor or substrate)."""
    sdf:   Path
    pdbqt: Path


@dataclass
class DockingRun:
    """A single Vina docking call: one receptor + one ligand → one output."""
    receptor:    Path         # .pdbqt used as receptor input
    ligand:      Path         # .pdbqt used as ligand input
    config_file: Path
    output_file: Path         # named e.g. asmt_sam.pdbqt or asmt_sam_harmol.pdbqt
    log_file:    Path
    grid_file:   Path
    label:       str          # human-readable description for logging

    @property
    def done(self) -> bool:
        return self.output_file.exists()


@dataclass
class ProteinDockingRunner:
    """
    Orchestrates the full docking pipeline for one protein × organism pair,
    handling cofactor pre-docking and per-reaction substrate docking.
    """
    protein:      str
    organism:     str
    ref_organism: bool

    # Resolved once in __post_init__
    protein_data: dict       = field(init=False)
    prot_base:    Path       = field(init=False)
    vina_base:    Path       = field(init=False)
    struct_base:  Path       = field(init=False)
    grid_file:    Path       = field(init=False)
    receptor:     ReceptorPaths = field(init=False)
    cofactors:    list[str]  = field(init=False)

    def __post_init__(self) -> None:
        self.protein_data = data.proteins[self.protein]
        uniprot_id        = self.protein_data["organisms"][self.organism]["uniprot_id"] if self.ref_organism else ""

        self.prot_base   = Path(data.PROTEINS_PATH) / self.protein
        self.vina_base   = self.prot_base / "analysis" / "structural" / "vina" / self.organism
        self.grid_file   = self.prot_base / "analysis" / "structural" / "p2rank" / self.organism / "docking_grid.json"
        self.cofactors   = self.protein_data.get("cofactors", [])

        if self.ref_organism:
            self.struct_base = self.prot_base / "reference" / self.organism / "structure"
        else:
            self.struct_base = self.prot_base / "analysis" / "structural" / "alphafold" / self.organism

        self.receptor = ReceptorPaths(
            pdb   = self.struct_base / f"{self.protein}_{uniprot_id}.pdb" if self.ref_organism else self.struct_base / "best_model.pdb",
            pdbqt = self.struct_base / f"{self.protein}_{uniprot_id}.pdbqt" if self.ref_organism else self.struct_base / "best_model.pdbqt",
        )

    @property
    def label(self) -> str:
        return f"{self.protein}/{self.organism}"

    def _ligand_paths(self, compound: str) -> LigandPaths:
        pubchem_id = data.compounds[compound]["pubchem_id"]
        cpd_base   = Path(data.COMPOUNDS_PATH) / compound / "structure"
        return LigandPaths(
            sdf   = cpd_base / f"{compound}_{pubchem_id}.sdf",
            pdbqt = cpd_base / f"{compound}_{pubchem_id}.pdbqt",
        )

    def _docking_run(self, receptor: Path, ligand_key: str, *name_parts: str) -> DockingRun:
        """Build a DockingRun with informative output filename."""
        stem        = "_".join(name_parts)          # e.g. asmt_sam or asmt_sam_harmol
        run_dir     = self.vina_base / stem
        ligand      = self._ligand_paths(ligand_key)
        return DockingRun(
            receptor    = receptor,
            ligand      = ligand.pdbqt,
            config_file = run_dir / "vina_config.txt",
            output_file = run_dir / f"{stem}.pdbqt",
            log_file    = run_dir / f"{stem}.log",
            grid_file   = self.grid_file,
            label       = f"{self.label} | {stem}",
        )

    def cofactor_dock(self) -> Optional[DockingRun]:
        """Return the cofactor docking run, or None if no cofactor."""
        if not self.cofactors:
            return None
        # Currently supports a single cofactor; extend list handling if needed
        cof = self.cofactors[0]
        return self._docking_run(self.receptor.pdbqt, cof, self.protein, cof)

    def complex_pdbqt(self) -> Optional[Path]:
        """Path for the receptor+cofactor complex .pdbqt."""
        if not self.cofactors:
            return None
        cof  = self.cofactors[0]
        stem = f"{self.protein}_{cof}"
        return self.vina_base / stem / f"{stem}_complex.pdbqt"

    def substrate_docks(self, receptor: Path) -> list[DockingRun]:
        """
        Return one DockingRun per unique non-cofactor substrate across all reactions.
        receptor is either the bare receptor or the complex .pdbqt.
        """
        seen: set[str] = set()
        runs: list[DockingRun] = []

        for rxn in self.protein_data["reactions"]:
            for substrate in rxn["substrates"]:
                if substrate in self.cofactors or substrate in seen:
                    continue
                seen.add(substrate)

                name_parts = (
                    [self.protein, self.cofactors[0], substrate]
                    if self.cofactors else
                    [self.protein, substrate]
                )
                runs.append(self._docking_run(receptor, substrate, *name_parts))

        return runs


# ──────────────────────────────────────────────────────────────────────────────
# Meeko preparation
# ──────────────────────────────────────────────────────────────────────────────

def prepare_receptor(runner: ProteinDockingRunner, *, dry_run: bool = False) -> None:
    """Convert receptor .pdb → .pdbqt. Skipped if .pdbqt already exists."""
    rec = runner.receptor

    if rec.pdbqt.exists():
        log.info("Receptor already prepared, skipping: %s", rec.pdbqt.name)
        return

    if not rec.pdb.exists():
        raise FileNotFoundError(f"Receptor PDB not found: {rec.pdb}")

    cmd = [
        MK_PREPARE_RECEPTOR,
        "--read_pdb", str(rec.pdb),
        "-o", str(rec.pdbqt).removesuffix(".pdbqt"),
        "-p",
    ]

    log.info("Preparing receptor: %s", rec.pdb.name)
    log.debug("Command: %s", " ".join(cmd))

    if dry_run:
        log.info("dry-run — skipping mk_prepare_receptor")
        return

    _run_cmd(cmd, label=f"mk_prepare_receptor ({runner.label})")
    log.info("Receptor prepared → %s", rec.pdbqt)


def prepare_ligand(ligand: LigandPaths, *, dry_run: bool = False) -> None:
    """Convert ligand .sdf → .pdbqt. Skipped if .pdbqt already exists."""
    if ligand.pdbqt.exists():
        log.info("Ligand already prepared, skipping: %s", ligand.pdbqt.name)
        return

    if not ligand.sdf.exists():
        raise FileNotFoundError(f"Ligand SDF not found: {ligand.sdf}")

    cmd = [
        MK_PREPARE_LIGAND,
        "-i", str(ligand.sdf),
        "-o", str(ligand.pdbqt),
    ]

    log.info("Preparing ligand: %s", ligand.sdf.name)
    log.debug("Command: %s", " ".join(cmd))

    if dry_run:
        log.info("dry-run — skipping mk_prepare_ligand")
        return

    _run_cmd(cmd, label=f"mk_prepare_ligand ({ligand.sdf.stem})")
    log.info("Ligand prepared → %s", ligand.pdbqt)


# ──────────────────────────────────────────────────────────────────────────────
# Complex building
# ──────────────────────────────────────────────────────────────────────────────

def build_complex(
    receptor_pdbqt: Path,
    cofactor_pose:  Path,
    complex_out:    Path,
    *,
    dry_run: bool = False,
) -> None:
    """
    Merge receptor .pdbqt and the best cofactor pose (first MODEL block from
    Vina output) into a single complex .pdbqt for use as the docking receptor.
    Skipped if complex already exists.
    """
    if complex_out.exists():
        log.info("Complex already built, skipping: %s", complex_out.name)
        return

    log.info("Building complex: %s + %s → %s",
             receptor_pdbqt.name, cofactor_pose.name, complex_out.name)

    if dry_run:
        log.info("dry-run — skipping complex build")
        return

    # Receptor: extract only ATOM/HETATM lines. Vina receptor format accepts no
    # other tags — REMARK, TORSDOF, ROOT, BRANCH, CONECT, etc. are ligand-only
    # records and cause "Unknown or inappropriate tag" errors.
    RECEPTOR_TAGS = {"ATOM", "HETATM"}
    receptor_lines = [
        line for line in receptor_pdbqt.read_text().splitlines(keepends=True)
        if line[:6].strip() in RECEPTOR_TAGS
    ]

    # Cofactor pose: extract only ATOM/HETATM lines from the first MODEL block.
    # Vina receptor format accepts no other tags — REMARK, TORSDOF, ROOT, BRANCH
    # etc. are ligand-only records and cause "Unknown or inappropriate tag" errors.
    pose_lines: list[str] = []
    inside = False
    for line in cofactor_pose.read_text().splitlines(keepends=True):
        if line.startswith("MODEL"):
            inside = True
            continue
        if line.startswith("ENDMDL"):
            break
        if inside and line[:6].strip() in RECEPTOR_TAGS:
            pose_lines.append(line)

    if not pose_lines:
        raise RuntimeError(f"No ATOM/HETATM lines found in cofactor pose: {cofactor_pose}")

    complex_out.parent.mkdir(parents=True, exist_ok=True)
    # Note: Do NOT add END record. Vina receptor files must contain only ATOM/HETATM lines.
    # The END tag would cause parse errors when this complex is used as a receptor.
    complex_out.write_text("".join(receptor_lines + pose_lines))
    log.info("Complex written → %s  (%d receptor + %d cofactor atoms)",
             complex_out, len(receptor_lines), len(pose_lines))


# ──────────────────────────────────────────────────────────────────────────────
# Vina config + docking
# ──────────────────────────────────────────────────────────────────────────────

def generate_config(
    run:            DockingRun,
    exhaustiveness: int = DEFAULT_EXHAUSTIVENESS,
    num_modes:      int = DEFAULT_NUM_MODES,
    energy_range:   int = DEFAULT_ENERGY_RANGE,
) -> None:
    """Write a Vina config file from the docking-grid JSON."""
    if not run.grid_file.exists():
        raise FileNotFoundError(f"Docking grid not found: {run.grid_file}")

    with run.grid_file.open() as fh:
        grid = json.load(fh)

    cx, cy, cz = (float(v) for v in grid["center"])
    sx, sy, sz = (float(v) for v in grid["size"])

    config_text = (
        f"receptor = {run.receptor}\n"
        f"ligand   = {run.ligand}\n"
        "\n"
        f"center_x = {cx}\n"
        f"center_y = {cy}\n"
        f"center_z = {cz}\n"
        "\n"
        f"size_x = {sx}\n"
        f"size_y = {sy}\n"
        f"size_z = {sz}\n"
        "\n"
        f"exhaustiveness = {exhaustiveness}\n"
        f"num_modes      = {num_modes}\n"
        f"energy_range   = {energy_range}\n"
    )

    run.config_file.parent.mkdir(parents=True, exist_ok=True)
    run.config_file.write_text(config_text)
    log.info("Config written → %s", run.config_file)


def run_docking(
    run:     DockingRun,
    *,
    dry_run: bool = False,
    exhaustiveness: int = DEFAULT_EXHAUSTIVENESS,
    num_modes:      int = DEFAULT_NUM_MODES,
    energy_range:   int = DEFAULT_ENERGY_RANGE,
) -> subprocess.CompletedProcess | None:
    """
    Generate config and run Vina for a single DockingRun.
    Skipped if output already exists.
    """
    if run.done:
        log.info("Already docked, skipping: %s", run.output_file.name)
        return None

    generate_config(run, exhaustiveness=exhaustiveness,
                    num_modes=num_modes, energy_range=energy_range)

    cmd = [
        VINA_PATH,
        "--config", str(run.config_file),
        "--out",    str(run.output_file),
        "--log",    str(run.log_file),
    ]

    log.info("Docking  %s", run.label)
    log.debug("Command: %s", " ".join(cmd))

    if dry_run:
        log.info("dry-run — skipping Vina subprocess")
        return None

    run.output_file.parent.mkdir(parents=True, exist_ok=True)
    result = _run_cmd(cmd, label=f"vina ({run.label})")
    log.info("Docking complete → %s", run.output_file)
    return result


# ──────────────────────────────────────────────────────────────────────────────
# Orchestration
# ──────────────────────────────────────────────────────────────────────────────

def run_protein(
    runner:  ProteinDockingRunner,
    *,
    dry_run: bool = False,
    exhaustiveness: int = DEFAULT_EXHAUSTIVENESS,
    num_modes:      int = DEFAULT_NUM_MODES,
    energy_range:   int = DEFAULT_ENERGY_RANGE,
) -> None:
    """
    Full pipeline for one protein × organism pair.

    With cofactor:
        prepare receptor → prepare cofactor → dock cofactor → build complex
        → for each unique non-cofactor substrate: prepare + dock

    Without cofactor:
        prepare receptor → for each unique substrate: prepare + dock
    """
    vina_kwargs = dict(dry_run=dry_run, exhaustiveness=exhaustiveness,
                       num_modes=num_modes, energy_range=energy_range)

    log.info("═══ Starting pipeline: %s ═══", runner.label)

    # ── 1. Receptor ───────────────────────────────────────────────────────────
    prepare_receptor(runner, dry_run=dry_run)

    # ── 2. Cofactor branch ────────────────────────────────────────────────────
    if runner.cofactors:
        cof_key   = runner.cofactors[0]
        cof_paths = runner._ligand_paths(cof_key)
        prepare_ligand(cof_paths, dry_run=dry_run)

        cof_dock    = runner.cofactor_dock()
        complex_out = runner.complex_pdbqt()

        run_docking(cof_dock, **vina_kwargs)
        build_complex(runner.receptor.pdbqt, cof_dock.output_file,
                      complex_out, dry_run=dry_run)

        substrate_receptor = complex_out
        log.info("Using complex as receptor for substrate docking: %s", complex_out.name)

    else:
        substrate_receptor = runner.receptor.pdbqt
        log.info("No cofactor — docking substrates against bare receptor")

    # ── 3. Substrate docking (one per unique non-cofactor substrate) ──────────
    substrate_runs = runner.substrate_docks(substrate_receptor)
    log.info("Substrate runs planned: %d", len(substrate_runs))

    for sub_run in substrate_runs:
        sub_key   = sub_run.ligand.stem.rsplit("_", 1)[0]   # reverse pubchem suffix
        # Re-derive the compound key from the ligand path name part
        sub_key   = next(
            k for k in data.compounds
            if sub_run.ligand.name.startswith(k)
        )
        sub_paths = runner._ligand_paths(sub_key)
        prepare_ligand(sub_paths, dry_run=dry_run)
        run_docking(sub_run, **vina_kwargs)

    log.info("═══ Pipeline complete: %s ═══", runner.label)


# ──────────────────────────────────────────────────────────────────────────────
# Internal helpers
# ──────────────────────────────────────────────────────────────────────────────

def _run_cmd(cmd: list[str], *, label: str) -> subprocess.CompletedProcess:
    """Run *cmd*, routing stdout/stderr to DEBUG logs. Raises on non-zero exit."""
    result = subprocess.run(cmd, capture_output=True, text=True)

    for line in result.stdout.splitlines():
        log.debug("[%s stdout] %s", label, line)
    for line in result.stderr.splitlines():
        log.debug("[%s stderr] %s", label, line)

    if result.returncode != 0:
        log.error("%s failed (exit %d):\n%s", label, result.returncode, result.stderr.strip())
        raise RuntimeError(f"{label} failed — see logs for details")

    return result


# ──────────────────────────────────────────────────────────────────────────────
# CLI
# ──────────────────────────────────────────────────────────────────────────────

def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Cofactor-aware Meeko + AutoDock Vina docking pipeline.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--protein",      required=True, help="Protein key (e.g. asmt)")
    p.add_argument("--organism",     required=True, help="Organism key (e.g. arabidopsis)")
    p.add_argument("--ref-organism", action="store_true",
                   help="Use reference structure path instead of AlphaFold")
    p.add_argument("--exhaustiveness", type=int, default=DEFAULT_EXHAUSTIVENESS)
    p.add_argument("--num-modes",      type=int, default=DEFAULT_NUM_MODES)
    p.add_argument("--energy-range",   type=int, default=DEFAULT_ENERGY_RANGE)
    p.add_argument("--dry-run", "-n",  action="store_true",
                   help="Prepare inputs and configs but skip all subprocess calls")
    p.add_argument("--verbose", "-v",  action="store_true")
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)

    if args.verbose:
        log.setLevel(logging.DEBUG)

    runner = ProteinDockingRunner(
        protein      = args.protein,
        organism     = args.organism,
        ref_organism = args.ref_organism,
    )

    try:
        run_protein(
            runner,
            dry_run        = args.dry_run,
            exhaustiveness = args.exhaustiveness,
            num_modes      = args.num_modes,
            energy_range   = args.energy_range,
        )
    except (FileNotFoundError, RuntimeError) as exc:
        log.error("%s", exc)
        return 1

    return 0


if __name__ == "__main__":
    # ── Quick test ────────────────────────────────────────────────────────────
    # Hardcoded runner for manual testing before the data.proteins loop is wired up.
    # Swap for sys.exit(main()) once the loop is in.
    
    for protein, protein_data in data.proteins.items():
    
        ptn_pockets_dir = f"{data.PROTEINS_PATH}/{protein}/analysis/structural/p2rank/"

        if not os.path.exists(ptn_pockets_dir): continue

        for organism in os.listdir(ptn_pockets_dir):
                        
            runner = ProteinDockingRunner(
                protein      = protein,
                organism     = organism,
                ref_organism = True if organism in protein_data['organisms'] else False,
            )
            run_protein(runner)

    # ── CLI (uncomment to switch) ─────────────────────────────────────────────
    # sys.exit(main())