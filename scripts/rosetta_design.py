#!/usr/bin/env python3
"""
PyRosetta Enzyme Design Pipeline
=================================
Usage:
    python rosetta_design.py --protein /path/to/alphafold.pdb --ligand /path/to/ligand.mol2
    python rosetta_design.py --protein /path/to/alphafold.pdb --ligand /path/to/ligand.sdf
    python rosetta_design.py --protein /path/to/alphafold.pdb --smiles "CCO"

    # Interactive mode — edit the inputs block and run:
    python rosetta_design.py

What it does:
  1. Auto-installs missing Python packages (biopython, rdkit)
  2. Converts/parameterises your ligand for PyRosetta
  3. Cleans the AlphaFold PDB
  4. Docks the ligand to the protein (baseline)
  5. Computes binding free energy (ddG in Rosetta Energy Units)
  6. Scans single-point mutations at the binding site
  7. Ranks mutants by ddG and writes ranked CSV
  8. Saves top-3 designed complexes as PDB files

Requirements:
  - PyRosetta (conda install -c conda-forge pyrosetta)  ← requires license
  - Everything else is auto-installed
"""

import argparse
import csv
import json
import os
import re
import shutil
import site
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path
from typing import Optional

# ======================================================================
# DEPENDENCY MANAGEMENT — auto-install missing packages
# ======================================================================

_MISSING: list[str] = []


def _import_or_queue(module: str, package: str = "") -> None:
    """Try importing *module*; if missing, queue *package* (or *module*) for install."""
    try:
        __import__(module)
    except ImportError:
        _MISSING.append(package or module)


def ensure_dependencies() -> None:
    """Install queued packages, then enforce PyRosetta is present."""
    _import_or_queue("Bio", "biopython")
    _import_or_queue("rdkit")

    if _MISSING:
        print(f"[deps] Installing missing packages: {', '.join(_MISSING)}")
        subprocess.check_call(
            [sys.executable, "-m", "pip", "install", * _MISSING,
             "--quiet", "--no-warn-script-location"]
        )
        print("[deps] Done  —  please re-run the script.")
        sys.exit(0)

    # PyRosetta cannot be pip-installed; give clear instructions.
    try:
        import pyrosetta   # noqa: F401
    except ImportError:
        print("[ERROR] PyRosetta is not installed.")
        print("  Install it with:  conda install -c conda-forge pyrosetta")
        print("  A license is required (free for academics): "
              "https://www.pyrosetta.org/")
        sys.exit(1)


ensure_dependencies()

# ======================================================================
# IMPORTS
# ======================================================================

from Bio.PDB import PDBParser, PDBIO, Select
from Bio.SeqUtils import seq1
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors

import pyrosetta
from pyrosetta import rosetta

# ======================================================================
# INPUTS — edit here OR use CLI flags
# ======================================================================

PROTEIN_PDB: str = ""
LIGAND_FILE: str = ""
SMILES: str = ""
LIGAND_OUT_NAME: str = "LIG"
BINDING_SITE_RES: list[int] = []
BINDING_SITE_CUTOFF: float = 8.0
ALLOWED_AA: str = "ARNDCQEGHILKMFPSTWYV"
EXCLUDE_WT: bool = True
OUTPUT_DIR: str = "enzyme_design_output"
TOP_N: int = 3
N_DOCK_POSES: int = 5
PYROSETTA_EXTRA: str = "-ex1 -ex2 -use_input_sc -ignore_unrecognized_res"


# ======================================================================
# LIGAND PREPARATION
# ======================================================================

def _find_molfile_to_params() -> Optional[Path]:
    """Locate the molfile_to_params.py script shipped with PyRosetta."""
    import pyrosetta
    root = Path(pyrosetta.__file__).parent
    for candidate in (
        root / "scripts" / "molfile_to_params.py",
        root / "tools" / "molfile_to_params.py",
        root.parent / "tools" / "molfile_to_params.py",
    ):
        if candidate.exists():
            return candidate
    return None


def _smiles_to_mol2(smiles: str, output: Path) -> None:
    """Generate a 3D mol2 file from a SMILES string via RDKit."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")
    mol = Chem.AddHs(mol)
    params = AllChem.EmbedMultipleConfs(mol, numConfs=50, randomSeed=42)
    if params == -1:
        params = AllChem.EmbedMultipleConfs(mol, numConfs=50, randomSeed=42,
                                            useRandomCoords=True)
    if params == -1:
        raise RuntimeError(f"Could not embed 3D conformer for: {smiles}")
    AllChem.UFFOptimizeMolecule(mol)
    Chem.MolToMolFile(mol, str(output.with_suffix(".sdf")))
    # RDKit can write mol2 via the Mol2 writer
    from rdkit.Chem import rdmolfiles
    rdkit_mol2_block = rdmolfiles.MolToMol2Block(mol)
    with open(str(output), "w") as f:
        f.write(rdkit_mol2_block)
    print(f"[lig]  Generated mol2 from SMILES: {output}")


def _convert_to_mol2(input_path: Path, output_mol2: Path) -> None:
    """Convert various ligand formats to MOL2 using RDKit + obabel fallback."""
    # Try RDKit first
    suppl = Chem.SDMolSupplier(str(input_path)) if input_path.suffix.lower() in (".sdf",) else None
    if suppl is not None:
        mols = [m for m in suppl if m is not None]
        if mols:
            mol = mols[0]
            from rdkit.Chem import rdmolfiles
            mol2_block = rdMolfiles.MolToMol2Block(mol)
            with open(str(output_mol2), "w") as f:
                f.write(mol2_block)
            print(f"[lig]  Converted {input_path.name} → {output_mol2.name} (RDKit)")
            return

    # Fallback to obabel
    if shutil.which("obabel"):
        subprocess.run(
            ["obabel", str(input_path), "-omol2", "-O", str(output_mol2)],
            check=True, capture_output=True
        )
        print(f"[lig]  Converted {input_path.name} → {output_mol2.name} (obabel)")
        return

    raise RuntimeError(
        f"Cannot convert {input_path.suffix} to mol2. "
        "Install rdkit or obabel."
    )


def _generate_params(mol2_file: Path, resname: str, working_dir: Path) -> Path:
    """Run molfile_to_params.py to create a .params file for PyRosetta."""
    script = _find_molfile_to_params()
    if script is None:
        raise RuntimeError(
            "Could not find molfile_to_params.py in PyRosetta installation. "
            "Try: conda install -c conda-forge pyrosetta"
        )

    params_out = working_dir / f"{resname}.params"

    subprocess.run(
        [sys.executable, str(script), "-n", resname, "--pdb",
         str(mol2_file), "--directory", str(working_dir)],
        check=True, capture_output=True, text=True,
        cwd=working_dir,
    )
    # molfile_to_params writes <resname>.params and <resname>_ligand.pdb
    if not params_out.exists():
        raise RuntimeError(f"molfile_to_params.py did not create {params_out}")
    print(f"[lig]  Params file: {params_out}")
    return params_out


def prepare_ligand(
    ligand_path: Optional[Path],
    smiles: str,
    working_dir: Path,
    resname: str = "LIG",
) -> tuple[Path, str, Optional[str]]:
    """
    Prepare a ligand for PyRosetta.

    Returns (params_file, resname, smiles_string).
    """
    working_dir.mkdir(parents=True, exist_ok=True)

    if smiles:
        mol2_path = working_dir / f"{resname}.mol2"
        _smiles_to_mol2(smiles, mol2_path)
        smiles_info = smiles
    elif ligand_path:
        mol2_path = working_dir / f"{resname}.mol2"
        ext = ligand_path.suffix.lower()
        if ext == ".mol2":
            shutil.copy2(ligand_path, mol2_path)
        elif ext == ".sdf":
            _convert_to_mol2(ligand_path, mol2_path)
        elif ext == ".pdb":
            # Use RDKit to perceive molecule from PDB
            mol = Chem.MolFromPDBFile(str(ligand_path))
            if mol is None:
                raise ValueError(f"Cannot read ligand PDB: {ligand_path}")
            from rdkit.Chem import rdMolfiles
            mol2_block = rdMolfiles.MolToMol2Block(mol)
            with open(str(mol2_path), "w") as f:
                f.write(mol2_block)
        else:
            _convert_to_mol2(ligand_path, mol2_path)
        smiles_info = _get_smiles_from_mol(mol2_path)
    else:
        raise ValueError("Provide either --ligand or --smiles")

    params_file = _generate_params(mol2_path, resname, working_dir)
    return params_file, resname, smiles_info


def _get_smiles_from_mol(mol2_path: Path) -> str:
    """Extract SMILES from a mol2 file via RDKit."""
    from rdkit.Chem import rdmolfiles
    mol = rdMolfiles.MolFromMol2File(str(mol2_path))
    if mol is None:
        return "unknown"
    return Chem.MolToSmiles(mol)


# ======================================================================
# PROTEIN PREPARATION
# ======================================================================

def clean_alphafold_pdb(input_pdb: Path, output_pdb: Path) -> None:
    """Remove non-standard residues, keep only first model, renumber cleanly."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", str(input_pdb))

    # Keep only the first model
    model = structure[0]

    class ProteinOnly(Select):
        def accept_residue(self, res):
            return res.get_id()[0] == " "  # standard residues only

    io = PDBIO()
    io.set_structure(model)
    io.save(str(output_pdb), select=ProteinOnly())
    print(f"[prot] Cleaned AlphaFold PDB → {output_pdb}")


# ======================================================================
# PYROSETTA SETUP
# ======================================================================

def setup_pyrosetta(params_file: Path, extra: str = "") -> None:
    opts = f"-extra_res_fa {params_file} {extra}".strip()
    if not pyrosetta._initialized:
        pyrosetta.init(extra_options=opts)
        print(f"[pyro] Initialised  |  options: {opts}")


def create_scorefxn() -> pyrosetta.ScoreFunction:
    sfxn = rosetta.core.scoring.get_score_function()
    print(f"[scf]  {sfxn.get_name()}")
    return sfxn


# ======================================================================
# POSE LOADING
# ======================================================================

def load_complex(protein_pdb: Path, ligand_params: Path) -> pyrosetta.Pose:
    pose = pyrosetta.pose_from_file(str(protein_pdb))
    print(f"[pose] {protein_pdb.name}: {pose.total_residue()} residues")
    return pose


def place_ligand_in_pose(
    pose: pyrosetta.Pose,
    params_file: Path,
    resname: str,
    chain: str = "X",
) -> list[int]:
    """Append ligand as a new chain and return its residue indices."""
    from pyrosetta.rosetta.core.chemical import ChemicalManager, ResidueTypeSet
    from pyrosetta.rosetta.core.conformation import ResidueFactory

    rts = ChemicalManager.Instance().get_residue_type_set("fa_standard")
    lig_restype = rts.name_map(resname)

    # Create the ligand residue
    lig_res = ResidueFactory.create_residue(lig_restype)

    # Append to pose via a jump (so it becomes a separate chain)
    pose.append_residue_by_jump(lig_res, pose.total_residue())

    lig_idx = pose.total_residue()
    # Set chain
    pose.pdb_info().chain(lig_idx, chain)

    print(f"[lig]  Placed {resname} at residue {lig_idx} (chain {chain})")
    return [lig_idx]


# ======================================================================
# BINDING SITE DETECTION
# ======================================================================

def get_binding_site(
    pose: pyrosetta.Pose,
    lig_res: list[int],
    cutoff: float,
) -> list[int]:
    """Return protein residue indices within *cutoff* of any ligand atom."""
    bs = set()
    lig_xyz = [
        pose.residue(ri).xyz(ai)
        for ri in lig_res
        for ai in range(1, pose.residue(ri).natoms() + 1)
    ]

    for ri in range(1, pose.total_residue() + 1):
        if ri in lig_res:
            continue
        res = pose.residue(ri)
        for ai in range(1, res.natoms() + 1):
            r_xyz = res.xyz(ai)
            for l_xyz in lig_xyz:
                if r_xyz.distance(l_xyz) < cutoff:
                    bs.add(ri)
                    break
            if ri in bs:
                break

    bs_list = sorted(bs)
    print(f"[bs]   {len(bs_list)} binding-site residues: {bs_list}")
    return bs_list


# ======================================================================
# DOCKING
# ======================================================================

def dock_ligand(
    pose: pyrosetta.Pose,
    lig_res: list[int],
    scorefxn: pyrosetta.ScoreFunction,
    n_poses: int,
) -> float:
    """Multi-trial ligand docking, keep the best-scoring pose."""
    from pyrosetta.rosetta.protocols.ligand_docking import LigandDockProtocol

    best_score = float("inf")
    best_pose = None

    dock = LigandDockProtocol()
    dock.set_ligand_residues(lig_res)
    dock.set_scorefxn(scorefxn)
    dock.set_high_res_scorefxn(scorefxn)

    for trial in range(n_poses):
        trial_pose = pose.clone()
        dock.apply(trial_pose)
        score = scorefxn(trial_pose)
        print(f"  [dock] trial {trial+1}/{n_poses}  score = {score:.3f} REU")
        if score < best_score:
            best_score = score
            best_pose = trial_pose

    if best_pose is not None:
        pose.assign(best_pose)
    print(f"[dock] Best total score: {best_score:.3f} REU")
    return best_score


# ======================================================================
# BINDING ENERGY (ddG)
# ======================================================================

def compute_ddg(
    pose: pyrosetta.Pose,
    lig_res: list[int],
    scorefxn: pyrosetta.ScoreFunction,
) -> float:
    """ddG = E_complex - (E_protein + E_ligand)  in REU."""
    E_complex = scorefxn(pose)

    # Ligand alone
    lig_pose = pyrosetta.Pose()
    lig_pose.append_residue_by_jump(pose.residue(lig_res[0]), 1)
    E_lig = scorefxn(lig_pose)

    # Protein alone
    prot_pose = pose.clone()
    for ri in sorted(lig_res, reverse=True):
        prot_pose.delete_residue_slow(ri)
    E_prot = scorefxn(prot_pose)

    ddG = E_complex - (E_prot + E_lig)
    print(f"[ddG]   {ddG:.3f} REU  (complex={E_complex:.1f}, "
          f"prot={E_prot:.1f}, lig={E_lig:.1f})")
    return ddG


# ======================================================================
# MUTAGENESIS
# ======================================================================

def mutate(pose: pyrosetta.Pose, position: int, aa_1letter: str) -> None:
    aa3 = rosetta.core.chemical.oneletter_to_three(aa_1letter)
    mut = rosetta.protocols.simple_moves.MutateResidue(position, aa3)
    mut.apply(pose)


def relax(pose: pyrosetta.Pose, scorefxn: pyrosetta.ScoreFunction,
          cycles: int = 2) -> None:
    relax_mover = rosetta.protocols.relax.FastRelax(scorefxn, cycles)
    relax_mover.apply(pose)


def scan_mutations(
    pose: pyrosetta.Pose,
    positions: list[int],
    lig_res: list[int],
    allowed_aa: str,
    scorefxn: pyrosetta.ScoreFunction,
    exclude_wt: bool,
) -> list[dict]:
    """Single-point mutation scan: mutate → relax → ddG."""
    results = []

    for pos in positions:
        wt = pose.residue(pos).name1()
        print(f"\n  [scan] Position {pos} (WT = {wt})")
        for aa in allowed_aa:
            if exclude_wt and aa == wt:
                continue
            trial = pose.clone()
            mutate(trial, pos, aa)
            relax(trial, scorefxn, lig_res)
            ddg = compute_ddg(trial, lig_res, scorefxn)
            results.append({
                "position": int(pos),
                "wt": wt,
                "mut": aa,
                "ddG": round(ddg, 3),
            })
            print(f"          {wt}{pos}{aa}  ddG = {ddg:.3f}")

    return results


# ======================================================================
# OUTPUT
# ======================================================================

def write_csv(results: list[dict], path: Path) -> None:
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["Rank", "Position", "WT", "Mutant", "ddG_REU"])
        for i, r in enumerate(results, 1):
            w.writerow([i, r["position"], r["wt"],
                        f"{r['wt']}{r['position']}{r['mut']}",
                        f"{r['ddG']:.3f}"])
    print(f"[out]  {path}")


# ======================================================================
# MAIN
# ======================================================================

def main() -> None:
    parser = argparse.ArgumentParser(
        description="PyRosetta Enzyme Design Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  python rosetta_design.py --protein model.pdb --ligand ligand.mol2\n"
            "  python rosetta_design.py --protein model.pdb --smiles 'CCO'\n"
            "  python rosetta_design.py   (uses INPUTS block)\n"
        ),
    )
    parser.add_argument("--protein", help="AlphaFold PDB file")
    parser.add_argument("--ligand", help="Ligand file (.mol2 / .sdf / .pdb)")
    parser.add_argument("--smiles", help="Ligand SMILES string")
    parser.add_argument("--resname", default="LIG",
                        help="Three-letter residue name for the ligand")
    parser.add_argument("--binding-residues", type=int, nargs="*",
                        help="Binding site residue numbers (auto-detect if empty)")
    parser.add_argument("--cutoff", type=float, default=8.0,
                        help="Binding site distance cutoff (Å)")
    parser.add_argument("--output", default="enzyme_design_output",
                        help="Output directory")
    parser.add_argument("--top-n", type=int, default=3,
                        help="Number of top designs to save")
    parser.add_argument("--dock-trials", type=int, default=5,
                        help="Docking trials per complex")
    args = parser.parse_args()

    # ── Resolve inputs: CLI > config block ──────────────────────────
    protein_pdb   = Path(args.protein or PROTEIN_PDB)
    ligand_file   = Path(args.ligand) if args.ligand else \
                    (Path(LIGAND_FILE) if LIGAND_FILE else None)
    smiles_str    = (args.smiles or SMILES) or ""
    resname       = args.resname or LIGAND_OUT_NAME
    bs_residues   = args.binding_residues or BINDING_SITE_RES
    bs_cutoff     = args.cutoff or BINDING_SITE_CUTOFF
    out_root      = Path(args.output or OUTPUT_DIR)
    top_n         = args.top_n or TOP_N
    dock_trials   = args.dock_trials or N_DOCK_POSES
    allowed_aa    = ALLOWED_AA
    exclude_wt    = EXCLUDE_WT
    pyro_extra    = PYROSETTA_EXTRA

    if not protein_pdb.exists():
        print(f"[ERROR] Protein PDB not found: {protein_pdb}")
        sys.exit(1)

    # ── Working directory ─────────────────────────────────────────────
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    work_dir = out_root / timestamp
    work_dir.mkdir(parents=True, exist_ok=True)
    ligand_dir = work_dir / "ligand"
    ligand_dir.mkdir(exist_ok=True)
    designs_dir = work_dir / "designs"
    designs_dir.mkdir(exist_ok=True)

    print(f"{'='*60}")
    print(f"  PyRosetta Enzyme Design Pipeline")
    print(f"  Output → {work_dir}")
    print(f"{'='*60}")

    # ── 1. Prepare the ligand ─────────────────────────────────────────
    print(f"\n{'─'*60}\n  1  Ligand preparation\n{'─'*60}")
    params_file, resname, resolved_smiles = prepare_ligand(
        ligand_file, smiles_str, ligand_dir, resname,
    )
    # Copy the ligand PDB for reference
    lig_pdb_src = ligand_dir / f"{resname}_ligand.pdb"
    if lig_pdb_src.exists():
        shutil.copy2(lig_pdb_src, work_dir / "ligand_reference.pdb")

    # ── 2. Clean protein structure ────────────────────────────────────
    print(f"\n{'─'*60}\n  2  Protein preparation\n{'─'*60}")
    clean_pdb = work_dir / "protein_clean.pdb"
    clean_alphafold_pdb(protein_pdb, clean_pdb)

    # ── 3. PyRosetta setup ────────────────────────────────────────────
    print(f"\n{'─'*60}\n  3  PyRosetta setup\n{'─'*60}")
    setup_pyrosetta(params_file, pyro_extra)
    scorefxn = create_scorefxn()

    # ── 4. Load complex ───────────────────────────────────────────────
    print(f"\n{'─'*60}\n  4  Load structure\n{'─'*60}")
    pose = load_complex(clean_pdb, params_file)
    lig_res = place_ligand_in_pose(pose, params_file, resname)
    ligand_chain = pose.pdb_info().chain(lig_res[0])

    # ── 5. Binding site ───────────────────────────────────────────────
    print(f"\n{'─'*60}\n  5  Binding site\n{'─'*60}")
    if bs_residues:
        positions = [r for r in bs_residues if r <= pose.total_residue()
                     and r not in lig_res]
        print(f"[bs]   Using user-specified residues: {positions}")
    else:
        positions = get_binding_site(pose, lig_res, bs_cutoff)

    if not positions:
        print("[WARN] No binding site residues detected. Docking skipped.")
        sys.exit(1)

    # ── 6. Baseline docking + ddG ─────────────────────────────────────
    print(f"\n{'─'*60}\n  6  Baseline docking & binding energy\n{'─'*60}")
    dock_ligand(pose, lig_res, scorefxn, dock_trials)
    baseline_ddg = compute_ddg(pose, lig_res, scorefxn)

    baseline_info = {
        "baseline_ddG_REU": round(baseline_ddg, 3),
        "protein_pdb": str(protein_pdb),
        "ligand_smiles": resolved_smiles or "",
        "ligand_params": str(params_file),
        "n_binding_site_residues": len(positions),
        "docking_trials": dock_trials,
        "score_function": scorefxn.get_name(),
    }
    with open(work_dir / "baseline.json", "w") as f:
        json.dump(baseline_info, f, indent=2)

    pose.dump_pdb(str(work_dir / "baseline_docked.pdb"))
    print(f"[save] Baseline docked pose → baseline_docked.pdb")

    # ── 7. Mutation scan ──────────────────────────────────────────────
    print(f"\n{'─'*60}\n  7  Mutation scan\n{'─'*60}")
    print(f"      {len(positions)} positions × {len(allowed_aa)} AAs "
          f"(~{len(positions) * (len(allowed_aa) - (1 if exclude_wt else 0))} "
          f"mutants)")
    all_results = scan_mutations(
        pose, positions, lig_res, allowed_aa, scorefxn, exclude_wt,
    )

    if not all_results:
        print("[WARN] No mutation results.")
        sys.exit(1)

    # ── 8. Rank ───────────────────────────────────────────────────────
    print(f"\n{'─'*60}\n  8  Ranking\n{'─'*60}")
    all_results.sort(key=lambda r: r["ddG"])
    write_csv(all_results, work_dir / "all_mutations_ranked.csv")

    print(f"\n  Top 10 mutants:")
    for rank, r in enumerate(all_results[:10], 1):
        delta = r["ddG"] - baseline_ddg
        print(f"    {rank:2d}. {r['wt']}{r['position']}{r['mut']}  "
              f"ddG = {r['ddG']:.3f}  (ΔΔG = {delta:+.3f})")

    # ── 9. Build top-N designs ────────────────────────────────────────
    print(f"\n{'─'*60}\n  9  Top {top_n} designs\n{'─'*60}")
    for rank in range(min(top_n, len(all_results))):
        r = all_results[rank]
        design = pose.clone()
        mutate(design, r["position"], r["mut"])
        relax(design, scorefxn, lig_res, cycles=3)
        fname = f"rank{rank+1}_{r['wt']}{r['position']}{r['mut']}_ddG{r['ddG']:.2f}.pdb"
        design.dump_pdb(str(designs_dir / fname))
        print(f"  [{rank+1}] {fname}")

    # ── Done ──────────────────────────────────────────────────────────
    print(f"\n{'='*60}")
    print(f"  DONE")
    print(f"  Baseline ddG:    {baseline_ddg:.3f} REU")
    print(f"  Best mutant:     "
          f"{all_results[0]['wt']}{all_results[0]['position']}"
          f"{all_results[0]['mut']}  "
          f"ddG = {all_results[0]['ddG']:.3f} REU")
    print(f"  All results:     {work_dir / 'all_mutations_ranked.csv'}")
    print(f"  Top designs:     {designs_dir}/")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
