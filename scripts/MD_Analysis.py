"""OpenMM MD pipeline for Boltz-predicted strictosidine synthase structures.

Based on the OpenMM quick-start tutorial. Converts a Boltz mmCIF prediction
to PDB (gemmi), solvates with TIP3P water (amber14 force field), minimizes,
equilibrates (NVT then NPT) and runs a short production MD, then analyzes the
trajectory (RMSD, RMSF, radius of gyration) with MDAnalysis.

Example:
    python MD_Analysis.py
    python MD_Analysis.py --input path/to/model.cif --prod-ns 2 --platform CUDA
"""

import argparse
import glob
import os
import sys

import gemmi
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import MDAnalysis as mda
import numpy as np
import pandas as pd

from openmm import LangevinMiddleIntegrator, MonteCarloBarostat, Platform
from openmm.app import (
    DCDReporter,
    ForceField,
    Modeller,
    PDBFile,
    Simulation,
    StateDataReporter,
)
from openmm.unit import atmosphere, bar, femtosecond, kelvin, nanometer, picosecond

DEFAULT_INPUT = (
    "/home/pedro/Desktop/projects/ayahuasca/pathways/harmine_biosynthesis/"
    "proteins/stricto/analysis/structural/boltz/outputs/"
    "boltz_results_tucunaca_stricto_TucunacaG00000060277.2_r1/predictions/"
    "tucunaca_stricto_TucunacaG00000060277.2_r1/"
    "tucunaca_stricto_TucunacaG00000060277.2_r1_model_0.cif"
)

FORCEFIELDS = ("amber14-all.xml", "amber14/tip3pfb.xml")
TEMPERATURE = 300.0 * kelvin
PRESSURE = 1.0 * atmosphere
DT = 2.0 * femtosecond
FRICTION = 1.0 / picosecond
PADDING = 1.0 * nanometer


def resolve_input(input_path, model_index):
    if os.path.isdir(input_path):
        cifs = sorted(glob.glob(os.path.join(input_path, "*.cif")))
        if not cifs:
            sys.exit(f"no .cif files found in {input_path}")
        return cifs[model_index]
    if not os.path.exists(input_path):
        sys.exit(f"input not found: {input_path}")
    return input_path


def cif_to_pdb(cif_path, pdb_path):
    os.makedirs(os.path.dirname(pdb_path), exist_ok=True)
    structure = gemmi.read_structure(cif_path)
    structure.write_pdb(pdb_path)


def build_system(pdb_path):
    pdb = PDBFile(pdb_path)
    forcefield = ForceField(*FORCEFIELDS)
    modeller = Modeller(pdb.topology, pdb.positions)
    modeller.addHydrogens(forcefield)
    modeller.addSolvent(forcefield, model="tip3p", padding=PADDING)
    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=forcefield.PME,
        nonbondedCutoff=1.0 * nanometer,
        constraints=forcefield.HBonds,
    )
    integrator = LangevinMiddleIntegrator(TEMPERATURE, FRICTION, DT)
    platform = Platform.getPlatformByName(args.platform)
    simulation = Simulation(modeller.topology, system, integrator, platform)
    simulation.context.setPositions(modeller.positions)
    simulation.minimizeEnergy()
    return simulation, modeller


def run_phase(simulation, modeller, duration_ps, name, out_prefix, report_every_ps, barostat=False):
    system = simulation.system
    if barostat:
        system.addForce(MonteCarloBarostat(PRESSURE, TEMPERATURE, 25))
        simulation.context.reinitialize(preserveState=True)
    simulation.reporters.append(
        StateDataReporter(
            out_prefix + f"_{name}_energies.csv",
            int(report_every_ps),
            step=True,
            time=True,
            temperature=True,
            potentialEnergy=True,
            kineticEnergy=True,
            totalEnergy=True,
            volume=True,
        )
    )
    simulation.reporters.append(DCDReporter(out_prefix + f"_{name}.dcd", int(report_every_ps)))
    n_steps = int(duration_ps / DT.value_in_unit(picosecond))
    print(f"running {name}: {duration_ps} ps ({n_steps} steps)")
    simulation.step(n_steps)
    simulation.reporters = []


def analyze(out_prefix, traj_name, solvated_pdb, args):
    u = mda.Universe(solvated_pdb, out_prefix + f"_{traj_name}.dcd")
    protein = u.select_atoms("protein")
    ca = u.select_atoms("protein and name CA")
    ref = u.select_atoms("protein and name CA").positions.copy()

    frames = []
    for ts in u.trajectory:
        pos = u.select_atoms("protein and name CA").positions
        r_fit, rmsd_val = mda.lib.algorithms.align.rotation_matrix(pos, ref)
        frames.append((ts.frame, ts.time, rmsd_val))
    rmsd_df = pd.DataFrame(frames, columns=["frame", "time_ps", "rmsd_Ca_nm"])

    rmsf = mda.analysis.rms.RMSF(ca).run()
    rmsf_df = pd.DataFrame(
        {
            "resid": ca.resids,
            "resname": ca.resnames,
            "rmsf_nm": rmsf.results.rmsf,
        }
    )

    rg = mda.analysis.rdf.RadiusOfGyration(protein).run()
    rg_df = pd.DataFrame({"frame": np.arange(len(rg.results.radius_of_gyration)), "rgyr_nm": rg.results.radius_of_gyration})

    fig, axes = plt.subplots(3, 1, figsize=(8, 10), sharex=False)
    axes[0].plot(rmsd_df["time_ps"], rmsd_df["rmsd_Ca_nm"])
    axes[0].set_title("Calpha RMSD vs first frame")
    axes[0].set_xlabel("time (ps)")
    axes[0].set_ylabel("RMSD (nm)")
    axes[1].plot(rmsf_df["resid"], rmsf_df["rmsf_nm"])
    axes[1].set_title("Calpha RMSF")
    axes[1].set_xlabel("residue")
    axes[1].set_ylabel("RMSF (nm)")
    axes[2].plot(rg_df["frame"], rg_df["rgyr_nm"])
    axes[2].set_title("Radius of gyration")
    axes[2].set_xlabel("frame")
    axes[2].set_ylabel("Rg (nm)")
    fig.tight_layout()
    fig.savefig(out_prefix + "_analysis.png", dpi=150)

    rmsd_df.to_csv(out_prefix + "_rmsd.csv", index=False)
    rmsf_df.to_csv(out_prefix + "_rmsf.csv", index=False)
    rg_df.to_csv(out_prefix + "_rgyr.csv", index=False)
    print(f"analysis figures/CSVs written to {os.path.dirname(out_prefix)}/")
    print(f"final mean Ca-RMSD: {rmsd_df['rmsd_Ca_nm'].mean():.3f} nm")
    print(f"residues with highest mobility: {rmsf_df.nlargest(5, 'rmsf_nm')[['resid', 'rmsf_nm']].to_string(index=False)}")


def main(args):
    cif_path = resolve_input(args.input, args.model)
    out_dir = args.output or os.path.dirname(cif_path)
    os.makedirs(out_dir, exist_ok=True)
    base = os.path.splitext(os.path.basename(cif_path))[0]
    out_prefix = os.path.join(out_dir, base)

    pdb_path = out_prefix + ".pdb"
    cif_to_pdb(cif_path, pdb_path)
    print(f"converted {cif_path} -> {pdb_path}")

    simulation, modeller = build_system(pdb_path)
    solvated_pdb = out_prefix + "_solvated.pdb"
    with open(solvated_pdb, "w") as f:
        PDBFile.writeFile(modeller.topology, simulation.context.getState(getPositions=True).getPositions(), f)
    print(f"solvated system written to {solvated_pdb}")

    run_phase(simulation, modeller, args.nvt_ps, "equil_nvt", out_prefix, 1, barostat=False)
    if args.npt_ps > 0:
        run_phase(simulation, modeller, args.npt_ps, "equil_npt", out_prefix, 2, barostat=True)
    run_phase(simulation, modeller, args.prod_ns * 1000, "production", out_prefix, args.report_ps, barostat=True)

    final_pdb = out_prefix + "_final.pdb"
    with open(final_pdb, "w") as f:
        PDBFile.writeFile(modeller.topology, simulation.context.getState(getPositions=True).getPositions(), f)
    print(f"final frame written to {final_pdb}")

    analyze(out_prefix, "production", solvated_pdb, args)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", default=DEFAULT_INPUT,
                        help="path to a Boltz .cif file or a directory of predictions")
    parser.add_argument("--model", type=int, default=0, help="model index to use when --input is a directory")
    parser.add_argument("--output", default=None, help="output directory (default: next to input)")
    parser.add_argument("--nvt-ps", type=float, default=100, help="NVT equilibration length in ps")
    parser.add_argument("--npt-ps", type=float, default=500, help="NPT equilibration length in ps")
    parser.add_argument("--prod-ns", type=float, default=1.0, help="production MD length in ns")
    parser.add_argument("--report-ps", type=float, default=5.0, help="write trajectory/state every N ps")
    parser.add_argument("--platform", default="CPU", help="OpenMM platform: CPU, CUDA, OpenCL, Reference")
    args = parser.parse_args()
    main(args)