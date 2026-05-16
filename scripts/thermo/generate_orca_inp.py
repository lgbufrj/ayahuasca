# Generate ORCA input file from compound structure SDF file
import os
from rdkit import Chem
from data import COMPOUNDS_PATH, compounds

max_iterations = 300
temperature = 298.15  # K
pressure = 1.0  # atm
scaling_factor = 0.985  # for M062X/def2-SVP


def orca_parameters_from_atoms(n_atoms: int):
    """
    ORCA parameters optimized for a 96-core / 1 TB RAM node.
    """
    if n_atoms <= 10:
        return {
            "method": "! M062X D3zero DEF2-TZVP OPT FREQ TightSCF DEFGRID3",
            "nprocs": 75,
            "maxcore": 12000,
        }

    elif n_atoms <= 35:
        return {
            "method": "! M062X D3zero DEF2-SVP OPT FREQ TightSCF DEFGRID2",
            "nprocs": 64,
            "maxcore": 14000,
        }

    elif n_atoms <= 50:
        return {
            "method": "! M062X D3zero DEF2-SVP RIJCOSX OPT FREQ TightSCF DEFGRID2",
            "nprocs": 32,
            "maxcore": 15000,
        }

    elif n_atoms <= 80:
        return {
            "method": "! M062X D3zero DEF2-SVP OPT FREQ TightSCF DEFGRID2",
            "nprocs": 24,
            "maxcore": 6000,
        }

    else:
        return {
            "method": "! M062X D3zero DEF2-SVP OPT FREQ TightSCF DEFGRID2",
            "nprocs": 8,
            "maxcore": 25000,
        }


def generate_orca_inp(compound_key):
    compound = compounds[compound_key]
    pubchem_id = compound["pubchem_id"]

    sdf_path = (
        f"{COMPOUNDS_PATH}/{compound_key}/structure/"
        f"{compound_key}_{pubchem_id}.sdf"
    )

    if not os.path.exists(sdf_path):
        raise FileNotFoundError(f"SDF file not found: {sdf_path}")

    mol = Chem.MolFromMolFile(sdf_path, removeHs=False)
    if mol is None:
        raise ValueError(f"Could not parse SDF: {sdf_path}")

    conf = mol.GetConformer()
    n_atoms = mol.GetNumAtoms()
    charge = Chem.GetFormalCharge(mol)

    params = orca_parameters_from_atoms(n_atoms)

    print(
        f"{compound_key}: {n_atoms} atoms → "
        f"nprocs={params['nprocs']}, "
        f"maxcore={params['maxcore']} MB, "
        f"{'RIJCOSX' if 'RIJCOSX' in params['method'] else 'RI-J'}"
    )

    orca_inp_content = f"""{params['method']}
%pal nprocs {params['nprocs']} end
%maxcore {params['maxcore']}

%geom
   MaxIter {max_iterations}
end

%freq
  Temp {temperature}
  Pressure {pressure}
  ScalFreq {scaling_factor}
end

* xyz {charge} 1
"""

    for atom in mol.GetAtoms():
        pos = conf.GetAtomPosition(atom.GetIdx())
        orca_inp_content += (
            f"{atom.GetSymbol():<3} "
            f"{pos.x:12.6f} {pos.y:12.6f} {pos.z:12.6f}\n"
        )

    orca_inp_content += "*\n"

    orca_inp_path = (
        f"{COMPOUNDS_PATH}/{compound_key}/quantum/"
        f"{compound_key}_{pubchem_id}.inp"
    )

    os.makedirs(os.path.dirname(orca_inp_path), exist_ok=True)
    with open(orca_inp_path, "w") as f:
        f.write(orca_inp_content)

    return orca_inp_path


if __name__ == "__main__":
    for compound_key in compounds:
        # Skip if INP already exists
        inp_path = (f"{COMPOUNDS_PATH}/{compound_key}/quantum/"
                    f"{compound_key}_{compounds[compound_key]['pubchem_id']}.inp")
        if os.path.exists(inp_path):
            print(f"ORCA input file already exists for {compound_key}, skipping...")
            continue
        generate_orca_inp(compound_key)
