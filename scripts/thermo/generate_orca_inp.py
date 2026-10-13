#!/usr/bin/env python3
"""
Generate a consistent set of ORCA 6.1.1 inputs for reaction free energies.
This script only writes inputs; the jobs run on a separate machine.

Every species gets the same two-step protocol:
  1. <tag>_optfreq.inp : M06-2X-D3(0)/def2-SVP  SMD(water)  TightOpt + NumFreq
  2. <tag>_sp.inp      : M06-2X-D3(0)/def2-TZVP SMD(water)  single point on
                         the optimized geometry (<tag>_optfreq.xyz)

G(species) = E_el(sp) + [G - E(el)](optfreq) + standard-state corrections

Run each _sp.inp from the same directory, after its _optfreq job finishes.
"""
import argparse
import os
from collections import Counter

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

from data import (COMPOUNDS_PATH, DEFAULT_PRESSURE, DEFAULT_TEMPERATURE,
                  compounds, proteins)

# --------------------------------------------------------------------------
# Level of theory: identical for every species. Never vary these by size.
# --------------------------------------------------------------------------
COMMON = "M062X D3ZERO RIJCOSX def2/J DEFGRID3 TightSCF CPCM(water)"
FREQ_KEYWORD = "NumFreq"  # works for M06-2X in any build; "Freq" if analytic is supported
OPTFREQ_KEYWORDS = f"! {COMMON} def2-SVP TightOpt {FREQ_KEYWORD}"
SP_KEYWORDS = f"! {COMMON} def2-TZVP"

SOLVATION_BLOCK = """%cpcm
  smd true
  SMDsolvent "water"
end
"""

TEMPERATURE = DEFAULT_TEMPERATURE  # K
PRESSURE = DEFAULT_PRESSURE        # atm; 1 M standard state applied afterwards
SCALING_FACTOR = 0.97              # verify for M06-2X/def2-SVP in Truhlar's database
GEOM_MAXITER = 300
SCF_MAXITER = 500
SLOWCONV_ABOVE = 50                # atoms; SlowConv changes the SCF path, not the energy

# --------------------------------------------------------------------------
# Resources for the machine that RUNS ORCA (not the one running this script):
# 2 x Xeon Gold 6438Y+ (2 NUMA nodes x 32 physical cores), 1 TB RAM.
# Each job must fit in one NUMA node, so nprocs never exceeds 32.
# --------------------------------------------------------------------------
CORES_PER_NUMA_NODE = 32
MAXCORE_MB = 8000  # per process; 64 busy cores x 8000 MB = 512 GB


def resources(n_atoms):
    """(nprocs, maxcore MB per process). Several jobs run side by side,
    each bound to one NUMA node by the runner."""
    if n_atoms <= 10:
        nprocs = 4
    elif n_atoms <= 40:
        nprocs = 16
    else:
        nprocs = 32
    return min(nprocs, CORES_PER_NUMA_NODE), MAXCORE_MB


# --------------------------------------------------------------------------
# Species settings. Keys must match the keys used in data.compounds/proteins.
# --------------------------------------------------------------------------
PROTON_KEYS = {"h+", "h_plus", "proton"}  # no calculation; literature value

MULTIPLICITY = {"o2": 3}  # everything else is a singlet

# Truncated NAD(P) models built from SMILES (keys must match data.compounds)
STRUCTURE_OVERRIDES = {
    "nadph": "CN1C=CCC(C(N)=O)=C1",   # 1-methyl-1,4-dihydronicotinamide
    "nadp+": "C[n+]1cccc(C(N)=O)c1",  # 1-methylnicotinamide cation
}
MODEL_SPECIES = set(STRUCTURE_OVERRIDES)  # SMARTS cross-check skipped for these

# pH 7 microspecies (also dominant at pH ~5-6), built from the PubChem SDF
# with stereochemistry kept. Each rule applies to every matching site.
PROTONATE_AMINES = {           # sp3 amines -> ammonium
    "tryptophan", "tryptamine", "serotonin", "pip", "sam", "sah",
    "tetrahydroharman", "tetrahydroharmol", "tetrahydroharmine", "strictosidine",
}
PROTONATE_IMINES = {           # C=N and pyridine-type N -> iminium / pyridinium
    "harmalan", "harmalol", "harmaline", "harman", "harmol", "harmine",
}
DEPROTONATE_ACIDS = {          # COOH -> COO- (remove for reactions below pH ~4.5)
    "tryptophan", "pip", "sah", "p2c", "iaa", "caffeic_acid", "ferulic_acid",
}

AMINE = Chem.MolFromSmarts(
    "[NX3;+0;!$(N-a);!$(N-[#6,#15,#16]=[#7,#8,#16]);!$(N-[#7,#8])]")
IMINE = Chem.MolFromSmarts("[$([NX2;+0]=[#6]),$([nX2;+0;H0;r6])]")
CARBOXYLIC_OH = Chem.MolFromSmarts("[OX2H1;+0][CX3]=O")


# --------------------------------------------------------------------------
# Reactions from data.proteins
# --------------------------------------------------------------------------
def collect_reactions():
    """Unique reactions across all proteins, labeled 'protein:id'.
    Protons in the YAML are ignored; they are inferred from the species."""
    unique = {}
    for ptn, pdata in proteins.items():
        for rxn in pdata.get("reactions", []):
            subs = [k for k in rxn["substrates"] if k not in PROTON_KEYS]
            prods = [k for k in rxn["products"] if k not in PROTON_KEYS]
            key = (tuple(sorted(subs)), tuple(sorted(prods)))
            label = f"{ptn}:{rxn.get('id', '?')}"
            if key in unique:
                unique[key]["labels"].append(label)
            else:
                unique[key] = {"labels": [label], "substrates": subs,
                               "products": prods, "smarts": rxn.get("smarts")}
    return list(unique.values())


def reaction_string(rxn):
    return f"{' + '.join(rxn['substrates'])} -> {' + '.join(rxn['products'])}"


# --------------------------------------------------------------------------
# Structures
# --------------------------------------------------------------------------
def species_tag(key):
    pubchem_id = compounds.get(key, {}).get("pubchem_id", "model")
    return f"{key}_{pubchem_id}"


def read_sdf(key):
    if key not in compounds:
        raise KeyError(f"'{key}' is neither in data.compounds nor in STRUCTURE_OVERRIDES")
    sdf_path = (f"{COMPOUNDS_PATH}/{key}/structure/"
                f"{key}_{compounds[key]['pubchem_id']}.sdf")
    if not os.path.exists(sdf_path):
        raise FileNotFoundError(f"SDF file not found: {sdf_path}")
    mol = Chem.MolFromMolFile(sdf_path, removeHs=False)
    if mol is None:
        raise ValueError(f"Could not parse SDF: {sdf_path}")
    if not mol.GetConformer().Is3D():
        raise ValueError(f"{sdf_path} has 2D coordinates; download the 3D SDF")
    if any(a.GetNumImplicitHs() for a in mol.GetAtoms()):
        raise ValueError(f"{sdf_path} lacks explicit hydrogens")
    return mol


def ph7_form(mol, key):
    """Applies the pH 7 edits for `key`; returns an H-free molecule with stereo kept."""
    Chem.AssignStereochemistryFrom3D(mol)
    mol = Chem.RWMol(Chem.RemoveHs(mol))
    rules = []
    if key in PROTONATE_AMINES:
        rules.append(("amine", AMINE, +1))
    if key in PROTONATE_IMINES:
        rules.append(("imine/pyridine N", IMINE, +1))
    if key in DEPROTONATE_ACIDS:
        rules.append(("carboxylic acid", CARBOXYLIC_OH, -1))
    for name, pattern, delta in rules:
        hits = [m[0] for m in mol.GetSubstructMatches(pattern)]
        if not hits:
            print(f"  note: {key}: no neutral {name} found (already charged in the SDF?)")
        for idx in hits:
            atom = mol.GetAtomWithIdx(idx)
            atom.SetNumExplicitHs(atom.GetTotalNumHs() + delta)
            atom.SetFormalCharge(atom.GetFormalCharge() + delta)
            atom.SetNoImplicit(True)
    Chem.SanitizeMol(mol)
    return mol.GetMol()


def embed(mol, n_confs=50, seed=42):
    """Lowest-MMFF-energy conformer out of n_confs ETKDG embeddings."""
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    conf_ids = list(AllChem.EmbedMultipleConfs(mol, numConfs=n_confs, params=params))
    if not conf_ids:
        raise ValueError(f"Embedding failed for {Chem.MolToSmiles(mol)}")
    results = AllChem.MMFFOptimizeMoleculeConfs(mol, maxIters=2000)
    energies = [e if status != -1 else float("inf") for status, e in results]
    best = conf_ids[min(range(len(conf_ids)), key=energies.__getitem__)]
    return mol, best


def load_species(key):
    if key in STRUCTURE_OVERRIDES:
        mol, conf_id = embed(Chem.MolFromSmiles(STRUCTURE_OVERRIDES[key]))
        source = "model, SMILES"
    else:
        mol = read_sdf(key)
        if key in PROTONATE_AMINES | PROTONATE_IMINES | DEPROTONATE_ACIDS:
            mol, conf_id = embed(ph7_form(mol, key))
            source = "PubChem SDF, pH 7 form"
        else:
            conf_id = mol.GetConformer().GetId()
            source = "PubChem SDF"

    charge = Chem.GetFormalCharge(mol)
    mult = MULTIPLICITY.get(key, 1)
    n_electrons = sum(a.GetAtomicNum() for a in mol.GetAtoms()) - charge
    if (n_electrons + mult - 1) % 2:
        raise ValueError(f"{key}: {n_electrons} electrons incompatible with multiplicity {mult}")

    return {
        "mol": mol,
        "conf_id": conf_id,
        "charge": charge,
        "mult": mult,
        "source": source,
        "formula": Counter(a.GetSymbol() for a in mol.GetAtoms()),
        "smiles": Chem.MolToSmiles(Chem.RemoveHs(mol)),
    }


# --------------------------------------------------------------------------
# Checks
# --------------------------------------------------------------------------
def proton_balance(rxn, species):
    """Infers the protons needed to balance a reaction.

    Returns (n_protons, problem). n_protons > 0 means H+ is consumed
    (substrate side), < 0 means released. problem is None when the
    reaction balances after adding protons, else a description.
    """
    d_atoms, d_charge = Counter(), 0
    for sign, side in ((-1, rxn["substrates"]), (1, rxn["products"])):
        for key in side:
            for el, n in species[key]["formula"].items():
                d_atoms[el] += sign * n
            d_charge += sign * species[key]["charge"]

    n_protons = d_charge  # each H+ on the substrate side cancels +1 charge and +1 H
    d_atoms["H"] -= n_protons
    residual = {el: n for el, n in d_atoms.items() if n}
    if residual:
        return n_protons, (f"atoms off even after {n_protons:+d} H+: "
                           f"products - substrates = {residual}")
    return n_protons, None


def heavy_atoms_from_smiles(side_smiles):
    RDLogger.DisableLog("rdApp.warning")
    try:
        mol = Chem.MolFromSmiles(side_smiles)
    finally:
        RDLogger.EnableLog("rdApp.warning")
    if mol is None:
        return None
    return Counter(a.GetSymbol() for a in mol.GetAtoms() if a.GetAtomicNum() > 1)


def check_smarts(rxn, species):
    """Compares heavy atoms in the SMARTS with the listed species.
    Returns a warning message, or None."""
    smarts = rxn.get("smarts")
    if not smarts or MODEL_SPECIES & set(rxn["substrates"] + rxn["products"]):
        return None
    try:
        lhs, rhs = smarts.split(">>")
    except ValueError:
        return "SMARTS has no '>>'"
    problems = []
    for name, side_smiles, side_keys in (("substrates", lhs, rxn["substrates"]),
                                         ("products", rhs, rxn["products"])):
        from_smarts = heavy_atoms_from_smiles(side_smiles)
        if from_smarts is None:
            problems.append(f"could not parse SMARTS {name}")
            continue
        from_species = Counter()
        for key in side_keys:
            from_species.update({el: n for el, n in species[key]["formula"].items()
                                 if el != "H"})
        if from_smarts != from_species:
            problems.append(f"{name}: SMARTS {dict(from_smarts)} vs "
                            f"species {dict(from_species)}")
    return "; ".join(problems) or None


# --------------------------------------------------------------------------
# ORCA inputs
# --------------------------------------------------------------------------
def write_inputs(key, s, force=False):
    mol = s["mol"]
    conf = mol.GetConformer(s["conf_id"])
    n_atoms = mol.GetNumAtoms()
    charge, mult = s["charge"], s["mult"]
    tag = species_tag(key)

    qdir = os.path.join(COMPOUNDS_PATH, key, "quantum")
    os.makedirs(qdir, exist_ok=True)
    optfreq_path = os.path.join(qdir, f"{tag}_optfreq.inp")
    sp_path = os.path.join(qdir, f"{tag}_sp.inp")

    if not force and os.path.exists(optfreq_path) and os.path.exists(sp_path):
        print(f"  {key}: inputs exist, skipping (use --force to overwrite)")
        return

    extra = []
    if mult > 1:
        extra.append("UKS")
    if n_atoms > SLOWCONV_ABOVE:
        extra.append("SlowConv")
    extra = (" " + " ".join(extra)) if extra else ""

    nprocs, maxcore = resources(n_atoms)
    blocks = (f"%pal nprocs {nprocs} end\n"
              f"%maxcore {maxcore}\n\n"
              f"%scf\n  MaxIter {SCF_MAXITER}\nend\n\n"
              f"{SOLVATION_BLOCK}")

    coords = ""
    for atom in mol.GetAtoms():
        p = conf.GetAtomPosition(atom.GetIdx())
        coords += f"{atom.GetSymbol():<3}{p.x:14.8f}{p.y:14.8f}{p.z:14.8f}\n"

    optfreq = f"""{OPTFREQ_KEYWORDS}{extra}
{blocks}
%geom
  MaxIter {GEOM_MAXITER}
end

%freq
  Temp {TEMPERATURE}
  Pressure {PRESSURE}
  ScalFreq {SCALING_FACTOR}
end

* xyz {charge} {mult}
{coords}*
"""

    single_point = f"""{SP_KEYWORDS}{extra}
{blocks}
* xyzfile {charge} {mult} {tag}_optfreq.xyz
"""

    with open(optfreq_path, "w") as f:
        f.write(optfreq)
    with open(sp_path, "w") as f:
        f.write(single_point)
    print(f"  {key}: wrote {os.path.basename(optfreq_path)}, {os.path.basename(sp_path)} "
          f"({n_atoms} atoms, nprocs={nprocs}, maxcore={maxcore})")


# --------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--force", action="store_true", help="overwrite existing inputs")
    parser.add_argument("--allow-unbalanced", action="store_true",
                        help="write inputs even if some reactions are unbalanced")
    args = parser.parse_args()

    reactions = collect_reactions()
    keys = sorted({k for r in reactions for k in r["substrates"] + r["products"]}
                  - PROTON_KEYS)

    print(f"Conditions: T = {TEMPERATURE} K, P = {PRESSURE} atm\n")
    print("Species:")
    species = {}
    for key in keys:
        species[key] = load_species(key)
        s = species[key]
        print(f"  {key:<20} {s['mol'].GetNumAtoms():>3} atoms  "
              f"charge {s['charge']:+d}  mult {s['mult']}  ({s['source']})\n"
              f"  {'':<20} {s['smiles']}")

    print("\nReaction checks:")
    unbalanced = False
    for rxn in reactions:
        labels = ", ".join(rxn["labels"])
        n_protons, problem = proton_balance(rxn, species)
        rxn["n_protons"] = n_protons  # keep for the ΔG step
        smarts = check_smarts(rxn, species)

        if problem:
            unbalanced = True
            status = "UNBALANCED"
        elif n_protons > 0:
            status = f"OK, {n_protons} H+ consumed"
        elif n_protons < 0:
            status = f"OK, {-n_protons} H+ released"
        else:
            status = "OK"

        print(f"  [{labels}] {reaction_string(rxn)}\n      {status}")
        if problem:
            print(f"      {problem}")
        if smarts:
            print(f"      SMARTS mismatch: {smarts}")

    if unbalanced and not args.allow_unbalanced:
        raise SystemExit("\nFix protonation states or reaction entries above "
                         "(or pass --allow-unbalanced).")

    print("\nWriting inputs:")
    for key in keys:
        write_inputs(key, species[key], force=args.force)


if __name__ == "__main__":
    main()