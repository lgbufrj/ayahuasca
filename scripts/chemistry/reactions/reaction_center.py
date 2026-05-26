# chemistry/reactions/reaction_center.py

from rdkit import Chem

from chemistry.helpers.rdkit_util import (
    collect_atoms,
    collect_bonds,
    hydrogen_count
)


# ============================================================
# PARSE MAPPED REACTION
# ============================================================

def parse_mapped_reaction(
    mapped_rxn,
    reaction_metadata=None
):
    """
    Parses mapped reaction while preserving:
        - molecule identity
        - substrate/cofactor/product role
        - atom provenance
    """

    reactants_smiles, products_smiles = mapped_rxn.split(">>")

    reactant_smiles = reactants_smiles.split(".")
    product_smiles  = products_smiles.split(".")

    reactants = []
    products  = []

    # --------------------------------------------------------
    # Metadata-aware parsing
    # --------------------------------------------------------

    substrates = []
    products_metadata = []

    if reaction_metadata is not None:

        substrates = reaction_metadata.get(
            "substrates",
            []
        )

        products_metadata = reaction_metadata.get(
            "products",
            []
        )

    # --------------------------------------------------------
    # Reactants
    # --------------------------------------------------------

    for i, smiles in enumerate(reactant_smiles):

        mol = Chem.MolFromSmiles(smiles)

        compound = (
            substrates[i]
            if i < len(substrates)
            else f"reactant_{i}"
        )

        role = infer_role(compound)

        reactants.append({

            "mol": mol,

            "compound": compound,

            "role": role,

            "index": i
        })

    # --------------------------------------------------------
    # Products
    # --------------------------------------------------------

    for i, smiles in enumerate(product_smiles):

        mol = Chem.MolFromSmiles(smiles)

        compound = (
            products_metadata[i]
            if i < len(products_metadata)
            else f"product_{i}"
        )

        role = infer_role(compound)

        products.append({

            "mol": mol,

            "compound": compound,

            "role": role,

            "index": i
        })

    return reactants, products


# ============================================================
# ROLE INFERENCE
# ============================================================

COFACTORS = {

    "nadph",
    "nadp+",

    "nadh",
    "nad+",

    "sam",
    "sah",

    "atp",
    "adp",

    "acetyl-coa",
    "coa"
}


def infer_role(compound):

    if compound.lower() in COFACTORS:
        return "cofactor"

    return "substrate"


# ============================================================
# ATOM PROVENANCE
# ============================================================

def build_atom_provenance(molecules):

    provenance = {}

    for entry in molecules:

        mol = entry["mol"]

        compound = entry["compound"]

        role = entry["role"]

        molecule_index = entry["index"]

        for atom in mol.GetAtoms():

            atom_map = atom.GetAtomMapNum()

            if atom_map == 0:
                continue

            provenance[atom_map] = {

                "compound": compound,

                "role": role,

                "molecule_index": molecule_index,

                "element": atom.GetSymbol()
            }

    return provenance


# ============================================================
# BOND CHANGES
# ============================================================

def detect_bond_changes(
    reactant_bonds,
    product_bonds,
    atom_provenance=None
):

    formed = []
    broken = []
    changed = []

    all_pairs = set(
        reactant_bonds.keys()
    ) | set(
        product_bonds.keys()
    )

    for pair in all_pairs:

        r = reactant_bonds.get(pair)
        p = product_bonds.get(pair)

        atom1, atom2 = pair

        provenance = None

        if atom_provenance is not None:

            provenance = {

                "atom1": atom_provenance.get(atom1),

                "atom2": atom_provenance.get(atom2)
            }

        # ----------------------------------------------------
        # formed
        # ----------------------------------------------------

        if r is None and p is not None:

            formed.append({

                "atoms": pair,

                "bond_order": p,

                "provenance": provenance
            })

        # ----------------------------------------------------
        # broken
        # ----------------------------------------------------

        elif r is not None and p is None:

            broken.append({

                "atoms": pair,

                "bond_order": r,

                "provenance": provenance
            })

        # ----------------------------------------------------
        # changed
        # ----------------------------------------------------

        elif r != p:

            changed.append({

                "atoms": pair,

                "reactant_order": r,

                "product_order": p,

                "provenance": provenance
            })

    return {

        "formed": formed,

        "broken": broken,

        "changed": changed
    }


# ============================================================
# HYDROGEN CHANGES
# ============================================================

def detect_hydrogen_changes(
    reactant_atoms,
    product_atoms,
    atom_provenance=None
):

    changes = []

    common_maps = set(
        reactant_atoms.keys()
    ) & set(
        product_atoms.keys()
    )

    for map_num in common_maps:

        r_atom = reactant_atoms[map_num]
        p_atom = product_atoms[map_num]

        r_h = hydrogen_count(r_atom)
        p_h = hydrogen_count(p_atom)

        if r_h != p_h:

            entry = {

                "atom_map": map_num,

                "element": r_atom.GetSymbol(),

                "reactant_h": r_h,

                "product_h": p_h,

                "delta_h": p_h - r_h
            }

            if atom_provenance is not None:

                entry["provenance"] = (
                    atom_provenance.get(map_num)
                )

            changes.append(entry)

    return changes


# ============================================================
# REACTION CENTER
# ============================================================

def infer_reaction_center(
    bond_changes,
    hydrogen_changes
):

    center_atoms = set()

    for category in [
        "formed",
        "broken",
        "changed"
    ]:

        for item in bond_changes[category]:

            center_atoms.update(
                item["atoms"]
            )

    for item in hydrogen_changes:

        center_atoms.add(
            item["atom_map"]
        )

    return sorted(center_atoms)


# ============================================================
# MAIN ANALYSIS
# ============================================================

def analyze_reaction_center(
    mapped_rxn,
    reaction_metadata=None
):

    reactants, products = parse_mapped_reaction(
        mapped_rxn,
        reaction_metadata
    )

    # --------------------------------------------------------
    # Raw molecules
    # --------------------------------------------------------

    reactant_mols = [
        x["mol"]
        for x in reactants
    ]

    product_mols = [
        x["mol"]
        for x in products
    ]

    # --------------------------------------------------------
    # Atom provenance
    # --------------------------------------------------------

    reactant_provenance = build_atom_provenance(
        reactants
    )

    product_provenance = build_atom_provenance(
        products
    )

    atom_provenance = {

        **reactant_provenance,
        **product_provenance
    }

    # --------------------------------------------------------
    # Collect atoms
    # --------------------------------------------------------

    reactant_atoms = collect_atoms(
        reactant_mols
    )

    product_atoms = collect_atoms(
        product_mols
    )

    # --------------------------------------------------------
    # Collect bonds
    # --------------------------------------------------------

    reactant_bonds = collect_bonds(
        reactant_mols
    )

    product_bonds = collect_bonds(
        product_mols
    )

    # --------------------------------------------------------
    # Detect changes
    # --------------------------------------------------------

    bond_changes = detect_bond_changes(
        reactant_bonds,
        product_bonds,
        atom_provenance
    )

    hydrogen_changes = detect_hydrogen_changes(
        reactant_atoms,
        product_atoms,
        atom_provenance
    )

    # --------------------------------------------------------
    # Reaction center
    # --------------------------------------------------------

    reaction_center = infer_reaction_center(
        bond_changes,
        hydrogen_changes
    )

    return {

        "reactants": reactants,

        "products": products,

        "atom_provenance": atom_provenance,

        "bond_changes": bond_changes,

        "hydrogen_changes": hydrogen_changes,

        "reaction_center": reaction_center,

        "reactant_atoms": reactant_atoms,

        "product_atoms": product_atoms
    }