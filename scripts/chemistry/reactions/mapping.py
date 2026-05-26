# chemistry/reactions/mapping.py

from rdkit import Chem

# ============================================================
# PARSE REACTION
# ============================================================

def parse_mapped_reaction(mapped_rxn):
    """
    Split mapped reaction SMILES into:
        reactants
        products
    """

    reactants_smiles, products_smiles = mapped_rxn.split(">>")

    reactants = [
        Chem.MolFromSmiles(x)
        for x in reactants_smiles.split(".")
    ]

    products = [
        Chem.MolFromSmiles(x)
        for x in products_smiles.split(".")
    ]

    return reactants, products


# ============================================================
# ATOM MAP DICT
# ============================================================

def atom_map_dict(mol):

    mapping = {}

    for atom in mol.GetAtoms():

        map_num = atom.GetAtomMapNum()

        if map_num != 0:

            mapping[map_num] = atom

    return mapping


# ============================================================
# BOND DICT
# ============================================================

def bond_dict(mol):

    bonds = {}

    for bond in mol.GetBonds():

        a1 = bond.GetBeginAtom().GetAtomMapNum()
        a2 = bond.GetEndAtom().GetAtomMapNum()

        if a1 == 0 or a2 == 0:
            continue

        pair = tuple(sorted((a1, a2)))

        bonds[pair] = bond.GetBondTypeAsDouble()

    return bonds


# ============================================================
# COLLECT ATOMS
# ============================================================

def collect_atoms(molecules):

    atoms = {}

    for mol in molecules:

        atoms.update(
            atom_map_dict(mol)
        )

    return atoms


# ============================================================
# COLLECT BONDS
# ============================================================

def collect_bonds(molecules):

    bonds = {}

    for mol in molecules:

        bonds.update(
            bond_dict(mol)
        )

    return bonds


# ============================================================
# HYDROGEN COUNT
# ============================================================

def hydrogen_count(atom):

    return (
        atom.GetNumExplicitHs()
        + atom.GetNumImplicitHs()
    )


# ============================================================
# DETECT BOND CHANGES
# ============================================================

def detect_bond_changes(
    reactant_bonds,
    product_bonds
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

        # ----------------------------------------------------
        # formed
        # ----------------------------------------------------

        if r is None and p is not None:

            formed.append({
                "atoms": pair,
                "bond_order": p
            })

        # ----------------------------------------------------
        # broken
        # ----------------------------------------------------

        elif r is not None and p is None:

            broken.append({
                "atoms": pair,
                "bond_order": r
            })

        # ----------------------------------------------------
        # changed
        # ----------------------------------------------------

        elif r != p:

            changed.append({
                "atoms": pair,
                "reactant_order": r,
                "product_order": p
            })

    return {
        "formed": formed,
        "broken": broken,
        "changed": changed
    }


# ============================================================
# DETECT HYDROGEN CHANGES
# ============================================================

def detect_hydrogen_changes(
    reactant_atoms,
    product_atoms
):

    changes = []

    common_maps = (
        set(reactant_atoms.keys())
        &
        set(product_atoms.keys())
    )

    for map_num in common_maps:

        r_atom = reactant_atoms[map_num]
        p_atom = product_atoms[map_num]

        r_h = hydrogen_count(r_atom)
        p_h = hydrogen_count(p_atom)

        if r_h != p_h:

            changes.append({

                "atom_map": map_num,

                "element": r_atom.GetSymbol(),

                "reactant_h": r_h,

                "product_h": p_h,

                "delta_h": p_h - r_h
            })

    return changes


# ============================================================
# FIND ATOM PROVENANCE
# ============================================================

def find_atom_provenance(
    atom_map,
    molecules,
    reaction_metadata,
    cofactors=[]
):
    """
    Find which compound/role an atom belongs to.
    """

    substrates = reaction_metadata.get(
        "substrates",
        []
    )

    for mol_idx, mol in enumerate(molecules):

        for atom in mol.GetAtoms():

            if atom.GetAtomMapNum() != atom_map:
                continue
            
            compound = substrates[mol_idx]

            role = (
                "cofactor"
                if compound in cofactors
                else "substrate"
            )

            return {

                "atom_map": atom_map,

                "compound": compound,

                "role": role,

                "molecule_index": mol_idx,

                "element": atom.GetSymbol()
            }

    return None


# ============================================================
# ASSIGN PROVENANCE
# ============================================================

def assign_provenance(
    bond_changes,
    hydrogen_changes,
    reactants,
    reaction_metadata,
    cofactors=[]
):

    # --------------------------------------------------------
    # Bond changes
    # --------------------------------------------------------

    for category in [
        "formed",
        "broken",
        "changed"
    ]:

        for item in bond_changes[category]:

            a1, a2 = item["atoms"]

            item["provenance"] = {

                "atom1": find_atom_provenance(
                    a1,
                    reactants,
                    reaction_metadata,
                    cofactors=cofactors
                ),

                "atom2": find_atom_provenance(
                    a2,
                    reactants,
                    reaction_metadata,
                    cofactors=cofactors
                )
            }

    # --------------------------------------------------------
    # Hydrogen changes
    # --------------------------------------------------------

    for item in hydrogen_changes:

        item["provenance"] = find_atom_provenance(
            item["atom_map"],
            reactants,
            reaction_metadata,
            cofactors=cofactors
        )