from rdkit import Chem


# ============================================================
# ATOM MAPS
# ============================================================

def atom_map_dict(mol):

    mapping = {}

    for atom in mol.GetAtoms():

        map_num = atom.GetAtomMapNum()

        if map_num != 0:
            mapping[map_num] = atom

    return mapping


# ============================================================
# BONDS
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
# HYDROGEN COUNT
# ============================================================

def hydrogen_count(atom):

    return atom.GetTotalNumHs()


# ============================================================
# COLLECTION
# ============================================================

def collect_atoms(molecules):

    atoms = {}

    for mol in molecules:

        atoms.update(
            atom_map_dict(mol)
        )

    return atoms



def collect_bonds(molecules):

    bonds = {}

    for mol in molecules:

        bonds.update(
            bond_dict(mol)
        )

    return bonds