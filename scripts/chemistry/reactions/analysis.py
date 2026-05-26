# chemistry/reactions/analysis.py

from chemistry.reactions.mapping import (
    parse_mapped_reaction,
    collect_atoms,
    collect_bonds,
    detect_bond_changes,
    detect_hydrogen_changes,
    assign_provenance,
)

from chemistry.reactions.reaction_center import (
    infer_reaction_center
)

from chemistry.reactions.mechanisms import (
    infer_mechanism,
    infer_transfer_vector,
)

# ============================================================
# MAIN ANALYSIS
# ============================================================

def analyze_reaction(
    mapped_rxn,
    reaction_metadata=None,
    cofactors=[]
):
    """
    Full reaction analysis pipeline.

    Parameters
    ----------
    mapped_rxn : str
        RXNMapper mapped reaction SMILES

    reaction_metadata : dict
        Original reaction config from proteins.yaml

    Returns
    -------
    dict
    """

    # --------------------------------------------------------
    # Parse reaction
    # --------------------------------------------------------

    reactants, products = parse_mapped_reaction(
        mapped_rxn
    )

    # --------------------------------------------------------
    # Collect atoms / bonds
    # --------------------------------------------------------

    reactant_atoms = collect_atoms(
        reactants
    )

    product_atoms = collect_atoms(
        products
    )

    reactant_bonds = collect_bonds(
        reactants
    )

    product_bonds = collect_bonds(
        products
    )

    # --------------------------------------------------------
    # Detect chemistry changes
    # --------------------------------------------------------

    bond_changes = detect_bond_changes(
        reactant_bonds,
        product_bonds
    )

    hydrogen_changes = detect_hydrogen_changes(
        reactant_atoms,
        product_atoms
    )

    # --------------------------------------------------------
    # Assign provenance
    # --------------------------------------------------------

    assign_provenance(
        bond_changes=bond_changes,
        hydrogen_changes=hydrogen_changes,
        reactants=reactants,
        reaction_metadata=reaction_metadata,
        cofactors=cofactors
    )

    # --------------------------------------------------------
    # Reaction center
    # --------------------------------------------------------

    reaction_center = infer_reaction_center(
        bond_changes,
        hydrogen_changes
    )

    # --------------------------------------------------------
    # Mechanism inference
    # --------------------------------------------------------

    mechanism = infer_mechanism(
        bond_changes=bond_changes,
        hydrogen_changes=hydrogen_changes,
        reaction_metadata=reaction_metadata
    )

    # --------------------------------------------------------
    # Transfer vector
    # --------------------------------------------------------

    transfer_vector = infer_transfer_vector(
        mechanism=mechanism,
        bond_changes=bond_changes,
        hydrogen_changes=hydrogen_changes,
        reaction_metadata=reaction_metadata
    )

    # --------------------------------------------------------
    # Final
    # --------------------------------------------------------

    return {

        "mechanism": mechanism,

        "reaction_center": reaction_center,

        "transfer_vector": transfer_vector,

        "bond_changes": bond_changes,

        "hydrogen_changes": hydrogen_changes,
    }