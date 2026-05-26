# ============================================================
# chemistry/reactions/rules.py
# ============================================================

MECHANISM_RULES = {

    # ========================================================
    # METHYL TRANSFER
    # ========================================================

    "methyl_transfer": {

        "type": "group_transfer",

        "transferred_element": "C",

        "donor_elements": ["S"],

        "acceptor_elements": ["O", "N"],

        "requires": {

            "formed_bonds": 1,
            "broken_bonds": 1
        }
    },

    # ========================================================
    # HYDRIDE TRANSFER
    # ========================================================

    "hydride_transfer": {

        "type": "redox",

        # ----------------------------------------
        # REQUIRED FEATURES
        # ----------------------------------------

        "requires": {

            "hydrogen_transfer": True,

            "hydrogen_loss": 1,

            "hydrogen_gain": 1,

            "bond_order_changes": True
        },

        # ----------------------------------------
        # TRANSFER PROPERTIES
        # ----------------------------------------

        "transferred_element": "H",

        "donor_elements": ["C"],

        "acceptor_elements": ["C", "N"],

        # ----------------------------------------
        # HEURISTICS / SCORING
        # ----------------------------------------

        "preferred_acceptor_compounds": [
            "nad",
            "nad+",
            "nadp",
            "nadp+",
        ],
        "preferred_donor_elements": ["C"],
        "preferred_acceptor_elements": ["C"]
    },

    # ========================================================
    # ACETYL TRANSFER
    # ========================================================

    "acetyl_transfer": {

        "type": "group_transfer",

        "transferred_fragment_size": 2,

        "donor_elements": ["S"],

        "acceptor_elements": ["O", "N"]
    }
}