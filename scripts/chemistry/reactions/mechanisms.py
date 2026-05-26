# chemistry/reactions/mechanisms.py

from chemistry.reactions.rules import (
    MECHANISM_RULES
)

# ============================================================
# INFER MECHANISM
# ============================================================

def infer_mechanism(
    bond_changes,
    hydrogen_changes,
    reaction_metadata=None
):

    formed = bond_changes["formed"]
    broken = bond_changes["broken"]
    changed = bond_changes["changed"]

    # ========================================================
    # GROUP TRANSFERS
    # ========================================================

    for mechanism, rule in MECHANISM_RULES.items():

        if rule["type"] != "group_transfer":
            continue

        required = rule.get(
            "requires",
            {}
        )

        if (
            len(formed)
            !=
            required.get(
                "formed_bonds",
                len(formed)
            )
        ):
            continue

        if (
            len(broken)
            !=
            required.get(
                "broken_bonds",
                len(broken)
            )
        ):
            continue

        transferred_element = rule.get(
            "transferred_element"
        )

        donor_elements = set(
            rule.get(
                "donor_elements",
                []
            )
        )

        acceptor_elements = set(
            rule.get(
                "acceptor_elements",
                []
            )
        )

        found = False

        for f in formed:

            fp = f["provenance"]

            if fp is None:
                continue

            a1 = fp["atom1"]
            a2 = fp["atom2"]

            if a1 is None or a2 is None:
                continue

            formed_atoms = [
                a1,
                a2
            ]

            transferred_atom = None
            acceptor_atom = None

            for atom in formed_atoms:

                if (
                    atom["element"]
                    ==
                    transferred_element
                ):
                    transferred_atom = atom

                elif (
                    atom["element"]
                    in
                    acceptor_elements
                ):
                    acceptor_atom = atom

            if (
                transferred_atom is None
                or acceptor_atom is None
            ):
                continue

            for b in broken:

                bp = b["provenance"]

                if bp is None:
                    continue

                b1 = bp["atom1"]
                b2 = bp["atom2"]

                if b1 is None or b2 is None:
                    continue

                broken_atoms = [
                    b1,
                    b2
                ]

                donor_atom = None

                for atom in broken_atoms:

                    if (
                        atom["element"]
                        in
                        donor_elements
                    ):
                        donor_atom = atom

                if donor_atom is None:
                    continue

                transferred_maps = {
                    b1["atom_map"],
                    b2["atom_map"]
                }

                if (
                    transferred_atom["atom_map"]
                    in
                    transferred_maps
                ):
                    found = True
                    break

            if found:
                return mechanism

    # ========================================================
    # REDOX
    # ========================================================

    has_hydrogen_transfer = (
        len(hydrogen_changes) > 0
    )

    has_bond_order_changes = (
        len(changed) > 0
    )

    has_cofactor_hydrogen = any(

        (
            h.get("provenance")
            and
            h["provenance"].get("role")
            == "cofactor"
        )

        for h in hydrogen_changes
    )

    if (
        has_hydrogen_transfer
        and has_bond_order_changes
        and has_cofactor_hydrogen
    ):
        return "hydride_transfer"

    return "unknown"


# ============================================================
# PRIORITY SCORING
# ============================================================

def score_transfer_candidate(
    mechanism,
    atom_data,
    role_preference
):
    """
    Scores transfer candidates dynamically.
    """

    score = 0

    role = atom_data["provenance"]["role"]
    element = atom_data["element"]

    # --------------------------------------------------------
    # Role preference
    # --------------------------------------------------------

    if role == role_preference:
        score += 10

    # --------------------------------------------------------
    # Mechanism-specific priorities
    # --------------------------------------------------------

    rules = MECHANISM_RULES.get(
        mechanism,
        {}
    )

    priorities = rules.get(
        "priority",
        {}
    )

    score += priorities.get(
        element,
        0
    )

    return score


# ============================================================
# INFER TRANSFER VECTOR
# ============================================================

def infer_transfer_vector(
    mechanism,
    bond_changes,
    hydrogen_changes,
    reaction_metadata=None
):

    if mechanism == "unknown":
        return None

    rule = MECHANISM_RULES.get(
        mechanism
    )

    if rule is None:
        return None

    # ========================================================
    # GROUP TRANSFER
    # ========================================================

    if rule["type"] == "group_transfer":

        transferred_element = rule.get(
            "transferred_element"
        )

        donor_elements = set(
            rule.get(
                "donor_elements",
                []
            )
        )

        acceptor_elements = set(
            rule.get(
                "acceptor_elements",
                []
            )
        )

        donor = None
        acceptor = None

        # ----------------------------------------------------
        # Find acceptor
        # ----------------------------------------------------

        for formed in bond_changes["formed"]:

            prov = formed["provenance"]

            if prov is None:
                continue

            atom1 = prov["atom1"]
            atom2 = prov["atom2"]

            if atom1 is None or atom2 is None:
                continue

            atoms = [atom1, atom2]

            transferred_atom = None
            acceptor_atom = None

            for atom in atoms:

                if (
                    atom["element"]
                    ==
                    transferred_element
                ):
                    transferred_atom = atom

                elif (
                    atom["element"]
                    in
                    acceptor_elements
                ):
                    acceptor_atom = atom

            if (
                transferred_atom is not None
                and
                acceptor_atom is not None
            ):

                donor = transferred_atom
                acceptor = acceptor_atom

                break

        # ----------------------------------------------------
        # Find donor attachment atom
        # ----------------------------------------------------

        donor_attachment = None

        if donor is not None:

            donor_map = donor["atom_map"]

            for broken in bond_changes["broken"]:

                prov = broken["provenance"]

                if prov is None:
                    continue

                atom1 = prov["atom1"]
                atom2 = prov["atom2"]

                if atom1 is None or atom2 is None:
                    continue

                atoms = [atom1, atom2]

                maps = {
                    atom1["atom_map"],
                    atom2["atom_map"]
                }

                if donor_map not in maps:
                    continue

                for atom in atoms:

                    if (
                        atom["atom_map"]
                        != donor_map
                        and
                        atom["element"]
                        in donor_elements
                    ):
                        donor_attachment = atom

        return {

            "donor": donor_attachment,

            "acceptor": acceptor,

            "transferred_group": donor
        }

    # ------------------------------------------------------------
    # REDOX / HYDRIDE TRANSFER
    # ------------------------------------------------------------

    elif rule["type"] == "redox":

        donor = None
        acceptor = None

        preferred_acceptors = set(
            x.lower()
            for x in rule.get(
                "preferred_acceptor_compounds",
                []
            )
        )

        # --------------------------------------------------------
        # classify H losses/gains
        # --------------------------------------------------------

        losses = [
            h for h in hydrogen_changes
            if h["delta_h"] < 0
        ]

        gains = [
            h for h in hydrogen_changes
            if h["delta_h"] > 0
        ]

        # --------------------------------------------------------
        # identify cofactor-side gain/loss
        # --------------------------------------------------------

        cofactor_gain = None
        substrate_gain = None

        cofactor_loss = None
        substrate_loss = None
        
        preferred_donor_elements = rule.get(
            "preferred_donor_elements",
            []
        )
        preferred_acceptor_elements = rule.get(
            "preferred_acceptor_elements",
            []
        )

        for h in gains:

            prov = h.get("provenance")

            if prov is None:
                continue

            compound = prov["compound"].lower()

            if compound in preferred_acceptors:
                if (
                    h["element"]
                    in preferred_acceptor_elements
                ):
                    cofactor_gain = h
            else:
                substrate_gain = h

        for h in losses:

            prov = h.get("provenance")

            if prov is None:
                continue

            compound = prov["compound"].lower()

            if compound in preferred_acceptors:
                if (
                    h["element"]
                    in preferred_donor_elements
                ):
                    cofactor_loss = h
            else:
                substrate_loss = h

        # --------------------------------------------------------
        # REDUCTION
        # substrate gains hydride from NADPH
        # --------------------------------------------------------

        if substrate_gain and cofactor_loss:

            donor = cofactor_loss
            acceptor = substrate_gain

        # --------------------------------------------------------
        # OXIDATION
        # substrate loses hydride to NADP+
        # --------------------------------------------------------

        elif substrate_loss and cofactor_gain:

            donor = substrate_loss
            acceptor = cofactor_gain

        # --------------------------------------------------------
        # fallback generic
        # --------------------------------------------------------

        else:

            donor = next(
                (
                    h for h in hydrogen_changes
                    if h["delta_h"] < 0
                ),
                None
            )

            acceptor = next(
                (
                    h for h in hydrogen_changes
                    if h["delta_h"] > 0
                ),
                None
            )

        return {
            "donor": donor,
            "acceptor": acceptor
        }

    return None