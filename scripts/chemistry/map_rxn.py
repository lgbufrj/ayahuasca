from rxnmapper import RXNMapper
from rdkit import Chem
from data import proteins


# ============================================================
# GLOBAL RULES (MECHANISM → NAC GEOMETRY)
# ============================================================
# NAC geometry is defined by THREE atoms: Donor (D), the
# transferred atom X, and Acceptor (A). The critical constraint
# is the D–X–A angle (near-linear for most enzymatic mechanisms)
# and the D–X / A–X distances separately, NOT the D–A distance.
#
# For hydride transfer:    C(donor)–H–C(acceptor)
# For proton transfer:     O/N(donor)–H–O/N(acceptor)
# For methyl transfer:     S/N(donor)–CH3–N/O(acceptor)
# For phosphoryl transfer: O(donor)–P–O(acceptor)
# For acyl transfer:       S/O(donor)–C=O–O/N(acceptor)

MECHANISM_RULES = {
    "hydride_transfer": {
        # H is the transferred atom. D–H and H–A cutoffs both tight.
        "donor_acceptor_cutoff": 3.5,   # C···C when H is midway
        "angle_min": 130,
        "angle_ideal": 165,
        "transferred_element": "H",
    },
    "proton_transfer": {
        # Shorter D–A than hydride because O/N are smaller than C.
        "donor_acceptor_cutoff": 3.2,
        "angle_min": 120,
        "angle_ideal": 175,
        "transferred_element": "H",
    },
    "methyl_transfer": {
        # SN2-like. D–C and C–A both ~2 Å at TS; D···A ~3.5–4 Å.
        "donor_acceptor_cutoff": 4.0,
        "angle_min": 150,
        "angle_ideal": 180,
        "transferred_element": "C",
    },
    "phosphoryl_transfer": {
        # Associative/dissociative. O–P–O angle near-linear.
        "donor_acceptor_cutoff": 4.5,
        "angle_min": 145,
        "angle_ideal": 175,
        "transferred_element": "P",
    },
    "acyl_transfer": {
        # Tetrahedral TS. Nucleophile attacks carbonyl C at ~105°.
        "donor_acceptor_cutoff": 3.5,
        "angle_min": 100,
        "angle_ideal": 109,
        "transferred_element": "C",
    },
    "oxidation": {
        # Generic: two-electron transfer with no single transferred atom.
        "donor_acceptor_cutoff": 4.5,
        "angle_min": 120,
        "angle_ideal": 180,
        "transferred_element": None,
    },
}


# ============================================================
# UTILITIES
# ============================================================

def atom_map_dict(mol):
    return {a.GetAtomMapNum(): a for a in mol.GetAtoms() if a.GetAtomMapNum()}


def bond_dict(mol):
    bonds = {}
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtom().GetAtomMapNum()
        a2 = bond.GetEndAtom().GetAtomMapNum()
        if a1 and a2:
            bonds[tuple(sorted((a1, a2)))] = bond.GetBondTypeAsDouble()
    return bonds


def hydrogen_count(atom):
    return atom.GetTotalNumHs()


# ============================================================
# PARSING
# ============================================================

def parse_mapped_reaction(mapped_rxn):
    r, p = mapped_rxn.split(">>")
    reactants = [Chem.MolFromSmiles(x) for x in r.split(".")]
    products  = [Chem.MolFromSmiles(x) for x in p.split(".")]
    return reactants, products


# ============================================================
# COLLECTION
# ============================================================

def collect_atoms(mols, heavy_only=False):
    atoms = {}
    for mol in mols:
        for atom in mol.GetAtoms():
            if heavy_only and atom.GetSymbol() == "H":
                continue
            m = atom.GetAtomMapNum()
            if m:
                atoms[m] = atom
    return atoms


def collect_bonds(mols):
    bonds = {}
    for mol in mols:
        bonds.update(bond_dict(mol))
    return bonds


# ============================================================
# BOND CHANGES
# ============================================================

def detect_bond_changes(r_bonds, p_bonds):
    formed, broken, changed = [], [], []
    for pair in set(r_bonds) | set(p_bonds):
        r = r_bonds.get(pair)
        p = p_bonds.get(pair)
        if   r is None: formed.append( {"atoms": pair, "bond_order": p})
        elif p is None: broken.append( {"atoms": pair, "bond_order": r})
        elif r != p:    changed.append({"atoms": pair, "reactant_order": r, "product_order": p})
    return {"formed": formed, "broken": broken, "changed": changed}


# ============================================================
# HYDROGEN CHANGES
# ============================================================

def detect_hydrogen_changes(r_atoms, p_atoms):
    changes = []
    for m in set(r_atoms) & set(p_atoms):
        delta = hydrogen_count(p_atoms[m]) - hydrogen_count(r_atoms[m])
        if delta:
            changes.append({
                "atom_map": m,
                "element":  r_atoms[m].GetSymbol(),
                "delta_h":  delta,
            })
    return changes


# ============================================================
# REACTION CENTER
# ============================================================

def reaction_center(bond_changes, h_changes):
    center = set()
    for cat in bond_changes.values():
        for x in cat:
            center.update(x["atoms"])
    for x in h_changes:
        center.add(x["atom_map"])
    return sorted(center)


# ============================================================
# MECHANISM CLASSIFIER
# ============================================================
# Priority order: phosphoryl > acyl > methyl > hydride >
# proton > oxidation (fallback). Each check looks at the
# *transferred atom's element* in the broken/formed bond pairs,
# not just whether delta_h exists.

def _involved_elements(bond_changes, r_atoms):
    """Return set of elements at atoms involved in bond changes."""
    elements = set()
    for cat in bond_changes.values():
        for x in cat:
            for atom_map in x["atoms"]:
                atom = r_atoms.get(atom_map)
                if atom:
                    elements.add(atom.GetSymbol())
    return elements


def _transferred_atom_element(bond_changes, r_atoms):
    """
    Find the element of the atom that appears in BOTH a broken bond
    and a formed bond — the true transferred atom at the reaction center.
    An atom that loses one bond and gains a new one is being transferred.
    """
    broken_atoms = set()
    formed_atoms = set()
    for x in bond_changes["broken"]:
        broken_atoms.update(x["atoms"])
    for x in bond_changes["formed"]:
        formed_atoms.update(x["atoms"])

    transferred = broken_atoms & formed_atoms
    elements = set()
    for m in transferred:
        a = r_atoms.get(m)
        if a:
            elements.add(a.GetSymbol())
    return elements


def classify_mechanism(bond_changes, h_changes, r_atoms):
    """
    Classify reaction mechanism from bond changes and H-count deltas.

    Decision logic (in priority order):
      1. Phosphoryl transfer: a phosphorus atom is the transferred atom.
      2. Acyl transfer: carbonyl carbon (C bonded to O with bond-order
         change) is transferred, leaving group is O or S.
      3. Methyl transfer: carbon atom is the transferred atom (SN2-like).
      4. Hydride transfer: H is transferred between two carbons.
         delta_h < 0 on a carbon AND delta_h > 0 on another carbon.
      5. Proton transfer: H is transferred to/from heteroatom (O, N, S).
      6. Oxidation (generic fallback): net electron change without a
         discrete transferred atom.
    """
    transferred = _transferred_atom_element(bond_changes, r_atoms)
    involved    = _involved_elements(bond_changes, r_atoms)

    # 1. Phosphoryl transfer
    if "P" in transferred or "P" in involved:
        return "phosphoryl_transfer"

    # 2. Acyl transfer — carbonyl C transferred, leaving group O or S
    if "C" in transferred and ("O" in involved or "S" in involved):
        for x in bond_changes.get("changed", []):
            atoms = [r_atoms.get(m) for m in x["atoms"]]
            syms  = {a.GetSymbol() for a in atoms if a}
            if "C" in syms and "O" in syms:
                return "acyl_transfer"

    # 3. Methyl transfer — C is transferred but not an acyl C
    if "C" in transferred:
        return "methyl_transfer"

    # 4. Hydride transfer — net H movement between two carbons
    donors_C    = [x for x in h_changes if x["element"] == "C" and x["delta_h"] < 0]
    acceptors_C = [x for x in h_changes if x["element"] == "C" and x["delta_h"] > 0]
    if donors_C and acceptors_C:
        return "hydride_transfer"

    # 5. Proton transfer — H movement involving at least one heteroatom
    heteroatoms   = {"O", "N", "S"}
    donors_het    = [x for x in h_changes if x["element"] in heteroatoms and x["delta_h"] < 0]
    acceptors_het = [x for x in h_changes if x["element"] in heteroatoms and x["delta_h"] > 0]
    if donors_het or acceptors_het:
        return "proton_transfer"

    return "oxidation"


# ============================================================
# DONOR / ACCEPTOR INFERENCE
# ============================================================
# Key fix from original: donor and acceptor must come from
# OPPOSITE sides of delta_h. Donor has delta_h < 0 (loses H or
# the transferred group). Acceptor has delta_h > 0 (gains it).
# We split the h_changes list first, THEN score within each half.

def _score_as_donor(atom, mechanism):
    """Higher = more likely to be the true donor for this mechanism."""
    sym = atom.GetSymbol()
    if sym == "H":
        return -999
    scores = {
        "hydride_transfer":    {"C": 100, "N": 20,  "O": 10},
        "proton_transfer":     {"O": 100, "N": 90,  "S": 60,  "C": 20},
        "methyl_transfer":     {"S": 100, "N": 80,  "O": 50},
        "phosphoryl_transfer": {"O": 100, "N": 60,  "S": 40},
        "acyl_transfer":       {"S": 100, "O": 80,  "N": 60},
        "oxidation":           {"C": 50,  "N": 50,  "O": 50,  "S": 50},
    }
    return scores.get(mechanism, {}).get(sym, 10)


def _score_as_acceptor(atom, mechanism):
    """Higher = more likely to be the true acceptor for this mechanism."""
    sym = atom.GetSymbol()
    if sym == "H":
        return -999
    scores = {
        "hydride_transfer":    {"C": 100, "N": 20,  "O": 10},
        "proton_transfer":     {"O": 100, "N": 90,  "S": 60,  "C": 20},
        "methyl_transfer":     {"N": 100, "O": 80,  "S": 60},
        "phosphoryl_transfer": {"O": 100, "N": 60,  "S": 40},
        "acyl_transfer":       {"O": 100, "N": 90,  "S": 60},
        "oxidation":           {"C": 50,  "N": 50,  "O": 50,  "S": 50},
    }
    return scores.get(mechanism, {}).get(sym, 10)


def infer_donor_acceptor(h_changes, r_atoms, bond_changes, mechanism):
    """
    Split h_changes into losers (delta_h < 0) and gainers (delta_h > 0),
    score each side separately, and return the best from each side.

    Falls back to bond_changes if h_changes is empty (e.g. for
    phosphoryl transfer where the transferred P has no H).
    """
    losers  = [x for x in h_changes if x["delta_h"] < 0]
    gainers = [x for x in h_changes if x["delta_h"] > 0]

    # For mechanisms where H isn't directly tracked (phosphoryl,
    # oxidation), use the bond-change atoms as a fallback.
    if not losers or not gainers:
        broken_atoms = [x["atoms"][0] for x in bond_changes["broken"]] if bond_changes["broken"] else []
        formed_atoms = [x["atoms"][1] for x in bond_changes["formed"]] if bond_changes["formed"] else []

        loser_map  = broken_atoms[0] if broken_atoms else None
        gainer_map = formed_atoms[0] if formed_atoms else None

        if not loser_map or not gainer_map:
            return None

        return {
            "donor":    {"atom_map": loser_map,  "element": r_atoms[loser_map].GetSymbol()  if loser_map  in r_atoms else "?"},
            "acceptor": {"atom_map": gainer_map, "element": r_atoms[gainer_map].GetSymbol() if gainer_map in r_atoms else "?"},
        }

    def best(candidates, score_fn):
        scored = [(score_fn(r_atoms[x["atom_map"]], mechanism), x)
                  for x in candidates if x["atom_map"] in r_atoms]
        if not scored:
            return candidates[0]
        return max(scored, key=lambda t: t[0])[1]

    return {
        "donor":    best(losers,  _score_as_donor),
        "acceptor": best(gainers, _score_as_acceptor),
    }


# ============================================================
# TRANSFERRED ATOM IDENTIFICATION
# ============================================================

def find_transferred_atom(bond_changes, r_atoms, mechanism):
    """
    Return the atom_map of the transferred atom X (e.g. H, CH3 carbon,
    phosphorus). Returns None for mechanisms like oxidation where no
    discrete atom is transferred.

    If multiple candidates exist, prefer the element expected for
    this mechanism.
    """
    expected = MECHANISM_RULES.get(mechanism, {}).get("transferred_element")
    if expected is None:
        return None

    broken_maps = set()
    for x in bond_changes["broken"]:
        broken_maps.update(x["atoms"])
    formed_maps = set()
    for x in bond_changes["formed"]:
        formed_maps.update(x["atoms"])

    candidates = broken_maps & formed_maps
    for m in candidates:
        a = r_atoms.get(m)
        if a and a.GetSymbol() == expected:
            return m

    # Fallback: return any candidate if element match fails
    # (can happen when RXNMapper assigns H atoms inconsistently)
    return next(iter(candidates), None)


# ============================================================
# NAC GENERATION (THREE-ATOM: DONOR – X – ACCEPTOR)
# ============================================================

def generate_nac(mechanism, donor, acceptor, transferred_atom_map):
    """
    Build NAC constraints using the three-atom representation.

    distance_pairs tracks:
      (a) donor    ↔ transferred atom  (bond being broken)
      (b) acceptor ↔ transferred atom  (bond being formed)

    The angle is:
      donor – transferred_atom – acceptor  (the attack angle)

    If no transferred atom is known (oxidation or RXNMapper could not
    assign explicit H), fall back to the two-atom donor↔acceptor
    representation with NO angle (a degenerate D–X–A angle where X
    equals one of the endpoints is meaningless geometrically).
    """
    rules = MECHANISM_RULES.get(mechanism, MECHANISM_RULES["oxidation"])

    if transferred_atom_map is None:
        # Two-atom fallback — no angle emitted
        return {
            "mechanism": mechanism,
            "distance_pairs": [{
                "atom1":  donor["atom_map"],
                "atom2":  acceptor["atom_map"],
                "cutoff": rules["donor_acceptor_cutoff"],
                "label":  "donor–acceptor",
            }],
            "angles": [],   # no vertex → no meaningful angle
        }

    return {
        "mechanism": mechanism,
        "distance_pairs": [
            {
                "atom1":  donor["atom_map"],
                "atom2":  transferred_atom_map,
                "cutoff": rules["donor_acceptor_cutoff"] * 0.55,
                "label":  "donor–transferred_atom (breaking bond)",
            },
            {
                "atom1":  acceptor["atom_map"],
                "atom2":  transferred_atom_map,
                "cutoff": rules["donor_acceptor_cutoff"] * 0.55,
                "label":  "acceptor–transferred_atom (forming bond)",
            },
        ],
        "angles": [{
            "atom1":   donor["atom_map"],
            "vertex":  transferred_atom_map,
            "atom2":   acceptor["atom_map"],
            "minimum": rules["angle_min"],
            "ideal":   rules["angle_ideal"],
            "label":   "donor–X–acceptor attack angle",
        }],
    }


# ============================================================
# FULL ANALYSIS PIPELINE
# ============================================================

def analyze(mapped_rxn):

    reactants, products = parse_mapped_reaction(mapped_rxn)

    # Chemistry layer: keep H for delta_h counting
    r_atoms_full = collect_atoms(reactants, heavy_only=False)
    p_atoms_full = collect_atoms(products,  heavy_only=False)

    r_bonds = collect_bonds(reactants)
    p_bonds = collect_bonds(products)

    bond_changes = detect_bond_changes(r_bonds, p_bonds)
    h_changes    = detect_hydrogen_changes(r_atoms_full, p_atoms_full)

    center    = reaction_center(bond_changes, h_changes)
    mechanism = classify_mechanism(bond_changes, h_changes, r_atoms_full)

    da = infer_donor_acceptor(h_changes, r_atoms_full, bond_changes, mechanism)
    if da is None:
        return {
            "mechanism":        mechanism,
            "reaction_center":  center,
            "donor_acceptor":   None,
            "nac":              None,
            "warning":          "Could not determine donor/acceptor atoms.",
        }

    transferred_map = find_transferred_atom(bond_changes, r_atoms_full, mechanism)
    nac             = generate_nac(mechanism, da["donor"], da["acceptor"], transferred_map)

    return {
        "mechanism":            mechanism,
        "reaction_center":      center,
        "donor_acceptor":       da,
        "transferred_atom_map": transferred_map,
        "nac":                  nac,
    }


# ============================================================
# MAIN
# ============================================================

if __name__ == "__main__":

    protein   = "asmt"
    rxn_mapper = RXNMapper()
    rxns       = [r["smiles"] for r in proteins[protein]["reactions"]]
    mapped_rxns = rxn_mapper.get_attention_guided_atom_maps(rxns)

    for r in mapped_rxns:
        results = analyze(r["mapped_rxn"])
        print(results)