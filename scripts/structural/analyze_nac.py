import json
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs, rdFMCS
from Bio.PDB import PDBParser
from Bio.PDB.MMCIFParser import MMCIFParser

from data import PROTEINS_PATH, COMPOUNDS_PATH, proteins, compounds
from chemistry.map_rxn import analyze as analyze_reaction, RXNMapper


# ============================================================
# COMPOUND FINGERPRINTS
# ============================================================

def _build_compound_fps(compounds):
    """
    Load each compound's SDF and compute a Morgan fingerprint.
    Built once at startup and reused across all reactions.
    """
    fps = {}
    for name in compounds:
        sdf_path = (
            f"{COMPOUNDS_PATH}/{name}/structure"
            f"/{name}_{compounds[name]['pubchem_id']}.sdf"
        )
        try:
            supplier = Chem.SDMolSupplier(str(sdf_path), removeHs=False)
            ref_mol  = supplier[0] if len(supplier) > 0 else None
        except Exception:
            ref_mol = None

        if ref_mol is not None:
            fps[name] = AllChem.GetMorganFingerprintAsBitVect(ref_mol, 2, 2048)

    return fps


# ============================================================
# FRAGMENT ROLE ASSIGNMENT
# ============================================================

def _assign_fragment_roles(fragments, protein_data, compound_fps, cofactor_names):
    """
    For each dot-separated reactant fragment, decide whether it is a
    "cofactor" or a "substrate" by fingerprint similarity against every
    known compound.

    A fragment is labelled "cofactor" if:
      - its best Tanimoto match (Morgan r=2) scores >= 0.85, AND
      - the matched compound name is in protein_data["cofactors"]

    Otherwise it is labelled "substrate".

    The 0.85 threshold is intentionally loose enough to survive heavy-H
    annotation differences between the mapped SMILES and the reference SDF,
    but tight enough to separate NAD+ from NADPH.

    Fragments that cannot be parsed (e.g. "[1H+]") fall through to
    "substrate" — they carry no coordinate atoms and will never appear
    in a NAC rule.
    """
    roles = []
    for frag_smiles in fragments:
        frag_mol = Chem.MolFromSmiles(frag_smiles)
        if frag_mol is None:
            roles.append("substrate")
            continue

        frag_fp = AllChem.GetMorganFingerprintAsBitVect(frag_mol, 2, 2048)

        best_name, best_sim = None, 0.0
        for name, ref_fp in compound_fps.items():
            sim = DataStructs.TanimotoSimilarity(frag_fp, ref_fp)
            if sim > best_sim:
                best_sim, best_name = sim, name

        if best_sim >= 0.85 and best_name in cofactor_names:
            roles.append("cofactor")
        else:
            roles.append("substrate")

    return roles


# ============================================================
# MCS-BASED ATOM MAP -> SDF INDEX RESOLUTION
# ============================================================

def _build_map_to_sdf_index(frag_smiles, sdf_mol):
    """
    Align a (mapped) reaction SMILES fragment onto a PubChem SDF molecule
    using Maximum Common Substructure, then return a dict:

        { atom_map_num: sdf_atom_index, ... }

    This is the only robust way to connect reaction atom map numbers
    (from RXNMapper) to positional indices in a PubChem SDF, because:
      - PubChem SDFs carry no atom map numbers
      - PubChem SDFs carry no atom names that RDKit exposes via _Name
      - Positional order between SMILES and SDF is not guaranteed
      - SMARTS built from the mapped SMILES fail on the SDF due to
        heavy-H annotation differences in neighbor connectivity

    The MCS is run on heavy atoms only (H stripped from both sides)
    with element and bond-order matching. The result is a bijection
    between SMILES atom indices and SDF atom indices for the matched
    subgraph; we then look up atom map numbers from the SMILES side.
    """
    smiles_mol = Chem.MolFromSmiles(frag_smiles)
    if smiles_mol is None or sdf_mol is None:
        return {}

    # Strip map numbers from a clean copy for MCS (they interfere with matching)
    smiles_clean = Chem.RWMol(Chem.RemoveHs(smiles_mol))
    for atom in smiles_clean.GetAtoms():
        atom.SetAtomMapNum(0)

    sdf_clean = Chem.RemoveHs(sdf_mol)

    result = rdFMCS.FindMCS(
        [smiles_clean, sdf_clean],
        atomCompare=rdFMCS.AtomCompare.CompareElements,
        bondCompare=rdFMCS.BondCompare.CompareOrder,
        matchValences=False,
        ringMatchesRingOnly=False,
        completeRingsOnly=False,
        timeout=10,
    )

    if result.canceled or result.numAtoms == 0:
        return {}

    query          = Chem.MolFromSmarts(result.smartsString)
    smiles_matches = smiles_clean.GetSubstructMatches(query)
    sdf_matches    = sdf_clean.GetSubstructMatches(query)

    if not smiles_matches or not sdf_matches:
        return {}

    # Use the first (largest) match from each side.
    # smiles_clean and smiles_mol have the same atom ordering (RWMol copy
    # preserves indices), so we can look up map numbers from the original.
    smiles_match = smiles_matches[0]
    sdf_match    = sdf_matches[0]

    map_to_sdf = {}
    for smiles_idx, sdf_idx in zip(smiles_match, sdf_match):
        map_num = smiles_mol.GetAtomWithIdx(smiles_idx).GetAtomMapNum()
        if map_num:
            map_to_sdf[map_num] = sdf_idx

    return map_to_sdf


# ============================================================
# REACTION CONFIG GENERATION
# ============================================================

def build_reaction_config(protein_key, rxn_mapper=None):
    """
    Run the rxnmapper pipeline for every reaction in proteins[protein_key]
    and return a geometry config dict keyed by reaction id.

    Fragment -> ligand role assignment uses fingerprint similarity against
    compound SDFs and the protein's cofactors list.

    Atom map numbers are resolved to SDF indices via MCS alignment at
    config-build time, so analyze_nac() needs only a simple dict lookup
    at runtime — no SMARTS, no map-number search in the SDF.
    """
    if rxn_mapper is None:
        rxn_mapper = RXNMapper()

    protein_data   = proteins[protein_key]
    cofactor_names = set(protein_data.get("cofactors", []))
    compound_fps   = _build_compound_fps(compounds)   # built once, reused

    smiles_list = [r["smiles"] for r in protein_data["reactions"]]
    mapped_list = rxn_mapper.get_attention_guided_atom_maps(smiles_list)

    config = {}

    for rxn_entry, mapped_entry in zip(protein_data["reactions"], mapped_list):

        rxn_id     = rxn_entry["id"]
        mapped_rxn = mapped_entry["mapped_rxn"]

        result = analyze_reaction(mapped_rxn)

        if result["nac"] is None:
            continue

        nac       = result["nac"]
        mechanism = result["mechanism"]

        # Split reactant SMILES into per-fragment strings
        reactant_smiles = mapped_rxn.split(">>")[0]
        fragments       = reactant_smiles.split(".")
        frag_mols       = [Chem.MolFromSmiles(f) for f in fragments]

        # Assign cofactor/substrate role to each fragment
        fragment_roles = _assign_fragment_roles(
            fragments, protein_data, compound_fps, cofactor_names
        )

        def find_fragment_for_map(atom_map_num):
            """Return (fragment_index, mol) for the fragment containing this map number."""
            for i, mol in enumerate(frag_mols):
                if mol is None:
                    continue
                for a in mol.GetAtoms():
                    if a.GetAtomMapNum() == atom_map_num:
                        return i, mol
            return None, None

        def map_to_role(frag_idx):
            if frag_idx is None:
                return "unknown"
            if frag_idx < len(fragment_roles):
                return fragment_roles[frag_idx]
            return f"ligand{frag_idx}"

        # --------------------------------------------------------
        # Pre-compute MCS alignment for each fragment against its
        # compound SDF, producing map_num -> sdf_index tables.
        # We load the SDF for whichever compound matches this fragment.
        # --------------------------------------------------------
        def sdf_mol_for_role(role):
            """Load the SDF molecule for the compound assigned to this role."""
            # Find which compound best matches fragments with this role
            for frag_smiles, frag_role in zip(fragments, fragment_roles):
                if frag_role != role:
                    continue
                frag_mol = Chem.MolFromSmiles(frag_smiles)
                if frag_mol is None:
                    continue
                frag_fp = AllChem.GetMorganFingerprintAsBitVect(frag_mol, 2, 2048)
                best_name, best_sim = None, 0.0
                for name, ref_fp in compound_fps.items():
                    sim = DataStructs.TanimotoSimilarity(frag_fp, ref_fp)
                    if sim > best_sim:
                        best_sim, best_name = sim, name
                if best_name and best_sim >= 0.5:
                    sdf_path = (
                        f"{COMPOUNDS_PATH}/{best_name}/structure"
                        f"/{best_name}_{compounds[best_name]['pubchem_id']}.sdf"
                    )
                    try:
                        supplier = Chem.SDMolSupplier(str(sdf_path), removeHs=True)
                        return supplier[0] if len(supplier) > 0 else None
                    except Exception:
                        return None
            return None

        # Build one MCS table per unique role used in this reaction's NAC
        roles_needed = set()
        for dp in nac["distance_pairs"]:
            for map_num in (dp["atom1"], dp["atom2"]):
                frag_idx, _ = find_fragment_for_map(map_num)
                roles_needed.add(map_to_role(frag_idx))
        for ag in nac["angles"]:
            for key in ("atom1", "vertex", "atom2"):
                map_num = ag.get(key)
                if map_num is not None:
                    frag_idx, _ = find_fragment_for_map(map_num)
                    roles_needed.add(map_to_role(frag_idx))

        mcs_tables = {}  # role -> { map_num: sdf_idx }
        for role in roles_needed:
            if role in ("unknown",):
                continue
            frag_idx = next(
                (i for i, r in enumerate(fragment_roles) if r == role), None
            )
            if frag_idx is None:
                continue
            sdf_mol = sdf_mol_for_role(role)
            if sdf_mol is None:
                mcs_tables[role] = {}
                continue
            mcs_tables[role] = _build_map_to_sdf_index(fragments[frag_idx], sdf_mol)

        def make_atom_entry(map_num):
            """Return the atom descriptor dict for a given atom map number."""
            frag_idx, _ = find_fragment_for_map(map_num)
            role = map_to_role(frag_idx)
            sdf_idx = mcs_tables.get(role, {}).get(map_num)
            return {
                "map_num": map_num,
                "sdf_idx": sdf_idx,   # None if MCS failed; fallback at runtime
                "role":    role,
            }

        # --------------------------------------------------------
        # distance_pairs
        # --------------------------------------------------------
        distance_pairs = []
        for dp in nac["distance_pairs"]:
            a1 = make_atom_entry(dp["atom1"])
            a2 = make_atom_entry(dp["atom2"])
            distance_pairs.append({
                "name":    dp.get("label", f"dist_{dp['atom1']}_{dp['atom2']}"),
                "ligand1": a1["role"],
                "atom1_map_num": a1["map_num"],
                "atom1_sdf_idx": a1["sdf_idx"],
                "ligand2": a2["role"],
                "atom2_map_num": a2["map_num"],
                "atom2_sdf_idx": a2["sdf_idx"],
                "cutoff":  dp["cutoff"],
            })

        # --------------------------------------------------------
        # angles — skip degenerate cases (two of three atoms identical)
        # --------------------------------------------------------
        angle_configs = []
        for ag in nac["angles"]:
            a1_map     = ag["atom1"]
            vertex_map = ag.get("vertex", ag.get("atom2"))
            a3_map     = ag["atom2"] if "vertex" in ag else ag.get("atom3", ag["atom2"])

            if len({a1_map, vertex_map, a3_map}) < 3:
                continue

            a1 = make_atom_entry(a1_map)
            av = make_atom_entry(vertex_map)
            a3 = make_atom_entry(a3_map)

            angle_configs.append({
                "name":    ag.get("label", f"angle_{a1_map}_{vertex_map}_{a3_map}"),
                "ligand1": a1["role"],
                "atom1_map_num": a1["map_num"],
                "atom1_sdf_idx": a1["sdf_idx"],
                "ligand2": av["role"],
                "atom2_map_num": av["map_num"],
                "atom2_sdf_idx": av["sdf_idx"],
                "ligand3": a3["role"],
                "atom3_map_num": a3["map_num"],
                "atom3_sdf_idx": a3["sdf_idx"],
                "minimum": ag["minimum"],
                "ideal":   ag["ideal"],
            })
            
        # store the fragment SMILES keyed by role, for map-number transfer
        reactant_smiles_by_role = {}
        for frag_smiles, role in zip(fragments, fragment_roles):
            if role in ("cofactor", "substrate"):
                reactant_smiles_by_role[role] = frag_smiles

        config[rxn_id] = {
            "description":     rxn_entry.get("description", mechanism),
            "mechanism":       mechanism,
            "reactant_smiles": reactant_smiles_by_role,   # ← new
            "geometry": {
                "distance_pairs": distance_pairs,
                "angles":         angle_configs,
            },
        }

    return config


# ============================================================
# STRUCTURE LOADING
# ============================================================

def load_structure(structure_path):
    lower = structure_path.lower()
    if lower.endswith(".cif") or lower.endswith(".mmcif"):
        parser = MMCIFParser(QUIET=True)
    elif lower.endswith(".pdb"):
        parser = PDBParser(QUIET=True)
    else:
        raise ValueError(f"Unsupported structure format: {structure_path}")
    return parser.get_structure("complex", structure_path)


# ============================================================
# MATH
# ============================================================

def distance(a, b):
    return float(np.linalg.norm(a - b))


def angle_degrees(a, b, c):
    """Angle at vertex b, in degrees."""
    ba = a - b
    bc = c - b
    cos_a = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
    return float(np.degrees(np.arccos(np.clip(cos_a, -1.0, 1.0))))


# ============================================================
# ATOM RESOLUTION
# ============================================================

def resolve_atom(mol, rule, prefix):
    """
    Resolve an atom to an RDKit index using, in order:
      1. Precomputed SDF index from MCS alignment  (rule has f"{prefix}_sdf_idx")
      2. Direct index                               (rule has f"{prefix}_index")  legacy

    The map-number and SMARTS paths from previous versions are removed:
    map numbers are absent from PubChem SDFs, and SMARTS derived from
    heavy-H-annotated SMILES fail on clean SDF molecules.
    """
    sdf_key   = f"{prefix}_sdf_idx"
    index_key = f"{prefix}_index"

    if sdf_key in rule and rule[sdf_key] is not None:
        return rule[sdf_key]

    if index_key in rule:
        return rule[index_key]

    raise ValueError(
        f"Cannot resolve atom for rule '{rule.get('name', '?')}': "
        f"no {sdf_key} or {index_key} present. "
        f"MCS alignment may have failed for this ligand."
    )


# ============================================================
# LIGAND EXTRACTION  (structure -> coordinate list)
# ============================================================

def extract_residue_atoms(structure, residue_name, residue_number=None):
    atoms = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.get_resname().strip() != residue_name:
                    continue
                if residue_number is not None and residue.id[1] != residue_number:
                    continue
                atoms.extend(residue.get_atoms())

    if not atoms:
        raise ValueError(f"Residue not found: {residue_name}")
    return atoms


# ============================================================
# ATOM MAPPING  (RDKit index -> BioPython atom -> 3D coord)
# ============================================================

def _transfer_map_numbers(rdkit_mol, mapped_smiles_mol):
    """
    Match mapped_smiles_mol onto rdkit_mol by substructure (MCS-free,
    uses direct substructure match which works when they are the same
    compound) and copy atom map numbers from the SMILES mol to rdkit_mol.

    Returns True if the transfer succeeded, False otherwise.
    """
    # useChirality=False so stereo annotation differences don't block the match
    match = rdkit_mol.GetSubstructMatch(
        mapped_smiles_mol,
        useChirality=False,
    )
    if not match:
        # Try the other direction — sometimes atom count differs by Hs
        match2 = mapped_smiles_mol.GetSubstructMatch(
            rdkit_mol,
            useChirality=False,
        )
        if not match2:
            return False
        # match2[i] = index in rdkit_mol that corresponds to atom i in mapped_smiles_mol
        for smiles_idx, sdf_idx in enumerate(match2):
            map_num = mapped_smiles_mol.GetAtomWithIdx(smiles_idx).GetAtomMapNum()
            rdkit_mol.GetAtomWithIdx(sdf_idx).SetAtomMapNum(map_num)
        return True

    # match[i] = index in rdkit_mol that corresponds to atom i in mapped_smiles_mol
    for smiles_idx, sdf_idx in enumerate(match):
        map_num = mapped_smiles_mol.GetAtomWithIdx(smiles_idx).GetAtomMapNum()
        rdkit_mol.GetAtomWithIdx(sdf_idx).SetAtomMapNum(map_num)
    return True


def build_atom_mapping(rdkit_mol, structure_atoms, mapped_smiles_mol=None):
    """
    Map RDKit atom index → BioPython atom by positional order.

    Positional order is reliable when the SDF was used to generate the
    CIF (which is the case for Boltz inputs built from the same SDF).
    If atom counts don't match, raise immediately rather than silently
    producing wrong coordinates.

    If mapped_smiles_mol is provided, map numbers are transferred onto
    rdkit_mol before returning so that resolve_atom() can use them.
    """
    if rdkit_mol.GetNumAtoms() != len(structure_atoms):
        raise ValueError(
            f"Atom count mismatch: "
            f"RDKit={rdkit_mol.GetNumAtoms()} "
            f"Structure={len(structure_atoms)}"
        )

    if mapped_smiles_mol is not None:
        ok = _transfer_map_numbers(rdkit_mol, mapped_smiles_mol)
        if not ok:
            print(
                "Warning: could not transfer atom map numbers onto SDF mol "
                "(substructure match failed). SMARTS fallback will be used."
            )

    return {i: a for i, a in enumerate(structure_atoms)}


# ============================================================
# NAC SCORING  (continuous, not binary ratio)
# ============================================================

def _distance_score(value, cutoff):
    """
    Continuous score in [0, 1].
    1.0 if value <= cutoff; exponential decay beyond cutoff
    (half-score at ~cutoff + 0.35 A).
    """
    if value <= cutoff:
        return 1.0
    return float(np.exp(-2.0 * (value - cutoff)))


def _angle_score(value, minimum, ideal):
    """
    Continuous score in [0, 1].
    1.0 at or above ideal; linear decay to 0 at minimum; 0 below minimum.
    """
    if value >= ideal:
        return 1.0
    if value < minimum:
        return 0.0
    return (value - minimum) / (ideal - minimum)


# ============================================================
# MAIN ANALYSIS
# ============================================================

def analyze_nac(
    structure_path,
    reaction_name,
    reaction_config,
    cofactor_sdf,
    substrate_sdf,
    cofactor_resname,
    substrate_resname,
    cofactor_residue_number=None,
    substrate_residue_number=None,
):
    config    = reaction_config[reaction_name]
    structure = load_structure(structure_path)

    # --------------------------------------------------------
    # Load ligand molecules from SDF (heavy atoms only)
    # --------------------------------------------------------
    cofactor_mol  = Chem.MolFromMolFile(cofactor_sdf,  removeHs=True)
    substrate_mol = Chem.MolFromMolFile(substrate_sdf, removeHs=True)

    if cofactor_mol is None:
        raise ValueError("Failed to parse cofactor SDF")
    if substrate_mol is None:
        raise ValueError("Failed to parse substrate SDF")

    # --------------------------------------------------------
    # Extract 3D atoms from structure (heavy atoms only via BioPython)
    # --------------------------------------------------------
    cofactor_struct_atoms  = extract_residue_atoms(
        structure, cofactor_resname,  cofactor_residue_number
    )
    substrate_struct_atoms = extract_residue_atoms(
        structure, substrate_resname, substrate_residue_number
    )

    # --------------------------------------------------------
    # Build RDKit-index -> BioPython atom mappings (positional)
    # --------------------------------------------------------
    # --------------------------------------------------------
    # Transfer atom map numbers from the reaction SMILES onto
    # the SDF mols so resolve_atom() can use map numbers directly.
    # --------------------------------------------------------
    reactant_smiles = config.get("reactant_smiles", {})

    cofactor_smiles_mol  = None
    substrate_smiles_mol = None

    if "cofactor"  in reactant_smiles:
        cofactor_smiles_mol  = Chem.MolFromSmiles(reactant_smiles["cofactor"])
    if "substrate" in reactant_smiles:
        substrate_smiles_mol = Chem.MolFromSmiles(reactant_smiles["substrate"])

    cofactor_map  = build_atom_mapping(
        cofactor_mol,  cofactor_struct_atoms,  cofactor_smiles_mol
    )
    substrate_map = build_atom_mapping(
        substrate_mol, substrate_struct_atoms, substrate_smiles_mol
    )

    ligands = {
        "cofactor":  {"mol": cofactor_mol,  "map": cofactor_map},
        "substrate": {"mol": substrate_mol, "map": substrate_map},
    }

    def get_coord(ligand_key, rule, prefix):
        lig    = ligands[ligand_key]
        idx    = resolve_atom(lig["mol"], rule, prefix)
        bpatom = lig["map"].get(idx)
        if bpatom is None:
            raise ValueError(
                f"Atom index {idx} has no structure coordinate "
                f"(rule: {rule.get('name', '?')})"
            )
        return np.array(bpatom.coord)

    # --------------------------------------------------------
    # Distances
    # --------------------------------------------------------
    distance_results = []

    for rule in config["geometry"]["distance_pairs"]:
        coord1 = get_coord(rule["ligand1"], rule, "atom1")
        coord2 = get_coord(rule["ligand2"], rule, "atom2")

        d      = distance(coord1, coord2)
        passes = d <= rule["cutoff"]
        score  = _distance_score(d, rule["cutoff"])

        distance_results.append({
            "name":   rule["name"],
            "value":  round(d, 3),
            "cutoff": rule["cutoff"],
            "passes": passes,
            "score":  round(score, 3),
        })

    # --------------------------------------------------------
    # Angles  (atom2 is always the vertex)
    # --------------------------------------------------------
    angle_results = []

    for rule in config["geometry"]["angles"]:
        coord1 = get_coord(rule["ligand1"], rule, "atom1")
        coord2 = get_coord(rule["ligand2"], rule, "atom2")   # vertex
        coord3 = get_coord(rule["ligand3"], rule, "atom3")

        a      = angle_degrees(coord1, coord2, coord3)
        passes = a >= rule["minimum"]
        score  = _angle_score(a, rule["minimum"], rule["ideal"])

        angle_results.append({
            "name":    rule["name"],
            "value":   round(a, 2),
            "minimum": rule["minimum"],
            "ideal":   rule["ideal"],
            "passes":  passes,
            "score":   round(score, 3),
        })

    # --------------------------------------------------------
    # NAC score: mean of all continuous per-rule scores
    # --------------------------------------------------------
    all_scores = [r["score"] for r in distance_results + angle_results]
    nac_score  = float(np.mean(all_scores)) if all_scores else 0.0

    if   nac_score >= 0.85: interpretation = "Excellent catalytic geometry. Very strong NAC evidence."
    elif nac_score >= 0.60: interpretation = "Good catalytic geometry. Catalysis is plausible."
    elif nac_score >= 0.30: interpretation = "Partial NAC compatibility. Geometry is suboptimal."
    else:                   interpretation = "Poor catalytic geometry. Catalysis is unlikely."

    return {
        "reaction":         reaction_name,
        "description":      config["description"],
        "mechanism":        config.get("mechanism", "unknown"),
        "distance_results": distance_results,
        "angle_results":    angle_results,
        "nac_score":        round(nac_score, 4),
        "interpretation":   interpretation,
    }


# ============================================================
# MAIN
# ============================================================

if __name__ == "__main__":

    protein      = "sard4"
    ref_organism = "tabaco"
    ooi          = "tucunaca"
    reaction     = "r1"
    substrate    = "tetrahydroharmol"
    cofactor     = "nadp+"

    rxn_mapper      = RXNMapper()
    reaction_config = build_reaction_config(protein, rxn_mapper)

    boltz_path         = f"{PROTEINS_PATH}/{protein}/analysis/structural/boltz"
    structure_filename = "sample_0_predicted_structure.cif"

    if ooi:
        structure_path = (
            f"{boltz_path}/{ooi}/{ref_organism}"
            f"/{protein}_{ooi}_{ref_organism}_{reaction}/{structure_filename}"
        )
    else:
        structure_path = (
            f"{boltz_path}/{ref_organism}"
            f"/{protein}_{ref_organism}_{reaction}/{structure_filename}"
        )

    results = analyze_nac(
        reaction_name=reaction,
        reaction_config=reaction_config,
        structure_path=structure_path,
        substrate_sdf=f"{COMPOUNDS_PATH}/{substrate}/structure/{substrate}_{compounds[substrate]['pubchem_id']}.sdf",
        cofactor_sdf=f"{COMPOUNDS_PATH}/{cofactor}/structure/{cofactor}_{compounds[cofactor]['pubchem_id']}.sdf",
        substrate_resname="LIG1",
        cofactor_resname="LIG2",
    )

    print(json.dumps(results, indent=4))