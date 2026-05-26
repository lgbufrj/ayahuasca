from typing import List, Dict, Any, Optional
import os
import json
from tqdm import tqdm

try:
    from rdkit import Chem
    RDKit_AVAILABLE = True
except Exception:
    RDKit_AVAILABLE = False

from rxnmapper import RXNMapper

from data import proteins
import data

from chemistry.reactions.analysis import analyze_reaction

from chemistry.structures.atom_map import build_atom_map_index
from chemistry.structures.pdb_coords import (
    load_structure_atoms,
    align_atommap_and_compute_distances,
)
from chemistry.structures.calculate_geometry import compute_transfer_geometry


def _parse_mapped_rxn_to_mols(mapped_rxn: str) -> List[Any]:
    mols = []
    if not RDKit_AVAILABLE:
        return mols
    text = mapped_rxn or ""
    fragments = []
    if '>>' in text:
        left, right = text.split('>>', 1)
        fragments += left.split('.')
        fragments += right.split('.')
    else:
        parts = text.split('>')
        for p in parts:
            fragments += p.split('.')

    for frag in fragments:
        s = frag.strip()
        if not s:
            continue
        try:
            m = Chem.MolFromSmiles(s)
            if m is not None:
                mols.append(m)
        except Exception:
            continue
    return mols


def _build_original_ligands_map(mapped_entries: List[Dict[str, Any]]) -> Optional[Dict[str, Any]]:
    if not RDKit_AVAILABLE:
        return None
    try:
        original_ligands = {}
        for me in mapped_entries:
            cname = me.get('compound')
            mol = me.get('mol')
            if mol is None:
                continue
            mcopy = Chem.Mol(mol)
            for a in mcopy.GetAtoms():
                try:
                    a.SetAtomMapNum(0)
                except Exception:
                    pass
            original_ligands[cname] = mcopy
        return original_ligands
    except Exception:
        return None


def run_pipeline(protein: str, outdir: Optional[str] = None, ref_organism: str = 'tabaco') -> List[str]:
    """Run pipeline for a single protein.

    Returns list of output file paths written.
    """
    rxn_mapper = RXNMapper()
    reactions = proteins[protein]["reactions"]
    rxns = [r["smiles"] for r in reactions]
    mapped_rxns = rxn_mapper.get_attention_guided_atom_maps(rxns)

    if outdir is None:
        outdir = os.path.join(os.getcwd(), 'outputs', protein)
    os.makedirs(outdir, exist_ok=True)

    out_paths = []
    for rxn_data, mapped in tqdm(list(zip(reactions, mapped_rxns)), desc=f"Processing {protein}", total=len(reactions)):
        rxn_id = rxn_data.get('id') or 'unknown'
        mapped_rxn = mapped.get('mapped_rxn') if isinstance(mapped, dict) else mapped

        results = analyze_reaction(
            mapped_rxn=mapped_rxn,
            reaction_metadata=rxn_data,
            cofactors=proteins[protein].get("cofactors")
        )

        atom_map_index = None
        pdb_mapping = None
        geom = None
        mapped_entries = []

        # Build atom_map_index if possible
        if RDKit_AVAILABLE and mapped_rxn:
            mols = _parse_mapped_rxn_to_mols(mapped_rxn)

            substrates = rxn_data.get('substrates', []) or []
            products = rxn_data.get('products', []) or []
            compound_order = list(substrates) + list(products)

            # filter tiny fragments
            filtered = []
            for mol in mols:
                try:
                    na = mol.GetNumAtoms()
                except Exception:
                    na = 0
                has_heavy = any(a.GetSymbol() != 'H' for a in mol.GetAtoms())
                if na <= 1 and not has_heavy:
                    continue
                filtered.append(mol)

            mapped_entries = []
            for mi, mol in enumerate(filtered):
                cname = compound_order[mi] if mi < len(compound_order) else f"mol_{mi}"
                mapped_entries.append({"compound": cname, "mol": mol})

            try:
                atom_map_index = build_atom_map_index(mapped_entries, substrates=substrates)
            except Exception as e:
                atom_map_index = {"error": str(e)}

        # Align to structure when available
        if atom_map_index:
            struct_path = rxn_data.get('structure_file') or rxn_data.get('structure') or proteins[protein].get('structure_file') or proteins[protein].get('structure')

            if not struct_path:
                base = getattr(data, 'PROTEINS_PATH', None)
                if base:
                    candidate = os.path.join(
                        base,
                        protein,
                        'analysis',
                        'structural',
                        'boltz',
                        ref_organism,
                        f"{protein}_{ref_organism}_{rxn_id}",
                        'sample_0_predicted_structure.cif',
                    )
                    if os.path.exists(candidate):
                        struct_path = candidate

            if struct_path and os.path.exists(struct_path):
                try:
                    pdb_idx = load_structure_atoms(struct_path)
                    original_ligands = _build_original_ligands_map(mapped_entries)
                    ligand_to_residue = rxn_data.get('ligand_to_residue') if isinstance(rxn_data, dict) else None

                    pdb_mapping = align_atommap_and_compute_distances(
                        atom_map_index,
                        pdb_idx,
                        transfer_vector=results.get('transfer_vector'),
                        original_ligands=original_ligands,
                        ligand_to_residue=ligand_to_residue,
                    )

                    try:
                        geom = compute_transfer_geometry(pdb_mapping.get('pdb_mapped', {}), results.get('transfer_vector'))
                    except Exception:
                        geom = None
                except Exception as e:
                    pdb_mapping = {"error": str(e)}
            else:
                pdb_mapping = None

        payload = {
            "protein": protein,
            "reaction_id": rxn_id,
            "mapped_rxn": mapped_rxn,
            "mechanism": results.get("mechanism"),
            "reaction_center": results.get("reaction_center"),
            "transfer_vector": results.get("transfer_vector"),
            "bond_changes": results.get("bond_changes"),
            "hydrogen_changes": results.get("hydrogen_changes"),
            "atom_map_index": atom_map_index,
            "pdb_mapping": pdb_mapping,
            "geometry": geom,
            "mapped_entries": None,
            "original_ligand_smiles": None,
        }

        try:
            if mapped_entries:
                out_mapped = []
                ols = {}
                if RDKit_AVAILABLE:
                    for me in mapped_entries:
                        cname = me.get('compound')
                        mol = me.get('mol')
                        s = None
                        try:
                            mcopy = Chem.Mol(mol)
                            for a in mcopy.GetAtoms():
                                try:
                                    a.SetAtomMapNum(0)
                                except Exception:
                                    pass
                            s = Chem.MolToSmiles(mcopy, isomericSmiles=True)
                        except Exception:
                            s = None
                        out_mapped.append({'compound': cname, 'smiles': s})
                        ols[cname] = s
                else:
                    for me in mapped_entries:
                        out_mapped.append({'compound': me.get('compound'), 'smiles': None})
                payload['mapped_entries'] = out_mapped
                payload['original_ligand_smiles'] = ols
        except Exception:
            pass

        outpath = os.path.join(outdir, f"{protein}_{rxn_id}.json")
        with open(outpath, 'w') as fh:
            json.dump(payload, fh, indent=2)
        out_paths.append(outpath)

    return out_paths
