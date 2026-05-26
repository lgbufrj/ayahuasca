"""Utility to build a stable atom-map -> original-atom index bridge.

The function implemented here is intentionally defensive: RXNMapper and similar
tools annotate RDKit atoms with temporary atom-map numbers (atom.GetAtomMapNum()).
During mapping we should preserve the original atom indices by setting an
integer atom property (for example via atom.SetIntProp("original_idx", idx)).

This module exposes `build_atom_map_index` which consumes a list of mapped
ligand molecules and (optionally) a lookup of original ligand molecules and
returns a dictionary keyed by atom-map number with the canonical bridge data.

Expected input shape for `mapped_molecules` is one of:

1) list of dicts: [{"compound": "nadp+", "mol": rdkit.Chem.Mol}, ...]
2) list of tuples: [("nadp+", rdkit.Chem.Mol), ...]
3) list of rdkit.Chem.Mol (in which case the caller must ensure atoms carry a
   "compound" property or the function will use a positional name)

The function will look for the following atom properties (in order) to find
the original atom index saved before mapping:
 - IntProp: "original_idx"
 - IntProp: "original_atom_idx"
 - IntProp: "source_atom_idx"
 - String prop with the same names (will be cast to int)

If none is present it falls back to the atom's current GetIdx().

Returned dict example:

{
    63: {
        "compound": "nadp+",
        "mapped_atom_idx": 5,
        "original_atom_idx": 14,
        "element": "C"
    },
    ...
}

This file intentionally avoids forcing heavy dependencies (like gemmi) and
assumes RDKit is available in the environment that uses it. If RDKit is not
installed an ImportError will be raised.

"""
from typing import List, Dict, Any, Union, Tuple

try:
    from rdkit import Chem
except Exception as e:
    raise ImportError("RDKit is required by chemistry.structures.atom_map: " + str(e))


def _extract_original_idx(atom: Chem.Atom) -> int:
    """Try several common property names and return an int index or None."""
    # Try RDKit IntProp first (faster / canonical API)
    for intprop in ("original_idx", "original_atom_idx", "source_atom_idx"):
        try:
            if atom.HasProp(intprop):
                v = atom.GetProp(intprop)
                try:
                    return int(v)
                except Exception:
                    # fall through to try GetIntProp
                    pass
        except Exception:
            # some RDKit versions may raise when prop missing; ignore
            pass

        try:
            # GetIntProp raises if prop missing; catch and continue
            return atom.GetIntProp(intprop)
        except Exception:
            pass

    # Nothing found
    return None


def build_atom_map_index(
    mapped_molecules: List[Union[Dict[str, Any], Tuple[str, Chem.Mol], Chem.Mol]],
    original_ligands: Dict[str, Chem.Mol] = None,
    substrates: List[str] = None,
) -> Dict[int, Dict[str, Any]]:
    """
    Build a mapping from reaction atom-map numbers to original ligand atom
    indices and basic atom metadata.

    Parameters
    ----------
    mapped_molecules
        A list containing mapped ligand molecules. Each item can be either:
        - a dict with keys "compound" and "mol" (recommended),
        - a tuple (compound_name, mol), or
        - a raw RDKit Mol (in which case the compound name will be "mol_{i}")

    original_ligands
        Optional dict mapping compound name -> original (unmapped) RDKit Mol.
        This is only used for reference; the preferred way to preserve the
        original atom index is to save it on each atom before mapping using
        atom.SetIntProp("original_idx", idx).

    substrates
        Optional list of compound names (strings). If provided, only entries
        whose compound name matches an entry in `substrates` (case-insensitive)
        will be included in the returned index. This lets callers restrict the
        mapping to substrates of interest and skip product molecules.

    Returns
    -------
    dict
        Keys are atom-map integers (atom.GetAtomMapNum()) and values are dicts
        with at least: compound, mapped_atom_idx, original_atom_idx, element.

    Notes
    -----
    This function is deliberately permissive: it will not fail if some atoms
    lack an atom-map number; those are skipped. If no saved original index is
    found, it falls back to the atom's current GetIdx().
    """

    # normalize substrate names for quick membership tests and remove
    # proton-like names (e.g. 'h+', 'h') which we don't want to map.
    if substrates:
        _proton_aliases = {"h", "h+", "h1+", "1h+", "[1h+]", "proton"}
        substrates_set = set()
        for s in substrates:
            try:
                s_norm = s.lower().strip()
            except Exception:
                continue
            if s_norm in _proton_aliases:
                # skip proton-equivalent substrate names
                continue
            substrates_set.add(s_norm)
    else:
        substrates_set = None

    index: Dict[int, Dict[str, Any]] = {}

    for i, entry in enumerate(mapped_molecules):
        if isinstance(entry, dict):
            compound = entry.get("compound")
            mol = entry.get("mol")
        elif isinstance(entry, tuple) and len(entry) == 2:
            compound, mol = entry
        else:
            # assume raw RDKit Mol
            mol = entry
            compound = f"mol_{i}"

        if mol is None:
            continue

        # if substrates filter provided, only include matching compounds
        comp_name = compound if compound is not None else f"mol_{i}"
        comp_norm = comp_name.lower().strip()
        if substrates_set is not None and comp_norm not in substrates_set:
            # skip products / any non-substrate
            continue

        if not isinstance(mol, Chem.Mol):
            raise TypeError(f"Expected RDKit Mol for compound {compound}, got {type(mol)}")

        for atom in mol.GetAtoms():
            mapnum = atom.GetAtomMapNum()
            if not mapnum:
                # unmapped atom, skip
                continue

            original_idx = _extract_original_idx(atom)

            if original_idx is None:
                # final fallback: use atom.GetIdx()
                original_idx = atom.GetIdx()

            index[mapnum] = {
                "compound": compound,
                "mapped_atom_idx": atom.GetIdx(),
                "original_atom_idx": original_idx,
                "element": atom.GetSymbol(),
            }

    return index
