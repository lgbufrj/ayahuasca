"""Simple PDB / mmCIF ligand atom coordinate parser.

This module exposes a small, dependency-light API to extract ligand atom
coordinates from PDB or mmCIF files. It tries to use gemmi if available for
robust parsing; otherwise it falls back to lightweight text parsers that
handle common mmCIF "atom_site" loops and PDB "HETATM"/"ATOM" records.

Provided functions
- load_structure_atoms(path) -> dict of residue-keys -> atom dicts
- find_ligand_atoms(atoms_index, resname, chain=None, resnum=None) -> dict

Residue-key format returned by load_structure_atoms is a tuple:
    (resname, chain_id, resnum, auth_seq_id)

Atom dicts map atom_name -> {
    "coord": (x, y, z),
    "element": element,
    "serial": int or None,
    "altloc": str or None,
}

The fallback mmCIF parser implemented here is minimal and will work for
standard mmCIF files that present atom_site as a simple loop with whitespace
separated columns. If your mmCIF files contain quoted fields or use
non-standard formatting, install `gemmi` (`pip install gemmi`) for robust
parsing.
"""
from typing import Dict, Tuple, Any, Optional
import os
import math
from collections import defaultdict

try:
    # RDKit optional dependency for degree/neighborhood heuristics
    from rdkit import Chem
    _RDKit_INSTALLED = True
except Exception:
    _RDKit_INSTALLED = False


def _try_gemmi_load(path: str) -> Optional[Dict[Tuple[str, str, str, str], Dict[str, Any]]]:
    try:
        import gemmi
    except Exception:
        return None

    struct = gemmi.read_structure(path)
    result: Dict[Tuple[str, str, str, str], Dict[str, Any]] = {}

    def _normalize_element(obj) -> Optional[str]:
        # obj may be a gemmi.Element, an object, or a string like "GEMMI.ELEMENT('N')"
        try:
            s = str(obj)
        except Exception:
            return None
        if not s:
            return None
        # common pattern GEMMI.ELEMENT('N') -> extract inside quotes
        import re
        m = re.search(r"'([A-Za-z]{1,2})'", s)
        if m:
            return m.group(1).upper()
        # fallback: take first element-like token
        m2 = re.search(r"([A-Z][a-z]?)", s)
        if m2:
            return m2.group(1).upper()
        return None

    # Traverse models -> chains -> residues -> atoms
    for model in struct:
        for chain in model:
            chain_id = chain.name
            for residue in chain:
                try:
                    resname = residue.name
                except Exception:
                    resname = ''
                try:
                    resnum = residue.seqid.num
                except Exception:
                    resnum = ''
                try:
                    auth_seq = residue.seqid.str()
                except Exception:
                    auth_seq = ''

                key = (resname or '', chain_id or '', str(resnum) if resnum is not None else '', auth_seq or '')
                atoms: Dict[str, Any] = {}
                for atom in residue:
                    try:
                        aname = atom.name
                        pos = atom.pos
                        coord = (pos.x, pos.y, pos.z)
                        # gemmi.Atom.element may be an Element object or a string
                        elt = _normalize_element(atom.element)
                        serial = int(atom.serial) if getattr(atom, 'serial', 0) else None
                        altloc = getattr(atom, 'altloc', None)
                        atoms[aname] = {
                            'coord': coord,
                            'element': elt,
                            'serial': serial,
                            'altloc': altloc,
                        }
                    except Exception:
                        # skip problematic atom entries
                        continue

                if atoms:
                    result[key] = atoms

    return result


def _parse_pdb(path: str) -> Dict[Tuple[str, str, str, str], Dict[str, Any]]:
    result: Dict[Tuple[str, str, str, str], Dict[str, Any]] = {}
    with open(path, 'r') as fh:
        for line in fh:
            if not (line.startswith('ATOM') or line.startswith('HETATM')):
                continue
            # PDB fixed-column fields (1-indexed in spec)
            # Columns: 7-11 serial, 13-16 name, 17 altLoc, 18-20 resName, 22 chainID,
            # 23-26 resSeq, 31-38 x, 39-46 y, 47-54 z, 77-78 element
            try:
                serial = int(line[6:11].strip())
            except Exception:
                serial = None
            name = line[12:16].strip()
            altloc = line[16].strip() or None
            resname = line[17:20].strip()
            chain = line[21].strip() or None
            resseq = line[22:26].strip()
            auth_seq = line[22:26].strip()
            try:
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
            except Exception:
                # malformed coordinate line, skip
                continue
            element = line[76:78].strip() or None

            key = (resname, chain or '', resseq or '', auth_seq or '')
            atoms = result.setdefault(key, {})
            atoms[name] = {
                "coord": (x, y, z),
                "element": element,
                "serial": serial,
                "altloc": altloc,
            }
    return result


def _parse_mmcif(path: str) -> Dict[Tuple[str, str, str, str], Dict[str, Any]]:
    # Minimal mmCIF atom_site loop parser. Works for simple, whitespace
    # separated atom_site loops. Not a full mmCIF implementation.
    fields = []
    data_started = False
    result: Dict[Tuple[str, str, str, str], Dict[str, Any]] = {}

    with open(path, 'r') as fh:
        for line in fh:
            s = line.strip()
            if not s:
                # blank line resets state
                if data_started:
                    break
                continue
            if s.startswith('loop_'):
                # start collecting headers
                fields = []
                data_started = False
                continue
            if s.startswith('_atom_site.'):
                fields.append(s)
                continue
            if fields and not s.startswith('_'):
                # first data line after headers
                data_started = True
            if data_started:
                # split respecting quoted tokens would be ideal; we attempt a
                # simple whitespace split which is fine for typical atom_site
                parts = s.split()
                if len(parts) < len(fields):
                    # possible multi-line or quoted fields; skip this entry
                    continue
                row = dict(zip(fields, parts))
                # extract common mmCIF atom_site fields
                resname = row.get('_atom_site.label_comp_id') or row.get('_atom_site.auth_comp_id')
                chain = row.get('_atom_site.label_asym_id') or row.get('_atom_site.auth_asym_id') or ''
                resseq = row.get('_atom_site.label_seq_id') or row.get('_atom_site.auth_seq_id') or ''
                aname = row.get('_atom_site.label_atom_id') or row.get('_atom_site.auth_atom_id')
                try:
                    x = float(row.get('_atom_site.Cartn_x', row.get('_atom_site.Cartn_x')))
                    y = float(row.get('_atom_site.Cartn_y', row.get('_atom_site.Cartn_y')))
                    z = float(row.get('_atom_site.Cartn_z', row.get('_atom_site.Cartn_z')))
                except Exception:
                    continue
                element = row.get('_atom_site.type_symbol')
                serial = None
                key = (resname or '', chain or '', resseq or '', resseq or '')
                atoms = result.setdefault(key, {})
                atoms[aname] = {
                    'coord': (x, y, z),
                    'element': element,
                    'serial': serial,
                    'altloc': None,
                }
    return result


def load_structure_atoms(path: str) -> Dict[Tuple[str, str, str, str], Dict[str, Any]]:
    """Load atom coordinates from a PDB or mmCIF file.

    Returns a dict keyed by (resname, chain, resnum, auth_seq_id) -> atom dicts
    mapping atom_name -> {coord, element, serial, altloc}.
    """
    if not os.path.exists(path):
        raise FileNotFoundError(path)

    # prefer gemmi if installed
    gemmi_result = _try_gemmi_load(path)
    if gemmi_result is not None:
        return gemmi_result

    # fallback based on extension
    lower = path.lower()
    if lower.endswith('.pdb'):
        return _parse_pdb(path)

    if lower.endswith('.cif') or lower.endswith('.mmcif'):
        parsed = _parse_mmcif(path)
        if parsed:
            return parsed
        raise RuntimeError('Failed to parse mmCIF file without gemmi; install gemmi for full support')

    # try PDB parser by default
    return _parse_pdb(path)


def find_ligand_atoms(atoms_index: Dict[Tuple[str, str, str, str], Dict[str, Any]],
                      resname: str,
                      chain: Optional[str] = None,
                      resnum: Optional[str] = None) -> Dict[str, Any]:
    """Return atoms for matching residue name (case-insensitive). If chain or
    resnum are provided they restrict the search.
    """
    resname_norm = (resname or '').strip().upper()
    out = {}
    for (rname, chain_id, rseq, auth_seq), atoms in atoms_index.items():
        if rname.upper() != resname_norm:
            continue
        if chain is not None and (chain_id or '').strip() != (chain or '').strip():
            continue
        if resnum is not None and (rseq or '').strip() != (resnum or '').strip():
            continue
        # merge atoms from matching residues; if multiple residues match the
        # returned dict may contain multiple residues' atoms (caller should
        # disambiguate by chain/resnum when possible)
        for k, v in atoms.items():
            out[k] = v
    return out


def _estimate_neighbor_counts(pdb_atoms: Dict[str, Any], threshold: float = 1.9) -> Dict[str, int]:
    """Estimate covalent neighbor counts for atoms in a PDB ligand by simple
    distance thresholding. Returns a dict atom_name -> neighbor_count.

    This is a heuristic used to disambiguate atoms by local connectivity when
    no explicit bond information is available in the coordinate file.
    """
    names = list(pdb_atoms.keys())
    coords = [pdb_atoms[n]['coord'] for n in names]
    neigh = {n: 0 for n in names}
    for i, a in enumerate(coords):
        for j, b in enumerate(coords):
            if i == j:
                continue
            dx = a[0] - b[0]
            dy = a[1] - b[1]
            dz = a[2] - b[2]
            d2 = dx * dx + dy * dy + dz * dz
            if d2 <= threshold * threshold:
                neigh[names[i]] += 1
    return neigh


def _residue_element_counts(pdb_atoms: Dict[str, Any]) -> Dict[str, int]:
    """Return a simple element -> count signature for a residue's atoms.

    Elements are upper-cased; missing element fields are ignored.
    """
    from collections import defaultdict
    counts = defaultdict(int)
    for info in pdb_atoms.values():
        elt = (info.get('element') or '').strip().upper()
        if not elt:
            continue
        counts[elt] += 1
    return dict(counts)


def _mol_element_counts(mol: 'Chem.Mol') -> Dict[str, int]:
    """Return element counts for an RDKit Mol (heavy atoms only).

    If RDKit is not available the function returns an empty dict.
    """
    if not _RDKit_INSTALLED or mol is None:
        return {}
    from collections import defaultdict
    counts = defaultdict(int)
    try:
        for a in mol.GetAtoms():
            sym = a.GetSymbol().upper()
            if sym == 'H':
                continue
            counts[sym] += 1
    except Exception:
        return {}
    return dict(counts)


def map_atommap_to_pdb_coords(
    atom_map_index: Dict[int, Dict[str, Any]],
    pdb_atoms_index: Dict[Tuple[str, str, str, str], Dict[str, Any]],
    original_ligands: Optional[Dict[str, 'Chem.Mol']] = None,
    ligand_to_residue: Optional[Dict[str, Tuple[str, Optional[str], Optional[str]]]] = None,
    distance_threshold: float = 1.9,
) -> Dict[int, Dict[str, Any]]:
    """Attempt to map reaction atom-map integers to PDB atom records (coords).

    Parameters
    - atom_map_index: mapping produced by `build_atom_map_index` (mapnum -> metadata)
    - pdb_atoms_index: output of `load_structure_atoms`
    - original_ligands: optional dict compound -> RDKit Mol (unmapped) used to
      compute RDKit atom degrees for better matching.
    - ligand_to_residue: optional mapping compound -> (resname, chain, resnum)
      to restrict which residue in the structure corresponds to the compound.

    Returns a dict keyed by atom_map integer with values containing at least:
    {"compound", "pdb_atom_name", "serial", "coord", "element"}

    Notes
    """
    result = {}

    # precompute per-compound element counts from atom_map_index as a
    # fallback when RDKit is not available
    comp_element_counts = {}
    for meta in atom_map_index.values():
        comp = meta.get('compound')
        elt = (meta.get('element') or '').strip().upper()
        if not comp:
            continue
        cc = comp_element_counts.setdefault(comp, {})
        if elt:
            cc[elt] = cc.get(elt, 0) + 1

    # Precompute a lookup of residues by simple keys for heuristic matching
    residues_by_name = defaultdict(list)
    for (resname, chain, rseq, auth_seq), atoms in pdb_atoms_index.items():
        residues_by_name[resname.upper()].append(((resname, chain, rseq, auth_seq), atoms))

    for mapnum, meta in atom_map_index.items():
        compound = meta.get('compound')
        element = (meta.get('element') or '').strip()
        orig_idx = meta.get('original_atom_idx')

        # find candidate residue(s) in pdb
        candidate_atoms = {}
        if ligand_to_residue and compound in ligand_to_residue:
            resname, chain, resnum = ligand_to_residue[compound]
            candidate_atoms = find_ligand_atoms(pdb_atoms_index, resname, chain, resnum)
        else:
            # try matching by residue name substring heuristics
            comp_norm = (compound or '').upper()
            # exact match first
            if comp_norm in residues_by_name:
                # choose first residue by default
                candidate_atoms = residues_by_name[comp_norm][0][1]
            else:
                # substring match
                for rname, lst in residues_by_name.items():
                    if comp_norm in rname or rname in comp_norm:
                        candidate_atoms = lst[0][1]
                        break

        # If no candidate atoms by name heuristics, try element-count fallback
        if not candidate_atoms:
            # preferred: use RDKit-derived molecule element counts
            target_counts = None
            if original_ligands and compound in original_ligands:
                target_counts = _mol_element_counts(original_ligands[compound])
            # fallback: use atom_map_index aggregated element counts
            if not target_counts:
                target_counts = comp_element_counts.get(compound) or {}

            if target_counts:
                # build residue signatures and pick the one with smallest
                # difference in counts (L1 distance)
                best = None
                best_score = None
                for (rkey, atoms) in pdb_atoms_index.items():
                    rc = _residue_element_counts(atoms)
                    # compute L1 distance over elements in target_counts
                    score = 0
                    for el, cnt in target_counts.items():
                        score += abs(cnt - rc.get(el, 0))
                    # penalize residues with many extra heavy atoms
                    extra = sum(v for k, v in rc.items() if k not in target_counts)
                    score += extra
                    if best_score is None or score < best_score:
                        best_score = score
                        best = atoms
                if best is not None:
                    candidate_atoms = best

        if not candidate_atoms:
            # nothing found, skip but register empty entry
            result[mapnum] = {
                'compound': compound,
                'pdb_atom_name': None,
                'serial': None,
                'coord': None,
                'element': element,
                'reason': 'no_candidate_residue'
            }
            continue

        # narrow candidates by element
        candidates = [ (name, info) for name, info in candidate_atoms.items() if (info.get('element') or '').strip().upper() == element.upper() ]
        if not candidates:
            # fallback: accept any element (some files lack element column)
            candidates = list(candidate_atoms.items())

        # if we have an RDKit mol for this compound and an original index, use
        # RDKit degree (heavy neighbor count) to disambiguate
        chosen = None
        if original_ligands and compound in original_ligands and orig_idx is not None and _RDKit_INSTALLED:
            try:
                rdm = original_ligands[compound]
                r_atom = rdm.GetAtomWithIdx(int(orig_idx))
                # heavy neighbor count in RDKit
                r_deg = sum(1 for n in r_atom.GetNeighbors() if n.GetSymbol() != 'H')
            except Exception:
                r_deg = None

            if r_deg is not None:
                # estimate neighbor counts in pdb candidates
                pdb_neigh = _estimate_neighbor_counts(candidate_atoms, threshold=distance_threshold)
                # filter candidates by neighbor count equal to r_deg
                cand_by_deg = [ (n, candidate_atoms[n]) for n in candidate_atoms.keys() if pdb_neigh.get(n) == r_deg and ((candidate_atoms[n].get('element') or '').strip().upper() == element.upper()) ]
                if len(cand_by_deg) == 1:
                    chosen = cand_by_deg[0]

        # if still ambiguous, try to use neighbor counts on filtered candidates
        if chosen is None:
            if len(candidates) == 1:
                chosen = candidates[0]
            else:
                # compute neighbor counts for only the candidate set
                sub_atoms = {n: candidate_atoms[n] for n, _ in candidates}
                pdb_neigh = _estimate_neighbor_counts(sub_atoms, threshold=distance_threshold)
                # try to find a unique candidate where neighbor count matches the
                # most common neighbor count for this element in the set.
                deg_counter = defaultdict(int)
                for n in sub_atoms:
                    deg_counter[pdb_neigh.get(n, 0)] += 1
                # choose candidate whose degree is rarest (heuristic)
                rare_deg = min(deg_counter.items(), key=lambda x: x[1])[0]
                rare_candidates = [n for n in sub_atoms.keys() if pdb_neigh.get(n, 0) == rare_deg]
                if len(rare_candidates) == 1:
                    chosen = (rare_candidates[0], sub_atoms[rare_candidates[0]])
                else:
                    # final fallback: pick first candidate
                    chosen = candidates[0] if candidates else None

        if chosen is None:
            result[mapnum] = {
                'compound': compound,
                'pdb_atom_name': None,
                'serial': None,
                'coord': None,
                'element': element,
                'reason': 'no_candidate_after_filters'
            }
            continue

        aname, info = chosen
        result[mapnum] = {
            'compound': compound,
            'pdb_atom_name': aname,
            'serial': info.get('serial'),
            'coord': tuple(info.get('coord')) if info.get('coord') else None,
            'element': info.get('element') or element,
        }

    return result


def align_atommap_and_compute_distances(
    atom_map_index: Dict[int, Dict[str, Any]],
    pdb_atoms_index: Dict[Tuple[str, str, str, str], Dict[str, Any]],
    transfer_vector: Dict[str, Any] = None,
    original_ligands: Optional[Dict[str, 'Chem.Mol']] = None,
    ligand_to_residue: Optional[Dict[str, Tuple[str, Optional[str], Optional[str]]]] = None,
    distance_threshold: float = 1.9,
) -> Dict[str, Any]:
    """Align atom-map -> PDB coords and compute donor-acceptor distance.

    Parameters
    - atom_map_index: output from build_atom_map_index
    - pdb_atoms_index: output from load_structure_atoms
    - transfer_vector: optional dict containing donor/acceptor atom_map ints
    - original_ligands: optional compound->RDKit Mol mapping for better matching
    - ligand_to_residue: optional compound->(resname,chain,resnum) mapping

    Returns a dict with keys:
    - pdb_mapped: mapnum -> pdb mapping (as returned by map_atommap_to_pdb_coords)
    - compounds: compound -> list of entries {mapnum, pdb_atom_name, coord, element}
    - transfer_distance: float or None
    - donor: donor mapping entry or None
    - acceptor: acceptor mapping entry or None
    - missing: list of atom_map ints with no pdb coord
    """
    mapped = map_atommap_to_pdb_coords(
        atom_map_index, pdb_atoms_index, original_ligands=original_ligands,
        ligand_to_residue=ligand_to_residue, distance_threshold=distance_threshold
    )

    compounds = defaultdict(list)
    missing = []
    for mapnum, info in mapped.items():
        comp = info.get('compound')
        coord = info.get('coord')
        entry = {
            'mapnum': mapnum,
            'compound': comp,
            'pdb_atom_name': info.get('pdb_atom_name'),
            'serial': info.get('serial'),
            'coord': coord,
            'element': info.get('element'),
        }
        if coord is None:
            missing.append(mapnum)
        compounds[comp].append(entry)

    donor_entry = None
    acceptor_entry = None
    transfer_distance = None
    if transfer_vector:
        try:
            dmap = transfer_vector.get('donor', {}).get('atom_map')
            amap = transfer_vector.get('acceptor', {}).get('atom_map')
        except Exception:
            dmap = amap = None

        if dmap is not None:
            donor_entry = mapped.get(dmap)
        if amap is not None:
            acceptor_entry = mapped.get(amap)

        # compute Euclidean distance if both coords available
        if donor_entry and acceptor_entry and donor_entry.get('coord') and acceptor_entry.get('coord'):
            a = donor_entry['coord']
            b = acceptor_entry['coord']
            try:
                dx = a[0] - b[0]
                dy = a[1] - b[1]
                dz = a[2] - b[2]
                transfer_distance = math.sqrt(dx * dx + dy * dy + dz * dz)
            except Exception:
                transfer_distance = None

    return {
        'pdb_mapped': mapped,
        'compounds': dict(compounds),
        'transfer_distance': transfer_distance,
        'donor': donor_entry,
        'acceptor': acceptor_entry,
        'missing': missing,
    }
