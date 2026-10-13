#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
boltz_reaction_distances.py
===========================

Distance between the two *reacting* atoms of docked substrates in Boltz-2 models.

Inputs
------
* ``--yaml``  pathway YAML (``compounds`` / ``organisms`` / ``proteins`` -> ``reactions``),
  where every reaction has a reaction SMILES/SMARTS (``smarts`` is used, ``smiles`` as fallback).
* ``--cif-dir`` folder with Boltz-2 mmCIF predictions (searched recursively). Substrates given to
  Boltz as SMILES appear as ligands ``LIG1``, ``LIG2``, ... (CCD ligands keep their CCD code).
* ``--csv``   table (CSV/TSV, or FASTA) with the enzyme sequences that were docked.

What it does
------------
1. Atom-maps every reaction in the YAML (iterative maximum-common-substructure between
   reactants and products, ignoring H) and finds the reacting atoms:
     tier 1  bond formed between two different substrates   (e.g. NAS phenol O + SAM CH3)
     tier 2  H / hydride transfer donor -> acceptor         (e.g. NADPH C4 -> P2C imine C)
     tier 3  most-changed atoms of two substrates (fallback)
   An explicit choice can be given per reaction with ``atoms_of_interest`` (see below).
2. For every CIF: extracts the protein sequence(s) and finds the enzyme in the CSV, identifies
   which ligand is which substrate (Boltz atom names -> heavy-atom graph -> 3D geometry),
   works out which reaction was docked, and measures the distance (minimum over
   symmetry-equivalent atoms, e.g. both O of O2, and over ligand copies).
3. Writes
     distance_atoms_of_interest.csv        one row per CIF model, all details
     distance_atoms_of_interest_wide.csv   one row per enzyme copy, <rid>_distance_A columns
     reaction_atoms_of_interest.csv        which atoms were chosen for each reaction and why
     aoi_depictions/                       2D drawings with the chosen atoms highlighted

Usage
-----
    python boltz_reaction_distances.py --yaml pathway.yaml --cif-dir boltz_results/ \\
        --csv enzymes.csv -o distance_atoms_of_interest.csv

    # check the chemistry first (no CIFs needed)
    python boltz_reaction_distances.py --yaml pathway.yaml --analyze-only -o out/

The CSV columns are auto-detected (sequence, ID, protein family, reaction id); override with
--seq-col / --id-col / --protein-col / --reaction-col / --sep. Each CIF is linked to its enzyme by
the protein sequence it contains; the protein family comes from the CSV, the file name, or
--protein; the reaction from the docked ligands (plus rr/r1/... tokens in the file name).

Optional per-reaction override in the YAML (exactly two selectors)::

    - id: rr
      ...
      atoms_of_interest:
        - {compound: tryptophan, smarts: '[NX3;H2:1][CX4]C(=O)O'}   # atom tagged :1
        - {residue: PLP, atom: C4A}                                # any CIF ligand atom

Requirements: python>=3.9, rdkit, numpy, pandas, pyyaml (all present in a Boltz-2 env).
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import html
import itertools
import json
import logging
import math
import re
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
import yaml
from rdkit import Chem, RDLogger
from rdkit.Chem import rdDepictor, rdFMCS
from rdkit.Chem.Draw import rdMolDraw2D

RDLogger.DisableLog("rdApp.*")
log = logging.getLogger("boltz_rxn_dist")

__version__ = "1.0"

# --------------------------------------------------------------------------------------
# constants
# --------------------------------------------------------------------------------------
AA3_TO_1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
    "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
    "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V", "SEC": "U",
    "PYL": "O", "MSE": "M", "UNK": "X", "ASX": "B", "GLX": "Z", "HID": "H", "HIE": "H",
    "HIP": "H", "CYX": "C",
}
NUCLEOTIDES = {"DA", "DC", "DG", "DT", "DU", "DN", "A", "C", "G", "U", "N"}
WATER_NAMES = {"HOH", "WAT", "DOD", "H2O", "TIP", "TIP3", "SOL"}
PT = Chem.GetPeriodicTable()
TIERS = {0: "override", 1: "bond_formed", 2: "H_transfer", 3: "reaction_center"}
GRAPH_TOLERANCES = (0.30, 0.45, 0.60)   # A added to the sum of covalent radii
AA_RE = re.compile(r"^[ACDEFGHIKLMNPQRSTVWYXBZUO*\-.]+$", re.I)


def _rcov(z: int) -> float:
    try:
        r = PT.GetRcovalent(int(z))
        return r if r > 0 else 0.77
    except Exception:  # pragma: no cover
        return 0.77


def _bo(bond) -> float:
    if bond is None:
        return 0.0
    if bond.GetIsAromatic():
        return 1.5
    return {
        Chem.BondType.SINGLE: 1.0, Chem.BondType.DOUBLE: 2.0,
        Chem.BondType.TRIPLE: 3.0, Chem.BondType.AROMATIC: 1.5,
    }.get(bond.GetBondType(), 1.0)


def _skeleton(mol: Chem.Mol, rings: bool = False) -> Chem.Mol:
    """Element-only graph: same atoms, all bonds single, no H, no charges."""
    rw = Chem.RWMol()
    for a in mol.GetAtoms():
        na = Chem.Atom(a.GetAtomicNum())
        na.SetNoImplicit(True)
        rw.AddAtom(na)
    for b in mol.GetBonds():
        rw.AddBond(b.GetBeginAtomIdx(), b.GetEndAtomIdx(), Chem.BondType.SINGLE)
    m = rw.GetMol()
    m.UpdatePropertyCache(strict=False)
    if rings:
        Chem.FastFindRings(m)
    return m


def _submol(skel: Chem.Mol, keep) -> tuple[Chem.Mol, list[int]]:
    keep = sorted(keep)
    pos = {old: new for new, old in enumerate(keep)}
    rw = Chem.RWMol()
    for old in keep:
        na = Chem.Atom(skel.GetAtomWithIdx(old).GetAtomicNum())
        na.SetNoImplicit(True)
        rw.AddAtom(na)
    for b in skel.GetBonds():
        u, v = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        if u in pos and v in pos:
            rw.AddBond(pos[u], pos[v], Chem.BondType.SINGLE)
    m = rw.GetMol()
    m.UpdatePropertyCache(strict=False)
    Chem.FastFindRings(m)
    return m, keep


def _formula_str(counter: Counter) -> str:
    if not counter:
        return ""
    items = sorted(counter.items(), key=lambda kv: (kv[0] != 6, kv[0] != 1, PT.GetElementSymbol(kv[0])))
    return "".join(f"{PT.GetElementSymbol(z)}{n if n > 1 else ''}" for z, n in items)


# --------------------------------------------------------------------------------------
# chemistry: one molecule of a reaction
# --------------------------------------------------------------------------------------
class Species:
    """One reactant or product fragment of a reaction string (heavy atoms only)."""

    def __init__(self, smiles: str, side: str, index: int, key: Optional[str] = None,
                 name: Optional[str] = None):
        self.smiles = smiles.strip()
        self.side = side            # "R" or "P"
        self.index = index
        self.key = key
        self.name = name or key or f"{'reactant' if side == 'R' else 'product'}{index + 1}"
        self.parse_warning = ""
        raw = Chem.MolFromSmiles(self.smiles)
        if raw is None:
            raw = Chem.MolFromSmiles(self.smiles, sanitize=False)
            if raw is not None:
                raw.UpdatePropertyCache(strict=False)
                Chem.FastFindRings(raw)
                self.parse_warning = "could not sanitize SMILES; H counts may be unreliable"
        if raw is None:
            raw = Chem.MolFromSmarts(self.smiles)
            if raw is None:
                raise ValueError(f"cannot parse fragment {self.smiles!r}")
            raw.UpdatePropertyCache(strict=False)
            self.parse_warning = "parsed as SMARTS; H counts may be unreliable"
        # Boltz-2 atom names: symbol.upper() + (canonical rank on the H-added molecule + 1)
        self.boltz_names: Optional[list[str]] = None
        try:
            mh = Chem.AddHs(raw)
            ranks = list(Chem.CanonicalRankAtoms(mh))
            self.boltz_names = [a.GetSymbol().upper() + str(r + 1)
                                for a, r in zip(mh.GetAtoms(), ranks) if a.GetAtomicNum() > 1]
        except Exception:
            self.boltz_names = None
        # heavy-atom molecule (H become implicit counts)
        try:
            mol = Chem.RemoveAllHs(raw, sanitize=not self.parse_warning)
        except Exception:
            mol = Chem.RemoveHs(raw, sanitize=False)
        if any(a.GetAtomicNum() <= 1 for a in mol.GetAtoms()):
            rw = Chem.RWMol(mol)
            for a in sorted((a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() <= 1), reverse=True):
                rw.RemoveAtom(a)
            mol = rw.GetMol()
            mol.UpdatePropertyCache(strict=False)
        self.mol = mol
        self.n_heavy = mol.GetNumAtoms()
        if self.boltz_names is not None and len(self.boltz_names) != self.n_heavy:
            self.boltz_names = None
        self.elements = [a.GetAtomicNum() for a in mol.GetAtoms()]
        self.h = [a.GetTotalNumHs() for a in mol.GetAtoms()]
        self.q = [a.GetFormalCharge() for a in mol.GetAtoms()]
        self.formula = Counter(self.elements)
        self.skel = _skeleton(mol) if self.n_heavy else None
        self.bonds = [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in mol.GetBonds()]
        try:
            self.canon = Chem.MolToSmiles(mol) if self.n_heavy else self.smiles
        except Exception:
            self.canon = self.smiles
        self._autos = None

    # automorphisms of the heavy-atom graph (used for symmetry-equivalent atoms)
    @property
    def autos(self) -> list[tuple]:
        if self._autos is None:
            if not self.n_heavy:
                self._autos = []
            else:
                ms = self.skel.GetSubstructMatches(self.skel, uniquify=False, useChirality=False,
                                                   maxMatches=20000)
                self._autos = list(ms) or [tuple(range(self.n_heavy))]
        return self._autos

    def atom_label(self, idx: int) -> str:
        """Expected Boltz atom name (e.g. C41); element#index if names are unavailable."""
        if self.boltz_names:
            return self.boltz_names[idx]
        return f"{PT.GetElementSymbol(self.elements[idx])}#{idx}"

    def boltz_name(self, idx: int) -> str:
        return self.boltz_names[idx] if self.boltz_names else ""

    def orbit(self, idx: int) -> tuple:
        """Topologically equivalent atoms (e.g. both O of O2)."""
        return tuple(sorted({a[idx] for a in self.autos})) if self.autos else (idx,)

    def env(self, idx: int) -> str:
        """Radius-1 environment as SMILES, the atom itself tagged with :1 (e.g. [CH3:1][S+])."""
        try:
            m = Chem.Mol(self.mol)
            for a in m.GetAtoms():
                a.SetAtomMapNum(0)
            at = m.GetAtomWithIdx(idx)
            at.SetAtomMapNum(1)
            atoms = [idx] + [n.GetIdx() for n in at.GetNeighbors()]
            bonds = [b.GetIdx() for b in m.GetBonds()
                     if b.GetBeginAtomIdx() in atoms and b.GetEndAtomIdx() in atoms]
            if not bonds:
                return Chem.MolFragmentToSmiles(m, atomsToUse=[idx], isomericSmiles=False)
            return Chem.MolFragmentToSmiles(m, atomsToUse=atoms, bondsToUse=bonds,
                                            isomericSmiles=False)
        except Exception:
            return ""


# --------------------------------------------------------------------------------------
# reaction atom mapping (iterative MCS on heavy-atom graphs)
# --------------------------------------------------------------------------------------
def _elem_overlap(a: Species, keep_a, b: Species, keep_b) -> int:
    ca = Counter(a.elements[i] for i in keep_a)
    cb = Counter(b.elements[i] for i in keep_b)
    return sum(min(n, cb.get(z, 0)) for z, n in ca.items())


def _link_bonus(i, j, pairs, R, P, amap) -> int:
    """Product bonds from newly mapped atoms to atoms already mapped from the same reactant
    that also exist in the reactant (rewards keeping molecules together)."""
    link = 0
    Pm, Rm = P[j].mol, R[i].mol
    for ra, pa in pairs:
        for nb in Pm.GetAtomWithIdx(pa).GetNeighbors():
            src = amap.get((j, nb.GetIdx()))
            if src is not None and src[0] == i and Rm.GetBondBetweenAtoms(ra, src[1]) is not None:
                link += 1
    return link


def _best_mcs_mapping(i, j, R, P, keepR, keepP, amap, timeout, max_combos):
    Ri, Pj = R[i], P[j]
    subR, idxR = _submol(Ri.skel, keepR)
    subP, idxP = _submol(Pj.skel, keepP)
    res = rdFMCS.FindMCS(
        [subR, subP],
        atomCompare=rdFMCS.AtomCompare.CompareElements,
        bondCompare=rdFMCS.BondCompare.CompareAny,
        ringMatchesRingOnly=False, completeRingsOnly=False,
        timeout=max(1, int(timeout)),
    )
    if res.numAtoms == 0:
        return None
    q = Chem.MolFromSmarts(res.smartsString)
    if q is None:
        return None
    mR = subR.GetSubstructMatches(q, uniquify=True, useChirality=False, maxMatches=200)
    mP = subP.GetSubstructMatches(q, uniquify=False, useChirality=False, maxMatches=10000)
    if not mR or not mP:
        return None
    qbonds = [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in q.GetBonds()]
    best = None
    n = 0
    for a in mR:
        for b in mP:
            n += 1
            if n > max_combos:
                break
            pairs = [(idxR[a[k]], idxP[b[k]]) for k in range(len(a))]
            cost = 0.0
            for ra, pa in pairs:
                cost += abs(Ri.h[ra] - Pj.h[pa]) + abs(Ri.q[ra] - Pj.q[pa])
            for k1, k2 in qbonds:
                cost += abs(_bo(Ri.mol.GetBondBetweenAtoms(pairs[k1][0], pairs[k2][0]))
                            - _bo(Pj.mol.GetBondBetweenAtoms(pairs[k1][1], pairs[k2][1])))
            link = _link_bonus(i, j, pairs, R, P, amap)
            score = cost - 2.0 * link
            if best is None or score < best[0] - 1e-9:
                best = (score, pairs)
        if n > max_combos:
            break
    return res.numAtoms, res.numBonds, best[0], best[1]


def map_reaction(R: list[Species], P: list[Species], timeout: float = 10.0,
                 max_combos: int = 20000) -> dict:
    """Greedy iterative MCS. Returns {(product_idx, atom): (reactant_idx, atom)}."""
    unR = {i: set(range(s.n_heavy)) for i, s in enumerate(R) if s.n_heavy}
    unP = {j: set(range(s.n_heavy)) for j, s in enumerate(P) if s.n_heavy}
    amap: dict = {}
    cache: dict = {}
    while True:
        best = None
        for i in sorted(unR):
            if not unR[i]:
                continue
            for j in sorted(unP):
                if not unP[j]:
                    continue
                ub = _elem_overlap(R[i], unR[i], P[j], unP[j])
                if ub == 0 or (best is not None and ub < best[0][0]):
                    continue
                ck = (i, j, frozenset(unR[i]), frozenset(unP[j]))
                if ck not in cache:
                    cache[ck] = _best_mcs_mapping(i, j, R, P, unR[i], unP[j], amap, timeout, max_combos)
                res = cache[ck]
                if res is None:
                    continue
                size, nbonds, score, pairs = res
                key = (size, nbonds, -score, -i, -j)
                if best is None or key > best[0]:
                    best = (key, i, j, pairs)
        if best is None:
            break
        _, i, j, pairs = best
        for ra, pa in pairs:
            amap[(j, pa)] = (i, ra)
            unR[i].discard(ra)
            unP[j].discard(pa)
    return amap


@dataclass
class Events:
    formed_inter: list = field(default_factory=list)   # [((i,a),(k,b), order)]
    formed_intra: list = field(default_factory=list)
    broken: list = field(default_factory=list)
    bo_change: list = field(default_factory=list)
    unknown_partner: dict = field(default_factory=dict)
    dH: dict = field(default_factory=dict)
    dq: dict = field(default_factory=dict)
    lost: set = field(default_factory=set)
    unmapped_products: list = field(default_factory=list)


def reaction_events(R: list[Species], P: list[Species], amap: dict) -> Events:
    ev = Events()
    rmap = {v: k for k, v in amap.items()}
    unknown = defaultdict(int)
    for j, Pj in enumerate(P):
        if not Pj.n_heavy:
            continue
        for b in Pj.mol.GetBonds():
            u, v = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
            ru, rv = amap.get((j, u)), amap.get((j, v))
            if ru is not None and rv is not None:
                if ru[0] != rv[0]:
                    ev.formed_inter.append((ru, rv, _bo(b)))
                else:
                    rb = R[ru[0]].mol.GetBondBetweenAtoms(ru[1], rv[1])
                    if rb is None:
                        ev.formed_intra.append((ru, rv, _bo(b)))
                    elif abs(_bo(rb) - _bo(b)) > 1e-6:
                        ev.bo_change.append((ru, rv, _bo(rb), _bo(b)))
            elif ru is not None or rv is not None:
                unknown[ru if ru is not None else rv] += 1
    ev.unknown_partner = dict(unknown)
    for i, Ri in enumerate(R):
        if not Ri.n_heavy:
            continue
        for u, v in Ri.bonds:
            pu, pv = rmap.get((i, u)), rmap.get((i, v))
            if (pu is None or pv is None or pu[0] != pv[0]
                    or P[pu[0]].mol.GetBondBetweenAtoms(pu[1], pv[1]) is None):
                ev.broken.append(((i, u), (i, v), _bo(Ri.mol.GetBondBetweenAtoms(u, v))))
        for a in range(Ri.n_heavy):
            p = rmap.get((i, a))
            if p is None:
                ev.lost.add((i, a))
                continue
            ev.dH[(i, a)] = P[p[0]].h[p[1]] - Ri.h[a]
            ev.dq[(i, a)] = P[p[0]].q[p[1]] - Ri.q[a]
    ev.unmapped_products = [(j, a) for j, Pj in enumerate(P) for a in range(Pj.n_heavy)
                            if (j, a) not in amap]
    return ev


@dataclass(frozen=True)
class Sel:
    """Atom selector: species (reactant index) + candidate atom indices, or a CIF residue/atom."""
    species: Optional[int] = None
    atoms: tuple = ()
    resname: Optional[str] = None
    atomname: Optional[str] = None


@dataclass
class Pair:
    tier: int
    sub: float
    a: Sel
    b: Sel

    @property
    def kind(self) -> str:
        return TIERS[self.tier]

    @property
    def rank(self) -> tuple:
        return (self.tier, self.sub)


def _change_score(ev: Events, x) -> int:
    s = 0
    nb = sum(1 for e in ev.formed_inter + ev.formed_intra + ev.broken if x in (e[0], e[1]))
    s += 4 * (nb + ev.unknown_partner.get(x, 0))
    s += 3 * abs(ev.dH.get(x, 0)) + 2 * abs(ev.dq.get(x, 0))
    s += sum(1 for e in ev.bo_change if x in (e[0], e[1]))
    s += 1 if x in ev.lost else 0
    return s


def candidate_pairs(R: list[Species], ev: Events) -> list[Pair]:
    """All candidate atom pairs, best first. Atoms are expanded to their symmetry orbit."""
    raw: list[tuple] = []
    # tier 1: bond formed between atoms of two different reactants
    #         (prefer bonds involving a heteroatom, e.g. amine N -> C=O before C-C ring closure)
    for ru, rv, _order in ev.formed_inter:
        a, b = sorted([ru, rv])
        het = sum(R[x[0]].elements[x[1]] != 6 for x in (a, b))
        raw.append((1, -het, a, b))
    # tier 2: H / hydride transfer (atom losing H in one reactant, gaining H in another);
    #         prefer C->C (hydride) over transfers involving N/O (proton shuffles)
    donors = sorted(x for x, d in ev.dH.items() if d < 0)
    acceptors = sorted(x for x, d in ev.dH.items() if d > 0)
    for d in donors:
        for a in acceptors:
            if d[0] != a[0]:
                ncarb = sum(R[x[0]].elements[x[1]] == 6 for x in (d, a))
                raw.append((2, -ncarb, d, a))
    # tier 3: most-changed atom(s) of each pair of reactants (fallback)
    mols = [i for i, s in enumerate(R) if s.n_heavy]
    top = {}
    for i in mols:
        sc = {a: _change_score(ev, (i, a)) for a in range(R[i].n_heavy)}
        mx = max(sc.values()) if sc else 0
        top[i] = [(a, s) for a, s in sc.items() if s == mx and s > 0]
    for i1, i2 in itertools.combinations(mols, 2):
        for a1, s1 in top[i1]:
            for a2, s2 in top[i2]:
                raw.append((3, -(s1 + s2), (i1, a1), (i2, a2)))
    pairs, seen = [], set()
    for tier, sub, x, y in raw:
        sa = Sel(x[0], R[x[0]].orbit(x[1]))
        sb = Sel(y[0], R[y[0]].orbit(y[1]))
        if tier != 2:  # unordered
            sa, sb = sorted([sa, sb], key=lambda s: (s.species, s.atoms))
        key = (tier, sa, sb)
        if key in seen:
            continue
        seen.add(key)
        pairs.append(Pair(tier, sub, sa, sb))
    pairs.sort(key=lambda p: (p.rank, p.a.species, p.a.atoms, p.b.species, p.b.atoms))
    return pairs


# --------------------------------------------------------------------------------------
# reaction analysis for the whole YAML
# --------------------------------------------------------------------------------------
def split_reaction(rxn: str) -> tuple[list[str], list[str]]:
    rxn = rxn.strip()
    if ">>" in rxn:
        lhs, rhs = rxn.split(">>", 1)
    else:
        parts = rxn.split(">")
        if len(parts) < 2:
            raise ValueError("no '>>' in reaction string")
        lhs, rhs = parts[0], parts[-1]
    frag = lambda side: [f for f in side.split(".") if f.strip()]  # noqa: E731
    return frag(lhs), frag(rhs)


@dataclass
class ReactionInfo:
    protein: str
    rid: str
    order: int
    rxn: dict
    reactants: list
    products: list
    amap: dict
    events: Events
    pairs: list
    notes: list
    source_field: str

    @property
    def uid(self) -> str:
        return f"{self.protein}/{self.rid}"

    @property
    def heavy_reactants(self) -> list[int]:
        return [i for i, s in enumerate(self.reactants) if s.n_heavy]

    @property
    def chem_key(self) -> tuple:
        return tuple(sorted(self.reactants[i].canon for i in self.heavy_reactants))

    def best_pairs(self) -> list[Pair]:
        if not self.pairs:
            return []
        r = min(p.rank for p in self.pairs)
        return [p for p in self.pairs if p.rank == r]

    def sel_label(self, sel: Sel, with_env: bool = True) -> str:
        """e.g. 'SAM:C41 [S+][CH3:1]' (atom = expected Boltz atom name; | = equivalent atoms)."""
        if sel.species is None:
            return f"{sel.resname}:{sel.atomname}"
        sp = self.reactants[sel.species]
        names = "|".join(sp.atom_label(a) for a in sel.atoms)
        env = f" {sp.env(sel.atoms[0])}" if with_env and sel.atoms else ""
        return f"{sp.name}:{names}{env}"


def _compound_name(compounds: dict, key: Optional[str]) -> Optional[str]:
    if key is None:
        return None
    c = compounds.get(key) or {}
    return str(c.get("name") or key)


def parse_override(rxn: dict, info: ReactionInfo) -> Optional[Pair]:
    ov = rxn.get("atoms_of_interest")
    if not ov:
        return None
    if not isinstance(ov, list) or len(ov) != 2:
        info.notes.append("atoms_of_interest ignored: needs exactly two selectors")
        return None
    sels = []
    for s in ov:
        if not isinstance(s, dict):
            info.notes.append(f"atoms_of_interest selector not understood: {s!r}")
            return None
        if "residue" in s:
            sels.append(Sel(resname=str(s["residue"]).upper(), atomname=str(s.get("atom", "")).upper()))
            continue
        comp = s.get("compound")
        idx = None
        for i, sp in enumerate(info.reactants):
            if comp is not None and (str(comp) == str(sp.key) or str(comp).lower() == str(sp.name).lower()
                                     or str(comp) == str(i)):
                idx = i
                break
        if idx is None:
            info.notes.append(f"atoms_of_interest: compound {comp!r} is not a reactant of this reaction")
            return None
        sp = info.reactants[idx]
        patt = Chem.MolFromSmarts(str(s.get("smarts", "")))
        if patt is None:
            info.notes.append(f"atoms_of_interest: bad SMARTS {s.get('smarts')!r}")
            return None
        tag = next((a.GetIdx() for a in patt.GetAtoms() if a.GetAtomMapNum() == 1), 0)
        hits = sp.mol.GetSubstructMatches(patt, uniquify=False)
        atoms = tuple(sorted({h[tag] for h in hits}))
        if not atoms:
            info.notes.append(f"atoms_of_interest: SMARTS {s.get('smarts')!r} not found in {sp.name}")
            return None
        sels.append(Sel(idx, atoms))
    return Pair(0, 0.0, sels[0], sels[1])


def analyze_reaction(pkey: str, rxn: dict, order: int, compounds: dict,
                     timeout: float) -> ReactionInfo:
    rid = str(rxn.get("id", f"rxn{order}"))
    notes: list[str] = []
    field_used = "smarts" if rxn.get("smarts") else "smiles"
    rstr = rxn.get("smarts") or rxn.get("smiles")
    if not rstr:
        raise ValueError("reaction has neither 'smarts' nor 'smiles'")
    lhs, rhs = split_reaction(str(rstr))
    subs = [str(x) for x in (rxn.get("substrates") or [])]
    prods = [str(x) for x in (rxn.get("products") or [])]
    if len(subs) != len(lhs):
        notes.append(f"{len(lhs)} reactant fragments in {field_used} but {len(subs)} substrates listed; "
                     f"fragments named by position")
    if len(prods) != len(rhs):
        notes.append(f"{len(rhs)} product fragments in {field_used} but {len(prods)} products listed")
    R = [Species(s, "R", i, subs[i] if len(subs) == len(lhs) else None,
                 _compound_name(compounds, subs[i]) if len(subs) == len(lhs) else None)
         for i, s in enumerate(lhs)]
    P = [Species(s, "P", j, prods[j] if len(prods) == len(rhs) else None,
                 _compound_name(compounds, prods[j]) if len(prods) == len(rhs) else None)
         for j, s in enumerate(rhs)]
    for sp in R + P:
        if sp.parse_warning:
            notes.append(f"{sp.name}: {sp.parse_warning}")
    fr = sum((s.formula for s in R), Counter())
    fp = sum((s.formula for s in P), Counter())
    if fr != fp:
        notes.append(f"heavy atoms not balanced: {_formula_str(fr)} -> {_formula_str(fp)}")
    amap = map_reaction(R, P, timeout=timeout)
    ev = reaction_events(R, P, amap)
    pairs = candidate_pairs(R, ev)
    info = ReactionInfo(pkey, rid, order, rxn, R, P, amap, ev, pairs, notes, field_used)
    ov = parse_override(rxn, info)
    if ov is not None:
        info.pairs.insert(0, ov)
    if len(info.heavy_reactants) < 2 and ov is None:
        notes.append("only one substrate with heavy atoms: no intermolecular atom pair "
                     "(add 'atoms_of_interest' to measure to a cofactor)")
    return info


def analyze_config(cfg: dict, timeout: float, only_proteins=None) -> dict:
    compounds = cfg.get("compounds") or {}
    out: dict = {}
    order = 0
    for pkey, prot in (cfg.get("proteins") or {}).items():
        pkey = str(pkey)
        if only_proteins and pkey not in only_proteins:
            continue
        for rxn in prot.get("reactions") or []:
            try:
                info = analyze_reaction(pkey, rxn, order, compounds, timeout)
            except Exception as e:  # keep going with the other reactions
                log.error("reaction %s/%s could not be analysed: %s", pkey, rxn.get("id"), e)
                order += 1
                continue
            out[info.uid] = info
            order += 1
    _consistency_checks(out)
    return out


def _consistency_checks(infos: dict) -> None:
    """Same compound key -> same structure across reactions, and vice versa."""
    key2canon = defaultdict(set)
    canon2keys = defaultdict(set)
    where = defaultdict(list)
    for info in infos.values():
        for sp in info.reactants + info.products:
            if sp.key is None or not sp.n_heavy:
                continue
            key2canon[sp.key].add(sp.canon)
            canon2keys[sp.canon].add(sp.key)
            where[sp.key].append(info.uid)
    for info in infos.values():
        for sp in info.reactants + info.products:
            if sp.key is None or not sp.n_heavy:
                continue
            others = sorted(canon2keys[sp.canon] - {sp.key})
            if others:
                used = sorted({u for o in others for u in where[o]})
                info.notes.append(f"'{sp.key}' has the same structure as '{'/'.join(others)}' "
                                  f"(used in {', '.join(used)}) - check the SMARTS")
            if len(key2canon[sp.key]) > 1:
                info.notes.append(f"'{sp.key}' is drawn with {len(key2canon[sp.key])} different "
                                  f"structures across reactions")
    for info in infos.values():  # de-duplicate
        info.notes = list(dict.fromkeys(info.notes))


# --------------------------------------------------------------------------------------
# minimal mmCIF reader (only _atom_site is needed)
# --------------------------------------------------------------------------------------
def _cif_split(line: str) -> list[str]:
    toks, i, n = [], 0, len(line)
    while i < n:
        c = line[i]
        if c.isspace():
            i += 1
            continue
        if c == "#":
            break
        if c in "'\"":
            j = i + 1
            while j < n and not (line[j] == c and (j + 1 == n or line[j + 1].isspace())):
                j += 1
            toks.append(line[i + 1:j])
            i = j + 1
        else:
            j = i
            while j < n and not line[j].isspace():
                j += 1
            toks.append(line[i:j])
            i = j
    return toks


def read_atom_site(path: Path) -> list[dict]:
    text = path.read_text(errors="replace").splitlines()
    tags: list[str] = []
    values: list[str] = []
    i, n = 0, len(text)
    while i < n:
        s = text[i].strip()
        if s == "loop_":
            j = i + 1
            loop_tags = []
            while j < n and text[j].strip().startswith("_"):
                loop_tags.append(text[j].strip().split()[0])
                j += 1
            if loop_tags and loop_tags[0].lower().startswith("_atom_site."):
                tags = [t.split(".", 1)[1] for t in loop_tags]
                k = j
                while k < n:
                    line = text[k]
                    st = line.strip()
                    if st.startswith(("loop_", "_", "data_", "save_")):
                        break
                    if line.startswith(";"):  # multi-line text field
                        buf = [line[1:]]
                        k += 1
                        while k < n and not text[k].startswith(";"):
                            buf.append(text[k])
                            k += 1
                        values.append("\n".join(buf).strip())
                        k += 1
                        continue
                    values.extend(_cif_split(line))
                    k += 1
                break
            i = j
            continue
        i += 1
    if not tags:
        return []
    nt = len(tags)
    if len(values) % nt:
        log.warning("%s: _atom_site has %d values for %d columns", path.name, len(values), nt)
    rows = []
    for r in range(len(values) // nt):
        rows.append(dict(zip(tags, values[r * nt:(r + 1) * nt])))
    return rows


def _elem_from(row: dict) -> int:
    sym = (row.get("type_symbol") or "").strip()
    if sym and sym not in ("?", ".", "X"):
        s = sym[0].upper() + sym[1:].lower()
        try:
            return PT.GetAtomicNumber(s)
        except Exception:
            pass
    name = re.sub(r"[^A-Za-z]", "", row.get("label_atom_id") or row.get("auth_atom_id") or "")
    for cand in (name[:2], name[:1]):
        if cand:
            try:
                return PT.GetAtomicNumber(cand[0].upper() + cand[1:].lower())
            except Exception:
                continue
    return 0


class LigandResidue:
    def __init__(self, chain: str, resname: str, resseq: str, names, elements, coords):
        keep = [k for k, z in enumerate(elements) if z > 1]
        self.chain = chain
        self.resname = resname
        self.resseq = resseq
        self.names = [names[k] for k in keep]
        self.elements = [elements[k] for k in keep]
        self.coords = np.asarray([coords[k] for k in keep], dtype=float).reshape(-1, 3)
        self.formula = Counter(self.elements)
        self.name_index = {nm: k for k, nm in enumerate(self.names)}
        self._D = None
        self._skels: dict = {}

    @property
    def label(self) -> str:
        return f"{self.chain}:{self.resname}"

    @property
    def D(self) -> np.ndarray:
        if self._D is None:
            d = self.coords[:, None, :] - self.coords[None, :, :]
            self._D = np.sqrt((d * d).sum(-1))
        return self._D

    def ideal(self, i: int, j: int) -> float:
        return _rcov(self.elements[i]) + _rcov(self.elements[j])

    def skeleton(self, tol: float) -> Chem.Mol:
        if tol not in self._skels:
            rw = Chem.RWMol()
            for z in self.elements:
                a = Chem.Atom(int(z))
                a.SetNoImplicit(True)
                rw.AddAtom(a)
            D, n = self.D, len(self.elements)
            for i in range(n):
                for j in range(i + 1, n):
                    if 0.4 < D[i, j] <= self.ideal(i, j) + tol:
                        rw.AddBond(i, j, Chem.BondType.SINGLE)
            m = rw.GetMol()
            m.UpdatePropertyCache(strict=False)
            self._skels[tol] = m
        return self._skels[tol]


@dataclass
class CifModel:
    model: str
    chains: dict            # chain -> one-letter sequence (protein chains)
    ligands: list           # [LigandResidue]


def parse_cif(path: Path) -> list[CifModel]:
    rows = read_atom_site(path)
    if not rows:
        raise ValueError("no _atom_site records")
    by_model = defaultdict(list)
    for r in rows:
        by_model[r.get("pdbx_PDB_model_num", "1")].append(r)
    models = []
    for mnum, mrows in by_model.items():
        residues = {}   # key -> dict
        order = []
        poly_asym = set()   # label_asym_ids (or chains) that contain ATOM records
        for r in mrows:
            alt = r.get("label_alt_id", ".")
            if alt not in (".", "?", "A", "1", ""):
                continue
            chain = r.get("auth_asym_id") or r.get("label_asym_id") or "?"
            asym = r.get("label_asym_id") or chain
            comp = (r.get("label_comp_id") or r.get("auth_comp_id") or "UNK").upper()
            seq = r.get("auth_seq_id") or r.get("label_seq_id") or "."
            ins = r.get("pdbx_PDB_ins_code", "?")
            group = (r.get("group_PDB") or "ATOM").upper()
            if group == "ATOM":
                poly_asym.add(asym)
            key = (chain, asym, seq, ins, comp)
            if key not in residues:
                residues[key] = {"group": group, "names": [], "elements": [], "coords": []}
                order.append(key)
            try:
                xyz = (float(r["Cartn_x"]), float(r["Cartn_y"]), float(r["Cartn_z"]))
            except (KeyError, ValueError):
                continue
            residues[key]["names"].append((r.get("label_atom_id") or r.get("auth_atom_id") or "").strip('"').upper())
            residues[key]["elements"].append(_elem_from(r))
            residues[key]["coords"].append(xyz)
        chains = defaultdict(list)
        ligands = []
        for key in order:
            chain, asym, seq, ins, comp = key
            res = residues[key]
            if comp in WATER_NAMES:
                continue
            # HETATM counts as polymer only for modified residues inside a polymer entity
            # (same label_asym_id as ATOM records); a docked TRP/SAM/... in its own asym is a ligand
            polymer = res["group"] == "ATOM" or (asym in poly_asym and (comp in AA3_TO_1 or comp in NUCLEOTIDES))
            if polymer:
                chains[chain].append(comp)
            elif any(z > 1 for z in res["elements"]):   # ignore lone H+/H atoms
                ligands.append(LigandResidue(chain, comp, seq, res["names"], res["elements"], res["coords"]))
        prot = {}
        for ch, comps in chains.items():
            n_aa = sum(c in AA3_TO_1 for c in comps)
            if n_aa >= max(1, 0.5 * len(comps)):
                prot[ch] = "".join(AA3_TO_1.get(c, "X") for c in comps)
        models.append(CifModel(str(mnum), prot, ligands))
    return models


# --------------------------------------------------------------------------------------
# ligand residue  <->  reaction species
# --------------------------------------------------------------------------------------
@dataclass
class ResidueMatch:
    residue: LigandResidue
    matches: list           # list of tuples: species atom -> residue atom (symmetry equivalents)
    method: str
    name_overlap: float
    geom_dev: Optional[float]
    quality: int            # 3 exact Boltz names, 2 partial names + graph, 1 graph only

    @property
    def label(self) -> str:
        return self.residue.label


def _geom_dev(sp: Species, res: LigandResidue, m: tuple) -> Optional[float]:
    """Sum over atoms with >=2 heavy neighbours of |mean bond angle - ideal(hybridisation)|."""
    tot, n = 0.0, 0
    for atom in sp.mol.GetAtoms():
        nb = [x.GetIdx() for x in atom.GetNeighbors()]
        if len(nb) < 2:
            continue
        hyb = atom.GetHybridization()
        if atom.GetIsAromatic() or hyb == Chem.HybridizationType.SP2:
            ideal = 120.0
        elif hyb == Chem.HybridizationType.SP3:
            ideal = 109.47
        elif hyb == Chem.HybridizationType.SP:
            ideal = 180.0
        else:
            continue
        c = res.coords[m[atom.GetIdx()]]
        vs = [res.coords[m[k]] - c for k in nb]
        angs = []
        for v1, v2 in itertools.combinations(vs, 2):
            den = np.linalg.norm(v1) * np.linalg.norm(v2)
            if den < 1e-6:
                continue
            angs.append(math.degrees(math.acos(max(-1.0, min(1.0, float(np.dot(v1, v2) / den))))))
        if angs:
            tot += abs(sum(angs) / len(angs) - ideal)
            n += 1
    return tot if n else None


def _name_overlap(sp: Species, res: LigandResidue) -> float:
    if not sp.boltz_names or not sp.n_heavy:
        return 0.0
    hit = 0
    for k, nm in enumerate(sp.boltz_names):
        idx = res.name_index.get(nm)
        if idx is not None and res.elements[idx] == sp.elements[k]:
            hit += 1
    return hit / sp.n_heavy


def match_species(sp: Species, res: LigandResidue) -> Optional[ResidueMatch]:
    if not sp.n_heavy or sp.formula != res.formula:
        return None
    overlap = _name_overlap(sp, res)
    # 1) exact Boltz atom names (validated against the geometry)
    if overlap == 1.0:
        m0 = tuple(res.name_index[nm] for nm in sp.boltz_names)
        ok = sum(1 for u, v in sp.bonds if res.D[m0[u], m0[v]] <= res.ideal(m0[u], m0[v]) + 0.8)
        if not sp.bonds or ok >= 0.9 * len(sp.bonds):
            eq = sorted({tuple(m0[a[k]] for k in range(sp.n_heavy)) for a in sp.autos})
            return ResidueMatch(res, eq, "atom_names", overlap, _geom_dev(sp, res, m0), 3)
    # 2) heavy-atom graph from coordinates (tightest tolerance that works)
    for tol in GRAPH_TOLERANCES:
        tgt = res.skeleton(tol)
        ms = tgt.GetSubstructMatches(sp.skel, uniquify=False, useChirality=False, maxMatches=20000)
        if not ms:
            continue
        groups = defaultdict(list)
        for m in ms:
            edges = frozenset(tuple(sorted((m[u], m[v]))) for u, v in sp.bonds)
            groups[edges].append(tuple(m))
        best_edges = min(groups, key=lambda e: sum(abs(res.D[i, j] - res.ideal(i, j)) for i, j in e))
        eq = sorted(set(groups[best_edges]))
        return ResidueMatch(res, eq, f"graph(tol={tol:.2f})", overlap, _geom_dev(sp, res, eq[0]),
                            2 if overlap >= 0.6 else 1)
    # 3) partial names as a last resort (e.g. badly distorted ligand)
    if overlap >= 0.9:
        m0 = [res.name_index.get(nm) for nm in sp.boltz_names]
        if all(x is not None for x in m0) and len(set(m0)) == sp.n_heavy:
            m0 = tuple(m0)
            eq = sorted({tuple(m0[a[k]] for k in range(sp.n_heavy)) for a in sp.autos})
            return ResidueMatch(res, eq, "atom_names(unvalidated)", overlap, _geom_dev(sp, res, m0), 1)
    return None


# --------------------------------------------------------------------------------------
# enzyme table (CSV / TSV / FASTA)
# --------------------------------------------------------------------------------------
def norm_seq(s: str) -> str:
    return re.sub(r"[^A-Za-z]", "", str(s)).upper()


def _read_fasta(path: Path) -> pd.DataFrame:
    recs, cur, buf = [], None, []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if cur is not None:
                recs.append((cur, "".join(buf)))
            cur, buf = line[1:].strip(), []
        elif line.strip():
            buf.append(line.strip())
    if cur is not None:
        recs.append((cur, "".join(buf)))
    return pd.DataFrame({"id": [h.split()[0] for h, _ in recs],
                         "description": [h for h, _ in recs],
                         "sequence": [s for _, s in recs]})


class EnzymeTable:
    SEQ_NAMES = ["sequence", "seq", "protein_sequence", "aa_sequence", "aa_seq", "prot_seq",
                 "protein_seq", "translation", "peptide", "aa", "seqs", "sequences"]
    ID_NAMES = ["id", "seq_id", "sequence_id", "protein_id", "enzyme_id", "gene_id", "name",
                "header", "accession", "transcript_id", "transcript", "gene", "locus", "entry",
                "query", "target", "identifier", "label"]

    def __init__(self, path: Path, cfg: dict, seq_col=None, id_col=None, protein_col=None,
                 reaction_col=None, sep=None):
        if path.suffix.lower() in (".fa", ".fasta", ".faa", ".fas", ".fna"):
            df = _read_fasta(path)
        else:
            sep = sep or self._guess_sep(path)
            df = pd.read_csv(path, sep=sep, dtype=str, keep_default_na=False, encoding="utf-8-sig")
        df.columns = [str(c).strip() for c in df.columns]
        self.df = df
        self.aliases = protein_aliases(cfg)
        self.rids = {str(r.get("id")).lower() for p in (cfg.get("proteins") or {}).values()
                     for r in (p.get("reactions") or [])}
        cols = list(df.columns)
        lower = {c.lower(): c for c in cols}
        # sequence column
        if seq_col:
            if seq_col not in cols:
                raise ValueError(f"--seq-col {seq_col!r} not in CSV columns {cols}")
            self.seq_col = seq_col
        else:
            self.seq_col = None
            scores = {c: self._seq_score(df[c]) for c in cols}
            for nm in self.SEQ_NAMES:
                if nm in lower and scores[lower[nm]] >= 0.5:
                    self.seq_col = lower[nm]
                    break
            if self.seq_col is None:
                best = max(scores, key=scores.get) if scores else None
                if best is not None and scores[best] >= 0.8:
                    self.seq_col = best
            if self.seq_col is None:
                raise ValueError(f"could not find a sequence column in {path.name} "
                                 f"(columns: {cols}); use --seq-col")
        rest = [c for c in cols if c != self.seq_col]
        # protein-family column
        self.protein_col = protein_col if protein_col in cols else None
        if protein_col and not self.protein_col:
            raise ValueError(f"--protein-col {protein_col!r} not in CSV columns {cols}")
        if self.protein_col is None:
            for c in rest:
                vals = [v.strip().lower() for v in df[c] if str(v).strip()]
                if vals and sum(v in self.aliases for v in vals) / len(vals) >= 0.6:
                    self.protein_col = c
                    break
        # reaction column
        self.reaction_col = reaction_col if reaction_col in cols else None
        if self.reaction_col is None:
            for c in rest:
                vals = [v.strip().lower() for v in df[c] if str(v).strip()]
                if vals and sum(v in self.rids for v in vals) / len(vals) >= 0.8:
                    self.reaction_col = c
                    break
        # id column: prefer a conventionally named column, but only if its values are unique
        if id_col:
            if id_col not in cols:
                raise ValueError(f"--id-col {id_col!r} not in CSV columns {cols}")
            self.id_col = id_col
        else:
            def unique(c):
                vals = [v for v in df[c] if str(v).strip()]
                return bool(vals) and len(set(vals)) >= 0.95 * len(df)
            self.id_col = None
            for nm in self.ID_NAMES:
                if nm in lower and lower[nm] != self.seq_col and unique(lower[nm]):
                    self.id_col = lower[nm]
                    break
            if self.id_col is None:
                for c in rest:
                    if c not in (self.protein_col, self.reaction_col) and unique(c):
                        self.id_col = c
                        break
        self.by_seq = defaultdict(list)
        for k, s in enumerate(df[self.seq_col]):
            ns = norm_seq(s)
            if ns:
                self.by_seq[ns].append(k)
        self.by_id = defaultdict(list)
        if self.id_col:
            for k, v in enumerate(df[self.id_col]):
                if str(v).strip():
                    self.by_id[str(v).strip()].append(k)
        log.info("enzyme table %s: %d rows | sequence column '%s' | id column '%s' | "
                 "protein column '%s' | reaction column '%s'", path.name, len(df), self.seq_col,
                 self.id_col, self.protein_col, self.reaction_col)

    @staticmethod
    def _guess_sep(path: Path) -> str:
        """Delimiter giving the most columns with a consistent column count over the first lines."""
        lines = [ln for ln in path.read_text(errors="replace").splitlines()[:30] if ln.strip()]
        best, best_key = ("\t" if path.suffix.lower() in (".tsv", ".tab") else ","), (0, 0.0)
        for d in ("\t", ";", ",", "|"):
            try:
                counts = [len(r) for r in csv.reader(lines, delimiter=d)]
            except csv.Error:
                continue
            if not counts or counts[0] < 2:
                continue
            consistent = sum(c == counts[0] for c in counts) / len(counts)
            key = (1 if consistent >= 0.9 else 0, counts[0] if consistent >= 0.9 else consistent)
            if key > best_key:
                best, best_key = d, key
        return best

    @staticmethod
    def _seq_score(col) -> float:
        vals = [str(v).strip().replace(" ", "") for v in col if str(v).strip()]
        if not vals:
            return 0.0
        nuc = set("ACGTUN-*.")
        ok = sum(1 for v in vals if len(v) >= 20 and AA_RE.match(v) and not set(v.upper()) <= nuc)
        return ok / len(vals)

    def lookup(self, seqs: list[str], stem: str, enzyme_hint: Optional[str] = None):
        """Return (row indices, how)."""
        for s in seqs:
            if s in self.by_seq:
                return list(self.by_seq[s]), "sequence_exact"
        best, best_len = None, 0
        for s in seqs:
            if len(s) < 20:
                continue
            for k, rows in self.by_seq.items():
                if (s in k or k in s) and min(len(s), len(k)) >= 0.8 * max(len(s), len(k)):
                    if min(len(s), len(k)) > best_len:
                        best, best_len = rows, min(len(s), len(k))
        if best:
            return list(best), "sequence_substring"
        if enzyme_hint and enzyme_hint in self.by_id:
            return list(self.by_id[enzyme_hint]), "id_from_regex"
        low = stem.lower()
        for ident in sorted(self.by_id, key=len, reverse=True):
            if len(ident) >= 3 and re.search(r"(?<![A-Za-z0-9])" + re.escape(ident.lower()) + r"(?![A-Za-z0-9])", low):
                return list(self.by_id[ident]), "id_in_filename"
        return [], "none"

    def row_family(self, k: int) -> Optional[str]:
        row = self.df.iloc[k]
        if self.protein_col:
            return self.aliases.get(str(row[self.protein_col]).strip().lower())
        fams = {self.aliases[str(v).strip().lower()] for c, v in row.items()
                if c != self.seq_col and str(v).strip().lower() in self.aliases}
        return fams.pop() if len(fams) == 1 else None

    def row_reaction(self, k: int) -> Optional[str]:
        if not self.reaction_col:
            return None
        v = str(self.df.iloc[k][self.reaction_col]).strip()
        return v or None

    def row_id(self, k: int) -> str:
        if self.id_col:
            return str(self.df.iloc[k][self.id_col])
        return f"row{k + 1}"

    def row_dict(self, k: int) -> dict:
        row = self.df.iloc[k]
        return {f"csv_{c}": row[c] for c in self.df.columns if c != self.seq_col}


def protein_aliases(cfg: dict) -> dict:
    """lower-case alias (key, abbreviation, gene name, name, UniProt id) -> protein key.
    The protein keys themselves always win; aliases shared by two proteins are dropped."""
    prots = cfg.get("proteins") or {}
    seen = defaultdict(set)
    for pk, p in prots.items():
        for v in [p.get("abbreviation"), p.get("gene_name"), p.get("name")]:
            if v:
                seen[str(v).strip().lower()].add(str(pk))
        for org in (p.get("organisms") or {}).values():
            if isinstance(org, dict) and org.get("uniprot_id"):
                seen[str(org["uniprot_id"]).strip().lower()].add(str(pk))
    al = {a: next(iter(ps)) for a, ps in seen.items() if len(ps) == 1}
    al.update({str(pk).strip().lower(): str(pk) for pk in prots})
    return al


# --------------------------------------------------------------------------------------
# per-CIF processing
# --------------------------------------------------------------------------------------
def _tokens(stem: str) -> list[str]:
    return [t for t in re.split(r"[^A-Za-z0-9+]+", stem.lower()) if t]


def _job_and_model(stem: str) -> tuple[str, Optional[int]]:
    m = re.match(r"^(?P<job>.+?)_model_(?P<model>\d+)$", stem)
    if m:
        return m.group("job"), int(m.group("model"))
    return stem, None


def _read_json(p: Path) -> dict:
    try:
        d = json.loads(p.read_text())
        return {f"boltz_{k}": v for k, v in d.items() if isinstance(v, (int, float, str, bool))}
    except Exception:
        return {}


class SpeciesMatcher:
    """Caches species <-> ligand residue matches for one CIF model."""

    def __init__(self, ligands: list[LigandResidue]):
        self.ligands = ligands
        self.cache: dict = {}

    def matches(self, sp: Species) -> list[ResidueMatch]:
        key = (sp.smiles, sp.canon)
        if key not in self.cache:
            out = []
            for lig in self.ligands:
                rm = match_species(sp, lig)
                if rm is not None:
                    out.append(rm)
            self.cache[key] = out
        return self.cache[key]


def score_reaction(info: ReactionInfo, sm: SpeciesMatcher, rid_hints: set) -> tuple:
    n_exact = n_any = 0
    geo = 0.0
    used = set()
    for i in info.heavy_reactants:
        ms = sm.matches(info.reactants[i])
        if not ms:
            continue
        best = max(ms, key=lambda m: (m.quality, m.name_overlap))
        n_any += 1
        n_exact += best.quality == 3
        geo += best.geom_dev or 0.0
        used.update(id(m.residue) for m in ms)
    hint = 1 if info.rid.lower() in rid_hints else 0
    extra = len(sm.ligands) - len(used)
    return (n_exact, hint, n_any, -extra, -round(geo, 1), -info.order)


def assign_species(info: ReactionInfo, sm: SpeciesMatcher) -> dict:
    """reactant index -> [ResidueMatch]; each ligand residue goes to its best-matching species."""
    best_for_res: dict = {}
    for i in info.heavy_reactants:
        for rm in sm.matches(info.reactants[i]):
            key = (rm.quality, rm.name_overlap, -(rm.geom_dev or 0.0))
            cur = best_for_res.get(id(rm.residue))
            if cur is None or key > cur[0]:
                best_for_res[id(rm.residue)] = (key, [i], rm)
            elif key == cur[0]:
                cur[1].append(i)
    assign = defaultdict(list)
    for key, idxs, _rm in best_for_res.values():
        for i in idxs:
            rm = next(m for m in sm.matches(info.reactants[i]) if m.residue is _rm.residue)
            assign[i].append(rm)
    # species with identical structure share residues
    for i in info.heavy_reactants:
        for k in info.heavy_reactants:
            if i != k and info.reactants[i].canon == info.reactants[k].canon and not assign.get(i) and assign.get(k):
                assign[i] = list(assign[k])
    return assign


def _candidates(sel: Sel, assign: dict, ligands: list) -> list:
    out = []
    if sel.species is not None:
        for rm in assign.get(sel.species, []):
            idx = sorted({m[a] for m in rm.matches for a in sel.atoms})
            out.extend((rm.residue, t) for t in idx)
    else:
        for lig in ligands:
            if lig.resname == sel.resname and sel.atomname in lig.name_index:
                out.append((lig, lig.name_index[sel.atomname]))
    return out


def measure_pair(pair: Pair, assign: dict, ligands: list):
    A = _candidates(pair.a, assign, ligands)
    B = _candidates(pair.b, assign, ligands)
    best = None
    for ra, ia in A:
        for rb, ib in B:
            if ra is rb:
                continue
            d = float(np.linalg.norm(ra.coords[ia] - rb.coords[ib]))
            if best is None or d < best[0]:
                best = (d, ra, ia, rb, ib)
    return best


def protein_meta(pkey: str, prot: dict) -> dict:
    orgs = prot.get("organisms") or {}
    refs = []
    for o, d in orgs.items():
        if isinstance(d, dict):
            ids = [str(d[k]) for k in ("uniprot_id", "rcsb_id", "rcsb") if d.get(k)]
            refs.append(f"{o}:{'/'.join(ids)}" if ids else str(o))
    return {
        "protein": pkey,
        "protein_name": prot.get("name", ""),
        "abbreviation": prot.get("abbreviation", ""),
        "gene_name": prot.get("gene_name", ""),
        "ec_number": prot.get("ec_number", ""),
        "compartment": prot.get("compartment", ""),
        "cofactors": ";".join(str(c) for c in (prot.get("cofactors") or [])),
        "oligomerization": str(prot.get("oligomerization", "") or ""),
        "oois": ";".join(str(o) for o in (prot.get("oois") or [])),
        "reference_organisms": "; ".join(refs),
    }


def reaction_meta(info: ReactionInfo, compounds: dict) -> dict:
    r = info.rxn
    subs = [str(x) for x in (r.get("substrates") or [])]
    prods = [str(x) for x in (r.get("products") or [])]
    heavy_smiles = ".".join(s.canon for s in info.reactants if s.n_heavy) + ">>" + \
        ".".join(s.canon for s in info.products if s.n_heavy)
    return {
        "reaction_id": info.rid,
        "reaction_ref": r.get("ref", ""),
        "substrates": ";".join(subs),
        "substrate_names": ";".join(_compound_name(compounds, s) or s for s in subs),
        "products": ";".join(prods),
        "product_names": ";".join(_compound_name(compounds, s) or s for s in prods),
        "reaction_smiles_heavy": heavy_smiles,
    }


def process_cif(path: Path, root: Path, cfg: dict, infos: dict, table: Optional[EnzymeTable],
                args) -> list[dict]:
    stem = path.stem
    job, model_idx = _job_and_model(stem)
    base = {"cif_file": str(path.relative_to(root)) if path.is_relative_to(root) else str(path),
            "cif_name": stem, "job_name": job, "boltz_model": model_idx}
    conf = _read_json(path.parent / f"confidence_{stem}.json")
    aff = _read_json(path.parent / f"affinity_{job}.json")
    try:
        models = parse_cif(path)
    except Exception as e:
        return [dict(base, status=f"error: {e}")]
    rx = re.search(args.name_regex, stem) if args.name_regex else None
    rgx = rx.groupdict() if rx else {}
    tokens = set(_tokens(stem))
    aliases = protein_aliases(cfg)
    all_rids = {i.rid.lower() for i in infos.values()}
    rows = []
    for mdl in models:
        row = dict(base)
        row["cif_model"] = mdl.model if len(models) > 1 else ""
        notes: list[str] = []
        seqs = [norm_seq(s) for s in mdl.chains.values()]
        row["protein_chains"] = ";".join(mdl.chains)
        row["protein_length"] = ";".join(str(len(s)) for s in mdl.chains.values())
        row["protein_seq_md5"] = hashlib.md5("|".join(sorted(set(seqs))).encode()).hexdigest()[:10] if seqs else ""
        # ---- enzyme lookup
        rows_idx, how = ([], "no_table")
        if table is not None:
            rows_idx, how = table.lookup(seqs, stem, rgx.get("enzyme"))
        row["enzyme_match"] = how
        # ---- protein family candidates
        fam = None
        fam_src = ""
        if rgx.get("protein"):
            fam = aliases.get(rgx["protein"].lower())
            fam_src = "name_regex"
        if fam is None and table is not None and rows_idx:
            fams = {table.row_family(k) for k in rows_idx} - {None}
            if len(fams) == 1:
                fam, fam_src = fams.pop(), "csv"
        if fam is None:
            hits = {aliases[t] for t in tokens if t in aliases}
            if len(hits) == 1:
                fam, fam_src = hits.pop(), "filename"
        cand = [i for i in infos.values() if (not args.protein or i.protein in args.protein)]
        if fam is not None:
            fc = [i for i in cand if i.protein == fam]
            if fc:
                cand = fc
            else:
                notes.append(f"protein family '{fam}' (from {fam_src}) has no analysable reaction")
                fam = None
        # ---- reaction hints
        rid_hints = set()
        if rgx.get("reaction"):
            rid_hints.add(rgx["reaction"].lower())
        rid_hints |= tokens & all_rids
        if table is not None and rows_idx:
            rid_hints |= {str(table.row_reaction(k)).lower() for k in rows_idx if table.row_reaction(k)}
        if not mdl.ligands:
            row["status"] = "no_ligands"
            rows.append(_finish_row(row, notes, conf, aff, table, rows_idx))
            continue
        sm = SpeciesMatcher(mdl.ligands)
        scored = sorted(((score_reaction(i, sm, rid_hints), i) for i in cand),
                        key=lambda t: t[0], reverse=True)
        if fam is not None and (not scored or scored[0][0][2] == 0):
            # docked ligands match no reaction of that family: try every protein in the YAML
            allc = [i for i in infos.values() if (not args.protein or i.protein in args.protein)]
            scored = sorted(((score_reaction(i, sm, rid_hints), i) for i in allc),
                            key=lambda t: t[0], reverse=True)
            if scored and scored[0][0][2] > 0:
                notes.append(f"ligands match no reaction of '{fam}' (from {fam_src}); "
                             f"using {scored[0][1].uid}")
                fam_src = "ligands(fallback)"
        if not scored or scored[0][0][2] == 0:
            row["status"] = "no_reaction_match"
            row["ligands_in_cif"] = "; ".join(f"{l.label}({_formula_str(l.formula)})" for l in mdl.ligands)
            rows.append(_finish_row(row, notes, conf, aff, table, rows_idx))
            continue
        top_score, info = scored[0]
        ties = [i for s, i in scored[1:] if s[:5] == top_score[:5]]
        if ties:
            same = [i.uid for i in ties if i.chem_key == info.chem_key]
            diff = [i.uid for i in ties if i.chem_key != info.chem_key]
            if same:
                notes.append(f"same chemistry also in {', '.join(same)}")
            if diff:
                notes.append(f"AMBIGUOUS reaction: also consistent with {', '.join(diff)} "
                             f"(add protein/reaction to the file name, the CSV, or use --name-regex)")
        geo_rivals = [(s, i) for s, i in scored[1:] if s[:4] == top_score[:4] and s[4] != top_score[4]
                      and i.chem_key != info.chem_key]
        if geo_rivals:
            s2, i2 = geo_rivals[0]
            notes.append(f"reaction chosen from 3D geometry (sp2/sp3 pattern) over {i2.uid}: angle "
                         f"deviation {-top_score[4]:.0f} vs {-s2[4]:.0f} deg - check if unexpected")
        if rid_hints and info.rid.lower() not in rid_hints:
            notes.append(f"reaction id hint {sorted(rid_hints)} not consistent with docked ligands")
        # narrow CSV rows with the reaction column if there were several
        if table is not None and len(rows_idx) > 1 and table.reaction_col:
            narrowed = [k for k in rows_idx if str(table.row_reaction(k)).lower() == info.rid.lower()]
            if narrowed:
                rows_idx = narrowed
        row.update(protein_meta(info.protein, cfg["proteins"][info.protein]))
        row.update(reaction_meta(info, cfg.get("compounds") or {}))
        row["protein_assignment"] = fam_src or "ligands"
        assign = assign_species(info, sm)
        lig_desc = []
        methods = []
        for lig in mdl.ligands:
            who = [(i, rm) for i, rms in assign.items() for rm in rms if rm.residue is lig]
            if who:
                i, rm = who[0]
                lig_desc.append(f"{lig.label}={info.reactants[i].name}[{rm.method}]")
                methods.append(rm.method)
            else:
                lig_desc.append(f"{lig.label}=?({_formula_str(lig.formula)})")
        row["ligands_in_cif"] = "; ".join(lig_desc)
        missing = [info.reactants[i].name for i in info.heavy_reactants if not assign.get(i)]
        row["substrates_not_docked"] = ";".join(missing)
        # ---- distances
        measured = []
        for p in info.pairs:
            res = measure_pair(p, assign, mdl.ligands)
            if res is not None:
                measured.append((p, res))
        if not info.pairs:
            row["status"] = "no_pair_defined"
            rows.append(_finish_row(row, notes + info.notes, conf, aff, table, rows_idx))
            continue
        if not measured:
            row["status"] = "partner_not_docked"
            rows.append(_finish_row(row, notes + info.notes, conf, aff, table, rows_idx))
            continue
        best_rank = min(p.rank for p, _ in measured)
        if best_rank != info.pairs[0].rank:
            notes.append(f"preferred pair type '{info.pairs[0].kind}' not measurable "
                         f"(substrate missing); using '{TIERS[best_rank[0]]}'")
        group = [(p, r) for p, r in measured if p.rank == best_rank]
        p, (d, ra, ia, rb, ib) = min(group, key=lambda t: t[1][0])
        row["status"] = "ok"
        row["pair_type"] = p.kind
        ends = []
        for tag, sel, res_, idx in (("atom1", p.a, ra, ia), ("atom2", p.b, rb, ib)):
            sp = info.reactants[sel.species] if sel.species is not None else None
            row[f"{tag}_compound"] = sp.key if sp is not None and sp.key else (sp.name if sp else sel.resname)
            row[f"{tag}_ligand"] = res_.resname
            row[f"{tag}_chain"] = res_.chain
            row[f"{tag}_cif_atom"] = res_.names[idx]
            row[f"{tag}_element"] = PT.GetElementSymbol(res_.elements[idx])
            row[f"{tag}_environment"] = sp.env(sel.atoms[0]) if sp is not None else ""
            ends.append(f"{sp.name if sp else res_.resname} {res_.names[idx]} ({res_.chain}:{res_.resname})")
        row["pair"] = "  <->  ".join(ends)
        row["distance_A"] = round(d, 3)
        row[f"within_{args.cutoff:g}A"] = bool(d <= args.cutoff)
        row["pymol_cmd"] = (f"distance aoi_{info.rid}, (chain {ra.chain} and name {ra.names[ia]}), "
                            f"(chain {rb.chain} and name {rb.names[ib]})")
        # other candidates: bond-formation and H-transfer pairs (+ fallback tier if that was used)
        listed = [(pp, r) for pp, r in measured if pp.tier <= 2 or pp.tier == p.tier]
        row["min_distance_any_candidate_A"] = round(min(r[0] for _, r in listed), 3)
        row["n_candidate_pairs"] = len(listed)
        cand_txt = []
        for pp, (dd, a1, i1, b1, i2) in sorted(listed, key=lambda t: (t[0].rank, t[1][0]))[:25]:
            cand_txt.append(f"{pp.kind}: {a1.resname}.{a1.names[i1]}-{b1.resname}.{b1.names[i2]}={dd:.2f}")
        row["candidate_pairs"] = " | ".join(cand_txt)
        row["mapping_method"] = ";".join(dict.fromkeys(methods))
        rows.append(_finish_row(row, notes + info.notes, conf, aff, table, rows_idx))
    return rows


def _finish_row(row, notes, conf, aff, table, rows_idx):
    if table is not None and rows_idx:
        row["enzyme_id"] = "|".join(dict.fromkeys(table.row_id(k) for k in rows_idx))
        row.update(table.row_dict(rows_idx[0]))
        if len(rows_idx) > 1:
            notes = notes + [f"{len(rows_idx)} CSV rows match this sequence"]
    row.update(conf)
    row.update(aff)
    row["notes"] = " | ".join(dict.fromkeys(n for n in notes if n))
    return row


# --------------------------------------------------------------------------------------
# outputs
# --------------------------------------------------------------------------------------
LONG_ORDER = [
    "cif_file", "cif_name", "job_name", "boltz_model", "cif_model", "status",
    "enzyme_id", "enzyme_match", "protein_chains", "protein_length", "protein_seq_md5",
    "protein", "protein_name", "abbreviation", "gene_name", "ec_number", "compartment",
    "cofactors", "oligomerization", "oois", "reference_organisms", "protein_assignment",
    "reaction_id", "reaction_ref", "substrates", "substrate_names", "products", "product_names",
    "ligands_in_cif", "substrates_not_docked", "pair_type", "pair",
    "atom1_compound", "atom1_ligand", "atom1_chain", "atom1_cif_atom", "atom1_element", "atom1_environment",
    "atom2_compound", "atom2_ligand", "atom2_chain", "atom2_cif_atom", "atom2_element", "atom2_environment",
    "distance_A", "WITHIN", "min_distance_any_candidate_A", "n_candidate_pairs", "candidate_pairs",
    "mapping_method", "pymol_cmd",
]


def order_columns(df: pd.DataFrame) -> pd.DataFrame:
    if "cif_model" in df.columns and not df["cif_model"].fillna("").astype(str).str.len().any():
        df = df.drop(columns=["cif_model"])
    within = [c for c in df.columns if c.startswith("within_")]
    order = []
    for c in LONG_ORDER:
        order.extend(within if c == "WITHIN" else [c])
    first = [c for c in order if c in df.columns]
    csvc = sorted(c for c in df.columns if c.startswith("csv_"))
    boltz = [c for c in df.columns if c.startswith("boltz_") and c != "boltz_model"]
    tail = [c for c in ("reaction_smiles_heavy", "notes") if c in df.columns]
    other = [c for c in df.columns if c not in set(first) | set(csvc) | set(boltz) | set(tail)]
    idx = first.index("protein") if "protein" in first else len(first)
    cols = first[:idx] + csvc + first[idx:] + boltz + other + tail
    return df[cols]


def wide_table(df: pd.DataFrame, cfg: dict) -> pd.DataFrame:
    if "protein" not in df.columns or "reaction_id" not in df.columns:
        return pd.DataFrame()
    ok = df[df["protein"].notna() & df["reaction_id"].notna()].copy()
    if ok.empty:
        return pd.DataFrame()
    # one row per distinct protein sequence (an "enzyme copy") and protein family
    blank = pd.Series("", index=ok.index)
    enz = ok["enzyme_id"].fillna("").astype(str) if "enzyme_id" in ok.columns else blank
    md5 = ok["protein_seq_md5"].fillna("").astype(str) if "protein_seq_md5" in ok.columns else blank
    ok["_enzyme"] = ("seq:" + md5).where(md5 != "", "id:" + enz.where(enz != "", ok["job_name"].astype(str)))
    ok["_enzyme_id"] = enz
    ok["_model"] = pd.to_numeric(ok.get("boltz_model"), errors="coerce").fillna(10**6)
    rid_order = []
    for p in (cfg.get("proteins") or {}).values():
        for r in p.get("reactions") or []:
            if str(r.get("id")) not in rid_order:
                rid_order.append(str(r.get("id")))
    multi = ok.groupby(["_enzyme", "protein", "reaction_id"]).size().max() > 1
    keep_cols = [c for c in ["enzyme_id", "protein_length"] if c in ok.columns] + \
        sorted(c for c in ok.columns if c.startswith("csv_")) + \
        [c for c in ["protein", "protein_name", "abbreviation", "gene_name", "ec_number",
                     "compartment", "cofactors", "oligomerization"] if c in ok.columns]
    out_rows = []
    for (enz, prot), g in ok.groupby(["_enzyme", "protein"], sort=True):
        g = g.sort_values(["_model", "cif_name"])
        ids = [x for x in dict.fromkeys(g["_enzyme_id"]) if x]
        rec = {"enzyme_id": "|".join(ids) if ids else enz}
        for c in keep_cols:
            if c == "enzyme_id":
                continue
            vals = [str(v) for v in dict.fromkeys(g[c].dropna()) if str(v) != ""]
            rec[c] = vals[0] if len(vals) == 1 else "|".join(vals)
        for rid in rid_order:
            gr = g[g["reaction_id"] == rid]
            if gr.empty:
                continue
            top = gr.iloc[0]
            rec[f"{rid}_substrates"] = top.get("substrates", "")
            rec[f"{rid}_pair"] = top.get("pair", "")
            rec[f"{rid}_distance_A"] = top.get("distance_A")
            if multi:
                dists = pd.to_numeric(gr.get("distance_A"), errors="coerce")
                rec[f"{rid}_min_distance_all_models_A"] = dists.min()
                rec[f"{rid}_n_models"] = int(dists.notna().sum())
            if "boltz_confidence_score" in gr.columns:
                rec[f"{rid}_confidence_score"] = top.get("boltz_confidence_score")
            if "boltz_affinity_pred_value" in gr.columns:
                rec[f"{rid}_affinity_pred_value"] = top.get("boltz_affinity_pred_value")
            rec[f"{rid}_status"] = top.get("status")
        out_rows.append(rec)
    wt = pd.DataFrame(out_rows)
    porder = {k: n for n, k in enumerate(cfg.get("proteins") or {})}
    wt = wt.assign(_p=wt["protein"].map(porder)).sort_values(["_p", "enzyme_id"]).drop(columns="_p")
    cols = ["enzyme_id"] + [c for c in keep_cols if c != "enzyme_id"]
    for rid in rid_order:
        cols += [c for c in wt.columns if c.startswith(f"{rid}_") and c not in cols]
    cols += [c for c in wt.columns if c not in cols]
    wt = wt[[c for c in cols if c in wt.columns]]
    for c in wt.columns:
        if c.endswith("_n_models"):
            wt[c] = wt[c].astype("Int64")
    return wt


def reactions_table(infos: dict, cfg: dict) -> pd.DataFrame:
    rows = []
    compounds = cfg.get("compounds") or {}
    for info in infos.values():
        best = info.best_pairs()
        prot = cfg["proteins"][info.protein]
        rec = {"protein": info.protein, "protein_name": prot.get("name", "")}
        rec.update(reaction_meta(info, compounds))
        rec["reaction_field_used"] = info.source_field
        rec["pair_type"] = best[0].kind if best else "none"
        rec["atoms_of_interest"] = " || ".join(
            f"{info.sel_label(p.a)}  <->  {info.sel_label(p.b)}" for p in best)
        for tag, attr in (("atom1", "a"), ("atom2", "b")):
            sels = [getattr(p, attr) for p in best]
            if not sels:
                continue
            s0 = sels[0]
            if s0.species is None:
                rec[f"{tag}_compound"] = s0.resname
                rec[f"{tag}_boltz_atom_names"] = s0.atomname
                continue
            sp = info.reactants[s0.species]
            atoms = sorted({a for s in sels if s.species == s0.species for a in s.atoms})
            rec[f"{tag}_compound"] = sp.key or sp.name
            rec[f"{tag}_boltz_atom_names"] = "|".join(sp.atom_label(a) for a in atoms)
            rec[f"{tag}_smiles_heavy_atom_index"] = "|".join(str(a) for a in atoms)
            rec[f"{tag}_environment"] = " ; ".join(dict.fromkeys(sp.env(a) for a in atoms))
        others = [p for p in info.pairs if p not in best and p.tier <= 2]
        rec["other_candidate_pairs"] = " || ".join(
            f"[{p.kind}] {info.sel_label(p.a, False)} <-> {info.sel_label(p.b, False)}" for p in others[:12])
        ev = info.events
        lab = lambda x: f"{info.reactants[x[0]].name}:{info.reactants[x[0]].atom_label(x[1])}"  # noqa: E731
        rec["bonds_formed"] = "; ".join(f"{lab(a)}-{lab(b)}" for a, b, _ in ev.formed_inter + ev.formed_intra)
        rec["bonds_broken"] = "; ".join(f"{lab(a)}-{lab(b)}" for a, b, _ in ev.broken)
        rec["H_changes"] = "; ".join(f"{lab(x)} {d:+d}H" for x, d in sorted(ev.dH.items()) if d)
        rec["notes"] = " | ".join(info.notes)
        rows.append(rec)
    df = pd.DataFrame(rows)
    front = ["protein", "protein_name", "reaction_id", "reaction_ref", "substrates", "substrate_names",
             "pair_type", "atoms_of_interest",
             "atom1_compound", "atom1_boltz_atom_names", "atom1_smiles_heavy_atom_index", "atom1_environment",
             "atom2_compound", "atom2_boltz_atom_names", "atom2_smiles_heavy_atom_index", "atom2_environment",
             "other_candidate_pairs", "bonds_formed", "bonds_broken", "H_changes",
             "products", "product_names", "reaction_smiles_heavy", "reaction_field_used", "notes"]
    return df[[c for c in front if c in df.columns] + [c for c in df.columns if c not in front]]


def draw_reaction(info: ReactionInfo, path: Path) -> None:
    best = info.best_pairs()
    hl = defaultdict(set)
    for p in best:
        for s in (p.a, p.b):
            if s.species is not None:
                hl[s.species].update(s.atoms)
    mols, legends, hls = [], [], []
    for i in info.heavy_reactants:
        sp = info.reactants[i]
        m = Chem.Mol(sp.mol)
        for a in hl.get(i, ()):
            bn = sp.boltz_name(a)
            m.GetAtomWithIdx(a).SetProp("atomNote", bn or sp.atom_label(a))
        rdDepictor.Compute2DCoords(m)
        mols.append(m)
        # RDKit's drawing font has no Greek/sub-/superscripts: fall back to the YAML key
        legends.append(sp.name if sp.name.isascii() else (sp.key or sp.name))
        hls.append(sorted(hl.get(i, ())))
    if not mols:
        return
    w, h = 380, 320
    d = rdMolDraw2D.MolDraw2DSVG(w * len(mols), h + 40, w, h)
    opts = d.drawOptions()
    opts.legendFontSize = 18
    opts.annotationFontScale = 0.9
    opts.highlightRadius = 0.45
    opts.fillHighlights = True
    colors = [{a: (1.0, 0.55, 0.1) for a in hh} for hh in hls]
    d.DrawMolecules(mols, legends=legends, highlightAtoms=hls, highlightAtomColors=colors,
                    highlightBonds=[[] for _ in mols])
    d.FinishDrawing()
    svg = d.GetDrawingText()
    title = (f"{info.uid}: {TIERS[best[0].tier] if best else 'no pair'} - " +
             " / ".join(f"{info.sel_label(p.a, False)} <-> {info.sel_label(p.b, False)}" for p in best))
    svg = svg.replace("</svg>", f"<text x='10' y='{h + 30}' font-family='sans-serif' font-size='15'>"
                      f"{html.escape(title)}</text></svg>")
    path.write_text(svg)


def write_depictions(infos: dict, outdir: Path) -> Optional[Path]:
    """One SVG per reaction in outdir/ plus a self-contained atoms_of_interest.html next to it."""
    outdir.mkdir(parents=True, exist_ok=True)
    blocks = []
    for info in infos.values():
        fn = outdir / f"{info.protein}_{info.rid}.svg"
        try:
            draw_reaction(info, fn)
        except Exception as e:  # drawing is optional
            log.warning("could not draw %s: %s", info.uid, e)
            continue
        svg = re.sub(r"^<\?xml[^>]*\?>\s*", "", fn.read_text())
        best = info.best_pairs()
        kind = best[0].kind if best else "no intermolecular pair"
        note = f"<p class='note'>{html.escape(' | '.join(info.notes))}</p>" if info.notes else ""
        blocks.append(f"<section><h3>{html.escape(info.uid)} <small>{html.escape(kind)}</small></h3>"
                      f"<div class='svg'>{svg}</div>{note}</section>")
    page = outdir.parent / "atoms_of_interest.html"
    page.write_text(
        "<!doctype html><html><head><meta charset='utf-8'><title>Atoms of interest</title>"
        "<style>body{font-family:system-ui,sans-serif;margin:16px;background:#fff;color:#111}"
        "section{border-top:1px solid #ddd;padding:8px 0}.svg{overflow-x:auto}"
        "small{color:#666;font-weight:normal}.note{color:#8a4b00;font-size:13px}</style></head><body>"
        "<h2>Atoms of interest per reaction</h2><p>Highlighted atoms are the pair measured in the "
        "CIFs; labels are the expected Boltz atom names (element + canonical rank).</p>"
        + "\n".join(blocks) + "</body></html>")
    return page


# --------------------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------------------
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[1] if __doc__ else "",
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--yaml", required=True, type=Path, help="pathway YAML with proteins/reactions")
    ap.add_argument("--cif-dir", type=Path, help="folder with Boltz-2 .cif files (recursive)")
    ap.add_argument("--csv", type=Path, help="CSV/TSV/FASTA with the docked enzyme sequences")
    ap.add_argument("-o", "--out", type=Path, default=Path("distance_atoms_of_interest.csv"),
                    help="output CSV (or a folder); default distance_atoms_of_interest.csv")
    ap.add_argument("--pattern", default="*.cif", help="glob for structure files (default *.cif)")
    ap.add_argument("--seq-col", help="CSV column with the protein sequence (auto-detected)")
    ap.add_argument("--id-col", help="CSV column with the enzyme/sequence ID (auto-detected)")
    ap.add_argument("--protein-col", help="CSV column with the protein family key, e.g. asmt")
    ap.add_argument("--reaction-col", help="CSV column with the reaction id (rr, r1, ...)")
    ap.add_argument("--sep", help="CSV delimiter (auto-detected; e.g. ';' or '\\t')")
    ap.add_argument("--protein", nargs="+", help="only consider these protein keys from the YAML")
    ap.add_argument("--name-regex", help="regex applied to the CIF file name; named groups "
                    "(?P<protein>..) (?P<reaction>..) (?P<enzyme>..) are used as hints")
    ap.add_argument("--cutoff", type=float, default=4.0,
                    help="flag poses with distance <= cutoff in a within_<cutoff>A column (default 4.0)")
    ap.add_argument("--mcs-timeout", type=float, default=10.0, help="seconds per MCS (default 10)")
    ap.add_argument("--analyze-only", action="store_true",
                    help="only analyse the reactions (write reaction_atoms_of_interest.csv)")
    ap.add_argument("--no-draw", action="store_true", help="do not write SVG depictions")
    ap.add_argument("-v", "--verbose", action="store_true")
    args = ap.parse_args(argv)
    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO,
                        format="%(levelname)s: %(message)s")

    cfg = yaml.safe_load(args.yaml.read_text())
    if not isinstance(cfg, dict) or "proteins" not in cfg:
        log.error("%s has no 'proteins' section", args.yaml)
        return 2
    cfg["proteins"] = {str(k): v for k, v in cfg["proteins"].items()}
    if args.protein:
        unknown = [p for p in args.protein if p not in cfg["proteins"]]
        if unknown:
            log.error("unknown --protein %s (known: %s)", unknown, list(cfg["proteins"]))
            return 2
    out = args.out
    if out.suffix.lower() != ".csv":
        out = out / "distance_atoms_of_interest.csv"
    out.parent.mkdir(parents=True, exist_ok=True)

    log.info("analysing reactions in %s ...", args.yaml.name)
    infos = analyze_config(cfg, args.mcs_timeout, args.protein)
    rt = reactions_table(infos, cfg)
    rt_path = out.parent / "reaction_atoms_of_interest.csv"
    rt.to_csv(rt_path, index=False)
    for info in infos.values():
        best = info.best_pairs()
        txt = " || ".join(f"{info.sel_label(p.a, False)} <-> {info.sel_label(p.b, False)}" for p in best) or "-"
        log.info("  %-14s %-15s %s", info.uid, best[0].kind if best else "none", txt)
    if not args.no_draw:
        write_depictions(infos, out.parent / "aoi_depictions")
    log.info("wrote %s", rt_path)
    if args.analyze_only:
        return 0
    if not args.cif_dir:
        log.error("--cif-dir is required unless --analyze-only")
        return 2

    table = None
    if args.csv:
        sep = args.sep.encode().decode("unicode_escape") if args.sep else None
        table = EnzymeTable(args.csv, cfg, args.seq_col, args.id_col, args.protein_col, args.reaction_col, sep)
    files = sorted(p for p in args.cif_dir.rglob(args.pattern) if p.is_file())
    if not files:
        log.error("no files matching %s under %s", args.pattern, args.cif_dir)
        return 2
    log.info("processing %d structure files ...", len(files))
    rows = []
    for k, f in enumerate(files, 1):
        try:
            rows.extend(process_cif(f, args.cif_dir, cfg, infos, table, args))
        except Exception as e:  # never stop the whole run for one file
            log.exception("failed on %s", f)
            rows.append({"cif_file": str(f), "cif_name": f.stem, "status": f"error: {e}"})
        if k % 200 == 0:
            log.info("  %d/%d", k, len(files))
    df = order_columns(pd.DataFrame(rows))
    for c in ("boltz_model", "n_candidate_pairs"):
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce").astype("Int64")
    df.to_csv(out, index=False)
    wt = wide_table(df, cfg)
    wide_path = out.with_name(out.stem + "_wide.csv")
    wt.to_csv(wide_path, index=False)
    counts = df["status"].value_counts().to_dict()
    log.info("status: %s", ", ".join(f"{k}={v}" for k, v in counts.items()))
    if table is not None and "enzyme_match" in df.columns:
        n_none = int((df["enzyme_match"] == "none").sum())
        if n_none:
            log.warning("%d structure(s) could not be matched to a row of %s (sequence/ID); "
                        "check --seq-col/--id-col", n_none, args.csv.name)
    log.info("wrote %s (%d rows) and %s (%d rows)", out, len(df), wide_path, len(wt))
    return 0


if __name__ == "__main__":
    sys.exit(main())