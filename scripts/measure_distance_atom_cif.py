"""
measure_distance_atom_cif.py
─────────────────────────────────────────────────────────────────────────────
Measures the distance (Ångströms) between ONE atom of a docked ligand and ONE
atom of its docked cofactor, across a whole folder of docked structures.

===============================================================================
  HOW TO USE — TWO STEPS, NOTHING TO TYPE ON THE COMMAND LINE
===============================================================================

  Everything you need to fill in sits in ONE clearly marked box near the top
  of this file, right below the imports:

      ██  INSERT YOUR DATA HERE  ██

  STEP 1   Set MODE = "show_atoms" in that box, then run:

               python measure_distance_atom_cif.py

           It prints the atom numbering of both SMILES together with the
           matching atom name found in your structures (e.g. which SMILES
           position is the sulfonium sulfur S1 of SAM). Read the two numbers
           you need from that table.

  STEP 2   Put those two numbers into LIGAND_ATOM and COFACTOR_ATOM in the
           box, set MODE = "measure", and run the same command again:

               python measure_distance_atom_cif.py

           A CSV with the distance in Ångströms is written to OUTPUT_DIR,
           one row per structure in CIF_FOLDER.

===============================================================================
  WHY NO ATOM NAMES ARE NEEDED
===============================================================================

Docked mmCIF/PDB files contain no bond connectivity, and the atom names are
program-generated (N21, C13, S1 …), so "the 10th atom of my SMILES" is not the
atom called C10 in the file. This script therefore locates each molecule by
*graph matching*: it rebuilds the molecule's bonds from the coordinates
(covalent radii) and finds the one arrangement that reproduces the element and
connectivity of your SMILES. You only ever refer to atoms by their position in
the SMILES you supply.

===============================================================================
  COMMAND LINE (optional — overrides the box)
===============================================================================

Every field can also be passed as a flag, which is handy for one-off runs:

    python measure_distance_atom_cif.py \\
        --cif-dir       /path/to/docked_structures \\
        --ligand-smiles '[H]c1c([H]...'  --ligand-atom 9 \\
        --cofactor-smiles '[H]O[C@@]1...' --cofactor-atom 18 \\
        --out-dir       /path/to/results

    --print-atoms   same as MODE = "show_atoms"

===============================================================================
  NOTES
===============================================================================

* Paste the SMILES exactly as you docked it. Explicit [H] atoms are stripped
  automatically, because docked coordinates are heavy-atom only.
* Formal charges are ignored while matching (they are absent from the
  coordinates), which is what makes sulfonium SAM match.
* A molecule that is too symmetric for a unique match is reported in the
  `status` column rather than silently guessed.
* Structures that do not contain one of the two molecules still get a row,
  with the reason, so the CSV always covers the whole folder.
─────────────────────────────────────────────────────────────────────────────
"""

from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass
from pathlib import Path

import gemmi
import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger

sys.path.insert(0, str(Path(__file__).resolve().parent))
from logger import build_logger  # noqa: E402

RDLogger.DisableLog("rdApp.*")
log = build_logger("atom_distance")

# Structure file extensions understood by --cif-dir
STRUCTURE_SUFFIXES = (".cif", ".mmcif", ".pdb", ".ent", ".cif.gz", ".pdb.gz")

# Covalent-radius slack (Å) used when rebuilding bonds from coordinates.
# 0.45 reproduces RDKit's DetermineConnectivity on the ligands in this project.
DEFAULT_BOND_TOLERANCE = 0.45

CSV_COLUMNS = [
    "structure",
    "ligand_resname", "ligand_chain", "ligand_seq",
    "ligand_atom", "ligand_element", "ligand_smiles_index",
    "cofactor_resname", "cofactor_chain", "cofactor_seq",
    "cofactor_atom", "cofactor_element", "cofactor_smiles_index",
    "distance_A", "status",
]


# ╔════════════════════════════════════════════════════════════════════════════╗
# ║                                                                            ║
# ║   ██  I N S E R T   Y O U R   D A T A   H E R E  ██                        ║
# ║                                                                            ║
# ║   Fill in the five numbered steps below, then just run:                      ║
# ║                                                                            ║
# ║        python measure_distance_atom_cif.py                                   ║
# ║                                                                            ║
# ╚════════════════════════════════════════════════════════════════════════════╝

# ── STEP 1 ─────────────────────────────────────────────────────────────────────
# Start with "show_atoms": the script prints the atom numbering of both SMILES
# and the matching atom name inside your structures (e.g. which SMILES position
# is the sulfonium sulfur S1 of SAM). Read the two numbers you need from that
# table, fill them into STEP 5, then come back and set this to "measure".
MODE = "show_atoms"                 # ← "show_atoms"  or  "measure"

# ── STEP 2 ─────────────────────────────────────────────────────────────────────
# The folder that holds your docked structures (.cif and/or .pdb files).
# Every structure inside it gets one row in the CSV.
CIF_FOLDER = r""                   # ← e.g. r"/home/pedro/docked_structures"

# ── STEP 3 ─────────────────────────────────────────────────────────────────────
# SMILES of the LIGAND, pasted exactly as you docked it.
# Keep the explicit [H] atoms if your docking input had them - they are
# stripped automatically.
LIGAND_SMILES = r""                # ← e.g. r"[H]c1c([H])c([H])c2c(c1[H])..."

# ── STEP 4 ─────────────────────────────────────────────────────────────────────
# SMILES of the COFACTOR, pasted exactly as you docked it.
COFACTOR_SMILES = r""              # ← e.g. r"[H]O[C@@]1([H])... [S@@+](...)"

# ── STEP 5 ─────────────────────────────────────────────────────────────────────
# Which atom of each SMILES to measure, counting 0, 1, 2, … through the heavy
# atoms from left to right in the SMILES above.
# Take these numbers from the table printed by MODE = "show_atoms".
LIGAND_ATOM = None                  # ← put the ligand atom number here
COFACTOR_ATOM = None                # ← put the cofactor atom number here

# ── STEP 6 ─────────────────────────────────────────────────────────────────────
# Where the resulting CSV is written. The folder is created if missing.
OUTPUT_DIR = r""                   # ← e.g. r"/home/pedro/results"
OUTPUT_FILENAME = "atom_distances.csv"

# ═════════════════════════════════════════════════════════════════════════════
#   OPTIONAL SETTINGS - safe to ignore, the defaults are what you normally want
# ═════════════════════════════════════════════════════════════════════════════
RECURSIVE = True          # also look inside sub-folders of CIF_FOLDER
ALL_PAIRS = False         # True = one row per ligand×cofactor copy, not just
                          #        the closest pair
MODEL_INDEX = 0           # which model to read from multi-model files
BOND_TOLERANCE = 0.45     # Å; nudge up if a ligand refuses to match

# ═════════════════════════════════════════════════════════════════════════════
#                       END OF USER DATA - NOTHING TO EDIT BELOW
# ═════════════════════════════════════════════════════════════════════════════


# ──────────────────────────────────────────────────────────────────────────────
# SMILES handling
# ──────────────────────────────────────────────────────────────────────────────

def parse_smiles(smiles: str, label: str) -> Chem.Mol:
    """Parse *smiles* and return a heavy-atom-only RDKit molecule.

    Explicit ``[H]`` atoms are removed, because docked coordinates are
    heavy-atom only. Atom indices therefore refer to the heavy atoms in the
    order they appear in the SMILES string.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Could not parse {label} SMILES: {smiles}")

    n_heavy = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1)
    if n_heavy == 0:
        raise ValueError(f"{label} SMILES has no heavy atoms: {smiles}")

    return Chem.RemoveHs(mol)


def select_smiles_atom(
    mol: Chem.Mol,
    label: str,
    index: int | None = None,
    map_num: int | None = None,
    atom_name: str | None = None,
    optional: bool = False,
) -> int | None:
    """Resolve which heavy atom of *mol* the user selected.

    Exactly one of *index* / *map_num* must be given, unless *atom_name* is
    used — then no SMILES atom is needed at all and ``None`` is returned.
    With *optional*, returning ``None`` is allowed (used by "show_atoms",
    where the numbers have not been chosen yet).
    """
    if atom_name is not None and index is None and map_num is None:
        return None

    if (index is None) == (map_num is None):
        if optional and index is None and map_num is None:
            return None
        raise ValueError(
            f"{label}: give exactly one of --{label}-atom or --{label}-atom-map"
        )

    if map_num is not None:
        hits = [
            atom.GetIdx()
            for atom in mol.GetAtoms()
            if atom.GetAtomMapNum() == map_num
        ]
        if not hits:
            raise ValueError(
                f"{label}: no atom carries the map number {map_num}. "
                f"Embed it in the SMILES as [C:{map_num}] / [N:{map_num}] and retry."
            )
        if len(hits) > 1:
            raise ValueError(
                f"{label}: map number {map_num} appears on {len(hits)} atoms; "
                "it must be unique."
            )
        return hits[0]

    if not 0 <= index < mol.GetNumAtoms():
        raise ValueError(
            f"{label}: atom index {index} is out of range 0-{mol.GetNumAtoms() - 1}. "
            f"Run with --print-atoms to list the valid indices."
        )
    return index


def indexed_smiles(mol: Chem.Mol) -> str:
    """Return the SMILES with every heavy atom labelled ``symbol:index+1``.

    RDKit treats map number 0 as "unset", so the printed number is the
    heavy-atom index plus one.
    """
    labelled = Chem.Mol(mol)
    for atom in labelled.GetAtoms():
        atom.SetAtomMapNum(atom.GetIdx() + 1)
    return Chem.MolToSmiles(labelled)


# ──────────────────────────────────────────────────────────────────────────────
# Structure handling
# ──────────────────────────────────────────────────────────────────────────────

@dataclass
class ResidueInstance:
    """One residue copy found in a structure file."""
    chain:     str
    resname:   str
    seq:       str
    atom_names: list[str]
    elements:  list[str]
    coords:    np.ndarray      # shape (n_atoms, 3)

    @property
    def signature(self) -> tuple[str, ...]:
        """Sorted element multiset — a cheap pre-filter before graph matching."""
        return tuple(sorted(self.elements))


def _is_hydrogen(atom: gemmi.Atom) -> bool:
    return atom.is_hydrogen()


def iter_residues(
    structure: gemmi.Structure,
    model_index: int = 0,
) -> list[ResidueInstance]:
    """Collect every residue of *model_index* as a :class:`ResidueInstance`.

    Waters and ions are kept (they are legitimate distance targets); altloc
    duplicates are collapsed onto the blank/'A' conformer.
    """
    if not len(structure):
        raise ValueError("structure contains no models")

    model = structure[min(model_index, len(structure) - 1)]
    residues: list[ResidueInstance] = []

    for chain in model:
        for residue in chain:
            names: list[str] = []
            elements: list[str] = []
            positions: list[tuple[float, float, float]] = []

            for atom in residue:
                altloc = atom.altloc
                if altloc not in ("\0", "", "A"):
                    continue
                if _is_hydrogen(atom):
                    continue

                element = atom.element
                if element.atomic_number == 0:
                    # Element not annotated (common in PDB files) — infer from name
                    element = gemmi.find_tabulated_residue(residue.name).find_atom(
                        atom.name, gemmi.ElType.UNKNOWN
                    )
                if element.atomic_number == 0:
                    element = _element_from_name(atom.name)

                pos = atom.pos
                names.append(atom.name)
                elements.append(element.name.capitalize())
                positions.append((pos.x, pos.y, pos.z))

            if not names:
                continue

            residues.append(ResidueInstance(
                chain=chain.name,
                resname=residue.name,
                seq=str(residue.seqid.num) if residue.seqid.num is not None else ".",
                atom_names=names,
                elements=elements,
                coords=np.asarray(positions, dtype=float),
            ))

    return residues


def _element_from_name(atom_name: str) -> gemmi.Element:
    """Best-effort element guess from a PDB atom name (e.g. ``CA1`` → C)."""
    for char in atom_name:
        if char.isalpha():
            try:
                return gemmi.Element(char.upper())
            except Exception:
                continue
    return gemmi.Element("C")


def build_connectivity(
    coords: np.ndarray,
    elements: list[str],
    tolerance: float = DEFAULT_BOND_TOLERANCE,
) -> np.ndarray:
    """Rebuild the bond list from coordinates using covalent radii.

    Two atoms are bonded when their separation is below the sum of their
    covalent radii plus *tolerance*. Returns a boolean adjacency matrix.
    """
    n = len(elements)
    if n == 1:
        return np.zeros((1, 1), dtype=bool)

    periodic = Chem.GetPeriodicTable()
    try:
        radii = np.array([
            periodic.GetRb0(periodic.GetAtomicNumber(el)) for el in elements
        ])
    except Exception:
        radii = np.full(n, 1.5)

    deltas = coords[:, None, :] - coords[None, :, :]
    distances = np.sqrt(np.einsum("ijk,ijk->ij", deltas, deltas))
    bonded = (distances < (radii[:, None] + radii[None, :] + tolerance))
    np.fill_diagonal(bonded, False)
    return bonded


def _neutralize_for_topology(mol: Chem.Mol) -> Chem.Mol:
    """Strip everything a graph match should ignore.

    Only element and connectivity identify an atom here, so formal charge
    (a sulfonium S+ in SAM, a carboxylate, …), isotope, explicit H counts,
    aromaticity and bond order must not influence the match — the SMILES
    carries those, the coordinates cannot.
    """
    for atom in mol.GetAtoms():
        atom.SetFormalCharge(0)
        atom.SetIsotope(0)
        atom.SetNumExplicitHs(0)
        atom.SetNoImplicit(True)
        atom.SetIsAromatic(False)
        atom.SetNumRadicalElectrons(0)
    for bond in mol.GetBonds():
        bond.SetBondType(Chem.BondType.SINGLE)
        bond.SetIsAromatic(False)
    return mol


def _graph_to_rdkit(adjacency: np.ndarray, elements: list[str]) -> Chem.Mol:
    """Turn a boolean adjacency matrix into a bond-order-free RDKit molecule."""
    editable = Chem.RWMol()
    for element in elements:
        editable.AddAtom(Chem.Atom(element))

    n = len(elements)
    for i in range(n):
        for j in range(i + 1, n):
            if adjacency[i, j]:
                editable.AddBond(i, j, Chem.BondType.SINGLE)

    return _neutralize_for_topology(editable.GetMol())


def _topology_query(mol: Chem.Mol) -> Chem.Mol:
    """Reduce the SMILES to a pure element-and-connectivity query."""
    return _neutralize_for_topology(Chem.Mol(mol))


def match_residue(
    query: Chem.Mol,
    residue: ResidueInstance,
    tolerance: float = DEFAULT_BOND_TOLERANCE,
) -> list[tuple[int, ...]] | None:
    """All ways *residue* can be superimposed onto *query* by topology.

    Returns a list of tuples where ``match[i]`` is the residue index that
    corresponds to query atom ``i``, or ``None`` if the residue is not the
    same molecule.
    """
    if residue.signature != tuple(sorted(a.GetSymbol() for a in query.GetAtoms())):
        return None

    adjacency = build_connectivity(residue.coords, residue.elements, tolerance)
    target = _graph_to_rdkit(adjacency, residue.elements)

    try:
        matches = target.GetSubstructMatches(
            _topology_query(query), uniquify=False, maxMatches=64
        )
    except Exception:
        return None

    return list(matches) or None


@dataclass
class Placement:
    """One located copy of a molecule, plus the atom the user selected."""
    residue:  ResidueInstance
    atom_idx: int              # index into residue.atom_names
    atom_name: str
    element: str
    smiles_index: int | None   # None when --*-atom-name was used
    ambiguous: bool


def locate_molecule(
    residues: list[ResidueInstance],
    query: Chem.Mol,
    smiles_index: int | None,
    label: str,
    tolerance: float = DEFAULT_BOND_TOLERANCE,
    resname: str | None = None,
    atom_name: str | None = None,
) -> list[Placement]:
    """Find every copy of *query* in *residues* and resolve the selected atom.

    *resname* restricts the search to one residue name. *atom_name* picks the
    atom directly (the SMILES is still used to confirm which residues are the
    right molecule). When *atom_name* is ``None``, the selected atom is located
    through the graph match, and a molecule that is symmetric enough for the
    match to be non-unique is flagged as ambiguous.
    """
    wanted_elements = tuple(sorted(a.GetSymbol() for a in query.GetAtoms()))
    placements: list[Placement] = []

    for residue in residues:
        if resname and residue.resname != resname:
            continue
        if residue.signature != wanted_elements:
            continue

        matches = match_residue(query, residue, tolerance)
        if not matches:
            continue

        ambiguous = False

        if atom_name is not None:
            if atom_name not in residue.atom_names:
                continue
            candidate = residue.atom_names.index(atom_name)
        else:
            if smiles_index is None:
                raise ValueError(
                    f"{label}: no atom selected — pass --{label}-atom, "
                    f"--{label}-atom-map or --{label}-atom-name"
                )
            # The selected atom can land on different positions when the
            # molecule is symmetric; every such position is chemically valid,
            # so take the lowest index and flag it for the user to verify.
            candidates = {match[smiles_index] for match in matches}
            if not candidates:
                continue
            candidate = min(candidates)
            ambiguous = len(candidates) > 1

        placements.append(Placement(
            residue=residue,
            atom_idx=candidate,
            atom_name=residue.atom_names[candidate],
            element=residue.elements[candidate],
            smiles_index=smiles_index,
            ambiguous=ambiguous,
        ))

    if not placements:
        where = f" with resname {resname}" if resname else ""
        raise ValueError(f"no copy of the {label} found{where}")

    return placements


# ──────────────────────────────────────────────────────────────────────────────
# Measurement
# ──────────────────────────────────────────────────────────────────────────────

def measure(placements_a: list[Placement], placements_b: list[Placement],
            all_pairs: bool = False) -> list[tuple[Placement, Placement, float]]:
    """Return ``(a, b, distance)`` triples in Ångströms.

    By default only the closest ligand/cofactor pair is kept, so a structure
    holding several binding sites still yields a single row. With
    *all_pairs* every combination is returned.
    """
    combinations = [
        (a, b, float(np.linalg.norm(
            a.residue.coords[a.atom_idx] - b.residue.coords[b.atom_idx]
        )))
        for a in placements_a
        for b in placements_b
    ]

    if not all_pairs:
        return [min(combinations, key=lambda item: item[2])]

    return sorted(combinations, key=lambda item: item[2])


def describe(placement: Placement) -> tuple[str, str, str, str, str, int | None]:
    """(resname, chain, seq, atom_name, element, smiles_index) for the CSV."""
    return (
        placement.residue.resname,
        placement.residue.chain,
        placement.residue.seq,
        placement.atom_name,
        placement.element,
        placement.smiles_index,
    )


# ──────────────────────────────────────────────────────────────────────────────
# File discovery
# ──────────────────────────────────────────────────────────────────────────────

def discover_structures(path: Path, recursive: bool) -> list[Path]:
    """Collect every structure file under *path*."""
    if path.is_file():
        return [path]

    if not path.is_dir():
        raise FileNotFoundError(f"not a file or folder: {path}")

    candidates = path.rglob("*") if recursive else path.glob("*")
    return sorted(
        candidate
        for candidate in candidates
        if candidate.is_file()
        and candidate.name.lower().endswith(STRUCTURE_SUFFIXES)
    )


# ──────────────────────────────────────────────────────────────────────────────
# Reporting
# ──────────────────────────────────────────────────────────────────────────────

def print_atom_table(
    ligand: Chem.Mol, ligand_index: int | None,
    cofactor: Chem.Mol, cofactor_index: int | None,
    structures: list[Path],
    settings: "Settings",
) -> None:
    """Print atom numbering so the user can fill in STEP 5 of the data block."""
    samples: list[tuple[str, str, ResidueInstance, list[tuple[int, ...]] | None]] = []

    if structures:
        first = structures[0]
        try:
            residues = iter_residues(read_structure(first), settings.model)
            for mol, label in ((ligand, "ligand"), (cofactor, "cofactor")):
                found = None
                for residue in residues:
                    matches = match_residue(mol, residue, settings.bond_tolerance)
                    if matches:
                        found = (residue, matches)
                        break
                samples.append((label, first.name,
                                found[0] if found else None,
                                found[1] if found else None))
        except Exception as exc:  # pragma: no cover - diagnostic only
            log.warning("could not read %s for the atom table: %s", first.name, exc)

    for mol, index, label in (
        (ligand, ligand_index, "ligand"),
        (cofactor, cofactor_index, "cofactor"),
    ):
        print(f"\n{label.upper()}  ({mol.GetNumAtoms()} heavy atoms)")
        print(f"  {indexed_smiles(mol)}")
        print("  (the number in the SMILES above is the STEP 5 number + 1; "
              "RDKit reserves 0)")

        sample = next((s for s in samples if s[0] == label), None)
        if sample is None or sample[2] is None:
            print("  (atom names unavailable - set CIF_FOLDER in STEP 2 to see "
                  "them resolved inside a real structure)")
            print(f"  {'STEP 5':>6}  {'el':<3}  note")
            for atom in mol.GetAtoms():
                marker = "  <-- selected" if atom.GetIdx() == index else ""
                note = f"map:{atom.GetAtomMapNum()}" if atom.GetAtomMapNum() else ""
                print(f"  {atom.GetIdx():>6}  {atom.GetSymbol():<3}  {note}{marker}")
            continue

        _, source, residue, matches = sample
        print(f"  resolved against {source} → residue "
              f"{residue.resname} chain {residue.chain} seq {residue.seq}")
        print(f"  {'STEP 5':>6}  {'el':<3}  {'atom':<6}  note")
        for atom in mol.GetAtoms():
            position = matches[0][atom.GetIdx()]
            marker = "  <-- selected" if atom.GetIdx() == index else ""
            note = f"map:{atom.GetAtomMapNum()}" if atom.GetAtomMapNum() else ""
            print(f"  {atom.GetIdx():>6}  {atom.GetSymbol():<3}  "
                  f"{residue.atom_names[position]:<6}  {note}{marker}")


def read_structure(path: Path) -> gemmi.Structure:
    """Read a mmCIF/PDB file with gemmi."""
    return gemmi.read_structure(str(path))


# ──────────────────────────────────────────────────────────────────────────────
# Settings: read the data block, then let the command line override it
# ──────────────────────────────────────────────────────────────────────────────

MODE_ALIASES = {
    "show_atoms": "show_atoms", "show": "show_atoms", "atoms": "show_atoms",
    "print": "show_atoms", "print_atoms": "show_atoms",
    "measure": "measure", "run": "measure", "distance": "measure",
}

# Every field the command line is allowed to override
OVERRIDABLE = (
    "cif_dir", "ligand_smiles", "cofactor_smiles",
    "ligand_atom", "cofactor_atom",
    "ligand_atom_map", "cofactor_atom_map",
    "ligand_atom_name", "cofactor_atom_name",
    "ligand_resname", "cofactor_resname",
    "all_pairs", "recursive", "model", "bond_tolerance",
)


@dataclass
class Settings:
    """Everything one run needs, after the data block and CLI are merged."""
    mode:               str
    cif_dir:            Path | None
    ligand_smiles:      str
    cofactor_smiles:    str
    ligand_atom:        int | None
    cofactor_atom:      int | None
    ligand_atom_map:    int | None
    cofactor_atom_map:  int | None
    ligand_atom_name:   str | None
    cofactor_atom_name: str | None
    ligand_resname:     str | None
    cofactor_resname:   str | None
    out_dir:            Path
    out_filename:       str
    output_dir_given:   bool     # False while OUTPUT_DIR is still blank in the block
    all_pairs:          bool
    recursive:          bool
    model:              int
    bond_tolerance:     float

    @property
    def out_path(self) -> Path:
        return self.out_dir / self.out_filename


def _clean_path(value: object) -> Path | None:
    """Turn a data-block string into a Path, treating blank as unset."""
    text = str(value).strip() if value is not None else ""
    return Path(text) if text else None


def settings_from_data_block() -> Settings:
    """Read the INSERT YOUR DATA HERE block at the top of this file."""
    mode_key = str(MODE).strip().lower()
    if mode_key not in MODE_ALIASES:
        sys.exit(
            f"error: MODE is {MODE!r}. Use \"show_atoms\" (to list the atom "
            f"numbering) or \"measure\" (to write the CSV)."
        )

    return Settings(
        mode=MODE_ALIASES[mode_key],
        cif_dir=_clean_path(CIF_FOLDER),
        ligand_smiles=str(LIGAND_SMILES or "").strip(),
        cofactor_smiles=str(COFACTOR_SMILES or "").strip(),
        ligand_atom=LIGAND_ATOM,
        cofactor_atom=COFACTOR_ATOM,
        ligand_atom_map=None,
        cofactor_atom_map=None,
        ligand_atom_name=None,
        cofactor_atom_name=None,
        ligand_resname=None,
        cofactor_resname=None,
        out_dir=_clean_path(OUTPUT_DIR) or Path.cwd(),
        out_filename=str(OUTPUT_FILENAME or "atom_distances.csv").strip(),
        output_dir_given=_clean_path(OUTPUT_DIR) is not None,
        all_pairs=bool(ALL_PAIRS),
        recursive=bool(RECURSIVE),
        model=int(MODEL_INDEX),
        bond_tolerance=float(BOND_TOLERANCE),
    )


def apply_cli(settings: Settings, argv: list[str] | None) -> Settings:
    """Let command-line flags override the data block."""
    p = argparse.ArgumentParser(
        description=(
            "Override any field of the 'INSERT YOUR DATA HERE' block. "
            "With no flags the script runs entirely from that block."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--cif-dir", type=Path, default=None,
                   help="Folder of docked structures, or a single structure file")
    p.add_argument("--ligand-smiles", default=None,
                   help="SMILES of the docked ligand, as used for docking")
    p.add_argument("--cofactor-smiles", default=None,
                   help="SMILES of the docked cofactor, as used for docking")
    p.add_argument("--ligand-atom", type=int, default=None,
                   help="0-based heavy-atom index inside the ligand SMILES")
    p.add_argument("--cofactor-atom", type=int, default=None,
                   help="0-based heavy-atom index inside the cofactor SMILES")
    p.add_argument("--ligand-atom-map", type=int, default=None,
                   help="Atom-map number embedded in the ligand SMILES, e.g. [N:1]")
    p.add_argument("--cofactor-atom-map", type=int, default=None,
                   help="Atom-map number embedded in the cofactor SMILES")
    p.add_argument("--ligand-atom-name", default=None,
                   help="Ligand atom name already in the file (skips matching)")
    p.add_argument("--cofactor-atom-name", default=None,
                   help="Cofactor atom name already in the file (skips matching)")
    p.add_argument("--ligand-resname", default=None,
                   help="Restrict the ligand search to this residue name (LIG1)")
    p.add_argument("--cofactor-resname", default=None,
                   help="Restrict the cofactor search to this residue name (LIG2)")
    p.add_argument("--out-dir", type=Path, default=None,
                   help="Folder the CSV is written to")
    p.add_argument("--out", type=Path, default=None,
                   help="Full CSV path, e.g. /tmp/results.csv (overrides --out-dir)")
    p.add_argument("--out-filename", default=None, help="CSV file name")
    p.add_argument("--all-pairs", action="store_true", default=None,
                   help="One row per ligand×cofactor copy, not just the closest")
    p.add_argument("--recursive", action="store_true", default=None,
                   help="Search the structure folder recursively")
    p.add_argument("--model", type=int, default=None,
                   help="Model index to read from multi-model files")
    p.add_argument("--bond-tolerance", type=float, default=None,
                   help="Covalent-radius slack (Å) when rebuilding bonds")
    p.add_argument("--print-atoms", action="store_true", default=None,
                   help="Print the atom numbering and exit (MODE = show_atoms)")

    overrides = p.parse_args(argv)

    applied = False
    for field in OVERRIDABLE:
        value = getattr(overrides, field)
        if value is not None:
            setattr(settings, field, value)
            applied = True

    if overrides.out_dir is not None:
        settings.out_dir = overrides.out_dir
        settings.output_dir_given = True
        applied = True
    if overrides.out_filename is not None:
        settings.out_filename = overrides.out_filename
        applied = True
    if overrides.out is not None:
        settings.out_dir = overrides.out.parent or Path(".")
        settings.out_filename = overrides.out.name
        settings.output_dir_given = True
        applied = True

    # Passing flags is an explicit request to do the work, so it implies
    # MODE = "measure" — otherwise a stale MODE in the data block would keep
    # the script printing tables. --print-atoms still wins if it is given.
    if overrides.print_atoms:
        settings.mode = "show_atoms"
    elif applied:
        settings.mode = "measure"

    return settings


def missing_data(settings: Settings) -> list[str]:
    """Return one clear line per field the user still has to fill in."""
    problems: list[str] = []

    if not settings.ligand_smiles:
        problems.append('LIGAND_SMILES is empty    → fill in STEP 3')
    if not settings.cofactor_smiles:
        problems.append('COFACTOR_SMILES is empty  → fill in STEP 4')

    if settings.mode != "measure":
        return problems

    if not settings.cif_dir:
        problems.append("CIF_FOLDER is empty       → fill in STEP 2")
    if (settings.ligand_atom is None
            and settings.ligand_atom_name is None
            and settings.ligand_atom_map is None):
        problems.append(
            "LIGAND_ATOM is still None  → fill in STEP 5 with a number from the "
            'table printed by MODE = "show_atoms"'
        )
    if (settings.cofactor_atom is None
            and settings.cofactor_atom_name is None
            and settings.cofactor_atom_map is None):
        problems.append(
            "COFACTOR_ATOM is still None → fill in STEP 5 with a number from the "
            'table printed by MODE = "show_atoms"'
        )
    if not settings.output_dir_given:
        problems.append("OUTPUT_DIR is empty       → fill in STEP 6")

    return problems


# ──────────────────────────────────────────────────────────────────────────────
# CLI
# ──────────────────────────────────────────────────────────────────────────────

def main(argv: list[str] | None = None) -> int:
    settings = apply_cli(settings_from_data_block(), argv)

    problems = missing_data(settings)
    if problems:
        print("\n  The INSERT YOUR DATA HERE block is not finished yet:\n")
        for problem in problems:
            print(f"    • {problem}")
        print(f"\n  Edit those lines near the top of "
              f"{Path(__file__).name} and run it again.\n")
        return 2

    ligand = parse_smiles(settings.ligand_smiles, "ligand")
    cofactor = parse_smiles(settings.cofactor_smiles, "cofactor")

    # In "show_atoms" mode the STEP 5 numbers are still None on purpose: that
    # is the whole point, the user is there to read them off the table.
    optional = settings.mode == "show_atoms"

    ligand_index = select_smiles_atom(
        ligand, "ligand", index=settings.ligand_atom,
        map_num=settings.ligand_atom_map, atom_name=settings.ligand_atom_name,
        optional=optional,
    )
    cofactor_index = select_smiles_atom(
        cofactor, "cofactor", index=settings.cofactor_atom,
        map_num=settings.cofactor_atom_map, atom_name=settings.cofactor_atom_name,
        optional=optional,
    )

    structures: list[Path] = []
    if settings.cif_dir:
        try:
            structures = discover_structures(settings.cif_dir, settings.recursive)
        except (FileNotFoundError, OSError) as exc:
            sys.exit(f"error: {exc}")

    if settings.mode == "show_atoms":
        if not structures:
            log.info("no CIF_FOLDER set yet - showing SMILES numbering only")
        else:
            log.info("resolved atom names against %s", structures[0].name)
        print_atom_table(ligand, ligand_index, cofactor, cofactor_index,
                         structures, settings)
        print("\n  Copy the two numbers you need into STEP 5 "
              "(LIGAND_ATOM / COFACTOR_ATOM),\n  then set MODE = \"measure\" "
              "and run the script again.\n")
        return 0

    if not structures:
        sys.exit(
            f"error: no .cif/.pdb files found in {settings.cif_dir} "
            f"(set RECURSIVE = True if they are in sub-folders)"
        )

    for mol, index, atom_name, label in (
        (ligand, ligand_index, settings.ligand_atom_name, "ligand"),
        (cofactor, cofactor_index, settings.cofactor_atom_name, "cofactor"),
    ):
        if atom_name is not None and index is None:
            log.info("%-8s : atom %s of %d heavy atoms",
                     label, atom_name, mol.GetNumAtoms())
        else:
            log.info("%-8s : SMILES atom %d (%s) of %d heavy atoms",
                     label, index, mol.GetAtomWithIdx(index).GetSymbol(),
                     mol.GetNumAtoms())

    log.info("reading %d structure(s) from %s", len(structures), settings.cif_dir)

    rows: list[dict] = []
    n_ok = 0

    for path in structures:
        row: dict = {column: "" for column in CSV_COLUMNS}
        row["structure"] = path.name
        row["ligand_smiles_index"] = ligand_index
        row["cofactor_smiles_index"] = cofactor_index

        try:
            residues = iter_residues(read_structure(path), settings.model)

            ligand_copies = locate_molecule(
                residues, ligand, ligand_index, "ligand",
                tolerance=settings.bond_tolerance,
                resname=settings.ligand_resname,
                atom_name=settings.ligand_atom_name,
            )
            cofactor_copies = locate_molecule(
                residues, cofactor, cofactor_index, "cofactor",
                tolerance=settings.bond_tolerance,
                resname=settings.cofactor_resname,
                atom_name=settings.cofactor_atom_name,
            )

            for lig, cof, distance in measure(
                ligand_copies, cofactor_copies, all_pairs=settings.all_pairs
            ):
                entry = dict(row)
                entry.update(dict(zip(
                    ("ligand_resname", "ligand_chain", "ligand_seq", "ligand_atom",
                     "ligand_element", "ligand_smiles_index"), describe(lig)
                )))
                entry.update(dict(zip(
                    ("cofactor_resname", "cofactor_chain", "cofactor_seq",
                     "cofactor_atom", "cofactor_element", "cofactor_smiles_index"),
                    describe(cof)
                )))
                entry["distance_A"] = round(distance, 3)

                if lig.ambiguous or cof.ambiguous:
                    entry["status"] = "ok (ambiguous mapping - verify atom names)"
                else:
                    entry["status"] = "ok"

                rows.append(entry)
                n_ok += 1

            if len(ligand_copies) > 1 or len(cofactor_copies) > 1:
                log.warning("%s: %d ligand / %d cofactor copies found; "
                            "closest pair reported (ALL_PAIRS = True for all)",
                            path.name, len(ligand_copies), len(cofactor_copies))

        except Exception as exc:
            row["status"] = f"error: {exc}"
            rows.append(row)
            log.error("%s: %s", path.name, exc)

    out_path = settings.out_path
    out_path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows, columns=CSV_COLUMNS).to_csv(out_path, index=False)

    log.info("wrote %d/%d row(s) → %s", n_ok, len(rows), out_path.resolve())
    return 0 if n_ok else 1


if __name__ == "__main__":
    try:
        sys.exit(main())
    except ValueError as exc:
        sys.exit(f"error: {exc}")
