import string
import yaml
from rdkit import Chem
from Bio import SeqIO
from data import proteins, PROTEINS_PATH, compounds, COMPOUNDS_PATH
import os
import sys
import glob

# ── CLI flags ──────────────────────────────────────────────────────────
SELECTED_PROTEINS: list[str] | None = None
SELECTED_OOIS: list[str] | None = None
EXPORT_DIR: str | None = None
OVERWRITE = "--overwrite" in sys.argv

for i, arg in enumerate(sys.argv):
    if arg in ("--protein", "-p") and i + 1 < len(sys.argv):
        SELECTED_PROTEINS = sys.argv[i + 1].split(",")
    elif arg in ("--ooi", "-o") and i + 1 < len(sys.argv):
        SELECTED_OOIS = sys.argv[i + 1].split(",")
    elif arg in ("--export", "-e") and i + 1 < len(sys.argv):
        EXPORT_DIR = sys.argv[i + 1]


# ============================================================
# Utility functions
# ============================================================

def read_fasta_sequence(fasta_path: str) -> str:
    record = next(SeqIO.parse(fasta_path, "fasta"))
    return str(record.seq)


def sdf_to_smiles(sdf_path: str) -> str:
    mol = Chem.MolFromMolFile(sdf_path, removeHs=False)
    if mol is None:
        raise ValueError(f"Could not parse SDF: {sdf_path}")
    return Chem.MolToSmiles(mol)


def chain_generator():
    alphabet = string.ascii_uppercase
    for letter in alphabet:
        yield letter
    for first in alphabet:
        for second in alphabet:
            yield first + second


# ============================================================
# Entity builders
# ============================================================

def build_protein_entity(
    chain_ids,
    fasta_path: str = None,
    sequence: str = None,
    msa_path: str = None,
    modifications: list = None,
    cyclic: bool = False,
):
    if sequence is None and fasta_path is not None:
        sequence = read_fasta_sequence(fasta_path)
    elif sequence is None:
        raise ValueError("Provide either fasta_path or sequence")

    protein_dict = {
        "id": chain_ids,
        "sequence": sequence,
    }
    if msa_path is not None:
        protein_dict["msa"] = msa_path
    if modifications is not None:
        protein_dict["modifications"] = modifications
    if cyclic:
        protein_dict["cyclic"] = True

    return {"protein": protein_dict}


def build_ligand_entity(
    chain_id: str,
    sdf_path: str = None,
    smiles: str = None,
    ccd: str = None,
):
    provided = [
        sdf_path is not None,
        smiles is not None,
        ccd is not None,
    ]
    if sum(provided) != 1:
        raise ValueError("Provide exactly ONE of: sdf_path, smiles, or ccd")

    ligand_dict = {"id": chain_id}
    if sdf_path is not None:
        ligand_dict["smiles"] = sdf_to_smiles(sdf_path)
    elif smiles is not None:
        ligand_dict["smiles"] = smiles
    elif ccd is not None:
        ligand_dict["ccd"] = ccd

    return {"ligand": ligand_dict}


# ============================================================
# YAML writer
# ============================================================

def write_boltz_yaml(
    protein_sequence: str,
    ligands: list[dict],
    output_yaml: str,
) -> None:
    chain_ids = chain_generator()
    sequences = []

    prot_chain = next(chain_ids)
    sequences.append(build_protein_entity(
        sequence=protein_sequence,
        chain_ids=prot_chain,
        msa_path="empty",
    ))

    for ligand in ligands:
        lig_chain = next(chain_ids)
        sequences.append(build_ligand_entity(
            chain_id=lig_chain,
            sdf_path=ligand["sdf"],
        ))

    boltz_dict = {
        "version": 1,
        "sequences": sequences,
    }

    if ligands:
        first_ligand_chain = sequences[1]["ligand"]["id"]
        boltz_dict["properties"] = [
            {"affinity": {"binder": first_ligand_chain}}
        ]

    os.makedirs(os.path.dirname(output_yaml), exist_ok=True)
    with open(output_yaml, "w") as f:
        yaml.dump(boltz_dict, f, default_flow_style=False, sort_keys=False)
    print(f"  Wrote {output_yaml}")


# ============================================================
# Ligand builder
# ============================================================

def build_ligands(reaction: dict) -> list[dict]:
    return [
        {
            "sdf": f"{COMPOUNDS_PATH}/{cpd}/structure/{cpd}_{compounds[cpd]['pubchem_id']}.sdf"
        }
        for cpd in reaction["substrates"]
    ]


def resolve_output(path_template: str, *parts: str) -> str:
    """Return export path or project path depending on EXPORT_DIR."""
    if EXPORT_DIR:
        fname = "_".join(p for p in parts if p) + ".yaml"
        return os.path.join(EXPORT_DIR, fname)
    return path_template


# ============================================================
# Main
# ============================================================

if __name__ == "__main__":
    protein_items = list(proteins.items())

    if SELECTED_PROTEINS:
        protein_items = [
            (name, data) for name, data in protein_items
            if name in SELECTED_PROTEINS
        ]

    for protein_name, protein_data in protein_items:
        print(f"\nProcessing {protein_name}...")

        for reaction in protein_data["reactions"]:
            ligands = build_ligands(reaction)

            # ── Reference organisms ───────────────────────────────────
            for ref_org, org_data in protein_data["organisms"].items():
                fasta_path = (
                    f"{PROTEINS_PATH}/{protein_name}/reference"
                    f"/{ref_org}/structure/{protein_name}_{org_data['uniprot_id']}.fasta"
                )
                if not os.path.isfile(fasta_path):
                    print(f"  Reference FASTA not found: {fasta_path}")
                    continue

                seq = read_fasta_sequence(fasta_path)
                proj_path = (
                    f"{PROTEINS_PATH}/{protein_name}/analysis/structural/boltz"
                    f"/{ref_org}/inputs/{protein_name}_{org_data['uniprot_id']}_{reaction['id']}.yaml"
                )
                out = resolve_output(
                    proj_path,
                    ref_org, protein_name, org_data["uniprot_id"], reaction["id"],
                )

                if os.path.exists(out) and not OVERWRITE:
                    print(f"  Skipping (exists): {out}")
                    continue

                write_boltz_yaml(seq, ligands, out)

            # ── Organisms of interest (all sequences from all_seqs.fasta) ──
            for ooi in protein_data["oois"]:
                if SELECTED_OOIS and ooi not in SELECTED_OOIS:
                    continue
                all_seqs_path = (
                    f"{PROTEINS_PATH}/{protein_name}/analysis/genomic/blast"
                    f"/{ooi}/all_seqs.fasta"
                )
                if not os.path.isfile(all_seqs_path):
                    print(f"  all_seqs.fasta not found for {ooi}: {all_seqs_path}")
                    continue

                records = [
                    rec for rec in SeqIO.parse(all_seqs_path, "fasta")
                    if rec.id and rec.id.strip()
                ]

                for rec in records:
                    proj_path = (
                        f"{PROTEINS_PATH}/{protein_name}/analysis/structural/boltz"
                        f"/{ooi}/inputs/{rec.id}_{reaction['id']}.yaml"
                    )
                    out = resolve_output(
                        proj_path,
                        ooi, protein_name, rec.id, reaction["id"],
                    )

                    if os.path.exists(out) and not OVERWRITE:
                        continue

                    write_boltz_yaml(str(rec.seq), ligands, out)

    print("\nDone!")
