import string
import yaml
from rdkit import Chem
from Bio import SeqIO
from data import proteins, PROTEINS_PATH, compounds, COMPOUNDS_PATH
import os
import pandas as pd

# ============================================================
# Utility functions
# ============================================================

def read_fasta_sequence(fasta_path: str) -> str:
    """
    Reads the first sequence from a FASTA file using Biopython.
    """

    record = next(SeqIO.parse(fasta_path, "fasta"))

    return str(record.seq)


def sdf_to_smiles(sdf_path: str) -> str:
    """
    Converts an SDF file to SMILES using RDKit.
    """

    mol = Chem.MolFromMolFile(sdf_path, removeHs=False)

    if mol is None:
        raise ValueError(f"Could not parse SDF: {sdf_path}")

    return Chem.MolToSmiles(mol)


def chain_generator():
    """
    Generates chain IDs:
    A, B, C ... Z, AA, AB, AC ...
    """

    alphabet = string.ascii_uppercase

    # Single-letter IDs
    for letter in alphabet:
        yield letter

    # Double-letter IDs
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
    """
    Builds a protein entity dictionary.
    """

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

    return {
        "protein": protein_dict
    }


def build_ligand_entity(
    chain_id: str,
    sdf_path: str = None,
    smiles: str = None,
    ccd: str = None,
):
    """
    Builds a ligand entity dictionary.
    """

    provided = [
        sdf_path is not None,
        smiles is not None,
        ccd is not None,
    ]

    if sum(provided) != 1:
        raise ValueError(
            "Provide exactly ONE of: sdf_path, smiles, or ccd"
        )

    ligand_dict = {
        "id": chain_id,
    }

    if sdf_path is not None:
        ligand_dict["smiles"] = sdf_to_smiles(sdf_path)

    elif smiles is not None:
        ligand_dict["smiles"] = smiles

    elif ccd is not None:
        ligand_dict["ccd"] = ccd

    return {
        "ligand": ligand_dict
    }


# ============================================================
# Main generator
# ============================================================

def generate_boltz_input(
    proteins,
    ligands=None,
    output_yaml="boltz_input.yaml",
):
    """
    Generates a Boltz YAML input.

    Parameters
    ----------
    proteins : list of dict

        Example:

        proteins = [
            {
                "fasta": "proteinA.fasta",
                "copies": 2
            },
            {
                "fasta": "proteinB.fasta",
                "copies": 1
            }
        ]

    ligands : list of dict

        Example:

        ligands = [
            {
                "sdf": "ligand.sdf"
            },
            {
                "ccd": "FAD"
            }
        ]
    """

    chain_ids = chain_generator()

    sequences = []

    # ========================================================
    # Proteins
    # ========================================================

    for protein in proteins:

        fasta_path = protein.get("fasta")
        
        if(fasta_path is not None):
            if not os.path.isfile(fasta_path):
                continue

        copies = protein.get("copies", 1)

        msa_path = protein.get("msa")

        modifications = protein.get("modifications")

        cyclic = protein.get("cyclic", False)

        protein_chain_ids = [
            next(chain_ids)
            for _ in range(copies)
        ]

        # Use scalar instead of list for monomers
        if len(protein_chain_ids) == 1:
            protein_chain_ids = protein_chain_ids[0]

        entity = build_protein_entity(
            fasta_path=protein.get("fasta"),
            sequence=protein.get("sequence"),
            chain_ids=protein_chain_ids,
            msa_path="empty",
        )

        sequences.append(entity)

    # ========================================================
    # Ligands / Cofactors
    # ========================================================

    if ligands is not None:

        for ligand in ligands:

            ligand_chain = next(chain_ids)

            entity = build_ligand_entity(
                chain_id=ligand_chain,
                sdf_path=ligand.get("sdf"),
                # smiles=ligand.get("smiles"),
                # ccd=ligand.get("ccd"),
            )

            sequences.append(entity)

    # ========================================================
    # Final YAML
    # ========================================================

    boltz_dict = {
        "version": 1,
        "sequences": sequences
    }
    
    # Add affinity property pointing to the first ligand's chain ID
    if ligands is not None and len(ligands) > 0:
        first_ligand_chain = sequences[len(proteins)]["ligand"]["id"]
        boltz_dict["properties"] = [
            {"affinity": {"binder": first_ligand_chain}}
        ]

    os.makedirs(os.path.dirname(output_yaml), exist_ok=True)

    with open(output_yaml, "w") as f:

        yaml.dump(
            boltz_dict,
            f,
            default_flow_style=False,
            sort_keys=False,
        )

    print(f"Boltz YAML written to: {output_yaml}")

def get_best_hit_fasta_sequence(blast_csv_path: str, blast_fasta_path: str) -> str:
    """
    Reads the best hit sequence ID from the first row of the BLAST CSV
    and extracts the matching sequence from the multifasta file.
    """
    if not os.path.isfile(blast_csv_path):
        return None
    df = pd.read_csv(blast_csv_path, sep=";")
    best_hit_id = str(df.iloc[0]["hit_id"])

    for record in SeqIO.parse(blast_fasta_path, "fasta"):
        if record.id == best_hit_id:
            return str(record.seq)

    # raise ValueError(f"Best hit ID '{best_hit_id}' not found in {blast_fasta_path}")

if __name__ == "__main__":
    
    for protein_name, protein_data in proteins.items():
        # if protein_name not in ["t5h"]: continue
        for reaction in protein_data["reactions"]:
            
            ligands = [
                {
                    "sdf": f"{COMPOUNDS_PATH}/{cpd}/structure/{cpd}_{compounds[cpd]['pubchem_id']}.sdf"
                } for cpd in reaction['substrates']
            ]

            for ref_organism, organism_data in protein_data['organisms'].items():
                
                input_protein = [{
                    "fasta": f"{PROTEINS_PATH}/{protein_name}/reference/{ref_organism}/structure/{protein_name}_{organism_data['uniprot_id']}.fasta",
                    "copies": 1
                }]
                

                output_file = f"{PROTEINS_PATH}/{protein_name}/analysis/structural/boltz/{ref_organism}/inputs/{protein_name}_{organism_data['uniprot_id']}_{reaction['id']}.yaml"
                # output_file = f"./ayahusca_boltz_inputs/{ref_organism}_{protein_name}_{organism_data['uniprot_id']}_{reaction['id']}.yaml"
                    
                generate_boltz_input(
                    proteins=input_protein,
                    ligands=ligands,
                    output_yaml=output_file
                )
                    
                for ooi in protein_data['oois']:
                    blast_csv  = f"{PROTEINS_PATH}/{protein_name}/analysis/genomic/blast/{ooi}/{ref_organism}/{protein_name}_{organism_data['uniprot_id']}_blast_filtered.csv"
                    blast_fasta = f"{PROTEINS_PATH}/{protein_name}/analysis/genomic/blast/{ooi}/{ref_organism}/{protein_name}_{organism_data['uniprot_id']}_blast_filtered.fasta"

                    best_hit_seq = get_best_hit_fasta_sequence(blast_csv, blast_fasta)
                    
                    if not best_hit_seq:
                        print(f"No valid best hit found for {protein_name} in {ref_organism} with OOI {ooi}. Skipping Boltz input generation.")
                        continue

                    input_protein = [{
                        "sequence": best_hit_seq,   # pass sequence directly
                        "copies": 1
                    }]
                    
                    output_file=f"{PROTEINS_PATH}/{protein_name}/analysis/structural/boltz/{ooi}/{ref_organism}/inputs/{protein_name}_{ooi}_{ref_organism}_{reaction['id']}.yaml"
                    # output_file = f"./ayahusca_boltz_inputs/{ooi}_{ref_organism}_{protein_name}_{organism_data['uniprot_id']}_{reaction['id']}.yaml"

                    generate_boltz_input(
                        proteins=input_protein,
                        ligands=ligands,
                        output_yaml=output_file
                    )