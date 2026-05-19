import os, subprocess
from data import PROTEINS_PATH, proteins
from Bio.PDB import PDBParser
import numpy as np
import pandas as pd
import json

PRANK_PATH = "/home/pedro/Desktop/Programas/p2rank/p2rank_2.5.1/prank"

def predict_pockets(input_file, prediction_output_path):

    os.makedirs(os.path.dirname(prediction_output_path), exist_ok=True)

    cmd = [
        PRANK_PATH,
        "predict",
        "-f", input_file,
        "-o", prediction_output_path
    ]

    subprocess.run(cmd)

def generate_docking_grid(input_file, prediction_output_path, final_output_file):

    prediction_output_file = prediction_output_path + f"/{input_file.split('/')[-1]}_predictions.csv"

    prediction_df = pd.read_csv(prediction_output_file)
    pocket_residues = set()
    for v in prediction_df[" residue_ids"][0].split(" "):
        if v != "":
            pocket_residues.add(int(str(v).split("_")[1]))

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", input_file)

    coords = []

    for residue in structure[0]["A"]:
        if residue.id[1] in pocket_residues:
            for atom in residue:
                coords.append(atom.coord)

    coords = np.array(coords)

    mins = coords.min(axis=0)
    maxs = coords.max(axis=0)

    padding = 10.0

    sizes = (maxs - mins) + padding
    center = coords.mean(axis=0)

    docking_grid = {
        "center": center.tolist(),
        "size": sizes.tolist()
    }

    os.makedirs(os.path.dirname(final_output_file), exist_ok=True)

    with open(final_output_file, "w") as f:
        json.dump(docking_grid, f, indent=4)

if __name__ == "__main__":
    
    protein = "asmt"

    protein_data = proteins[protein]

    organism = "tucunaca"
    # organism_data = protein_data['organisms'][organism]

    # input_file = PROTEINS_PATH + f"/{protein}/reference/{organism}/structure/{protein}_{organism_data['uniprot_id']}.pdb"
    input_file = PROTEINS_PATH + f"/{protein}/analysis/structural/alphafold/{organism}/best_model.pdb"
        
    prediction_output_path = PROTEINS_PATH + f"/{protein}/analysis/structural/p2rank/{organism}/output"
    final_output_file = PROTEINS_PATH + f"/{protein}/analysis/structural/p2rank/{organism}/docking_grid.json"

    predict_pockets(input_file, prediction_output_path)
    generate_docking_grid(input_file, prediction_output_path, final_output_file)