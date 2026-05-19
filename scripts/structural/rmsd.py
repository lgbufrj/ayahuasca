import os, subprocess
from data import proteins, PROTEINS_PATH
import re
import pandas as pd

USALIGN_PATH = "/home/pedro/Desktop/Programas/USalign/USalign"

def usalign(pdb1, pdb2):
    """
    TM-score has values in (0,1] with 1 indicating an identical structure match, where a TM-score ≥0.5 (or 0.45) means the structures share the same global topology for proteins (or RNAs).
    """

    result = subprocess.run(
        [USALIGN_PATH, pdb1, pdb2],
        capture_output=True,
        text=True
    )

    text = result.stdout

    rmsd = re.search(r"RMSD=\s+([\d.]+)", text)
    tm   = re.search(r"TM-score=\s+([\d.]+)", text)

    return {
        "rmsd": float(rmsd.group(1)),
        "tm_score": float(tm.group(1))
    }

if __name__ == "__main__":
    
    for protein, protein_data in proteins.items():
        
        results = []
        
        for ooi in protein_data["oois"]:
        
            ooi_pdb = PROTEINS_PATH + f"/{protein}/analysis/structural/alphafold/{ooi}/best_model.pdb"
        
            if not os.path.exists(ooi_pdb):
                # print(f"Mobile PDB not found for {protein} ({ooi}): {ooi_pdb}")
                continue
        
            for organism, organism_data in protein_data['organisms'].items():
            
                reference_pdb = PROTEINS_PATH + f"/{protein}/reference/{organism}/structure/{protein}_{organism_data['uniprot_id']}.pdb"           
                
                if not os.path.exists(reference_pdb):
                    # print(f"Reference PDB not found for {protein} ({organism}): {reference_pdb}")
                    continue

                alignment = usalign(reference_pdb, ooi_pdb)

                # print(f"RMSD for {protein} ({ooi} x {organism}): {alignment['rmsd']}\nTM-score for {protein} ({ooi} x {organism}): {alignment['tm_score']}")
    
                results.append({
                    "protein": protein,
                    "ooi": ooi,
                    "reference_organism": organism,
                    "rmsd": alignment['rmsd'],
                    "tm_score": alignment['tm_score']
                })

        # Save results to a CSV
        output_path = PROTEINS_PATH + f"/{protein}/analysis/structural/alignment"
        os.makedirs(output_path, exist_ok=True)
        df = pd.DataFrame(results)
        df.to_csv(f"{output_path}/structural_alignments.csv", index=False)

        print(f"Results saved to {protein} | structural_alignments.csv")