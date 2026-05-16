import subprocess
import data
import os

# Function to get protein names and their corresponding fasta names, and combine them
def get_protein_organisms(proteins):
    result = []
    for prot_name, info in proteins.items():
        for species, org_info in info["organisms"].items():
            uniprot = org_info["uniprot_id"]
            rcsb = org_info.get("rcsb_id")
            fasta_name = f"{prot_name}_{uniprot}"
            rcsb_id = f"{rcsb}"
            result.append((prot_name, species, fasta_name, rcsb_id, uniprot))
    return result

prot = get_protein_organisms(data.proteins)

# command to run BLASTP for each protein
for prot_name, species, fasta_name, rcsb_id, uniprot in prot:
    
    if os.path.exists(f"./proteins/{prot_name}/reference/{species}/structure/{rcsb_id}.pdbqt"):
        continue
    input = [f"./proteins/{prot_name}/reference/{species}/structure/x_ray_{rcsb_id}.pdb"]
    
    output = [f"./proteins/{prot_name}/reference/{species}/structure/x_ray_{rcsb_id}.pdbqt"]

    if not os.path.exists(input[0]):    
        input = [f"./proteins/{prot_name}/reference/{species}/structure/{uniprot}.pdb"]
        
        output = [f"./proteins/{prot_name}/reference/{species}/structure/{uniprot}.pdbqt"]


    for input, output in zip(input, output):
        cmd = ["obabel", input, "-O", output] 
        
        with open(output, "w") as out:
            subprocess.run(cmd, stdout=out, check=True)

print("Convertion finished!")
