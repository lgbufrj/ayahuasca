import subprocess
import data

# Function to get protein names and their corresponding fasta names, and combine them
def get_protein_organisms(proteins):
    result = []
    for prot_name, info in proteins.items():
        for species, org_info in info["organisms"].items():
            uniprot = org_info["uniprot_id"]
            fasta_name = f"{prot_name}_{uniprot}"
            result.append((prot_name, species, fasta_name))
    return result

def run_mafft(input_fasta, output_fasta):
    command = f"mafft --auto --preservecase {input_fasta} > {output_fasta}"
    try:
        subprocess.run(command, shell=True, check=True)
        print(f"MAFFT alignment completed: {output_fasta}")
    except subprocess.CalledProcessError as e:
        print(f"Error running MAFFT: {e}")

prot = get_protein_organisms(data.proteins)

# command to run BLASTP for each protein
for prot_name, species, fasta_name in prot:
    for organism in [org if data.organisms[org]["ref"]==False else None for org in data.organisms]:
        if organism is None: continue
        
        input = [f"{data.PROTEINS_PATH}/{prot_name}/analysis/genomic/blast/{organism}/{species}/{fasta_name}_blast_filtered.fasta", 
                f"{data.PROTEINS_PATH}/{prot_name}/analysis/genomic/blast/{organism}/{species}/reference_and_{fasta_name}_blast_filtered.fasta"
                ]
        output = [f"{data.PROTEINS_PATH}/{prot_name}/analysis/genomic/blast/{organism}/{species}/{fasta_name}_aligned.fasta",
                f"{data.PROTEINS_PATH}/{prot_name}/analysis/genomic/blast/{organism}/{species}/reference_and_{fasta_name}_aligned.fasta"
                ]

        for input_fasta, output_fasta in zip(input, output):
            run_mafft(input_fasta, output_fasta)

print("Alignment finished!")
