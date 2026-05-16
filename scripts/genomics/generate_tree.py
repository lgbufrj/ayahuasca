import os
import subprocess
from Bio import SeqIO, AlignIO
import re
import shutil

PATH_TO_IQTREE = "/home/pedro/Desktop/Programas/iqtree-3.0.1-Linux/bin/iqtree3"

ptn_name = "sard4"
uniprot_id = "Q9FLY0"
reference_species = "arabidopsis"

aligned_file =      f"../proteins/{ptn_name}/analysis/genomic/blast/{reference_species}/reference_and_{ptn_name}_{uniprot_id}_aligned.fasta"
parsed_fasta =       f"../proteins/{ptn_name}/analysis/phylogenetic/{reference_species}/reference_and_{ptn_name}_{uniprot_id}_aligned_parsed.fasta"
output_file =       f"../proteins/{ptn_name}/analysis/phylogenetic/{reference_species}/tree/{ptn_name}_{uniprot_id}"
# mb_output_file =    f"../proteins/{ptn_name}/analysis/phylogenetic/{reference_species}/{ptn_name}_{uniprot_id}.nex"

# Definir outgroup e raiz como a sequencia de referencia
outgroups = ["sp|Q9FLY0|SARD4_ARATH"]
root = "sp|Q9FLY0|SARD4_ARATH"

# print(f"Outgroups: {outgroups}")

def generate_parsed_fasta():
    genes = []
    for record in SeqIO.parse(aligned_file, "fasta"):
        desc = record.description
        match = re.match(
                r"([^|]+)\|([^|]+).*?transmembrane domain:(\w+)\|signal peptide:(\w+).*?\bchr([0-9]+[a-zA-Z]?)\b",
                desc,
            )
        if match:
            id_ = match.group(1)
            name = match.group(2)
            td = "y" if match.group(3).lower() == "yes" else "n"
            sp = "y" if match.group(4).lower() == "yes" else "n"
            chr_ = match.group(5)
            record.id = f"{id_} | {ptn_name} | td {td} | sp {sp} | {chr_}".replace(" ", "_")
            record.description = ""   
        else:
            print(f"Warning: Description format not recognized for record {record.id}. Skipping.")
        genes.append(record)

    # Create directory if it doesn't exist
    os.makedirs(os.path.dirname(parsed_fasta), exist_ok=True)
    SeqIO.write(genes, parsed_fasta, "fasta")

def generate_phylogenetic_tree(parsed_fasta, output_file):
    if shutil.which(PATH_TO_IQTREE) is None:
        raise FileNotFoundError(f"IQ-TREE not found in PATH (tried '{PATH_TO_IQTREE}').")

    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    command = [
        PATH_TO_IQTREE, "-s", parsed_fasta, "-m", "MFP",
        "-alrt", "1000", "-B", "1000",
        "--redo", "--prefix", output_file
    ]

    result = subprocess.run(command, capture_output=True, text=True)
    print("STDOUT:", result.stdout)
    print("STDERR:", result.stderr)

    if result.returncode == 0:
        print("Phylogenetic tree generated successfully.")
    else:
        raise RuntimeError(f"IQ-TREE failed with code {result.returncode}")

# def setUpAndRunMrBayes(file):

#     # Convert aligned FASTA file to NEXUS format
#     AlignIO.convert(aligned_file, "fasta", file, "nexus", "DNA")

#     # Check if the conversion was successful
#     if not os.path.isfile(file):
#         raise FileNotFoundError(f"Failed to convert {aligned_file} to NEXUS format.")

#     # Set up the MrBayes configuration
#     with open(file, "a") as mb_file:
#         mb_file.write("Begin mrbayes;\n")
#         mb_file.write("  set autoclose=yes;\n")
#         mb_file.write("  lset nst=6 rates=gamma;\n")
#         mb_file.write("  mcmc ngen=1000000 printfreq=1000 samplefreq=100 nchains=4 savebrlens=yes;\n")
#         mb_file.write("  sumt burnin=250;\n")
#         mb_file.write("End;\n")

#     # Run MrBayes
#     subprocess.run(["mb", file])
#     # Clean up
#     # os.remove("mb.nex")
#     # os.remove("mb.nex.bak")
#     # os.remove("mb.nex.log")
#     # os.remove("mb.nex.tre")
#     # os.remove("mb.nex.p") 

if __name__ == "__main__":
    # Check if the input file exists
    if not os.path.isfile(aligned_file):
        raise FileNotFoundError(f"Input file '{aligned_file}' not found.")
    
    # Generate the parsed FASTA file
    generate_parsed_fasta()
    
    # Generate the phylogenetic tree
    generate_phylogenetic_tree(parsed_fasta, output_file)

    # Replace every '_' with ' ' in the output file
    with open(output_file + ".treefile", "r") as file:
        tree_data = file.read()

    tree_data = tree_data.replace('_', ' ')

    with open(output_file + ".treefile", "w") as file:
        file.write(tree_data)

    print("Done!")

    # setUpAndRunMrBayes(mb_output_file)

