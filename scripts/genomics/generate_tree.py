import os
import subprocess
from Bio import SeqIO, AlignIO
import re
import shutil
import pandas as pd

from pietree import PieTree

from data import PROTEINS_PATH, proteins

PATH_TO_IQTREE = "/home/pedro/Desktop/Programas/iqtree-3.0.1-Linux/bin/iqtree3"

ptn = "stricto"
ptn_data = proteins[ptn]

ooi = "tucunaca"
ref_org = "tabaco"
uniprot_id = ptn_data["organisms"][ref_org]["uniprot_id"]
# Definir outgroup e raiz como a sequencia de referencia
outgroups = [f"{ref_org}_{ptn}_{uniprot_id}"]
root = outgroups[0]

metadata_file =      f"{PROTEINS_PATH}/{ptn}/analysis/genomic/blast/{ooi}/{ref_org}/{ptn}_{uniprot_id}_blast_filtered.csv"
aligned_file =      f"{PROTEINS_PATH}/{ptn}/analysis/genomic/blast/{ooi}/all_seqs_aligned.fasta"
# aligned_file =      f"{PROTEINS_PATH}/{ptn}/analysis/genomic/blast/{ooi}/{ref_org}/reference_and_{ptn}_{uniprot_id}_aligned.fasta"
# output_file =       f"{PROTEINS_PATH}/{ptn}/analysis/phylo/{ooi}/{ref_org}/tree/{ptn}_{uniprot_id}"
output_file =       f"{PROTEINS_PATH}/{ptn}/analysis/phylo/{ooi}/tree/{ptn}"
# mb_output_file =    f"{PROTEINS_PATH}/{ptn}/analysis/phylo/{ooi}/{ref_org}/{ptn}_{uniprot_id}.nex"

# print(f"Outgroups: {outgroups}")

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


def draw_tree(treefile, output_file):
    tree = PieTree.from_newick(path=treefile, support_format="{bootstrap}/{alrt}")
    tree.annotate(pd.read_csv(metadata_file, sep=';'), on="hit_id")
    tree.tips.rename("{hit_description}")
    # tree.clade(tree.find_nodes(
    #     lambda n: n.metadata.get("expect", 1) == 0
    # )).highlight(label="Expect = 0", label_position="center_right")
    tree.metadata("chr").highlight(allow_single_tip=True, label_position="center_right")
    tree.savefig(output_file, size=(1800,1200))

if __name__ == "__main__":
    # Check if the input file exists
    if not os.path.isfile(aligned_file):
        raise FileNotFoundError(f"Input file '{aligned_file}' not found.")
    
    # Generate the phylogenetic tree
    # generate_phylogenetic_tree(aligned_file, output_file)
    draw_tree(f"{output_file}.treefile", f"{output_file}.svg")

    print("Done!")

    # setUpAndRunMrBayes(mb_output_file)

