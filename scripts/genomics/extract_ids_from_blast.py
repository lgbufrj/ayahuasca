from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pandas import read_csv
import re
import sys, os
from data import proteins, organisms, PROTEINS_PATH, GENOME_PATH
from helper_functions import format_fasta_description

# Function to get protein names, species, and fasta names
def get_protein_organisms(proteins):
    result = []
    for prot_name, info in proteins.items():
        for species, org_info in info["organisms"].items():
            uniprot = org_info["uniprot_id"]
            fasta_name = f"{prot_name}_{uniprot}"
            result.append((prot_name, species, fasta_name))
    return result

# Function to merge two fasta files into one
def merge_fastas(file_1, file_2, output_file):
    
    with open(output_file, 'w') as outfile:
        for fasta_file in [file_1, file_2]:
            
            with open(fasta_file, 'r') as infile:
                for line in infile:
                    outfile.write(line)

#function to extract IDs from a CSV file
def extract_ids_from_csv(list_file):
    df = read_csv(list_file, sep=';')
    return df['hit_id'].to_list()

print("Extracting ids...")

prot = get_protein_organisms(proteins)

# Loop through each protein, species, and fasta name
for prot_name, species, fasta_name in prot:
    
    for ooi in [org if organisms[org]["ref"]==False else None for org in organisms]:
        if ooi is None: continue
        
        blast_res_path = f"{PROTEINS_PATH}/{prot_name}/analysis/genomic/blast/{ooi}/{species}/{fasta_name}_blast_filtered.fasta"
    
        if os.path.exists(blast_res_path):
            print(f"Skipping {ooi} {species} {fasta_name}...")
            continue
        
        input_csv_path = f"{PROTEINS_PATH}/{prot_name}/analysis/genomic/blast/{ooi}/{species}/{fasta_name}_blast_filtered.csv"

        if not os.path.exists(input_csv_path):
            print(f"Skipping {ooi} {species} {fasta_name}...")
            continue

        print(f"Extracting IDs for {prot_name} in {ooi}...")
        ids = extract_ids_from_csv(input_csv_path)

        records = []

        for seq_record in SeqIO.parse(f"{GENOME_PATH}/{ooi}/phased/prot.fasta", "fasta"):
            if seq_record.id.split("|")[0] in ids:
                seq_record.id = format_fasta_description(seq_record.description, description_only=True)
                seq_record.description = ""
                records.append(seq_record)
        
        SeqIO.write(records, blast_res_path, "fasta")

        ref_ptn_path = f"{PROTEINS_PATH}/{prot_name}/reference/{species}/structure/{fasta_name}.fasta"
        merged_out_path = f"{PROTEINS_PATH}/{prot_name}/analysis/genomic/blast/{ooi}/{species}/reference_and_{fasta_name}_blast_filtered.fasta"

        merge_fastas(ref_ptn_path, blast_res_path, merged_out_path)

print("Extraction completed!")