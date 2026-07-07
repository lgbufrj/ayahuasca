from Bio import SeqIO
from pandas import read_csv
from collections import defaultdict
import os
from data import proteins, organisms, PROTEINS_PATH, GENOME_PATH
from helper_functions import format_fasta_description

MAX_RESULTS = 10

def get_protein_organisms(proteins):
    result = []
    for prot_name, info in proteins.items():
        for species, org_info in info["organisms"].items():
            uniprot = org_info["uniprot_id"]
            fasta_name = f"{prot_name}_{uniprot}"
            result.append((prot_name, species, fasta_name))
    return result

def extract_ids_from_csv(list_file):
    df = read_csv(list_file, sep=';')
    return df['hit_id'].to_list()

print("Extracting ids...")

prot = get_protein_organisms(proteins)
combined = defaultdict(lambda: {"seen_ids": set(), "records": []})

for prot_name, species, fasta_name in prot:
    for ooi in [org if organisms[org]["ref"]==False else None for org in organisms]:
        if ooi is None: continue

        blast_res_path = f"{PROTEINS_PATH}/{prot_name}/analysis/genomic/blast/{ooi}/{species}/{fasta_name}_blast_filtered.fasta"

        if os.path.exists(blast_res_path):
            print(f"Skipping {ooi} {species} {fasta_name}...")
            blast_recs = list(SeqIO.parse(blast_res_path, "fasta"))
        else:
            input_csv_path = f"{PROTEINS_PATH}/{prot_name}/analysis/genomic/blast/{ooi}/{species}/{fasta_name}_blast_filtered.csv"

            if not os.path.exists(input_csv_path):
                print(f"Skipping {ooi} {species} {fasta_name}...")
                continue

            print(f"Extracting seqs for {prot_name} in {ooi} x {species}...")
            ids = extract_ids_from_csv(input_csv_path)[:MAX_RESULTS]

            seq_by_id = {
                rec.id.split("|")[0]: rec
                for rec in SeqIO.parse(f"{GENOME_PATH}/{ooi}/phased/prot.fasta", "fasta")
            }

            blast_recs = []
            for seq_id in ids:
                seq_record = seq_by_id.get(seq_id)
                if seq_record:
                    seq_record.id = format_fasta_description(seq_record.description, description_only=True)
                    seq_record.description = ""
                    blast_recs.append(seq_record)

            SeqIO.write(blast_recs, blast_res_path, "fasta")

        ref_ptn_path = f"{PROTEINS_PATH}/{prot_name}/reference/{species}/structure/{fasta_name}.fasta"
        merged_out_path = f"{PROTEINS_PATH}/{prot_name}/analysis/genomic/blast/{ooi}/{species}/reference_and_{fasta_name}_blast_filtered.fasta"

        ref_rec = next(SeqIO.parse(ref_ptn_path, "fasta"))
        ref_rec.id = f"{species}_{fasta_name}"
        ref_rec.description = ""
        SeqIO.write([ref_rec, *blast_recs], merged_out_path, "fasta")

        key = (prot_name, ooi)
        if ref_rec.id not in combined[key]["seen_ids"]:
            combined[key]["seen_ids"].add(ref_rec.id)
            combined[key]["records"].append(ref_rec)
        for rec in blast_recs:
            if rec.id not in combined[key]["seen_ids"]:
                combined[key]["seen_ids"].add(rec.id)
                combined[key]["records"].append(rec)

print("Writing combined files...")
for (prot_name, ooi), data in combined.items():
    out_dir = f"{PROTEINS_PATH}/{prot_name}/analysis/genomic/blast/{ooi}"
    combined_path = f"{out_dir}/all_seqs.fasta"
    SeqIO.write(data["records"], combined_path, "fasta")
    print(f"Combined FASTA for {prot_name} / {ooi}: {combined_path}")

print("Extraction completed!")
