from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from pandas import read_csv
from collections import defaultdict
import os
import sys
from data import proteins, organisms, PROTEINS_PATH, GENOME_PATH
from helper_functions import format_fasta_description, sequence_key

OVERWRITE = "--overwrite" in sys.argv
GENOME_TYPE = "non_phased"
for i, arg in enumerate(sys.argv):
    if arg in ("--genome-type", "-g") and i + 1 < len(sys.argv):
        GENOME_TYPE = sys.argv[i + 1]

MAX_RESULTS = 100

def get_protein_organisms(proteins):
    result = []
    for prot_name, info in proteins.items():
        for species, org_info in info["organisms"].items():
            uniprot = org_info["uniprot_id"]
            fasta_name = f"{prot_name}_{uniprot}"
            result.append((prot_name, species, fasta_name))
    return result

def extract_ids_from_csv(list_file):
    try:
        df = read_csv(list_file, sep=';')
        return df['hit_id'].to_list()
    except:
        return []

def read_cached_records(path):
    records = list(SeqIO.parse(path, "fasta"))
    if any(not rec.id for rec in records):
        return None
    return records

print("Extracting ids...")

prot = get_protein_organisms(proteins)
combined = defaultdict(lambda: {"seen_ids": set(), "records": []})

for prot_name, species, fasta_name in prot:
    for ooi in [org if organisms[org]["ref"]==False else None for org in organisms]:
        if ooi is None: continue

        blast_res_path = f"{PROTEINS_PATH}/{prot_name}/analysis/genomic/blast/{ooi}/{species}/{fasta_name}_blast_filtered.fasta"

        if os.path.exists(blast_res_path) and not OVERWRITE:
            cached_recs = read_cached_records(blast_res_path)
            if cached_recs is not None:
                print(f"Skipping {ooi} {species} {fasta_name}...")
                blast_recs = cached_recs
            else:
                print(f"Cached FASTA for {ooi} {species} {fasta_name} has records without ids, regenerating...")
        else:
            cached_recs = None

        if cached_recs is None:
            input_csv_path = f"{PROTEINS_PATH}/{prot_name}/analysis/genomic/blast/{ooi}/{species}/{fasta_name}_blast_filtered.csv"

            if not os.path.exists(input_csv_path):
                print(f"Skipping {ooi} {species} {fasta_name}...")
                continue

            genome_prot_path = f"{GENOME_PATH}/{ooi}/{GENOME_TYPE}/prot.fasta"
            if not os.path.exists(genome_prot_path):
                print(f"Skipping {ooi} {species} {fasta_name}: {genome_prot_path} not found")
                continue

            print(f"Extracting seqs for {prot_name} in {ooi} x {species}...")
            ids = [sequence_key(hit_id) for hit_id in extract_ids_from_csv(input_csv_path)[:MAX_RESULTS]]
            ids = [hit_id for hit_id in ids if hit_id]
            if ids == []: continue

            seq_by_id = {}
            for rec in SeqIO.parse(genome_prot_path, "fasta"):
                key = sequence_key(rec.id)
                if key:
                    seq_by_id.setdefault(key, rec)

            blast_recs = []
            seen = set()
            for seq_id in ids:
                seq_record = seq_by_id.get(seq_id)
                if seq_record is None or seq_id in seen:
                    continue
                seen.add(seq_id)
                _, new_id, _name, _td, _sp, _chr_ = format_fasta_description(seq_record.description)
                if not new_id:
                    print(f"  Warning: skipping {seq_id}, no usable id in '{seq_record.description}'")
                    continue
                blast_recs.append(SeqRecord(seq_record.seq, id=new_id, description=""))

            if not blast_recs:
                print(f"  Warning: no genome sequence matched {len(ids)} hit id(s) for {fasta_name}; regenerate the CSV with filter_blast_results.py --overwrite")

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
