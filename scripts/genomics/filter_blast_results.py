from Bio.Blast import NCBIXML
from Bio import SeqIO
import re
import sys, os
from data import proteins, organisms, PROTEINS_PATH, GENOME_PATH
from helper_functions import format_fasta_description

EXPECT_CUTOFF = 1e-5

# Function to get protein names and their corresponding fasta names, and combine them
def get_protein_organisms(proteins):
    result = []
    for prot_name, info in proteins.items():
        for species, org_info in info["organisms"].items():
            uniprot = org_info["uniprot_id"]
            fasta_name = f"{prot_name}_{uniprot}"
            result.append((prot_name, species, fasta_name))
    return result

prots = get_protein_organisms(proteins)

for prot_name, ref_organism, fasta_name in prots:
    
    # print(prot_name, ref_organism, fasta_name)

    for ooi in [org if organisms[org]["ref"]==False else None for org in organisms]:
        if ooi is None: continue

        output_file = f"{PROTEINS_PATH}/{prot_name}/analysis/genomic/blast/{ooi}/{ref_organism}/{fasta_name}_blast_filtered.csv"

        if os.path.exists(output_file):
            print(f"Output file for {ref_organism} {prot_name} against {ooi} already exists, skipping...")
            continue

        input_file = f"{PROTEINS_PATH}/{prot_name}/analysis/genomic/blast/{ooi}/{ref_organism}/{fasta_name}.xml"

        if not os.path.exists(input_file):
            print(f"Input file for {ref_organism} {prot_name} against {ooi} does not exist, skipping...")
            continue

        with open(input_file) as result:

            blast_records = NCBIXML.parse(result)  # Lê o arquivo XML, NCBIXML.parse()retorna um iterador. Em termos simples, um iterador permite percorrer a saída do BLAST, recuperando registros BLAST um por um para cada resultado de pesquisa do BLAST:

            for blast_record in blast_records:  # cada query no FASTA
                
                length_cutoff = (blast_record.query_length * 0.5)

                with open(output_file, "w") as out_handle:
                    
                    out_handle.write(f"{'enzyme'};{'hit_id'};{'expect'};{'hit_len'};{'bit_score'};{'identity'};{'similarity'};{'gaps'};{'td'};{'sp'};{'chr'};{'hit_description'}\n") #Cabeçalho do arquivo de saída
                    
                    for alignment in blast_record.alignments: # Acessa cada alinhamento dentro do registro BLAST
                        for hsp in alignment.hsps: # Acessa cada HSP (High-scoring Segment Pair) dentro do alinhamento e pode ter múltiplos HSPs
                            if hsp.expect < EXPECT_CUTOFF: # hsp.expect representa o valor E-value e define o corte para menores que 0.001
                                if alignment.length >= length_cutoff: # Define o corte para o comprimento do alinhamento
                                    
                                    new_description, id_, name, td, sp, chr_ = format_fasta_description(alignment.hit_def)
                                    
                                    enzyme = blast_record.query.split()[0]
                                    hit_id = id_
                                    expect = hsp.expect
                                    hit_len = alignment.length
                                    bit_score = hsp.bits
                                    identity = hsp.identities
                                    similarity = hsp.positives
                                    gaps = hsp.gaps
                                    hit_description = new_description
                                    
                                    out_handle.write(f"{enzyme};{hit_id};{expect};{hit_len};{bit_score};{identity};{similarity};{gaps};{td};{sp};{chr_};{hit_description}\n")


        print(f"CSV archives created for {ref_organism} {prot_name} against {ooi}!")