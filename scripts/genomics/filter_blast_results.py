from Bio.Blast import NCBIXML
from Bio import SeqIO
import pandas as pd
import re
import sys, os
from data import proteins, organisms, PROTEINS_PATH, GENOME_PATH
from helper_functions import format_fasta_description

EXPECT_CUTOFF = 1e-20
COVERAGE_CUTOFF = 50
MAX_LENGTH_MULTIPLIER = 2
OVERWRITE = "--overwrite" in sys.argv

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

        if os.path.exists(output_file) and not OVERWRITE:
            print(f"Output file for {ref_organism} {prot_name} against {ooi} already exists, skipping...")
            continue

        input_file = f"{PROTEINS_PATH}/{prot_name}/analysis/genomic/blast/{ooi}/{ref_organism}/{fasta_name}.xml"

        if not os.path.exists(input_file):
            print(f"Input file for {ref_organism} {prot_name} against {ooi} does not exist, skipping...")
            continue

        with open(input_file) as result:

            blast_records = NCBIXML.parse(result)  # Lê o arquivo XML, NCBIXML.parse()retorna um iterador. Em termos simples, um iterador permite percorrer a saída do BLAST, recuperando registros BLAST um por um para cada resultado de pesquisa do BLAST:

            for blast_record in blast_records:  # cada query no FASTA

                rows = []

                for alignment in blast_record.alignments: # Acessa cada alinhamento dentro do registro BLAST
                    for hsp in alignment.hsps: # Acessa cada HSP (High-scoring Segment Pair) dentro do alinhamento e pode ter múltiplos HSPs
                        if hsp.expect < EXPECT_CUTOFF: # hsp.expect representa o valor E-value e define o corte para menores que 0.001
                            coverage = (hsp.query_end - hsp.query_start + 1) / blast_record.query_length * 100
                            if coverage < COVERAGE_CUTOFF:
                                continue

                            if alignment.length > blast_record.query_length * MAX_LENGTH_MULTIPLIER:
                                continue

                            new_description, id_, name, td, sp, chr_ = format_fasta_description(alignment.hit_def)

                            rows.append({
                                'enzyme': blast_record.query.split()[0],
                                'hit_id': id_,
                                'expect': hsp.expect,
                                'hit_len': alignment.length,
                                'bit_score': hsp.bits,
                                'coverage': round(min(coverage, 100.0), 2),
                                'identity': round((hsp.identities / alignment.length) * 100, 2),
                                'similarity': round((hsp.positives / alignment.length) * 100, 2),
                                'gaps': hsp.gaps,
                                'td': td,
                                'sp': sp,
                                'chr': chr_,
                                'hit_description': new_description,
                            })

                df = pd.DataFrame(rows)
                df.to_csv(output_file, sep=';', index=False)


        print(f"CSV archives created for {ref_organism} {prot_name} against {ooi}!")