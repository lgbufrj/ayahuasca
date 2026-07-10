from Bio.Blast import NCBIXML
from Bio import SeqIO
import pandas as pd
import re
import csv
import sys, os
from data import (
    proteins,
    organisms,
    PROTEINS_PATH,
    GENOME_PATH,
    TRANSCRIPTOME_PATH,
)
from helper_functions import format_fasta_description

EXPECT_CUTOFF = 1e-20
COVERAGE_CUTOFF = 50
MAX_LENGTH_MULTIPLIER = 2
OVERWRITE = "--overwrite" in sys.argv

# ── Abundance filter configuration ─────────────────────────────────────────
ABUNDANCE_FILTER = "--abundance-filter" in sys.argv
MIN_FPKM = 1.0
MIN_TPM = 1.0
MIN_COVERAGE = 0.0
GENOME_TYPE = "non_phased"

for i, arg in enumerate(sys.argv):
    if arg == "--min-fpkm" and i + 1 < len(sys.argv):
        MIN_FPKM = float(sys.argv[i + 1])
    elif arg == "--min-tpm" and i + 1 < len(sys.argv):
        MIN_TPM = float(sys.argv[i + 1])
    elif arg == "--min-coverage" and i + 1 < len(sys.argv):
        MIN_COVERAGE = float(sys.argv[i + 1])
    elif arg in ("--genome-type", "-g") and i + 1 < len(sys.argv):
        GENOME_TYPE = sys.argv[i + 1]

TRANSCRIPTOME_FILENAMES: dict[str, str] = {
    "tucunaca": "tucunaca.tab",
    "caupuri": "caupuri.tab",
}


def load_abundances(organism: str) -> dict[str, dict[str, float]]:
    filename = TRANSCRIPTOME_FILENAMES.get(organism)
    if not filename:
        return {}
    ab_path = f"{TRANSCRIPTOME_PATH}/{organism}/{filename}"
    if not os.path.exists(ab_path):
        return {}
    abundances: dict[str, dict[str, float]] = {}
    with open(ab_path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            gene_id = row["Gene ID"].strip()
            try:
                abundances[gene_id] = {
                    "coverage": float(row.get("Coverage", 0)),
                    "fpkm": float(row.get("FPKM", 0)),
                    "tpm": float(row.get("TPM", 0)),
                }
            except (ValueError, KeyError):
                pass
    return abundances


def load_mapping(map_csv: str) -> dict[str, str]:
    mapping: dict[str, str] = {}
    if not os.path.exists(map_csv):
        return mapping
    with open(map_csv) as f:
        reader = csv.DictReader(f)
        for row in reader:
            mapping[row["phased_id"]] = row["non_phased_id"]
    return mapping

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

        # ── Load abundance data for this organism ──────────────────────────
        abundances: dict[str, dict[str, float]] = {}
        mapping: dict[str, str] = {}
        if ABUNDANCE_FILTER:
            abundances = load_abundances(ooi)
            if GENOME_TYPE == "phased":
                map_csv = (
                    f"{PROTEINS_PATH}/{prot_name}/analysis/genomic/blast"
                    f"/{ooi}/{ref_organism}/gene_id_map.csv"
                )
                mapping = load_mapping(map_csv)
                if not mapping:
                    print(
                        f"  Warning: no gene_id_map.csv for phased genome, "
                        f"abundance filter will be incomplete for {prot_name} / {ooi}"
                    )

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

                            # ── Abundance filter ────────────────────────────────
                            if ABUNDANCE_FILTER and abundances and id_:
                                if GENOME_TYPE == "phased":
                                    non_phased_id = mapping.get(id_, "")
                                    lookup_id = non_phased_id.rsplit(".", 1)[0] if non_phased_id else ""
                                else:
                                    lookup_id = id_.rsplit(".", 1)[0]

                                ab = abundances.get(lookup_id) if lookup_id else None
                                if ab is not None:
                                    if (
                                        ab["fpkm"] < MIN_FPKM
                                        or ab["tpm"] < MIN_TPM
                                        or ab["coverage"] < MIN_COVERAGE
                                    ):
                                        continue

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