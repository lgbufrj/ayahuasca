import subprocess
import data

def run_mafft(input_fasta, output_fasta):
    command = f"mafft --auto --preservecase {input_fasta} > {output_fasta}"
    try:
        subprocess.run(command, shell=True, check=True)
        print(f"MAFFT alignment completed: {output_fasta}")
    except subprocess.CalledProcessError as e:
        print(f"Error running MAFFT: {e}")

for prot_name, info in data.proteins.items():
    ref_orgs = [(species, f"{prot_name}_{org_info['uniprot_id']}")
                for species, org_info in info["organisms"].items()]

    for ooi in [org for org in data.organisms if not data.organisms[org]["ref"]]:

        for species, fasta_name in ref_orgs:
            input_fastas = [
                f"{data.PROTEINS_PATH}/{prot_name}/analysis/genomic/blast/{ooi}/{species}/{fasta_name}_blast_filtered.fasta",
                f"{data.PROTEINS_PATH}/{prot_name}/analysis/genomic/blast/{ooi}/{species}/reference_and_{fasta_name}_blast_filtered.fasta",
            ]
            output_fastas = [
                f"{data.PROTEINS_PATH}/{prot_name}/analysis/genomic/blast/{ooi}/{species}/{fasta_name}_aligned.fasta",
                f"{data.PROTEINS_PATH}/{prot_name}/analysis/genomic/blast/{ooi}/{species}/reference_and_{fasta_name}_aligned.fasta",
            ]

            for input_fasta, output_fasta in zip(input_fastas, output_fastas):
                run_mafft(input_fasta, output_fasta)

        combined_fasta = f"{data.PROTEINS_PATH}/{prot_name}/analysis/genomic/blast/{ooi}/all_seqs.fasta"
        combined_aligned = f"{data.PROTEINS_PATH}/{prot_name}/analysis/genomic/blast/{ooi}/all_seqs_aligned.fasta"
        run_mafft(combined_fasta, combined_aligned)

print("Alignment finished!")
