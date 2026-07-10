import os
import subprocess
from Bio import SeqIO
import shutil
import pandas as pd

from pietree import PieTree

from data import PROTEINS_PATH, PAPER_PATH, proteins

PATH_TO_IQTREE = "/home/pedro/Desktop/Programas/iqtree-3.0.1-Linux/bin/iqtree3"
PTNS = ["stricto", "t5h"]


def clean_alignment(input_fasta: str, output_fasta: str) -> int:
    records = list(SeqIO.parse(input_fasta, "fasta"))
    cleaned = [r for r in records if r.id and r.id.strip()]
    removed = len(records) - len(cleaned)
    if removed:
        print(f"  Removed {removed} sequence(s) with empty names")
    SeqIO.write(cleaned, output_fasta, "fasta")
    return len(cleaned)


def generate_phylogenetic_tree(parsed_fasta: str, output_prefix: str, outgroup_names: list[str]) -> None:
    if shutil.which(PATH_TO_IQTREE) is None:
        raise FileNotFoundError(f"IQ-TREE not found at '{PATH_TO_IQTREE}'.")

    os.makedirs(os.path.dirname(output_prefix), exist_ok=True)

    command = [
        PATH_TO_IQTREE, "-s", parsed_fasta, "-m", "MFP",
        "-alrt", "1000", "-B", "1000",
        "--redo", "--prefix", output_prefix,
    ]
    if outgroup_names:
        command += ["-o", ",".join(outgroup_names)]

    result = subprocess.run(command, capture_output=True, text=True)

    tail = result.stdout[-1500:] if len(result.stdout) > 1500 else result.stdout
    print(tail)

    if result.returncode != 0:
        raise RuntimeError(f"IQ-TREE failed with code {result.returncode}\n{result.stderr}")

    print("Phylogenetic tree generated successfully.")


def draw_tree(treefile: str, metadata_file: str, output_svg: str) -> None:
    tree = PieTree.from_newick(path=treefile, support_format="{bootstrap}/{alrt}")
    tree.annotate(pd.read_csv(metadata_file, sep=','), on="hit_id")
    tree.tips.rename("{hit_description}")
    tree.metadata("chr").highlight(allow_single_tip=True, label_position="center_right")
    tree.savefig(output_svg, size=(1800, 1200))


if __name__ == "__main__":
    for ptn in PTNS:
        ptn_data = proteins[ptn]

        for ooi in ptn_data["oois"]:
            aligned_file = f"{PROTEINS_PATH}/{ptn}/analysis/genomic/blast/{ooi}/all_seqs_aligned.fasta"

            if not os.path.isfile(aligned_file):
                print(f"  Skipping {ptn}/{ooi}: aligned file not found")
                continue

            print(f"\nProcessing {ptn} / {ooi}...")

            ref_orgs = list(ptn_data["organisms"].keys())
            outgroups = [
                f"{org}_{ptn}_{ptn_data['organisms'][org]['uniprot_id']}"
                for org in ref_orgs
            ]

            phylo_dir = f"{PROTEINS_PATH}/{ptn}/analysis/phylo/{ooi}"
            os.makedirs(phylo_dir, exist_ok=True)

            clean_fasta = f"{phylo_dir}/clean_alignment.fasta"
            num_seqs = clean_alignment(aligned_file, clean_fasta)

            if num_seqs < 4:
                print(f"  Skipping {ptn}/{ooi}: only {num_seqs} sequences (need at least 4)")
                continue

            output_prefix = f"{phylo_dir}/tree/{ptn}"
            generate_phylogenetic_tree(clean_fasta, output_prefix, outgroups)

            metadata_file = f"{PAPER_PATH}/tables/{ooi}_blast_abundance.csv"
            if not os.path.exists(metadata_file):
                print(f"  Warning: metadata not found at {metadata_file}, skipping tree drawing")
                continue

            draw_tree(f"{output_prefix}.treefile", metadata_file, f"{output_prefix}.svg")
            print(f"  Done! Tree: {output_prefix}.svg")
