from Bio import SeqIO
import re
import subprocess

ENZYME = "stricto"

genes = []

for record in SeqIO.parse("./faseado/cds.fasta", "fasta"):
	desc = record.description

	length = len(record.seq)

	## Stricto (best-hit: TucunacaPG00000159720.1)
	if (length >= 350) and (("strictosidine" in desc) or ("STRICTOSIDINE" in desc)) and not ("hypothetical" in desc) and not ("uncharacterized" in desc):

	## SARD4 (best-hit: )
	# if ("protein SAR DEFICIENT 4" in desc):

	## T5H (best-hit: )
	# if ("tryptamine 5-hydroxylase-like" in desc):

	## ASMT (best-hit: TucunacaPG00000082609.1)
	# if ("Methyltransf 2 domain-containing protein" in desc) and ("class I-like SAM-binding methyltransferase" in desc):
	
	## Peroxidase (best-hit: TucunacaPG00000159720.1)
	# if ("Peroxidase" in desc) and ("Removal of H[2]O[2] oxidation of toxic reductants biosynthesis and degradation of lignin suberization auxin catabolism response to environmental stresses such as wounding pathogen attack and oxidative stress" in desc):

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
			record.id = f"{id_} | {name} | td {td} | sp {sp} | {chr_}".replace(" ", "_")
			record.description = ""

		genes += [record]

print(len(genes))

SeqIO.write(genes, f"./faseado/phylo/{ENZYME}/genes.fasta", "fasta")

# mafft --auto ./faseado/phylo/t5h/genes.fasta > ./faseado/phylo/t5h/genes_aligned.fasta
subprocess.run(
	["mafft", "--auto", f"./faseado/phylo/{ENZYME}/genes.fasta"],
	stdout=open(f"./faseado/phylo/{ENZYME}/genes_aligned.fasta", "w"),
	check=True
)

# /Applications/iqtree_2.3.6/bin/iqtree2 -s ~/local_desktop/artigos/ayahuasca/faseado/phylo/stricto/genes_aligned.fasta -m MFP -alrt 1000 -B 1000