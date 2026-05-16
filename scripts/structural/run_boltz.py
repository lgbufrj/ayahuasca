import os
import subprocess
from data import PROTEINS_PATH, proteins

# boltz predict ../pathways/harmine_biosynthesis/proteins/asmt/reference/arabidopsis/structure/asmt_Q9T003.yaml --out_dir ../pathways/harmine_biosynthesis/proteins/asmt/reference/arabidopsis/structure/ --use_msa_server


input_path = f"{PROTEINS_PATH}/asmt/reference/arabidopsis/structure/asmt_Q9T003.yaml"
output_path = f"{PROTEINS_PATH}/asmt/analysis/structural/boltz/arabidopsis/"

os.makedirs(output_path, exist_ok=True)

subprocess.run([
    "boltz", "predict",
    input_path,
    "--out_dir", output_path,
    "--use_msa_server",
    "--accelerator", "gpu"
])