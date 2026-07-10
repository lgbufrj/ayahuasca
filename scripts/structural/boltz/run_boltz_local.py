import os
import subprocess
from tqdm import tqdm

input_dir = "./boltz_inputs"
output_dir = "./boltz_results"

os.makedirs(output_dir, exist_ok=True)

yaml_files = [f for f in os.listdir(input_dir) if f.endswith(".yaml")]

for yaml_file in tqdm(yaml_files, desc="Running Boltz predictions"):
    input_path = os.path.join(input_dir, yaml_file)
    subprocess.run([
        "boltz", "predict",
        input_path,
        "--out_dir", output_dir,
        "--use_msa_server",
        "--recycling_steps", 10, # AlphaFold: 10
        "--diffusion_samples", 25 # AlphaFold: 25
    ])