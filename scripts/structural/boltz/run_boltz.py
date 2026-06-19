import os
import subprocess
from data import PROTEINS_PATH, proteins

# boltz predict {input_path} --out_dir {output_dir} --use_msa_server


input_path = f"{PROTEINS_PATH}/asmt/analysis/structural/boltz/tucunaca/tabaco/inputs/asmt_tucunaca_tabaco_rr.yaml"
output_dir = f"{PROTEINS_PATH}/asmt/analysis/structural/boltz/tucunaca/tabaco/"

os.makedirs(output_dir, exist_ok=True)

# subprocess.run([
#     "boltz", "predict",
#     input_path,
#     "--out_dir", output_dir,
#     "--use_msa_server",
#     "--accelerator", "gpu"
# ])

from boltz_api import Boltz

client = Boltz(api_key="sk_bc_ws_live_MbYfb_DnYObGSLLniFYHV4DRQ2BoffmFyLY2haV0RSA")

run_dir = client.experiments.run_structure_and_binding(
    entities=[
        {
        "type": "protein", 
        "value": "MNEIQTNSNLINRDDDEEEAQADIEIWDYVFGFVKMAVVKCAIELGISEAIENHGGPISLSELAASLNCDPSGLHRIMRFLIHYRFFKETVDGYVHTALSRRLLLKVPNSMADIILMESSHVMLEPWHQLSSYLLDSKKPPFERAHGIDLWKFCSLNPSYSKLIDDAMACDARLAVKAVIQGCPEIFKGIGTMVDVGGGNGTALNMFVKAFPWIQGINFDLPHVVEVAPKLDGVKHVGGDMFQSVPKADAAYLMKVLHDWSDDESIQILRRCREAIEESKGKVIIVESVLEKDEDCDRLEFVRLMLDMVMLAHTSKGKERTLKQWDYVLHQAGFSSYDIKAIDTYHSIIIAVP", 
        "chain_ids": ["A"]
        },
        {
            "type": "ligand_smiles",
            "value": "[H]Oc1c([H])c([H])c2c(c1[H])c(C([H])([H])C([H])([H])N([H])C(=O)C([H])([H])[H])c([H])n2[H]",
            "chain_ids": ["B"],
        },
        {
            "type": "ligand_smiles",
            "value": "[H]O[C@@]1([H])[C@@]([H])(O[H])[C@]([H])(n2c([H])nc3c(N([H])[H])nc([H])nc32)O[C@]1([H])C([H])([H])[S@@+](C([H])([H])[H])C([H])([H])C([H])([H])[C@@]([H])(C(=O)[O-])N([H])[H]",
            "chain_ids": ["C"],
        },
    ],
    model="boltz-2.1",
    name="asmt_tucunaca_tabaco_rr",
    properties={
        "affinity": {
            "binder": "B"
        }   
    }
)

print(run_dir)