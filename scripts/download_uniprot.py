import requests
import os
from data import PROTEINS_PATH, proteins, organisms, databases

OUTPUT_PATH = PROTEINS_PATH+"/{ptn_name}/reference/{species}/structure/{ptn_name}_{uniprot_id}.fasta"

def download_uniprot():
    for ptn_name, ptn_data in proteins.items():
        for species, org_data in ptn_data.get("organisms", {}).items():
            if org_data and "uniprot_id" in org_data:
                
                db = "uniparc" if org_data['uniprot_id'].startswith("UPI") else "uniprot"

                url = f"{databases['uniprot']['url']}/{db}/{org_data['uniprot_id']}.fasta"
                filepath = OUTPUT_PATH.format(ptn_name=ptn_name, species=species, uniprot_id=org_data["uniprot_id"])

                if os.path.exists(filepath):
                    print(f"{species} {ptn_name} ({org_data['uniprot_id']}) already exists, skipping download.")
                    continue

                print(f"Downloading {species} {ptn_name} ({org_data['uniprot_id']})...")
                response = requests.get(url)

                # print(ptn_name, species, url)

                if response.status_code == 200:
                    os.makedirs(os.path.dirname(filepath), exist_ok=True)  # create directories if missing

                    with open(filepath, "w") as file:
                        file.write(response.text)
                    print(f"Downloaded {ptn_name} ({org_data['uniprot_id']}) → {filepath}")
            else:
                print(f"Failed to download {ptn_name} ({org_data['uniprot_id']})")

if __name__ == "__main__":
    download_uniprot()
