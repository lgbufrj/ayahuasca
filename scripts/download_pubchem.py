import os
import requests
# import pubchempy as pcp
from data import pathways, compounds, databases, COMPOUNDS_PATH

FILE_PATH = "{COMPOUNDS_PATH}/{cname}/structure/{cname}_{cid}.sdf"

def download_sdf(cid, file_path, record_type="3d"):
    url = f"{databases['pubchem']['url']}/rest/pug/compound/cid/{cid}/SDF?record_type={record_type}"
    response = requests.get(url)
    response.raise_for_status()  # raise error if request failed
    with open(file_path, "wb") as f:
        f.write(response.content)

for cname, cdata in compounds.items():
    cid = cdata["pubchem_id"]
    file_path = FILE_PATH.format(COMPOUNDS_PATH=COMPOUNDS_PATH, cname=cname, cid=cid)
    
    # Skip if file already exists
    if os.path.exists(file_path):
        print(f"SDF file already exists for {cname} (CID: {cid}), skipping...")
        continue

    # make sure folder exists
    os.makedirs(os.path.dirname(file_path), exist_ok=True)

    try:
        download_sdf(cid, file_path)
        print(f"Downloaded 3D {cname} (CID: {cid}) to {file_path}")
    except Exception as e:
        print(f"Error downloading 3D {cname} (CID: {cid}): {e}")
        print(f"Trying to download 2D {cname} (CID: {cid})...")
        try:
            download_sdf(cid, file_path, record_type="2d")
            print(f"Downloaded 2D {cname} (CID: {cid}) to {file_path}")
        except Exception as e:
            print(f"Error downloading 2D {cname} (CID: {cid}): {e}")
