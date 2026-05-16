from data import CURRENT_PROJECT, COMPOUNDS_PATH, pathways, compounds
import subprocess

def convert_compounds():
    for compound in compounds:
        
        compound_file = COMPOUNDS_PATH + f"/{compound}/structure/{compound}_{compounds[compound]['pubchem_id']}" + ".sdf"
        output_file = COMPOUNDS_PATH + f"/{compound}/structure/{compound}_{compounds[compound]['pubchem_id']}" + ".pdbqt"

        cmd = ["obabel", "-isdf", compound_file, "-opdbqt", "-O", output_file]
        subprocess.run(cmd, check=True)
        
if __name__ == "__main__":
    convert_compounds()