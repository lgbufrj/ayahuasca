from data import COMPOUNDS_PATH, compounds
import subprocess

def convert_compound(compound, input_type, output_type):
    compound_file = COMPOUNDS_PATH + f"/{compound}/structure/{compound}_{compounds[compound]['pubchem_id']}" + f".{input_type}"
    output_file = COMPOUNDS_PATH + f"/{compound}/structure/{compound}_{compounds[compound]['pubchem_id']}" + f".{output_type}"

    # cmd = ["obabel", "-isdf", compound_file, "-opdbqt", "-O", output_file]
    cmd = ["obabel", f"-i{input_type}", compound_file, f"-o{output_type}", "-O", output_file]
    subprocess.run(cmd, check=True)
        
if __name__ == "__main__":
    for compound in compounds:
        convert_compound(compound, "sdf", "pdb")