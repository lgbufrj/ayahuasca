import os
import re
import json
from data import COMPOUNDS_PATH, compounds

# Extract Gibbs value from ORCA output file
def parse_orca_output(output_file_path):
    """
    Parse the ORCA output file to extract the Gibbs free energy value.
    The section is:
    -------------------
    GIBBS FREE ENERGY
    -------------------

    The Gibbs free energy is G = H - T*S

    Total enthalpy                    ...   -151.20782510 Eh 
    Total entropy correction          ...     -0.01930703 Eh    -12.12 kcal/mol
    -----------------------------------------------------------------------
    Final Gibbs free energy         ...   -151.22713212 Eh

    For completeness - the Gibbs free energy minus the electronic energy
    G-E(el)                           ...     -0.00698568 Eh     -4.38 kcal/mol


    Maximum memory used throughout the entire PROP-calculation: 26.2 MB
    
    We want to extract the total enthalpy, total entropy correction, final Gibbs free energy and G-E(el) values.
    """
    gibbs_data = {}
    with open(output_file_path, "r") as f:
        content = f.read()
        
        gibbs_section = re.search(r"GIBBS FREE ENERGY(.*?)Maximum memory used", content, re.DOTALL)
        if gibbs_section:
            section_text = gibbs_section.group(1)
            
            enthalpy_match = re.search(r"Total enthalpy\s+\.\.\.\s+([-+]?\d*\.\d+|\d+)\s+Eh", section_text)
            entropy_match = re.search(r"Total entropy correction\s+\.\.\.\s+([-+]?\d*\.\d+|\d+)\s+Eh\s+([-+]?\d*\.\d+|\d+)\s+kcal/mol", section_text)
            final_gibbs_match = re.search(r"Final Gibbs free energy\s+\.\.\.\s+([-+]?\d*\.\d+|\d+)\s+Eh", section_text)
            g_minus_eel_match = re.search(r"G-E\(el\)\s+\.\.\.\s+([-+]?\d*\.\d+|\d+)\s+Eh\s+([-+]?\d*\.\d+|\d+)\s+kcal/mol", section_text)
            
            if enthalpy_match:
                # Total Enthalpy (Eh)
                gibbs_data["total_enthalpy_eh"] = float(enthalpy_match.group(1))
            if entropy_match:
                # Total Entropy Correction (Eh)
                gibbs_data["total_entropy_correction_eh"] = float(entropy_match.group(1))
                # Total Entropy Correction (kcal/mol)
                gibbs_data["total_entropy_correction_kcalpmol"] = float(entropy_match.group(2))
            if final_gibbs_match:
                # Final Gibbs Free Energy (Eh)
                gibbs_data["final_gibbs_eh"] = float(final_gibbs_match.group(1))
            if g_minus_eel_match:
                # Gibbs Free Energy minus Electronic Energy (Eh)
                gibbs_data["g_minus_eel_eh"] = float(g_minus_eel_match.group(1))
                # Gibbs Free Energy minus Electronic Energy (kcal/mol)
                gibbs_data["g_minus_eel_kcalpmol"] = float(g_minus_eel_match.group(2))
    
    return gibbs_data

def generate_json_report(gibbs_data, output_json_path):
    os.makedirs(os.path.dirname(output_json_path), exist_ok=True)
    with open(output_json_path, "w") as json_file:
        json.dump(gibbs_data, json_file, indent=4)
    print(f"Gibbs free energy data saved to {output_json_path}")

if __name__ == "__main__":

    for cpd, info in compounds.items():
        cpd_id = info["pubchem_id"]
    
        output_file = f"{COMPOUNDS_PATH}/{cpd}/quantum/{cpd}_{cpd_id}.out"
        json_output_file = f"{COMPOUNDS_PATH}/{cpd}/thermo/{cpd}_{cpd_id}.json"
        
        if os.path.exists(output_file):
            
            # Skip if JSON already exists
            if os.path.exists(json_output_file):
                print(f"JSON report already exists for {cpd} (PubChem ID: {cpd_id}), skipping...")
                continue
            
            gibbs_info = parse_orca_output(output_file)
            print(f"Gibbs free energy data for {cpd} (PubChem ID: {cpd_id}):")
            for key, value in gibbs_info.items():
                print(f"  {key}: {value}")
                
            generate_json_report(gibbs_info, json_output_file)
            
        else:
            print(f"Output file not found: {output_file}")