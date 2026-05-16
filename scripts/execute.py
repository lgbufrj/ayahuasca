import subprocess
import data

execute = {
    "RUN_BLAST": {
        "path": f"{data.SCRIPTS_PATH}/run_blast.py", 
        "log": "Running BLAST searches...",
        "run": False
    },
    "FILTER_BLAST_RESULTS": {
        "path": f"{data.SCRIPTS_PATH}/filter_blast_results.py", 
        "log": "Filtering BLAST results...",
        "run": True
    },
    "EXTRACT_IDS_FROM_BLAST": {
        "path": f"{data.SCRIPTS_PATH}/extract_ids_from_blast.py", 
        "log": "Extracting IDs from BLAST results...",
        "run": False
    },
    "ALIGNMENT": {
        "path": f"{data.SCRIPTS_PATH}/alignment.py", 
        "log": "Performing sequence alignments...",
        "run": False
    },
}

for script in execute:
    if(execute[script]["run"]):
        print(execute[script]["log"])
        script_command = f"python3 {execute[script]['path']}"
        subprocess.run(script_command, shell=True, check=True)
    
print("All scripts executed!")