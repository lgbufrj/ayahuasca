"""
    - Criar um banco de dados para o BLAST a partir do FASTA com as sequencias
        makeblastdb -in genome.cds_maririT.fasta -dbtype nucl -out cds_database
"""

import subprocess
import data
import os

GENOME_TYPES = ["phased", "non_phased"]
SEQUENCE_TYPES = ["cds", "prot"]

organisms = [org if data.organisms[org]["ref"]==False else None for org in data.organisms]

for organism in organisms:
    
    if organism is None: continue
    
    for genome_type in GENOME_TYPES:
        for sequence_type in SEQUENCE_TYPES:
            
            # Paths
            genome_fasta_path = f"{data.GENOME_PATH}/{organism}/{genome_type}/{sequence_type}.fasta"
            blast_db_path = f"{data.GENOME_PATH}/{organism}/{genome_type}/blast/{sequence_type}_{genome_type}_db"

            # Create necessary directories
            os.makedirs(os.path.dirname(blast_db_path), exist_ok=True)    
            
            db_type = "nucl" if sequence_type == "cds" else "prot"
            
            # Command to create BLAST database
            cmd = [
                "makeblastdb",
                "-in", genome_fasta_path,
                "-dbtype", db_type,
                "-out", blast_db_path
            ]
            # Execute the command
            subprocess.run(cmd, check=True)

print("BLAST database created successfully!")
