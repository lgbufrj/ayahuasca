import os
from data import INPUTS_PATH


TEMPLATE_TEXT = """# Template for {pathway_name} input file
# ============================================
# Description: Input parameters for the {pathway_name} pathway
# Instructions: Replace placeholder values (those with capitalized names) with actual data
# IMPORTANT: Maintain YAML formatting (indentation with 2 spaces, proper list syntax)

# COMPOUNDS SECTION
# -----------------
# List all chemical compounds involved in the pathway
# Required fields: name, pubchem_id
# Find PubChem IDs at https://pubchem.ncbi.nlm.nih.gov/
compounds:
  cpd1:
    name: "Compound 1"  # Replace with actual compound name
    pubchem_id: "000000"  # Replace with PubChem ID (numerical value as string)
  cpd2:
    name: "Compound 2"  # Replace with actual compound name
    pubchem_id: "000001"  # Replace with PubChem ID

# ORGANISMS SECTION
# -----------------
# Define organisms used in the analysis
# ref: true = reference organism (known sequences available)
# ref: false = organism of interest (may need to download genome)
organisms:
  org1:
    species: "Organism one"  # Replace with scientific name (e.g., "Homo sapiens")
    ncbi_id: "000000"  # Replace with NCBI Taxonomy ID (find at https://www.ncbi.nlm.nih.gov/taxonomy)
    ref: false  # Set to false if this is your organism of interest
    genome_files:
      prot:  # Protein sequences file
        url: ""  # Provide URL to download protein sequences (leave empty if local file)
        compressed: true  # Whether the file is compressed (gz, zip, etc.)
        final_format: "fasta"  # Expected format after decompression
        download: true  # Whether to download from URL
  org2:
    species: "Organism two"  # Replace with scientific name for reference organism
    ncbi_id: "000001"  # Replace with NCBI Taxonomy ID
    ref: true  # Set to true if this is a reference organism

# PROTEINS SECTION
# ----------------
# Define all proteins/enzymes involved in the pathway
# Each protein must have reactions that connect compounds
proteins:
  ptn1:
    name: "Protein 1"  # Full protein name
    ec_number: "0.0.0.0"  # EC number (Enzyme Commission) - find at https://www.expasy.org/enzyme/
    gene_name: "ptn1"  # Gene name or symbol
    abbreviation: "PTN 1"  # Short abbreviation for the protein
    oois: ["org1"]  # Organism(s) of interest where this protein is found
    reference_organisms:
      org2:  # Reference organism key (must match an organism defined above)
        uniprot_id: "X0XXX0"  # UniProt ID for reference sequence (find at https://www.uniprot.org/)
        rcsb_id: "0XX0"  # PDB structure ID if available (find at https://www.rcsb.org/)
    cofactors: "COFACTOR 1"  # Cofactor name (e.g., "NAD+", "Mg2+", etc.) - use null if none
    reactions:  # List of reactions catalyzed by this protein
      - id: "rr"  # Reaction ID (unique identifier)
        ref: "true"  # Whether this is a reference reaction with known kinetics
        substrates: ["cpd1"]  # List of substrate compound IDs
        products: ["cpd2"]  # List of product compound IDs
  ptn2:
    name: "Protein 2"  # Full protein name
    ec_number: "0.0.0.0"  # EC number
    gene_name: "ptn2"  # Gene name
    abbreviation: "PTN 2"  # Abbreviation
    oois: ["org1"]  # Organism(s) of interest
    reference_organisms:
      org2:  # Reference organism
        uniprot_id: "X0XXX0"  # UniProt ID
        rcsb_id: "0XX0"  # PDB structure ID
    cofactors: null  # Use null if no cofactors (or use a string for cofactor name)
    reactions:
      - id: "r1"  # Reaction ID
        ref: "true"  # Reference reaction
        substrates: ["cpd1"]  # Input compounds
        products: ["cpd2"]  # Output compounds

# NOTES:
# ------
# - All IDs (compound, organism, protein, reaction) must be unique within their sections
# - Use list notation with hyphens for arrays (e.g., ["cpd1", "cpd2"])
# - Keep indentation consistent (2 spaces per level)
# - Use quotes around values to ensure proper YAML parsing
# - For null values, use 'null' or leave empty (e.g., cofactors: null)
"""


def generate_template_input(pathway_name):
    """
    Generate a template input YAML file for a given pathway.
    
    Parameters
    ----------
    pathway_name : str
        Name of the pathway (will be used as filename)
    
    Returns
    -------
    str
        Path to the generated template file
    
    Raises
    ------
    FileExistsError
        If a template already exists for this pathway
    """
    # Ensure inputs directory exists
    os.makedirs(INPUTS_PATH, exist_ok=True)
    
    template_file = os.path.join(INPUTS_PATH, f"{pathway_name}.yaml")
    
    # Check if file already exists
    if os.path.exists(template_file):
        raise FileExistsError(
            f"Template file already exists at {template_file}. "
            f"Please delete it first or choose a different pathway name."
        )
    
    # Write template file
    with open(template_file, "w") as f:
        f.write(TEMPLATE_TEXT.format(pathway_name=pathway_name))
    
    print(f"✓ Template input file created: {template_file}")
    print(f"✓ Please edit this file and replace all placeholder values")
    
    return template_file


if __name__ == "__main__":
    import sys
    
    if len(sys.argv) < 2:
        print("Usage: python generate_template_input.py <pathway_name>")
        print("Example: python generate_template_input.py example_pathway")
        sys.exit(1)
    
    pathway_name = sys.argv[1]
    try:
        generate_template_input(pathway_name)
    except FileExistsError as e:
        print(f"Error: {e}")
        sys.exit(1)