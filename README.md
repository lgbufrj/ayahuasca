# Ayahuasca Biosynthesis Analysis

A comprehensive computational pipeline for analyzing the biosynthetic pathways of psychoactive alkaloid compounds found in *Psychotria viridis* and related species, with a focus on harmine, mescaline, and yuremamine biosynthesis.

## Overview

This project combines genomic, transcriptomic, and structural data to:

- Identify and characterize enzymes involved in alkaloid biosynthesis
- Perform homology searches against multiple plant genomes
- Conduct sequence alignments and phylogenetic analysis
- Generate structural models for computational studies
- Compute thermodynamic properties of pathway intermediates
- Produce publication-ready figures and tables

## Project Structure

```
├── documentation.md              # Detailed project documentation (PT-BR)
├── instrucoes.txt               # Step-by-step instructions (PT-BR)
├── inputs/                      # Input pathway definitions (YAML format)
│   ├── harmine_biosynthesis.yaml
│   ├── mescaline_biosynthesis.yaml
│   └── yuremamine_biosynthesis.yaml
├── pathways/                    # Analysis results organized by pathway
│   ├── harmine_biosynthesis/
│   │   ├── compounds/           # Chemical compound data
│   │   ├── genome/              # Genomic sequences & BLAST databases
│   │   ├── proteins/            # Per-protein analysis results
│   │   ├── transcriptome/       # Transcriptomic data
│   │   └── paper/               # Publication figures and tables
│   ├── mescaline_biosynthesis/
│   └── yuremamine_biosynthesis/
└── scripts/                     # Python analysis pipeline
    ├── genomics/                # Genome analysis (BLAST, alignment, etc.)
    ├── structural/              # Molecular structure processing
    ├── thermo/                  # Thermodynamic calculations (ORCA)
    ├── util/                    # Utility functions
    ├── data.py                  # Central configuration & IDs
    ├── execute.py               # Master execution script
    ├── requirements.txt         # Python dependencies
    └── pathway_schema.json      # YAML schema validation
```

## Getting Started

### Prerequisites

- Python 3.8+
- Virtual environment manager (venv, conda, etc.)
- External tools (optional, for full pipeline):
  - BLAST/DIAMOND (for homology searches)
  - MAFFT (for sequence alignment)
  - ORCA (for quantum chemistry calculations)
  - AlphaFold/OmegaFold (for structure prediction)

### Installation

1. **Clone or set up the repository:**
   ```bash
   cd /home/pedro/Desktop/projects/ayahuasca
   ```

2. **Create a Python virtual environment:**
   ```bash
   python -m venv .venv
   source .venv/bin/activate  # On Windows: .venv\Scripts\activate
   ```

3. **Install dependencies:**
   ```bash
   pip install -r scripts/requirements.txt
   ```

## Usage

### Basic Workflow

All scripts should be executed from the `scripts/` directory:

```bash
cd scripts/
python execute.py          # Run the main pipeline
```

### Running Individual Scripts

```bash
# Genomics pipeline
python genomics/run_blast.py           # Search for protein homologs
python genomics/alignment.py           # Perform sequence alignments
python genomics/generate_tree.py       # Build phylogenetic trees

# Structural analysis
python structural/convert_sdf_to_pdbqt.py   # Prepare structures for docking
python structural/run_boltz.py              # Molecular dynamics simulations

# Thermodynamics
python thermo/generate_orca_inp.py     # Prepare quantum chemistry inputs
python thermo/run_orca.py              # Execute ORCA calculations
python thermo/parse_orca_output.py     # Parse and analyze results

# Data generation
python generate_tables.py               # Create publication tables
```

### Configuration

Edit `scripts/data.py` to configure:

- Current project pathway: `CURRENT_PROJECT`
- Base paths and directories
- Protein identifiers (UniProt, local names)
- Compound identifiers (PubChem, SMILES)
- External database URLs

## Pipeline Overview

### 1. Genomic Data Preparation
- Download or provide genomic FASTA sequences
- Build BLAST databases: `genomics/make_blast_db.py`

### 2. Homology Searches
- Run BLAST/DIAMOND against target genomes
- Filter results by E-value and identity thresholds
- Extract matching sequences

### 3. Sequence Analysis
- Multiple sequence alignment (MAFFT)
- Phylogenetic tree construction
- Domain annotation

### 4. Structural Modeling
- Download structures from PDB
- Convert between molecular formats (SDF, PDB, PDBQT)
- Prepare for molecular dynamics or docking

### 5. Thermodynamic Analysis
- Generate quantum chemistry input files
- Run ORCA DFT calculations
- Calculate ΔG, ΔH, ΔS for pathway intermediates

### 6. Publication
- Generate summary tables
- Create publication-quality figures
- Organize results in `paper/` directory

## Important Notes

### File Naming Conventions

- Use lowercase English names with underscores (no spaces)
- Examples: `tryptamine_oxidase.fasta`, `protein_name/`
- Maintain consistent IDs in `scripts/data.py`

### Directory Organization

Each protein should have its own folder in `proteins/<protein_name>/`:
```
proteins/tdc/
├── reference/
│   ├── sequence.fasta
│   └── structure/
├── results/
│   ├── blast_results.txt
│   ├── alignment.fasta
│   └── tree.nwk
└── annotations/
    ├── domains.tsv
    └── predictions.gff
```

### Executing from Different Locations

Scripts use relative paths from the `scripts/` directory. If running from elsewhere, ensure proper path configuration in your shell environment.

## Key Components

### `genomics/`
- `run_blast.py` - Execute BLAST searches
- `alignment.py` - Multiple sequence alignment
- `generate_tree.py` - Phylogenetic analysis
- `extract_ids_from_blast.py` - Parse BLAST results
- `translate_genome.py` - Convert DNA → Protein sequences

### `structural/`
- `convert_sdf_to_pdbqt.py` - Prepare ligands for docking
- `convert_pdb_to_pdbqt.py` - Prepare receptors
- `run_boltz.py` - Molecular dynamics

### `thermo/`
- `generate_orca_inp.py` - Create quantum chemistry inputs
- `run_orca.py` - Execute DFT calculations
- `compute_delta_G.py` - Analyze thermodynamic results

## Dependencies

Key Python packages:
- **BioPython** - Sequence analysis and BLAST parsing
- **pandas** - Data manipulation and tables
- **matplotlib/seaborn** - Visualization
- **rdkit** - Cheminformatics
- **numpy/scipy** - Numerical computing

See `scripts/requirements.txt` for complete list and versions.

## Output Structure

Results are organized by pathway and analysis type:

```
pathways/harmine_biosynthesis/
├── compounds/          # Chemical structures and properties
├── genome/             # Genomic sequences and BLAST results
├── proteins/           # Per-protein BLAST hits and alignments
├── transcriptome/      # Transcript sequences
└── paper/              # Publication-ready outputs
    ├── tables/         # Summary statistics
    └── plots/          # Figures
```

## Troubleshooting

| Issue | Solution |
|-------|----------|
| Import errors | Ensure virtual environment is activated and `requirements.txt` installed |
| File not found | Verify you're running scripts from `scripts/` directory |
| BLAST/MAFFT not found | Install external tools and add to PATH |
| Permission denied | Check file permissions and ensure write access to `pathways/` |

## Contributing

When adding new features:
1. Create scripts in appropriate subdirectory (`genomics/`, `structural/`, `thermo/`)
2. Update `scripts/data.py` for new IDs or configurations
3. Document dependencies in `requirements.txt`
4. Add entry to `execute.py` if part of main pipeline
5. Update this README

## References

### Documentation
- [BioPython Tutorial](https://biopython.org/wiki/Documentation)
- [BLAST Documentation](https://www.ncbi.nlm.nih.gov/books/NBK279690/)
- [MAFFT Manual](https://mafft.cbrc.jp/alignment/software/)
- [RDKit Docs](https://www.rdkit.org/docs/)
- [ORCA Documentation](https://www.orcasmp.de/)

### External Tools
- NCBI BLAST: https://blast.ncbi.nlm.nih.gov/
- UniProt: https://www.uniprot.org/
- PubChem: https://pubchem.ncbi.nlm.nih.gov/
- PDB: https://www.rcsb.org/

## License

[Add license information here]

## Contact

[Add contact information here]

---

For detailed information, see `documentation.md` (Portuguese) or `instrucoes.txt` for step-by-step instructions.
