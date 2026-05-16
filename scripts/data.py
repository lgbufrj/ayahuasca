import os
import json
from glob import glob
try:
    import yaml
except Exception:
    yaml = None
try:
    import jsonschema
except Exception:
    jsonschema = None


CURRENT_PROJECT = "harmine_biosynthesis"

BASE_PATH = "/home/pedro/Desktop/projects/ayahuasca"

SCRIPTS_PATH = f"{BASE_PATH}/scripts"

SCHEMA_PATH = os.path.join(SCRIPTS_PATH, "pathway_schema.json")

INPUTS_PATH = os.path.join(BASE_PATH, "inputs")

PATHWAYS_PATH = f"{BASE_PATH}/pathways"

PROJECT_PATH = f"{PATHWAYS_PATH}/{CURRENT_PROJECT}"

PAPER_PATH = f"{PROJECT_PATH}/paper"
COMPOUNDS_PATH = f"{PROJECT_PATH}/compounds"
GENOME_PATH = f"{PROJECT_PATH}/genome"
TRANSCRIPTOME_PATH = f"{PROJECT_PATH}/transcriptome"
PROTEINS_PATH = f"{PROJECT_PATH}/proteins"
ORGANISMS_PATH = f"{PROJECT_PATH}/organisms"


def _load_pathway_file(path):
    if path.endswith((".yml", ".yaml")):
        if yaml is None:
            raise RuntimeError("PyYAML is required to load YAML pathway files (pip install pyyaml)")
        with open(path, "r", encoding="utf-8") as fh:
            return yaml.safe_load(fh)
    elif path.endswith(".json"):
        with open(path, "r", encoding="utf-8") as fh:
            return json.load(fh)
    else:
        return None


# Load all pathway files from inputs/
pathways = {}
for p in glob(os.path.join(INPUTS_PATH, "*.*")):
    name = os.path.splitext(os.path.basename(p))[0]
    try:
        data = _load_pathway_file(p)
    except Exception as e:
        # skip invalid files but warn
        print(f"Warning: failed to load pathway file {p}: {e}")
        data = None
    if data:
        pathways[name] = data


# Optional validation against a JSON Schema in scripts/
_schema = None
if jsonschema is not None and os.path.exists(SCHEMA_PATH):
    try:
        with open(SCHEMA_PATH, "r", encoding="utf-8") as fh:
            _schema = json.load(fh)
    except Exception:
        _schema = None

if _schema is not None and jsonschema is not None:
    for pname, pdata in pathways.items():
        try:
            jsonschema.validate(instance=pdata, schema=_schema)
        except Exception as e:
            print(f"Validation error for pathway {pname}: {e}")


# EXPORT (keep the same top-level exports used across the project)
if CURRENT_PROJECT not in pathways:
    raise RuntimeError(f"CURRENT_PROJECT '{CURRENT_PROJECT}' not found in {INPUTS_PATH}. Available: {list(pathways.keys())}")


# Databases for compounds, proteins, and organisms

databases = {
    "pubchem": {"name": "PubChem", "url": "https://pubchem.ncbi.nlm.nih.gov"},
    "uniprot": {"name": "UniProt", "url": "https://www.uniprot.org"},
    "rcsb": {"name": "RCSB PDB", "url": "https://www.rcsb.org"},
    "ncbi": {"name": "NCBI", "url": "https://www.ncbi.nlm.nih.gov"},
    "kegg": {"name": "KEGG", "url": "https://www.genome.jp/kegg"},
    "metacyc": {"name": "MetaCyc", "url": "https://metacyc.org"}
}

# Intracellular compartments
compartments = {
    "cytoplasm": {"name": "Cytoplasm", "pH": 7.2, "ionic_strength": 0.15, "temperature": 298.15},
    "vacuole": {"name": "Vacuole", "pH": 5.5, "ionic_strength": 0.1, "temperature": 298.15},
    "apoplast": {"name": "Apoplast", "pH": 5.0, "ionic_strength": 0.1, "temperature": 298.15},
    "chloroplast": {"name": "Chloroplast", "pH": 8.0, "ionic_strength": 0.1, "temperature": 298.15},
    "mitochondrion": {"name": "Mitochondrion", "pH": 7.8, "ionic_strength": 0.15, "temperature": 298.15},
}

compounds = pathways[CURRENT_PROJECT].get("compounds", {})
organisms = pathways[CURRENT_PROJECT].get("organisms", {})
proteins = pathways[CURRENT_PROJECT].get("proteins", {})