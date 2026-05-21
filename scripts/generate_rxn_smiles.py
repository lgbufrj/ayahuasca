#!/usr/bin/env python3

from pathlib import Path

import yaml
from rdkit import Chem

from data import INPUTS_PATH, PATHWAYS_PATH

# ============================================================
# PATHS
# ============================================================

INPUTS_DIR = Path(INPUTS_PATH)
PATHWAYS_DIR = Path(PATHWAYS_PATH)

# ============================================================
# HELPERS
# ============================================================

SPECIAL_SMILES = {
    "h+": "[H+]",
    "proton": "[H+]",
    "water": "O",
    "h2o": "O",
}


def load_smiles_from_sdf(sdf_path):
    """
    Read first molecule from SDF and return canonical SMILES.
    """

    supplier = Chem.SDMolSupplier(str(sdf_path), removeHs=False)

    if len(supplier) == 0 or supplier[0] is None:
        raise ValueError(f"Invalid SDF: {sdf_path}")

    mol = supplier[0]

    return Chem.MolToSmiles(mol, canonical=True)


def build_smiles_map(compounds_dir):
    """
    Build:
        {
            "harmol": "...",
            "harmine": "...",
            ...
        }

    Searches recursively for SDF files.
    """

    smiles_map = dict(SPECIAL_SMILES)

    if not compounds_dir.exists():
        raise FileNotFoundError(
            f"Compounds directory not found: {compounds_dir}"
        )

    sdf_files = compounds_dir.rglob("*.sdf")

    found_any = False

    for sdf_file in sdf_files:

        found_any = True

        # Example:
        # harmol_68094.sdf -> harmol
        compound_name = sdf_file.stem.split("_")[0].lower()

        try:

            smiles = load_smiles_from_sdf(sdf_file)

            smiles_map[compound_name] = smiles

            print(
                f"[SDF] {compound_name} <- {sdf_file.name}"
            )

        except Exception as e:

            print(
                f"[WARNING] Failed reading "
                f"{sdf_file}: {e}"
            )

    if not found_any:
        print(
            f"[WARNING] No SDF files found in "
            f"{compounds_dir}"
        )

    return smiles_map

def reaction_to_smiles(substrates, products, smiles_map):

    substrate_smiles = []
    product_smiles = []

    for compound in substrates:

        compound_key = compound.lower()

        if compound_key not in smiles_map:
            raise KeyError(f"Missing substrate SDF: {compound}")

        substrate_smiles.append(smiles_map[compound_key])

    for compound in products:

        compound_key = compound.lower()

        if compound_key not in smiles_map:
            raise KeyError(f"Missing product SDF: {compound}")

        product_smiles.append(smiles_map[compound_key])

    return (
        ".".join(substrate_smiles)
        + ">>"
        + ".".join(product_smiles)
    )


# ============================================================
# MAIN
# ============================================================

yaml_files = sorted(INPUTS_DIR.glob("*.yaml"))

if not yaml_files:
    print("No YAML files found.")
    raise SystemExit

total_updated = 0
total_skipped = 0

for yaml_file in yaml_files:

    pathway_name = yaml_file.stem

    compounds_dir = (
        PATHWAYS_DIR
        / pathway_name
        / "compounds"
    )

    print(f"\n=== Processing {yaml_file.name} ===")

    try:
        smiles_map = build_smiles_map(compounds_dir)

    except Exception as e:
        print(f"[ERROR] {e}")
        continue

    with open(yaml_file, "r") as f:
        data = yaml.safe_load(f)

    updated = 0
    skipped = 0

    proteins = data.get("proteins", {})

    for protein_id, protein_data in proteins.items():

        reactions = protein_data.get("reactions", [])

        for reaction in reactions:

            reaction_id = reaction.get("id", "unknown")

            # Skip existing smiles
            if "smiles" in reaction:
                skipped += 1
                continue

            try:

                rxn_smiles = reaction_to_smiles(
                    reaction.get("substrates", []),
                    reaction.get("products", []),
                    smiles_map
                )

                reaction["smiles"] = rxn_smiles

                updated += 1

                print(
                    f"[OK] {protein_id}:{reaction_id}"
                )

            except Exception as e:

                print(
                    f"[ERROR] {protein_id}:{reaction_id}: {e}"
                )

    # overwrite IN PLACE
    with open(yaml_file, "w") as f:
        yaml.dump(
            data,
            f,
            sort_keys=False,
            default_flow_style=False
        )

    print(
        f"Updated: {updated} | Skipped: {skipped}"
    )

    total_updated += updated
    total_skipped += skipped

print("\n===================================")
print(f"Total updated: {total_updated}")
print(f"Total skipped: {total_skipped}")
print("Done.")