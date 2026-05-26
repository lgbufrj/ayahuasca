from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
from data import COMPOUNDS_PATH, compounds

def draw_molecule_with_atom_numbers(
    mol,
    output_path="molecule_atoms.png",
    size=(1200, 1200)
):
    """
    Draws molecule with atom indices.
    """

    drawer = rdMolDraw2D.MolDraw2DCairo(
        size[0],
        size[1]
    )

    options = drawer.drawOptions()

    # Label atoms with indices
    for atom in mol.GetAtoms():

        options.atomLabels[
            atom.GetIdx()
        ] = f"{atom.GetSymbol()}{atom.GetIdx()}"

    rdMolDraw2D.PrepareAndDrawMolecule(
        drawer,
        mol
    )

    drawer.FinishDrawing()

    png = drawer.GetDrawingText()

    with open(output_path, "wb") as f:
        f.write(png)

    print(f"Saved atom-numbered molecule to: {output_path}")
    
if __name__ == "__main__":
    
    for compound in compounds:
        compound_dir = f"{COMPOUNDS_PATH}/{compound}/structure"
        compound_sdf = f"{compound_dir}/{compound}_{compounds[compound]['pubchem_id']}.sdf"
        compound_png = f"{compound_dir}/{compound}_{compounds[compound]['pubchem_id']}.png"
        
        mol = Chem.MolFromMolFile(
            compound_sdf,
            removeHs=True
        )
        draw_molecule_with_atom_numbers(
            mol,
            compound_png
        )