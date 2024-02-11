import sys
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np

def center_molecule_coordinates(sdf_file):
    """
    Centers the coordinates of a molecule in an SDF file and writes the centered molecule
    to a new SDF file with '_centered' appended to the original file name.

    Parameters:
    - sdf_file: Path to the original SDF file containing the molecule.
    """
    # Load the molecule
    suppl = Chem.SDMolSupplier(sdf_file, removeHs=False, sanitize=False)
    mol = next(suppl)  # Assuming there's only one molecule per file

    if mol is not None:
        # Get the conformer & its coordinates
        conf = mol.GetConformer()
        coords = conf.GetPositions()
        coords = np.array(coords)

        # Calculate mean coordinates and center the molecule
        mean_coords = np.mean(coords, axis=0)
        centered_coords = coords - mean_coords

        # Update the coordinates of the conformer
        for i, atom_coords in enumerate(centered_coords):
            conf.SetAtomPosition(i, atom_coords)

        # Prepare the output file name
        output_file = sdf_file.rsplit('.', 1)[0] + '_centered.sdf'

        # Write the centered molecule to a new SDF file
        writer = Chem.SDWriter(output_file)
        writer.write(mol)
        writer.close()
        
        print(f"Centered molecule written to {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python center_molecule.py <path/to/your/molecule.sdf>")
    else:
        sdf_file = sys.argv[1]
        center_molecule_coordinates(sdf_file)
