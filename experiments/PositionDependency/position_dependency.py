# Script to show the similarity value dependency on the position the same change occurs 

import numpy as np  
from hsr.pre_processing import *
from hsr.pca_transform import * 
from hsr.fingerprint import *
from hsr.similarity import *
from hsr.utils import *
import os 
import sys
from copy import deepcopy
sys.path.append(os.path.abspath('../'))
from perturbations import *

np.set_printoptions(precision=4, suppress=True)

def generate_all_single_deuterium_variants(molecule):
    hydrogens = [atom for atom in molecule.GetAtoms() if atom.GetSymbol() == "H"]
    modified_molecules = []

    for hydrogen in hydrogens:
        modified_molecule = deepcopy(molecule)
        hydrogen_atom = modified_molecule.GetAtomWithIdx(hydrogen.GetIdx())  # Get the corresponding hydrogen in the copied molecule
        hydrogen_atom.SetIsotope(2)  # Set isotope number to 2 for deuterium
        modified_molecules.append(modified_molecule)

    return modified_molecules

cwd = os.getcwd()

print(f'Position dependency of feature change')

original_molecule = load_molecules_from_sdf(f'{cwd}/molecule.sdf', removeHs=False, sanitize=False)[0]

modified_molecules = generate_all_single_deuterium_variants(original_molecule)

# Print the similarities between the original molecule and its perturbed versions
for i,modified_molecule in enumerate(modified_molecules):
    similarity = compute_similarity(original_molecule, modified_molecule, DEFAULT_FEATURES, scaling='matrix', chirality=False)
    print(f'Hydrogen {i+1}: {similarity:.4f}') 