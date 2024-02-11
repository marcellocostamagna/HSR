# Script to show the similarity value dependency on the position the same change occurs 

import numpy as np  
from hsr.pre_processing import *
from hsr.pca_transform import * 
from hsr.fingerprint import *
from hsr.similarity import *
from hsr.utils import *
from trials.perturbations import *
import os 
from copy import deepcopy

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

# Load the molecule
molecules = load_molecules_from_sdf(f'{cwd}/sd_data/change_position_dependency.sdf', removeHs=False, sanitize=False)

original_molecule = molecules[0]

modified_molecules = generate_all_single_deuterium_variants(original_molecule)

# Print the similarities between the original molecule and its perturbed versions
# create a disctionary to store the similarity values for each hydrogen id 
similarity_values = []
for idx,modified_molecule in enumerate(modified_molecules):
    similarity = compute_similarity(original_molecule, modified_molecule, DEFAULT_FEATURES, scaling='matrix', chirality=False)
    similarity_values.append((idx,similarity))
    
# print all (and only) the idxs
print(similarity_values)
print(similarity_values[0][:])
print(similarity_values[0][1])



import pymol
from pymol import cmd

# Initialize PyMOL
pymol.finish_launching()

# Load the molecule
molecule_path = f'{cwd}/sd_data/change_position_dependency.sdf'
cmd.load(molecule_path, 'molecule')

def amplify_differences(value, min_val, max_val):
    normalized = (value - min_val) / (max_val - min_val)
    amplified = 1 / (1 + np.exp(-10 * (normalized - 0.5)))
    return amplified

def similarity_to_color(value, min_val, max_val):
    amplified = amplify_differences(value, min_val, max_val)
    return (1 - amplified, amplified, 0)  # Red to Green gradient

# Calculate min and max for scaling
min_val = min([val for _, val in similarity_values])
max_val = max([val for _, val in similarity_values])

model = 'molecule'
state = 1

# Apply color based on similarity
for idx, similarity in similarity_values:
    color = similarity_to_color(similarity, min_val, max_val)
    color_name = f'hydrogen_color_{idx}'
    cmd.set_color("test_color", [1.0, 0.0, 0.0])
    cmd.set_color(color_name, color)
    # cmd.color("test_color", f"id {idx}")
    cmd.color(color_name, f'/molecule//H`{idx+1}')  # Hydrogen atom selection in PyMOL might need adjustment

cmd.show('spheres', 'elem H')

# Save session
cmd.save(f'{cwd}/experiments/Figure_8/visualization_session.pse')



# change color of the molecule
# remove all hydrogens (hide (hydro))
# unlabel all atoms
# re-scale pseudos spheres
