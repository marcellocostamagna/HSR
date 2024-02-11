# Scripts for calculating between all pairs of molecules in a list of molecules.

import numpy as np  
from nd_sim.pre_processing import *
from nd_sim.pca_transform import * 
from nd_sim.fingerprint import *
from nd_sim.similarity import *
from nd_sim.utils import *
from trials.perturbations import *
import os 

np.set_printoptions(precision=4, suppress=True)

cwd = os.getcwd()
# PRE-PROCESSING
# List of molecules from SDF file
molecules = []
molecule1 = load_molecules_from_sdf(f'{cwd}/sd_data/optoiso_test/ORCA/rocs/1/rocs_mol1.sdf', removeHs=False, sanitize=False)
molecule2 = load_molecules_from_sdf(f'{cwd}/sd_data/optoiso_test/ORCA/rocs/2/rocs_mol2.sdf', removeHs=False, sanitize=False)
molecules_new = [molecule1[0], molecule2[0]]

# molecules_original = load_molecules_from_sdf(f'{cwd}/sd_data/optoiso_test/rocs_problematic_mols.sdf', removeHs=False, sanitize=False)

molecules = molecules_new #+ molecules_original

print(molecules)

# ### ROTATE MOLECULES ###
rotated_molecules = []
for molecule in molecules:
    angle1 = np.random.randint(0, 360)
    angle2 = np.random.randint(0, 360)
    angle3 = np.random.randint(0, 360)
    mol = rotate_molecule(molecule, angle1, angle2, angle3)
    rotated_molecules.append(mol)

fingerprints = [generate_fingerprint_from_molecule(molecule, features=None, scaling='matrix') for molecule in rotated_molecules]

n_molecules = len(fingerprints)
print(f'Shape Similarity:')
for i in range(n_molecules):
    for j in range(i+1, n_molecules):
        similarity = compute_similarity_score(fingerprints[i], fingerprints[j])
        print(f"{i+1}-{j+1}: {similarity:.4f}")
        
fingerprints = [generate_fingerprint_from_molecule(molecule, features=DEFAULT_FEATURES, scaling='matrix') for molecule in rotated_molecules]        
    
print(f'\nSimilarity with features:')
for i in range(n_molecules):
    for j in range(i+1, n_molecules):
        similarity = compute_similarity_score(fingerprints[i], fingerprints[j])
        print(f"{i+1}-{j+1}: {similarity:.4f}")