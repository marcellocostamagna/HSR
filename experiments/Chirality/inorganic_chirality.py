# Scrpits collectiing examples of chirality and isomerism

import numpy as np  
from hsr.pre_processing import *
from hsr.pca_transform import * 
from hsr.fingerprint import *
from hsr.similarity import *
from hsr.utils import *
from trials.perturbations import *
import os 

cwd = os.getcwd()

### INORGANIC OPTICAL ISOMERS ###
molecules = load_molecules_from_sdf(f'{cwd}/sd_data/lamda_delta_isomerism_1.sdf', removeHs=False, sanitize=False)

### ROTATE MOLECULES ###
rotated_molecules = []
for molecule in molecules:
    angle1 = np.random.randint(0, 360)
    angle2 = np.random.randint(0, 360)
    angle3 = np.random.randint(0, 360)
    mol = rotate_molecule(molecule, angle1, angle2, angle3)
    rotated_molecules.append(mol)
    
fingerprints = [generate_fingerprint_from_molecule(molecule, DEFAULT_FEATURES, scaling='matrix', chirality=True)[0] for molecule in rotated_molecules]

# COMPARE ALL PAIRS OF MOLECULES
# Compute similarity between all pairs of fingerprints
n_molecules = len(fingerprints)
for i in range(n_molecules):
    for j in range(i+1, n_molecules):
        similarity = compute_similarity_score(fingerprints[i], fingerprints[j])
        print(f"{i+1}-{j+1}: {similarity:.4f}")