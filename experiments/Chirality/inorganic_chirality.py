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
molecules = load_molecules_from_sdf(f'{cwd}/experiments/Chirality/lambda_delta_isomers.sdf', removeHs=False, sanitize=False)

### ROTATE MOLECULES ###
rotated_molecules = []
for molecule in molecules:
    angle1 = np.random.randint(0, 360)
    angle2 = np.random.randint(0, 360)
    angle3 = np.random.randint(0, 360)
    mol = rotate_molecule(molecule, angle1, angle2, angle3)
    rotated_molecules.append(mol)
    
representations = [None, DEFAULT_FEATURES]

print(f'\nHSR similarity of inorganic enantiomers: \n')
for representation in representations:
    fingerprints = [generate_fingerprint_from_molecule(molecule, representation, scaling='matrix', chirality=True)[0] for molecule in rotated_molecules]

    similarity = compute_similarity_score(fingerprints[0], fingerprints[1])
    if representation is None:
        print(f"HSR (3D): {similarity:.4f}")
    else:
        print(f"HSR (6D): {similarity:.4f}")
