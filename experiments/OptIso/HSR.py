import numpy as np  
from hsr.pre_processing import *
from hsr.pca_transform import * 
from hsr.fingerprint import *
from hsr.similarity import *
from hsr.utils import *
from trials.perturbations import *
import os 

np.set_printoptions(precision=4, suppress=True)

cwd = os.getcwd()

print(f'\nContinuity (HSR): \n')
sorted_files = sorted(os.listdir(f'{cwd}/experiments/OptIso/molecules'))
for file in sorted_files:
    if file.endswith('.sdf'):
        molecules = load_molecules_from_sdf(f'{cwd}/experiments/OptIso/molecules/{file}', removeHs=False, sanitize=False)
    
    ### ROTATE MOLECULES ###
    rotated_molecules = []
    for molecule in molecules:
        angle1 = np.random.randint(0, 360)
        angle2 = np.random.randint(0, 360)
        angle3 = np.random.randint(0, 360)
        mol = rotate_molecule(molecule, angle1, angle2, angle3)
        rotated_molecules.append(mol)

    fingerprints = [generate_fingerprint_from_molecule(molecule, features=None, scaling='matrix') for molecule in rotated_molecules]

    n_molecules = len(fingerprints)
    for i in range(n_molecules):
        for j in range(i+1, n_molecules):
            similarity = compute_similarity_score(fingerprints[i], fingerprints[j])
            print(f"{file[:-4]}: {similarity:.4f}")
            