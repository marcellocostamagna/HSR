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

print(f'\nHSR similarity of organic enantiomers: \n')
sorted_files = sorted(os.listdir(f'{cwd}/experiments/Chirality/organic_molecules'), key=lambda x: int(x.split('-')[0]))
for file in sorted_files:
    if file.endswith('.sdf'):
        molecules = load_molecules_from_sdf(f'{cwd}/experiments/Chirality/organic_molecules/{file}', removeHs=False, sanitize=False)
    
    ### ROTATE MOLECULES ###
    rotated_molecules = []
    for molecule in molecules:
        angle1 = np.random.randint(0, 360)
        angle2 = np.random.randint(0, 360)
        angle3 = np.random.randint(0, 360)
        mol = rotate_molecule(molecule, angle1, angle2, angle3)
        rotated_molecules.append(mol)
        
    fingerprints = [generate_fingerprint_from_molecule(molecule, None, scaling='matrix', chirality=True)[0] for molecule in rotated_molecules]

    # COMPARE MOLECULES
    similarity = compute_similarity_score(fingerprints[0], fingerprints[1])
    print(f"{file[:-4]}: {similarity:.4f}\n")