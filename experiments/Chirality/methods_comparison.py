# Scrpits collectiing examples of chirality and isomerism

import numpy as np  
from hsr.pre_processing import *
from hsr.pca_transform import * 
from hsr.fingerprint import *
from hsr.similarity import *
from hsr.utils import *
from trials.perturbations import *
from experiments.Chirality import USR_OptIso
from experiments.Chirality import CSR
import rdkit
from rdkit.Chem import AllChem
from rdkit.Geometry.rdGeometry import Point3D


import os 

cwd = os.getcwd()

### INORGANIC OPTICAL ISOMERS ###
molecules = load_molecules_from_sdf(f'{cwd}/sd_data/lamda_delta_isomerism_2.sdf', removeHs=False, sanitize=False)

### ROTATE MOLECULES ###
rotated_molecules = []
for molecule in molecules:
    angle1 = np.random.randint(0, 360)
    angle2 = np.random.randint(0, 360)
    angle3 = np.random.randint(0, 360)
    mol = rotate_molecule(molecule, angle1, angle2, angle3)
    rotated_molecules.append(mol)

fingerprints_HSR = [generate_fingerprint_from_molecule(molecule, DEFAULT_FEATURES, scaling='matrix', chirality=True)[0] for molecule in rotated_molecules]
# fingerprints_HSR = [generate_fingerprint_from_molecule(molecule, None, scaling='matrix', chirality=True)[0] for molecule in rotated_molecules]

# COMPARE ALL PAIRS OF MOLECULES
# Compute similarity between all pairs of fingerprints
n_molecules = len(fingerprints_HSR)
for i in range(n_molecules):
    for j in range(i+1, n_molecules):
        similarity_HSR = compute_similarity_score(fingerprints_HSR[i], fingerprints_HSR[j])
        similarity_CSR = CSR.compute_similarity(rotated_molecules[i], rotated_molecules[j])
        similarity_USR_Optoiso = USR_OptIso.compute_similarity(rotated_molecules[i], rotated_molecules[j])
        # print(f"{i+1}-{j+1}: HSR:{similarity_HSR:.4f}" + f" CSR:{similarity_CSR:.4f}" + f" USR_OptoIso:{similarity_USR_Optoiso:.4f}")
        print(f"{i+1}-{j+1}: HSR:{similarity_HSR}" + f" CSR:{similarity_CSR}" + f" USR_OptIso:{similarity_USR_Optoiso}")