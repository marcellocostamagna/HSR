# Scrpits collectiing examples of chirality and isomerism

import numpy as np  
from hsr.pre_processing import *
from trials.perturbations import *
from experiments.usr import *
from rdkit.Chem.rdMolDescriptors import GetUSRScore, GetUSR, GetUSRCAT
import os 

cwd = os.getcwd()

# PRE-PROCESSING
# List of molecules from SDF file
# molecules = load_molecules_from_sdf(f'{cwd}/sd_data/coordination_isomerism_3d.sdf', removeHs=False, sanitize=False)
# molecules = load_molecules_from_sdf(f'{cwd}/sd_data/linkage_isomerism.sdf', removeHs=False, sanitize=False)
# molecules = load_molecules_from_sdf(f'{cwd}/sd_data/linkage_isomerism_avo.sdf', removeHs=False, sanitize=False)
# molecules = load_molecules_from_sdf(f'{cwd}/sd_data/fac_mer_isomerism.sdf', removeHs=False, sanitize=False)
# molecules = load_molecules_from_sdf(f'{cwd}/sd_data/cis_trans_isomerism.sdf', removeHs=False, sanitize=False)
# molecules = load_molecules_from_sdf(f'{cwd}/sd_data/cis_trans_isomerism_planar.sdf', removeHs=False, sanitize=False)
# molecules = load_molecules_from_sdf(f'{cwd}/sd_data/cis_trans_isomerisms_planar_substituted.sdf', removeHs=False, sanitize=False)
molecules = load_molecules_from_sdf(f'{cwd}/sd_data/agostic/agostic_mols.sdf', removeHs=False, sanitize=False)

### ROTATE MOLECULES ###
rotated_molecules = []
for molecule in molecules:
    angle1 = np.random.randint(0, 360)
    angle2 = np.random.randint(0, 360)
    angle3 = np.random.randint(0, 360)
    mol = rotate_molecule(molecule, angle1, angle2, angle3)
    rotated_molecules.append(mol)
    

# COMPARE ALL PAIRS OF MOLECULES
n_molecules = len(molecules)
print(f'(in-house) USR Similarity:')
for i in range(n_molecules):
    for j in range(i+1, n_molecules):
        similarity = compute_similarity(molecules[i], molecules[j])
        # similarity = compute_similarity(molecules[i], molecules[j])
        print(f"{i+1}-{j+1}: {similarity:.4f}")

usrs = [GetUSR(mol) for mol in rotated_molecules]
print(f'USR Similarity with rdkit:')
for i in range(n_molecules):
    for j in range(i+1, n_molecules):
        similarity = GetUSRScore(usrs[i], usrs[j])
        print(f"{i+1}-{j+1}: {similarity:.4f}")
        
usrcats = [GetUSRCAT(mol) for mol in rotated_molecules]
print(f'USRCAT Similarity with rdkit:')
for i in range(n_molecules):
    for j in range(i+1, n_molecules):
        similarity = GetUSRScore(usrcats[i], usrcats[j])
        print(f"{i+1}-{j+1}: {similarity:.4f}")
