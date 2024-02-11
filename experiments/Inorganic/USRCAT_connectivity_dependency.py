import numpy as np  
import os 
import sys
from hsr.pre_processing import *

from rdkit.Chem.rdMolDescriptors import GetUSRScore, GetUSR, GetUSRCAT
sys.path.append(os.path.abspath('../'))
from perturbations import *
from usr import *

cwd = os.getcwd()

# PRE-PROCESSING
molecules = load_molecules_from_sdf(f'{cwd}/connectivity_dependency.sdf', removeHs=False, sanitize=False)
# mol1 --> Dative bonds = single bonds
# mol2 --> Remove dative bond betwenn Co and O
# mol3 --> Remove dative bond betwenn Co and O and Co and Ns

### ROTATE MOLECULES ###
rotated_molecules = []
for molecule in molecules:
    angle1 = np.random.randint(0, 360)
    angle2 = np.random.randint(0, 360)
    angle3 = np.random.randint(0, 360)
    mol = rotate_molecule(molecule, angle1, angle2, angle3)
    rotated_molecules.append(mol)
    
n_molecules = len(molecules)
usrcats = [GetUSRCAT(mol) for mol in rotated_molecules]
print(f'USRCAT Similarity with rdkit:')
for i in range(n_molecules):
    for j in range(i+1, n_molecules):
        similarity = GetUSRScore(usrcats[i], usrcats[j])
        print(f"{i+1}-{j+1}: {similarity:.4f}")
