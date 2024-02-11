# Scripts for creating (and saving) isomers from a molecule by reflecting it for each axis.

import numpy as np  
from hsr.pre_processing import *
from hsr.pca_transform import * 
from hsr.fingerprint import *
from hsr.similarity import *
from hsr.utils import *
from trials.perturbations import *
import os 
from rdkit import Chem

def print_3d_coordinates(mol):
    print(f"3D coordinates \n")
    conf = mol.GetConformer()
    for atom in mol.GetAtoms():
        pos = conf.GetAtomPosition(atom.GetIdx())
        print(f"({pos.x:.2f}, {pos.y:.2f}, {pos.z:.2f})")


np.set_printoptions(precision=4, suppress=True)

cwd = os.getcwd()
# PRE-PROCESSING
# List of molecules from SDF file
molecules = load_molecules_from_sdf(f'{cwd}/sd_data/Co_diamine_2.sdf', removeHs=False, sanitize=False)


## SIMILARITY BETWEEN ISOMERS CONSTRUCTED BY REFLECTION FOR EACH AXIS ##
axis = ['x', 'y', 'z']
## Reflect molecule to construct an optical isomer ##
mol = molecules[0]


# Open an SDF writer and specify the output SDF file
with Chem.SDWriter(f'{cwd}/sd_data/lamda_delta_isomerism_2.sdf') as writer:
    # Write the original molecule
    writer.write(mol)

    # Reflect the molecule for each axis and compute similarity
    for ax in axis:
        mol1 = reflect_molecule_coordinate(mol, coordinate=ax)
        writer.write(mol1)

writer.close()



