import rdkit
import numpy as np
import os
from nd_sim.pre_processing import *
from rdkit.Chem import AllChem
from rdkit.Chem.rdMolDescriptors import GetUSRScore, GetUSR
from oddt import toolkit
from oddt import shape

np.set_printoptions(precision=4, suppress=True)

def read_molecules_from_file(file_path):
    mols = []
    for mol in toolkit.readfile(os.path.splitext(file_path)[1][1:], file_path):
        mols.append(mol)
    return mols

cwd = os.getcwd()
# PRE-PROCESSING
# List of molecules from SDF file
molecules_rdkit = load_molecules_from_sdf(f'{cwd}/sd_data/optoiso_test/mols.sdf', removeHs=False, sanitize=False)

# Compute USR similarity (rdkit function) between all pairs of molecules
usrs = [GetUSR(mol) for mol in molecules_rdkit]
n_molecules = len(usrs)
print(f'USR Similarity with rdkit:')
for i in range(n_molecules-2):
    for j in range(i+1, n_molecules):
        similarity = GetUSRScore(usrs[i], usrs[j])
        print(f"{i+1}-{j+1}: {similarity:.4f}")

        