import numpy as np  
from hsr.pre_processing import *
from trials.perturbations import *
from experiments.usr import *
from rdkit.Chem.rdMolDescriptors import GetUSRScore, GetUSR, GetUSRCAT
import os 

cwd = os.getcwd()

print(f'\nUSR and USRCAT similarity of inorganic compounds: \n')
sorted_files = sorted(os.listdir(f'{cwd}/experiments/Inorganic/molecules'), key=lambda x: int(x.split('-')[0]))
for file in sorted_files:
    if file.endswith('.sdf'):
        molecules = load_molecules_from_sdf(f'{cwd}/experiments/Inorganic/molecules/{file}', removeHs=False, sanitize=False)
    
    ### ROTATE MOLECULES ###
    rotated_molecules = []
    for molecule in molecules:
        angle1 = np.random.randint(0, 360)
        angle2 = np.random.randint(0, 360)
        angle3 = np.random.randint(0, 360)
        mol = rotate_molecule(molecule, angle1, angle2, angle3)
        rotated_molecules.append(mol)
    
    print(f'{file[:-4]}:')
    # COMPARE ALL PAIRS OF MOLECULES
    n_molecules = len(molecules)
    for i in range(n_molecules):
        for j in range(i+1, n_molecules):
            similarity = compute_similarity(molecules[i], molecules[j])
            print(f"(in-house) USR:: {similarity:.4f}")

    usrs = [GetUSR(mol) for mol in rotated_molecules]
    for i in range(n_molecules):
        for j in range(i+1, n_molecules):
            similarity = GetUSRScore(usrs[i], usrs[j])
            print(f"(rdkit) USR:: {similarity:.4f}")
            
    usrcats = [GetUSRCAT(mol) for mol in rotated_molecules]
    for i in range(n_molecules):
        for j in range(i+1, n_molecules):
            similarity = GetUSRScore(usrcats[i], usrcats[j])
            print(f"(rdkit) USRCAT:: {similarity:.4f} \n")
