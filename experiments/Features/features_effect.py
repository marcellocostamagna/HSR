import os
import sys
import numpy as np  
from hsr.pre_processing import *
from hsr.pca_transform import * 
from hsr.fingerprint import *
from hsr.similarity import *
from hsr.utils import *
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import re
from rdkit.Chem import Descriptors
sys.path.append(os.path.abspath('../'))
from perturbations import *

PROTON_FEATURES = {
    'protons' : extract_proton_number,
    }

CHARGE_FEATURES = {
    'protons' : extract_proton_number,
    'charge' : extract_formal_charge,
    }

NEUTRON_FEATURES = {   
    'protons' : extract_proton_number,
    'neutrons' : extract_neutron_difference_from_common_isotope,
    }    
    

np.set_printoptions(precision=4, suppress=True)
directory_path = f'{os.getcwd()}/amines' 
pattern = re.compile(r"^(?:[1-9]|1[0-4])-.+\.sdf$")  # Pattern to match files like "1-name.sdf" 

file_names = []
for file in os.listdir(directory_path):
    if pattern.match(file):
        file_names.append(file[:-4])

file_names = sorted(file_names, key=lambda x: int(os.path.basename(x).split('-')[0]))
features_set_charge = [None, PROTON_FEATURES, DEFAULT_FEATURES]
feature_set_neutrons = [None, PROTON_FEATURES, NEUTRON_FEATURES]
features_labels_charge = ['3D:[x,y,z]', '4D:[x,y,z,√p]','5D:[x,y,z,√p,q]']
features_labels_neutrons = ['3D:[x,y,z]', '4D:[x,y,z,√p]','5D:[x,y,z,√p,√n]']

similarity_scores_H = {feature: [] for feature in features_labels_charge}
similarity_scores_C13 = {feature: [] for feature in features_labels_neutrons}

distances_H_list = {feature: [] for feature in features_labels_charge}
distances_C13_list = {feature: [] for feature in features_labels_neutrons}


# Initialize a DataFrame for storing the results
similarity_table_H = pd.DataFrame(index=file_names, columns=['NumAtoms', 'MolWeight'] + features_labels_charge)
similarity_table_C13 = pd.DataFrame(index=file_names, columns=['NumAtoms', 'MolWeight'] + features_labels_neutrons)

# for file_path in file_names:
for file_name in file_names:
    file_path = os.path.join(directory_path, file_name + '.sdf')
    molecules = load_molecules_from_sdf(file_path, removeHs=False, sanitize=False)
    
    # Extract the number of atoms and molecular weight from the first molecule
    num_atoms = molecules[0].GetNumAtoms()
    mol_weight = Descriptors.MolWt(molecules[0])
    similarity_table_H.at[file_name, 'NumAtoms'] = num_atoms
    similarity_table_H.at[file_name, 'MolWeight'] = mol_weight
    similarity_table_C13.at[file_name, 'NumAtoms'] = num_atoms
    similarity_table_C13.at[file_name, 'MolWeight'] = mol_weight

    ### ROTATE MOLECULES ###
    rotated_molecules = []
    for molecule in molecules:
        angle1 = np.random.randint(0, 360)
        angle2 = np.random.randint(0, 360)
        angle3 = np.random.randint(0, 360)
        mol = rotate_molecule(molecule, angle1, angle2, angle3)
        rotated_molecules.append(mol)

    similarity_scores_H_list_dil = []
    similarity_scores_C13_list_dil = []

    for i in range(3):
        fingerprint1_charge =  generate_fingerprint_from_molecule(molecules[0], features=features_set_charge[i], scaling='matrix')
        fingerprint1_neutrons = generate_fingerprint_from_molecule(molecules[0], features=feature_set_neutrons[i], scaling='matrix')
        fingerprint2_charge = generate_fingerprint_from_molecule(molecules[1], features=features_set_charge[i], scaling='matrix')
        fingerprint3_neutrons = generate_fingerprint_from_molecule(molecules[2], features=feature_set_neutrons[i], scaling='matrix')

        similarity_H = compute_similarity_score(fingerprint1_charge, fingerprint2_charge)
        similarity_C13 = compute_similarity_score(fingerprint1_neutrons, fingerprint3_neutrons) 
            
        similarity_scores_H[features_labels_charge[i]].append(similarity_H)
        similarity_scores_C13[features_labels_neutrons[i]].append(similarity_C13)
        
        similarity_table_H.at[file_name, features_labels_charge[i]] = similarity_H
        similarity_table_C13.at[file_name, features_labels_neutrons[i]] = similarity_C13
 
    for i in range(3):
        fingerprint1_charge =  generate_fingerprint_from_molecule(molecules[0], features=features_set_charge[i], scaling='matrix')
        fingerprint1_neutrons = generate_fingerprint_from_molecule(molecules[0], features=feature_set_neutrons[i], scaling='matrix')
        fingerprint2_charge = generate_fingerprint_from_molecule(molecules[1], features=features_set_charge[i], scaling='matrix')
        fingerprint3_neutrons = generate_fingerprint_from_molecule(molecules[2], features=feature_set_neutrons[i], scaling='matrix')

        similarity_H = compute_similarity_score(fingerprint1_charge, fingerprint2_charge)
        similarity_C13 = compute_similarity_score(fingerprint1_neutrons, fingerprint3_neutrons)
        
        distance_H = ((1 - similarity_H)* len(fingerprint1_charge)) / similarity_H
        distance_C13 = ((1 - similarity_C13)* len(fingerprint1_neutrons)) / similarity_C13

        distances_H_list[features_labels_charge[i]].append(distance_H)
        distances_C13_list[features_labels_neutrons[i]].append(distance_C13)
        
        
similarity_table_H.reset_index(inplace=True)
similarity_table_H.rename(columns={'index': 'Molecules'}, inplace=True)

similarity_table_C13.reset_index(inplace=True)
similarity_table_C13.rename(columns={'index': 'Molecules'}, inplace=True)


# # Save the results to a text file with formatted string alignment
with open(f'{os.getcwd()}/similarity_scores.txt', 'w') as f:
    f.write(f'Similarity results for protonated amines (constant dilution)\n\n')
    f.write(similarity_table_H.to_string(index=False))
    f.write('\n\n')
    f.write(f'Similarity results for deuterated (C_13) amines (constant dilution)\n\n')
    f.write(similarity_table_C13.to_string(index=False))
    f.write('\n\n')

# Setting up plot aesthetics
matplotlib.rcParams['font.size'] = 10
matplotlib.rcParams['font.family'] = 'Arial'

fig, ((ax3, ax4), (ax5, ax6)) = plt.subplots(2, 2, figsize=(6.6, 6), sharex='col')

# Define a list of markers to cycle through
markers = ['o', 's', '^']

# Plotting for ax3 and ax5 (left column) with different markers
for i, label in enumerate(features_labels_charge):
    marker = markers[i % len(markers)] 
    ax3.plot(range(1, len(file_names) + 1), similarity_scores_H[label], label=label, marker=marker)
ax3.set_ylabel('Similarity Score')
ax3.legend()

for i, label in enumerate(features_labels_charge):
    marker = markers[i % len(markers)]  
    ax5.plot(range(1, len(file_names) + 1), distances_H_list[label], label=label, marker=marker)
ax5.set_xlabel('Number of Carbons')
ax5.set_ylabel('Distance')
ax5.legend()

# Plotting for ax4 and ax6 (right column) with different markers
for i, label in enumerate(features_labels_neutrons):
    marker = markers[i % len(markers)] 
    ax4.plot(range(1, len(file_names) + 1), similarity_scores_C13[label], label=label, marker=marker)
ax4.set_ylabel('Similarity Score')
ax4.legend()

for i, label in enumerate(features_labels_neutrons):
    marker = markers[i % len(markers)] 
    ax6.plot(range(1, len(file_names) + 1), distances_C13_list[label], label=label, marker=marker)
ax6.set_xlabel('Number of Carbons')
ax6.set_ylabel('Distance')
ax6.legend()

# Explicitly setting x-ticks to ensure consistency across all plots
x_ticks = range(1, len(file_names) + 1)
ax3.set_xticks(x_ticks)
ax4.set_xticks(x_ticks)
ax5.set_xticks(x_ticks)
ax6.set_xticks(x_ticks)

plt.tight_layout(pad=1.0, h_pad=0.5, w_pad=0.5) 

plt.savefig(f'{os.getcwd()}/similarity_&_distances.svg', format='svg')
plt.show()



