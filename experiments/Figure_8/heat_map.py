from rdkit import Chem
import numpy as np
import os
from hsr.pre_processing import *
from hsr.similarity import *
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors


def generate_all_single_deuterium_variants_and_compute_similarity(molecule, output_path):
    hydrogens = [atom for atom in molecule.GetAtoms() if atom.GetSymbol() == "H"]
    similarity_values = {}

    for i, hydrogen in enumerate(hydrogens):
        modified_molecule = Chem.Mol(molecule)  # Deep copy of the molecule
        hydrogen_atom = modified_molecule.GetAtomWithIdx(hydrogen.GetIdx())
        hydrogen_atom.SetIsotope(2)  # Change to deuterium for the variant

        similarity = compute_similarity(molecule, modified_molecule)
        
        # Tagging the original hydrogen in the molecule
        molecule.GetAtomWithIdx(hydrogen.GetIdx()).SetProp("hydrogenID", f"H_{i+1}")
        # get atom coordinates
        pos = molecule.GetConformer().GetAtomPosition(hydrogen.GetIdx())
        
        molecule.SetProp(f"H_{i+1}", f"{pos.x},{pos.y},{pos.z}")
        
        # Store the similarity score with this tag
        similarity_values[f"H_{i+1}"] = similarity

    # Save the original molecule with hydrogen tags to an SDF file
    writer = Chem.SDWriter(output_path)
    writer.write(molecule)
    writer.close()

    return similarity_values

# Load your molecule here
cwd = os.getcwd()
sdf_file_path = os.path.join(cwd, f'{cwd}/original_molecule.sdf')
molecules = Chem.SDMolSupplier(sdf_file_path, removeHs=False, sanitize=False)
original_molecule = molecules[0]

# Assuming the output path for the tagged molecule
sdf_output_path = os.path.join(cwd, f'{cwd}/tagged_molecule.sdf')

# Generate deuterium variants, compute similarities, and tag hydrogens
similarity_values = generate_all_single_deuterium_variants_and_compute_similarity(original_molecule, sdf_output_path)

# print(f"Similarity values: {similarity_values}")

import pymol
from pymol import cmd
from pymol.cgo import *

# Initialize PyMOL
pymol.finish_launching()

def similarity_to_color(similarity_values_dict):
    # Convert similarity values from dict values to a list
    # round the similarity values to 4 decimal places
    for key, value in similarity_values_dict.items():
        similarity_values_dict[key] = round(value, 4)
    similarity_values = list(similarity_values_dict.values())
    # print(f'similarity_values: {similarity_values}')
    # print(len(similarity_values))
    unique_similarity_values = sorted(set(similarity_values))
    
    # print(f'unique_similarity_values: {unique_similarity_values}')
    # print(len(unique_similarity_values))
    
    # Normalize the similarity values to a [0, 1] range for color mapping
    min_val, max_val = min(similarity_values), max(similarity_values)
    normalized_unique_values = [(val - min_val) / (max_val - min_val) for val in unique_similarity_values]
    # print(len(normalized_unique_values))
    # Apply a non-linear transformation to emphasize differences
    emphasized_unique_values = [np.power(val, 2) for val in normalized_unique_values]  # Squaring to spread values
    
    # Create a color map from light blue to dark blue
    cmap = mcolors.LinearSegmentedColormap.from_list("grad", ["lightblue", "darkblue"])
    
    # Get corresponding color for each emphasized value
    unique_colors = [cmap(value) for value in emphasized_unique_values]
    # print(len(unique_colors))
    
     # Map each original similarity value to its corresponding color
    value_to_color = {val: cmap(np.power((val - min_val) / (max_val - min_val), 2)) for val in unique_similarity_values}
    color_mapping = [value_to_color[val] for val in similarity_values]
    
    # Create the legend figure
    fig, ax = plt.subplots(figsize=(10, 1))
    
    # Plot colors with positions
    for i, color in enumerate(unique_colors):
        ax.plot([i, i+1], [0, 0], color=color, linewidth=10, solid_capstyle='butt')

    # Adjust x-ticks to match the position of unique similarity values
    ax.set_xticks(np.arange(0.5, len(unique_similarity_values), 1))
    ax.set_xticklabels([f"{val:.4f}" for val in unique_similarity_values], rotation=45, ha="right", fontsize=10, font='Arial')
    
    # Hide y-axis and spines
    ax.set_yticks([])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    
    plt.tight_layout()  # Adjust layout to make room for the rotated x-tick labels
    # save figure as svg
    plt.savefig(f'{os.getcwd()}/similarity_legend.svg', format='svg')
    plt.show()

    return color_mapping

# Load the tagged molecule in PyMOL
cmd.load(sdf_output_path, 'tagged_molecule')

# Load the same molecule in RDKit
rdkit_mol = Chem.SDMolSupplier(sdf_output_path)[0]

coord_to_tag = {}
for i in range(13):
    tag = f"H_{i+1}"
    pos = [float(x) for x in rdkit_mol.GetProp(tag).split(',')]
    coord_to_tag[i] = [tag, pos]
# print(f'coord_to_tag: {coord_to_tag}')

# MAKING SPHERES FOR THE HYDROGENS
# generate the colors for the hydrogen atoms
colors = similarity_to_color(similarity_values)

# for i in range(13):
#     print(f'color: {colors[i][0], colors[i][1], colors[i][2]}')
#     print(f'pos: {[coord_to_tag[i][1][0], coord_to_tag[i][1][1], coord_to_tag[i][1][2]]}')

# create a list pf spheres with the correct coordinates and colors, and radius 1
d = 0.3
obj =  [COLOR, colors[0][0], colors[0][1], colors[0][2],  SPHERE, coord_to_tag[0][1][0], coord_to_tag[0][1][1], coord_to_tag[0][1][2], d,
        COLOR, colors[1][0], colors[1][1], colors[1][2],  SPHERE, coord_to_tag[1][1][0], coord_to_tag[1][1][1], coord_to_tag[1][1][2], d,
        COLOR, colors[2][0], colors[2][1], colors[2][2],  SPHERE, coord_to_tag[2][1][0], coord_to_tag[2][1][1], coord_to_tag[2][1][2], d,
        COLOR, colors[3][0], colors[3][1], colors[3][2],  SPHERE, coord_to_tag[3][1][0], coord_to_tag[3][1][1], coord_to_tag[3][1][2], d,
        COLOR, colors[4][0], colors[4][1], colors[4][2],  SPHERE, coord_to_tag[4][1][0], coord_to_tag[4][1][1], coord_to_tag[4][1][2], d,
        COLOR, colors[5][0], colors[5][1], colors[5][2],  SPHERE, coord_to_tag[5][1][0], coord_to_tag[5][1][1], coord_to_tag[5][1][2], d,
        COLOR, colors[6][0], colors[6][1], colors[6][2],  SPHERE, coord_to_tag[6][1][0], coord_to_tag[6][1][1], coord_to_tag[6][1][2], d,
        COLOR, colors[7][0], colors[7][1], colors[7][2],  SPHERE, coord_to_tag[7][1][0], coord_to_tag[7][1][1], coord_to_tag[7][1][2], d,
        COLOR, colors[8][0], colors[8][1], colors[8][2],  SPHERE, coord_to_tag[8][1][0], coord_to_tag[8][1][1], coord_to_tag[8][1][2], d,
        COLOR, colors[9][0], colors[9][1], colors[9][2],  SPHERE, coord_to_tag[9][1][0], coord_to_tag[9][1][1], coord_to_tag[9][1][2], d,
        COLOR, colors[10][0], colors[10][1], colors[10][2],  SPHERE, coord_to_tag[10][1][0], coord_to_tag[10][1][1], coord_to_tag[10][1][2], d,
        COLOR, colors[11][0], colors[11][1], colors[11][2],  SPHERE, coord_to_tag[11][1][0], coord_to_tag[11][1][1], coord_to_tag[11][1][2], d,
        COLOR, colors[12][0], colors[12][1], colors[12][2],  SPHERE, coord_to_tag[12][1][0], coord_to_tag[12][1][1], coord_to_tag[12][1][2], d]

# show the spheres
cmd.load_cgo(obj, 'gradients')

# modify the transparency of the spheres objects
cmd.set("cgo_transparency", 0.0 , "gradients")

