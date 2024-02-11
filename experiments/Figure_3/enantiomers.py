import os
import sys
from hsr.pre_processing import *
from hsr.pca_transform import *
import pymol
from pymol import cmd
from pymol.cgo import *
sys.path.append(os.path.abspath('../'))
from perturbations import *
from usr import get_geometrical_center

# MOLECULE SELCTION AND CHIRALITY
# Select one or the other to visualize an enantiomer at the time
# file_name = 'helicene_M.sdf'
file_name = 'helicene_P.sdf'

# Select chirality to be True or False
# chirality = True
chirality = False

cwd = os.getcwd()

molecule = load_molecules_from_sdf(f'{cwd}/{file_name}', removeHs=False, sanitize=False)[0]

molecule_coordinates = []
for atom in molecule.GetAtoms():
    position = molecule.GetConformer().GetAtomPosition(atom.GetIdx())
    molecule_coordinates.append([position.x, position.y, position.z])

ctd = get_geometrical_center(molecule_coordinates)
print(f'Geometrical center: {ctd}')

molecule_array = molecule_to_ndarray(molecule, features=None)

if chirality:
    _, _, PCs = compute_pca_using_covariance(molecule_array, chirality=chirality, return_axes=True)
else:
    _, PCs = compute_pca_using_covariance(molecule_array, return_axes=True)

# extract the three eigenvectors
pc1 = PCs[:, 0]
pc2 = PCs[:, 1]
pc3 = PCs[:, 2]

norm_pc1 = pc1 / np.linalg.norm(pc1)
norm_pc2 = pc2 / np.linalg.norm(pc2)
norm_pc3 = pc3 / np.linalg.norm(pc3)

# Initialize PyMOL
pymol.finish_launching()

# Load the molecule
molecule_path = f'{cwd}/{file_name}'

cmd.load(molecule_path, 'molecule')



w = 0.06 # cylinder width 
h = 0.25 # cone hight
d = w * 1.618 # cone base diameter

l1=5
l2=5
l3=5
cone_tip_pc1 = l1*pc1 + norm_pc1 * h
cone_tip_pc2 = l2*pc2 + norm_pc2 * h
cone_tip_pc3 = l3*pc3 + norm_pc3 * h

obj = [CYLINDER, ctd[0], ctd[1], ctd[2], l1*pc1[0], l2*pc1[1], l3*pc1[2], w, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0,
       CYLINDER, ctd[0], ctd[1], ctd[2], l2*pc2[0], l2*pc2[1], l2*pc2[2], w, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0,
       CYLINDER, ctd[0], ctd[1], ctd[2], l3*pc3[0], l3*pc3[1], l3*pc3[2], w, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0,
       CONE, l1*pc1[0], l1*pc1[1], l1*pc1[2], cone_tip_pc1[0], cone_tip_pc1[1], cone_tip_pc1[2], d, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 
       CONE, l2*pc2[0], l2*pc2[1], l2*pc2[2], cone_tip_pc2[0], cone_tip_pc2[1], cone_tip_pc2[2], d, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 
       CONE, l3*pc3[0], l3*pc3[1], l3*pc3[2], cone_tip_pc3[0], cone_tip_pc3[1], cone_tip_pc3[2], d, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0]

label_offset =  0.1
cmd.pseudoatom('label_x', pos=[cone_tip_pc1[0] + label_offset, cone_tip_pc1[1], cone_tip_pc1[2]], label='PC1')
cmd.pseudoatom('label_y', pos=[cone_tip_pc2[0], cone_tip_pc2[1] + label_offset, cone_tip_pc2[2]], label='PC2')
cmd.pseudoatom('label_z', pos=[cone_tip_pc3[0], cone_tip_pc3[1], cone_tip_pc3[2] + label_offset], label='PC3')

cmd.load_cgo(obj, 'axes')
