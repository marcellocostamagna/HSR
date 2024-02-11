import os
from hsr.pre_processing import *
from hsr.pca_transform import *
from experiments.usr import get_geometrical_center
from experiments.perturbations import *
import pymol
from pymol import cmd
from pymol.cgo import *
from rdkit.Geometry import Point3D

cwd = os.getcwd()

# Select one or the other to visualize an enantiomer at the time
file_name = 'helicene_M_centered.sdf'
# file_name = 'helicene_P_centered.sdf'

molecule = load_molecules_from_sdf(f'{cwd}/experiments/Figure_3/{file_name}', removeHs=False, sanitize=False)[0]

molecule_coordinates = []
for atom in molecule.GetAtoms():
    # print atom symbol
    position = molecule.GetConformer().GetAtomPosition(atom.GetIdx())
    molecule_coordinates.append([position.x, position.y, position.z])

ctd = get_geometrical_center(molecule_coordinates)
print(f'Geometrical center: {ctd}')

# # rotate the molecule
# angle1 = 0
# angle2 = 0
# angle3 = 0
# mol = rotate_molecule(molecule, angle1, angle2, angle3)
# # translate the molecule with rdkit
# mol = translate_molecule(mol,0,0,0)

molecule_array = molecule_to_ndarray(molecule, features=None)

# print(f'{molecule_array}')


# print(f'{ctd}')

_, PCs = compute_pca_using_covariance(molecule_array, chirality=True, return_axes=True)

# print(f'{PCs}')

# extract the three eigenvectors (one for each coloumn)
pc1 = PCs[:, 0]
pc2 = PCs[:, 1]
pc3 = PCs[:, 2]

norm_pc1 = pc1 / np.linalg.norm(pc1)
norm_pc2 = pc2 / np.linalg.norm(pc2)
norm_pc3 = pc3 / np.linalg.norm(pc3)



# new_molecule = Chem.Mol(molecule)

# conf = new_molecule.GetConformer()
# for i, pos in enumerate(molecule_array):
#     conf.SetAtomPosition(i, Point3D(pos[0], pos[1], pos[2]))

# writer = Chem.SDWriter(f'{cwd}/experiments/Figure_3/helicene_M_centered.sdf')
# writer.write(molecule)
# writer.close()    

# Initialize PyMOL
pymol.finish_launching()

# Load the molecule
molecule_path = f'{cwd}/experiments/Figure_3/{file_name}'

cmd.load(molecule_path, 'molecule')

# # Iterate over atoms in the molecule and print their coordinates
# coords = []
# cmd.iterate_state(1, molecule, 'stored.coords.append((x,y,z))')

# # Print the coordinates
# for i, coord in enumerate(coords):
#     print(f'Atom {i+1}: {coord}')

w = 0.06 # cylinder width 
# x = 3 # cylinder length
# y = 4 # cylinder length
# z = 5 # cylinder length
h = 0.25 # cone hight
d = w * 1.618 # cone base diameter

# Calculate the cone tip positions by extending the vectors a bit further than the cone base
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

# x_offset =  h + 0.1  # Adjust label position beyond the tip of the cones
# y_offset =  h + 0.1  # Adjust label position beyond the tip of the cones
# z_offset =  h + 0.1  # Adjust label position beyond the tip of the cones
label_offset =  0.1
cmd.pseudoatom('label_x', pos=[cone_tip_pc1[0] + label_offset, cone_tip_pc1[1], cone_tip_pc1[2]], label='PC1')
cmd.pseudoatom('label_y', pos=[cone_tip_pc2[0], cone_tip_pc2[1] + label_offset, cone_tip_pc2[2]], label='PC2')
cmd.pseudoatom('label_z', pos=[cone_tip_pc3[0], cone_tip_pc3[1], cone_tip_pc3[2] + label_offset], label='PC3')
# cmd.pseudoatom('label_x', pos=[x_offset, 0, 0], label='PC1')
# cmd.pseudoatom('label_y', pos=[0, y_offset, 0], label='PC2')
# cmd.pseudoatom('label_z', pos=[0, 0, z_offset], label='PC3')


cmd.load_cgo(obj, 'axes')

# Save the session if needed
cmd.save(f'{cwd}/experiments/Figure_3/visualization_session.pse')

