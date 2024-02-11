import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.stats import skew
import math as m
from hsr import *
from hsr.pre_processing import *
from hsr.pca_transform import *
from rdkit.Geometry import Point3D


# load molecule from SDF file
molecule = load_molecules_from_sdf('sd_data/Figure_2_molecule.sdf', removeHs=False, sanitize=False)[0]

# 3D plot of the molecule's coordinates
fig = plt.figure(figsize=(10, 8))
ax1 = fig.add_subplot(1,2,1, projection='3d')
coords1 = []
for atom in molecule.GetAtoms():
    pos = molecule.GetConformer().GetAtomPosition(atom.GetIdx())
    coords1.append([pos.x, pos.y, pos.z])
    print(f"({pos.x:.2f}, {pos.y:.2f}, {pos.z:.2f})")
    ax1.scatter(pos.x, pos.y, pos.z) #, color='black') # label=atom.GetSymbol())    
# plot the point (0,0,0)
ax1.scatter(0, 0, 0, color='r', marker='x')
ax1.set_xlabel('X')
ax1.set_ylabel('Y')
ax1.set_zlabel('Z')

coords1 = np.array(coords1)

# get the molecule's coordinates transformed by PCA
molecule_array = molecule_to_ndarray(molecule, features=None, removeHs=False)
molecule_transformed = compute_pca_using_covariance(molecule_array, chirality=False)

# 3D plot of the molecule's coordinates transformed by PCA
ax2 = fig.add_subplot(1,2,2, projection='3d')
coords2 = []
for atom in molecule_transformed:
    coords2.append([atom[0], atom[1], atom[2]])
    ax2.scatter(atom[0], atom[1], atom[2])

coords2 = np.array(coords2)
# plot the point (0,0,0)
ax2.scatter(0, 0, 0, color='r', marker='x')
ax2.set_xlabel('X')
ax2.set_ylabel('Y')
ax2.set_zlabel('Z')

# Assuming arbitrary lengths for the axes for demonstration
axis_length = 5.0  # Adjust based on your data scale
# Plotting axes from the origin
for direction in np.eye(3):  # Standard basis vectors for 3D
    ax2.quiver(0, 0, 0, direction[0]*axis_length, direction[1]*axis_length, direction[2]*axis_length, color='k')

# Adjust axis limits to make units appear uniform
def set_uniform_axes(ax, coords):
    max_range = np.array([coords[:,0].max()-coords[:,0].min(), coords[:,1].max()-coords[:,1].min(), coords[:,2].max()-coords[:,2].min()]).max() / 2.0
    mid_x = (coords[:,0].max()+coords[:,0].min()) * 0.5
    mid_y = (coords[:,1].max()+coords[:,1].min()) * 0.5
    mid_z = (coords[:,2].max()+coords[:,2].min()) * 0.5
    ax.set_xlim(mid_x - max_range, mid_x + max_range)
    ax.set_ylim(mid_y - max_range, mid_y + max_range)
    ax.set_zlim(mid_z - max_range, mid_z + max_range)

set_uniform_axes(ax1, coords1)
set_uniform_axes(ax2, coords2)

plt.show()

new_molecule = Chem.Mol(molecule)

conf = new_molecule.GetConformer()
for i, pos in enumerate(molecule_transformed):
    conf.SetAtomPosition(i, Point3D(pos[0], pos[1], pos[2]))
    
# save the molecule to an SDF file
writer = Chem.SDWriter('sd_data/Figure_2_molecule_transformed.sdf')
writer.write(new_molecule)
writer.close()