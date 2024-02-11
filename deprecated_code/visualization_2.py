import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from hsr import *
from hsr.pre_processing import *
from hsr.pca_transform import *

# Function to adjust axis limits for uniform scaling
def set_uniform_axes(ax, coords):
    max_range = np.array([coords[:,0].max()-coords[:,0].min(), 
                          coords[:,1].max()-coords[:,1].min(), 
                          coords[:,2].max()-coords[:,2].min()]).max() / 2.0
    mid_x = (coords[:,0].max()+coords[:,0].min()) * 0.5
    mid_y = (coords[:,1].max()+coords[:,1].min()) * 0.5
    mid_z = (coords[:,2].max()+coords[:,2].min()) * 0.5
    ax.set_xlim(mid_x - max_range, mid_x + max_range)
    ax.set_ylim(mid_y - max_range, mid_y + max_range)
    ax.set_zlim(mid_z - max_range, mid_z + max_range)

# Function to plot molecule coordinates
def plot_molecule(ax, molecule, title):
    coords = []
    for atom in molecule.GetAtoms():
        pos = molecule.GetConformer().GetAtomPosition(atom.GetIdx())
        coords.append([pos.x, pos.y, pos.z])
        ax.scatter(pos.x, pos.y, pos.z, color='blue')  # Original coordinates in blue
    ax.scatter(0, 0, 0, color='r', marker='x')  # Plot the point (0,0,0)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title(title)
    set_uniform_axes(ax, np.array(coords))

# Function to plot PCA-transformed molecule coordinates
def plot_transformed_molecule(ax, molecule_transformed, title):
    for atom in molecule_transformed:
        ax.scatter(atom[0], atom[1], atom[2], color='green')  # Transformed coordinates in green
    ax.scatter(0, 0, 0, color='r', marker='x')  # Plot the point (0,0,0)
    
    # Assuming arbitrary lengths for the axes for demonstration
    axis_length = 5.0  # Adjust based on your data scale
    # Plotting axes from the origin
    for direction in np.eye(3):  # Standard basis vectors for 3D
        ax.quiver(0, 0, 0, direction[0]*axis_length, direction[1]*axis_length, direction[2]*axis_length, color='k')

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title(title)
    set_uniform_axes(ax, np.array(molecule_transformed))
    

# Load molecules from SDF files
molecule_1 = load_molecules_from_sdf('sd_data/helicene_L.sdf', removeHs=False, sanitize=False)[0]
molecule_2 = load_molecules_from_sdf('sd_data/helicene_R.sdf', removeHs=False, sanitize=False)[0]

# Create a figure for the plots
fig = plt.figure(figsize=(20, 10))

# Plot original coordinates of molecule_1
ax1 = fig.add_subplot(2, 2, 1, projection='3d')
plot_molecule(ax1, molecule_1, 'Molecule 1 Original')

# Plot PCA-transformed coordinates of molecule_1
ax2 = fig.add_subplot(2, 2, 2, projection='3d')
molecule_1_array = molecule_to_ndarray(molecule_1, features=None, removeHs=False)
molecule_1_transformed = compute_pca_using_covariance(molecule_1_array, chirality=False)
plot_transformed_molecule(ax2, molecule_1_transformed, 'Molecule 1 PCA Transformed')

# Plot original coordinates of molecule_2
ax3 = fig.add_subplot(2, 2, 3, projection='3d')
plot_molecule(ax3, molecule_2, 'Molecule 2 Original')

# Plot PCA-transformed coordinates of molecule_2
ax4 = fig.add_subplot(2, 2, 4, projection='3d')
molecule_2_array = molecule_to_ndarray(molecule_2, features=None, removeHs=False)
molecule_2_transformed = compute_pca_using_covariance(molecule_2_array, chirality=False)
plot_transformed_molecule(ax4, molecule_2_transformed, 'Molecule 2 PCA Transformed')

plt.show()

