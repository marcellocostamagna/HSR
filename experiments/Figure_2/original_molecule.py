import os
from hsr.pre_processing import *
from experiments.perturbations import *
import pymol
from pymol import cmd
from pymol.cgo import *

cwd = os.getcwd()

molecule = load_molecules_from_sdf(f'{cwd}/sd_data/Figure_1_molecule.sdf', removeHs=False, sanitize=False)[0]

# rotate the molecule
angle1 = 0
angle2 = 30
angle3 = -35
mol = rotate_molecule(molecule, angle1, angle2, angle3)

# translate the molecule with rdkit
mol = translate_molecule(mol,0,0,2)

# Save the molecule
# save the molecule to an SDF file
writer = Chem.SDWriter(f'{cwd}/experiments/Figure_2/original_molecule.sdf')
writer.write(mol)
writer.close()

# Initialize PyMOL
pymol.finish_launching()

# Load the molecule
molecule_path = f'{cwd}/experiments/Figure_2/original_molecule.sdf'

cmd.load(molecule_path, 'molecule')

w = 0.06 # cylinder width 
x = 3 # cylinder length
y = 4 # cylinder length
z = 5 # cylinder length
h = 0.25 # cone hight
d = w * 1.618 # cone base diameter

obj = [CYLINDER, 0.0, 0.0, 0.0,   x, 0.0, 0.0, w, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0,
       CYLINDER, 0.0, 0.0, 0.0, 0.0,   y, 0.0, w, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0,
       CYLINDER, 0.0, 0.0, 0.0, 0.0, 0.0,   z, w, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0,
       CONE,   x, 0.0, 0.0, h+x, 0.0, 0.0, d, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 
       CONE, 0.0,   y, 0.0, 0.0, h+y, 0.0, d, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 
       CONE, 0.0, 0.0,   z, 0.0, 0.0, h+z, d, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0]

x_offset = x + h + 0.1  
y_offset = y + h + 0.1  
z_offset = z + h + 0.1  
cmd.pseudoatom('label_x', pos=[x_offset, 0, 0], label='X')
cmd.pseudoatom('label_y', pos=[0, y_offset, 0], label='Y')
cmd.pseudoatom('label_z', pos=[0, 0, z_offset], label='Z')


cmd.load_cgo(obj, 'axes')


