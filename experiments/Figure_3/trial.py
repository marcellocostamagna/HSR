
import numpy as np
import os
from hsr.pre_processing import *
from experiments.usr import *
import pymol
from pymol import cmd
import rdkit
from rdkit.Chem import AllChem

cwd = os.getcwd()

np.set_printoptions(precision=4, suppress=True)

molecule = load_molecules_from_sdf(f'{cwd}/sd_data/helicene_isomerism.sdf', removeHs=False, sanitize=False)[0]
# molecule = load_molecules_from_sdf(f'{cwd}/sd_data/helicene_isomerism.sdf', removeHs=False, sanitize=False)[1]

molecule_coords = get_molecule_coordinates(molecule)

ctd = get_geometrical_center(molecule_coords)
print(f'Geometrical center: {ctd}')

# find the closest atom to the geometrical center of the molecule
ctc_coords = find_closest_atom(molecule_coords, ctd)

# find the furthest atom to the geometrical center of the molecule
ftc_coords = find_furthest_atom(molecule_coords, ctd)

# find the furthest atom from the furthest atom of the molecule
ftf_coords = find_furthest_atom_from_furthest_atom(molecule_coords, ftc_coords)

import pymol
from pymol import cmd

# Initialize PyMOL
pymol.finish_launching()

# Load the molecule
molecule_path = f'{cwd}/sd_data/helicene_isomerism.sdf'
cmd.load(molecule_path, 'molecule')
    
# Create pseudatoms at each point of interest
cmd.pseudoatom('geom_center', pos=ctd.tolist(), color='cyan', label='Geom Center', )
cmd.pseudoatom('closest_to_center', pos=ctc_coords.tolist(), color='yellow', label='Closest to Center')
cmd.pseudoatom('furthest_to_center', pos=ftc_coords.tolist(), color='magenta', label='Furthest to Center')
cmd.pseudoatom('furthest_to_furthest', pos=ftf_coords.tolist(), color='green', label='Furthest to Furthest')

# Markers settings

cmd.show('spheres', 'geom_center')
cmd.show('spheres', 'closest_to_center')
cmd.show('spheres', 'furthest_to_center')
cmd.show('spheres', 'furthest_to_furthest')

cmd.set('sphere_scale', 0.5, 'geom_center')
cmd.set('sphere_scale', 0.5, 'closest_to_center')
cmd.set('sphere_scale', 0.5, 'furthest_to_center')
cmd.set('sphere_scale', 0.5, 'furthest_to_furthest')

# remove all hydrogens (hide (hydro))
# cmd.hide('everything', 'hydro')

cmd.group("pseudos", "geom_center closest_to_center furthest_to_center furthest_to_furthest")

# # Optionally, adjust the view and representation
# cmd.zoom('all', buffer=20)  # Adjust zoom and buffer to fit all objects nicely
# cmd.show('spheres', 'pseudos')  # Show pseudatoms as spheres
# cmd.set('sphere_scale', 0.5, 'pseudos')  # Adjust sphere size

# Save the session if needed
cmd.save(f'{cwd}/experiments/Figure_1/visualization_session.pse')


# change color of the molecule
# remove all hydrogens (hide (hydro))
# unlabel all atoms
# re-scale pseudos spheres
