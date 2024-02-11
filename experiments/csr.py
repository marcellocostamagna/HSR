import numpy as np
from scipy.spatial import distance
from scipy.stats import skew
import os
from hsr.pre_processing import *
from numpy.linalg import norm

cwd = os.getcwd()

def get_molecule_coordinates(molecule, removeHs=False):
    """
    Get the coordinates of a molecule
    """
    molecule_coordinates = []
    for atom in molecule.GetAtoms():
        if removeHs and atom.GetSymbol() == "H":
            continue
        # print atom symbol
        position = molecule.GetConformer().GetAtomPosition(atom.GetIdx())
        molecule_coordinates.append([position.x, position.y, position.z])
    
    return np.array(molecule_coordinates)

def get_geometrical_center(coordinates):
    """
    Get the geometrical center of a molecule
    """
    return np.mean(coordinates, axis=0)

def find_furthest_atom(coordinates, geometrical_center):
    """
    Find the furthest atom to the geometrical center of a molecule.
    """
    distances = np.linalg.norm(coordinates - geometrical_center, axis=1)
    max_distance = np.max(distances)
    # get the indices of the atoms with the maximum distance
    furthest_atoms_indices = np.where(distances == max_distance)[0]
    # get the coordinates of the furthest atoms
    furthest_atom_coordinates = coordinates[furthest_atoms_indices[0]]
    # print(f"Furthest atom number: {furthest_atoms_indices[0]+1}")  # Print the atom number
    return furthest_atom_coordinates

def find_furthest_atom_from_furthest_atom(coordinates, furthest_atoms_coordinates):
    """
    Find the furthest atom from the furthest atom of a molecule.
    """
    distances = np.linalg.norm(coordinates - furthest_atoms_coordinates, axis=1)
    max_distance = np.max(distances)
    # get the indices of the atom with the maximum distance
    furthest_atom_indices = np.where(distances == max_distance)[0]
    # get the coordinates of the furthest atom
    furthest_atom_coordinates = coordinates[furthest_atom_indices[0]]
    # print(f"Furthest atom from the furthest atom number: {furthest_atom_indices[0]+1}")  # Print the atom number
    return furthest_atom_coordinates

def cross_product_point(ctd, ftd, ftf):
    """
    Computes the cross product of the vectors ftd-ctd and ftf-ctd
    """
    v_a = ftd-ctd
    v_b = ftf-ctd
    v_c = ((norm(v_a) / (2 * norm(np.cross(v_a, v_b)))) * np.cross(v_a, v_b))
    # print if the cross product is positive or negative
    # if np.dot(v_c, v_b) > 0:
    #     print(f"Positive cross product: {v_c}")
    # else:
    #     print(f"Negative cross product: {v_c}")
    ref_point_chiral = ctd + v_c
    return ref_point_chiral

def compute_distances(molecule_coordinates, reference_points):
    """
    Computes the Euclidean distance of each point from each reference point
    """
    distances = np.empty((molecule_coordinates.shape[0], len(reference_points)))
    for i, point in enumerate(molecule_coordinates):
        for j, ref_point in enumerate(reference_points):
            distances[i, j] = distance.euclidean(point, ref_point)
    return distances

def compute_statistics(distances):
    means = np.mean(distances, axis=0)
    std_devs = np.std(distances, axis=0)
    skewness = np.cbrt(skew(distances, axis=0))

    # check if skewness is nan
    skewness[np.isnan(skewness)] = 0
    
    statistics_matrix = np.vstack((means, std_devs, skewness)).T 
    # add all rows to a list   
    statistics_list = [element for row in statistics_matrix for element in row]

    return statistics_list  

def calculate_partial_score(moments1: list, moments2:list):
    partial_score = 0
    for i in range(len(moments1)):
        partial_score += abs(moments1[i] - moments2[i])
    return partial_score / len(moments1)

def get_similarity_measure(partial_score):
    return 1/(1 + partial_score)

def compute_similarity(molecule1, molecule2, removeHs=False):
    # get the coordinates of the molecules
    molecule1_coordinates = get_molecule_coordinates(molecule1, removeHs)
    molecule2_coordinates = get_molecule_coordinates(molecule2, removeHs)
    # get the geometrical center of the molecules
    molecule1_geometrical_center = get_geometrical_center(molecule1_coordinates)
    molecule2_geometrical_center = get_geometrical_center(molecule2_coordinates)
    # find the furthest atom to the geometrical center of the molecules
    molecule1_furthest_atom = find_furthest_atom(molecule1_coordinates, molecule1_geometrical_center)
    molecule2_furthest_atom = find_furthest_atom(molecule2_coordinates, molecule2_geometrical_center)
    # find the furthest atom from the furthest atom of the molecules
    molecule1_furthest_atom_from_furthest_atom = find_furthest_atom_from_furthest_atom(molecule1_coordinates, molecule1_furthest_atom)
    molecule2_furthest_atom_from_furthest_atom = find_furthest_atom_from_furthest_atom(molecule2_coordinates, molecule2_furthest_atom)
    # find the chiral reference point of the molecules
    molecule1_chiral_reference_point = cross_product_point(molecule1_geometrical_center, molecule1_furthest_atom, molecule1_furthest_atom_from_furthest_atom)
    molecule2_chiral_reference_point = cross_product_point(molecule2_geometrical_center, molecule2_furthest_atom, molecule2_furthest_atom_from_furthest_atom)
    # generate reference points in two lists
    molecule1_reference_points = [molecule1_geometrical_center,
                                  molecule1_furthest_atom,
                                  molecule1_furthest_atom_from_furthest_atom,
                                  molecule1_chiral_reference_point]
    molecule2_reference_points = [molecule2_geometrical_center,
                                  molecule2_furthest_atom, 
                                  molecule2_furthest_atom_from_furthest_atom,
                                  molecule2_chiral_reference_point]
    # compute the distances of each point from each reference point
    molecule1_distances = compute_distances(molecule1_coordinates, molecule1_reference_points)
    molecule2_distances = compute_distances(molecule2_coordinates, molecule2_reference_points)
    # compute the statistics of the distances
    molecule1_statistics = compute_statistics(molecule1_distances)
    # print(f'Molecule 1 fingerprint: {molecule1_statistics}')
    molecule2_statistics = compute_statistics(molecule2_distances)
    # print(f'Molecule 2 fingerprint: {molecule2_statistics}')
    # compute the partial score
    partial_score = calculate_partial_score(molecule1_statistics, molecule2_statistics)
    # compute the similarity measure
    similarity = get_similarity_measure(partial_score)
    return similarity