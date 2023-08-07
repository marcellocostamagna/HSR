# Script that provides the fingerprints of the molecules represented as a 6D array

import numpy as np
from scipy.spatial import distance
from scipy.stats import skew

def get_reference_points(dimensionality):
    """
    Returns the reference points
    """
    centroid = np.zeros(dimensionality)
    axis_points = np.eye(dimensionality)
    reference_points = np.vstack((centroid, axis_points))
    return reference_points

def compute_distances(molecule_data: np.ndarray):
    """
    Computes the Euclidean distance of each point from each reference point
    """
    reference_points = get_reference_points(molecule_data.shape[1])
    distances = np.empty((molecule_data.shape[0], len(reference_points)))
    for i, point in enumerate(molecule_data):
        for j, ref_point in enumerate(reference_points):
            distances[i, j] = distance.euclidean(point, ref_point)
    return distances

def compute_statistics(distances):
    means = np.mean(distances, axis=1)
    std_devs = np.std(distances, axis=1)
    skewness = skew(distances, axis=1)
    # check if skewness is nan
    skewness[np.isnan(skewness)] = 0
    
    statistics_matrix = np.vstack((means, std_devs, skewness)).T 
    # add all rows to a list   
    statistics_list = [element for row in statistics_matrix for element in row]

    return statistics_list  

def get_fingerprint(molecule_data: np.ndarray):
    """Computes the fingerprint of a given data set."""
    # Compute the Euclidean distance of each point from each reference point (which are fixed)
    distances = compute_distances(molecule_data)
    # Compute the statistics of the distances (mean, std_dev, skewness)
    fingerprint = compute_statistics(distances.T)
    
    return fingerprint
