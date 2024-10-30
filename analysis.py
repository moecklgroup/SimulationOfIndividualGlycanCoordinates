########################################################
#### File with functions to analyse the coordinates ####
########################################################

import numpy as np
from scipy.spatial import KDTree
from datetime import datetime
from scipy.signal import find_peaks
from scipy.stats import gaussian_kde

import coordinates


def find_nN(coords, neighborhood=np.inf, k=1):
    """ Finds the kth nearest neighbor distances between the given coordinates and return them in a list.

           Parameters
           -----------
           coords:          2d numpy array
                            points which you want to analyse
           neighborhood:    float
                            maximal distance a point can have in order to still be considered a neighbor; infinity is default
           k:               integer
                            neighbor considered; k=1 is default -> nearest neighbor

           Returns
           ----------
           nN_distances:    array of floats
                            distances to the k-th nearest neighbors
       """

    tree = KDTree(coords)   # transform the coordinates into a KD Tree

    ds, index = tree.query(coords, k+1, distance_upper_bound=neighborhood)  # finds nearest neighbor
    nN_distances = ds[:, k]

    return nN_distances


def labeling_nN_all(coords, labeffs, neighborhood=np.inf,  k=1, repeat=1):
    """ Creates one array with the labeling efficiencies and another one with the corresponding nearest neighbor distances.

           Parameters
           -----------
           coords:          2d numpy array
                            points which you want to analyse
           neighborhood:    float
                            maximal distance a point can have in order to still be considered a neighbor; infinity is default
           k:               integer
                            neighbor considered; k=1 is default -> nearest neighbor
           labeffs:         numpy array
                            labeling efficiencies that should be analysed
           repeat:          integer
                            value for how often the array should be randomly reduced for each labeling efficiency

           Returns
           ----------
           labeling:       numpy array
                           array containing the labeling efficiencies
           nN:             2d numpy array
                           containing the nearest neighbor distances for all labeling efficiencies and repeats,
                           ready to put it into the violin plot

    """
    labeling = np.empty(len(labeffs))         # create a new array for the labeling efficiencies

    nN = []
    index = 0                   # take the highest as the starting index

    for i in labeffs:     # for all labeling efficiencies
        start_analysis = datetime.now()
        new_points = int(np.round(i*len(coords)))   # calculate how many points the new array has
        labeling[index] = i
        if new_points < 2:                         # if there are too little points, take inf as the nN distance
            nN.append(np.array([np.inf]))
            index += 1
        else:
            for j in np.arange(repeat):             # take repeat different versions of reduced coordinates
                print('reducing starts')
                start_reducing = datetime.now()
                new_coords = coordinates.coordinates_reduce(coords, new_points)
                end_reducing = datetime.now()
                print('reducing finished and took ' + str(end_reducing - start_reducing))
                additional_current_nN_distances = find_nN(new_coords, neighborhood, k)
                if j == 0:
                    current_nN_distances = additional_current_nN_distances
                else:
                    current_nN_distances = np.concatenate((current_nN_distances, additional_current_nN_distances))

            nN.append(current_nN_distances)
            index += 1
        end_analysis = datetime.now()
        print('Analysis of 1 labeling efficiency took ' + str(end_analysis - start_analysis))

    return labeling, nN


def find_maxima(labeling, nN):
    """ finds local maxima of nearest neighbor distributions. x and kde represent the whole nN distribution,
        peak_distances and peak_y the nN distances and the corresponding frequencies of the local maxima.

          Parameters
          -----------
          labeling:             numpy array
                                numpy array with investigated labeling efficiencies
          nN:                   2d numpy array
                                nearest neighbor distances

          Returns
          ----------
           x_values_all:        2d numpy array
                                all nearest neighbor distances
           kde_values_all:      2d numpy array
                                all kde values
           peak_distances_all:  2d numpy array
                                nN distances of local maxima
           peaks_y_all:         2d numpy array
                                kde values of local maxima
    """

    x_values_all = []
    kde_values_all = []
    peak_distances_all = []
    peaks_y_all = []

    for le in range(len(labeling)):
        kde = gaussian_kde(nN[le])
        x_values = np.linspace(min(nN[le]), max(nN[le]), 100)
        kde_values = kde(x_values)
        peaks, _ = find_peaks(kde_values, height=0)
        peaks_y = kde_values[peaks]
        peak_distances = x_values[peaks]

        x_values_all.append(x_values)
        kde_values_all.append(kde_values)
        peak_distances_all.append(peak_distances)
        peaks_y_all.append(peaks_y)

    return x_values_all, kde_values_all, peak_distances_all, peaks_y_all
