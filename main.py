#####################################################
### Main: to define routines and acts as a script ###
#####################################################

# import packages
import numpy as np
import importlib
import shutil
from datetime import datetime
import os
import time
import csv

# import other files
import variables as v
import checks
import plot
import coordinates
import analysis


def create_analyse_configuration(name, boolean_plot_coords, boolean_plot_analysis, boolean_coords_add, boolean_analysis,
                                 labeffs, repeat_reducing, boolean_save_coordstxt,
                                 boolean_save_coordsplots, boolean_save_analysistxt, boolean_save_analysisplots, k=1):
    """ Creates coordinates and depending on the boolean_value saves, plots and analyses them.

            Parameters
            -----------
            name:                       String
                                        name to put into the name for saved files and heading of plotting
            boolean_plot_coords:        boolean
                                        to decide whether the coordinates should be plotted
            boolean_plot_analysis:      boolean
                                        to decide whether the analysis should be plotted
            boolean_coords_add:         boolean
                                        to decide whether the additional coordinates
                                        should be taken into account for plotting and saving
            boolean_analysis:           boolean
                                        to decide whether the coordinates should be analysed
            labeffs:                    numpy array
                                        array containing all the labeling efficiencies for the analysis
            repeat_reducing:            integer
                                        how often should the array be randomly reduced for each labeling efficiency
            boolean_save_coordstxt:     boolean
                                        to decide whether the coordinates should be saved
            boolean_save_coordsplots:   boolean
                                        to decide whether the plot of coordinates should be saved
            boolean_save_analysistxt:   boolean
                                        to decide whether the file with nN distances should be saved
            boolean_save_analysisplots: boolean
                                        to decide whether the plot with the analysis should be saved
            k:                          integer
                                        to specify the which k. neighbor distance should be analysed. default: 1


            Returns
            ----------
            None.

    """

    # Look at proteins
    nr_prots = int(np.round((v.prot_dens * (v.FOV_x * v.FOV_y))))  # number of proteins within the FOV

    if nr_prots == 0:
        coords_prots = None
        coords_sias_prots = None
    else:
        # coordinates of proteins
        coords_prots = coordinates.random_coordinates(nr_prots, v.FOV_x, v.FOV_y)

        # coordinates of sias on proteins
        coords_sias_prots = coordinates.coordinates_grouped_around(coords_prots, v.prob_x_glycoprots,
                                                                   v.prob_x_sias_prots,
                                                                   v.FOV_x, v.FOV_y, v.mu_sias_glycoprots,
                                                                   v.sigma_sias_glycoprots,
                                                                   v.mu_sias, v.sigma_sias)

    # Look at Lipids
    nr_glycolips = int(np.round(
        (v.lip_dens * (v.FOV_x * v.FOV_y)) * (1 - v.prob_x_glycolips[0])))  # number of glycolipids within the FOV
    nr_sias_lips = int(
        np.round(nr_glycolips * v.prob_x_sias_lips[1]))  # number of sialin acids of lipids within the FOV

    # coordinates of sias on glycolipids
    coords_sias_lips = coordinates.random_coordinates(nr_sias_lips, v.FOV_x, v.FOV_y)
    #
    #
    # combine coordinates of sias from lipids and sias from proteins
    if coords_sias_prots is None:
        coords = coords_sias_lips
    else:
        coords = np.concatenate((coords_sias_prots, coords_sias_lips), axis=0)

    # Do we plot additional coordinates?
    if boolean_coords_add:
        if coords_prots is None:
            coords_add = coords_sias_lips
        elif coords_sias_lips is None:
            coords_add = coords_prots
        else:
            coords_add = np.concatenate((coords_prots, coords_sias_lips))
    else:
        coords_add = np.empty(shape=(1, 2))

    # Save the coordinates if wanted
    if boolean_save_coordstxt:
        if not os.path.exists(v.path + r'/coordinates_' + name + '.txt'):
            np.savetxt(v.path + r'/coordinates_' + name + '.txt', coords)
        else:
            print(fr"ERROR: {v.path}/coordinates_{name}.txt was not saved, since it already existed!!!")
        if boolean_coords_add and not os.path.exists(v.path + r'/coordinates_add_' + name + '.txt'):
            np.savetxt(v.path + r'/coordinates_add_' + name + '.txt', coords_add)
        else:
            print(fr"ERROR: {v.path}/coordinates_add_{name}.txt was not saved, since it already existed!!!")
    #
    #
    if boolean_plot_coords:
        # plot the coordinates
        plot_coordinates(name, coords, coords_add, boolean_save_coordsplots)

    if boolean_analysis:
        nN_analysis(name, coords, labeffs, repeat_reducing, boolean_plot_analysis, boolean_save_analysistxt,
                    boolean_save_analysisplots, k=k)


def create_analyse_configuration_3d(name, boolean_plot_coords, boolean_plot_analysis, boolean_coords_add,
                                    boolean_analysis, labeffs, repeat_reducing, boolean_save_coordstxt,
                                    boolean_save_coordsplots, boolean_save_analysistxt, boolean_save_analysisplots,
                                    k=1):
    """ Creates 3D coordinates and depending on the boolean_value saves, plots and analyses them.

            Parameters
            -----------
            name:                       String
                                        name to put into the name for saved files and heading of plotting
            boolean_plot_coords:        boolean
                                        to decide whether the coordinates should be plotted
            boolean_plot_analysis:      boolean
                                        to decide whether the analysis should be plotted
            boolean_coords_add:         boolean
                                        to decide whether the additional coordinates
                                        should be taken into account for plotting and saving
            boolean_analysis:           boolean
                                        to decide whether the coordinates should be analysed
            labeffs:                    numpy array
                                        array containing all the labeling efficiencies for the analysis
            repeat_reducing:            integer
                                        how often should the array be randomly reduced for each labeling efficiency
            boolean_save_coordstxt:     boolean
                                        to decide whether the coordinates should be saved
            boolean_save_coordsplots:   boolean
                                        to decide whether the plot of coordinates should be saved
            boolean_save_analysistxt:   boolean
                                        to decide whether the file with nN distances should be saved
            boolean_save_analysisplots: boolean
                                        to decide whether the plot with the analysis should be saved
            k:                          integer
                                        to specify the which k. neighbor distance should be analysed. default: 1


            Returns
            ----------
            None.

    """

    # Look at proteins
    nr_prots = int(np.round((v.prot_dens * (v.FOV_x * v.FOV_y))))  # number of proteins within the FOV

    if nr_prots == 0:
        coords_prots = None
        coords_sias_prots = None
    else:
        # coordinates of glycoproteins
        coords_prots = coordinates.random_coordinates_3d(nr_prots, v.FOV_x, v.FOV_y, v.mu_glycan_height,
                                                         v.alpha_glycan_height, v.scale_glycan_height,
                                                         v.lower_glycan_height, v.upper_glycan_height)
        # coordinates of sias on proteins
        coords_sias_prots = coordinates.coordinates_grouped_around_3d(coords_prots, v.prob_x_glycoprots,
                                                                      v.prob_x_sias_prots, v.FOV_x, v.FOV_y,
                                                                      v.mu_sias_glycoprots, v.sigma_sias_glycoprots,
                                                                      v.mu_sias, v.sigma_sias, v.forbidden_angle_glycans)

    # Look at Lipids
    nr_glycolips = int(np.round(
        (v.lip_dens * (v.FOV_x * v.FOV_y)) * (1 - v.prob_x_glycolips[0])))  # number of glycolipids within the FOV
    nr_sias_lips = int(
        np.round(nr_glycolips * v.prob_x_sias_lips[1]))  # number of sialin acids of lipids within the FOV

    # coordinates of sias on glycolipids
    coords_sias_lips = coordinates.random_coordinates_3d(nr_sias_lips, v.FOV_x, v.FOV_y, 0, 0, 0,
                                                         0, 0)

    # combine coordinates of sias from lipids and sias from proteins
    if coords_sias_prots is None:
        coords = coords_sias_lips
    else:
        coords = np.concatenate((coords_sias_prots, coords_sias_lips), axis=0)

    # Do we plot additional coordinates?
    if boolean_coords_add:
        if coords_prots is None:
            coords_add = coords_sias_lips
        elif coords_sias_lips is None:
            coords_add = coords_prots
        else:
            coords_add = np.concatenate((coords_prots, coords_sias_lips), axis=0)
    else:
        coords_add = np.empty(shape=(1, 3))

    # Save the coordinates if wanted
    if boolean_save_coordstxt:
        if not os.path.exists(v.path + r'/coordinates3d_' + name + '.txt'):
            np.savetxt(v.path + r'/coordinates3d_' + name + '.txt', coords)
        else:
            print(fr"ERROR: {v.path}/coordinates3d_{name}.txt was not saved, since it already existed!!!")
        if boolean_coords_add and not os.path.exists(v.path + r'/coordinates_add3d_' + name + '.txt'):
            np.savetxt(v.path + r'/coordinates_add3d_' + name + '.txt', coords_add)
        else:
            print(fr"ERROR: {v.path}/coordinates_add3d_{name}.txt was not saved, since it already existed!!!")
    #
    #
    if boolean_plot_coords:
        # plot the coordinates
        plot_coordinates3d(name, coords, coords_add, v.FOV_x, v.FOV_y, boolean_save_coordsplots)

    if boolean_analysis:
        nN_analysis(name, coords, labeffs, repeat_reducing, boolean_plot_analysis, boolean_save_analysistxt,
                    boolean_save_analysisplots, k=k)


def load_coords(filename_coords, filename_coords_add=''):
    """ Load already saved coordinates.

           Parameters
           -----------
           filename_coords:         String
                                    path to coordinates file
           filename_coords_add:     String
                                    path to additional coordinates file

           Returns
           ----------
           coords:              2D numpy array
                                array with coordinates inside
           coords_add:          2D numpy array
                                array with additional coordinates inside
           boolean_coords_add:  boolean
                                states whether there were additional coordinates loaded
    """
    # load coordinates
    coords = np.loadtxt(filename_coords)
    if filename_coords_add != '':
        boolean_coords_add = True
        coords_add = np.loadtxt(filename_coords_add)
    else:
        boolean_coords_add = False
        coords_add = np.empty(shape=(1, 2))
    return coords, coords_add, boolean_coords_add


def load_coords3d(filename_coords, filename_coords_add=''):
    """ Load already saved coordinates.

           Parameters
           -----------
           filename_coords:         String
                                    path to coordinates file
           filename_coords_add:     String
                                    path to additional coordinates file

           Returns
           ----------
           coords:              2D numpy array
                                array with coordinates inside
           coords_add:          2D numpy array
                                array with additional coordinates inside
           boolean_coords_add:  boolean
                                states whether there were additional coordinates loaded
    """
    # load coordinates
    coords = np.loadtxt(filename_coords)
    if filename_coords_add != '':
        boolean_coords_add = True
        coords_add = np.loadtxt(filename_coords_add)
    else:
        boolean_coords_add = False
        coords_add = np.empty(shape=(1, 3))
    return coords, coords_add, boolean_coords_add


def nN_analysis(name, coords, labeffs, repeat_reducing, boolean_plot_analysis,
                boolean_save_analysistxt, boolean_save_analysisplots, k=1):
    """ Starts the nearest neighbour analysis.

           Parameters
           -----------
           name:                        String
                                        name to put into the name for saved files and heading of plotting
           coords:                      2D numpy array
                                        coordinates that should be analysed
           labeffs:                     numpy array
                                        array containing all the labeling efficiencies for the analysis
           repeat_reducing:             integer
                                        value for how often the array should be randomly reduced for
                                        each labeling efficiency
           boolean_plot_analysis:       boolean
                                        to decide whether the analysis should be plotted
           boolean_save_analysistxt:    boolean
                                        to decide whether the file with nN distances should be saved
           boolean_save_analysisplots:  boolean
                                        to decide whether the plot with the analysis should be saved
           k:                           integer
                                        to specify the which k. neighbor distance should be analysed. default: 1

           Returns
           ----------
           None.

    """
    labeling, nN = analysis.labeling_nN_all(coords, labeffs, repeat=repeat_reducing, k=k)

    min_max_med = np.empty((len(labeling), 4))
    for i in range(len(labeling)):
        min_max_med[i, 0] = labeling[i]
        min_max_med[i, 1] = min(nN[i])
        min_max_med[i, 2] = max(nN[i])
        min_max_med[i, 3] = np.median(np.array(nN[i]))

    if boolean_plot_analysis:
        plot.plot_labeling_nN_violin('Nearest neighbor analysis ' + name, nN, labeling, path=v.path,
                                     save=boolean_save_analysisplots)

    if boolean_save_analysistxt:
        if not os.path.exists(v.path + r'/labeling_vio_' + name + '.txt'):
            np.savetxt(v.path + r'/labeling_vio_' + name + '.txt', labeling)
        else:
            print(fr"ERROR: {v.path}/labeling_vio_{name}.txt was not saved, since it already existed!!!")

        if not os.path.exists(v.path + r'/nN_vio_' + name + '.txt'):
            with open(v.path + r'/nN_vio_' + name + '.txt', 'w', newline='') as file:
                csv_writer = csv.writer(file)
                csv_writer.writerows(nN)
        else:
            print(fr"ERROR: {v.path}/nN_vio_{name}.txt was not saved, since it already existed!!!")

        if not os.path.exists(v.path + r'/min_max_med_' + name + '.txt'):
            np.savetxt(v.path + fr"/min_max_med_" + name + '.txt', min_max_med)
        else:
            print(fr"ERROR: {v.path}/min_max_med_{name}.txt was not saved, since it already existed!!!")


def nN_analysis_only(name, filename_coords, labeffs, repeat_reducing, boolean_plot_analysis,
                     boolean_save_analysistxt, boolean_save_analysisplots, k=1):
    """ Load already saved coordinates and start the nearest neighbour analysis.

           Parameters
           -----------
           name:                        String
                                        name to put into the name for saved files and heading of plotting
           filename_coords:             String
                                        location and name of saved coordinates
           labeffs:                     numpy array
                                        array containing all the labeling efficiencies for the analysis
           repeat_reducing:             integer
                                        value for how often the array should be randomly reduced for each labeling efficiency
           boolean_plot_analysis:       boolean
                                        to decide whether to plot the analysis
           boolean_save_analysistxt:    boolean
                                        to decide whether the file with nN distances should be saved
           boolean_save_analysisplots:  boolean
                                        to decide whether the plot with the analysis should be saved
           k:                           integer
                                        to specify the which k. neighbor distance should be analysed. default: 1

           Returns
           ----------
           None.
    """
    # load coordinates
    coords, coords_add, boolean_coords_add = load_coords(filename_coords)

    # analyse them
    nN_analysis(name, coords, labeffs, repeat_reducing, boolean_plot_analysis, boolean_save_analysistxt,
                boolean_save_analysisplots, k)


def nN_analysis_only3d(name, filename_coords, labeffs, repeat_reducing, boolean_plot_analysis,
                       boolean_save_analysistxt, boolean_save_analysisplots, k=1):
    """ Load already saved  3D coordinates and start the nearest neighbour analysis.

           Parameters
           -----------
           name:                        String
                                        name to put into the name for saved files and heading of plotting
           filename_coords:             String
                                        location and name of saved coordinates
           labeffs:                     numpy array
                                        array containing all the labeling efficiencies for the analysis
           repeat_reducing:             integer
                                        value for how often the array should be randomly reduced for each labeling efficiency
           boolean_plot_analysis:       boolean
                                        to decide whether to plot the analysis
           boolean_save_analysistxt:    boolean
                                        to decide whether the file with nN distances should be saved
           boolean_save_analysisplots:  boolean
                                        to decide whether the plot with the analysis should be saved
           k:                           integer
                                        to specify the which k. neighbor distance should be analysed. default: 1

           Returns
           ----------
           None.
    """
    # load coordinates
    coords, coords_add, boolean_coords_add = load_coords3d(filename_coords)
    print('coordinatates loaded')

    # analyse them
    nN_analysis(name, coords, labeffs, repeat_reducing, boolean_plot_analysis, boolean_save_analysistxt,
                boolean_save_analysisplots, k)


def plot_nN_analysis_only(name, filepath_labeling, filepath_nN, path_save=v.path, y_lim=0.0):
    """ Plot already analyzed nN distributions.

               Parameters
               -----------
               name:                        String
                                            name to put into the name for saved files and heading of plotting
               filepath_labeling:           String
                                            location and name of labeling file
               filepath_nN:                 String
                                            location and name of file with nN distances
               path_save:                   String
                                            location where the plot should be saved
               y_lim:                       float
                                            maximum nN distance displayed; if 0 (default), all are depicted

               Returns
               ----------
               None.
        """
    labeling = np.loadtxt(filepath_labeling)
    print('labeling loaded')

    with open(filepath_nN, 'r') as file:
        csv_reader = csv.reader(file, delimiter=',')
        nN = [list(map(float, row)) for row in csv_reader]

    print('nN distances loaded')

    try:
        len(labeling)
        plot.plot_labeling_nN_violin(name, nN, labeling, path=path_save, save=True, ylim=y_lim)
    except:
        plot.plot_labeling_nN_violin(name, nN, [labeling], path=path_save, save=True, ylim=y_lim)


def plot_coordinates(name, coords, coords_add, boolean_save_coordsplots, FOV_x=v.FOV_x, FOV_y=v.FOV_y):
    """ Plots 2D coordinates.

           Parameters
           -----------
           name:                        String
                                        name to put into the name for saved files and heading of plotting
           coords:                      2D numpy array
                                        coordinates that should be plotted
           coords_add:                  2D numpy array
                                        coordinates that could be additionally plotted for a sanity check
           boolean_save_coordsplots:    boolean
                                        to decide whether the plot should be saved
           FOV_x:                       float
                                        FOV in x-direction (in [µm])
           FOV_y:                       float
                                        FOV in y-direction (in [µm])
           Returns
           ----------
            None.

    """

    plot.plot_single(name, coords, coords_add=coords_add, FOV_x=FOV_x, FOV_y=FOV_y, path=v.path,
                         save=boolean_save_coordsplots)  # plot the possible distribution of sialin acids


def plot_coordinates3d(name, coords, coords_add, FOV_x, FOV_y, boolean_save_coordsplots):
    """ Plots 3D coordinates.

           Parameters
           -----------
           name:                        String
                                        name to put into the name for saved files and heading of plotting
           coords:                      2D numpy array
                                        coordinates that should be plotted
           coords_add:                  2D numpy array
                                        coordinates that could be additionally plotted for a sanity check
           FOV_x:                       float
                                        FOV in x-direction (in [µm])
           FOV_y:                       float
                                        FOV in y-direction (in [µm])
           boolean_save_coordsplots:    boolean
                                        to decide whether the plot should be saved

           Returns
           ----------
            None.

    """

    plot.plot_single_3d(name, coords, v.path, coords_add, FOV_x, FOV_y,
                        boolean_save_coordsplots)  # plot the possible distribution of sialin acids


def plot_coordinates_only(name, filename_coords, filename_coords_add='', FOV_x=v.FOV_x, FOV_y=v.FOV_y,
                          boolean_save_coordsplots=False):
    """ loads and plots coordinates.

           Parameters
           -----------
           name:                        String
                                        name to put into the name for saved files and heading of plotting
           filename_coords:             String
                                        location and name of saved coordinates
           filename_coords_add:         String
                                        location and name of saved additional coordinates
           FOV_x:                       float
                                        FOV in x-direction (in [µm])
           FOV_y:                       float
                                        FOV in y-direction (in [µm])
           boolean_save_coordsplots:    boolean
                                        to decide whether the plot should be saved

           Returns
           ----------
            None.

    """
    # load coordinates
    coords, coords_add, boolean_coords_add = load_coords(filename_coords, filename_coords_add)

    # plot the coordinates

    plot.plot_single(name, coords, coords_add=coords_add, path=v.path, FOV_x=FOV_x, FOV_y=FOV_y,
                         save=boolean_save_coordsplots)  # plot the possible distribution of sialin acids


def plot_coordinates_only3d(name, filename_coords, filename_coords_add='', FOV_x=v.FOV_x, FOV_y=v.FOV_y,
                            boolean_save_coordsplots=False):
    """ loads and plots coordinates.

           Parameters
           -----------
           name:                        String
                                        name to put into the name for saved files and heading of plotting
           filename_coords:             String
                                        location and name of saved coordinates
           filename_coords_add:         String
                                        location and name of saved additional coordinates
           FOV_x:                       float
                                        FOV in x-direction (in [µm])
           FOV_y:                       float
                                        FOV in y-direction (in [µm])
           boolean_save_coordsplots:    boolean
                                        to decide whether the plot should be saved

           Returns
           ----------
            None.

    """
    # load coordinates
    coords, coords_add, boolean_coords_add = load_coords3d(filename_coords, filename_coords_add)

    # plot the coordinates
    plot.plot_single_3d(name, coords, coords_add=coords_add, path=v.path, FOV_x=FOV_x, FOV_y=FOV_y,
                        save=boolean_save_coordsplots)  # plot the possible distribution of sialin acids


def plot_coordinates_scatter_only3d(name, filename_coords, filename_coords_add='', FOV_x=v.FOV_x, FOV_y=v.FOV_y,
                                    boolean_save_coordsplots=False):
    """ loads and plots coordinates.

           Parameters
           -----------
           name:                        String
                                        name to put into the name for saved files and heading of plotting
           filename_coords:             String
                                        location and name of saved coordinates
           filename_coords_add:         String
                                        location and name of saved additional coordinates
           FOV_x:                       float
                                        FOV in x-direction (in [µm])
           FOV_y:                       float
                                        FOV in y-direction (in [µm])
           boolean_save_coordsplots:    boolean
                                        to decide whether the plot should be saved

           Returns
           ----------
            None.

    """
    # load coordinates
    coords, coords_add, boolean_coords_add = load_coords3d(filename_coords, filename_coords_add)

    # plot the coordinates
    plot.plot_scatter_3d(name, coords, coords_add=coords_add, path=v.path, FOV_x=FOV_x, FOV_y=FOV_y,
                         save=boolean_save_coordsplots)  # plot the possible distribution of sialin acids


def find_local_maxima(folderpath, name):
    """ finds local maxima of nearest neighbor distributions. x and kde represent the whole nN distribution, peakDis and
        peaky the nN distances and the corresponding kde values of the local maxima.

              Parameters
              -----------
              folderpath:           String
                                    path to folder where the nN distances are saved and the results get saved
              name:                 String
                                    name of the simulation to analyze

              Returns
              ----------
               None.
       """
    labeling = np.loadtxt(folderpath + fr"\labeling_vio_" + name + fr"_0.txt")

    with open(folderpath + fr"\nN_vio_" + name + "_0.txt", 'r') as file:
        csv_reader = csv.reader(file, delimiter=',')
        nN = [list(map(float, row)) for row in csv_reader]

    x, kde, peak_dis, peaks_y = analysis.find_maxima(labeling, nN)
    if not os.path.exists(v.path + r'/localMax_x_' + name + '_0.txt'):
        with open(v.path + r'/localMax_x_' + name + '_0.txt', 'w', newline='') as file:
            csv_writer = csv.writer(file)
            csv_writer.writerows(x)
    else:
        print(fr"ERROR: {v.path}/localMax_x_{name}_0.txt was not saved, since it already existed!!!")

    if not os.path.exists(v.path + r'/localMax_kde_' + name + '_0.txt'):
        with open(v.path + r'/localMax_kde_' + name + '_0.txt', 'w', newline='') as file:
            csv_writer = csv.writer(file)
            csv_writer.writerows(kde)
    else:
        print(fr"ERROR: {v.path}/localMax_kde_{name}_0.txt was not saved, since it already existed!!!")

    if not os.path.exists(v.path + r'/localMax_peakDis_' + name + '_0.txt'):
        with open(v.path + r'/localMax_peakDis_' + name + '_0.txt', 'w', newline='') as file:
            csv_writer = csv.writer(file)
            csv_writer.writerows(peak_dis)
    else:
        print(fr"ERROR: {v.path}/localMax_peakDis_{name}_0.txt was not saved, since it already existed!!!")

    if not os.path.exists(v.path + r'/localMax_peaky_' + name + '_0.txt'):
        with open(v.path + r'/localMax_peaky_' + name + '_0.txt', 'w', newline='') as file:
            csv_writer = csv.writer(file)
            csv_writer.writerows(peaks_y)
    else:
        print(fr"ERROR: {v.path}/localMax_peaky_{name}_0.txt was not saved, since it already existed!!!")


def routine(boolean_plot_coords, boolean_plot_analysis, nr_rounds, boolean_coords_add, boolean_analysis, labeffs,
            repeat_reducing, boolean_save_coordstxt, boolean_save_coordsplots,
            boolean_save_analysistxt, boolean_save_analysisplots, k=1):
    """ Starts routine for plotting and analysing.

              Parameters
              -----------
              boolean_plot_coords:          boolean
                                            to decide whether coordinates should be plotted
              boolean_plot_analysis:        boolean
                                            to decide whether analysis should be plotted
              nr_rounds:                    integer
                                            number of times completely new random coordinates are generated
              boolean_coords_add:           boolean
                                            to decide whether additional coordinates should be considered as a sanity check
              boolean_analysis:             boolean
                                            to decide whether coordinates should be analysed
              labeffs:                      numpy array
                                            array containing all the labeling efficiencies for the analysis
              repeat_reducing:              integer
                                            value for how often the array should be randomly reduced
                                            for each labeling efficiency
              boolean_save_coordstxt:       boolean
                                            to decide whether the coordinates should be saved
              boolean_save_coordsplots:     boolean
                                            to decide whether the plot of coordinates should be saved
              boolean_save_analysistxt:     boolean
                                            to decide whether the file with nN distances should be saved
              boolean_save_analysisplots:   boolean
                                            to decide whether the plot with the analysis should be saved
              k:                            integer
                                            to specify the which k. neighbor distance should be analysed. default: 1

              Returns
              ----------
               None.
       """

    # check variables
    checks.check_variables()
    print(v.path)

    os.mkdir(v.path)

    # save variables
    shutil.copy2(r'variables.py', v.path)

    for rnd in range(nr_rounds):
        start = datetime.now()  # start measuring the time
        name = v.global_name + '_' + str(rnd)
        create_analyse_configuration(name, boolean_plot_coords, boolean_plot_analysis, boolean_coords_add,
                                     boolean_analysis, labeffs, repeat_reducing,
                                     boolean_save_coordstxt, boolean_save_coordsplots, boolean_save_analysistxt,
                                     boolean_save_analysisplots, k)
        end = datetime.now()
        print('Time needed: ' + str(end - start))


def routine3d(boolean_plot_coords, boolean_plot_analysis, nr_rounds, boolean_coords_add, boolean_analysis, labeffs,
              repeat_reducing, boolean_save_coordstxt, boolean_save_coordsplots,
              boolean_save_analysistxt, boolean_save_analysisplots, k=1):
    """ Starts routine for plotting and analysing.

              Parameters
              -----------
              boolean_plot_coords:          boolean
                                            to decide whether coordinates should be plotted
              boolean_plot_analysis:        boolean
                                            to decide whether analysis should be plotted
              nr_rounds:                    integer
                                            number of times completely new random coordinates are generated
              boolean_coords_add:           boolean
                                            to decide whether additional coordinates should be considered as a sanity check
              boolean_analysis:             boolean
                                            to decide whether coordinates should be analysed
              labeffs:                      numpy array
                                            array containing all the labeling efficiencies for the analysis
              repeat_reducing:              integer
                                            value for how often the array should be randomly reduced
                                            for each labeling efficiency
              boolean_save_coordstxt:       boolean
                                            to decide whether the coordinates should be saved
              boolean_save_coordsplots:     boolean
                                            to decide whether the plot of coordinates should be saved
              boolean_save_analysistxt:     boolean
                                            to decide whether the file with nN distances should be saved
              boolean_save_analysisplots:   boolean
                                            to decide whether the plot with the analysis should be saved
              k:                            integer
                                            to specify the which k. neighbor distance should be analysed. default: 1

              Returns
              ----------
               None.
       """

    # check variables
    checks.check_variables()
    print(v.path)

    os.mkdir(v.path)

    # save variables
    shutil.copy2(r'variables.py', v.path)

    for rnd in range(nr_rounds):
        start = datetime.now()  # start measuring the time
        name = v.global_name + '_' + str(rnd)
        create_analyse_configuration_3d(name, boolean_plot_coords, boolean_plot_analysis, boolean_coords_add,
                                        boolean_analysis, labeffs, repeat_reducing, boolean_save_coordstxt,
                                        boolean_save_coordsplots, boolean_save_analysistxt, boolean_save_analysisplots,
                                        k)
        end = datetime.now()
        print('Time needed: ' + str(end - start))


def change_variable(new_name, new_path, variable_to_change='None', new_value=None):
    """ changes variables and simulation name.

          Parameters
          -----------
          new_name:                 String
                                    new v.global_name
          new_path:                 String
                                    new v.path
          variable_to_change:       String
                                    name of variable that should be changed
          new_value:                any
                                    value of variable that should be changed

          Returns
          ----------
           None.
    """

    # read the file content
    with open(v.variables_path, 'r') as file:
        content = file.read()

    # Find the line containing the new variable assignments and change the values
    lines = content.split('\n')
    for i, line in enumerate(lines):
        if v.global_name in line:
            lines[i] = f"global_name = '{new_name}'"
        if v.path in line:
            lines[i] = f"path = r'{new_path}'"
        if variable_to_change in line:
            lines[i] = f"{variable_to_change} = {new_value}"

    # combine new file content again
    new_content = '\n'.join(lines)

    # Write the modified content back to the file
    with open(v.variables_path, 'w') as file:
        file.write(new_content)

    # reload the module to actually save the values in the storage
    importlib.reload(v)
    time.sleep(5)


routine3d(False, True, 1, True, True, [1],
        1, True, True, True,
        True, 1)

'''
change_prob_x_glycoprots = [np.array([0.3, 0, 0.7]), np.array([0.3, 0.7])]

change_prob_x_sias_prots = [np.array([0.3, 0, 0.7]), np.array([0.3, 0.7])]

for index1, array1 in enumerate(change_prob_x_glycoprots):
        array_string1 = 'np.' + repr(array1)
        change_variable(f"PeakFinder_glycoprots_{array1}_sias[]", rf'F:\3D\PeakFinderG_glycoprots_{array1}_sias[]',
                            'prob_x_glycoprots', array_string1)

        for index2, array2 in enumerate(change_prob_x_sias_prots):
            array_string2 = 'np.' + repr(array2)

            change_variable(f"PeakFinder_glycoprots_{array1}_sias_{array2}", rf'F:\3D\PeakFinderG_glycoprots_{array1}_sias_{array2}',
                            'prob_x_sias_prots', array_string2)

            routine3d(False, True, 1, True, True,
                      [1], 1, True, True, True,
                      True, 1)
'''