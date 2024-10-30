#####################################################
#### File with functions to plot the coordinates ####
#####################################################
import variables as v
import numpy as np
import matplotlib.pyplot as plt
import os
import seaborn as sns

# official color-code
colour_Neu5Ac = (0.65, 0.26, 0.6)


def plot_single(plot_name, coordinates, path='', coords_add=np.array([[np.nan, np.nan]]), FOV_x=v.FOV_x, FOV_y=v.FOV_y, save=False):
    """ plots coordinates as sialic acids in a coordinate system with the heading of plot_name

        Parameters
        -----------
        plot_name:      String
                        heading of the plot and name under which it will be saved
        coordinates:    2D numpy-array
                        pairs of coordinates in an array, that should be plotted
        path:           String
                        path to which the final plot should be saved if save = True
        coords_add:     2D numpy-array
                        pairs of coordinates in an array that should be plotted additional for sanity check
        FOV_x:          float
                        x-dimension of plotted FOV
        FOV_y:          float
                        y-dimension of plotted FOV
        save:           boolean
                        if true, the plot is saved

        Returns
        ----------
        None.
        Plots the coordinates and saves the plot.
    """
    fig, axis = plt.subplots(figsize=((FOV_x/FOV_y)*8, (FOV_x/FOV_y)*8))
    axis.plot(coordinates[:, 0], coordinates[:, 1], "D", markersize=12, color=colour_Neu5Ac)      # plot coordinates with a red diamond
    axis.plot(coords_add[:, 0], coords_add[:, 1], '.', markersize=12, color="black")
    axis.set(xlabel=r'x-axis in [$\mu \mathrm{m}$]', ylabel=r'y-axis in [$\mu \mathrm{m}$]')     # set axis labels
    axis.set_title(plot_name)  # set a title
    axis.set_xlim([0, FOV_x])   # set limits to the FOV
    axis.set_ylim([0, FOV_y])
    plt.tight_layout()                                                                         # use the tight layout
    if save:
        if not os.path.exists(path + '\\' + plot_name + '.pdf'):
            plt.savefig(path + '\\' + plot_name + '.pdf', bbox_inches='tight')      # save the figure
            plt.close()
        else:
            print(fr"ERROR: {path}\{plot_name}.pdf was not saved, since it already existed!!!")
    else:
        plt.show()   # show the plot


def plot_single_3d(plot_name, coordinates, path='', coords_add=np.array([[np.nan, np.nan, np.nan]]), FOV_x=v.FOV_y, FOV_y=v.FOV_y, save=False):
    """ plots coordinates as sialic acids in a coordinate system with the heading of plot_name

        Parameters
        -----------
        plot_name:      String
                        heading of the plot and name under which it will be saved
        coordinates:    2D numpy-array
                        tupel of coordinates(x, y, z) in an array, that should be plotted
        path:           String
                        path to which the final plot should be saved if save = True
        coords_add:     2D numpy-array
                        tupel of coordinates in an array that should be plotted additional for sanity check
        FOV_x:          float
                        x-dimension of plotted FOV
        FOV_y:          float
                        y-dimension of plotted FOV
        save:           boolean
                        if true, the plot is saved

        Returns
        ----------
        None.
        Plots the coordinates and saves the plot.
    """
    fig, axis = plt.subplots(figsize=(10, 8))

    z_min = 0
    z_max = max(np.max(coordinates[:, 2]), np.max(coords_add[:, 2]))

    # Scatter plot with color representing the third dimension
    sc1 = axis.scatter(coordinates[:, 0], coordinates[:, 1], c=coordinates[:, 2], marker='D', s=150, cmap='viridis', vmin=z_min, vmax=z_max)

    # Additional coordinates for sanity check
    axis.scatter(coords_add[:, 0], coords_add[:, 1], c=coords_add[:, 2], marker='.', s=150, cmap='viridis', vmin=z_min, vmax=z_max)

    plt.colorbar(sc1, ax=axis, label=r'z-axis in [$\mu \mathrm{m}$]', pad=0.1)

    axis.set_title(plot_name)  # set a title
    axis.set_xlim([0, FOV_x])   # set limits to the FOV
    axis.set_ylim([0, FOV_y])
    plt.tight_layout()                                                                         # use the tight layout
    if save:
        if not os.path.exists(path + '\\' + plot_name + '.pdf'):
            plt.savefig(path + '\\' + plot_name + '.pdf', bbox_inches='tight')      # save the figure
            plt.close()
        else:
            print(fr"ERROR: {path}\{plot_name}.pdf was not saved, since it already existed!!!")
    else:
        plt.show()   # show the plot


def plot_scatter_3d(plot_name, coordinates, coords_add=np.array([np.nan, np.nan, np.nan]), FOV_x=v.FOV_x,
                    FOV_y=v.FOV_y, path='', save=False):
    """ plots coordinates as sialic acids in a coordinate system with the heading of plot_name

        Parameters
        -----------
        plot_name:      String
                        heading of the plot and name under which it will be saved
        coordinates:    2D numpy-array
                        tupel of coordinates(x, y, z) in an array, that should be plotted
        coords_add:     2D numpy-array
                        tupel of added coordinates(x, y, z) in an array, that should be plotted
        FOV_x:          float
                        x-dimension of plotted FOV
        FOV_y:          float
                        y-dimension of plotted FOV
        path:           String
                        path to which the final plot should be saved if save = True
        coords_add:     2D numpy-array
                        tupel of coordinates in an array that should be plotted additional for sanity check
        save:           boolean
                        if true, the plot is saved

        Returns
        ----------
        None.
        Plots the coordinates and saves the plot.
    """
    fig = plt.figure()
    axis = fig.add_subplot(projection='3d')
    axis.scatter(coordinates[:, 0], coordinates[:, 1], coordinates[:, 2])
    axis.scatter(coords_add[:, 0], coords_add[:, 1], coords_add[:, 2])

    axis.set_title(plot_name)  # set a title
    axis.set_xlim([0, FOV_x])   # set limits to the FOV
    axis.set_ylim([0, FOV_y])
    plt.tight_layout()                                                                         # use the tight layout
    if save:
        if not os.path.exists(path + '\\' + plot_name + '_new.pdf'):
            plt.savefig(path + '\\' + plot_name + '_new.pdf', bbox_inches='tight')      # save the figure
            plt.close()
        else:
            print(fr"ERROR: {path}\{plot_name}.pdf was not saved, since it already existed!!!")
    else:
        plt.show()   # show the plot


def plot_labeling_nN_violin(plot_name, nN, labeling, path='', save=False, ylim=0):
    """ plots violin plot of mean nN distances over labeling efficiencies.

            Parameters
            -----------
            plot_name:      String
                            heading of the plot and name under which it will be saved
            nN:             2D numpy-array
                            with nN distances
            labeling:       numpy array
                            array containing the labeling efficiencies
            path:           String
                            path to which the final plot should be saved if save = True
            save:           Boolean
                            if true, the plot is saved
            ylim:           float
                            If unequal to 0, determines the upper y limit of the plot

            Returns
            ----------
            None.
            Plots mean nN distance over labeling efficiency and saves the plot.
        """
    fig, axis = plt.subplots()
    sns.violinplot(data=nN, cut=0, density_norm='count', color='lightgrey')
    axis.set_xticks(range(len(labeling)))
    axis.set_xticklabels(labeling)
    axis.set(xlabel=r'labeling efficiency', ylabel=r'nearest neighbor distances in $\mu \mathrm{m}$')  # set axis labels
    if ylim != 0:
        axis.set_ylim(0-1/30*ylim, ylim)
    axis.set_title(plot_name)  # set a title
    plt.tight_layout()                                                                             # use tight layout
    if save:
        if not os.path.exists(path + rf'\nN_violin_' + plot_name + '.pdf'):  # check whether the path already exists
            plt.savefig(path + rf'\nN_violin_' + plot_name + '.pdf', bbox_inches='tight')    # save the figure
            plt.close()
        else:
            print(fr"ERROR: {path}\nN_violin_{plot_name}.pdf was not saved, since it already existed!!!")
    else:
        plt.show()      # show the plot


