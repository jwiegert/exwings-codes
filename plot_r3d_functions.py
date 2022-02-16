# Functions for plotting various input and output data for Radmc3d
# ------------------------------------------------------------ #
# Useful packages
import os
import csv
import numpy as np

import matplotlib.pyplot as plt

# Basic definitions
AUcm = 1.49598e13

# ------------------------------------------------------------ #

def set_figurefonts():
    # Import required packages
    from matplotlib import rc
    rc('font',**{'family':'serif','serif':['serif']})
    rc('text', usetex=True)

    # Set size of ticks and tick fonts, labels and scales.
    rc('xtick.major',size=8)
    rc('xtick.minor',size=4)
    rc('ytick.major',size=8)
    rc('ytick.minor',size=4)

    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)


# This function sets standard plot settings
def set_ticklabels(xlabel,ylabel,xscale,yscale):

    # Import required packages
    from matplotlib import rc
    rc('font',**{'family':'serif','serif':['serif']})
    rc('text', usetex=True)

    # Set size of ticks and tick fonts, labels and scales.
    rc('xtick.major',size=8)
    rc('xtick.minor',size=4)
    rc('ytick.major',size=8)
    rc('ytick.minor',size=4)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.xlabel(xlabel,fontsize=18)
    plt.ylabel(ylabel,fontsize=18)
    if xscale == "log":
        plt.xscale('log')
    if yscale == "log":
        plt.yscale('log')

# This function plots the specified grid
def plot_grid(
        gridpath:str='../grid_distances.csv',
        amrpath:str='../amr_grid.inp',
        nbins:int=10
    ):
    """
    Loads and plots the current grid.

    Inputs
    ------
    gridpath: path to grid_distances.csv
    amrpath: path to amr_grid.inp
    nbins: number of bins for radial histogram plot
    """

    # Check if both files exist
    if os.path.exists(gridpath) == True and os.path.exists(amrpath) == True:

        # Extract necessary info from amr_grid
        with open(amrpath, 'r') as f:
            for nn,line in enumerate(f.readlines()):

                # nleafs is the number of cells in the grid
                # This is supposed to be the second number on line ten
                if nn == 9:
                    nleafs = int(line.split(' ')[1])

        # Create griddistances array
        griddistances = np.zeros((nleafs,4))


        # Load distances to cells
        with open(gridpath, 'r') as f:
            griddistancesfile = csv.reader(f)

            # Set index o 0
            nn = 0
            for row in griddistancesfile:

                # Skip the header and add radial, x,y,z distances to array
                if row[0][-1].isdigit() == True:
                    
                    # Save radial, x, y, z distances (in this order)
                    griddistances[nn,:] = [float(row[0])/AUcm,float(row[1])/AUcm,float(row[2])/AUcm,float(row[3])/AUcm]
                    # Increase index
                    nn += 1

        # Plots (need a nice way to set fontsizes and fonts?)
        fig, ax = plt.subplots(2,2)

        # Plot radial distances
        ax[0,0].hist(griddistances[:,0],bins=nbins)
        ax[0,0].set(
            title='Radial distances', 
            xlabel='Number of cells', 
            ylabel='Radial to centrum of grid (AU)'
        )

        # Plot coordinates in each plane
        ax[0,1].plot(griddistances[:,1],griddistances[:,2],'b.')
        ax[0,1].set(
            title='Cells in X-Y-plane', 
            xlabel='X coord (AU)', 
            ylabel='Y coord (AU)'
        )

        ax[1,0].plot(griddistances[:,1],griddistances[:,3],'b.')
        ax[1,0].set(
            title='Cells in X-Z-plane', 
            xlabel='X coord (AU)', 
            ylabel='Z coord (AU)'
        )

        ax[1,1].plot(griddistances[:,2],griddistances[:,3],'b.')
        ax[1,1].set(
            title='Cells in Y-Z-plane', 
            xlabel='Y coord (AU)', 
            ylabel='Z coord (AU)'
        )

        # Better spacing between figures
        fig.tight_layout()
        fig.show()

    else:
        print(f'ERROR: plot_grid can not find {gridpath} and/or {amrpath}.')
