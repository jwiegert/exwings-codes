# Functions for plotting various input and output data for Radmc3d
# ------------------------------------------------------------ #
# Useful packages


import matplotlib.pyplot as plt

import create_r3d_functions as c3d

import numpy as np

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


# This function sets standard plot settings TODO: REMOVE THIS?
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

# ------------------------------------------------------------ #

# This function plots the specified grid
def plot_grid(
        gridpath:str='../grid_distances.csv',
        sizepath:str='../grid_cellsizes.csv',
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
    AUcm = 1.49598e13

    # Load data
    griddistances = c3d.load_griddistances(gridpath,amrpath)
    gridsizes = c3d.load_cellsizes(sizepath,amrpath)
    
    # Load some grid props
    nleafs = c3d.load_gridprops()[2]
    ncellsizes = np.size(np.unique(gridsizes))

    # Change units to AU
    for nn in range(nleafs):
        gridsizes[nn] = gridsizes[nn]/AUcm
        for nx in range(4):
            griddistances[nn,nx] = griddistances[nn,nx]/AUcm

    # Plots (need a nice way to set fontsizes and fonts?)
    fig, ax = plt.subplots(3,2)

    # Plot radial distances
    ax[0,0].hist(griddistances[:,0],bins=nbins)
    ax[0,0].set(
        title='Radial distances', 
        ylabel='Number of cells', 
        xlabel='Radial to centrum of grid (AU)'
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

    ax[2,0].hist(gridsizes,bins=ncellsizes)
    ax[2,0].set(
        title='Cell sizes',
        xlabel='Cell size (AU)',
        ylabel='Number of cells'
    )

    # Better spacing between figures
    fig.tight_layout()
    fig.show()

# Plot the densities as function of radial distance from centrum of star
# def plot_radialdensity()

