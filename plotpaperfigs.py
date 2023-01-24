# Plots various figures for first co5bold-r3d-paper
import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np

import analyze_co5bold_functions as a5d
import analyze_r3d_functions as a3d


# Figure settings
rc('font',**{'family':'serif','serif':['serif']})
rc('text', usetex=True)
rc('xtick.major',size=8)
rc('xtick.minor',size=4)
rc('ytick.major',size=8)
rc('ytick.minor',size=4)

# Constants
AUcm = 1.49598e13 # cm


# Plot choices
plot_coboldgrid = 'n'
plot_opticalthickness = 'n'
plot_grainsizehist_all = 'n'
plot_grainsizehist_one = 'n'
plot_absscat = 'y'

# ----------------------------------------------------------------
# FIG Cut through of CO5BOLD grid for st28gm06n052 with cell 
# sizes against distance from centre of the grid.

if plot_coboldgrid == 'y':

    c5dgrid,cellcourners,cellsize = a5d.load_grid_properties(
        savpath='../co5bold_data/dst28gm06n052/st28gm06n052_186.sav'
    )
    Ncourners = cellcourners[:,0].size - 1
    cellsizes = np.zeros(Ncourners)
    cellcoords = np.zeros(Ncourners)

    for nn in range(Ncourners):
        cellsizes[nn] = (cellcourners[nn+1,0] - cellcourners[nn,0])/AUcm
        cellcoords[nn] = (0.5*(cellcourners[nn+1,0] + cellcourners[nn,0]))/AUcm


    fig, ax = plt.figure('grid-cutthrough', figsize=(6, 4)), plt.axes()

    ax.plot(cellcoords,cellsizes,'.',markersize=2)

    ax.set(ylim=[0,0.14])
    ax.set_xlabel(r'Distance from centre along one axis (au)',fontsize=18)
    ax.set_ylabel(r'Cell size (au)',fontsize=18,)
    ax.tick_params(axis='both', which='major', labelsize=15)

    fig.tight_layout()
    fig.savefig('figs/grid-cutthrough.pdf', facecolor="white")
    fig.show()

# ----------------------------------------------------------------
#
# Plot Figure with average LOS optical thickness of gasmodel

if plot_opticalthickness == 'y':

    Reffective,fig,ax = a3d.plot_opticalthick(
        path='../r3dresults/st28gm06n056/140/'
    )
    print(f'{Reffective / 1.49598e13} AU')

    ax.set_xlabel(r'Distance along LOS (AU)',fontsize=18)
    ax.set_ylabel(r'Optical thickness, $\tau$',fontsize=18)
    ax.tick_params(axis='both', which='major', labelsize=15)

    fig.tight_layout()
    fig.savefig('figs/optthick_los.pdf', facecolor="white")
    fig.show()

# ----------------------------------------------------------------
#
# Plot Figure with Mass-vs-grain size histogram
#
# Warning: slow because of loading lots dust densities

plot_grainsizehist_all = 'n'
if plot_grainsizehist_all == 'y':

    fig,ax = a5d.plot_grainsizemass_histogram(
        model='st28gm06n052',
        phases=[186,190,198]
    )

    ax[0].set_ylabel(r'Dust mass (g)',fontsize=18)    
    ax[1].set_ylabel(r'Dust mass (g)',fontsize=18)
    ax[2].set_ylabel(r'Dust mass (g)',fontsize=18)

    ax[2].set_xlabel(r'Grain size ($\mu$m)',fontsize=18)

    ax[0].tick_params(axis='both', which='major', labelsize=15)
    ax[1].tick_params(axis='both', which='major', labelsize=15)
    ax[2].tick_params(axis='both', which='major', labelsize=15)

    fig.tight_layout()
    fig.savefig('figs/grainsize_hist.pdf', facecolor="white")

    fig.show()

if plot_grainsizehist_one == 'y':

    fig,ax = a5d.plot_grainsizemass_histogram(
        model='st28gm06n052',
        phases=[186,190,198]
    )

    ax.set_ylabel(r'Dust mass (g)',fontsize=18)    
    ax.set_xlabel(r'Grain size ($\mu$m)',fontsize=18)
    ax.tick_params(axis='both', which='major', labelsize=15)

    fig.tight_layout()
    fig.savefig('figs/grainsize_hist.pdf', facecolor="white")

    fig.show()

# ----------------------------------------------------------------
#
# Plot Figure absorption and scattering, and angles

if plot_absscat == 'y':

    fig,ax = a3d.plot_allkappa(
        path='../r3dresults/opacities_st28gm06n052/'
    )

    ax[0].set_ylabel(r'$\kappa_{\rm abs}$ (cm$^2$/g)',fontsize=18)    
    ax[1].set_ylabel(r'$\kappa_{\rm scat}$ (cm$^2$/g)',fontsize=18)
    ax[2].set_ylabel(r'$\kappa_{\rm abs} + \kappa_{\rm scat}$ (cm$^2$/g)',fontsize=18)
    ax[3].set_ylabel(r'$\left< \cos \theta \right>$',fontsize=18)

    ax[3].set_xlabel(r'Wavelength ($\mu$m)',fontsize=18)

    ax[0].tick_params(axis='both', which='major', labelsize=15)
    ax[1].tick_params(axis='both', which='major', labelsize=15)
    ax[2].tick_params(axis='both', which='major', labelsize=15)
    ax[3].tick_params(axis='both', which='major', labelsize=15)

    fig.tight_layout()
    fig.savefig('figs/abs_scat_angle.pdf', facecolor="white")

    fig.show()
