# Plots various figures for first co5bold-r3d-paper
import matplotlib.pyplot as plt
from scipy import ndimage
import numpy as np
import scipy.ndimage
import re
import os

from matplotlib import rc
from mpl_toolkits.axes_grid1 import make_axes_locatable

import analyze_co5bold_functions as a5d
import analyze_r3d_functions as a3d
import create_r3d_functions as c3d

# Figure settings
rc('font',**{'family':'serif','serif':['serif']})
rc('text', usetex=True)
rc('xtick.major',size=8)
rc('xtick.minor',size=4)
rc('ytick.major',size=8)
rc('ytick.minor',size=4)

# Constants
c = 2.998e8 # m/s
AUcm = 1.49598e13 # cm
Lsol = 3.828e26 # W
Rsolmeter = 6.955e8 # m
cubesize = 222757675648155.62/AUcm # Yes, hardcoded to large grid, change if needed for 
                                   # other grids, ie, from amr_grid, first coordinate 
                                   # is courner coordinate of first base cell
radian = 206264800 # milliasec
baselineVLTI = 201.92 # metres
baselineChandra = 300 # TODO metres

# Plot-list

plot_dustmass = 'y'
plot_075grainsize = 'n'





# Plot dust mass vs time for all three model
if plot_dustmass == 'y':
    # Load data from *052 model
    path052 = '../r3dresults/st28gm06n052_timedep_nospikes/'
    phasetimes052 = np.loadtxt(path052+'snapshot_yr.dat')[:,1]
    dustmass052 = np.loadtxt(path052+'dustmass.dat')[:,1]
    time052 = phasetimes052 - phasetimes052[0]

    # Load data from *074 model
    path074 = '../r3dresults/st28gm06n074_nospikes/'
    phasetimes074 = np.loadtxt(path074+'snapshot_yr.dat')[:,1]
    dustmass074 = np.loadtxt(path074+'dustmass.dat')[:,1]
    time074 = phasetimes074 - phasetimes074[0]

    # Load data from *075 model
    path075 = '../r3dresults/st28gm06n075_nospikes/'
    phasetimes075 = np.loadtxt(path075+'snapshot_yr.dat')[:,1]
    dustmass075 = np.loadtxt(path075+'dustmass.dat')[:,1]
    time075 = phasetimes075 - phasetimes075[0]

    # Iniatiate figure object
    fig,ax = plt.subplots(
        2,1, 
        figsize=(6,6),
        gridspec_kw={'height_ratios': [2, 1]}
    )

    # Plot all dust masses on top of eachother
    ax[0].plot(time075,dustmass075,'g')
    ax[0].plot(time074,dustmass074,'r')
    ax[0].plot(time052,dustmass052,'b')
    # And 075 in its own plot
    ax[1].plot(time075,dustmass075,'g')

    ax[0].set_ylabel(r'Dust mass ($M_\odot $)',fontsize=18)
    ax[1].set_xlabel('Time (yrs)',fontsize=18)
    ax[0].tick_params(axis='both', which='major', labelsize=15)
    ax[1].tick_params(axis='both', which='major', labelsize=15)
    ax[0].set_xlim([time052[0],time052[-1]])
    ax[1].set_xlim([time052[0],time052[-1]])

    # TODO
    # add labels in a box
    # save figure and put in paper



    # Save figure
    fig.tight_layout()
    fig.show()





# Plot max grain size bin vs time of 075
if plot_075grainsize == 'y':
    # Load data from *075 model
    path075 = '../r3dresults/st28gm06n075_nospikes/'
    phasetimes075 = np.loadtxt(path075+'snapshot_yr.dat')[:,1]
    grainsize075 = np.loadtxt(path075+'maxgrainsize.dat')[:,1]

    # Extract some numbers
    print(f'Max grain size. Median: {np.median(grainsize075)}. Average {np.mean(grainsize075)}')

    # Plot and save figure
    fig,ax = plt.figure(figsize=(6,4)), plt.axes()
    ax.plot(phasetimes075,grainsize075)
    ax.set_xlabel('Time (yrs)',fontsize=18)
    ax.set_ylabel(r'Max grain size ($\mu$m)',fontsize=18)
    ax.tick_params(axis='both', which='major', labelsize=15)
    fig.tight_layout()
    fig.savefig('../r3dplots/075_maxgrainsize.pdf', dpi=300)


