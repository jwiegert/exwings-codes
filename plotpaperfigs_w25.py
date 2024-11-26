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
import analyze_timedep_functions as atf

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

# Set LOS-angles
angles = [
    'i000_phi000',
    'i090_phi000',
    'i090_phi090',
    'i090_phi270',
    'i180_phi000',
    'i270_phi000',
]



# Plot-list

plot_dustmass = 'n'
plot_075grainsize = 'n'

plot_allseds = 'n'
plot_luminosities = 'y'


# Plots below ----------------------------------------------------------------#
#
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
    legendlist = [
        'st28gm06n052','st28gm06n074','st28gm06n075'
    ]
    # Plot all dust masses on top of eachother
    ax[0].plot(time052,dustmass052,'b', label=legendlist[0])
    ax[0].plot(time074,dustmass074,'r', label=legendlist[1])
    ax[0].plot(time075,dustmass075,'g', label=legendlist[2])

    # And 075 in its own plot
    ax[1].plot(time075,dustmass075,'g', label=legendlist[2])

    ax[0].set_ylabel(r'Dust mass ($M_\odot $)',fontsize=18)
    ax[1].set_xlabel('Time (yrs)',fontsize=18)
    ax[0].tick_params(axis='both', which='major', labelsize=15)
    ax[1].tick_params(axis='both', which='major', labelsize=15)
    ax[0].set_xlim([time052[0],time052[-1]])
    ax[1].set_xlim([time052[0],time052[-1]])
    ax[0].legend(fontsize=14)
    ax[1].legend(fontsize=14)

    # Save figure
    fig.tight_layout()
    ax[0].yaxis.set_label_coords(-0.09,0.09) # Moves ylabel to centre of vertical
    fig.savefig(
        f'figs/052_074_075_dustmasscompare.pdf', dpi=300, facecolor="white"
    )
    fig.show()
#
#####################################################################################
# Plot max grain size bin vs time of 075 
if plot_075grainsize == 'y':
    # Load data from *075 model
    path075 = '../r3dresults/st28gm06n075_nospikes/'
    phasetimes075 = np.loadtxt(path075+'snapshot_yr.dat')[:,1]
    grainsize075 = np.loadtxt(path075+'maxgrainsize.dat')[:,1]

    # Extract some numbers
    print(f'Max grain size. Median: {np.median(grainsize075)}. Average {np.mean(grainsize075)}')

    # Plot and save figure
    fig,ax = plt.figure(figsize=(6,3)), plt.axes()
    ax.plot(phasetimes075,grainsize075)
    ax.set_xlabel('Time (yrs)',fontsize=18)
    ax.set_ylabel(r'Max grain size ($\mu$m)',fontsize=18)
    ax.tick_params(axis='both', which='major', labelsize=15)
    fig.tight_layout()
    fig.savefig('figs/075_maxgrainsize.pdf', dpi=300)
#
#####################################################################################
# Plot all SEDs in three subplots 
if plot_allseds == 'y':
    
    # Set figure objects    
    fig,ax = plt.subplots(
        1,3,
        figsize=(12,4.5)
    )
    # Set paths and loop through them
    paths = [
        '../r3dresults/st28gm06n052_timedep_nospikes/',
        '../r3dresults/st28gm06n074_nospikes/',
        '../r3dresults/st28gm06n075_nospikes/'
    ]
    for nmodel,path in enumerate(paths):
        # Extract all available phase numbers
        phases = [int(filename) for filename in os.listdir(path) if os.path.isdir(path+filename)]
        phases.sort()

        # Extract average, std and minmax SEDs
        flux_average,flux_std,flux_min,flux_max = atf.extract_averageseds(
            path=path,
            phases=phases,
            angles=angles,
            save_datafile='n',
        )
        # Plot all
        ax[nmodel].set(
            xscale='log',
            yscale='log',
            xlim=[0.5,40],
            ylim=[1e5,2e8]
        )
        ax[nmodel].tick_params(axis='both', which='major', labelsize=15)
        if nmodel == 0:
            linecolour = 'darkblue'
            allcolour = 'deepskyblue'
        if nmodel == 1:
            linecolour = 'darkred'
            allcolour = 'salmon' #lr tomato
        if nmodel == 2:
            linecolour = 'darkgreen'
            allcolour = 'lawngreen'

        for phase in phases:
            for angle in angles:

                # Loop over all SEDs
                wavelength, fluxes = a3d.load_spectrum(
                    path=f'{path}{phase}/spectrum_{angle}.out'
                )
                # Plot all SEDs as a lighter shade
                ax[nmodel].plot(
                    wavelength,fluxes,
                    color=allcolour
                )
        
        # Plot average, stdrange of all
        ax[nmodel].plot(
            wavelength,flux_average,
            color=linecolour, linewidth=2
        )
        ax[nmodel].plot(
            wavelength,flux_average-flux_std,
            color=linecolour,linestyle='--', linewidth=2
        )
        ax[nmodel].plot(
            wavelength,flux_average+flux_std,
            color=linecolour,linestyle='--', linewidth=2
        )

    # Set final axis-settings
    ax[0].set_ylabel(r'$F_\nu$ (Jy at 1 pc)', fontsize=18)
    ax[1].set_xlabel(r'Wavelength ($\mu$m)',fontsize=18)
    ax[0].set_title('st28gm06n052', fontsize=15)
    ax[1].set_title('st28gm06n074', fontsize=15)
    ax[2].set_title('st28gm06n075', fontsize=15)

    # Save figure
    fig.tight_layout()
    fig.savefig(
        'figs/all_seds.pdf', 
        facecolor='white',
        dpi=300
    )
    fig.show()
#
#####################################################################################
# Plot all luminosities in three subplots 
if plot_luminosities == 'y':

    # TODO
    #  flatten at bit
    # large tick-fonts



    # Set figure objects
    fig,ax = plt.subplots(
        1,3,
        figsize=(12,4)
    )
    # Set paths and loop through them
    paths = [
        '../r3dresults/st28gm06n052_timedep_nospikes/',
        '../r3dresults/st28gm06n074_nospikes/',
        '../r3dresults/st28gm06n075_nospikes/'
    ]
    paths_nodust = [
        '../r3dresults/st28gm06n052_timedep_nodust/',
        '../r3dresults/st28gm06n074_nodust/',
        '../r3dresults/st28gm06n075_nodust/'
    ]
    models = [
        'st28gm06n052',
        'st28gm06n074',
        'st28gm06n075'
    ]
    Nangles = 6

    for nmodel,path in enumerate(paths):

        # Load time axis and lum-data
        phasetimes = np.loadtxt(path+'snapshot_yr.dat')
        luminosities = np.loadtxt(path+'luminosity.dat')
        luminosity_nodust = np.loadtxt(paths_nodust[nmodel]+'luminosity.dat')[:,-1]

        # Set colours and tick font size
        if nmodel == 0:
            linecolour = 'darkblue'
            allcolour = 'skyblue'
        if nmodel == 1:
            linecolour = 'darkred'
            allcolour = 'mistyrose'
        if nmodel == 2:
            linecolour = 'darkgreen'
            allcolour = 'palegreen'
        ax[nmodel].set_title(models[nmodel])
        ax[nmodel].tick_params(axis='both', which='major', labelsize=15)

        # Plot all angles with dust
        for nangle in range(Nangles):
            ax[nmodel].plot(
                phasetimes[:,1],luminosities[:,1+nangle],
                color=allcolour
            )
        # Plot average with dust
        ax[nmodel].plot(
            phasetimes[:,1],luminosities[:,-1],
            color=linecolour,
            linewidth=2
        )
        # Plot average without dust
        ax[nmodel].plot(
            phasetimes[:,1],luminosity_nodust,
            linestyle='--',
            color='black',
        )

    # Set final axis-settings
    ax[0].set_ylabel(r'$L$ ($L_\odot$)', fontsize=18)
    ax[1].set_xlabel(r'Simulation time (yrs)',fontsize=18)

    fig.tight_layout()
    fig.savefig(
        'figs/luminosities.pdf', 
        facecolor='white',
        dpi=300
    )