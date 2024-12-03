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
angles_label = [
    '0-0',
    '90-0',
    '90-90',
    '90-270',
    '180-0',
    '270-0'
]
# Set a list of model labels
models_label = [
    'st28gm06n052',
    'st28gm06n074',
    'st28gm06n075',
]



# Plot-list

plot_dustmass = 'n'
plot_075grainsize = 'n'

plot_allseds = 'n'
plot_luminosities = 'n'

plot_rsourceevents = 'y'



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
    legendlist = models_label
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
    models = models_label
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
#
#####################################################################################
# Plot average period of Rsource events per model and angle
if plot_rsourceevents == 'y':

    Nmodels = len(models_label)
    Nangles = len(angles)


    # Set figure objects    
    fig,ax = plt.figure(
        figsize=(6,4)
    ), plt.axes()

    ax.tick_params(axis='both', which='major', labelsize=15)
    ax.set_ylim([0,23])    
    ax.set_xlim([-0.5,5.5])
    # Combine each style to one set for each model, plot in one figure
    model052 = [
        [
            6.971389087364163,
            5.809490906136802,
            4.979563633831545,
            5.809490906136802,
            5.809490906136802,
            5.809490906136802
        ],[
            0.95050,
            1.7690,
            3.3041,
            1.7426,
            2.5875,
            1.3729,
        ]
    ]
    model074 = [
        [
            -1,
            4.9795575666192935,
            -1,
            5.809483827722509,
            5.809483827722509,
            11.618967655445019,
        ],[
            -1,
            0.47541,
            -1,
            1.5319,
            2.1658,
            2.3242,
        ]
    ]
    model075 = [
        [
            8.71422145903946,
            11.618961945385948,
            34.85688583615784,
            17.42844291807892,
            6.9713771672315685 ,
            2.4897775597255603,
        ],[
            0.43570,
            0.15844,
            0.158437231328179,
            0.15844,
            0.22181,
            0.23766,
        ]
    ]
    # Plot 052
    # x-numbering must be -1 to get correct labels....
    for nangle in range(Nangles):
        ax.plot(nangle-0.05,model052[0][nangle],'bo',markersize=6)
        ax.plot(
            [nangle-0.05,nangle-0.05],
            [
                model052[0][nangle]-0.5*model052[1][nangle],
                model052[0][nangle]+0.5*model052[1][nangle]
            ]
            ,'b'
        )
        # Plot 074
        ax.plot(nangle,model074[0][nangle],'rs',markersize=6)
        ax.plot(
            [nangle,nangle],
            [
                model074[0][nangle]*model074[1][nangle],
                model074[0][nangle]*model074[1][nangle]
            ]
            ,'r'
        )


        # Plot 075
        ax.plot(nangle+0.05,model075[0][nangle],'gd',markersize=6)
        ax.plot(
            [nangle+0.05,nangle+0.05],
            [
                model075[0][nangle]+0.5*model075[1][nangle],
                model075[0][nangle]+0.5*model075[1][nangle]
            ]
            ,'g'
        )
    # List the labels so that theres 1 per model.
    ax.plot(-1,-1,'bo',markersize=6,label=models_label[0])
    ax.plot(-1,-1,'rs',markersize=6,label=models_label[1])
    ax.plot(-1,-1,'gd',markersize=6,label=models_label[2])
    ax.legend(
        loc='upper left',
        fontsize=14
    )
    ax.set_title(
        r'$R_{\rm source}$',
        fontsize=18
    )
    # Modify xticklabels
    ax.set_xticks([0,1,2,3,4,5]) 
    ax.set_xticklabels(
        angles_label
    ) 
    # Set axislabels
    ax.set_ylabel(r'Av. period \& event length (yrs)', fontsize=18)
    ax.set_xlabel(r'LOS-angle',fontsize=18)


    # TODO
    # ha med F2 och F10-statistiken i subplots under här??




    # ylabel: average period & event length (yr)
    # xlabel: LOS-angle

    fig.tight_layout()
    fig.savefig(
        'figs/periods_rsource.pdf', 
        facecolor='white',
        dpi=300
    )

