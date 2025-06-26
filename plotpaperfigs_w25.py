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
plot_052fluxdensity = 'n'
plot_052exampleimages = 'y'

plot_LOSevents = 'n'
plot_fluxvariations = 'n'
plot_datacompare = 'n'

# For vr-prop
plotvr_exampleimages = 'n'
plotvr_radiusplot = 'n'
plotvr_datacompare = 'n'


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
    ax[1].set_xlabel('Synchronised time (yrs)',fontsize=18)
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
    ax.set_xlabel('Simulation time (yrs)',fontsize=18)
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
    fig,ax2d = plt.subplots(
        2,2,
        figsize=(10,7)
    )
    ax = ax2d.ravel()
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
    Nangles = len(angles)

    # List all luminosities of these data
    #    052_nodust - average lum: 7623.889592760182 Lsol
    #    074_nodust - average lum: 7632.387339366515 Lsol
    #    075_nodust - average lum: 7603.132628959276 Lsol
    #
    #
    #    052_nospikes - max-max lum: 11216.202 Lsol
    #    074_nospikes - max-max lum: 11017.651 Lsol
    #    075_nospikes - max-max lum: 10948.281 Lsol
    #
    #    052_nospikes - average max lum: 7948.66 Lsol
    #    074_nospikes - average max lum: 8237.554 Lsol
    #    075_nospikes - average max lum: 9101.957 Lsol
    #
    #    052_nospikes - average lum: 5325.214027149322 Lsol
    #    074_nospikes - average lum: 6146.020846153845 Lsol
    #    075_nospikes - average lum: 7314.732800904978 Lsol
    #
    #    052_nospikes - average min lum: 3217.095 Lsol
    #    074_nospikes - average min lum: 4403.186 Lsol
    #    075_nospikes - average min lum: 5911.108 Lsol
    #
    #    052_nospikes - min-min lum: 583.556 Lsol
    #    074_nospikes - min-min lum: 766.543 Lsol
    #    075_nospikes - min-min lum: 3815.841 Lsol
    #
    # Each sublist is for each model, then I add in them in the same order
    # as above.
    minmaxlums = [
        [7623.890, 11216.202, 7948.660, 5325.214, 3217.095,  583.556],
        [7632.387, 11017.651, 8237.554, 6146.021, 4403.186,  766.543],
        [7603.133, 10948.281, 9101.957, 7314.733, 5911.108, 3815.841]
    ]

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
        ax[nmodel].set_ylim(0,11500)
    # Plot min-max-average luminosoties of each model
    ax[3].set_title(r'$L$ ranges')
    # max-max
    ax[3].fill_between(
        [0,1,2], 
        [minmaxlums[0][1], minmaxlums[1][1], minmaxlums[2][1]],
        [minmaxlums[0][5], minmaxlums[1][5], minmaxlums[2][5]],
        color='khaki'
    )
    # average-max
    ax[3].fill_between(
        [0,1,2], 
        [minmaxlums[0][2], minmaxlums[1][2], minmaxlums[2][2]],
        [minmaxlums[0][4], minmaxlums[1][4], minmaxlums[2][4]],
        color='gold',linestyle='-'
    )
    # average
    ax[3].plot(
        [0,1,2], 
        [minmaxlums[0][3], minmaxlums[1][3], minmaxlums[2][3]]
        ,color='darkgoldenrod',linestyle='-'
    )
    # Last: stars luminosity
    ax[3].plot(
        [0,1,2], 
        [minmaxlums[0][0], minmaxlums[1][0], minmaxlums[2][0]]
        ,'k--'
    )
    # Set xticklabels for last plot
    ax[3].set_xticks([0,1,2]) 
    ax[3].set_xticklabels(models) 
    ax[3].set_ylim(0,11500)
    ax[3].set_xlim(0,2)
    ax[3].tick_params(axis='both', which='major', labelsize=15)

    # Set final general axis-settings
    ax[0].set_ylabel(r'Luminosity ($L_\odot$)', fontsize=18)
    ax[2].set_ylabel(r'Luminosity ($L_\odot$)', fontsize=18)
    ax[2].set_xlabel(r'Simulation time (yrs)',fontsize=18)

    # Plot and save
    fig.tight_layout()
    fig.savefig(
        'figs/luminosities.pdf', 
        facecolor='white',
        dpi=300
    )
    fig.show()
#
####################################################################################
# Plot flux density at 2um over time for 052
#
if plot_052fluxdensity == 'y':

    # Wavelength to plot
    wavelength = 2

    # Path to model of choice
    path = '../r3dresults/st28gm06n052_timedep_nospikes/'
    model = '052'

    # Load all phases
    phases = [int(filename) for filename in os.listdir(path) if os.path.isdir(path+filename)]
    phases.sort()
    Nphases = len(phases)
    phasetimes = np.loadtxt(path+'snapshot_yr.dat')[:Nphases,1]

    # Load wavelengthgrid and extract index for wavelength
    wavelengths,spectrum = a3d.load_spectrum(
        path = f'{path}{phases[0]}/spectrum_{angles[0]}.out'
    )
    wavelengths = np.array(wavelengths)
    wavelengthindex = int(np.where(wavelengths >= wavelength)[0][0]-1)
    fluxall = np.zeros((len(angles),Nphases))

    # Load average Flux densities
    # extract flux at wavelength in MJy
    # and declare array for relative fluxes vs time
    average_seds = np.loadtxt(f'{path}average_sed.dat')
    average_flux = average_seds[wavelengthindex,1] * 1e-6
    relative_flux = np.zeros((len(angles),Nphases))

    # Set figure settings
    #fig, ax = plt.figure(num=f'F(t) of {wavelength} um', figsize=(6, 4)), plt.axes()
    fig,ax = plt.subplots(
        2,1,
        figsize=(6, 6)
    )
    ax[0].set_ylabel(rf'$F(2\,\mu$m$)$, MJy at 1\,pc',fontsize=18)
    ax[1].set_ylabel(r'$F(2\,\mu$m$) / F_{\rm average}$',fontsize=18)
    ax[0].set_yscale('log')
    ax[1].set_ylim(0,2)
    ax[1].set_xlabel(r'Simulation time (yrs)',fontsize=18)
    for nn in range(2):
        ax[nn].tick_params(axis='both', which='major', labelsize=15)
        ax[nn].set_xlim(phasetimes[0],phasetimes[-1])
    linecolour = 'darkblue'
    allcolour = 'deepskyblue'

    # Loop through all phases and extract all flux densities at wavelength
    # And the relative flux density (to average flux density)
    for nangle,angle in enumerate(angles):
        for nphase,phase in enumerate(phases):

            wavelengths,spectrum = a3d.load_spectrum(
                path = f'{path}{phase}/spectrum_{angle}.out'
            )
            fluxall[nangle,nphase] = spectrum[wavelengthindex]*1e-6

            relative_flux[nangle,nphase] = fluxall[nangle,nphase]/average_flux

        # Plot each angles time dependent F in MJy
        ax[0].plot(phasetimes,fluxall[nangle,:],allcolour)

        # Plot relative flux densities of each angle and dashed lines for the limits
        ax[1].plot(phasetimes,relative_flux[nangle,:],allcolour)
        ax[1].plot(
            [phasetimes[0],phasetimes[-1]],
            [1,1],
            'k-'
        )
        ax[1].plot(
            [phasetimes[0],phasetimes[-1]],
            [0.8,0.8],
            'k--'
        )
        ax[1].plot(
            [phasetimes[0],phasetimes[-1]],
            [0.6,0.6],
            'k:'
        )

    # Save average of each angles flux density at each time
    fluxaverage = []
    for nphase in range(Nphases):
        fluxaverage.append(
            np.mean(fluxall[:,nphase])
        )
    # and plot average flux density
    ax[0].plot(phasetimes,fluxaverage,linecolour)

    # and save figure (and show)
    fig.tight_layout()
    plt.savefig(f'figs/{model}_fluxtime_{wavelength}um.pdf', dpi=300)
    fig.show()
#
#####################################################################################
# Plot 6 snapshots at 2 and 10um to showcase what we're looking at
if plot_052exampleimages == 'y':

    # Set paths
    path = '../r3dresults/st28gm06n052_timedep_nospikes/'
    modelabbreviation = '052'

    # Define wavelenvths
    wavelengths = [
        '02','10'
    ]
    # Chose snapshots
    snapshots = [
        290, 292, 294, 296, 298
    ]
    snapshots_times = []

    # Extract corresponding snapshot-times
    snapshot_file = np.loadtxt(
        path+'snapshot_yr.dat'
    )
    for snaptime in snapshot_file:
        for snapshot in snapshots:
            if snapshot == snaptime[0]:
                snapshots_times.append(snaptime[1])

    # Create image objects and fill subplots
    fig,ax = plt.subplots(
        len(snapshots),len(wavelengths),
        figsize=(5, 17)
    )
    # Loop through wavelengths
    for ncolumn,wavelength in enumerate(wavelengths):
        imagefilename = f'image_i000_phi000_{wavelength}um.out'

        # Loop through snapshots
        for ntime,snapshot in enumerate(snapshots):

            # Load image
            image2d,image2dlog,flux,axisplot = a3d.load_images(
                path=f'{path}{snapshot}/',
                image=imagefilename,
                distance=1
            )
            # Change to MJy
            image2d = image2d/1e6
            # Set scales for the two wavelengths
            if ncolumn == 0:
                scale = [1e-2,7]
            if ncolumn == 1:
                scale = [1e-2,1]
            # Plot image and save colourbar info
            imbar = ax[ntime,ncolumn].imshow(
                image2d, 
                origin='lower', 
                extent=axisplot, 
                vmin=scale[0],
                vmax=scale[1],
                cmap=plt.get_cmap('hot')
            )
            # write time and flux on top of each left column
            if ncolumn == 0:
                ax[ntime,0].set_title(
                    f'{snapshots_times[ntime]:.2f} yrs, {flux*1e-6:.3f} MJy',
                    fontsize = 10
                )
            if ncolumn == 1:
                # and only flux on second column
                ax[ntime,1].set_title(
                    f'{flux*1e-6:.3f} MJy',
                    fontsize = 10
                )
                # Remove ylabelticks for second column except bottom
                ax[ntime,1].axes.yaxis.set_ticklabels([])
            # And xticklabels for all rows except bottom
            if ntime < len(snapshots)-1:
                ax[ntime,ncolumn].axes.xaxis.set_ticklabels([])

            # Offset AU on outer plots
            if ncolumn == 0:
                ax[ntime,ncolumn].set_ylabel('Offset (au)', fontsize = 14)
            # And colour bar on bottom bar
            if ntime == len(snapshots)-1:
                divider = make_axes_locatable(ax[ntime,ncolumn])
                cax = divider.append_axes(
                    'bottom', 
                    size='6%', 
                    pad=0.7,
                )
                cb0 = plt.colorbar(
                    imbar, 
                    cax=cax, 
                    orientation = 'horizontal', 
                )
                cb0.ax.tick_params(labelsize=12)
                # Only label on one, TODO move this to middle of plot
                if ncolumn == 0:
                    cb0.set_label(
                        label = 'Flux density (MJy at 1pc)', fontsize=12
                    )

    # NOTE
    # final layout of image is fixed manually in GIMP

    fig.tight_layout()
    fig.savefig(
        'figs/052exampleimages.png', 
        facecolor='white',
        dpi=300
    )
    #fig.show()

#
#####################################################################################
# Plot average period of Rsource, F2um, F10um events per model and angle
if plot_LOSevents == 'y':

    Nmodels = len(models_label)
    Nangles = len(angles)

    fig,ax = plt.subplots(
        1,3,
        figsize=(13,4.5)
    )
    for nmodel in range(Nmodels):
        ax[nmodel].tick_params(axis='both', which='major', labelsize=15)
        ax[nmodel].set_ylim([0,28])    
        ax[nmodel].set_xlim([-0.5,5.5])
    # Combine each style to one set for each model, plot in one figure
    #
    # model[0][:] - perioder per LOS
    # model[1][:] - medeleventlength per LOS
    #
    rsource052 = [
        [
            6.971389087364163,
            5.809490906136802,
            4.979563633831545,
            5.809490906136802,
            5.809490906136802,
            5.809490906136802
        ],[
            0.950497746248196,
            1.7689819166285872,
            3.304111213148491,
            1.742579201455026,
            2.5874660870089783,
            1.3729411890251721,
        ]
    ]
    rsource074 = [
        [
            -10,
            4.9795575666192935,
            -10,
            5.809483827722509,
            5.809483827722509,
            11.618967655445019
        ],[
            -10,
            0.47541265629394047,
            -10,
            1.5318852258360305,
            3.301476779819031,
            2.3242396529925977,
        ]
    ]
    rsource075 = [
        [
            8.71422145903946,
            11.618961945385948,
            -20,
            17.42844291807892,
            6.9713771672315685,
            2.4897775597255603,
        ],[
            0.43570238615249224,
            0.158437231328179,
            0.158437231328179,
            0.158437231328179,
            0.2218121238594506,
            0.2376558469922685,
        ]
    ]
    # 2um-flux density metric
    f2um052 = [
        [
            4.357118179602602,
            6.971389087364163,
            4.357118179602602,
            6.971389087364163,
            5.809490906136802,
            5.809490906136802,
        ],
        [
            0.63367,
            1.7109,
            2.7129,
            1.5525,
            2.2178,
            0.97690,
        ]
    ]
    f2um074 = [
        [
            -20,
            5.809483827722509,
            -10,
            3.4856902966335057,
            4.357112870791882,
            6.971380593267011,
        ],[
            0.1584708854313135,
            0.6867071702023585,
            -10,
            1.140990375105457,
            2.674196191653415,
            1.679791385571923,
        ]
    ]
    f2um075 = [
        [
            2.4897775597255603,
            4.9795551194511205,
            8.71422145903946,
            8.71422145903946,
            5.809480972692974,
            17.42844291807892,
        ],[
            0.5318964194588867,
            0.49794558417427687,
            0.475311693984537,
            0.2376558469922685,
            0.42249928354181066,
            0.158437231328179,
        ]
    ]
    # 10um flux density matric
    f10um052 = [
        [
            3.872993937424535,
            6.971389087364163,
            17.428472718410408,
            4.979563633831545,
            6.971389087364163,
            4.979563633831545,
        ],[
            0.8096832653225374,
            1.077230779081289,
            1.267330328330928,
            0.5431415692846835,
            0.3801990984992784,
            0.8599741513674155,
        ]
    ]
    f10um074 = [
        [
            6.971380593267011,
            5.809483827722509,
            4.9795575666192935,
            17.428451483167528,
            -10,
            6.971380593267011,
        ],[
            1.2043787292779826,
            1.0300607553035377,
            0.656522239644013,
            0.47541265629394047,
            -10,
            0.728966072984042,
        ]
    ]
    f10um075 = [
        [
            -10,
            17.42844291807892,
            -20,
            -10,
            -20,
            17.42844291807892,
        ],[
            -10,
            0.2376558469922685,
            0.475311693984537,
            -10,
            0.316874462656358,
            0.316874462656358,
        ]
    ]
    # Declare variables for average periods
    plot_average052 = 0
    plot_average074 = 0
    plot_average075 = 0
    N052_zeros = 0
    N074_zeros = 0
    N075_zeros = 0
    #
    # Plot Rsource metric
    for nangle in range(Nangles):

        # Plot 052
        # rsource
        ax[0].plot(nangle-0.05,rsource052[0][nangle],'go',markersize=6)
        ax[0].plot(
            [nangle-0.05,nangle-0.05],
            [
                rsource052[0][nangle]-0.5*rsource052[1][nangle],
                rsource052[0][nangle]+0.5*rsource052[1][nangle]
            ]
            ,'g'
        )
        # F2
        ax[0].plot(nangle,f2um052[0][nangle],'bs',markersize=6)
        ax[0].plot(
            [nangle,nangle],
            [
                f2um052[0][nangle]-0.5*f2um052[1][nangle],
                f2um052[0][nangle]+0.5*f2um052[1][nangle]
            ]
            ,'b'
        )
        # F10
        ax[0].plot(nangle+0.05,f10um052[0][nangle],'rd',markersize=6)
        ax[0].plot(
            [nangle+0.05,nangle+0.05],
            [
                f10um052[0][nangle]-0.5*f10um052[1][nangle],
                f10um052[0][nangle]+0.5*f10um052[1][nangle]
            ]
            ,'r'
        )
        # Plot 074
        # rsource
        ax[1].plot(nangle-0.05,rsource074[0][nangle],'go',markersize=6)
        ax[1].plot(
            [nangle-0.05,nangle-0.05],
            [
                rsource074[0][nangle]-0.5*rsource074[1][nangle],
                rsource074[0][nangle]+0.5*rsource074[1][nangle]
            ]
            ,'g'
        )
        # F2
        ax[1].plot(nangle,f2um074[0][nangle],'bs',markersize=6)
        ax[1].plot(
            [nangle,nangle],
            [
                f2um074[0][nangle]-0.5*f2um074[1][nangle],
                f2um074[0][nangle]+0.5*f2um074[1][nangle]
            ]
            ,'b'
        )
        #
        # F10
        ax[1].plot(nangle+0.05,f10um074[0][nangle],'rd',markersize=6)
        ax[1].plot(
            [nangle+0.05,nangle+0.05],
            [
                f10um074[0][nangle]-0.5*f10um074[1][nangle],
                f10um074[0][nangle]+0.5*f10um074[1][nangle]
            ]
            ,'r'
        )
        #
        # Plot 075
        # rsource
        ax[2].plot(nangle-0.05,rsource075[0][nangle],'go',markersize=6)
        ax[2].plot(
            [nangle-0.05,nangle-0.05],
            [
                rsource075[0][nangle]-0.5*rsource075[1][nangle],
                rsource075[0][nangle]+0.5*rsource075[1][nangle]
            ]
            ,'g'
        )
        #
        # F2
        ax[2].plot(nangle,f2um075[0][nangle],'bs',markersize=6)
        ax[2].plot(
            [nangle,nangle],
            [
                f2um075[0][nangle]-0.5*f2um075[1][nangle],
                f2um075[0][nangle]+0.5*f2um075[1][nangle]
            ]
            ,'b'
        )
        # F10
        ax[2].plot(nangle+0.05,f10um075[0][nangle],'rd',markersize=6)
        ax[2].plot(
            [nangle+0.05,nangle+0.05],
            [
                f10um075[0][nangle]-0.5*f10um075[1][nangle],
                f10um075[0][nangle]+0.5*f10um075[1][nangle]
            ]
            ,'r'
        )
        # Sum together all periods of all metrics for this model that are
        # not upper limits
        # Ie remove all negative numbers, ie alot of if statements
        # And count number of zeros for each model
        # NOTE this MUST be in the end of this loop
        if rsource052[0][nangle] < 0:
            rsource052[0][nangle] = 0
            N052_zeros += 1
        if f2um052[0][nangle] < 0:
            f2um052[0][nangle] = 0
            N052_zeros += 1
        if f10um052[0][nangle] < 0:
            f10um052[0][nangle] = 0
            N052_zeros += 1
        #
        if rsource074[0][nangle] < 0:
            rsource074[0][nangle] = 0
            N074_zeros += 1
        if f2um074[0][nangle] < 0:
            f2um074[0][nangle] = 0
            N074_zeros += 1
        if f10um074[0][nangle] < 0:
            f10um074[0][nangle] = 0
            N074_zeros += 1
        #
        if rsource075[0][nangle] < 0:
            rsource075[0][nangle] = 0
            N075_zeros += 1
        if f2um075[0][nangle] < 0:
            f2um075[0][nangle] = 0
            N075_zeros += 1
        if f10um075[0][nangle] < 0:
            f10um075[0][nangle] = 0
            N075_zeros += 1
        #
        # And add them all together
        plot_average052 += rsource052[0][nangle] + f2um052[0][nangle] + f10um052[0][nangle]
        plot_average074 += rsource074[0][nangle] + f2um074[0][nangle] + f10um074[0][nangle]
        plot_average075 += rsource075[0][nangle] + f2um075[0][nangle] + f10um075[0][nangle]
    #
    # And take average of them
    # 3*6 numbers per model (6 angles, 3 metrics) = 18 minus all zeros!
    plot_average052 /= 18 - N052_zeros
    plot_average074 /= 18 - N074_zeros
    plot_average075 /= 18 - N075_zeros
    #
    # And plot these averages
    ax[0].plot(
        [-1,6],[plot_average052,plot_average052],'k--'
    )
    ax[1].plot(
        [-1,6],[plot_average074,plot_average074],'k--'
    )
    ax[2].plot(
        [-1,6],[plot_average075,plot_average075],'k--'
    )
    print('Plot averages')
    print(f'  052: {plot_average052} yrs')
    print(f'  074: {plot_average074} yrs')
    print(f'  075: {plot_average075} yrs')

    # Plot those with not period with a line from zero to show the event length
    #
    ax[1].plot([0,0,],[0,f2um074[1][0]],'b',linewidth = 3)
    #
    ax[2].plot([2,2],[0,rsource075[1][2]],'g',linewidth = 3)
    ax[2].plot([2,2],[0,f10um075[1][2]],'r',linewidth = 3)
    ax[2].plot([4,4],[0,f10um075[1][4]],'r',linewidth = 3)
    #
    # Plot averaged period of all period for each model to compare with table
    # 052:
    # (5.8648 + 10.844 + 34.856945436820816 + 5.7127 + 12.297 + 7.5339)/6
    # 074:
    # (16.322 + 30.016 + 34.856945436820816 + 15.056 + 18.743 + 12.836)/6
    # 75:
    # (30.500 + 34.856945436820816 + 18.037 + 8.0226 + 34.85688583615784 + 29.047)/6
    average_period052 = 12.85155757280347
    average_period074 = 21.304990906136798
    average_period075 = 25.886738545496442
    # Plots    
    ax[0].plot(
        [-1,6],[average_period052,average_period052],'k:'
    )
    ax[1].plot(
        [-1,6],[average_period074,average_period074],'k:'
    )
    ax[2].plot(
        [-1,6],[average_period075,average_period075],'k:'
    )
    # List the labels so that theres 1 per model.
    labelpanel = 0
    ax[labelpanel].plot(-1,-1,'go',markersize=6,label=r'$R_{\rm source}$')
    ax[labelpanel].plot(-1,-1,'bs',markersize=6,label=r'$F(2\,\mu$m$)$')
    ax[labelpanel].plot(-1,-1,'rd',markersize=6,label=r'$F(10\,\mu$m$)$')
    ax[labelpanel].legend(
        #loc='upper left',
        fontsize=14
    )
    # Modify xticklabels and title
    for nmodel in range(Nmodels):
        ax[nmodel].set_title(models_label[nmodel],fontsize=14)
        ax[nmodel].set_xticks([0,1,2,3,4,5]) 
        ax[nmodel].set_xticklabels(angles_label) 
    # Set axislabels
    ax[0].set_ylabel(r'Average period \& event length (yrs)', fontsize=18)
    ax[1].set_xlabel(r'LOS-angle',fontsize=18)

    # Save figure
    fig.tight_layout()
    fig.savefig(
        'figs/periods_allmetrics.pdf', 
        facecolor='white',
        dpi=300
    )
    fig.show()
#
#####################################################################################
# Plot statistics on F2um och F10um for each model.
#
# %                          052           074            075
# %       2um:
# %           Fmean      :   68.73         80.89          95.14
# %           Fst        :   30.18         25.98          12.77
# %           Fmin       :    2.17          3.69          54.38
# %           Fmin/Fmean :    0.0316        0.0456         0.5716
# %           Fmin-proc  :    3.16%         4.56%         57.16%
# %
# %       10um
# %           Fmean      :   27.65         25.48          21.63
# %           Fst        :    5.01          4.23           1.94
# %           Fmax       :   39.28         35.79          27.23
# %           Fmax/Fmean :    1.4206        1.4049         1.2588
# %           Fmax-proc  :   42.06%        40.49%         25.88%
#
if plot_fluxvariations == 'y':
    #
    # Set path settings
    paths = [
        '../r3dresults/st28gm06n052_timedep_nospikes/',
        '../r3dresults/st28gm06n074_nospikes/',
        '../r3dresults/st28gm06n075_nospikes/'
    ]
    pathsstellar = [
        '../r3dresults/st28gm06n052_timedep_nodust/',
        '../r3dresults/st28gm06n074_nodust/',
        '../r3dresults/st28gm06n075_nodust/'
    ]
    modelabbreviations = [
        '052',
        '074',
        '075'
    ]
    Nangles = len(angles)
    # Set wavelengths
    obswavelengths = [
        2,10
    ]
    # Initiate subplot objects
    fig,ax = plt.subplots(
        len(obswavelengths),1,
        figsize=(6,6)
    )
    # Loop over models
    sed_average_2um = []
    sed_std_2um = []
    sed_min_2um = []

    sed_average_10um = []
    sed_std_10um = []
    sed_max_10um = []

    star_average_2um = []
    star_std_2um = []

    star_average_10um = []
    star_std_10um = []

    for path,pathstar in zip(paths,pathsstellar):

        # Extract list of folders from the path folder, and make them to numbers, and sort them!
        phases = [int(filename) for filename in os.listdir(path) if os.path.isdir(path+filename)]
        phases.sort()
        Nphases = len(phases)
        # And load corresponding time
        snapshot_times = np.loadtxt(path+'snapshot_yr.dat')[:,1]
        # Load average SED-data
        # 0           1         2    3     4
        # Wavelength  AverageF  STD  MaxF  MinF
        average_seds = np.loadtxt(f'{path}average_sed.dat')
        average_stars = np.loadtxt(f'{pathstar}average_sed.dat')
        wavelengths = average_seds[:,0]


        # Loop over obswavelengths
        for wavelength in obswavelengths:

            # Extract index of chosen wavelength and average flux density at wavelength
            waveindex = np.argwhere(wavelengths >= wavelength)[0][0] - 1

            if wavelength == 2:
                # Save data for each
                sed_average_2um.append(average_seds[waveindex,1]*1e-6)
                sed_std_2um.append(average_seds[waveindex,2]*1e-6)
                sed_min_2um.append(average_seds[waveindex,4]*1e-6)

                # Save stellar fluxes
                star_average_2um.append(average_stars[waveindex,1]*1e-6)
                star_std_2um.append(average_stars[waveindex,2]*1e-6)

            if wavelength == 10:
                # Save data for each
                sed_average_10um.append(average_seds[waveindex,1]*1e-6)
                sed_std_10um.append(average_seds[waveindex,2]*1e-6)
                sed_max_10um.append(average_seds[waveindex,3]*1e-6)

                # Save stellar fluxes
                star_average_10um.append(average_stars[waveindex,1]*1e-6)
                star_std_10um.append(average_stars[waveindex,2]*1e-6)


    # Add fields with stellar fluxes
    # 2um
    ax[0].fill_between(
        [0,1,2], 
        [
            star_average_2um[0]-star_std_2um[0],
            star_average_2um[1]-star_std_2um[1],
            star_average_2um[2]-star_std_2um[2]
        ],
        [
            star_average_2um[0]+star_std_2um[0],
            star_average_2um[1]+star_std_2um[1],
            star_average_2um[2]+star_std_2um[2]
        ],
        color='khaki',linestyle='-'
    )
    ax[0].plot(
        [0,1,2],
        star_average_2um,'darkgoldenrod'
    )
    # 10um
    ax[1].fill_between(
        [0,1,2], 
        [
            star_average_10um[0]-star_std_10um[0],
            star_average_10um[1]-star_std_10um[1],
            star_average_10um[2]-star_std_10um[2]
        ],
        [
            star_average_10um[0]+star_std_10um[0],
            star_average_10um[1]+star_std_10um[1],
            star_average_10um[2]+star_std_10um[2]
        ],
        color='khaki',linestyle='-'
    )
    ax[1].plot(
        [0,1,2],
        star_average_10um,'darkgoldenrod'
    )

    # And then plot dusty flux densities
    for nmodel in range(len(paths)):
        # Set some plot colours
        if nmodel == 0:
            linecolour = 'darkblue'
            allcolour = 'deepskyblue'
        if nmodel == 1:
            linecolour = 'darkred'
            allcolour = 'salmon' #lr tomato
        if nmodel == 2:
            linecolour = 'darkgreen'
            allcolour = 'lawngreen'

        # Plot 2um
        # First mininum
        ax[0].plot(
            [nmodel-0.05,nmodel+0.05],
            [sed_min_2um[nmodel],sed_min_2um[nmodel]],
            color=allcolour
        )
        ax[0].plot(
            [nmodel,nmodel],
            [sed_min_2um[nmodel],sed_average_2um[nmodel]],':',
            color=allcolour
        )
        # And then average-std
        ax[0].plot(
            nmodel,sed_average_2um[nmodel],'.',markersize=10,color=linecolour
        )
        ax[0].plot(
            [nmodel,nmodel],
            [sed_average_2um[nmodel]-sed_std_2um[nmodel],sed_average_2um[nmodel]+sed_std_2um[nmodel]],
            color=linecolour
        )

        # Plot 10um
        # First maximum
        ax[1].plot(
            [nmodel-0.05,nmodel+0.05],
            [sed_max_10um[nmodel],sed_max_10um[nmodel]],
            color=allcolour
        )
        ax[1].plot(
            [nmodel,nmodel],
            [sed_max_10um[nmodel],sed_average_10um[nmodel]],':',
            color=allcolour
        )
        # And then average-std
        ax[1].plot(
            nmodel,sed_average_10um[nmodel],'.',markersize=10,color=linecolour
        )
        ax[1].plot(
            [nmodel,nmodel],
            [sed_average_10um[nmodel]-sed_std_10um[nmodel],sed_average_10um[nmodel]+sed_std_10um[nmodel]],
            color=linecolour
        )
    # Final settings for the plots

    # Set xticklabels
    ax[0].set_xticks([0,1,2])
    ax[0].set_xticklabels(['','','']) 
    ax[1].set_xticks([0,1,2]) 
    ax[1].set_xticklabels(models_label) 
    ax[0].set_xlim(-0.3,2.3)
    ax[1].set_xlim(-0.3,2.3)

    # Set final general axis-settings
    ax[0].set_ylabel(r'$F(2\,\mu$m$)$, MJy at 1\,pc', fontsize=18)
    ax[1].set_ylabel(r'$F(10\,\mu$m$)$, MJy at 1\,pc', fontsize=18)
    ax[0].tick_params(axis='both', which='major', labelsize=15)
    ax[1].tick_params(axis='both', which='major', labelsize=15)
    ax[0].set_ylim(0,120)
    ax[1].set_ylim(15,45)

    fig.tight_layout()
    fig.savefig(
        'figs/dustfluxes.pdf',
        facecolor='white',
        dpi=300
    )
    fig.show()

# Plot comparisons with data as found at 
# suh2021, smiths2001, aavso.rg
# 
if plot_datacompare == 'y':
    # First panel plot with colour statistics
    # K[2.2](2mass) - W3[12] of O-rich AGB stars from Suh2021
    # both WISE and IRAS catalogue, a total of 4038(IRAS)+5253(WISE) AGB stars
    #
    # My model colours are computed with appropiate filters
    # as downloaded from SVO-filter databse.
    # http://svo2.cab.inta-csic.es/svo/theory/fps3/index.php
    #
    # All data:
    #                Suh+2021 
    #                fig9 IRAS      fig9 WISE
    #               K[2.2]-W3[12]   K[2.2]-W3[12]    n052       n074        n075
    #   maxK-minW3:  -               -               0.747      0.553       0.757
    #          Max: 0.9             1                1.727      1.580       1.317
    #          STD: 2.18            0.91             0.313      0.187       0.041
    #        Medel: 4.39            2.15             1.942      1.704       1.416
    #          STD: 2.18            0.91             0.313      0.187       0.041
    #          Min: 10-11           6                4.742      4.145       1.655
    #   minK-maxW3:  -               -               5.722      5.171       2.214
    #
    # X-axis to have models
    modelnames = [
        'st28gm06n052',
        'st28gm06n074',
        'st28gm06n075',
    ]
    modelsymbol = [
        'd','s','o'
    ]
    # Y-axel, div färger
    # modellerna med "error-bar-aktiga" linjer

    # Y-axis with model data in errorbar like
    # max-min, avr+std, avr, avr-std, , min-max
    # 052,074,075
    # Numbers from jupyternotebook with explore_timedep.ipynb
    modeldata = [
        [
            0.747 ,  1.381 , 1.942 , 2.685 , 5.722
        ],
        [
            0.553 , 1.251 , 1.704 , 2.240 , 5.171 
        ],
        [
            0.757 , 1.188 , 1.416 , 1.652 , 2.214
        ]
    ]
    linestyles = [
        ':','--','-'
    ]
    # Set up figure object
    fig,ax = plt.subplots(
        2,1,
        figsize=(6,9)
    )
    # For VR-prop:
    #fig,ax = plt.subplots(
    #    1,2,
    #    figsize=(10.5,4.5)
    #)


    # Plot fields for observed statistics from Suh21
    #
    # From Fig9, 2 right panels two catalogues, 
    #   first IRAScat-minmax
    #   then  WISEcat-minmax
    #   In middle, first: IRAS-mean, second: WISE-mean


    suhdata = [
        0.34,   # iras min
        0.58,   # wise min
        4.39,   # iras mean
        2.15,   # wise mean
        8.02,   # wise max
        13.6,   # iras max
    ]
    suhstd = [
        2.18,   # iras std
        0.91,   # wise std
    ]
    suhcolours = [
        'khaki','paleturquoise','red','blue'
    ]
    modelaxis = [0,4]
    modelpositions = [0.8,2,3.2]
    for nfield in range(4):
        if nfield < 2:
            ax[0].fill_between(
                modelaxis,
                [suhdata[nfield],suhdata[nfield]],
                [suhdata[-nfield-1],suhdata[-nfield-1]],
                color=suhcolours[nfield]
            )
        if nfield > 1:
            # Plot average
            ax[0].plot(
                modelaxis,
                [suhdata[nfield],suhdata[nfield]],
                color=suhcolours[nfield],linewidth=3
            )
            # and std
            ax[0].plot(
                modelaxis,
                [suhdata[nfield] + suhstd[nfield-2],suhdata[nfield] + suhstd[nfield-2]],
                '--',color=suhcolours[nfield],linewidth=2
            )
            ax[0].plot(
                modelaxis,
                [suhdata[nfield] - suhstd[nfield-2],suhdata[nfield] - suhstd[nfield-2]],
                '--',color=suhcolours[nfield],zorder=3,linewidth=2
            )

    # Plot model colour ranges
    for nmodel,modeldat in enumerate(modeldata):
        # error bars
        for nn in range(2):
            ax[0].plot(
                [modelpositions[nmodel],modelpositions[nmodel]],
                [modeldat[nn],modeldat[-nn-1]],
                'k',linestyle=linestyles[nn],zorder=4
            )
            ax[0].plot(
                [modelpositions[nmodel]-0.2,modelpositions[nmodel]+0.2],
                [modeldat[nn],modeldat[nn]],
                'k',zorder=4
            )
            ax[0].plot(
                [modelpositions[nmodel]-0.2,modelpositions[nmodel]+0.2],
                [modeldat[-nn-1],modeldat[-nn-1]],
                'k',zorder=4
            )
        ax[0].plot(
            [modelpositions[nmodel],modelpositions[nmodel]],
            [modeldat[2],modeldat[2]],
            modelsymbol[nmodel],color='k',markersize=12,zorder=4
        )
    # Set xlabels and tick settings
    ax[0].set_xlim(0,4)
    ax[0].set_xticks(modelpositions) 
    ax[0].set_xticklabels(modelnames) 
    ax[0].set_ylabel('K[2.2]$-$W3[12]', fontsize=18)
    ax[0].tick_params(axis='both', which='major', labelsize=15)
    #
    #
    # NEXT PANEL
    # compare with data one specific sources
    #
    # Y-axis: minimum/average flux ratio
    # X-axis: total time of included data
    #
    # Model data:
    modeltime =  34.857
    # each row: model
    # each column: Vband, Cobe2.2, Cobe3.5
    modelfluxratios = np.array([
        [0.00089898,0.04972 , 0.20425],
        [0.0033807 ,0.06202 , 0.23182],
        [0.058986  ,0.60118 , 0.69670],
    ])
    # Observed data:
    # Save all band data separately
    visualdata = np.array([
        [  8.06 ,0.083179  ],
        [ 14.47 ,0.15850   ],
        [131.18 ,0.012023  ],
        [143.96 ,0.044055  ],
        [144.73 ,0.039811  ],
        [171.37 ,0.0091201 ],
    ])
    Kbanddata = np.array([
        [0.68 , 0.700],
        [0.71 , 0.720],
        [0.71 , 0.667],
        [0.74 , 0.647],
        [0.74 , 0.558],
    ])
    W1banddata = np.array([
        [ 0.68, 0.857  ],
        [ 0.71, 0.674  ],
        [ 0.71, 0.765  ],
        [ 0.74, 0.75   ],
        [ 0.74, 0.664  ],
        [10.4 , 0.39810],
        [10.4 , 0.30199],
    ])
    allobsdata = [
        visualdata,
        Kbanddata,
        W1banddata
    ]
    # Start plotting
    #
    # Settings for subplotinset:
    # x0,y0,width,height
    axin = ax[1].inset_axes(
        bounds=[
            0.09,0.4,0.3,0.39
        ]
    )

    # Colours for each band
    wavecolour = [
        'b','orange','r'
    ]
    for nwave,obsdata in enumerate(allobsdata):

        obstime = obsdata[:,0]
        fluxratios = obsdata[:,1]

        ax[1].plot(obstime,fluxratios,'.',color=wavecolour[nwave])
        ax[1].plot(obstime,fluxratios,'--',color=wavecolour[nwave])
        axin.plot(obstime,fluxratios,'.',color=wavecolour[nwave])
        axin.plot(obstime,fluxratios,'--',color=wavecolour[nwave])

        # Plot model results
        for nmodel in range(3):
            ax[1].plot(modeltime,modelfluxratios[nmodel,nwave],modelsymbol[nmodel],color=wavecolour[nwave])
            #axin.plot(modeltime,modelfluxratios[nmodel,nwave],modelsymbol[nmodel],color=wavecolour[nwave])

    # Add model labels
    for nmodel,modelname in enumerate(modelnames):
        ax[1].plot(-1,-1,modelsymbol[nmodel],color='k',markersize=6,label=modelname)

    # Set axis settings
    ax[1].set_xscale('log')
    ax[1].set_yscale('log')
    axin.set_xlim(0.65,0.8)
    axin.set_ylim(0.5,0.9)
    #ax[1].indicate_inset_zoom(axin)
    ax[1].set_xlabel('Time series length (years)', fontsize=18)
    ax[1].tick_params(axis='both', which='major', labelsize=15)
    ax[1].set_ylabel(r'$F_\nu ($min$)/F_\nu ($average$)$', fontsize=18)
    ax[1].legend(
        #loc='upper right',
        fontsize=14
    )
    # Final figure fixes, and save fig.
    
    fig.tight_layout()
    fig.savefig(
        'figs/data_compare.pdf',
        facecolor='white',
        dpi=300
    )
    fig.show()

#######################################################################
# Plot for VR-prop

if plotvr_exampleimages == 'y':

    # Set paths
    path = '../r3dresults/st28gm06n052_timedep_nospikes/'
    modelabbreviation = '052'

    # Define wavelenvths
    wavelengths = [
        '02','10'
    ]
    # Chose snapshots
    snapshots = [
        289, 292, 295, 298
    ]
    snapshots_times = []

    # Extract corresponding snapshot-times
    snapshot_file = np.loadtxt(
        path+'snapshot_yr.dat'
    )
    for snaptime in snapshot_file:
        for snapshot in snapshots:
            if snapshot == snaptime[0]:
                snapshots_times.append(snaptime[1])

    # Create image objects and fill subplots
    fig,ax = plt.subplots(
        len(wavelengths),len(snapshots),
        figsize=(10, 5)
    )
    # Loop through wavelengths
    for nrow,wavelength in enumerate(wavelengths):
        imagefilename = f'image_i000_phi000_{wavelength}um.out'

        # Loop through snapshots
        for ntime,snapshot in enumerate(snapshots):

            # Load image
            image2d,image2dlog,flux,axisplot = a3d.load_images(
                path=f'{path}{snapshot}/',
                image=imagefilename,
                distance=1
            )
            # Compute and apply gamma function of each image
            scale = np.max(image2d)-np.min(image2d)
            gamma = 0.3*np.log(image2d.max())/np.log(image2d.mean())
            imageplot = ((image2d / scale)**gamma) * scale
            # Plot image
            ax[nrow,ntime].imshow(
                imageplot, origin='lower', extent=axisplot, cmap=plt.get_cmap('hot')
            )
            # write time and flux on top of each column
            if nrow == 0:
                ax[0,ntime].set_title(
                    f'{snapshots_times[ntime]:.2f} yrs, {flux*1e-6:.3f} MJy',
                    fontsize = 10
                )
            # and only flux on second row
            if nrow == 1:
                ax[1,ntime].set_title(
                    f'{flux*1e-6:.3f} MJy',
                    fontsize = 10
                )

    # Offset AU på yttre plots
    ax[0,0].set_ylabel('Offset (au)', fontsize = 14)
    ax[1,0].set_ylabel('Offset (au)', fontsize = 14)
    for ntime in range(len(snapshots)):
        ax[1,ntime].set_xlabel('Offset (au)', fontsize = 14)

    fig.tight_layout()
    fig.savefig(
        'figs/vrexampleimages.pdf', 
        facecolor='white',
        dpi=300
    )
    #fig.show()

if plotvr_radiusplot == 'y':
    # Plot radial variations of models compared to
    # This original stellar radius
    Rstar = 1.65

    # Chose model
    models = [
        #'st28gm06n052_timedep_nodust',
        #'st28gm06n052_timedep_nospikes',
        'st28gm06n074_nodust',
        'st28gm06n074_nospikes',
        #'st28gm06n075_nodust',
        #'st28gm06n075_nospikes',
    ]
    # Chose wavelength (in um)
    wavelengths = [
        #'01',
        '02',
        #'10'
    ]
    # Set LOS-angles
    angles = [
        'i000_phi000',
        'i090_phi000',
        'i090_phi090',
        'i090_phi270',
        'i180_phi000',
        'i270_phi000',
    ]
    # Set up figure
    fig, ax = plt.subplots(
        1,len(models),
        figsize=(11,4)
    )
    # Plot all
    for nmodel,modelname in enumerate(models):
        for wavelength in wavelengths:
            
            path = f'../r3dresults/{modelname}/'

            # Load and plot radii vs time
            source_radii = np.loadtxt(f'{path}/source_radius_{wavelength}um.dat')

            # Load corresponding time
            snapshot_times = np.loadtxt(path+'snapshot_yr.dat')[:,1]

            # Set figure object
            ax[nmodel].set_xlabel(r'Sim. time (yrs)',fontsize=18)
            ax[nmodel].tick_params(axis='both', which='major', labelsize=15)
            ax[nmodel].set_xlim([
                snapshot_times[0],snapshot_times[-1]
            ])
            ax[nmodel].set_ylim([
                0,1.25
            ])

            # First plot angles
            for nangle in range(len(angles)):
                ax[nmodel].plot(snapshot_times,source_radii[:,nangle+2]/Rstar,'lightgrey')
            # Then angle-averaged
            ax[nmodel].plot(snapshot_times,source_radii[:,1]/Rstar,'k')

            # Plot "table radius" and average radius
            source_radii_average = np.mean(source_radii[:,1])
            ax[nmodel].plot(
                [snapshot_times[0],snapshot_times[-1]],
                [1,1],'r:'
            )
            ax[nmodel].plot(
                [snapshot_times[0],snapshot_times[-1]],
                [source_radii_average/Rstar,source_radii_average/Rstar],'--',color='mediumblue'
            )
    # Final settings
    ax[0].set_ylabel(r'Radius ($R_\star$)', fontsize=18)
    # Save figure
    fig.tight_layout()
    fig.savefig(f"figs/{models[0].split('_')[0]}_sourceradius_{wavelength}um.pdf", dpi=300, facecolor="white")

# Plot special data compare for VR prop

# Plot comparisons with data as found at 
# suh2021, smiths2001, aavso.rg
# 
if plotvr_datacompare == 'y':
    # First panel plot with colour statistics
    # K[2.2](2mass) - W3[12] of O-rich AGB stars from Suh2021
    # both WISE and IRAS catalogue, a total of 4038(IRAS)+5253(WISE) AGB stars
    #
    # My model colours are computed with appropiate filters
    # as downloaded from SVO-filter databse.
    # http://svo2.cab.inta-csic.es/svo/theory/fps3/index.php
    #
    # All data:
    #                Suh+2021 
    #                fig9 IRAS      fig9 WISE
    #               K[2.2]-W3[12]   K[2.2]-W3[12]    n052       n074        n075
    #   maxK-minW3:  -               -               0.747      0.553       0.757
    #          Max: 0.9             1                1.727      1.580       1.317
    #          STD: 2.18            0.91             0.313      0.187       0.041
    #        Medel: 4.39            2.15             1.942      1.704       1.416
    #          STD: 2.18            0.91             0.313      0.187       0.041
    #          Min: 10-11           6                4.742      4.145       1.655
    #   minK-maxW3:  -               -               5.722      5.171       2.214
    #
    # X-axis to have models
    modelnames = [
        'st28gm06n052',
        'st28gm06n074',
        'st28gm06n075',
    ]
    modelsymbol = [
        'd','s','o'
    ]
    # Y-axel, div färger
    # modellerna med "error-bar-aktiga" linjer

    # Y-axis with model data in errorbar like
    # max-min, avr+std, avr, avr-std, , min-max
    # 052,074,075
    # Numbers from jupyternotebook with explore_timedep.ipynb
    modeldata = [
        [
            0.747 ,  1.381 , 1.942 , 2.685 , 5.722
        ],
        [
            0.553 , 1.251 , 1.704 , 2.240 , 5.171 
        ],
        [
            0.757 , 1.188 , 1.416 , 1.652 , 2.214
        ]
    ]
    linestyles = [
        ':','--','-'
    ]
    # Set up figure object
    fig,ax = plt.figure(figsize=(7,5)), plt.axes()
    #
    # Plot fields for observed statistics from Suh21
    #
    # From Fig9, 2 right panels two catalogues, 
    #   first IRAScat-minmax
    #   then  WISEcat-minmax
    #   In middle, first: IRAS-mean, second: WISE-mean
    #
    suhdata = [
        0.34,   # iras min
        0.58,   # wise min
        4.39,   # iras mean
        2.15,   # wise mean
        8.02,   # wise max
        13.6,   # iras max
    ]
    suhstd = [
        2.18,   # iras std
        0.91,   # wise std
    ]
    suhcolours = [
        'khaki','paleturquoise','red','blue'
    ]
    modelaxis = [0,4]
    modelpositions = [0.8,2,3.2]
    for nfield in range(4):
        if nfield < 2:
            ax.fill_between(
                modelaxis,
                [suhdata[nfield],suhdata[nfield]],
                [suhdata[-nfield-1],suhdata[-nfield-1]],
                color=suhcolours[nfield]
            )
        if nfield > 1:
            # Plot average
            ax.plot(
                modelaxis,
                [suhdata[nfield],suhdata[nfield]],
                color=suhcolours[nfield],linewidth=3
            )
            # and std
            ax.plot(
                modelaxis,
                [suhdata[nfield] + suhstd[nfield-2],suhdata[nfield] + suhstd[nfield-2]],
                '--',color=suhcolours[nfield],linewidth=2
            )
            ax.plot(
                modelaxis,
                [suhdata[nfield] - suhstd[nfield-2],suhdata[nfield] - suhstd[nfield-2]],
                '--',color=suhcolours[nfield],zorder=3,linewidth=2
            )
    #
    # Plot model colour ranges
    for nmodel,modeldat in enumerate(modeldata):
        # error bars
        for nn in range(2):
            ax.plot(
                [modelpositions[nmodel],modelpositions[nmodel]],
                [modeldat[nn],modeldat[-nn-1]],
                'k',linestyle=linestyles[nn],zorder=4,
            )
            ax.plot(
                [modelpositions[nmodel]-0.2,modelpositions[nmodel]+0.2],
                [modeldat[nn],modeldat[nn]],
                'k',zorder=4
            )
            ax.plot(
                [modelpositions[nmodel]-0.2,modelpositions[nmodel]+0.2],
                [modeldat[-nn-1],modeldat[-nn-1]],
                'k',zorder=4
            )
        ax.plot(
            [modelpositions[nmodel],modelpositions[nmodel]],
            [modeldat[2],modeldat[2]],
            modelsymbol[nmodel],color='k',markersize=12,zorder=4
        )
    # Set xlabels and tick settings
    ax.set_xlim(0,4)
    ax.set_xticks(modelpositions) 
    ax.set_xticklabels(modelnames) 
    ax.set_ylabel('K[2.2]$-$W3[12]', fontsize=18)
    ax.tick_params(axis='both', which='major', labelsize=15)
    #
    # Final figure fixes, and save fig.
    fig.tight_layout()
    fig.savefig(
        'figs/vrdata_compare.pdf',
        facecolor='white',
        dpi=300
    )
    fig.show()


