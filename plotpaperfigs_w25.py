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
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
})
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

# Set various general numbers
Rstar = 1.65                            # Stellar radius in AU
Rin = 2*Rstar                           # Annalus radii
Rout = 6*Rstar                          # Annalus radii
annulus_area = np.pi*(Rout-Rin)**2      # Annalus area in AU2
pix_area = (876562500000.0/AUcm)**2     # Pixel area in AU2 (521**2 & 30**2AU**2 images)



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
Nangles = len(angles)
# Set a list of model and model labels
models = [
    'st28gm06n052_timedep',
    'st28gm06n074',
    'st28gm06n075',
]
models_label = [
    r'Model-A: $\alpha_{\rm stick} = 1$',
    r'Model-B: $\alpha_{\rm stick} = 0.1$',
    r'Model-C: $\alpha_{\rm stick} = 0.01$',
]
Nmodels = len(models)

# Plot-list
plot_dustmass = 'n'             # Plots dust masses vs time of all 3 models
plot_075grainsize = 'n'         # Plots grain size vs time of 075
plot_052exampleimages = 'n'     # Plots a number of example images of 052
plot_numbclouds = 'n'           # Plots number of clouds per time for each angle&model
plot_datacompare = 'n'          # Plots colour comparisons for each model with data
plot_LOSevents = 'n'            # Plots angle-dependent cloud-periods and probabilities
plot_cloudareas = 'n'           # Plots histogram of N clouds per area size
plot_bestrandomsample = 'n'     # Plots three example figures and cloud sizes
plot_allrandomsample = 'y'      # Plots all 24 random images of all models


# For vr-prop
plotvr_exampleimages = 'n'
plotvr_radiusplot = 'n'
plotvr_datacompare = 'n'

# For whyAGB5-talk
# TODO



# SKIP THESE
plot_allseds = 'n'              # SKIP
plot_luminosities = 'n'         # SKIP
plot_052fluxdensity = 'n'       # SKIP
plot_fluxvariations = 'n'       # SKIP



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

    # Print average and median dust masses
    print(f'052\n  mean: {np.mean(dustmass052)} Msol\n  median: {np.median(dustmass052)} Msol\n')
    print(f'074\n  mean: {np.mean(dustmass074)} Msol\n  median: {np.median(dustmass074)} Msol\n')
    print(f'075\n  mean: {np.mean(dustmass075)} Msol\n  median: {np.median(dustmass075)} Msol\n')

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
    ax[0].set_xlim([0,np.max([time052[-1],time074[-1],time075[-1]])])
    ax[0].set_ylim([0,3.5e-7])
    ax[1].set_xlim([0,time075[-1]])
    ax[1].set_ylim([0,1.6e-8])
    ax[0].legend(fontsize=14)
    ax[1].legend(fontsize=14)

    # Save figure
    fig.tight_layout()
    ax[0].yaxis.set_label_coords(-0.09,0.09) # Moves ylabel to centre of vertical
    fig.savefig(
        f'figs/052_074_075_dustmasscompare.pdf', 
        dpi=300, 
        facecolor="white"
    )
    #fig.show()
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
    ax.set_xlim([phasetimes075[0],phasetimes075[-1]])
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
    #fig.show()
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
    #fig.show()
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
        path = f'{path}{phases[0]:03d}/spectrum_{angles[0]}.out'
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
                path = f'{path}{phase:03d}/spectrum_{angle}.out'
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
    plt.savefig(
        f'figs/{model}_fluxtime_{wavelength}um.pdf', 
        facecolor='white',
        dpi=300
    )
    #fig.show()
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
            # Change to MJy/asec2 or au2 at 1pc
            image2d = image2d/1e6
            # Set scales for the two wavelengths
            if ncolumn == 0:
                scale = [1e-2,5]
            if ncolumn == 1:
                scale = [1e-2,9e-1]
            # Plot image and save colourbar info
            imbar = ax[ntime,ncolumn].imshow(
                image2d, 
                origin='lower', 
                extent=axisplot, 
                vmin=scale[0],
                vmax=scale[1],
                cmap=plt.get_cmap('hot')
            )
            # write time on top of each left column
            # flux: , {flux*1e-6:.3f} MJy
            if ncolumn == 0:
                ax[ntime,0].set_title(
                    f'{snapshots_times[ntime]:.2f} yrs',
                    fontsize = 10,
                    loc='left'
                )
            if ncolumn == 1:
                # and only flux on second column
                #ax[ntime,1].set_title(
                #    f'{flux*1e-6:.3f} MJy',
                #    fontsize = 10
                #)
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
                #divider = make_axes_locatable(ax[ntime,ncolumn])
                #cax = divider.append_axes(
                #    'bottom', 
                #    size='6%', 
                #    pad=0.7,
                #)
                cb0 = plt.colorbar(
                    imbar, 
                    #cax=cax, 
                    orientation = 'horizontal', 
                )
                cb0.ax.tick_params(labelsize=12)
                # Only label on one, TODO move this to middle of plot (manual labour in gimp)
                if ncolumn == 0:
                    cb0.set_label(
                        label = r'Flux density (MJy au$^{-2}$, at 1pc)', fontsize=12
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
    """
    All data for this plot

        Total included time: 379.1480984769462 time units
    MODEL  Area limit    N clouds   N clouds/LOS  Total port.   Aver period (yrs)   FFT-period (yrs)
    052     0.1          168        28.00         0.4023        2.26                4.21
            0.2          142        23.67         0.2580        2.67                4.21
            0.3          110        18.33         0.1697        3.45                4.21
            0.4           81        13.50         0.1184        4.68               12.64
            0.5           64        10.67         0.0775        5.92               12.23
            0.6           45         7.50         0.0509        8.43               12.64

        Total included time: 426.6804631133246 time units
    074     0.1          133        22.17         0.2812        3.21                5.93
            0.2           95        15.83         0.1489        4.49                6.46
            0.3           63        10.50         0.0848        6.77                7.90
            0.4           48         8.00         0.0515        8.89               11.85
            0.5           32         5.33         0.0326       13.33               11.85
            0.6           17         2.83         0.0182       25.10               11.85

        Total included time: 419.07515235389167 time units
    075     0.1           80        13.33         0.1475        5.24                3.68
            0.2           42         7.00         0.0551        9.98                9.98
            0.3           20         3.33         0.0230       20.95               23.28
            0.4            8         1.33         0.0094       52.38               34.92
            0.5            4         0.66         0.0034      104.77                N/A
            0.6            1         0.17         0.0008      419.08                N/A


    FRACT_AREA: 0.1 periods
    Model          0-0      90-0     90-90     90-270     180-0     270-0
    Model-A:052    1.98     2.18     2.43      2.53       2.11      2.43
    Model-B:074    3.39     3.23     3.23      3.09       3.09      3.23
    Model-C:075    5.82     5.82     4.37      4.37       5.37      6.35
                    Portion of total time
    Model-A:052    0.5187   0.4887   0.4009    0.3358     0.3433    0.3308
    Model-B:074    0.3252   0.2851   0.3029    0.2717     0.2673    0.2383
    Model-C:075    0.1609   0.1382   0.1405    0.1836     0.1337    0.1292


    FRACT_AREA: 0.2 periods
    Model          0-0      90-0     90-90     90-270     180-0     270-0
    Model-A:052    2.26     2.26     2.75      3.01       2.87      3.16
    Model-B:074    3.56     4.18     4.74      5.47       4.74      4.74
    Model-C:075    7.76    11.65     7.76     11.65       9.98     13.97
                    Portion of total time
    Model-A:052    0.3533   0.3007   0.2706    0.2055     0.2080    0.2130
    Model-B:074    0.1715   0.1782   0.1715    0.1403     0.1069    0.1270
    Model-C:075    0.0725   0.0567   0.0385    0.0431     0.0567    0.0635


    FRACT_AREA: 0.3 periods
    Model          0-0      90-0     90-90     90-270     180-0     270-0
    Model-A:052    2.87     3.16     3.01      4.21       3.95      3.95
    Model-B:074    5.08     5.47     7.90      7.11       8.89      7.90
    Model-C:075   17.47    23.29    34.94     17.47      34.94     13.97
                    Length
    Model-A:052    0.2456   0.2030   0.1779    0.1328     0.1253    0.1353
    Model-B:074    0.0846   0.1047   0.1136    0.0802     0.0512    0.0757
    Model-C:075    0.0227   0.0340   0.0091    0.0136     0.0249    0.0340


    FRACT_AREA: 0.4 periods
    Model          0-0      90-0     90-90     90-270     180-0     270-0
    Model-A:052    4.21     3.72     4.21      5.27       7.02      4.86
    Model-B:074    7.90     7.11     7.11      7.90      17.78     11.86
    Model-C:075   34.94    34.94     N/A      69.87      69.87     34.94
                    Portion of total time
    Model-A:052    0.1729   0.1504   0.1078    0.1002     0.0802    0.1002
    Model-B:074    0.0490   0.0690   0.0802    0.0445     0.0223    0.0445
    Model-C:075    0.0113   0.0204   N/A       0.0023     0.0136    0.0091


    FRACT_AREA: 0.5 periods
    Model          0-0      90-0     90-90     90-270     180-0     270-0
    Model-A:052    4.86     4.86     6.32      5.75       9.03      6.32
    Model-B:074   14.23    10.16    10.16     11.86      23.71     17.78
    Model-C:075    N/A     69.87     N/A       N/A       34.94     69.87
                    Length
    Model-A:052    0.1328   0.0977   0.0526    0.0526     0.0626    0.0677
    Model-B:074    0.0356   0.0423   0.0490    0.0312     0.0089    0.0290
    Model-C:075    N/A      0.0091   N/A       N/A        0.0091    0.0023


    FRACT_AREA: 0.6 periods
    Model          0-0      90-0     90-90     90-270     180-0     270-0
    Model-A:052    5.75     7.02    12.64     10.54      10.54      7.90
    Model-B:074   23.71    17.78    14.28     17.78       N/A      71.14
    Model-C:075    N/A     69.87     N/A       N/A        N/A       N/A
                    Portion of total time
    Model-A:052    0.0952   0.0702   0.0301    0.0301     0.0426    0.0376
    Model-B:074    0.0200   0.0267   0.0334    0.0134     N/A       0.0156
    Model-C:075    N/A      0.0045   N/A       N/A        N/A       N/A

    LOWER LIMITS (PERIODS)
    052:
        NONE
    074:
        0.6: 180-0, 270-0
    075:
        0.4: 90-90
        0.5: 0-0, 90-90, 90-270
        0.6: 0-0, 90-90, 90-270, 180-0, 270-0
    """;

    fig,ax = plt.subplots(
        2,3,
        figsize=(13,7)
    )
    for nrow in range(2):
        for nmodel in range(Nmodels):
            ax[nrow][nmodel].tick_params(axis='both', which='major', labelsize=15)
            ax[nrow][nmodel].set_xlim([-0.5,5.5])
            ax[1][nmodel].set_ylim([0,1])
    ax[0][0].set_ylim([0,14])
    ax[0][1].set_ylim([0,30])
    ax[0][2].set_ylim([0,77])
    farea_colours = ['g','c','b','m','r','darkorange']
    farea_markers = ['o','p','s','h','d','8']
    farea_label = ['0.1','0.2','0.3','0.4','0.5','0.6']
    # Combine each style to one set for each model, plot in one figure
    #
    # model[0][:] - perioder per LOS
    # model[1][:] - Andel av tiden med moln per LOS, sannolikhet för detektion
    # model_average[0] - medelvärde från stortabell på medelperiod
    # model_average[1] - medelvärde från stortabell på medelandel av tot tid
    #
    # Period per LOS in same LOS-order as listed in angle-labels
    # Fract_area = 0.1
    #
    farea01_052 = [
        [
            1.98, 2.18, 2.43, 2.53, 2.11, 2.43
        ],[
            0.5187, 0.4887, 0.4009, 0.3358, 0.3433, 0.3308
        ]
    ]
    farea01_052_average = [2.26,0.4023]
    #
    farea01_074 = [
        [
            3.39, 3.23, 3.23, 3.09, 3.09, 3.23
        ],[
            0.3252, 0.2851, 0.3029, 0.2717, 0.2673, 0.2383
        ]
    ]
    farea01_074_average = [3.21, 0.2812]
    #
    farea01_075 = [
        [
            5.82, 5.82, 4.37, 4.37, 5.37, 6.35
        ],[
            0.1609, 0.1382, 0.1405, 0.1836, 0.1337, 0.1292
        ]
    ]
    farea01_075_average = [5.24, 0.1475]
    #
    # Fract_area = 0.2
    #
    farea02_052 = [
        [
            2.26, 2.26, 2.75, 3.01, 2.87, 3.16
        ],[
            0.3533, 0.3007, 0.2706, 0.2055, 0.2080, 0.2130
        ]
    ]
    farea02_052_average = [2.67, 0.2580]
    #
    farea02_074 = [
        [
            3.56, 4.18, 4.74, 5.47, 4.74, 4.74
        ],[
            0.1715, 0.1782, 0.1715, 0.1403, 0.1069, 0.1270
        ]
    ]
    farea02_074_average = [4.49, 0.1489]
    #
    farea02_075 = [
        [
            7.76, 11.65, 7.76, 11.65, 9.98, 13.97
        ],[
            0.0725, 0.0567, 0.0385, 0.0431, 0.0567, 0.0635
        ]
    ]
    farea02_075_average = [9.98, 0.0551]
    #
    # Fract_area = 0.3
    #
    farea03_052 = [
        [
            2.87, 3.16, 3.01, 4.21, 3.95, 3.95
        ],[
            0.2456, 0.2030, 0.1779, 0.1328, 0.1253, 0.1353
        ]
    ]
    farea03_052_average = [3.45, 0.1697]
    #
    farea03_074 = [
        [
            5.08, 5.47, 7.90, 7.11, 8.89, 7.90
        ],[
            0.0846, 0.1047, 0.1136, 0.0802, 0.0512, 0.0757
        ]
    ]
    farea03_074_average = [6.77, 0.0848]
    #
    farea03_075 = [
        [
            17.47, 23.29, 34.94, 17.47, 34.94, 13.97
        ],[
            0.0227, 0.0340, 0.0091, 0.0136, 0.0249, 0.0340
        ]
    ]
    farea03_075_average = [20.95, 0.0230]
    #
    # Fract_area = 0.4
    #
    farea04_052 = [
        [
            4.21, 3.72, 4.21, 5.27, 7.02, 4.86
        ],[
            0.1729, 0.1504, 0.1078, 0.1002, 0.0802, 0.1002
        ]
    ]
    farea04_052_average = [4.68, 0.1184]
    #
    farea04_074 = [
        [
            7.90, 7.11, 7.11, 7.90, 17.78, 11.86
        ],[
            0.0490, 0.0690, 0.0802, 0.0445, 0.0223, 0.0445
        ]
    ]
    farea04_074_average = [8.89, 0.0515]
    #
    farea04_075 = [
        [
            34.94, 34.94, -10, 69.87, 69.87, 34.94
        ],[
            0.0113, 0.0204, -10, 0.0023, 0.0136, 0.0091
        ]
    ]
    farea04_075_average = [52.38, 0.0094]
    #
    # Fract_area = 0.5
    #
    farea05_052 = [
        [
            4.86, 4.86, 6.32, 5.75, 9.03, 6.32
        ],[
            0.1328, 0.0977, 0.0526, 0.0526, 0.0626, 0.0677
        ]
    ]
    farea05_052_average = [5.92, 0.0775]
    #
    farea05_074 = [
        [
            14.23, 10.16, 10.16, 11.86, 23.71, 17.78
        ],[
            0.0356, 0.0423, 0.0490, 0.0312, 0.0089, 0.0290
        ]
    ]
    farea05_074_average = [13.33, 0.0326]
    #
    farea05_075 = [
        [
            -10, 69.87, -10, -10, 34.94, 69.87
        ],[
            -10, 0.0091, -10, -10, 0.0091, 0.0023
        ]
    ]
    farea05_075_average =  [104.77, 0.0034]
    #
    # Fract_area = 0.6
    #
    farea06_052 = [
        [
            5.75, 7.02, 12.64, 10.54, 10.54, 7.90
        ],[
            0.0952, 0.0702, 0.0301, 0.0301, 0.0426, 0.0376
        ]
    ]
    farea06_052_average = [8.43, 0.0509]
    #
    farea06_074 = [
        [
            23.71, 17.78, 14.28, 17.78, -10, 71.14
        ],[
            0.0200, 0.0267, 0.0334, 0.0134, -10, 0.0156
        ]
    ]
    farea06_074_average = [25.10, 0.0182]
    #
    farea06_075 = [
        [
            -10, 69.87, -10, -10, -10, -10
        ],[
            -10, 0.0045, -10, -10, -10, -10
        ]
    ]
    farea06_075_average = [419.08, 0.0008]
    #
    # Plot farea
    #
    for nangle in range(Nangles):

        # Plot 052
        # Plot farea0.1
        ax[0][0].plot(
            nangle,
            farea01_052[0][nangle],
            color=farea_colours[0],
            marker=farea_markers[0],
            markersize=6
        )
        # Plot farea0.2
        ax[0][0].plot(
            nangle,
            farea02_052[0][nangle],
            color=farea_colours[1],
            marker=farea_markers[1],
            markersize=6
        )
        # farea0.3
        ax[0][0].plot(
            nangle,
            farea03_052[0][nangle],
            color=farea_colours[2],
            marker=farea_markers[2],
            markersize=6
        )
        # farea0.4
        ax[0][0].plot(
            nangle,
            farea04_052[0][nangle],
            color=farea_colours[3],
            marker=farea_markers[3],
            markersize=6
        )
        # farea0.5
        ax[0][0].plot(
            nangle,
            farea05_052[0][nangle],
            color=farea_colours[4],
            marker=farea_markers[4],
            markersize=6
        )
        # farea0.6
        ax[0][0].plot(
            nangle,
            farea06_052[0][nangle],
            color=farea_colours[5],
            marker=farea_markers[5],
            markersize=6
        )
        # Plot 074
        # Plot farea0.1
        ax[0][1].plot(
            nangle,
            farea01_074[0][nangle],
            color=farea_colours[0],
            marker=farea_markers[0],
            markersize=6
        )
        # Plot farea0.2
        ax[0][1].plot(
            nangle,
            farea02_074[0][nangle],
            color=farea_colours[1],
            marker=farea_markers[1],
            markersize=6
        )
        # farea0.3
        ax[0][1].plot(
            nangle,
            farea03_074[0][nangle],
            color=farea_colours[2],
            marker=farea_markers[2],
            markersize=6
        )
        # farea0.4
        ax[0][1].plot(
            nangle,
            farea04_074[0][nangle],
            color=farea_colours[3],
            marker=farea_markers[3],
            markersize=6
        )
        # farea0.5
        ax[0][1].plot(
            nangle,
            farea05_074[0][nangle],
            color=farea_colours[4],
            marker=farea_markers[4],
            markersize=6
        )
        # farea0.6
        ax[0][1].plot(
            nangle,
            farea06_074[0][nangle],
            color=farea_colours[5],
            marker=farea_markers[5],
            markersize=6
        )
        # Plot 075
        # Plot farea0.1
        ax[0][2].plot(
            nangle,
            farea01_075[0][nangle],
            color=farea_colours[0],
            marker=farea_markers[0],
            markersize=6
        )
        # Plot farea0.2
        ax[0][2].plot(
            nangle,
            farea02_075[0][nangle],
            color=farea_colours[1],
            marker=farea_markers[1],
            markersize=6
        )
        # farea0.3
        ax[0][2].plot(
            nangle,
            farea03_075[0][nangle],
            color=farea_colours[2],
            marker=farea_markers[2],
            markersize=6
        )
        # farea0.4
        ax[0][2].plot(
            nangle,
            farea04_075[0][nangle],
            color=farea_colours[3],
            marker=farea_markers[3],
            markersize=6
        )
        # farea0.5
        ax[0][2].plot(
            nangle,
            farea05_075[0][nangle],
            color=farea_colours[4],
            marker=farea_markers[4],
            markersize=6
        )
        # farea0.6
        ax[0][2].plot(
            nangle,
            farea06_075[0][nangle],
            color=farea_colours[5],
            marker=farea_markers[5],
            markersize=6
        )
    #
    # And plot averages
    # For 052
    #
    ax[0][0].plot(
        [-1,6],[farea01_052_average[0],farea01_052_average[0]],'--',color=farea_colours[0]
    )
    ax[0][0].plot(
        [-1,6],[farea02_052_average[0],farea02_052_average[0]],'--',color=farea_colours[1]
    )
    ax[0][0].plot(
        [-1,6],[farea03_052_average[0],farea03_052_average[0]],'--',color=farea_colours[2]
    )
    ax[0][0].plot(
        [-1,6],[farea04_052_average[0],farea04_052_average[0]],'--',color=farea_colours[3]
    )
    ax[0][0].plot(
        [-1,6],[farea05_052_average[0],farea05_052_average[0]],'--',color=farea_colours[4]
    )
    ax[0][0].plot(
        [-1,6],[farea06_052_average[0],farea06_052_average[0]],'--',color=farea_colours[5]
    )
    #
    # For 074
    #
    ax[0][1].plot(
        [-1,6],[farea01_074_average[0],farea01_074_average[0]],'--',color=farea_colours[0]
    )
    ax[0][1].plot(
        [-1,6],[farea02_074_average[0],farea02_074_average[0]],'--',color=farea_colours[1]
    )
    ax[0][1].plot(
        [-1,6],[farea03_074_average[0],farea03_074_average[0]],'--',color=farea_colours[2]
    )
    ax[0][1].plot(
        [-1,6],[farea04_074_average[0],farea04_074_average[0]],'--',color=farea_colours[3]
    )
    ax[0][1].plot(
        [-1,6],[farea05_074_average[0],farea05_074_average[0]],'--',color=farea_colours[4]
    )
    ax[0][1].plot(
        [-1,6],[farea06_074_average[0],farea06_074_average[0]],'--',color=farea_colours[5]
    )
    #
    # For 075
    #
    ax[0][2].plot(
        [-1,6],[farea01_075_average[0],farea01_075_average[0]],'--',color=farea_colours[0]
    )
    ax[0][2].plot(
        [-1,6],[farea02_075_average[0],farea02_075_average[0]],'--',color=farea_colours[1]
    )
    ax[0][2].plot(
        [-1,6],[farea03_075_average[0],farea03_075_average[0]],'--',color=farea_colours[2]
    )
    ax[0][2].plot(
        [-1,6],[farea04_075_average[0],farea04_075_average[0]],'--',color=farea_colours[3]
    )
    ax[0][2].plot(
        [-1,6],[farea05_075_average[0],farea05_075_average[0]],'--',color=farea_colours[4]
    )
    ax[0][2].plot(
        [-1,6],[farea06_075_average[0],farea06_075_average[0]],'--',color=farea_colours[5]
    )
    #
    # PLOT SECOND ROW: cloud-portions
    #
    for nangle in range(Nangles):
        #
        # Cloud portions 052
        # Plot farea0.1
        ax[1][0].plot(
            nangle,
            farea01_052[1][nangle],
            color=farea_colours[0],
            marker=farea_markers[0],
            markersize=6
        )
        # Plot farea0.2
        ax[1][0].plot(
            nangle,
            farea02_052[1][nangle],
            color=farea_colours[1],
            marker=farea_markers[1],
            markersize=6
        )
        # farea0.3
        ax[1][0].plot(
            nangle,
            farea03_052[1][nangle],
            color=farea_colours[2],
            marker=farea_markers[2],
            markersize=6
        )
        # farea0.4
        ax[1][0].plot(
            nangle,
            farea04_052[1][nangle],
            color=farea_colours[3],
            marker=farea_markers[3],
            markersize=6
        )
        # farea0.5
        ax[1][0].plot(
            nangle,
            farea05_052[1][nangle],
            color=farea_colours[4],
            marker=farea_markers[4],
            markersize=6
        )
        # farea0.6
        ax[1][0].plot(
            nangle,
            farea06_052[1][nangle],
            color=farea_colours[5],
            marker=farea_markers[5],
            markersize=6
        )
        #
        # Cloud portions 074
        # Plot farea0.1
        ax[1][1].plot(
            nangle,
            farea01_074[1][nangle],
            color=farea_colours[0],
            marker=farea_markers[0],
            markersize=6
        )
        # Plot farea0.2
        ax[1][1].plot(
            nangle,
            farea02_074[1][nangle],
            color=farea_colours[1],
            marker=farea_markers[1],
            markersize=6
        )
        # farea0.3
        ax[1][1].plot(
            nangle,
            farea03_074[1][nangle],
            color=farea_colours[2],
            marker=farea_markers[2],
            markersize=6
        )
        # farea0.4
        ax[1][1].plot(
            nangle,
            farea04_074[1][nangle],
            color=farea_colours[3],
            marker=farea_markers[3],
            markersize=6
        )
        # farea0.5
        ax[1][1].plot(
            nangle,
            farea05_074[1][nangle],
            color=farea_colours[4],
            marker=farea_markers[4],
            markersize=6
        )
        # farea0.6
        ax[1][1].plot(
            nangle,
            farea06_074[1][nangle],
            color=farea_colours[5],
            marker=farea_markers[5],
            markersize=6
        )
        #
        # Cloud portions 075
        # Plot farea0.1
        ax[1][2].plot(
            nangle,
            farea01_075[1][nangle],
            color=farea_colours[0],
            marker=farea_markers[0],
            markersize=6
        )
        # Plot farea0.2
        ax[1][2].plot(
            nangle,
            farea02_075[1][nangle],
            color=farea_colours[1],
            marker=farea_markers[1],
            markersize=6
        )
        # farea0.3
        ax[1][2].plot(
            nangle,
            farea03_075[1][nangle],
            color=farea_colours[2],
            marker=farea_markers[2],
            markersize=6
        )
        # farea0.4
        ax[1][2].plot(
            nangle,
            farea04_075[1][nangle],
            color=farea_colours[3],
            marker=farea_markers[3],
            markersize=6
        )
        # farea0.5
        ax[1][2].plot(
            nangle,
            farea05_075[1][nangle],
            color=farea_colours[4],
            marker=farea_markers[4],
            markersize=6
        )
        # farea0.6
        ax[1][2].plot(
            nangle,
            farea06_075[1][nangle],
            color=farea_colours[5],
            marker=farea_markers[5],
            markersize=6
        )
    #
    # And plot DETECTION PROBABILITY averages
    # For 052
    #
    ax[1][0].plot(
        [-1,6],[farea01_052_average[1],farea01_052_average[1]],'--',color=farea_colours[0]
    )
    ax[1][0].plot(
        [-1,6],[farea02_052_average[1],farea02_052_average[1]],'--',color=farea_colours[1]
    )
    ax[1][0].plot(
        [-1,6],[farea03_052_average[1],farea03_052_average[1]],'--',color=farea_colours[2]
    )
    ax[1][0].plot(
        [-1,6],[farea04_052_average[1],farea04_052_average[1]],'--',color=farea_colours[3]
    )
    ax[1][0].plot(
        [-1,6],[farea05_052_average[1],farea05_052_average[1]],'--',color=farea_colours[4]
    )
    ax[1][0].plot(
        [-1,6],[farea06_052_average[1],farea06_052_average[1]],'--',color=farea_colours[5]
    )
    #
    # For 074
    #
    ax[1][1].plot(
        [-1,6],[farea01_074_average[1],farea01_074_average[1]],'--',color=farea_colours[0]
    )
    ax[1][1].plot(
        [-1,6],[farea02_074_average[1],farea02_074_average[1]],'--',color=farea_colours[1]
    )
    ax[1][1].plot(
        [-1,6],[farea03_074_average[1],farea03_074_average[1]],'--',color=farea_colours[2]
    )
    ax[1][1].plot(
        [-1,6],[farea04_074_average[1],farea04_074_average[1]],'--',color=farea_colours[3]
    )
    ax[1][1].plot(
        [-1,6],[farea05_074_average[1],farea05_074_average[1]],'--',color=farea_colours[4]
    )
    ax[1][1].plot(
        [-1,6],[farea06_074_average[1],farea06_074_average[1]],'--',color=farea_colours[5]
    )
    #
    # For 075
    #
    ax[1][2].plot(
        [-1,6],[farea01_075_average[1],farea01_075_average[1]],'--',color=farea_colours[0]
    )
    ax[1][2].plot(
        [-1,6],[farea02_075_average[1],farea02_075_average[1]],'--',color=farea_colours[1]
    )
    ax[1][2].plot(
        [-1,6],[farea03_075_average[1],farea03_075_average[1]],'--',color=farea_colours[2]
    )
    ax[1][2].plot(
        [-1,6],[farea04_075_average[1],farea04_075_average[1]],'--',color=farea_colours[3]
    )
    ax[1][2].plot(
        [-1,6],[farea05_075_average[1],farea05_075_average[1]],'--',color=farea_colours[4]
    )
    ax[1][2].plot(
        [-1,6],[farea06_075_average[1],farea06_075_average[1]],'--',color=farea_colours[5]
    )
    #
    # PLOT LOWER LIMITS FOR THOSE WITH NO DETECTIONS
    #
    #    074:
    #    0.6: 180-0 (4), 270-0 (5)
    ax[0][1].plot(
        4,
        29.4,
        color=farea_colours[5],
        marker='^',
        markersize=6
    )
    ax[0][1].plot(
        5,
        29.4,
        color=farea_colours[5],
        marker='^',
        markersize=6
    )
    #
    #    075:
    #    0.4: 90-90
    ax[0][2].plot(
        2,
        75,
        color=farea_colours[3],
        marker='^',
        markersize=6
    )
    #    0.5: 0-0, (90-90), 90-270
    ax[0][2].plot(
        0,
        75,
        color=farea_colours[4],
        marker='^',
        markersize=6
    )
    ax[0][2].plot(
        3,
        75,
        color=farea_colours[4],
        marker='^',
        markersize=6
    )
    #    0.6: (0-0), (90-90), (90-270), 180-0, 270-0
    ax[0][2].plot(
        4,
        75,
        color=farea_colours[5],
        marker='^',
        markersize=6
    )
    ax[0][2].plot(
        5,
        75,
        color=farea_colours[5],
        marker='^',
        markersize=6
    )




    #
    # List the labels so that theres 1 per model.
    #
    labelpanel = 2
    for nn in range(len(farea_label)):
        ax[1][labelpanel].plot(
            -1,-1,
            color=farea_colours[nn],
            marker=farea_markers[nn],
            markersize=6,
            label=farea_label[nn]
        )
    ax[1][labelpanel].legend(
        #loc='upper left',
        title='Cloud area limit',
        fontsize=12
    )
    #
    # Modify xticklabels and title
    #
    for nmodel in range(Nmodels):
        for nrow in range(2):
            ax[nrow][nmodel].set_xticks([0,1,2,3,4,5]) 
            ax[nrow][nmodel].set_xticklabels(angles_label) 
        ax[0][nmodel].set_title(models_label[nmodel],fontsize=14)
        #ax[0][nmodel].set_xticklabels(['','','','','','']) 
        
    #
    # Set axislabels
    #
    ax[0][0].set_ylabel(r'Average period (yrs)', fontsize=18)
    ax[1][0].set_ylabel(r'Detection probability', fontsize=18)
    ax[1][1].set_xlabel(r'LOS-angle',fontsize=18)
    #
    # Save figure
    #
    fig.tight_layout()
    fig.savefig(
        'figs/periods_allmetrics.pdf', 
        facecolor='white',
        dpi=300
    )
    #fig.show()
#
# 
# Plot number of clouds per model and LOS angle
#
if plot_numbclouds == 'y':

    # Set various settings
    fig,ax = plt.subplots(
        Nangles,Nmodels,
        figsize=(13,14)
    )
    max_flux_contrast = 0.01
    fract_starareas = [0.1,0.3,0.5]
    fract_stararea_colours = ['r','g','b']
    #
    # Too messy with all farea and legend is too large. Doesnt say more to 
    # have all, might as well use the subset
    #fract_starareas = [0.1,0.2,0.3,0.4,0.5,0.6]
    #fract_stararea_colours = ['g','c','b','m','r','darkorange']
    #
    # Loop over models
    for nmodel,model in enumerate(models):

        # Load snapshot times for each model
        phasetimes = np.loadtxt(f'../r3dresults/{model}_nospikes/snapshot_yr.dat')[:,1]

        # Loop over fraction star area and load nblobs for all angles
        for nstararea,stararea in enumerate(fract_starareas):

            angles,nsnaps,nblobs,blob_areas = atf.load_imageblob_files(
                filepath=f'../r3dresults/{model}_nospikes/',
                fract_stararea = stararea
            )
            Nsnaps = len(nsnaps)
            Nangles = len(angles)

            # Plot for each model and each angle
            for nangle,angle in enumerate(angles_label):
                ax[nangle,nmodel].step(
                    phasetimes,
                    nblobs[:,nangle],
                    where='mid',
                    label=f'{stararea}',
                    color=fract_stararea_colours[nstararea]
                )
                # And set plot settings
                ax[nangle,nmodel].set_xlim(phasetimes[0],phasetimes[-1])
                ax[nangle,nmodel].set_ylim(0,5)
                ax[nangle,nmodel].tick_params(axis='both', which='major', labelsize=15)
                ax[nangle,nmodel].tick_params(axis='both', which='major', labelsize=15)
                ax[nangle,0].set_ylabel(
                    f'N clouds at {angle}',
                    fontsize=18
                )
        #
        ax[0,nmodel].set_title(
            f'{models_label[nmodel]}',
            fontsize=14
        )
    ax[-1,1].set_xlabel(
        'Simulation time (yrs)',
        fontsize=18
    )
    ax[-1,2].legend(
        #loc='upper left',
        title='Cloud area limit',
        fontsize=14
    )
    # 4 blobs is max

    fig.tight_layout()
    fig.savefig(
        'figs/nblobs_allangles.pdf',
        facecolor='white',
        dpi=300
    )
    #plt.show()


# Plot histograms over cloud
if plot_cloudareas == 'y':

    # Set various settings
    fig,ax = plt.subplots(
        1,Nmodels,
        figsize=(13,4)
    )
    Rstar = 1.65
    Rin = 2*Rstar
    Rout = 6*Rstar
    max_flux_contrast = 0.01
    fract_starareas = [0.1,0.2,0.3,0.4,0.5,0.6]
    fract_stararea_colours = ['g','c','b','m','r','darkorange']
    # Loop over models
    for nmodel,model in enumerate(models):

        # Declare/reset blob area list
        blob_areas_all = []

        # Load cloud areas
        angles,nsnaps,nblobs,blob_areas = atf.load_imageblob_files(
            filepath = f'../r3dresults/{model}_nospikes/',
            max_flux_contrast = 0.01,
            load_blobareas = 'y'
        )
        # Loop through all blob areas and put them in one huge 1D-list
        for blob_area in blob_areas:
            for blob in blob_area:
                blob_areas_all.append(blob)


        binareas, binranges, patches = ax[nmodel].hist(
            blob_areas_all,
            bins=int(np.round(np.max(blob_areas_all)*5)),
            histtype='step',
            align='mid',
            color='k',
            bottom=0.001,
            log=True
        )

        # Plot fractional area limit
        for nfract_stararea,fract_stararea in enumerate(fract_starareas):
            ax[nmodel].plot(
                [fract_stararea * np.pi*Rstar**2 , fract_stararea * np.pi*Rstar**2],
                [0.001,np.max(binareas)+10000],
                color=fract_stararea_colours[nfract_stararea],
                label = f'{fract_stararea}'
            )
        # Set settings
        ax[nmodel].tick_params(axis='both', which='major', labelsize=15)
        ax[nmodel].tick_params(axis='both', which='major', labelsize=15)
        ax[nmodel].set_ylim(0.5,np.max(binareas)+2000)
        ax[nmodel].set_xlim(-0.2,15)
        ax[nmodel].set_title(models_label[nmodel], fontsize=14)
    ax[-1].legend(
        title='Cloud area limit',
        fontsize=14
    )
    ax[1].set_xlabel(r'Cloud area (au$^2$)', fontsize=18)
    ax[0].set_ylabel('Number of clouds', fontsize=18)



    fig.tight_layout()
    fig.savefig(
        'figs/blobareas_histogram.pdf',
        facecolor='white',
        dpi=300
    )

# Plot "best of" random sample figures
if plot_bestrandomsample == 'y':
    # ie those with at least 1 cloud >0.1Astar and most clouds of those
    # 
    #    052
    #        297 : 46.898 yrs at i090_phi000 : 3.0   0.0   0.0   0.0   0.0   0.0
    #        297                 i090_phi000 : [1.18792849 0.74159698 0.01373328 1.08492891 1.60679345 0.81713 0.04806647 0.01029996 0.00343332 0.03089987]
    #    074
    #        413 : 109.641 yrs at i090_phi090 : 1.0   0.0   0.0   0.0   0.0   0.0
    #        413                 i090_phi090 : [0.67636391 0.03089987 1.0299958  0.61799748 0.08583298 0.05493311]
    #    075
    #        136 : 65.753 yrs at i090_phi270 : 1.0   0.0   0.0   0.0   0.0   0.0
    #        136                 i090_phi270 : [0.04806647 1.11582878 0.02059992 0.00343332 0.13733277]
    #
    imagepaths = [
        '../r3dresults/st28gm06n052_timedep_nospikes/297/',
        '../r3dresults/st28gm06n074_nospikes/413/',
        '../r3dresults/st28gm06n075_nospikes/136/',
    ]
    imagenames = [
        'image_i090_phi000_10um.out',
        'image_i090_phi090_10um.out',
        'image_i090_phi270_10um.out',
    ]
    imagetimes = [
        'Model-A: 46.90 yr',
        'Model-B: 109.64 yr',
        'Model-C: 65.75 yr',
    ]
    # Overwrite labels for shorter version
    models_label = [
        'Model-A',
        'Model-B',
        'Model-C',
    ]
    # Chosen snapshots from random sample
    nsnapshots = [
        297,413,136
    ]
    # Angles of chosen snapshots
    nangles = [
        1,2,3
    ]
    # Max number of snapshots for each model
    models_snapshots = [
        400,
        450,
        442,
    ]
    model_colours = [
        'b','r','g'
    ]
    model_symbols = [
        'd','s','o'
    ]
    # Set up image object and image scale
    fig,ax = plt.subplots(
        1,4,
        figsize=(13,3.5),
    )
    scale = [
        np.sqrt(1e-2),
        np.sqrt(9e-1),
    ]
    # And numbers for contour lines in images
    StarArea = np.pi * Rstar**2
    max_flux_contrast = 0.01
    exposure_limit = np.sqrt(240000*1e-6)
    #
    # Loop through the 3 different model selections
    for nimage in range(3):
        # Set all parameters
        model = models[nimage]
        nangle = nangles[nimage]
        angle_label = angles_label[nangle]
        nsnapshot = nsnapshots[nimage]
        model_snapshot = models_snapshots[nimage]
        imagepath = imagepaths[nimage]
        imagename = imagenames[nimage]
        imagetime = imagetimes[nimage]

        # Load stellar flux to extract stellar flux limits
        # in MJy and sqrt
        nodustpath = imagepath.replace('nospikes', 'nodust')
        temp1,temp2,stellarflux,temp3 = a3d.load_images(
            path=nodustpath,
            image=imagename
        )
        stellarflux /= StarArea
        Fstarlimit = np.sqrt(max_flux_contrast * stellarflux * 1e-6)

        # Create first 3 panels with the chosen images
        # Load image
        image2d,image2dlog,flux,axisplot = a3d.load_images(
            path=imagepath,
            image=imagename,
            distance=1
        )
        # Change to MJy and scale (same as original example images)
        # but with sqrt as well
        image2d = np.sqrt(
            image2d/1e6
        )
        # Plot image and save colourbar info
        ax[nimage].imshow(
            image2d, 
            origin='lower', 
            extent=axisplot, 
            vmin=scale[0],
            vmax=scale[1],
            cmap=plt.get_cmap('hot')
        )
        # Exposure limit
        ax[nimage].contour(
            image2d, 
            origin='lower', 
            extent=axisplot, 
            levels=[exposure_limit],
            linewidths=1,
            colors='cyan'
        )
        # Contrast limit line
        ax[nimage].contour(
            image2d, 
            origin='lower', 
            extent=axisplot, 
            levels=[Fstarlimit],
            linewidths=1,
            colors='lightgrey'
        )
        # Write time for each observation
        ax[nimage].set_title(
            f'{imagetime}, {angle_label}',
            fontsize = 12,
            loc='left'
        )
        # Create final panel with cloud sizes
        # load cloud sizes
        temp1,temp2,temp3,blob_areas,temp4 = atf.load_imageblob_files(
            filepath=f'../r3dresults/{model}_nospikes/'
        )
        # Extract chosen snapshots and angles
        row_counter = 1
        angle_counter = 0
        for blob_area in blob_areas:
            if angle_counter == nangle and row_counter == nsnapshot:
                # Print sizes of clouds for each model-selection
                print(
                    f'{model} {row_counter:03d} {angles[angle_counter]} : {blob_area}'
                )
                # Sort blob areas and plot, largest first
                blob_area.sort()
                blob_area = blob_area[::-1]
                nblobs = len(blob_area)
                # Save largest number of clouds
                if nimage == 0:
                    # Set initial max nblobs
                    nblobs_max = nblobs
                if nblobs > nblobs_max:
                    # Change max nblobs if there are more nblobs in later models
                    nblobs_max = nblobs
                # Plot cloud sizes in final panel
                ax[-1].plot(
                    list(range(1,nblobs+1)),blob_area,
                    model_symbols[nimage],
                    color=model_colours[nimage],
                    markersize=10,
                    label=models_label[nimage],
                )
            # Update counters for cloud size extraction
            row_counter += 1
            if row_counter > model_snapshot:
                row_counter = 1
                angle_counter += 1

        # Set general axis settings
        ax[nimage].tick_params(axis='both', which='major', labelsize=15)
        ax[nimage].tick_params(axis='both', which='major', labelsize=15)

    # Plot fractional limits, pixel size and largest possible size
    fract_starareas = [
        0.1,0.2,0.3,0.4
    ]
    for fract_limit in fract_starareas:
        ax[-1].plot(
            [0,nblobs_max+1],[fract_limit*np.pi*Rstar**2,fract_limit*np.pi*Rstar**2],
            'k:'
        )
    # And a line at zero
    ax[-1].plot(
        [0,nblobs_max+1],[0,0],
        'k--'
    )

    # Set figure settings
    ax[0].set_ylabel('Offset (au)', fontsize = 14)
    ax[1].set_xlabel('Offset (au)', fontsize = 14)
    ax[-1].set_ylabel(r'Cloud area (au$^2$)', fontsize = 14)
    ax[-1].set_xlabel('Cloud number', fontsize = 14)
    ax[-1].legend(
        loc='upper right',
        fontsize=12
    )
    ax[-1].tick_params(axis='both', which='major', labelsize=15)
    ax[-1].tick_params(axis='both', which='major', labelsize=15)
    ax[-1].set_xticks(list(range(1,nblobs_max+1)))
    ax[-1].set_xlim([0.8,nblobs_max+0.2])
    ax[-1].set_ylim([-0.2,2])
    ax[-1].set_box_aspect(1)

    # Save figure
    fig.tight_layout()
    fig.savefig(
        'figs/best_randomsample.pdf',
        facecolor='white',
        dpi=300
    )


# Plot 24*3 images for appendix, ie whole random sample
if plot_allrandomsample == 'y':

    # Set seed and random-generator-object
    rng = np.random.default_rng(
        seed=42
    )
    nsnap_start = 60
    models_snapshots = [
        400,
        450,
        442,
    ]
    # Loop over models here
    for nmodel in range(Nmodels):
        model_snapshot = models_snapshots[nmodel]
        model = models[nmodel]

        # Define snapshots and angles
        # Observational snapshots
        obs_snapshots = np.sort(rng.integers(
            nsnap_start,high=model_snapshot,size=24
        ))
        # Randomized obs-angle for each snapshot as well
        obs_angles = rng.integers(
            0,high=6,size=24
        )
        # Fix so that no snapshot-angle-combo appears twice
        # Predefine first step
        old_combo = [obs_snapshots[0],obs_angles[0]]
        counter = 0
        for nsnap,nangle in zip(obs_snapshots,obs_angles):
            # Compare next snapshot with previous one
            if counter > 0:
                if old_combo == [nsnap,nangle]:
                    # If theres an similarity just add 1... (works with seed=42)
                    obs_snapshots[counter] += 1
                # Update previous step
                old_combo = [nsnap,nangle]
            # Update step counter
            counter += 1

        # Set url settings
        path = f'../r3dresults/{model}_nospikes/'

        # Load snapshot-times for this model
        snapshot_file = np.loadtxt(
            f'{path}snapshot_yr.dat'
        )
        # Set up image object
        fig,ax = plt.subplots(
            6,4,
            figsize=(9,14),
        )
        # reset snapshot-number-counter
        snap_counter_lin = 0
        snap_counterX = 0
        snap_counterY = 0

        # Loop over observation snapshots and angles
        for n_obs,n_angle in zip(obs_snapshots,obs_angles):

            # Change to ints
            n_obs = n_obs.astype(int)
            n_angle = n_angle.astype(int)

            # Set angle name
            angle = angles[n_angle]

            # Extract correct time for each panel
            for snaptime in snapshot_file:
                if obs_snapshots[snap_counter_lin] == snaptime[0]:
                    obs_time = snaptime[1]


            # Load image
            image2d,image2dlog,flux,axisplot = a3d.load_images(
                path=f'{path}{n_obs:03d}/',
                image=f'image_{angle}_10um.out',
                distance=1
            )
            # Testplot image
            #
            # Change to MJy and to square-root stype
            image2d = np.sqrt(
                image2d/1e6
            )
            # Set scale same as for original example images for the two wavelengths
            scale = [
                np.sqrt(1e-2),np.sqrt(9e-1)
            ]
            # Plot image and save colourbar info
            ax[snap_counterY,snap_counterX].imshow(
                image2d, 
                origin='lower', 
                extent=axisplot, 
                vmin=scale[0],
                vmax=scale[1],
                cmap=plt.get_cmap('hot')
            )
            # Write time for each observation
            ax[snap_counterY,snap_counterX].set_title(
                f'{obs_time:.2f} yrs, {angles_label[n_angle]}',
                fontsize = 12,
                loc='left'
            )
            # Remove all ticks and tick labels
            ax[snap_counterY,snap_counterX].set_xticks([])
            ax[snap_counterY,snap_counterX].set_yticks([])
            #ax[snap_counterY,snap_counterX].axes.yaxis.set_ticklabels([])
            #ax[snap_counterY,snap_counterX].axes.xaxis.set_ticklabels([])

            # Update snapshot counters
            snap_counter_lin += 1
            snap_counterX += 1
            if snap_counterX == 4:
                snap_counterX = 0
                snap_counterY += 1


        # Save figure
        fig.tight_layout()
        fig.savefig(
            f'figs/all_randomsample_{model}.pdf',
            facecolor='white',
            dpi=300
        )






#####################################################################################
# OLD DISCARDED FIGS BELOW
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
    #fig.show()

# Plot comparisons with data as found at 
# suh2021, smiths2001, aavso.rg
# 
if plot_datacompare == 'y':
    # Plot with colour statistics
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
    # Y-axis to have colour ranges
    modelsymbol = [
        'd','s','o'
    ]
    #
    # Y-axis with model data in errorbars
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
    fig,ax = plt.figure(figsize=(6,4)), plt.axes()
    plotaxis = [0,6]
    datapositions = [1,2,3,4,5]
    #
    # Plot fields for observed statistics from Suh21
    # From Suh21-Fig9, 2 right panels two catalogues, 
    #
    irasdata = [
        0.34,           # min
        4.39-2.18,      # avr-std
        4.39,           # Average
        4.39+2.18,      # avr+std
        13.6,           # Max
    ]
    wisedata = [
        0.58,           # min
        2.15-0.91,      # avr-std
        2.15,           # Average
        2.15+0.91,      # avr+std
        8.02,           # Max
    ]
    # Plot model colour ranges
    for nmodel,modeldat in enumerate(modeldata):
        for nn in range(2):
            ax.plot(
                [datapositions[nmodel+2],datapositions[nmodel+2]],
                [modeldat[nn],modeldat[-nn-1]],
                'k',linestyle=linestyles[nn]
            )
            ax.plot(
                [datapositions[nmodel+2]-0.3,datapositions[nmodel+2]+0.3],
                [modeldat[nn],modeldat[nn]],
                'k'
            )
            ax.plot(
                [datapositions[nmodel+2]-0.3,datapositions[nmodel+2]+0.3],
                [modeldat[-nn-1],modeldat[-nn-1]],
                'k'
            )
        ax.plot(
            [datapositions[nmodel+2],datapositions[nmodel+2]],
            [modeldat[2],modeldat[2]],
            modelsymbol[nmodel],color='k',markersize=12
        )
    #
    # Plot Suh21-data
    # IRAS
    for nn in range(2):
        ax.plot(
            [datapositions[0],datapositions[0]],
            [irasdata[nn],irasdata[-nn-1]],
            'r',linestyle=linestyles[nn]
        )
        ax.plot(
            [datapositions[0]-0.3,datapositions[0]+0.3],
            [irasdata[nn],irasdata[nn]],
            'r'
        )
        ax.plot(
            [datapositions[0]-0.3,datapositions[0]+0.3],
            [irasdata[-nn-1],irasdata[-nn-1]],
            'r'
        )
    ax.plot(
        [datapositions[0],datapositions[0]],
        [irasdata[2],irasdata[2]],
        'p',color='r',markersize=12
    )
    ax.plot(
        [0,6],
        [irasdata[2],irasdata[2]],
        'r:'
    )
    # WISE
    for nn in range(2):
        ax.plot(
            [datapositions[1],datapositions[1]],
            [wisedata[nn],wisedata[-nn-1]],
            'b',linestyle=linestyles[nn]
        )
        ax.plot(
            [datapositions[1]-0.3,datapositions[1]+0.3],
            [wisedata[nn],wisedata[nn]],
            'b'
        )
        ax.plot(
            [datapositions[1]-0.3,datapositions[1]+0.3],
            [wisedata[-nn-1],wisedata[-nn-1]],
            'b'
        )
    ax.plot(
        [datapositions[1],datapositions[1]],
        [wisedata[2],wisedata[2]],
        'h',color='b',markersize=12
    )
    ax.plot(
        [0,6],
        [wisedata[2],wisedata[2]],
        'b:'
    )




    # Set xlabels and tick settings
    ax.set_xlim(0.5,5.5)
    ax.set_xticks(datapositions) 
    ax.set_xticklabels(
        ['IRAS', 'WISE', 'Model-A', 'Model-B', 'Model-C']
    ) 
    ax.set_ylabel('K[2.2]$-$W3[12]', fontsize=18)
    ax.tick_params(axis='both', which='major', labelsize=15)
    #
    # Final fig fix and save
    #
    fig.tight_layout()
    fig.savefig(
        'figs/data_compare.pdf',
        facecolor='white',
        dpi=300
    )
    #fig.show()

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
    fig.savefig(
        f"figs/{models[0].split('_')[0]}_sourceradius_{wavelength}um.pdf", 
        dpi=300, 
        facecolor="white"
    )

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
    #fig.show()


