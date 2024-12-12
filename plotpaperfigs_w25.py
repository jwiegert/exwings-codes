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
plot_052fluxdensity = 'y'

plot_rsourceevents = 'n'



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
####################################################################################
# Plot flux density at 2um over time for 052
#
if plot_052fluxdensity == 'y':

    # TODO ?





    # Set figure settings
    #fig, ax = plt.figure(num=f'F(t) of {wavelength} um', figsize=(6, 4)), plt.axes()
    fig,ax = plt.subplots(
        1,2,
        figsize=(12, 4)
    )
    ax[0].set_ylabel(rf'$F({wavelength}\,\mu$m$)$, MJy at 1\,pc',fontsize=18)
    for nn in range(2):
        ax[nn].set_xlabel(r'Sim. time (yrs)',fontsize=18)
        ax[nn].tick_params(axis='both', which='major', labelsize=15)
        ax[nn].set_yscale('log')



    # Path to model of choice
    path = '../r3dresults/st28gm06n052_timedep_nospikes/'
    model = '052'

    # Wavelength to plot
    wavelength = 2

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
    # extract flux at wavelength
    # and declare array for relative fluxes vs time
    average_seds = np.loadtxt(f'{path}average_sed.dat')
    average_flux = average_seds[wavelengthindex,1]
    relative_flux = np.zeros((len(angles),Nphases))



    # Loop through all phases and extract all flux densities at wavelength
    # And the relative flux density (to average flux density)
    for nangles,angle in enumerate(angles):
        for nphase,phase in enumerate(phases):

            wavelengths,spectrum = a3d.load_spectrum(
                path = f'{path}{phase}/spectrum_{angle}.out'
            )
            fluxall[nangles,nphase] = spectrum[wavelengthindex]*1e-6

            relative_flux[nangles,nphase] = fluxall[nangles,nphase]/average_flux


        # TODO
        # byt till samma färg som för 052 i SED-plotten


        # Plot each angles time dependent F in MJy
        ax[0].plot(phasetimes,fluxall[nangles,:],'lightgrey')

        # TODO
        # Plot relative flux densities of each angle and dashed lines for the limits
        # 1.0: 'k-'
        # 0.8: 'k--'
        # 0.6: 'k:'



    # Save average of each angles flux density at each time
    fluxaverage = []
    for nphase in range(Nphases):
        fluxaverage.append(np.mean(fluxall[:,nphase]))

    # and plot average flux density
    ax.plot(phasetimes,fluxaverage,'k')

    # and save figure
    fig.tight_layout()
    #plt.savefig(f'../r3dplots/{model}_fluxtime_{wavelength}um.pdf', dpi=300)
    fig.show()










#
#####################################################################################
# Plot average period of Rsource events per model and angle
if plot_rsourceevents == 'y':

    Nmodels = len(models_label)
    Nangles = len(angles)


    # Set figure objects    
    #fig,ax = plt.figure(
    #    figsize=(6,4)
    #), plt.axes()
    fig,ax = plt.subplots(
        3,1,
        figsize=(6,12)
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
    ax[1].set_ylabel(r'Average period \& event length (yrs)', fontsize=18)
    ax[2].set_xlabel(r'LOS-angle',fontsize=18)

    # Save figure
    fig.tight_layout()
    fig.savefig(
        'figs/periods_allmetrics.pdf', 
        facecolor='white',
        dpi=300
    )
    fig.show()
