# Plots various figures for first co5bold-r3d-paper
import matplotlib.pyplot as plt
import numpy as np
import re

from matplotlib import rc
from mpl_toolkits.axes_grid1 import make_axes_locatable

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

# 
plot_coboldgrid = 'n'

#
plot_opticalthickness = 'n'

# Grain properties
plot_grainsizehist_all = 'n'
plot_grainsizehist_one = 'n'
plot_grainsizeradius = 'n'
plot_absscat = 'n'

#
plot_temperatureradial = 'n' # Only cobold-T, no comparison, not used
plot_temperaturecompare = 'n'

# Plot SEDs
plot_seds_cobolddarwin = 'n'
plot_seds_point = 'n'
plot_represenativeseds = 'n'




# Plot various images
plot_images_obscured = 'y'



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
    fig.savefig('figs/grid-cutthrough.pdf', dpi=300, facecolor="white")
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
    fig.savefig('figs/optthick_los.pdf', dpi=300, facecolor="white")
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
    fig.savefig('figs/grainsize_hist.pdf', dpi=300, facecolor="white")

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
    fig.savefig('figs/grainsize_hist.pdf', dpi=300, facecolor="white")

    fig.show()


# Plot figure with grain sizes vs R
if plot_grainsizeradius == 'y':

    phases = ['186','190','198']
    fig, ax = plt.subplots(
        len(phases),1,
        figsize=(6,13)
    )
    a5d.plot_grainsizeradius(phases=phases,ax=ax)


    #size_path='../r3dresults/st28gm06n052_staranddust_1/grid_cellsizes.csv',


    ax[0].set_ylabel(fr'Grain sizes ($\mu$m), phase: {phases[0]}',fontsize=18)
    ax[1].set_ylabel(fr'Grain sizes ($\mu$m), phase: {phases[1]}',fontsize=18)
    ax[2].set_ylabel(fr'Grain sizes ($\mu$m), phase: {phases[2]}',fontsize=18)
    ax[0].set_xlabel(r'',fontsize=18)
    ax[1].set_xlabel(r'',fontsize=18)
    ax[2].set_xlabel(r'Distance (AU)',fontsize=18)
    ax[0].tick_params(axis='both', which='major', labelsize=15)
    ax[1].tick_params(axis='both', which='major', labelsize=15)
    ax[2].tick_params(axis='both', which='major', labelsize=15)

    fig.tight_layout()
    fig.savefig('figs/grainsizes_radius.pdf', dpi=300, facecolor="white")
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
    fig.savefig('figs/abs_scat_angle.pdf', dpi=300, facecolor="white")

    fig.show()

# ----------------------------------------------------------------
#
# Plot Figure with radial dependence of temperatures

if plot_temperatureradial == 'y':
    fig,ax, Tbins,Tstd,Tminmax,radial_range = a3d.plot_temperaturebins_radius(
        temperature_path='../r3dresults/st28gm06n052_staranddust/186/dust_temperature.dat',
        grid_path='../r3dresults/st28gm06n052_staranddust/grid_distances.csv',
        amr_path='../r3dresults/st28gm06n052_staranddust/amr_grid.inp',
        numb_specie=1
    )

    ax.set_ylabel(r'Gas \&  dust temperature (K)',fontsize=18)
    ax.set_xlabel(r'Distance (AU)',fontsize=18)

    ax.tick_params(axis='both', which='major', labelsize=15)
    fig.tight_layout()
    fig.savefig('figs/temperature_radial.pdf', dpi=300, facecolor="white")

    fig.show()


# Plot comparisons of temperature between co5bold and radmc3d
if plot_temperaturecompare == 'y':

    # Set up subplots
    fig,ax = plt.subplots(
        3,1,
        figsize=(6,13)
    )

    phase = 186


    # Load star's radius here
    Mstar,Rstar,Lstar = a5d.load_star_information(
        savpath='../co5bold_data/dst28gm06n052/st28gm06n052_186.sav',
        printoutput='n'
    )
    Rstar /= AUcm


    # Load cobold-T and create subplot
    T_c5d, Tstd_c5d, Tmax_c5d, Tmin_c5d, radial_range = a3d.plot_temperaturebins_radius(
        temperature_path=f'../r3dresults/st28gm06n052_staranddust_1/{phase}/dust_temperature.dat',
        grid_path='../r3dresults/st28gm06n052_staranddust_1/grid_distances.csv',
        amr_path='../r3dresults/st28gm06n052_staranddust_1/amr_grid.inp',
        numb_specie = 1,
        ax=ax[0]
    )
    ax[0].set_ylabel(r'$T_{\rm CO5BOLD}$ (K)',fontsize=18)
    ax[0].set_xlabel(r'')
    ax[0].set(xlim=(0,26))
    ax[0].tick_params(axis='both', which='major', labelsize=15)



    # Load R3d-temperatures and bin each grain size's tmeperature
    # and take averages, max-min, stds of these
    # Re-use radial range from c5d-plot since the grid is identical

    # Load coordintaes of R3D-cells and change to AU
    # extract size of grid cube (round up because this is the centre of cell coord)
    # and create binned radius-array
    cellcoords = a3d.load_griddistances(
        gridpath='../r3dresults/st28gm06n052_pointtemperature/grid_distances.csv',
        amrpath='../r3dresults/st28gm06n052_pointtemperature/amr_grid.inp'
    )
    cubesize = np.ceil(cellcoords[:,1].max()/AUcm )
    radii = cellcoords[:,0]/AUcm
    Nbins = 100
    radial_bins = np.linspace(0,radii.max(),Nbins+1)


    # Create temporary temperature arrays
    Nspecies = 10
    temperatures_avr = np.zeros((Nbins,Nspecies))
    temperatures_max = np.zeros((Nbins,Nspecies))
    temperatures_min = np.zeros((Nbins,Nspecies))
    temperatures_std = np.zeros((Nbins,Nspecies))

    # Bin each grain size temperature
    for nspecie in range(Nspecies):

        # Load each specie's R3D-temperature
        Ncells, Ntemp, temperatures = a3d.load_temperature(
            path=f'../r3dresults/st28gm06n052_pointtemperature/{phase}/dust_temperature.dat',
            numb_specie=nspecie+1
        )

        # Bin the temperatures and save in a Nbins*Nspecies arrays
        for nn in range(Nbins):
            ncells = np.where((radii >= radial_bins[nn]) & (radii < radial_bins[nn+1]))[0]
            temperatures_avr[nn,nspecie] = temperatures[ncells].mean()
            temperatures_max[nn,nspecie] = temperatures[ncells].max()
            temperatures_min[nn,nspecie] = temperatures[ncells].min()
            temperatures_std[nn,nspecie] = temperatures[ncells].std()


    # Save average/max-min/std's of each specie 
    # This saves the maximum standard deviation of eahc specie at each radial bin
    Tr3d_avr = np.zeros(Nbins)
    Tr3d_max = np.zeros(Nbins)
    Tr3d_min = np.zeros(Nbins)
    Tr3d_std = np.zeros(Nbins)

    for nn in range(Nbins):
        Tr3d_avr[nn] = temperatures_avr[nn,:].mean()
        Tr3d_max[nn] = temperatures_max[nn,:].max()
        Tr3d_min[nn] = temperatures_min[nn,:].min()
        Tr3d_std[nn] = temperatures_std[nn,:].max()


    ax[1].plot(radial_range,Tr3d_avr,'k')

    ax[1].fill_between(
        radial_range,
        Tr3d_min,
        Tr3d_max,
        color='b',
        alpha=0.2
    )

    ax[1].fill_between(
        radial_range,
        Tr3d_avr-Tr3d_std,
        Tr3d_avr+Tr3d_std,
        color='b',
        alpha=0.4
    )


    ax[1].set_ylabel(r'$T_{\rm RADMC-3D}$ (K)',fontsize=18)
    ax[1].set_xlabel(r'')
    ax[1].set(xlim=(0,26), ylim=(0,4000))
    ax[1].tick_params(axis='both', which='major', labelsize=15)

    # A plot with Tc5d / Tr3d
    ax[2].plot(radial_range,T_c5d/Tr3d_avr,'b')

    ax[2].fill_between(
        radial_range,
        (T_c5d-Tstd_c5d)/(Tr3d_avr-Tr3d_std),
        (T_c5d+Tstd_c5d)/(Tr3d_avr+Tr3d_std),
        color='b',
        alpha=0.4
    )
    ax[2].plot([1.65,1.65],[0,10],'r:')
    ax[2].set_ylabel(r'$T_{\rm CO5BOLD}$ / $T_{\rm RADMC-3D}$',fontsize=18)
    ax[2].set_xlabel(r'Distance (AU)',fontsize=18)
    ax[2].set(
        ylim=(1,3),
        xlim=(0,26)
    )
    ax[2].tick_params(axis='both', which='major', labelsize=15)

    # A vertical line at limit of grid
    for nn in range(3):
        # Grid cube
        ax[nn].plot([cubesize,cubesize],[0,4100],'k:')
        # And stellar radius
        ax[nn].plot([Rstar,Rstar],[0,4100],'r:',linewidth=1)


    # Show all plots
    plt.tight_layout()
    fig.show()

    #Save figure
    fig.savefig(f'figs/temperatures_{phase}.pdf', dpi=300, facecolor="white")




# ----------------------------------------------------------------
#
# Plot Figure subplots of co5bold-darwin Seds

if plot_seds_cobolddarwin == 'y':

    paths = [
        '../r3dresults/st28gm06n052_staranddust_nospikes/',
        '../r3dresults/st28gm06n052_darwinsource/',
    ]
    
    
    phases = ['186','190','198']
    spectra = [
        '/spectrum_i000_phi000.out',
        '/spectrum_i090_phi000.out',
        '/spectrum_i090_phi090.out',
        '/spectrum_i090_phi270.out',
        '/spectrum_i180_phi000.out',
        '/spectrum_i270_phi000.out'
    ]

    legendlist = [
        '0, 0',
        '90, 0',
        '90, 90',
        '90, 270',
        '180, 0',
        '270, 0',
    ]
    fig,ax = plt.subplots(len(paths)+1,len(phases), figsize = (12, 10))

    
    sedscobold = []
    sedcoboldnumbers = []
    
    sedsdarwin = []
    sedsdarwinnumbers = []

    # Extract SEDs and fill figure-axes-objects
    for nphase,phase in enumerate(phases):
        for npath,path in enumerate(paths):
            for nangle,spectrum in enumerate(spectra):

                wavelength,sed = a3d.load_spectrum(
                    path = path+phase+spectrum
                )
                
                if npath == 0:
                    sedscobold.append(sed)
                    sedcoboldnumbers.append(f'{nphase} {nangle}')
                    
                if npath == 1:
                    sedsdarwin.append(sed)
                    sedsdarwinnumbers.append(f'{nphase} {nangle}')


                ax[npath][nphase].plot(wavelength,sed,label = legendlist[nangle])

            ax[npath][nphase].set(xscale='log',yscale='log')
            ax[npath][nphase].set_xlim(5e-1,6e1)
            ax[npath][nphase].set_ylim(1e6,1.3e8)
            ax[npath][nphase].tick_params(axis='both', which='major', labelsize=15)



    # Compute chi2-numbers of each SED and add to plots
    for nn in range(len(sedcoboldnumbers)):
        if sedcoboldnumbers[nn] == sedsdarwinnumbers[nn]:
            
            #chisq_array, chiaq_reduced = a3d.compute_chisquare(
            #    simulation = np.array(sedsdarwin[nn]),
            #    observation = np.array(sedscobold[nn]),
            #    obssigma = np.array(sedscobold[nn])
            #)

            chisq_array = (np.array(sedsdarwin[nn]) - np.array(sedscobold[nn]))/np.array(sedscobold[nn])
            
            # Plot line (angle order should be the same as above)
            nphase = int(sedcoboldnumbers[nn][0])
            ax[2][nphase].plot(wavelength,chisq_array)
            ax[2][nphase].set(xscale='log',yscale='linear')
            ax[2][nphase].set_xlim(5e-1,6e1)
            #ax[2][nphase].set_ylim(0,3)
            ax[2][nphase].tick_params(axis='both', which='major', labelsize=15)
    for nn in range(3):
        ax[2][nn].plot([wavelength[0],wavelength[-1]],[0,0],'k:')
    
    ax[0][0].set_ylabel(r'$F({\rm CO5BOLD})$, Jy at 1 pc', fontsize=18)
    ax[1][0].set_ylabel(r'$F({\rm DARWIN})$, Jy at 1 pc', fontsize=18)
    ax[2][0].set_ylabel(
        r'$\frac{F({\rm DARWIN}) - F({\rm CO5BOLD})}{F({\rm CO5BOLD}}$', 
        fontsize=18
    )
    for nn in range(3):
        ax[-1][nn].set_xlabel(r'Wavelength ($\mu$m)', fontsize=18)
    ax[0][2].legend(title=r'$i$, $\phi$')

    fig.tight_layout()
    fig.show()
    
    #Save figure
    fig.savefig(f'figs/seds_all_cobold_darwin.pdf', dpi=300, facecolor="white")


# Plots SEDs of point-source-temperature
if plot_seds_point == 'y':

    paths = [
        '../r3dresults/st28gm06n052_pointtemperature/',
        '../r3dresults/st28gm06n052_staranddust_nospikes/',
    ]
    
    
    phases = ['186','190','198']
    spectra = [
        '/spectrum_i000_phi000.out',
        '/spectrum_i090_phi000.out',
        '/spectrum_i090_phi090.out',
        '/spectrum_i090_phi270.out',
        '/spectrum_i180_phi000.out',
        '/spectrum_i270_phi000.out'
    ]

    legendlist = [
        '0, 0',
        '90, 0',
        '90, 90',
        '90, 270',
        '180, 0',
        '270, 0',
    ]
    fig,ax = plt.subplots(len(paths),len(phases), figsize = (12, 6))


    # Loop through phases
    for nphase,phase in enumerate(phases):
        # Loop through spectra
        for nangle,spectrum in enumerate(spectra):

            # Load all point source seds
            wavelength,sed = a3d.load_spectrum(
                path = paths[0]+phase+spectrum
            )
            
            # Plot on first row
            ax[0][nphase].plot(wavelength,sed,label = legendlist[nangle])
        
        # Set settings for first row
        ax[0][nphase].set(xscale='log',yscale='log')
        ax[0][nphase].set_xlim(5e-1,6e1)
        ax[0][nphase].set_ylim(1e6,1.3e8)
        ax[0][nphase].tick_params(axis='both', which='major', labelsize=15)


    # Compute chi2-numbers of each SED and add to plots
    # (Yes I load the point-seds twice, whatev)
    for nphase,phase in enumerate(phases):
        for spectrum in spectra:

            # Load each spectra
            wavelength,sedpoint = a3d.load_spectrum(
                path = paths[0]+phase+spectrum
            )
            wavelength,sedcobold = a3d.load_spectrum(
                path = paths[1]+phase+spectrum
            )

            # Compute comparisons
            #chisq_array, chiaq_reduced = a3d.compute_chisquare(
            #    simulation = np.array(sedpoint),
            #    observation = np.array(sedcobold),
            #    obssigma = np.array(sedcobold)
            #)
            chisq_array = (np.array(sedpoint) - np.array(sedcobold))/np.array(sedcobold)
            
            # Plot line (angle order should be the same as above)
            ax[1][nphase].plot(wavelength,chisq_array)
        ax[1][nphase].plot([wavelength[0],wavelength[-1]],[0,0],'k:')
        ax[1][nphase].set(xscale='log',yscale='linear')
        ax[1][nphase].set_xlim(5e-1,6e1)
        ax[1][nphase].tick_params(axis='both', which='major', labelsize=15)

    ax[0][0].set_ylabel(r'$F({\rm point})$, Jy at 1 pc', fontsize=18)
    ax[1][0].set_ylabel(
        r'$\frac{F({\rm point}) - F({\rm CO5BOLD})}{F({\rm CO5BOLD})}$', 
        fontsize=18
    )
    for nn in range(3):
        ax[-1][nn].set_xlabel(r'Wavelength ($\mu$m)', fontsize=18)
    ax[0][2].legend(title=r'$i$, $\phi$')

    fig.tight_layout()
    fig.show()
    
    #Save figure
    fig.savefig(f'figs/seds_all_pointtemperature.pdf', dpi=300, facecolor="white")


# Plot 3 representative SEDer vid obskuration
if plot_represenativeseds == 'y':


    # Set up settings for plots
    fig,ax = plt.subplots(1,4, figsize = (12, 4))
    linestyles = [
        'b-',
        'k:',
        'k--'
    ]
    for nn in range(len(ax)):
        ax[nn].set(xscale='log',yscale='log')
        ax[nn].set_xlim(3e-1,6e1)
        ax[nn].set_ylim(5e5,1.3e8)
        ax[nn].tick_params(axis='both', which='major', labelsize=15)
        ax[nn].set_xlabel(r'Wavelength ($\mu$m)', fontsize=18)
    ax[0].set_ylabel(r'$F$ (Jy at 1 pc)', fontsize=18)
    

    # Links to files with "normal" SED
    imagefiles = [
        '../r3dresults/st28gm06n052_staranddust_nospikes/198/spectrum_i180_phi000.out',
        '../r3dresults/st28gm06n052_nodust/198/spectrum_i180_phi000.out',
        '../r3dresults/st28gm06n052_nostar/198/spectrum_i180_phi000.out'
    ]
    # Load and plot SEDs and wavelength grid
    for nn,imagefile in enumerate(imagefiles):
        wavelength,sed = a3d.load_spectrum(
            path = imagefile
        )
        ax[0].plot(wavelength,sed,linestyles[nn])
    ax[0].text(x=5e0,y=7e7,s=r'\noindent CO5BOLD\\$i = 180$, $\phi = 0$')

    # Links to obscured cobold-sed
    imagefiles = [
        '../r3dresults/st28gm06n052_staranddust_nospikes/198/spectrum_i090_phi090.out',
        '../r3dresults/st28gm06n052_nodust/198/spectrum_i090_phi090.out',
        '../r3dresults/st28gm06n052_nostar/198/spectrum_i090_phi090.out'
    ]
    # Load and plot SEDs and wavelength grid
    for nn,imagefile in enumerate(imagefiles):
        wavelength,sed = a3d.load_spectrum(
            path = imagefile
        )
        ax[1].plot(wavelength,sed,linestyles[nn])
    ax[1].text(x=5e0,y=7e7,s=r'\noindent CO5BOLD\\$i = 90$, $\phi = 90$')


    # Links to DARWIN-version of same LOS
    imagefiles = [
        '../r3dresults/st28gm06n052_darwinsource/198/spectrum_i090_phi090.out',
        '../r3dresults/st28gm06n052_darwinnodust/10/spectrum_i000_phi000.out',
        '../r3dresults/st28gm06n052_nostar/198/spectrum_i090_phi090.out'
    ]
    # Load and plot SEDs and wavelength grid
    for nn,imagefile in enumerate(imagefiles):
        wavelength,sed = a3d.load_spectrum(
            path = imagefile
        )
        ax[2].plot(wavelength,sed,linestyles[nn])
    ax[2].text(x=5e0,y=7e7,s=r'\noindent DARWIN\\$i = 90$, $\phi = 90$')


    # Links to r3d-temperature-version of same LOS
    imagefiles = [
        '../r3dresults/st28gm06n052_pointtemperature/198/spectrum_i090_phi090.out',
        '../r3dresults/st28gm06n052_pointtemperature/198/spectrum_i090_phi090_nodust.out',
        '../r3dresults/st28gm06n052_pointtemperature/198/spectrum_i090_phi090_nostar.out',
    ]
    # Load and plot SEDs and wavelength grid
    for nn,imagefile in enumerate(imagefiles):
        wavelength,sed = a3d.load_spectrum(
            path = imagefile
        )
        ax[3].plot(wavelength,sed,linestyles[nn])
    ax[3].text(x=5e0,y=7e7,s=r'\noindent Point source\\$i = 90$, $\phi = 90$')



    fig.tight_layout()
    fig.show()

    #Save figure
    fig.savefig(f'figs/seds_obscuredexamples.pdf', dpi=300, facecolor="white")


# ----------------------------------------------------------------------
# Plot various images

if plot_images_obscured == 'y':

    imagelist=[
        '../r3dresults/st28gm06n052_staranddust_1/198/image_i090_phi090_1um.out',
        '../r3dresults/st28gm06n052_staranddust_1/198/image_i090_phi090_2um.out',
        '../r3dresults/st28gm06n052_staranddust_1/198/image_i090_phi090_5um.out',
    ]
    distance=1

    # Number of plots
    Nplots = len(imagelist)

    # Set fig-and-axis settings for subplots
    fig, ax = plt.subplots(
        Nplots,1,
        figsize = (5,13)
    )
    #        dpi = 300, 

    # Load image data and save in various lists
    for nn,image in enumerate(imagelist):

        # Extract path and imagename from image
        imagestrings = re.split('/', image)
        path = f'{imagestrings[0]}/{imagestrings[1]}/'
        modelname = imagestrings[2]
        phase = imagestrings[3]
        imagefilename = imagestrings[4]

        # extract inclination and wavelength
        imagestrings = re.split('_', imagefilename)
        incl = imagestrings[1][1:]
        phi = imagestrings[2][3:]
        wavelengthum = imagestrings[3][:-6]

        # Remove 0 if first character in incl and phi
        incl = (incl[1:] if incl.startswith('0') else incl)
        phi = (phi[1:] if phi.startswith('0') else phi)

        # Load data
        image2d,image2dlog,flux,axisplot = a3d.load_images(
            path=f'{path}{modelname}/{phase}',
            image=imagefilename,
            distance=distance
        )

        # Save linear data in list (list of arrays)
        imagedata = image2d

        # Plot image at spot nn, set title and axis labels
        im0 = ax[nn].imshow(
            imagedata, 
            origin='lower', extent=axisplot, 
            cmap=plt.get_cmap('hot')
        )
        ax[nn].set(
            xlim=(-7,7),
            ylim=(-7,7)
        )
        ax[nn].set_ylabel(rf'Offset (AU), $\lambda =$ {wavelengthum} $\mu$m', fontsize=18)
        ax[nn].tick_params(axis='both', which='major', labelsize=15)
    ax[-1].set_xlabel('Offset (AU)', fontsize=18)
    


    # Change figure size
    fig.tight_layout()
    #fig.show()

    #Save figure
    fig.savefig(f'figs/images_obscuredexamples.pdf', dpi=300, facecolor="white")



