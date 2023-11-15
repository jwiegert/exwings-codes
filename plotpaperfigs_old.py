# Plots various figures for first co5bold-r3d-paper
import matplotlib.pyplot as plt
import numpy as np
import scipy.ndimage
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
Lsol = 3.828e26 # W
Rsolmeter = 6.955e8 # m
cubesize = 222757675648155.62/AUcm # Yes, hardcoded to large grid, change if needed for 
                                   # other grids, ie, from amr_grid, first coordinate 
                                   # is courner coordinate of first base cell
radian = 206264800 # milliasec
#baselineVLTI = 130.23 # metres
baselineVLTI = 201.92 # metres
diameterJWST = 6.5 # metres

# Dust model snapshots labels
phasetimes = [
    '29.31 yrs',
    '29.95 yrs',
    '31.21 yrs'
]


# Plot choices

# 
plot_coboldgrid = 'n'

#
plot_opticalthickness = 'n'

# Grain properties
plot_grainsizehist_one = 'n'
plot_grainsizeradius = 'n'
plot_absscat = 'y'

#
plot_temperatureradial = 'n' # Only cobold-T, no comparison, not used
plot_temperaturecompare = 'n'

# Plot SEDs
plot_seds_cobolddarwin = 'n'    # Old fig
plot_seds_point = 'n'           # Old fig
plot_seds_obscured = 'n'
plot_seds_cobold = 'n'
plot_seds_darwinpoint = 'n'




# Plot various images
plot_images_examples = 'n'
plot_images_darwinpoint = 'n'
plot_images_obscured = 'n'
plot_images_convolved_jwst = 'n'
plot_images_convolved_vlti = 'n'



# Observables
compute_luminosities = 'n'
compute_tenmicronfluxdensities = 'n'
measuredustcloudflux = 'n'
plot_resolutiondistance = 'n'
check_smoothedimage_radius = 'y'
plot_smoothedimage_radius = 'n'

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

    ax.set_xlabel(r'Distance along LOS (au)',fontsize=18)
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


# This plots histogram of mass of grain size bins
if plot_grainsizehist_one == 'y':

    fig,ax = a5d.plot_grainsizemass_histogram(
        model='st28gm06n052_staranddust_1',
        phases=[186,190,198],
        phaselabels=phasetimes
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
    for nn in range(len(phases)):
        ax[nn].set_ylabel(fr'Grain sizes ($\mu$m), $t_{nn+1} = $\,{phasetimes[nn]}',fontsize=18)
        ax[nn].tick_params(axis='both', which='major', labelsize=15)
    ax[0].set_xlabel(r'',fontsize=18)
    ax[1].set_xlabel(r'',fontsize=18)
    ax[2].set_xlabel(r'Distance (au)',fontsize=18)

    fig.tight_layout()
    fig.savefig('figs/grainsizes_radius.pdf', dpi=300, facecolor="white")
    fig.show()


# ----------------------------------------------------------------
#
# Plot Figure absorption and scattering, and angles

if plot_absscat == 'y':

    fig,ax = a3d.plot_allkappa(
        path='../r3dresults/opacities_mie_st28gm06n052/'
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
    ax.set_xlabel(r'Distance (au)',fontsize=18)

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
    ax[0].set_ylabel(r'$T({\rm CO5BOLD})$, K',fontsize=18)
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
        color='r',
        alpha=0.2
    )

    ax[1].fill_between(
        radial_range,
        Tr3d_avr-Tr3d_std,
        Tr3d_avr+Tr3d_std,
        color='r',
        alpha=0.4
    )


    ax[1].set_ylabel(r'$T({\rm point})$, K',fontsize=18)
    ax[1].set_xlabel(r'')
    ax[1].set(xlim=(0,26), ylim=(0,4000))
    ax[1].tick_params(axis='both', which='major', labelsize=15)

    # A plot with Tc5d or Tr3d divided by T_theory
    Teff = 2800   # Table1
    pindex = -0.9 # Bladh 2012
    temperature_theory = Teff * (Rstar/(2*radial_range))**(2/(4+pindex))

    legendlist = [
        r'$T($CO5BOLD$)$',
        r'$T($point$)$'
    ]

    ax[2].plot(radial_range,T_c5d/temperature_theory,'b',linewidth=2,label = legendlist[0])
    ax[2].plot(radial_range,Tr3d_avr/temperature_theory,'r',linewidth=2,label = legendlist[1])
    ax[2].legend(fontsize=13)

    ax[2].fill_between(
        radial_range,
        (T_c5d-Tstd_c5d)/temperature_theory,
        (T_c5d+Tstd_c5d)/temperature_theory,
        color='b',
        alpha=0.4
    )

    ax[2].fill_between(
        radial_range,
        (Tr3d_avr-Tr3d_std)/temperature_theory,
        (Tr3d_avr+Tr3d_std)/temperature_theory,
        color='r',
        alpha=0.4
    )

    ax[2].set_ylabel(r'$T({\rm simulated})$ / $T({\rm theory})$',fontsize=18)
    ax[2].set_xlabel(r'Distance (au)',fontsize=18)
    ax[2].set(
        ylim=(0.5,2),
        xlim=(0,26)
    )
    ax[2].tick_params(axis='both', which='major', labelsize=15)


    # A vertical line at limit of grid for all subplots
    for nn in range(3):
        # Grid cube
        ax[nn].plot([cubesize,cubesize],[0,4100],'k:',linewidth=1)
        # And stellar radius
        ax[nn].plot([Rstar,Rstar],[0,4100],'r:',linewidth=1)

    # Horisontal line for 1
    ax[2].plot([radial_range[0],radial_range[-1]],[1,1],'g:')

    # Show all plots
    plt.tight_layout()
    fig.show()

    #Save figure
    fig.savefig(f'figs/temperatures_{phase}.pdf', dpi=300, facecolor="white")

    # Also plot Ttheory by itself for comparison
    #plt.figure('Ttheory')
    #plt.plot(radial_range,temperature_theory)
    #plt.show()



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
            print(path)
            for nangle,spectrum in enumerate(spectra):

                wavelength,sed = a3d.load_spectrum(
                    path = path+phase+spectrum
                )
                
                # Print all sed-luminosities
                luminosity = a3d.compute_luminosity(
                    wavelengths=wavelength,
                    spectrum=sed
                )
                print(f'{phase}{spectrum}: {luminosity/Lsol}')


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
        # Line for the 1-1 correspondance
        ax[2][nn].plot([wavelength[0],wavelength[-1]],[0,0],'k:')
        # Phase number as title on top-plots
        ax[0][nn].set_title(rf'$t_{nn+1} = $\,{phasetimes[nn]}', fontsize=16)

    ax[0][0].set_ylabel(r'$F_\nu({\rm CO5BOLD})$, Jy at 1 pc', fontsize=18)
    ax[1][0].set_ylabel(r'$F_\nu({\rm DARWIN})$, Jy at 1 pc', fontsize=18)
    ax[2][0].set_ylabel(
        r'$\frac{F({\rm DARWIN}) - F({\rm CO5BOLD})}{F({\rm CO5BOLD}}$', 
        fontsize=18
    )
    for nn in range(3):
        ax[-1][nn].set_xlabel(r'Wavelength ($\mu$m)', fontsize=18)
    ax[0][2].legend(title=r'$i$, $\phi$', fontsize=13)

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

    print(paths[0])
    # Loop through phases
    for nphase,phase in enumerate(phases):
        # Loop through spectra
        for nangle,spectrum in enumerate(spectra):

            # Load all point source seds
            wavelength,sed = a3d.load_spectrum(
                path = paths[0]+phase+spectrum
            )
            
            # Print all sed-luminosities
            luminosity = a3d.compute_luminosity(
                wavelengths=wavelength,
                spectrum=sed
            )
            print(f'{phase}{spectrum}: {luminosity/Lsol}')

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

    ax[0][0].set_ylabel(r'$F_\nu({\rm point})$, Jy at 1 pc', fontsize=18)
    ax[1][0].set_ylabel(
        r'$\frac{F({\rm point}) - F({\rm CO5BOLD})}{F({\rm CO5BOLD})}$', 
        fontsize=18
    )
    for nn in range(3):
        ax[-1][nn].set_xlabel(r'Wavelength ($\mu$m)', fontsize=18)
        ax[0][nn].set_title(rf'$t_{nn+1} = $\,{phasetimes[nn]}', fontsize=16)
    ax[0][2].legend(title=r'$i$, $\phi$', fontsize=13)

    fig.tight_layout()
    fig.show()
    
    #Save figure
    fig.savefig(f'figs/seds_all_pointtemperature.pdf', dpi=300, facecolor="white")


# Plot 3 representative SEDer vid obskuration
if plot_seds_obscured == 'y':


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
    ax[0].set_ylabel(r'$F_\nu$ (Jy at 1 pc)', fontsize=18)
    

    # Links to files with SED from opposite direction
    imagefiles = [
        '../r3dresults/st28gm06n052_staranddust_nospikes/198/spectrum_i090_phi270.out',
        '../r3dresults/st28gm06n052_nodust/198/spectrum_i090_phi270.out',
        '../r3dresults/st28gm06n052_nostar/198/spectrum_i090_phi270.out'
    ]
    ax[0].text(x=5e0,y=7e7,s=r'\noindent CO5BOLD\\$i = 90$, $\phi = 270$')


    # Load and plot SEDs and wavelength grid
    for nn,imagefile in enumerate(imagefiles):
        wavelength,sed = a3d.load_spectrum(
            path = imagefile
        )
        ax[0].plot(wavelength,sed,linestyles[nn])

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


# Plot all SEDs of CO5BOLD data
if plot_seds_cobold == 'y':

    path = '../r3dresults/st28gm06n052_staranddust_nospikes/'
    phases = [
        '186','190','198'
    ]
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
    fig,ax = plt.subplots(1,len(phases), figsize = (13, 4.5))

    sedscobold = []
    sedcoboldnumbers = []
    
    # Extract SEDs and fill figure-axes-objects
    for nphase,phase in enumerate(phases):
        for nangle,spectrum in enumerate(spectra):

            wavelength,sed = a3d.load_spectrum(
                path = path+phase+spectrum
            )
            
            # Print all sed-luminosities
            luminosity = a3d.compute_luminosity(
                wavelengths=wavelength,
                spectrum=sed
            )
            print(f'{phase}{spectrum}: {luminosity/Lsol}')

            sedscobold.append(sed)
            sedcoboldnumbers.append(f'{nphase} {nangle}')
                
            ax[nphase].plot(wavelength,sed,label = legendlist[nangle])

        ax[nphase].set(xscale='log',yscale='log')
        ax[nphase].set_xlim(5e-1,6e1)
        ax[nphase].set_ylim(1e6,1.3e8)
        ax[nphase].tick_params(axis='both', which='major', labelsize=15)
        ax[nphase].set_xlabel(r'Wavelength ($\mu$m)', fontsize=18)    
        ax[nphase].set_title(rf'$t_{nphase+1} = $\,{phasetimes[nphase]}', fontsize=16)

    # Ylabel and legend with angles
    ax[0].set_ylabel(r'$F_\nu({\rm CO5BOLD})$, Jy at 1 pc', fontsize=18)
    ax[2].legend(title=r'$i$, $\phi$', fontsize=13)

    fig.tight_layout()
    fig.show()
    
    #Save figure
    fig.savefig(f'figs/seds_all_cobold.pdf', dpi=300, facecolor="white")










# Plot t3-SEDs of Darwin and Point-source-data
if plot_seds_darwinpoint == 'y':
    
    paths = [
        '../r3dresults/st28gm06n052_darwinsource/',
        '../r3dresults/st28gm06n052_pointtemperature/'
    ]
    phase = '198'
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
    fig,ax = plt.subplots(len(paths),1, figsize = (5, 8))
    
    sedsdarwin = []
    sedspoint = []


    # Extract SEDs and fill figure-axes-objects
    for npath,path in enumerate(paths):
        print(path)
        for nangle,spectrum in enumerate(spectra):

            wavelength,sed = a3d.load_spectrum(
                path = path+phase+spectrum
            )
            
            # Print all sed-luminosities
            luminosity = a3d.compute_luminosity(
                wavelengths=wavelength,
                spectrum=sed
            )
            print(f'{phase}{spectrum}: {luminosity/Lsol}')

            ax[npath].plot(wavelength,sed,label = legendlist[nangle])

        ax[npath].set(xscale='log',yscale='log')
        ax[npath].set_xlim(5e-1,6e1)
        ax[npath].set_ylim(1e6,1.3e8)
        ax[npath].tick_params(axis='both', which='major', labelsize=15)

    # Time stamp above
    ax[0].set_title(rf'$t_3 = $\,{phasetimes[2]}', fontsize=16)


    ax[0].set_ylabel(r'$F_\nu({\rm DARWIN})$, Jy at 1 pc', fontsize=18)
    ax[1].set_ylabel(r'$F_\nu({\rm point})$, Jy at 1 pc', fontsize=18)
    ax[1].set_xlabel(r'Wavelength ($\mu$m)', fontsize=18)
    ax[1].legend(title=r'$i$, $\phi$', fontsize=13)

    fig.tight_layout()
    fig.show()
    
    #Save figure
    fig.savefig(f'figs/seds_darwin_point.pdf', dpi=300, facecolor="white")








# ----------------------------------------------------------------------
# Plot various images

if plot_images_examples == 'y':
    # Plot 3 images at one walanvength, one per phase, lin and log

    # Chose wavelength here
    wavelengthum = 10

    distance = 1
    imagelist = [
        f'../r3dresults/st28gm06n052_staranddust_1/186/image_i000_phi000_{wavelengthum}um.out',
        f'../r3dresults/st28gm06n052_staranddust_1/190/image_i000_phi000_{wavelengthum}um.out',
        f'../r3dresults/st28gm06n052_staranddust_2/198/image_i000_phi000_{wavelengthum}um.out'
    ]
    
    # Number of plots
    Nplots = len(imagelist)

    # Set fig-and-axis settings for subplots
    fig, ax = plt.subplots(
        2,Nplots,
        figsize = (12,7.3),
    )

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

        # Load data
        image2d,image2dlog,flux,axisplot = a3d.load_images(
            path=f'{path}{modelname}/{phase}',
            image=imagefilename,
            distance=distance
        )
        # Change to MJy per asec2
        image2d = image2d*1e-6

        # Limits per wavelength suggestions:
        #  1um: 0 -> 3   (lin), -10 -> 6 (log)
        # 10um: 0 -> 1.5 (lin),  1  -> 6 (log)


        imlin = ax[0][nn].imshow(
            image2d, 
            origin='lower', extent=axisplot, 
            cmap=plt.get_cmap('hot'),
            vmin=0,vmax=1.5
        )
        ax[0][nn].set_title(rf'$t_{nn+1} = $\,{phasetimes[nn]}', fontsize=15)
        ax[0][nn].tick_params(axis='both', which='major', labelsize=15)

        imlog = ax[1][nn].imshow(
            image2dlog, 
            origin='lower', extent=axisplot, 
            cmap=plt.get_cmap('hot'),
            vmin=1,vmax=6
        )
        ax[1][nn].set_xlabel('Offset (au)',fontsize=18)
        ax[1][nn].tick_params(axis='both', which='major', labelsize=15)

    ax[0][0].set_ylabel('Offset (au)',fontsize=18)
    ax[1][0].set_ylabel('Offset (au)',fontsize=18)

    # In case I want to remove "inbetween" axis labels
    #ax[0][1].axes.get_yaxis().set_visible(False)
    #ax[0][2].axes.get_yaxis().set_visible(False)

    #ax[0][0].axes.get_xaxis().set_visible(False)    
    #ax[0][1].axes.get_xaxis().set_visible(False)
    #ax[0][2].axes.get_xaxis().set_visible(False)

    #ax[1][1].axes.get_yaxis().set_visible(False)
    #ax[1][2].axes.get_yaxis().set_visible(False)




    # Set colour bar settings and label, first row (linear)
    divider = make_axes_locatable(ax[0][-1])
    cax = divider.append_axes(
        'right', 
        size='5%', 
        pad=0.05
    )
    cb0 = plt.colorbar(
        imlin, 
        cax=cax, 
        orientation = 'vertical',
    )
    cb0.set_label(
        label = rf'$F_\nu$(MJy/asec$^2$) at {wavelengthum} $\mu$m \& {distance} pc',fontsize= 15
    )
    cb0.ax.tick_params(labelsize=15)

    # Set colour bar settings and label, second row (logarithmic)
    divider = make_axes_locatable(ax[1][-1])
    cax = divider.append_axes(
        'right', 
        size='5%', 
        pad=0.05
    )
    cb0 = plt.colorbar(
        imlog, 
        cax=cax, 
        orientation = 'vertical'
    )
    cb0.set_label(
        label = rf'$\log F$(Jy/asec$^2$) at {wavelengthum} $\mu$m \& {distance} pc',fontsize= 15
    )
    cb0.ax.tick_params(labelsize=15)

    fig.tight_layout()
    #fig.show()

    #Save figure
    fig.savefig(f'figs/images_{wavelengthum}umexamples.pdf', dpi=300, facecolor="white")



if plot_images_darwinpoint == 'y':
    # Plot Darwin-and-pointsource example images at 10um, of 190-phase

    imagelist=[
        '../r3dresults/st28gm06n052_darwinsource/190/image_i000_phi000_10um.out',
        '../r3dresults/st28gm06n052_pointtemperature/190/image_i000_phi000_10um.out'
    ]
    distance=1

    # Number of plots
    Nplots = len(imagelist)

    # Set fig-and-axis settings for subplots
    fig, ax = plt.subplots(
        Nplots,1,
        figsize = (5.5,8.3),
    )

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

        # Load data
        image2d,image2dlog,flux,axisplot = a3d.load_images(
            path=f'{path}{modelname}/{phase}',
            image=imagefilename,
            distance=distance
        )
        # Change to MJy per asec2
        image2d = image2d*1e-6

        imlin = ax[nn].imshow(
            image2d, 
            origin='lower', extent=axisplot, 
            cmap=plt.get_cmap('hot'),
            vmin=0,vmax=0.7
        )
        ax[nn].tick_params(axis='both', which='major', labelsize=15)
    
    ax[0].set_title(r'M2n315u6 \&\ dust', fontsize=15)
    ax[1].set_title(r'Point source \&\ dust', fontsize=15)
    ax[1].set_xlabel('Offset (au)', fontsize=18)
    ax[0].set_ylabel('Offset (au)', fontsize=18)
    ax[1].set_ylabel('Offset (au)', fontsize=18)


    # Set colour bar settings and label
    #divider = make_axes_locatable(ax[1])
    #cax = divider.append_axes(
    #    'right', 
    #    size='5%', 
    #    pad=0.05,
    #)
    # [Xcoord, Ycoord, Width, Height]
    cax = fig.add_axes([0.83, 0.17, 0.03, 0.7])
    cb0 = plt.colorbar(
        imlin, 
        cax=cax, 
        orientation = 'vertical', 
    )
    #shrink=0.6,pad=0.15
    cb0.set_label(
        label = rf'$F_\nu$(MJy/asec$^2$) at {wavelengthum} $\mu$m \& {distance} pc', fontsize=15
    )
    cb0.ax.tick_params(labelsize=15)

    fig.tight_layout()
    #fig.show()

    # Save figure
    fig.savefig(f'figs/images_10umdarwinpoint.pdf', dpi=300, facecolor="white")



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
        # Use default vertical scale
        ax[nn].imshow(
            imagedata, 
            origin='lower', 
            vmin= 0,
            extent=axisplot, 
            cmap=plt.get_cmap('hot')
        )
        ax[nn].set(
            xlim=(-7,7),
            ylim=(-7,7)
        )
        ax[nn].set_ylabel(rf'Offset (au), $\lambda =$ {wavelengthum} $\mu$m', fontsize=18)
        ax[nn].tick_params(axis='both', which='major', labelsize=15)
    ax[-1].set_xlabel('Offset (au)', fontsize=18)
    

    fig.tight_layout()
    #fig.show()

    #Save figure
    fig.savefig(f'figs/images_obscuredexamples.pdf', dpi=300, facecolor="white")


if plot_images_convolved_jwst == 'y':
    # Plot figure with convolved images, ie as observed
    # This one with 3 images with JWST, all phases, one wavelength

    distanceJWST = 40 # parsec


    imagelist_1um = [
        '../r3dresults/st28gm06n052_staranddust_1/186/image_i000_phi000_1um.out',
        '../r3dresults/st28gm06n052_staranddust_1/190/image_i000_phi000_1um.out',
        '../r3dresults/st28gm06n052_staranddust_1/198/image_i000_phi000_1um.out'
    ]


    # JWST vid 1um och 40pc
    #
    # Full-width at half-maximum: 1.22*lambda/D
    fwhmJWST_1um = 1.22 * 1e-6 / diameterJWST * radian
    sigmaJWST_1um = fwhmJWST_1um/2.355
    
    # Initialise fig-ax
    figJWST_1um, axJWST_1um = plt.subplots(
        1,3, 
        figsize = (12,5.8),
        num='JWST at 1um 40pc'
    )

    # Loop through and plot each image
    for nn,image in enumerate(imagelist_1um):
        distance = distanceJWST

        # Extract path and imagename from image
        imagestrings = re.split('/', image)

        modelname = imagestrings[2]
        imagefilename = imagestrings[4]
        phase = imagestrings[3]

        path = f'{imagestrings[0]}/{imagestrings[1]}/{modelname}/{phase}'

        # Load image
        image2d,image2dlog,flux,axisplot = a3d.load_images(
            path=path,
            image=imagefilename,
            distance=distance
        )
        # Change to Jy/mas2
        image2d = image2d * 1e-6

        # Extract props and compute distance-dependant scales
        Npix = np.shape(image2d)[0] # Number of pixels along one side
        auperpixel = 2*axisplot[0]/Npix  # number of au per pixel
        masperpixel = auperpixel/distance * 1000  # number of mas per pixel
        size_milliasec = masperpixel * Npix # size of image in mas

        # Image axis-limits
        axisplotmilliasec = [0.5*size_milliasec,-0.5*size_milliasec,-0.5*size_milliasec,0.5*size_milliasec]

        # Change sigma to order in number of pixels
        sigmaJWST_1um_pixel = sigmaJWST_1um / masperpixel

        # Convolve with Gaussian filter
        image2d = scipy.ndimage.gaussian_filter(
            image2d,
            sigma = sigmaJWST_1um_pixel
        )

        # Plot image
        imJWST_1um = axJWST_1um[nn].imshow(
            image2d, 
            origin='lower', extent=axisplotmilliasec, 
            cmap=plt.get_cmap('hot'),
            vmin=0, vmax=1.5
        )
        axJWST_1um[nn].tick_params(axis='both', which='major', labelsize=15)
        axJWST_1um[nn].set_xlabel('Offset (mas)',fontsize=18)

        # FWHM-cirle to show beam
        axJWST_1um[nn].add_patch(
            plt.Circle((250,-250), 0.5*fwhmJWST_1um, color='cyan', fill=False)
        )
        # Star-radius-circle to show size of stellar disc
        axJWST_1um[nn].add_patch(
            plt.Circle((0,0), 1.65/distance*1000, color='lime', fill=False, linestyle=':')
        )
        # Phase in titles
        axJWST_1um[nn].set_title(f'{phase}', fontsize=15)


    # Set ylabel for first plot only
    axJWST_1um[0].set_ylabel('Offset (mas)',fontsize=18)

    # Remove label texts for yaxis on the rest
    axJWST_1um[1].axes.get_yaxis().set_visible(False)
    axJWST_1um[2].axes.get_yaxis().set_visible(False)

    # Tight layout before colourbar
    figJWST_1um.tight_layout()

    # Set colour bar settings and label
    cax = figJWST_1um.add_axes([0.33, 0.95, 0.4, 0.03])
    cb0 = plt.colorbar(
        imJWST_1um, 
        cax=cax, 
        orientation = 'horizontal', 
    )
    cb0.set_label(
        label = rf'$F_\nu$(Jy/mas$^2$) at 1\,$\mu$m \& {distanceJWST}\,pc', fontsize=15
    )
    cb0.ax.tick_params(labelsize=15)

    # Show and save figure
    #figJWST_1um.show()
    figJWST_1um.savefig(f'figs/images_JWST_1um40pc.pdf', dpi=300, facecolor="white")



if plot_images_convolved_vlti == 'y':
    #
    # Plot figure with convolved images, ie as observed
    # This one with 6 images with VLTI, all phases, three wavelength
    # VLTI at 10um and 200pc
    # VLTI at 3.5um and 200pc
    # VLTI at  1um and 200pc
    #
    distanceVLTI = 200 # parsec
    shortwavelength = 1.625 # um
    wavelengths = [1.6,3.5,10]

    imagelist_1um = [
        f'../r3dresults/st28gm06n052_staranddust_1/186/image_i000_phi000_{shortwavelength}um.out',
        f'../r3dresults/st28gm06n052_staranddust_1/190/image_i000_phi000_{shortwavelength}um.out',
        f'../r3dresults/st28gm06n052_staranddust_1/198/image_i000_phi000_{shortwavelength}um.out'
    ]

    imagelist_3um = [
        '../r3dresults/st28gm06n052_staranddust_1/186/image_i000_phi000_3.5um.out',
        '../r3dresults/st28gm06n052_staranddust_1/190/image_i000_phi000_3.5um.out',
        '../r3dresults/st28gm06n052_staranddust_1/198/image_i000_phi000_3.5um.out',
    ]

    imagelist_10um = [
        '../r3dresults/st28gm06n052_staranddust_1/186/image_i000_phi000_10um.out',
        '../r3dresults/st28gm06n052_staranddust_1/190/image_i000_phi000_10um.out',
        '../r3dresults/st28gm06n052_staranddust_2/198/image_i000_phi000_10um.out',
    ]
    # 198 at 10um has a spike at default seed
    imagelists = [
        imagelist_1um,imagelist_3um,imagelist_10um
    ]
    # Maximum vertical scales
    vmaxima = [
        4,4,1
    ]
    # Contour levels
    starlevels = [
        [0.2,1.2,2.2,3.2,4.2,5.2,6.2],
        [0.2,1.2,2.2,3.2,4.2,5.2],
        [0.2,0.4,0.6,0.8,1.0,1.2,1.4]
    ]
    backlevels = [
        [0.0000001,0.00001,0.001,0.1],
        [0.00001,0.0001,0.001,0.01,0.1],
        [0.00001,0.0001,0.001,0.01,0.1]
    ]
    # Initialise fig-ax for normal images
    figVLTI, axVLTI = plt.subplots(
        3,3, 
        figsize = (11,10),
        num='VLTI at 90pc'
    )
    #         figsize = (12,7.3)  ifall jag har 3*2.bilder

    # Create a new figure with contour-plots of convolved images
    figcontour, axcontour = plt.subplots(
        3,3, 
        figsize = (11,10),
        num='Contours at 90pc'
    )



    distance = distanceVLTI

    # Loop through wavelengths/images
    for nwave,imagelist in enumerate(imagelists):

        # Full-width at half-maximum: lambda/(2*baseline)
        fwhmVLTI = wavelengths[nwave]*1e-6 / (2*baselineVLTI) * radian
        sigmaVLTI = fwhmVLTI/2.355

        # For reference, print FWHM:
        print(f'VLTI {wavelengths[nwave]}um FWHM: {fwhmVLTI} mas')

        # Loop through images/phases
        for nn,image in enumerate(imagelist):

            # Extract path and imagename from image
            imagestrings = re.split('/', image)

            modelname = imagestrings[2]
            imagefilename = imagestrings[4]
            phase = imagestrings[3]
            path = f'{imagestrings[0]}/{imagestrings[1]}/{modelname}/{phase}'

            # Load image
            image2d,image2dlog,flux,axisplot = a3d.load_images(
                path=path,
                image=imagefilename,
                distance=distance
            )
            # Change to Jy/mas2
            image2d = image2d * 1e-6

            # Extract props and compute distance-dependant scales
            Npix = np.shape(image2d)[0] # Number of pixels along one side
            auperpixel = 2*axisplot[0]/Npix  # number of au per pixel
            masperpixel = auperpixel/distance * 1000  # number of mas per pixel
            size_milliasec = masperpixel * Npix # size of image in mas

            # Image axis-limits
            axisplotmilliasec = [0.5*size_milliasec,-0.5*size_milliasec,-0.5*size_milliasec,0.5*size_milliasec]

            # Change sigma to order in number of pixels
            sigmaVLTI_pixels = sigmaVLTI / masperpixel

            # Convolve with Gaussian filter
            image2d = scipy.ndimage.gaussian_filter(
                image2d,
                sigma = sigmaVLTI_pixels
            )

            # Plot image
            imVLTI = axVLTI[nwave][nn].imshow(
                image2d, 
                origin='lower', extent=axisplotmilliasec, 
                cmap=plt.get_cmap('hot'),
                vmin=0, vmax=vmaxima[nwave]
            )
            axVLTI[nwave][nn].tick_params(axis='both', which='major', labelsize=15)

            # FWHM-cirle to show beam
            axVLTI[nwave][nn].add_patch(
                plt.Circle(
                    (0.75*axisplotmilliasec[0],0.75*axisplotmilliasec[1]), 
                    radius=0.5*fwhmVLTI, color='cyan', fill=True)
            )
            # Star-radius-circle to show size of stellar disc (two circles to make it more visible)
            axVLTI[nwave][nn].add_patch(
                plt.Circle(
                    (0,0), 
                    radius=1.65/distance*1000, 
                    color='lime', fill=False, linestyle=':', linewidth=2
                )
            )
            axVLTI[nwave][nn].add_patch(
                plt.Circle(
                    (0,0), 
                    radius=1.65/distance*1000, 
                    color='b', fill=False, linestyle=':', linewidth=1
                )
            )

            # Plot contours in separate figure, first with background
            axcontour[nwave][nn].contour(
                image2d,
                origin='lower', extent=axisplotmilliasec,
                colors='k',
                levels=backlevels[nwave],
                linewidths=1
            )
            # Then contours for the star
            axcontour[nwave][nn].contour(
                image2d,
                origin='lower', extent=axisplotmilliasec,
                colors='r',
                levels=starlevels[nwave],
                linewidths=1
            )
            # Flip x-axis
            axcontour[nwave][nn].invert_xaxis()

            # FWHM-cirle to show beam
            axcontour[nwave][nn].add_patch(
                plt.Circle(
                    (0.75*axisplotmilliasec[0],0.75*axisplotmilliasec[1]), 
                    radius=0.5*fwhmVLTI, color='green', fill=True, zorder=10)
            )
            # Star-radius-circle to show size of stellar disc
            axcontour[nwave][nn].add_patch(
                plt.Circle(
                    (0,0), 
                    radius=1.65/distance*1000, 
                    color='b', fill=False, linestyle='--', linewidth=2, 
                    zorder=10
                )
            )
            axcontour[nwave][nn].tick_params(axis='both', which='major', labelsize=15)


            # Set titles for first row for both figs
            if nwave == 0:
                axVLTI[nwave][nn].set_title(rf'$t_{nn+1} = $\,{phasetimes[nn]}', fontsize=15)
                axcontour[nwave][nn].set_title(rf'$t_{nn+1} = $\,{phasetimes[nn]}', fontsize=15)
            # Set xlabels for final row for both figs
            if nwave == 2:
                axVLTI[nwave][nn].set_xlabel('Offset (mas)',fontsize=18)
                axcontour[nwave][nn].set_xlabel('Offset (mas)',fontsize=18)


        axVLTI[nwave][0].set_ylabel('Offset (mas)',fontsize=18)
        axcontour[nwave][0].set_ylabel('Offset (mas)',fontsize=18)

        # Add text in figure with wavelength
        axcontour[nwave][0].text(
            x = -10, y = 60,
            s = rf'\noindent $\lambda = ${wavelengths[nwave]}\,$\mu$m',
            fontsize = 13
        )

        # Set colour bar settings and label
        divider = make_axes_locatable(axVLTI[nwave][-1])
        cax = divider.append_axes(
            'right', 
            size='5%', 
            pad=0.05
        )
        cb0 = plt.colorbar(
            imVLTI, 
            cax=cax, 
            orientation = 'vertical'
        )
        cb0.set_label(
            label = rf'$F_\nu$(Jy/mas$^2$) at {wavelengths[nwave]}\,$\mu$m',fontsize= 15
        )
        cb0.ax.tick_params(labelsize=15)

    # Add patches for the areas I take flux densities from
    # 186:
    # x: -5.5 till -4     -> mas  -27.5 -> -20
    # y: -4.5 till -3     -> mas  -22.5 -> -15
    axcontour[-1][0].add_patch(
        plt.Rectangle((-27.5,-22.5), 7.5, 7.5, color='orange', fill=False, zorder=10)
    )
    # 190:
    # x: 1 till 4       -> mas  5 -> 20  (15)
    # y: -2.5 till 0.5  -> mas  -12.5 -> 2.5 (15)
    axcontour[-1][1].add_patch(
        plt.Rectangle((5,-12.5), 15, 15, color='orange', fill=False, zorder=10)
    )
    # 198:
    # x: 1 - 6      -> mas  5 -> 30 (25)
    # y: -6.5 - -5  -> mas  -32.5 -> -25 (7.5)
    axcontour[-1][2].add_patch(
        plt.Rectangle((5,-32.5), 25, 7.5, color='orange', fill=False, zorder=10)
    )


    # Final settings for figures, save and show if you want to
    figVLTI.tight_layout()
    figcontour.tight_layout()

    #Save figure
    figVLTI.savefig(f'figs/images_VLTI_{distanceVLTI}pc.pdf', dpi=300, facecolor="white")
    figcontour.savefig(f'figs/contours_VLTI_{distanceVLTI}pc.pdf', dpi=300, facecolor="white")

    #figVLTI.show()
    #figcontour.show()


# -------------------------------------------------------------------------------
# Compute various observables
#

if compute_luminosities == 'y':
    # Compute luminosities of all SEDs without dust

    # Load st28gm06n052_nodust
    #   186, 190, 198
    # Load st28gm06n056
    #   140, 141, 142
    #
    paths = [
        '../r3dresults/st28gm06n052_nodust/',
        '../r3dresults/st28gm06n056/',
    ]
    phaseses = [
        ['186','190','198'],
        ['140','141','142']
    ]
    spectra = [
        '/spectrum_i000_phi000.out',
        '/spectrum_i090_phi000.out',
        '/spectrum_i090_phi090.out',
        '/spectrum_i090_phi270.out',
        '/spectrum_i180_phi000.out',
        '/spectrum_i270_phi000.out'
    ]

    for npath, path in enumerate(paths):
        print(path)
        phases = phaseses[npath]
        for phase in phases:

            lumsum = 0
            
            for spectrum in spectra:

                wavelength,sed = a3d.load_spectrum(
                    path = path+phase+spectrum
                )

                # Print all sed-luminosities
                luminosity = a3d.compute_luminosity(
                    wavelengths=wavelength,
                    spectrum=sed
                )
                luminosity = luminosity/Lsol
                print(f'{phase}{spectrum}: {luminosity}')

                lumsum += luminosity
            
            print(f'    Average: {lumsum/len(spectra)}')


if compute_tenmicronfluxdensities == 'y':

    # Computes total stellar flux density in various directions at only 10um
    # and totalt dust flux density in various directions (and wavelengths?)
    # lists these and Fdust/Fstar-ratios
    #
    # In debris discs we used Ldust/Lstar-ratios also...

    # Load stellar flux denisty:
    # r3dresults/st28gm06n052_nodust

    # Load dust flux density:
    # r3dresults/st28gm06n052_nostar

    models = [
        'st28gm06n052_staranddust_nospikes',
        'st28gm06n052_nodust',
        'st28gm06n052_nostar'
    ]
    contents = [
        'Total flux'
        'Stellar flux',
        'Dust flux'
    ]
    phases = [
        '186',
        '190',
        '198'
    ]
    spectra = [
        'spectrum_i000_phi000.out',
        'spectrum_i090_phi000.out',
        'spectrum_i090_phi090.out',
        'spectrum_i090_phi270.out',
        'spectrum_i180_phi000.out',
        'spectrum_i270_phi000.out'
    ]
    angles = [
        '0-0',
        '90-0',
        '90-90',
        '90-270',
        '180-0',
        '270-0'
    ]

    # Loop through models and print lsits for tables
    print('Total flux                   Stellar Flux                   Dust flux                   Contrast (Fd/F*)                Ratio (Ld/L*)')
    for phase in phases:
        print(f'  {phase}                            {phase}                     {phase}')

        for nangle,spectrum in enumerate(spectra):

            # Reset numbers
            totalflux = 0
            starflux = 0
            dustflux = 0

            # Load total fluxes, at 200pc
            # Load stellar fluxes, at 200pc, comparing to images I show
            wavelengths,totalfluxes = a3d.load_spectrum(
                path = f'../r3dresults/st28gm06n052_staranddust_nospikes/{phase}/{spectrum}',
                distance = 200
            )
            tenmicron = np.where(np.array(wavelengths) > 10)[0][0]
            totalflux = 0.5 * (
                totalfluxes[tenmicron-1] + totalfluxes[tenmicron]
            )



            # Load stellar fluxes, at 200pc, comparing to images I show
            wavelengths,starfluxes = a3d.load_spectrum(
                path = f'../r3dresults/st28gm06n052_nodust/{phase}/{spectrum}',
                distance = 200
            )
            tenmicron = np.where(np.array(wavelengths) > 10)[0][0]
            starflux = 0.5 * (
                starfluxes[tenmicron-1] + starfluxes[tenmicron]
            )
            # Also load stellar luminosity
            starlum = a3d.compute_sed_luminosity(
                path = f'../r3dresults/st28gm06n052_nodust/{phase}/{spectrum}',
                distance = 200
            )



            # Load dust fluxes
            wavelengths,dustfluxes = a3d.load_spectrum(
                path = f'../r3dresults/st28gm06n052_nostar/{phase}/{spectrum}',
                distance = 200
            )
            tenmicron = np.where(np.array(wavelengths) > 10)[0][0]
            dustflux = 0.5 * (
                dustfluxes[tenmicron-1] + dustfluxes[tenmicron]
            )
            # Also load dust luminosity
            dustlum = a3d.compute_sed_luminosity(
                path = f'../r3dresults/st28gm06n052_nostar/{phase}/{spectrum}',
                distance = 200
            )

            # Print everything
            print(f'    {angles[nangle]}:{totalflux}           {starflux}           {dustflux}           {dustflux/starflux}           {dustlum/starlum}')



if measuredustcloudflux == 'y':
    # Measures dust flux densities in patches that I define
    # Also plots the patches for reference

    fig, ax = a3d.plot_imagesubplots(
        imagelist = [
            '../r3dresults/st28gm06n052_staranddust_1/186/image_i000_phi000_10um.out',
            '../r3dresults/st28gm06n052_staranddust_1/190/image_i000_phi000_10um.out',
            '../r3dresults/st28gm06n052_staranddust_2/198/image_i000_phi000_10um.out',
        ],
        distance = 1,
        scale = 'lin'
    )

    # 186:
    # x: -5.5 till -4
    # y: -4.5 till -3
    ax[0].add_patch(
        plt.Rectangle((-5.5,-4.5), 1.5, 1.5, color='cyan', fill=False)
    )

    # 190:
    # x: 1 till 4
    # y: -2.5 till 0.5
    ax[1].add_patch(
        plt.Rectangle((1,-2.5), 3, 3, color='cyan', fill=False)
    )

    # 198:
    # x: 1 - 6
    # y: -6.5 - -5
    ax[2].add_patch(
        plt.Rectangle((1,-6.5), 5, 1.5, color='cyan', fill=False)
    )

    # Extract fluxes of these areas
    phases = [186,190,198]

    xranges = [
        [-5.5,-4],
        [1,4],
        [1,6]
    ]
    yranges = [
        [-4.5,-3],
        [-2.5,0.5],
        [-6.5,-5]
    ]

    # Initilize a second figure to check that I actually measure the correct areas
    fig2, ax2 = plt.subplots(1,3)

    # Loop over phases
    for nphase,phase in enumerate(phases):

        # Load images themselves
        image2d,image2dlog,totalflux,axisplot = a3d.load_images(
            path = f'../r3dresults/st28gm06n052_staranddust_1/{phase}/',
            image = 'image_i000_phi000_10um.out',
            distance = 200
        )
        # Change to Jy/mas2
        image2d = image2d * 1e-6

        # Create axis-array
        axscale = np.linspace(axisplot[0],axisplot[1],image2d.shape[0])

        # rectangle-ranges:
        xrange = np.where((axscale <= xranges[nphase][1]) & (axscale >= xranges[nphase][0]))[0]
        yrange = np.where((-axscale >= yranges[nphase][0]) & (-axscale <= yranges[nphase][1]))[0]

        # Fluxes within
        fluxes = image2d[yrange[0]:yrange[-1],xrange[0]:xrange[-1]]

        print(f'{phase}, x: {xranges[nphase][0]} to {xranges[nphase][1]}')
        print(f'     y: {yranges[nphase][0]} to {yranges[nphase][1]}')
        print(f'     max: {fluxes.max()}     mean: {fluxes.mean()} Jy mas-2')

        ax2[nphase].imshow(fluxes, origin='lower', cmap=plt.get_cmap('hot'))


    fig.show()
    fig2.show()



if plot_resolutiondistance == 'y':
    # Plots image estimate of within which distance may
    # one resolive surface features?
    # papers to reference for this?


    # TODO
    # REDO THIS PLOT
    # remove this plot? Not used any longer

    # Full-width at half-maximum in mas:  1.22*lambda/D 
    # or for VLTI, lambda/2baseline
    fwhmVLTI_10um = 1e-5 / (2*baselineVLTI) * radian
    fwhmJWST_10um =  1.22 * 1e-5 / diameterJWST * radian
    fwhmVLTI_1um = 1.625*1e-6 / (2*baselineVLTI) * radian
    fwhmJWST_1um =  1.22 * 1e-6 / diameterJWST * radian

    print('FWHM')
    print(f'{fwhmVLTI_10um} mas, VLTI at 10um')
    print(f'{fwhmVLTI_1um} mas, VLTI at 1.625um')
    print(f'{fwhmJWST_10um} mas, JWST at 10 um')
    print(f'{fwhmJWST_1um} mas, JWST at 1 um')
    print('')


    # Beam area in mas2
    areaVLTI_10um = np.pi * (0.5*fwhmVLTI_10um)**2
    areaVLTI_1um = np.pi * (0.5*fwhmVLTI_1um)**2
    areaJWST_10um = np.pi * (0.5*fwhmJWST_10um)**2
    areaJWST_1um = np.pi * (0.5*fwhmJWST_1um)**2

    print('Beam area (fwhm)')
    print(f'{areaVLTI_10um} mas2, VLTI at 10um')
    print(f'{areaVLTI_1um} mas2, VLTI at 1.625um')
    print(f'{areaJWST_10um} mas2, JWST at 10um')
    print(f'{areaJWST_1um} mas2, JWST at 1um')
    print('')

    # Stellar disc area as function of distance
    distance = np.logspace(0,np.log10(2e3),1000)  # distance array in pc
    detect_limit = np.linspace(4,4,1000)  # detect limit is area-ratio of 4
    Rstar = 1.65 / distance * 1000  # Star radius in mas (au/pc)
    Astar = np.pi*Rstar**2  # Star surface area in mas2

    # Ratio between stellar disc and beam surfaces, or R^2
    area_ratio_VLTI_10um = Astar / areaVLTI_10um
    area_ratio_VLTI_1um = Astar / areaVLTI_1um
    area_ratio_JWST_10um = Astar / areaJWST_10um
    area_ratio_JWST_1um = Astar / areaJWST_1um

    # Indeces of "resolution limits"
    indeces_detected_VLTI_10um = np.where(area_ratio_VLTI_10um >= 4)[0]
    indeces_detected_VLTI_1um = np.where(area_ratio_VLTI_1um >= 4)[0]
    indeces_detected_JWST_10um = np.where(area_ratio_JWST_10um >= 4)[0]
    indeces_detected_JWST_1um = np.where(area_ratio_JWST_1um >= 4)[0]

    print(f'VLTI at  1.625um: <{distance[indeces_detected_VLTI_1um[-1]]} pc')
    print(f'VLTI at 10    um: <{distance[indeces_detected_VLTI_10um[-1]]} pc')
    print(f'JWST at  1    um: <{distance[indeces_detected_JWST_1um[-1]]} pc')
    print(f'JWST at 10    um: <{distance[indeces_detected_JWST_10um[-1]]} pc')


    # Plot figure
    fig,ax = plt.figure(), plt.axes()

    # First fillbetween gave bugs in pdf, mitigated by plotting from index 300
    ax.plot(distance,area_ratio_VLTI_1um,'b',linewidth = 2)
    ax.plot(distance,area_ratio_VLTI_1um,'k',linewidth = 1)
    ax.fill_between(
        distance[indeces_detected_VLTI_1um[300:]],
        detect_limit[indeces_detected_VLTI_1um[300:]],
        area_ratio_VLTI_1um[indeces_detected_VLTI_1um[300:]],
        color='b',
        alpha=1
    )
    ax.plot(distance,area_ratio_VLTI_10um,'c',linewidth = 2)
    ax.plot(distance,area_ratio_VLTI_10um,'k--',linewidth = 1)
    ax.fill_between(
        distance[indeces_detected_VLTI_10um[300:]],
        detect_limit[indeces_detected_VLTI_10um[300:]],
        area_ratio_VLTI_10um[indeces_detected_VLTI_10um[300:]],
        color='c',
        alpha=1
    )

    ax.plot(distance,area_ratio_JWST_1um,'r',linewidth = 2)
    ax.plot(distance,area_ratio_JWST_1um,'k-.',linewidth = 1)
    ax.fill_between(
        distance[indeces_detected_JWST_1um],
        detect_limit[indeces_detected_JWST_1um],
        area_ratio_JWST_1um[indeces_detected_JWST_1um],
        color='r',
        alpha=1
    )
    ax.plot(distance,area_ratio_JWST_10um,'orange',linewidth = 2)
    ax.plot(distance,area_ratio_JWST_10um,'k:',linewidth = 1)
    ax.fill_between(
        distance[indeces_detected_JWST_10um],
        detect_limit[indeces_detected_JWST_10um],
        area_ratio_JWST_10um[indeces_detected_JWST_10um],
        color='orange',
        alpha=1
    )

    ax.text(
        x=1.5,y=6.5,
        s=r'JWST: 10 \&\ 1$\mu$m',
        backgroundcolor='white',
        fontsize=15
    )
    ax.text(
        x=40,y=7.5,
        s=r'VLTI: 10 \&\ 1.625$\mu$m',
        backgroundcolor='white',
        fontsize=15
    )

    ax.plot(distance,detect_limit,'k--')
    ax.set_ylim(1,10)
    ax.set_xlim(1,2000)
    ax.set_xscale('log')
    ax.set_xlabel('Distance (pc)', fontsize=18)
    ax.set_ylabel(r'$R_\star^2 / R_{\rm Beam}^2$', fontsize=18)
    ax.tick_params(axis='both', which='major', labelsize=15)

    fig.tight_layout()
    fig.show()

    fig.savefig(f'figs/plot_resolution.pdf', dpi=300, facecolor="white")


# Check "radius" of source at various images
if check_smoothedimage_radius == 'y':

    # Standard images
    #
    #    '../r3dresults/st28gm06n052_staranddust_1/186/image_i000_phi000_1.625um.out'
    #    '../r3dresults/st28gm06n052_staranddust_1/190/image_i000_phi000_1.625um.out'
    #    '../r3dresults/st28gm06n052_staranddust_1/198/image_i000_phi000_1.625um.out'
    #
    #    '../r3dresults/st28gm06n052_staranddust_1/186/image_i000_phi000_3.5um.out'
    #    '../r3dresults/st28gm06n052_staranddust_1/190/image_i000_phi000_3.5um.out'
    #    '../r3dresults/st28gm06n052_staranddust_1/198/image_i000_phi000_3.5um.out'
    #
    #    '../r3dresults/st28gm06n052_staranddust_1/186/image_i000_phi000_10um.out'
    #    '../r3dresults/st28gm06n052_staranddust_1/190/image_i000_phi000_10um.out'
    #    '../r3dresults/st28gm06n052_staranddust_2/198/image_i000_phi000_10um.out'

    distance = 200 # parsec
    relativelimit = 0.1
    phase = 198
    wavelength = 10 # um

    # Load correct seed also
    if phase == 198 and wavelength == 10:
        seed = 2
    else:
        seed = 1
    image = f'../r3dresults/st28gm06n052_staranddust_{seed}/{phase}/image_i000_phi000_{wavelength}um.out'

    # Angular resolution of "beam" is then approx:
    fwhmVLTI = wavelength*1e-6 / (2*baselineVLTI) * radian
    sigmaVLTI = fwhmVLTI/2.355


    # Load star's radius here (for comparison plot)
    # Mstar: gram
    # Rstar: cm
    # Lstar: W
    #Mstar,Rstar,Lstar = a5d.load_star_information(
    #    savpath='../co5bold_data/dst28gm06n052/st28gm06n052_186.sav',
    #    printoutput='n'
    #)
    #Rstar /= AUcm
    Rstar = 1.65*AUcm


    # Extract path and imagename from image
    imagestrings = re.split('/', image)
    modelname = imagestrings[2]
    imagefilename = imagestrings[4]
    phase = imagestrings[3]
    path = f'{imagestrings[0]}/{imagestrings[1]}/{modelname}/{phase}'

    # Load image
    image2d,image2dlog,flux,axisplot = a3d.load_images(
        path=path,
        image=imagefilename,
        distance=distance
    )
    # Change to Jy/mas2
    image2d = image2d * 1e-6


    # Extract props and compute distance-dependant scales
    Npix = np.shape(image2d)[0] # Number of pixels along one side
    auperpixel = 2*axisplot[0]/Npix  # number of au per pixel
    masperpixel = auperpixel/distance * 1000  # number of mas per pixel
    size_au = auperpixel * Npix # Size of image in au


    # Change sigma to order in number of pixels
    sigmaVLTI_pixels = sigmaVLTI / masperpixel

    # Convolve with Gaussian filter
    image2d = scipy.ndimage.gaussian_filter(
        image2d,
        sigma = sigmaVLTI_pixels
    )


    # Number of circular annulii are a quarter to make 
    # sure there are enough pixels in each annulus
    # And range is half the size of the image
    radial_range = np.linspace(auperpixel,size_au*0.5,Npix*0.25)
    radial_fluxes = np.zeros(int(Npix*0.25))
    radial_npixels = np.zeros(int(Npix*0.25))


    # List of pixel numbers with 0 in centrum, radially and along each axis
    rpix = list(range(0,int(Npix*0.5),2))
    xpix = list(range(int(-Npix*0.5),int(Npix*0.5)))
    ypix = list(range(int(-Npix*0.5),int(Npix*0.5)))


    # Loop through annulii, add fluxes of each
    print('Looping through image')
    for nn,nr in enumerate(rpix):
        for nx in xpix:
            for ny in ypix:
                if np.sqrt(nx**2 + ny**2) > nr and \
                    np.sqrt(nx**2 + ny**2) <= nr+2:
                    
                    imagex = int(nx+Npix*0.5)
                    imagey = int(ny+Npix*0.5)

                    # Add fluxes of each annulus and number of pixels
                    radial_fluxes[nn] += image2d[imagey,imagex]
                    radial_npixels[nn] += 1
    print('    done')

    # Average flux per annulus
    radial_fluxes /= radial_npixels

    # At what radius is half maximum? Or 25%? Or 10%?
    fluxlimit = relativelimit*image2d.max()
    if fluxlimit < radial_fluxes.max():
        radius_fluxlimit = radial_range[np.where(radial_fluxes >= fluxlimit)[0].max()]
    else:
        radius_fluxlimit = 0
    print('')
    print(f'Max flux of smoothed image: {image2d.max()} Jy')
    print(f'Max-flux ratio (image2d/annulus): {image2d.max()/radial_fluxes.max()}')
    print(f'{modelname}_{phase}, {imagefilename}')
    print(f'Obs-radius: {radius_fluxlimit} AU')

    # Plot to check
    plt.plot(radial_range,radial_fluxes)
    plt.plot([Rstar,Rstar],[0,radial_fluxes.max()])
    plt.plot([radius_fluxlimit,radius_fluxlimit],[0,fluxlimit])
    plt.xlim(0,4)
    plt.show()

    plt.plot(radial_range,radial_npixels,'.')


if plot_smoothedimage_radius == 'y':

    # To plot source-radius vs wavelength
    wavelengths = [1.6,3.5,10.0]

    t1_50perc = [0,2.0587,2.4116]
    t2_50perc = [np.nan,np.nan,1.7057]
    t3_50perc = [0,1.4704,2.0071]

    t1_25perc = [1.9410,2.6469,2.9999]
    t2_25perc = [1.2351,1.9410,2.4116]  # This
    t3_25perc = [0,2.1763,2.6469]

    t1_10perc = [2.5293,2.9999,3.3528]
    t2_10perc = [1.9410,2.4116,3.5881]  # This
    t3_10perc = [2.0587,2.6469,3.1175]

    # Manually change these to whatever you want to plot
    plot_these_lists = [t2_50perc,t2_25perc,t2_10perc]
    legendlist = [
        r'50\,\%',
        r'25\,\%',
        r'10\,\%'
    ]

    # Initialise plot (add more styles if needed)
    fig, ax = plt.figure(figsize=(6, 4)), plt.axes()
    linestyles = ['-','--','-.',':']
    markerstyles = ['v','o','*','^']

    for nn,sourceradius in enumerate(plot_these_lists):

        ax.plot(
            wavelengths,sourceradius,
            linestyle=linestyles[nn],marker=markerstyles[nn],
            markersize = 8,
            label = legendlist[nn]
        )

    ax.plot([1,11],[1.65,1.65],'k:')
    ax.legend(fontsize=13)

    ax.set_xlabel(r'Wavelength ($\mu$m)',fontsize=18)
    ax.set_ylabel(r'Source radius (au)',fontsize=18,)
    ax.set_xlim([1,11])
    ax.tick_params(axis='both', which='major', labelsize=15)

    fig.tight_layout()

    fig.savefig('figs/source_radius.pdf', dpi=300, facecolor="white")
    fig.show()

