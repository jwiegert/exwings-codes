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
plot_grainsizeradius = 'n'

plot_absscat = 'n'

plot_temperatureradial = 'n'
plot_temperaturecompare = 'n'

# Make these as just one humongous figure?
# that would be 3*4 SEDs in one figure. 
# My main results?
plot_manyseds = 'y'
#plot_seds_cobold = 'n'
#plot_seds_darwin = 'n'
#plot_seds_points = 'n'
#plot_seds_pointstemper = 'n'

# plot co5bold-pictures
# plot one darwin and one point-source-picture, 190 I think




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

    phase = 198

    # Load cobold-T and create subplot
    T_c5d,Tstd_c5d,Tminmax_c5d,radius_c5d = a3d.plot_temperaturebins_radius(
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


    # Load r3d-T
    T_r3d,Tstd_r3d,Tminmax_r3d,radius_r3d = a3d.plot_temperaturebins_radius(
        temperature_path=f'../r3dresults/st28gm06n052_pointtemperature/{phase}/dust_temperature.dat',
        grid_path='../r3dresults/st28gm06n052_pointtemperature/grid_distances.csv',
        amr_path='../r3dresults/st28gm06n052_pointtemperature/amr_grid.inp',
        numb_specie = 1,
        ax=ax[1]
    )
    ax[1].set_ylabel(r'$T_{\rm RADMC-3D}$ (K)',fontsize=18)
    ax[1].set_xlabel(r'')
    ax[1].set(xlim=(0,26))
    ax[1].tick_params(axis='both', which='major', labelsize=15)

    # A plot with Tc5d / Tr3d
    ax[2].plot(radius_r3d,T_c5d/T_r3d,'b')

    ax[2].fill_between(
        radius_r3d,
        (T_c5d-Tstd_c5d)/(T_r3d-Tstd_r3d),
        (T_c5d+Tstd_c5d)/(T_r3d+Tstd_r3d),
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

    # Show all plots
    plt.tight_layout()
    fig.show()

    #Save figure
    fig.savefig(f'figs/temperatures_{phase}.pdf', dpi=300, facecolor="white")




# ----------------------------------------------------------------
#
# Plot Figure subplots of ALL Seds

if plot_manyseds == 'y':

    paths = [
        '../r3dresults/st28gm06n052_staranddust_nospikes/',
        '../r3dresults/st28gm06n052_darwinsource/',
    ]
#        '../r3dresults/st28gm06n052_pointsource/',
#        '../r3dresults/st28gm06n052_pointtemperature/',
    
    
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
        r'0, 0',
        r'90, 0',
        r'90, 90',
        r'90, 270',
        r'180, 0',
        r'270, 0',
    ]


    fig,ax = plt.subplots(len(paths)+1,len(phases))


    # Spara SEDerna och kör över som jämförelser till en tredje rad    
    # Kör en slags jämförelse mot 10% av c5d-model istället för ena kurvan
    # ie, chi2 = 1 -> kurvorna är 0.1* varandra
    # chi2 > 1 -> mer än 0.1 från andrandra
    #
    # kanske köra att "sigma" är 100% av c5d-model istället?
    # det innebär att chi=1 -> 2*model = c5dmodel?
    
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


                ax[npath][nphase].plot(wavelength,sed)

            ax[npath][nphase].set(
                xscale='log',yscale='log'
            )


    # Compute chi2-numbers of each SED and add to plots
    for nn in range(len(sedcoboldnumbers)):
        if sedcoboldnumbers[nn] == sedsdarwinnumbers[nn]:
            #print(f'{sedcoboldnumbers[nn]} | {sedsdarwinnumbers[nn]}')
            chisq_array, chiaq_reduced = a3d.compute_chisquare(
                simulation = np.array(sedsdarwin[nn]),
                observation = np.array(sedscobold[nn]),
                obssigma = 0.1*np.array(sedscobold[nn])
            )
            
            # Plot line (angle order should be the same as above)
            nphase = int(sedcoboldnumbers[nn][0])
            ax[2][nphase].plot(wavelength,chisq_array)
            ax[2][nphase].set(
                xscale='log',yscale='log'
            )

    # Add line at sigma-limit of chi2-plots
    # dvs chi2 = 1 -> D = 1.1C eller 0.9C, om sigma = 0.1C
    ax[2]



    # TODO
    # Limit axis to something useful

    # ax.set dessa endast på yttre plots
    #                 ylabel=f'Flux density (Jy at 1 pc)',
    #                xlabel=r'Wavelength ($\mu$m)',
    ax[0][0].set(
        ylabel=r'Flux density (Jy at 1 pc)',
    )
    ax[1][0].set(
        ylabel=r'Flux density (Jy at 1 pc)',
    )
    ax[2][0].set(
        ylabel=r'$\chi^2(\lambda)$',
    )

    
    
    ax[-1][1].set(
        xlabel=r'Wavelength ($\mu$m)',
    )
    ax[1][-1].legend(legendlist)




    fig.tight_layout()
    fig.show()
    



