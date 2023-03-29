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
plot_absscat = 'n'
plot_temperatureradial = 'n'
plot_temperaturecompare = 'y'

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

if plot_temperaturecompare == 'y':
    # Load cobold-T
    fig1,ax1, T_c5d,Tstd_c5d,Tminmax_c5d,radius_c5d = a3d.plot_temperaturebins_radius(
        temperature_path='../r3dresults/st28gm06n052_staranddust/186/dust_temperature.dat',
        grid_path='../r3dresults/st28gm06n052_staranddust/grid_distances.csv',
        amr_path='../r3dresults/st28gm06n052_staranddust/amr_grid.inp',
        numb_specie = 1
    )

    # Load r3d-T
    fig2,ax2, T_r3d,Tstd_r3d,Tminmax_r3d,radius_r3d = a3d.plot_temperaturebins_radius(
        temperature_path='../r3dresults/st28gm06n052_pointtemperature/186/dust_temperature.dat',
        grid_path='../r3dresults/st28gm06n052_pointtemperature/grid_distances.csv',
        amr_path='../r3dresults/st28gm06n052_pointtemperature/amr_grid.inp',
        numb_specie = 1
    )


    # Compute and plot chi2
    Tsigma = Tstd_c5d + Tstd_r3d
    Tminmax = Tminmax_c5d + Tminmax_r3d

    # With std
    fig3,ax3,c2,c2red = a3d.plot_chisquare(
        T_c5d,T_r3d,Tsigma,radius_c5d
    )
    ax3.set(xscale='lin')

    # with minmax
    fig4,ax4,c2,c2red = a3d.plot_chisquare(
        T_c5d,T_r3d,Tminmax,radius_c5d
    )
    ax4.set(xscale='lin')


    # A plot with Tc5d / Tr3d
    fig5,ax5 = plt.figure(figsize=(6, 4)), plt.axes()
    ax5.plot(radius_r3d,T_c5d/T_r3d,'b')

    ax5.fill_between(
        radius_r3d,
        (T_c5d-Tminmax_c5d)/(T_r3d-Tminmax_r3d),
        (T_c5d+Tminmax_c5d)/(T_r3d+Tminmax_r3d),
        color='b',
        alpha=0.2
    )

    ax5.fill_between(
        radius_r3d,
        (T_c5d-Tstd_c5d)/(T_r3d-Tstd_r3d),
        (T_c5d+Tstd_c5d)/(T_r3d+Tstd_r3d),
        color='b',
        alpha=0.4
    )

    ax5.plot([1.65,1.65],[0,10],'r:')
    ax5.set(
        ylim=(0,5)
    )

    # Show all plots
    fig1.show()
    fig2.show()
    fig3.show()
    fig4.show()
    fig5.show()




