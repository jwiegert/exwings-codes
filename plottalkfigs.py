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

plot_coboldseds = 'n'
plot_cobolddarwin_compare = 'y'


# ----------------------------------------------------------------
#
# Plot Figures with subplots of Seds

if plot_coboldseds == 'y':

    #paths = ['../r3dresults/st28gm06n052_staranddust_nospikes/']
    #paths = ['../r3dresults/st28gm06n052_darwinsource/']
    #paths = ['../r3dresults/st28gm06n052_pointsource/']
    paths = ['../r3dresults/st28gm06n052_pointtemperature/']
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
    fig,ax = plt.subplots(len(paths),len(phases), figsize = (12, 4))


    # Extract SEDs and fill figure-axes-objects
    for nphase,phase in enumerate(phases):
        for npath,path in enumerate(paths):
            for nangle,spectrum in enumerate(spectra):

                wavelength,sed = a3d.load_spectrum(
                    path = path+phase+spectrum
                )

                ax[nphase].plot(
                    wavelength,sed,
                    label = legendlist[nangle],
                )

            ax[nphase].set(
                xscale='log',yscale='log'
            )
            ax[nphase].title.set_text(phase)

        ax[nphase].set_xlim(5e-1,6e1)
        ax[nphase].set_ylim(1e6,1.3e8)
        ax[nphase].tick_params(axis='both', which='major', labelsize=15)


    # Axis, labels and so on settings    
    ax[0].set_ylabel(r'Flux density (Jy at 1 pc)', fontsize=18)
    ax[1].set_xlabel(r'Wavelength ($\mu$m)',fontsize=18)
    ax[2].legend(title=r'$i$, $\phi$')
    
    fig.tight_layout()
    #fig.savefig('figs/pointtemperature_seds.pdf', dpi=300, facecolor="white")
    fig.show()


if plot_cobolddarwin_compare == 'y':

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


    fig,ax = plt.subplots(1,len(phases), figsize = (12, 4))


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

        # Plot line at sigma-limit of chi2-plots
        ax[nphase].plot(
            [wavelength[0],wavelength[-1]],[1,1],'r:'
        )
        # Set to log-scale
        ax[nphase].set(
            xscale='log',yscale='linear'
        )
        ax[nphase].set_ylim(0,30)
        ax[nphase].tick_params(axis='both', which='major', labelsize=15)


    # Compute chi2-numbers of each SED and add to plots
    for nn in range(len(sedcoboldnumbers)):
        if sedcoboldnumbers[nn] == sedsdarwinnumbers[nn]:


            chisq_array, chiaq_reduced = a3d.compute_chisquare(
                simulation = np.array(sedsdarwin[nn]),
                observation = np.array(sedscobold[nn]),
                obssigma = 0.1*np.array(sedscobold[nn])
            )

            # Plot line (angle order should be the same as above)
            nphase = int(sedcoboldnumbers[nn][0])
            nangle = int(sedcoboldnumbers[nn][-1])
            ax[nphase].plot(
                wavelength,np.sqrt(chisq_array),
                label = legendlist[nangle]
            )




    # TODO
    # Limit axis to something useful

    # ax.set dessa endast på yttre plots
    #                 ylabel=f'Flux density (Jy at 1 pc)',
    #                xlabel=r'Wavelength ($\mu$m)',
    ax[0].set(
        ylabel=r'$|\chi(\lambda)|$',
    )
    ax[1].set(
        xlabel=r'Wavelength ($\mu$m)',
    )
    ax[2].legend(title=r'$i$, $\phi$')

    ax[0].tick_params(axis='both', which='major', labelsize=15)
    ax[1].tick_params(axis='both', which='major', labelsize=15)
    ax[2].tick_params(axis='both', which='major', labelsize=15)


    fig.tight_layout()
    fig.show()
    










