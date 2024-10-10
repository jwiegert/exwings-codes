# Functions for analyzing time dependent R3D-output
# import as a3t
# ------------------------------------------------------------ #
# Various useful packages
import analyze_r3d_functions as a3d
import matplotlib.pyplot as plt
import numpy as np

# My own functions
import analyze_r3d_functions as a3d

# Basic definitions
AUcm = 1.49598e13 # cm
Msol = 1.989e33 # g
Rsol = 6.955e10 # cm


# ------------------------------------------------------------ #
# List of functions
#
# plot_allseds()
# extract_averageseds()
#
# TODO
# dustmass(time)-plotter
# F_lambda(time)-plotter
# F_lambda-variations-to-average-extracter
#
# ------------------------------------------------------------ #
#
# Functions below
#
# ------------------------------------------------------------ #
#
# Plot functions

# Plot all SEDs and overlay with average, std, and min-max SEDs
def plot_allseds(
        path:str='',
        phases:list=[],
        angles = [
            'i000_phi000',
            'i090_phi000',
            'i090_phi090',
            'i090_phi270',
            'i180_phi000',
            'i270_phi000',
        ]
    ):
    """
    Plots all listed seds in one plot within path/phases folders, with
    specified angles.

    ARGUMENTS
      path: string with path to folder with all phase-subfolders
      phases: list of phase numbers (sub folder names)


    RETURNS
      fig, ax : figure and axis matplotlib-objects
    """

    # Set settings for the plot
    fig, ax = plt.figure(figsize=(6, 4)), plt.axes()
    ax.set(
        xscale='log',
        yscale='log',
        xlim=[0.5,40],
        ylim=[1e5,2e8]
    )
    ax.set_ylabel(r'$F_\nu$, Jy at 1 pc', fontsize=18)
    ax.set_xlabel(r'Wavelength ($\mu$m)',fontsize=18)
    ax.tick_params(axis='both', which='major', labelsize=15)

    # Extract length of wavelengthgrid
    wavelength, fluxes = a3d.load_spectrum(
        path=f'{path}{phases[0]}/spectrum_{angles[0]}.out'
    )
    Nwave = len(wavelength)

    # And total number of SEDs
    Nseds = len(phases) * len(angles)

    # Array with all SEDs
    allfluxes = np.zeros((Nwave,Nseds))

    # Pre-set SED number counter
    sedcounter = -1


    for phase in phases:
        for angle in angles:
            sedcounter += 1

            # Loop over all SEDs
            wavelength, fluxes = a3d.load_spectrum(
                path=f'{path}{phase}/spectrum_{angle}.out'
            )
            # Extract and save all flux densities at each wavelength
            allfluxes[:,sedcounter] = fluxes
            
            # Plot all SEDs as grey
            ax.plot(wavelength,fluxes,'lightgrey')

    # Extract average, std and minmax SEDs
    flux_average,flux_std,flux_min,flux_max = extract_averageseds(
        path=path,
        phases=phases,
        angles=angles,
        save_datafile='n',
    )

    # Plot average, stdrange and minmax
    ax.plot(wavelength,flux_average,'black')
    ax.plot(wavelength,flux_average-flux_std,'k--')
    ax.plot(wavelength,flux_average+flux_std,'k--')
    ax.plot(wavelength,flux_min,'k:')
    ax.plot(wavelength,flux_max,'k:')

    fig.tight_layout()

    return fig,ax



#---------------------------------------------------------------#
# Extract various data-functions

# Save average, std, and minmax-SEDs
def extract_averageseds(
        path:str = '',
        phases:list = [],
        angles:list = [
            'i000_phi000',
            'i090_phi000',
            'i090_phi090',
            'i090_phi270',
            'i180_phi000',
            'i270_phi000',
        ],
        save_datafile:str = 'y',
    ):
    """
    Extracts and returns arrays with average, std and minmax SEDs of
    input SED files.

    ARGUMENTS
      path: string with path to folder with all phase-subfolders
      phases: list of phase numbers (sub folder names)
      angles: list of angles in style 'i000_phi000'
      save_datafile: string, 'y' to save data to dat file, 'n' to return as arrays

    RETURNS
      Either dat file or arrays with
        average,std,min,max flux densities (in this order)
    """

    # Extract length of wavelengthgrid
    wavelength, fluxes = a3d.load_spectrum(
        path=f'{path}{phases[0]}/spectrum_{angles[0]}.out'
    )
    Nwave = len(wavelength)

    # And total number of SEDs
    Nseds = len(phases) * len(angles)

    # Array with all SEDs
    allfluxes = np.zeros((Nwave,Nseds))

    # Number the SEDs
    sedcounter = -1

    for phase in phases:
        for angle in angles:
            sedcounter += 1

            # Loop over all SEDs
            wavelength, fluxes = a3d.load_spectrum(
                path=f'{path}{phase}/spectrum_{angle}.out'
            )
            # Extract and save all flux densities at each wavelength
            allfluxes[:,sedcounter] = fluxes

    if save_datafile == 'y':
        # Create output file and save SEDs in it
        with open(f'{path}average_sed.dat', 'w') as fsed:
            # Declare/reset variables
            flux_average = 0.0
            flux_std = 0.0
            flux_min = 0.0
            flux_max = 0.0

            # Write header
            fsed.write(f'# Average SEDs of model {path.split("/")[-2]} snapshots {phases[0]} - {phases[-1]}, all in Jy at 1 pc\n')
            fsed.write('# Wavelength       AverageF              STD                  MaxF                 MinF\n')

            # Loop through wavelengths and save average, std, and minmax values at each wavelength
            for nwave in range(Nwave):
                flux_average = np.mean(allfluxes[nwave,:])
                flux_std = np.std(allfluxes[nwave,:])
                flux_min = np.min(allfluxes[nwave,:])
                flux_max = np.max(allfluxes[nwave,:])

                # Save fluxes
                fsed.write(f'{wavelength[nwave]}    {flux_average}    {flux_std}    {flux_max}    {flux_min}\n')
        return f'Extract_averageseds: written file {path}average_sed.dat'

    else:
        # Just return arrays with data
        # Declare arrays
        flux_average = np.zeros(Nwave)
        flux_std = np.zeros(Nwave)
        flux_min = np.zeros(Nwave)
        flux_max = np.zeros(Nwave)

        # Loop through wavelengths and save average, std, and minmax values at each wavelength
        for nwave in range(Nwave):
            flux_average[nwave] = np.mean(allfluxes[nwave,:])
            flux_std[nwave] = np.std(allfluxes[nwave,:])
            flux_min[nwave] = np.min(allfluxes[nwave,:])
            flux_max[nwave] = np.max(allfluxes[nwave,:])

        return flux_average,flux_std,flux_min,flux_max
