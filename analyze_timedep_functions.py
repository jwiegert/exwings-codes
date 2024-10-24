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
# extract_sourcesize()
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


# Extract radii of "main source" at included snapshots and at various
# angles
def extract_sourcesize(
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
        wavelength:str = '01',
        relativelimit:float = 0.12,
        save_datafile:str = 'y',
    ):
    """
    Extracts a corresponding source radius (in au) based on the total area of all
    pixels that contain a flux density above a certain fraction of the maximum
    flux density of the image.

    ARGUMENTS
      path: a string to the path to your model data
      phases: list with strings to the subfolder (snapshotfolders) of your model
      angles: list of angles to investigate, see default setting for syntax
              'i000_phi000'
      wavelength: a string with wavelength (um) of image to look at, syntax: '01' for 1um
      relativelimit: a float with fraction of image max flux to base area on
                     Default value (0.12) is based on tests with only star of models
                     st28gm06n052, st28gm06n074, st28gm06n075
      save_datafile: string, 'y' or 'n' for yes or no on saving a dat-file with all output

    RETURNS
      stellar_radius_average: len(phases)-long array with angle-averaged radii in au
      stellar_radii: len(phases)*len(angles)-shaped array with all radii for all phases
                     and angles.
      (optional): dat file with all these data named '{path}source_radius_{wavelength}um.dat'
    """

    # Extract some number
    Nphases = len(phases)
    Nangles = len(angles)

    # Declare Radius(time and angle)-array
    stellar_radii = np.zeros((Nphases,Nangles))
    stellar_radius_average = np.zeros(Nphases)


    # Loop through angles and snapshots
    for nphase,phase in enumerate(phases):
        for nangle,angle in enumerate(angles):
            #
            # Load image (image2d is in Jy/asec2)
            image2d,image2dlog,flux,axisplot = a3d.load_images(
                path=f'{path}/{phase}/',
                image=f'image_{angle}_{wavelength}um.out',
                distance=1
            )
            # Extract props and compute distance-dependant scales
            Npix = np.shape(image2d)[0]              # Number of pixels along one side
            Npixhalf = int(Npix*0.5)                 # Half size
            auperpixel = np.abs(2*axisplot[0]/Npix)  # Number of au per pixel
            fluxlimit = relativelimit*image2d.max()  # Flux density limit of source area

            # Loop through image and count pixels >fluxlimit
            Npixels = 0
            for nx in range(Npix):
                for ny in range(Npix):
                    if image2d[ny,nx] >= fluxlimit:
                        Npixels += 1

            # Compute corresponding radius at each angle and snapshot
            stellar_radii[nphase,nangle] = np.sqrt(
                auperpixel*auperpixel*Npixels/np.pi
            )
        # Extract angle-averaged radius for each snapshot
        stellar_radius_average[nphase] = np.mean(stellar_radii[nphase,:])

    # Save data if requested
    if save_datafile == 'y':
        with open(f'{path}source_radius_{wavelength}um.dat', 'w') as fradius:
            #
            # Write header (and all included angles) 
            fradius.write(f'# Average radius of main source at {wavelength} in au\n')
            fradius.write('# Snapshot  R_average')
            for angle in angles:
                fradius.write(f'   R_{angle}')
            fradius.write('\n#\n')
            #
            # Write radii for each time
            for nphase,phase in enumerate(phases):
                fradius.write(f'  {phase}       {stellar_radius_average[nphase]:.5f}')
                # And included angle
                for nangle,angle in enumerate(angles):
                    fradius.write(f'     {stellar_radii[nphase,nangle]:.5f}    ')
                # Add new line at each snapshotline
                fradius.write('\n')
        print(f'extract_sourcesize: Done\n  Output written to {path}source_radius_{wavelength}um.dat')
    else:
        # Else just return data
        return stellar_radius_average, stellar_radii

