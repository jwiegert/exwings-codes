# Functions for analyzing time dependent R3D-output
# import as a3t
# ------------------------------------------------------------ #
# Various useful packages
import matplotlib.pyplot as plt
import numpy as np
import scipy


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
# extract_surfacetemp()
# extract_imageblobs()
# load_imageblob_files()
# compute_period()
# extract_events()
# 
#
# TODO?
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

    for phase in phases:
        for angle in angles:

            # Loop over all SEDs
            wavelength, fluxes = a3d.load_spectrum(
                path=f'{path}{phase:03d}/spectrum_{angle}.out'
            )
            
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
        path=f'{path}{phases[0]:03d}/spectrum_{angles[0]}.out'
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
                path=f'{path}{phase:03d}/spectrum_{angle}.out'
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
    # Automatically add / to end of path if it's missing
    if path[-1] != '/':
        path += '/'

    # Extract some number
    Nphases = len(phases)
    Nangles = len(angles)

    # Declare Radius(time and angle)-array
    stellar_radii = np.zeros((Nphases,Nangles))
    stellar_radius_average = np.zeros(Nphases)

    # Average flux densities of nodust models are:
    #
    #    average_nodust
    #    1.00791589806 um: 42398912.69025754 Jy
    #    2.010669926771 um: 97403587.81062156 Jy
    #    10.07270043879 um: 20601793.91693989 Jy
    #
    #
    # Define flux limits depending on wavelength-choice
    if int(wavelength) == 1:
        fluxlimit = 0.035*42398912.69025754
    elif int(wavelength) == 2:
        fluxlimit = 0.045*97403587.81062156
    elif int(wavelength) == 10:
        fluxlimit = 0.047*20601793.91693989
    else:
        raise ValueError('  ERROR: wrong wavelength, only 1, 2, and 10 are available.')

    # Loop through angles and snapshots
    for nphase,phase in enumerate(phases):
        for nangle,angle in enumerate(angles):
            #
            # Load image (image2d is in Jy/asec2)
            image2d,image2dlog,flux,axisplot = a3d.load_images(
                path=f'{path}{phase:03d}/',
                image=f'image_{angle}_{wavelength}um.out',
                distance=1
            )
            # Extract props and compute distance-dependant scales
            Npix = np.shape(image2d)[0]              # Number of pixels along one side
            auperpixel = np.abs(2*axisplot[0]/Npix)  # Number of au per pixel

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
            fradius.write(f'# Using a flux limit of {fluxlimit} Jy.\n')
            fradius.write('# Snapshot  R_average')
            for angle in angles:
                fradius.write(f'   R_{angle}')
            fradius.write('\n#\n')
            #
            # Write radii for each time
            for nphase,phase in enumerate(phases):
                fradius.write(f'  {phase:03d}       {stellar_radius_average[nphase]:.5f}')
                # And included angle
                for nangle,angle in enumerate(angles):
                    fradius.write(f'     {stellar_radii[nphase,nangle]:.5f}    ')
                # Add new line at each snapshotline
                fradius.write('\n')
        print(f'extract_sourcesize: Done\n  Output written to {path}source_radius_{wavelength}um.dat')
    else:
        # Else just return data
        return stellar_radius_average, stellar_radii


# Extract temperatures of a shell around the radius listed in the files created by
# extract_sourcesize().
def extract_surfacetemp(
        path_model:str='../r3dresults/',
        path_radii:str='../r3dresults/source_radius_01um.dat',
    ):
    """
    TODO write some info, tar alla snapshots i radius-filen o kör på
    Tar bara temperaturen från stjärnan! viktigt!
    
    ARGUMENTS
      path_radii: string with path to folder containing radius as func of time
      path_model: string to folder containing all snapshot data
    
    RETURNS
      surface_temperature.dat: file with average temperatures of shell at Radius
                               for each included snapshot/phasenumber
    """
    print('Running atf.extract_surfacetemp()')

    # Automatically add / to end of path if it's missing
    if path_model[-1] != '/':
        path_model += '/'

    # Extract modelname from model path
    modelname = path_model.split('/')[-2]

    # Load radii and phase numbers
    source_radii = np.loadtxt(path_radii)
    source_radius_au = source_radii[:,1]
    source_radius_cm = source_radii[:,1] * AUcm
    phases = source_radii[:,0]

    # Load grid cell sizes and positions
    nbasecells,gridside,nrefinements,nleafs,cellsizes,gridref_in,gridref_out = a3d.load_grid_information(
        gridinfo_path=f'{path_model}/grid_info.txt'
    )
    # put all refinements limits and courner of box in one array
    gridref = np.concatenate((gridref_in,gridref_out[::-1],[np.sqrt(3)*0.5*gridside]))

    # Load radial grid distances to all cells
    griddistances = a3d.load_griddistances(
        gridpath=f'{path_model}/grid_distances.csv',
        amrpath=f'{path_model}/amr_grid.inp'
    )
    griddistances = griddistances[:,0]

    # Form an array to save temperatures in
    surface_temperatures = np.zeros(len(phases))

    # Loop through all snapshots
    for nphase,phase in enumerate(phases):
        # Remove float-status on phase-number to work with folder names
        phase = int(phase)

        # Load cellsize of cells at this phase's radius
        refindex = np.argwhere(gridref >= source_radius_au[nphase])[0][0]-1

        if refindex == -1 or refindex == 7:
            # Base cells are in grid centrum and outside the refinements
            cellsize = cellsizes[0]
        elif refindex < 4:
            # Inner refinements (within to surface of star)
            cellsize = cellsizes[refindex+1]
        elif refindex > 3 and refindex < 7:
            # Outer refinements (outside the star)
            refindex = 2-refindex
            cellsize = cellsizes[refindex]

        # Load indeces of grid cells around Rstar plus/minus cellsize, in cm
        cellsize *= AUcm
        radius_range_cm = [
            source_radius_cm[nphase]-cellsize,
            source_radius_cm[nphase]+cellsize,
        ]
        shellindeces = np.argwhere(
            (griddistances >= radius_range_cm[0]) & (griddistances <= radius_range_cm[1])
        )
        if len(shellindeces) == 0:
            raise ValueError('  ERROR: no shellindeces found!')

        # Extract all temperatures at this radius (with these indeces) and save average
        # temperature for each snapshot
        Ncells,Nspecies,temperatures = a3d.load_temperature(
            path=f'{path_model}/{phase:03d}/dust_temperature_onestar_smoothed.dat'
        )
        surface_temperatures[nphase] = np.mean(temperatures[shellindeces])

    # Write dat-file with temperatures
    with open(f'{path_model}/surface_temperature.dat', 'w') as ft:
        # Write header
        ft.write(f'# Approximate stellar surface temperature as based on model {modelname}\n')
        #
        # Write data
        for phase,Tsurface in zip(phases,surface_temperatures):
            ft.write(f'    {int(phase):03d}    {Tsurface:.3f}\n')

    print('Extract approximate surface temperature: Done')


# Function to find number of and sizes of blobs in images around the star
def extract_imageblobs(
        imagepath:str='../image_i090_phi090_10um.out',
        Rstar:float=1.65,
        max_flux_contrast:float=0.01,
        Flimit_exposure:float=240e3,
        fract_stararea:float=0.5,
    ):
    """
    Identify and return statistics on dust blobs in images with a star in the
    middle of the image. Blobs are outside Rin and inside Rout from middle of
    image (in AU). Blobs pixels must be higher than a low cut-off, for example
    pixel flux density must be stronger than the contrast to the star that is
    allowed by the telescope. In the case of VLTI/MATISSE, more 1 per cent 
    (private conversation, Claudia Paladini, ESO, Chile).

    ARGUMENTS
      imagepath:str = path to imagefile
      Rstar:float = radius of star in AU
      max_flux_contrast:float = A number <1 that defines the allowed contrast
                                to the central star's flux density for dust
                                clouds to be observable
      Flimit_exposure:float = Flimit based on exposure time. 240kJy at 1pc is
                              based on 32min with VLTI/MATISSE at 10um reaches 
                              SNR~3.3 for 6Jy, 1AU-source at 200pc.
      fract_stararea:float = A number <1 that defines smallest allowed size
                             in multiples of the stellar disc's surface rea
                             in the image plane.

    RETURNS
      Nlargeblobs : the number of blobs with image surface area larger than
                    fract_stararea x star_area and with pixels brighter than
                    max_flux_contrast x stellar flux density
      np.max(blobareas) : area in AU^2 of the largest dust blob in the image
      stellarflux_info : String with info on how stellar flux is estimated
    """
    # Set annulus limits
    Rin = 2*Rstar
    Rout = 6*Rstar

    # Define fractional stellar area-limit in AU^2
    StarArea = np.pi * Rstar**2
    FracStarArea = fract_stararea * StarArea


    # Extract path to folder and to image
    path = '/'.join(imagepath.split('/')[:-1])+'/'
    imagename = imagepath.split('/')[-1]

    # Load image, pixel units are in Jy/asec2
    image2d,image2dlog,totalflux,axisplot = a3d.load_images(
        path=path,
        image=imagename
    )
    # Extract stellar flux density of image with no dust for Flimit-number
    # and change its unit to Jy/AU2 by division with constant stellar disc
    # area.
    if 'nospikes' in path:
        nodustpath = path.replace('nospikes', 'nodust')
        image2dstar,image2dlog,stellarflux,axisplotstar = a3d.load_images(
            path=nodustpath,
            image=imagename
        )
        stellarflux /= StarArea
        stellarflux_info = f'Nodust data exist; uses stellar flux and/or exposure limit {Flimit_exposure} Jy/AU2'
    else:
        # Or if no star data exist, use the max flux of the image but in
        # units of AU2 , ie approximate stellar flux density with max-pixel
        # averaged over star's disc surface area. It's normally approximately
        # the same as the star's flux density divided by stellar area.
        stellarflux = np.max(image2d)/StarArea
        stellarflux_info = ' !NO! nodust-data does not exist; uses max pixel flux density averaged over stellar surface and/or exposure limit {Flimit_exposure} Jy/AU2'

    # Extract image props
    # Pixel-resolution of image
    Npix = np.shape(image2d)[0] 
    halfNpix = Npix/2
    # Size of each pixel in AU
    pixsize_au = axisplot[-1]/halfNpix
    # And create new image array 
    image2dmod = np.zeros(np.shape(image2d))
    # Set pixel flux limit
    Flimit = max_flux_contrast * stellarflux

    # Loop through image remove pixels outside radial ranges or below flux limit
    for xpix in range(Npix):
        for ypix in range(Npix):

            # Current pixel-distance to centre in AU
            Rpix = np.sqrt(
                (xpix - halfNpix)**2 + (ypix - halfNpix)**2
            )*pixsize_au

            # Set pixels outside R-limits to zero
            if Rpix > Rout or Rpix < Rin:
                image2dmod[ypix,xpix] = 0

            # Set pixels below exposure limit to zero
            elif image2d[ypix,xpix] < Flimit_exposure:
                image2dmod[ypix,xpix] = 0

            # Set pixels below contrast-to-star-limit to zero
            elif image2d[ypix,xpix] < Flimit:
                image2dmod[ypix,xpix] = 0

            # And the rest to 1
            else:
                image2dmod[ypix,xpix] = 1

    # Label blobs in binary image

    # Array defining neighbouring pixels for each blob identification
    structure = np.ones((3, 3), dtype=int)
    # Extract connected components
    image2dlabeled, ncomponents = scipy.ndimage.label(image2dmod, structure)

    if ncomponents > 0:
        # Extract size of each component
        # Array to fill with blob-areas in AU^2
        blobareas = np.zeros(ncomponents)
        # And counter of large blobs
        Nlargeblobs = 0

        for nblob in range(1,ncomponents+1):

            # Number of pixels per blob time pixel au^2 size
            blobareas[nblob-1] = np.where(image2dlabeled.ravel()==nblob)[0].size * pixsize_au**2

            # Number of blobs larger than Fractional stellar area
            if blobareas[nblob-1] >= FracStarArea:
                Nlargeblobs += 1

        # Return Number of large blobs and array with area of all blobs
        return Nlargeblobs, blobareas, stellarflux_info

    else:
        # If there are no components, return just zeros and warning
        return 0,np.zeros(1),'WARNING: No components found'



# function that loads data from ascii files with information on 
def load_imageblob_files(
        filepath:str='../r3dresults/st28gm06n052_timedep_nospikes/',
        max_flux_contrast:float=0.01,
        fract_stararea:float=0.1,
        load_blobareas:str='y'
    ):
    """
    Loads and returns arrays with pre-computed numbers on dusty blobs
    in Radmc3d-simulated images.

    ARGUMENTS
      filepath:str = path to folder that contains files. Filenames are
                     with this syntax:
                     imageblobs_maxarea_Flim{max_flux_contrast}.dat
                     imageblobs_numb_Flim{max_flux_contrast}_Farea{fract_stararea}.dat
      max_flux_contrast:float = Flux contrast-limit, contrast between blob-pixels
                                and stellar flux density.
      fract_stararea:float = Smallest allowed size of cloud in image-plane. Set
                             as fractions of stellar disc area.
    RETURNS
      angles: 1D array of angles (strings) in data
      nsnaps: 1D array of snap shot numbers included in data (ints)
      nblobs: 2D array with number of dust clouds per snapshot and angle
      blob_areas: 2D array with area (AU^2) of largest dust cloud per snapshot and angle
    """
    # Automatically add / to end of path if it's missing
    if filepath[-1] != '/':
        filepath += '/'

    # Declare various lists
    nsnaps = []
    nblobs = []
    blob_areas = []

    # Load all number of blobs
    with open(f'{filepath}imageblobs_numb_Flim{max_flux_contrast}_Farea{fract_stararea}.dat', 'r') as fnblobs:
        for nblobs_line in fnblobs.readlines():
            # Extract angles
            if '#  Phase  i' in nblobs_line:
                angles = nblobs_line.split()[2:]
            # Extract all number of blobs
            if nblobs_line[0] != '#':
                nblobs.append(nblobs_line.split()[1:])
                nsnaps.append(int(nblobs_line.split()[0]))
    # Change from strings in lists to arrays with numbers
    nblobs = np.array(nblobs).astype('int')
    nsnaps = np.array(nsnaps).astype('int')


    # Load all blog areas (if set to yes)
    if load_blobareas == 'y' or load_blobareas == 'Y' or load_blobareas == 'yes':
        # Changes so that each line in list is an array with floats
        # Each float is the area of each cloud at that specific nsnap and nangle.
        # Number of snapshots and angles are given above.
        with open(f'{filepath}imageblobs_blobareas_Flim{max_flux_contrast}.dat', 'r') as farea:
            for area_line in farea.readlines():
                if area_line[0] != '#':
                    area_line = area_line.split('[')[1][:-2]
                    area_line = area_line.split(', ')
                    blob_areas.append(np.array(area_line).astype('float'))
    else:
        # If not loading all areas, return list with zeros
        blob_areas = [0]


    # return all angles, snapshotnumbers, numb of blobs and largest blob areas
    return angles,nsnaps,nblobs,blob_areas




# Uses FFT to compute period of signal
def compute_period(
        signal:list=[0,1,0,-1,0,1,0,-1,0,1,0,-1,0],
        timeaxis:list=[0,1,2,3,4,5,6,7,8,9,10,11,12],
        plot_spec:str='n'
    ):
    """
    Input list or array with some signal and time axis and get a main period
    and a plot of the power spectrum of the signal.
    """
    from scipy.signal import find_peaks

    # Normalise signal amplitude to zero
    signal_zeroed = signal - np.mean(signal)

    # Take fourier transform and abs to get power spectrum of signal
    signal_fft = np.abs(np.fft.rfft(signal_zeroed, norm="ortho"))

    # Get corresponding frequencies for the power spectrum
    Ntimesteps = len(timeaxis)
    freqs = np.fft.fftfreq(Ntimesteps)[:len(signal_fft)]

    # And two main periods
    delta_timestep = (timeaxis[-1] - timeaxis[0])/Ntimesteps
    period_axis = 1/freqs * delta_timestep
    peakcoords = np.argpartition(signal_fft, -2)[-4:]
    periods = period_axis[peakcoords[::-1]]

    # Save for emergencies
    #main_period = period_axis[np.argmax(signal_fft)]

    # Plot to check
    if plot_spec == 'y':
        plt.figure(num='Fourier power spectrum')
        plt.plot(period_axis,signal_fft)

    return periods


# Extract events statistics
def extract_events(
        eventdata:list=[],
        timerange:list=[],
        relative_limit:float=0.1,
        show_data:str='y'
    ):
    """
    Extracts spike-statistics from input data array/list. Spikes, or events, are defined
    as when the data reaches above the input relative limit (times max-data) above the
    median data.

    ARGUMENTS
      eventdata: list or array with spikey data
      timerange: list of array with corresponding time (x) axis
      relative_limit: events are when data > median(data) + relative_limit * max(data)
      show_data: str, 'y' or 'n', 'y' will plot the data and the limits to check if they
                 are as expected.
    
    RETURNS
      Prints statistics in stdout. Plots figure if wanted.
    """
    # Reset timearray to start at zero
    timerange = timerange - timerange[0]

    # Extract useful statistics
    data_median = np.mean(eventdata)
    data_max = np.max(eventdata)

    # That gives event limits
    event_limit = data_median + relative_limit*(data_max-data_median)

    # Extract event statistics, number of and length of events
    event_sum = 0    # Length of all events
    event_count = 0  # Number of events

    prev_data = 0    
    for data in eventdata:
        if data > event_limit:
            event_sum += 1
            if prev_data <= event_limit:
                event_count += 1
        prev_data = data
    
    # Translate event_sum to length in time units
    event_sum *= timerange[1]

    # Print output
    # TODO write periodicitet
    # = totala tiden delat på antal events -> tid mellan events i medel
    print(f'  Total included time: {timerange[-1]} time units')
    print(f'  Relative limit: {relative_limit}')
    print(f'  Number of events: {event_count}')
    print(f'  Average pediod between events: {timerange[-1]/event_count} time units')
    print(f'  Total length of events: {event_sum} time units')
    print(f'  Average length of events: {event_sum/event_count} time units')
    print(f'  Portion of events: {event_sum/timerange[-1]}')

    if show_data == 'y':
        # Plot the data and events in case its wanted
        fig,ax = plt.figure(), plt.axes()
        ax.plot(
            timerange,eventdata
        )
        ax.plot(
            [timerange[0],timerange[-1]],
            [data_median,data_median],
            'k--'
        )
        ax.plot(
            [timerange[0],timerange[-1]],
            [event_limit,event_limit],
            'k:'
        )
        fig.show()

