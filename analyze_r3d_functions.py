# Various functions useful for analyzing in and out-data of RADMC3D
# ------------------------------------------------------------ #
# Useful packages
import os
import csv
import re
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm, xlabel, xscale, yscale
from mpl_toolkits.axes_grid1 import make_axes_locatable

# My own libraries
import analyze_co5bold_functions as a5d

# Useful numbers
c = 2.998e8 # speed of light in m/s
pc = 30.857e15 # 1 parsec in m
AUcm = 1.49598e13 # AU in cm
kboltz = 1.3806503e-23 # Boltzmann constant in Si
hplanck = 6.626068e-34 # Planck constant in Si


# ------------------------------------------------------------ #
# List of functions
#
#
# Functions that load various r3d input data
# ------------------------------------------
#
# load_grid_information()
#    TODO fill info
#
# load_grid_properties(
#    amrpath:str='../amr_grid.inp'
# )
#
# load_griddistances(
#    gridpath:str='../grid_distances.csv',
#    amrpath:str='../amr_grid.inp',
# )
#
# load_cellsizes(
#    sizepath:str='../grid_cellsizes.csv',
#    amrpath:str='../amr_grid.inp',
# )
#
# load_wavelengthgrid(
#    path:str='../wavelength_micron.inp'
# )
#
# load_dustdensity(
#    path:str = '../dust_density.inp',
#    numb_specie:int = 1
# )
#
# load_grainsizes(
#    grainsize_path:str='../grain_sizes_186.dat'
# )
#
# load_onekappa(
#    specie_name:str='',
#    specie_number:int=0,
#    path:str='../'
# )
#
#
# Load output data from R3D
# -------------------------
#
# load_temperature(
#    path:str='../dust_temperature.dat',
#    numb_specie:int=1
# )
#
# load_spectrum(
#    path:str='../r3dsims/spectrum.out',
#    distance:float=1
# )
#
# load_images(
#    path:str='../',
#    image:str='image.out',
#    distance:float=1
# )
#
# extract_dustmasses(
#    log_files = [],
#    Nspecies = 1,
#    first_species = 1
# )
#
# load_photocentre_file(
#   file_path:str='../r3dresults/photocentre.dat'
# )
#
# Functions to set Plot-settings
# ------------------------------
#
# TODO implement or remove these?
#
# set_figurefonts():
#
# set_ticklabels(
#    xlabel:str,
#    ylabel:str,
#    xscale:str,
#    yscale:str
# )
#
#
# Plot various details of input and output data of R3D
# ----------------------------------------------------
#
# plot_grid(
#    gridpath:str = '../grid_distances.csv',
#    sizepath:str = '../grid_cellsizes.csv',
#    amrpath:str = '../amr_grid.inp',
#    nbins:int = 10
# )
#
# plot_onedensity_radius(
#    density_path:str='../dust_density.inp',
#    grid_path:str='../grid_distances.csv',
#    amr_path:str='../amr_grid.inp',
#    numb_specie:int=1
# )
#
# plot_alldensity_radius(
#    path:str='../'
# )
#
# These two plots all cell's temperatures
# plot_onetemperature_radius(
#    temperature_path:str='../dust_temperature.dat',
#    grid_path:str='../grid_distances.csv',
#    amr_path:str='../amr_grid.inp',
#    numb_specie:int=1
# )
#
# plot_alltemperature_radius(
#    path:str='../'
# )
#
# This plots binned averages of the temperatures instead
# plot_temperaturebins_radius(
#    temperature_path:str='../dust_temperature.dat',
#    grid_path:str='../grid_distances.csv',
#    amr_path:str='../amr_grid.inp',
#    numb_specie:int=1
# )
#
#
# plot_grainsize_radius(
#    gridpath:str='../r3dresults/st28gm06n052/grid_distances.csv',
#    amrpath:str='../r3dresults/st28gm06n052/amr_grid.inp',
#    grainsizepath:str='../grain_sizes_186.dat'
# )
#
# plot_onekappa(
#    specie_name:str = '',
#    specie_number:int = 0,
#    path:str = '../'
# )
#
# plot_allkappa(
#    path:str='../'
# )
#
# plot_sed(
#    path:str='../spectrum.out',
#    distance:float=1
# )
#
# plot_sedsmany(
#    pathlist:list=['../spectrum.out'],
#    legendlist:list=['spectrum1'],
#    distance:float=1
# )
#
# plot_images(
#    path:str='../',
#    images:list=['image.out'],
#    distance:float=1
# )
#
# plot_imagesubplots(
#    imagelist:list=['../image.out'],
#    distance:float=1,
#    scale:str='lin'
# )
#
# plot_imagecrosssections(
#    path:str='../',
#    images:list=['image.out'],
#    distance:float=1
# )
#
# plot_opticalthick(
#    path:str = '../',
# )
#
#
# Compute different quantities
# ----------------------------
#
# compute_sed_luminosity(
#    path:str = '../r3dsims/spectrum.out',
#    distance:float = 1
# )
#
# compute_luminosity(
#    wavelengths:list,
#    spectrum:list,
#    distance:float=1
# )
#
# translate_filter(
#    filterfile
#    wavelengthr3d
#)
#
# remove_sedspikes()
#    path:str = '../r3dsims/spectrum.out',
#    distance:float = 1
# )
#
# remove_imagespikes(
#    folders:list=['../'],
#    imagefilename:str='image.out'
# )
#
# compare_seds(
#    sedpaths:list=[
#       '../spectrum0.out',
#       '../spectrum1.out'
#    ]
# )
#
#
# compare_images(
#    TODO
# )
#
# compute_blackbody(
#   peak_flux
#   peak_freq
#   wavelengths
# )
#
# compute_photocentre(
#   TODO
#)
#
#
# ------------------------------------------------------------ #
# Functions that load various r3d input data

# Load grid information
def load_grid_information(
        gridinfo_path:str='../grid_info.txt'
    ):
    """
    Loads basic grid information.

    ARGUMENTS
      gridinfo_path: path to grid_info.txt-file as generated by grid
                     radmc3d-creation-functions.

    RETURNS
      nbasecells: int with number of base cells along one cube side
      gridside: float with length of one cube side in AU
      nrefinements: int with number of grid refinements
      nleafs: int with number of cells in whole grid
      cellsizes: nrefinements+1 long np.array with all unique cell sizes in AU
      gridref_in: nrefinements long np.array with inner radial refinements limints in AU
      gridref_out: nrefinements long np.array with outer radial refinements limints in AU    
    """

    # Check if file exists
    if os.path.exists(gridinfo_path) == True:

        # Extract basic numbers
        with open(gridinfo_path, 'r') as f:
            for nn,line in enumerate(f.readlines()):

                # Extract number of base cells along one side
                if nn == 3:
                    nbasecells = int(re.findall('\d+', line)[0])
                # Extract size of grid side
                if nn == 4:
                    gridside = float(re.findall('\d+.\d+', line)[0])
                # Extract number of refinements
                if nn == 5:
                    nrefinements = int(re.findall('\d', line)[0])
                # Extact number of cells (nleafs)
                if nn == 6:
                    nleafs = int(re.findall('\d+', line)[0])

        # Extract cell sizes and distances to refinements (in and out)
        cellsizes = np.zeros(nrefinements+1)
        gridref_in = np.zeros(nrefinements)
        gridref_out = np.zeros(nrefinements)

        with open(gridinfo_path, 'r') as f:
            for nn,line in enumerate(f.readlines()):
                
                # Extract cell sizes
                if nn == 2:
                    # Base cells
                    cellsizes[0] = float(re.findall('\d+\.\d+', line)[0])

                # Extract refinement limits
                # Line 8 is nref = 1
                if line[:34] == f'Radial distances to refinement {nn-7}: ':
                    numbers = re.findall('\d+\.\d+', line)
                    gridref_in[nn-8] = float(numbers[0])
                    gridref_out[nn-8] = float(numbers[1])

                # Child cell sizes
                # Line 13 is nref 1
                if line[:19] == f'Child cell size {nn-12}: ':
                    cellsizes[nn-12] = float(re.findall('\d+\.\d+', line)[0])

        return nbasecells,gridside,nrefinements,nleafs,cellsizes,gridref_in,gridref_out

    else:
        return f'ERROR: load_grid_information can not find {gridinfo_path}.'




# Load grid properties
def load_grid_properties(
        amrpath:str='../amr_grid.inp'
    ):
    """
    Loads basic proparties from amr_grid.inp

    ARGUMENTS
      amrpath: str, path to amr_grid.inp

    RETURNS
      (In this order)
      nxyz,nrefines,nleafs,nbranch,gridedge
      - Number of base cells
      - Number of grid refinements
      - Number of grid cells (nleafs)
      - Number of branches
      - Size of the whole grid cube
    """

    # Check if file exists
    if os.path.exists(amrpath) == True:

        with open(amrpath, 'r') as f:
            for nn,line in enumerate(f.readlines()):

                # Number of base cells along one side
                if nn == 7:
                    nxyz = int(line.split(' ')[0])

                # Number of refinements,
                # number of cells (nleafs),
                # number of branches
                if nn == 9:
                    nrefines = int(line.split(' ')[0])
                    nleafs = int(line.split(' ')[1])
                    nbranch = int(line.split(' ')[2])
                
                # Total size of the grid in cm
                if nn == 11:
                    gridedge = np.abs(float(line))
        gridedge *= 2
        
        # Return everything!
        return nxyz,nrefines,nleafs,nbranch,gridedge
    
    else:
        return f'ERROR: load_grid_properties can not find {amrpath}.'


# Load griddistances
def load_griddistances(
        gridpath:str='../grid_distances.csv',
        amrpath:str='../amr_grid.inp',
    ):
    """
    Loads array of distances to cells of the current grid.
    Distances are from centrum of grid (radial,x,y, and z) in centimeters.

    ARGUMENTS
      gridpath: path to grid_distances.csv
      amrpath: path to amr_grid.inp

    RETURNS
      griddistances: np.array of size nleafs*4. 
        Each column is distances in cm
        [:,0] : radial distances
        [:,1] 2] and 3] : x,y, and z distances
    """

    # Check if file exists
    if os.path.exists(gridpath) == True:

        # Extract necessary info from amr_grid
        nleafs = load_grid_properties(amrpath)[2]

        # Create griddistances array
        griddistances = np.zeros((nleafs,4))

        # Load distances to cells
        with open(gridpath, 'r') as f:
            griddistancesfile = csv.reader(f)

            # Set index o 0
            nn = 0
            for row in griddistancesfile:

                # Skip the header and add radial, x,y,z distances to array
                if row[0][-1].isdigit() == True:
                    
                    # Save radial, x, y, z distances (in this order)
                    griddistances[nn,:] = [float(row[0]),float(row[1]),float(row[2]),float(row[3])]
                    # Increase index
                    nn += 1

        # Return griddistances array    
        return griddistances
    
    else:
        raise TypeError(f'ERROR: load_griddistances can not find {gridpath}.')


# Load list of grid cell sizes
def load_cellsizes(
        sizepath:str='../grid_cellsizes.csv',
        amrpath:str='../amr_grid.inp',
    ):
    """
    Loads array of grid cell sizes in centimeters
    
    ARGUMENTS
      sizepath: path to grid_cellsizes.csv
      amrpath: path to amr_grid.inp

    RETURNS
      gridsizes: np.array of size nleafs*1 with cell sizes in same order as in dust_density
    """

    # Check if file exists
    if os.path.exists(sizepath) == True:

        # Extract necessary info from amr_grid
        nleafs = load_grid_properties(amrpath)[2]

        # Create griddistances array
        gridsizes = np.zeros(nleafs)

        # Load distances to cells
        with open(sizepath, 'r') as f:
            gridsizefile = csv.reader(f)

            # Set index o 0
            nn = 0
            for row in gridsizefile:

                # Skip the header and add radial, x,y,z distances to array
                if row[0][-1].isdigit() == True:
                    
                    # Save radial, x, y, z distances (in this order)
                    gridsizes[nn] = float(row[0])
                    # Increase index
                    nn += 1

        # Return griddistances array    
        return gridsizes
    
    else:
        return f'ERROR: load_cellsizes can not find {sizepath}.'


# Load R3d wavelength grid
def load_wavelengthgrid(
        path:str='../wavelength_micron.inp'
    ):
    """
    Ancillary function. Loads and extracts wavelengths from wavelength grid.
    
    Input
    -----
    path: str, path to wavelength_micron.inp

    Output
    ------
    wavelengths: list of wavelengths in micron
    nwave: number of wavelength grid points    
    """

    wavelengths = []

    with open(path,'r') as f:

        for nn,line in enumerate(f.readlines()):
            
            # Extract number of wavelengths
            if nn == 0:
                nwave = int(line)
            
            # Extract wavelengths
            else:
                wavelengths.append(float(line))
            
    return wavelengths,nwave


# Load one density
def load_dustdensity(
        path:str = '../dust_density.inp',
        numb_specie:int = 1
    ):
    """
    Load and extracts densities of one dust species of dust_density.inp

    INPUT
    -----
    path: string with path and filename of density file
    numb_specie: integer with number of the species you want to load (default=1)
    
    OUTPUT
    ------
    Ncells: Number of cells in grid (nleafs)
    Nspec: Number of species in dust_density file
    dust_densities: np.array containing densities of all cells for species numb_specie

    """

    if numb_specie == 0:
        print(f'ERROR number of species cant be 0')
    else:
        # Read header
        with open(path,'r') as f:
            for nn,line in enumerate(f.readlines()):

                # Number of cells
                if nn == 1:
                    Ncells = int(line)

                # Total number of species
                if nn == 2:
                    Nspecies = int(line)

        # Check that the chosen number of species exists
        if numb_specie > Nspecies:
            print('\nERROR\nThis dust species doesnt exist.\n')

        else:
            # Reduce species number by 1 (index starts at 0)
            numb_specie -= 1

            # Create density np.array
            dust_densities = np.zeros(Ncells)

            # Extract dust densities
            with open(path,'r') as f:
                for nn,line in enumerate(f.readlines()):
                    if nn > 2+numb_specie*Ncells and nn <= 2+(numb_specie+1)*Ncells:
                        dust_densities[nn-3-Ncells*numb_specie] = float(line)

            return Ncells,Nspecies,dust_densities


# Load grainsizes
def load_grainsizes(
        grainsize_path:str='../grain_sizes_186.dat'
    ):
    """
    TODO info
    Loads already extracted-to-R3D grain sizes
    ARGUMENTS
    RETURNS
      sizes : array with all grain sizes per grid cell
      Nleafs : number of cells
    """
    sizes = []
    if os.path.exists(grainsize_path) == True:
        with open(grainsize_path, 'r') as fsizes:
            for line in fsizes.readlines():
                if line[0] != '#':
                    sizes.append(float(line.strip('\n')))
    else:
        return f'ERROR: bin_grainsizes can not find {grainsize_path}.'
    sizes = np.array(sizes)
    Nleafs = np.size(sizes)

    return sizes,Nleafs


# Load absorptionscattering data
def load_onekappa(
        specie_name:str='',
        specie_number:int=0,
        path:str='../'
    ):
    """
    Loads absorption/scattering/scattering angles of one dust specie
    
    ARGUMENTS
      specie_name: a string
      OR
      specie_number: an integer
      path: path to folder containing dustkappa and dustopac-files

    RETURNS
      specie_name: a string
      wavelenths: a list of wavelengths in micro meters
      kappadata: a list with lsits containing 
          [0] absorption in cm^2/g
          [1] scattering in cm^2/g
          [2] total extinction (abs + scat)
          [3] mean scattering angles in <cos theta>
    """

    # Automatically add / to end of path if it's missing
    if path[-1] != '/':
        path += '/'

    # Extract name from dustopac-file if a number is given
    if specie_number != 0:
        with open(f'{path}dustopac.inp', 'r') as f:
            for nn,line in enumerate(f.readlines()):
                # Specie-names are on line numbers 1+specie_number*4
                if nn == 1+specie_number*4:
                    specie_name = line.strip()
        print(f'Extracting species number {specie_number}: {specie_name}')

    # Extract species by given name if name is given
    elif len(specie_name) != 0:
        print(f'Extracting species {specie_name}')
    
    # If nothing is given, print error
    else:
        print('ERROR: no species number nor name is given.')

    # Extract absorption, scattering and scattering angles
    
    # First we need to automatically adapt according to available filename
    if os.path.exists(f'{path}dustkappa_{specie_name}.inp') == True:
        dustkappa_filename = f'{path}dustkappa_{specie_name}.inp'

    elif os.path.exists(f'{path}dustkapscatmat_{specie_name}.inp') == True:
        dustkappa_filename = f'{path}dustkapscatmat_{specie_name}.inp'

    else:
        print(f'ERROR: kappa files are neither dustkappa_X.inp not dustkapscatmat.inp! Need to add additional file-name-style in a3d.load_onekappa()!')
        return 'STOPPING'

    # Declare lists and counters
    wavelengths = []
    absorption = []
    scattering = []
    scattangle = []
    nn = 0

    # Open and read dustkappa-files
    with open(dustkappa_filename, 'r') as f:
        for line in f.readlines():

            # Skip comments and header
            if line[0] != '#':
                # second line is number of wavelengths (not used)
                #if nn == 1:
                #    Nwave = int(line)
                
                # If line contains four elements: then it's wavelength, abs, scat, angle
                # If it's longer or just 1 then it's angles and polarisation matrix elements
                
                templine = line.split()
                Ndata = len(templine) - 1

                if Ndata > 0 and Ndata < 4:

                    wavelengths.append(float(templine[0]))
                    absorption.append(float(templine[1]))

                    if Ndata > 1:
                        scattering.append(float(templine[2]))
                    if Ndata > 2:
                        scattangle.append(float(templine[3]))
                    
                    # Save this Ndata also
                    finalNdata = Ndata
                
                # Increase line counter
                nn += 1
    
    if finalNdata == 1:
        return specie_name,wavelengths,absorption
    if finalNdata == 2:
        return specie_name,wavelengths,[absorption,scattering,np.array(absorption)+np.array(scattering)]
    if finalNdata == 3:
        return specie_name,wavelengths,[absorption,scattering,np.array(absorption)+np.array(scattering),scattangle]


# ------------------------------------------------------------ #
# Load output data from R3D

# Load dust_temperature-file
def load_temperature(
        path:str='../dust_temperature.dat',
        numb_specie:int=1
    ):
    """
    Load and extracts temperatures in output file dust_temperature.dat

    ARGUMENTS
      path: string with path and filename of temperature file
      nspecie: number of species you want to load
    
    RETURNS
      Ncells: Number of cells in grid (nleafs)
      Nspecies: number of species in data
      dust_temperatures: np.array containing temperatures in grid
    """

    # Read header
    with open(path,'r') as f:
        for nn,line in enumerate(f.readlines()):

            # Number of cells
            if nn == 1:
                Ncells = int(line)
            
            # Number of species
            if nn == 2:
                Nspecies = int(line)
    
    # Check that species number is correct
    if numb_specie <= Nspecies and numb_specie > 0:
        # Create density np.array
        dust_temperatures = np.zeros(Ncells)

        # Extract dust densities (reduce species number by 1, since python starts at index=0)
        ncounter = 0
        numb_specie -= 1

        with open(path,'r') as f:
            for nn,line in enumerate(f.readlines()):
                if (2 + Ncells*numb_specie) < nn < (2+ Ncells*(numb_specie+1)):
                    dust_temperatures[ncounter] = float(line)
                    ncounter += 1

        return Ncells,Nspecies,dust_temperatures
    
    else:
        # Otherwise, return error message
        print(f'ERROR: number of species in data is {Nspecies}, your specie-number was {numb_specie}')


# Load SED
def load_spectrum(
        path:str='../spectrum.out',
        distance:float=1
    ):
    """
    Loads and returns SED and wavelength.

    ARGUMENTS
      path: str with path and filename of spectrum.out file
      distance: float, distance to star in pc

    RETURNS
      wavelengths: list with wavelength grid in microns
      spectrum: list with SED in Jy normalised to indicated distance
    """

    # Declare lists
    wavelengths = []
    spectrum = []

    # Load data
    with open(path, 'r') as f:

        # Data starts at row 4, so first data at nn == 3, where each row has two numbers
        for nn,line in enumerate(f.readlines()):
            if nn > 2:

                # Extract wavelengths (in um)
                wavelengths.append(float(line.split('   ')[0]))

                # Extract and recompute flux density to Jy (at 1pc for now)
                spectrum.append(float(line.split('   ')[1]) * 1e23 / distance**2)
    
    # Return data
    return wavelengths,spectrum


# Prepare image data for plots
def load_images(
        path:str='../',
        image:str='image.out',
        distance:float=1
    ):
    """
    Loads and outputs array of one image in linear and log-scale, and some numbers.
    
    ARGUMENTS
      path: string with path to folder containing r3d-image files
      image: string with image file name
      distance: set distance to the source in pc, default 1pc
    
    RETURNS
      image2d: 2d-array with image
      image2dlog: 2d-array with image with logarithmic flux scale
      totalflux: total flux density of image in Jy
      axisplot: list with axis sizes, for matplotlib-plotting
    """

    # Automatically add / to end of path if it's missing
    if path[-1] != '/':
        path += '/'

    # Declare lists
    image1d = []

    # Load images
    with open(path+image, 'r') as f:
        for nl,line in enumerate(f.readlines()):
            
            # row 1: pixels by pixels
            if nl == 1:
                npixels_x = int(line.split()[0])
                npixels_y = int(line.split()[1])
                npixels = max([npixels_x,npixels_y])            
            # row 3: pixel size is in cm, divide by AUcm for AU
            if nl == 3:
                pixelsize_au = float(line.split()[0])/AUcm
            
            # NOTE might be useful later, commented out for now
            # row 4: wavelenght in um
            #if nl == 4:
            #    wavelength = float(line)

            # row 6 onward: pixels
            # image.out's pixels has unit
            # erg s-1 cm-2 Hz-1 ster-1
            if nl > 5:
                # Some rows are empty (and they contain space and \n, so strip them lines)
                if len(line.strip()) > 0:
                    image1d.append(float(line.strip()))

    # Extract some useful quantities
    # pixel size in asec (pixelsize in au and distance in pc gives distance in asec)
    pixelsize_as = pixelsize_au / distance
    
    # Size of whole image in AU and image-axis-scales
    size_au = pixelsize_au * npixels
    axisplot  = [-0.5*size_au,0.5*size_au,-0.5*size_au,0.5*size_au]

    # Total flux density of the image in Jy
    # Transoform to Jy/pix and sum all
    # 1 Jy = 1e23 erg/(s cm2 Hz)
    # 1 asec = 1/(180/pi * 3600)^2 ster = 2.35044305391e-11 ster
    # 1 pixel = pixelsize_as^2  asec^2
    totalflux = sum(image1d) * 1.e23 * 2.35044305391e-11 * pixelsize_as**2

    # Create 2D arrays
    image2d = np.zeros((npixels,npixels))
    image2dlog = np.zeros((npixels,npixels))
    nx,ny = 0,0

    for flux in image1d:
        # Convert image1d to 2d and change unit to Jy/asec2
        image2d[nx,ny] = flux * 1.e23 * 2.35044305391e-11
        
        # Remove zeros before logging (1e-6 is small enough, probably)
        # Smallest number >0 is 1.1, so log10 is just larger than 0
        if image2d[nx,ny] == 0:
            image2dlog[nx,ny] = -6
        else:
            image2dlog[nx,ny] = np.log10(flux * 1.e23 * 2.35044305391e-11)

        # Move nx and ny
        nx = nx + 1
        if nx == npixels_x:
            nx = 0
            ny = ny + 1

    return image2d,image2dlog,totalflux,axisplot


# Extract total dust mass from r3d-log files
def extract_dustmasses(
        log_files = [],
        Nspecies = 1,
        first_species = 1
    ):
    """
    Extract and sum dust masses from saved r3d-log file's headers

    ARGUMENTS
      log_files: list of file names of r3d-log files
      Nspecies: total number of species in log file
      first_species: numbering on the first dust species to include in sum
    
    RETURNS
      Prints file name and sum of specified dust species:
      first_species to Nspecies
    """
    for file in log_files:
        # Reset dust mass
        totalmass = 0

        # Load beginning of log-file
        with open(file, 'r') as flog:
            # Loop through beginning of log file
            # Line 37 contains first dust species
            for nline,line in enumerate(flog.readlines()):
                if nline >= 36+first_species and nline <= 36+Nspecies:

                    # Extract lines with dust masses
                    if line.split()[:4] == ['Dust', 'mass', 'in', 'species']:
                        # Dust mass is 6th "word of the line"
                        masses = float(line.split()[6])
                        totalmass += masses

                if nline > 36+Nspecies:
                    break
        # Output total mass and phase number
        print(f'{file}    {totalmass}')


# Load data from photocentre-files
def load_photocentre_file(
        file_path:str='../r3dresults/photocentre.dat'
    ):
    """
    Extract X-Y-R coordinates of photocentre out of photocentre.dat-files as created with 
    function a3d.write_photocentre_files from RADMC3D-images.

    ARGUMENTS
      file_path: string with path to (including) photocentre.dat-file
    
    RETURNS
      angles: list of included angles
      snapshots: list of included snapshot numbers
      coordinatelist: multidimensional list with arrays containing all coordinates for all angles
                      and snapshots:     coordinatelist[angle][snapshot , X/Y/R]
                      e.g. coordinatelist[0][snapshot , X/Y/R]  is Nsnapshots x 3 in size
    """

    coordinatelist = []
    snapshots = []
    line_counter = 0

    with open(file_path, 'r') as fphotoc:

        # Loop through file
        for nline,line in enumerate(fphotoc.readlines()):

            # Skip all comments
            if line[0] != '#':
                
                # Extract number of snapshots
                if line.split('=')[0] == 'Nsnapshots':
                    Nsnapshots = int(line.split('=')[1])
                
                # Extract angles
                if line[:4] == '    ':
                    # save angles in list
                    angles = line.split('    ')[1:]
                    # remove \n on last angle
                    angles[-1] = angles[-1].rstrip()
                    # Save number of angles
                    Nangles = len(angles)

                    # Create a list containing arrays for each coordinate
                    # for each angle
                    for _ in angles:
                        coordinatelist.append(np.zeros((Nsnapshots,3)))

                # Extract snapshot numbers and photocentre positions
                # for each angle
                if nline > 5:

                    # Separate out snapshotnumber and angle-dependent XYR-coords
                    linedata = line.split('    ')

                    # Save snapshot numbers
                    snapshots.append(int(linedata[0]))

                    # Save XYR-coords for each angle
                    for nangle in range(Nangles+1):
                        if nangle > 0:
                        
                            # Extract coordinates from list
                            pcX = linedata[nangle].split('  ')[0]
                            pcY = linedata[nangle].split('  ')[1]
                            pcR = linedata[nangle].split('  ')[2]

                            # Write to coordinate list
                            coordinatelist[nangle-1][line_counter,0] = pcX
                            coordinatelist[nangle-1][line_counter,1] = pcY
                            coordinatelist[nangle-1][line_counter,2] = pcR

                    # Data line counter (should be at the end, should start
                    # at zero)
                    line_counter += 1

    return angles,snapshots,coordinatelist

# ------------------------------------------------------------ #
# Functions to set Plot-settings

# TODO implement or remove these?

def set_figurefonts():

    # Import required packages
    from matplotlib import rc
    rc('font',**{'family':'serif','serif':['serif']})
    rc('text', usetex=True)

    # Set size of ticks and tick fonts, labels and scales.
    rc('xtick.major',size=8)
    rc('xtick.minor',size=4)
    rc('ytick.major',size=8)
    rc('ytick.minor',size=4)

    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)


# This function sets standard plot settings TODO: REMOVE THIS?
def set_ticklabels(
        xlabel:str,
        ylabel:str,
        xscale:str,
        yscale:str
    ):
    """
    TODO add info
    """

    # Import required packages
    from matplotlib import rc
    rc('font',**{'family':'serif','serif':['serif']})
    rc('text', usetex=True)

    # Set size of ticks and tick fonts, labels and scales.
    rc('xtick.major',size=8)
    rc('xtick.minor',size=4)
    rc('ytick.major',size=8)
    rc('ytick.minor',size=4)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.xlabel(xlabel,fontsize=18)
    plt.ylabel(ylabel,fontsize=18)
    if xscale == "log":
        plt.xscale('log')
    if yscale == "log":
        plt.yscale('log')

# ------------------------------------------------------------ #
# Plot various details of input and output data of R3D

# Plots the specified grid
def plot_grid(
        gridpath:str = '../grid_distances.csv',
        sizepath:str = '../grid_cellsizes.csv',
        amrpath:str = '../amr_grid.inp',
        nbins:int = 10
    ):
    """
    Loads and plots the current grid.

    ARGUMENTS
      gridpath: path to grid_distances.csv
      amrpath: path to amr_grid.inp
      nbins: number of bins for radial histogram plot
    """

    # Load data
    griddistances = load_griddistances(gridpath,amrpath)
    gridsizes = load_cellsizes(sizepath,amrpath)
    
    # Load some grid props
    nleafs = load_grid_properties()[2]
    ncellsizes = np.size(np.unique(gridsizes))

    # Change units to AU
    for nn in range(nleafs):
        gridsizes[nn] = gridsizes[nn]/AUcm
        for nx in range(4):
            griddistances[nn,nx] = griddistances[nn,nx]/AUcm

    # Plots (need a nice way to set fontsizes and fonts?)
    fig, ax = plt.subplots(3,2)

    # Plot radial distances
    ax[0,0].hist(griddistances[:,0],bins=nbins)
    ax[0,0].set(
        title='Radial distances', 
        ylabel='Number of cells', 
        xlabel='Radial to centrum of grid (AU)'
    )

    # Plot coordinates in each plane
    ax[0,1].plot(griddistances[:,1],griddistances[:,2],'b.')
    ax[0,1].set(
        title='Cells in X-Y-plane', 
        xlabel='X coord (AU)', 
        ylabel='Y coord (AU)'
    )

    ax[1,0].plot(griddistances[:,1],griddistances[:,3],'b.')
    ax[1,0].set(
        title='Cells in X-Z-plane', 
        xlabel='X coord (AU)', 
        ylabel='Z coord (AU)'
    )

    ax[1,1].plot(griddistances[:,2],griddistances[:,3],'b.')
    ax[1,1].set(
        title='Cells in Y-Z-plane', 
        xlabel='Y coord (AU)', 
        ylabel='Z coord (AU)'
    )

    ax[2,0].hist(gridsizes,bins=ncellsizes)
    ax[2,0].set(
        title='Cell sizes',
        xlabel='Cell size (AU)',
        ylabel='Number of cells'
    )

    # Better spacing between figures
    fig.tight_layout()
    fig.show()


def plot_onedensity_radius(
        density_path:str='../dust_density.inp',
        grid_path:str='../grid_distances.csv',
        amr_path:str='../amr_grid.inp',
        numb_specie:int=1
    ):
    """
    Plots one figure with radial density distribution of one dust species of a 
    dust_density-file of your choice.

    ARGUMENTS
      density_path: path to density.inp-file
      grid_path: path to grid_distances.csv'
      amr_path: path to amr_grid.inp
      numb_specie: number of the dust species you want to plot

    RETURNS
      Shows figure
    """

    # Load griddistances
    griddistances = load_griddistances(
        gridpath=grid_path,
        amrpath=amr_path
    )/AUcm

    # load dust_density
    density = load_dustdensity(
        path=density_path,
        numb_specie=numb_specie
    )[2]

    # Load and plots r3d density data for ONE dust specie
    fig, ax = plt.figure(), plt.axes()
    ax.plot(
        griddistances[:,0],density,
        linestyle='',marker='.',markersize=1
    )
    ax.set(
        ylabel=r'Density (g\,cm$^{-3}$)',
        xlabel=r'Distance (AU)',
        title=f'Dust species {numb_specie}'
    )
    fig.show()


def plot_alldensity_radius(
        path:str='../'
    ):
    """
    Plots one figure with radial density distribution of all species in dust_density.inp

    ARGUMENTS
      path: path to density.inp-file's folder

    RETURNS
      Shows figures
    """

    # Automatically add / to end of path if it's missing
    if path[-1] != '/':
        path += '/'

    # Get paths to necessary files
    density_path = path+'dust_density.inp'
    grid_path = path+'../grid_distances.csv'
    amr_path = path+'amr_grid.inp'
    
    # Get info on number of dust species in 

    # Load griddistances
    griddistances = load_griddistances(
        gridpath=grid_path,
        amrpath=amr_path
    )/AUcm

    # load first dust_density
    Ncells,Nspec,density = load_dustdensity(
        path=density_path,
        numb_specie=1
    )
    densities = [density]

    # Load the rest of the dust densities and put them in a list of np.arrays
    if Nspec > 1:
        for nn in range(2,Nspec+1):
            densities.append(
                load_dustdensity(
                    path=density_path,
                    numb_specie=nn
                )[2]
            )
        
    # Control colours of each density distribution
    colour = cm.rainbow(np.linspace(0, 1, Nspec))

    # Set objects for plot with all in the same figure
    fig, ax = plt.figure(), plt.axes()
    
    for nn, c in enumerate(colour):

        density = np.array(densities[nn])

        ax.plot(
            griddistances[np.where(density > 0)[0],0],density[np.where(density > 0)[0]],
            markeredgecolor=c,linestyle='',marker='.',markersize=1
        )
    ax.set(
        ylabel=r'Density (g cm$^{-3}$)',
        xlabel=r'Distance (AU)',
        title=f'Dust species 1 to {Nspec}'
    )
    fig.tight_layout()
    fig.show()

    # Plot of all species in same figure logarithmic yscale
    fig, ax = plt.figure(), plt.axes()
    
    for nn, c in enumerate(colour):

        density = np.array(densities[nn])

        ax.plot(
            griddistances[np.where(density > 0)[0],0],density[np.where(density > 0)[0]],
            markeredgecolor=c,linestyle='',marker='.',markersize=1
        )
    ax.set(
        ylabel=r'Density (g cm$^{-3}$)',
        xlabel=r'Distance (AU)',
        yscale='log',
        title=f'Dust species 1 to {Nspec}'
    )
    fig.tight_layout()
    fig.show()

    # Set objects for subplots (two columns, increasing number of rows)
    # ax[rows,columns]
    fig,ax = plt.subplots((-(-Nspec//2)),2)

    for nn, c in enumerate(colour):

        density = np.array(densities[nn])

        ax.ravel()[nn].plot(
            griddistances[np.where(density > 0)[0],0],density[np.where(density > 0)[0]],
            markeredgecolor=c,linestyle='',marker='.',markersize=1
        )
        ax.ravel()[nn].set(
            ylabel=r'Density (g cm$^{-3}$)',
            xlabel=r'Distance (AU)',
            title=f'Dust species {nn+1}'
        )
    fig.tight_layout()
    fig.show()


def plot_onetemperature_radius(
        temperature_path:str='../dust_temperature.dat',
        grid_path:str='../grid_distances.csv',
        amr_path:str='../amr_grid.inp',
        numb_specie:int=1
    ):
    """
    Plots one figure with radial temperature distribution of one specie.

    ARGUMENTS
      temperature_path: path to dust_temperature.dat
      grid_path: path to grid_distances.csv'
      amr_path: path to amr_grid.inp

    RETURNS
      figure-object
    """

    # Load griddistances
    griddistances = load_griddistances(
        gridpath=grid_path,
        amrpath=amr_path
    )/AUcm

    # load dust_temperature
    temperatures = load_temperature(
        path=temperature_path,
        numb_specie=numb_specie
    )[2]

    # Load and plots r3d density data for ONE dust specie
    fig, ax = plt.figure(), plt.axes()
    ax.plot(
        griddistances[:,0],temperatures,
        linestyle='',marker='.',markersize=1
    )
    ax.set(
        ylabel=r'Temperature (K)',
        xlabel=r'Distance (AU)',
        title='Grid cell temperatures'
    )
    return fig


# Plot all temperatures
def plot_alltemperature_radius(
        path:str='../'
    ):
    """
    Plots figures with radial temperature distribution of all species within temperature file.

    ARGUMENTS
      temperature_path: path to dust_temperature.dat
      grid_path: path to grid_distances.csv'
      amr_path: path to amr_grid.inp

    RETURNS
      figure-object
    """

    # Automatically add / to end of path if it's missing
    if path[-1] != '/':
        path += '/'

    # Get paths to necessary files
    temperature_path = path+'dust_temperature.dat'
    grid_path = path+'../grid_distances.csv'
    amr_path = path+'amr_grid.inp'

    # Load griddistances
    griddistances = load_griddistances(
        gridpath=grid_path,
        amrpath=amr_path
    )/AUcm

    # Load first dust_temperature
    Ncells,Nspec,temperature = load_temperature(
        path=temperature_path,
        numb_specie=1
    )
    temperatures = [temperature]

    # Load rest of temperatures
    if Nspec > 1:
        for numb_specie in range(2,Nspec+1):
            temperatures.append(
                load_temperature(
                    path=temperature_path,
                    numb_specie=numb_specie
                )[2]
            )
    
    # Control colours of each density distribution
    colour = cm.rainbow(np.linspace(0, 1, Nspec))

    # Set objects for plot with all in the same figure
    fig1, ax1 = plt.figure(), plt.axes()
    
    for nn, c in enumerate(colour):

        temperature = np.array(temperatures[nn])

        ax1.plot(
            griddistances[
                np.where(temperature > 0)[0],0
            ],temperature[
                np.where(temperature > 0)[0]
            ],
            markeredgecolor=c,linestyle='',marker='.',markersize=1
        )
    ax1.plot(
        1.65,2800,'r*'
    )
    ax1.set(
        ylabel=r'Cell temperature (K)',
        xlabel=r'Distance (AU)',
        title=f'Dust species 1 to {Nspec}'
    )
    fig1.tight_layout()
    fig1.show()

    # Set objects for subplots (two columns, increasing number of rows)
    # ax[rows,columns]
    fig2,ax2 = plt.subplots((-(-Nspec//2)),2)

    for nn, c in enumerate(colour):

        temperature = np.array(temperatures[nn])

        ax2.ravel()[nn].plot(
            griddistances[
                np.where(temperature > 0)[0],0
            ],temperature[
                np.where(temperature > 0)[0]
            ],
            markeredgecolor=c,linestyle='',marker='.',markersize=1
        )
        ax2.ravel()[nn].set(
            ylabel=r'Cell temperature (K)',
            xlabel=r'Radius (AU)',
            title=f'Dust species {nn+1}'
        )
    fig2.tight_layout()
    return fig1, fig2


# Plots temperature of radial bins (for paper-usage)
def plot_temperaturebins_radius(
        temperature_path:str='../dust_temperature.dat',
        grid_path:str='../grid_distances.csv',
        amr_path:str='../amr_grid.inp',
        numb_specie:int=1,
        ax=0
    ):
    """
    Plots average temperature of 100 spherical shells, 
    and max-min-values, and STD of each shell.
    Can take a few minutes when loading the larger data sets.
    Only returns fig and ax objects if input-ax = 0, is ax and axis-object
    this function also returns fig and ax-objects (see below RETURNS)

    ARGUMENTS
      temperature_path:str='../dust_temperature.dat',
      grid_path:str='../grid_distances.csv',
      amr_path:str='../amr_grid.inp',
      numb_specie:int=1
      ax: ax object to plot into, if none, don't change
          if you only want the binned arrays, write 'no' here.

    RETURNS
      (fig),(ax), temperature_bins,temperature_std,minmax_bins,radial_range
        figure-object (if ax-object did not exist before running func)
        axes-object (if ax-object did not exist before running func)
        Radial temperature-array - in K
        and standard deviation
        and Max-min
        Radial range array - in AU
    """

    # Load R3D-temperature-file
    Ncells, Nspecies, temperatures = load_temperature(
        path=temperature_path,
        numb_specie=numb_specie
    )

    # Load coordintaes of R3D-cells and change to AU
    cellcoords = load_griddistances(
        gridpath=grid_path,
        amrpath=amr_path
    )
    radii = cellcoords[:,0]/AUcm

    # Load star's radius here
    Mstar,Rstar,Lstar,Tstar = a5d.load_star_information(
        savpath='../../exwings_archivedata/co5bold_data/dst28gm06n052/st28gm06n052_186.sav',
        printoutput='n'
    )
    Rstar /= AUcm


    # Set up bins and arrays for binned data
    Nbins = 100
    radial_bins = np.linspace(0,radii.max(),Nbins+1)
    radial_range = np.zeros(Nbins)
    temperature_bins = np.zeros(Nbins)
    temperature_max = np.zeros(Nbins)
    temperature_min = np.zeros(Nbins)
    temperature_std = np.zeros(Nbins)

    # Loop through bins and save binned data
    for nn in range(Nbins):
        ncells = np.where((radii >= radial_bins[nn]) & (radii < radial_bins[nn+1]))[0]
        temperature_bins[nn] = temperatures[ncells].mean()
        temperature_max[nn] = temperatures[ncells].max()
        temperature_min[nn] = temperatures[ncells].min()
        temperature_std[nn] = temperatures[ncells].std()
        radial_range[nn] = radial_bins[nn] + 0.5 * radial_bins[1]

    # Check if input ax-object exists (!= 0)
    if ax == 0:
        # Create figure-ax-objects
        fig, ax = plt.figure(figsize=(6, 4)), plt.axes()
        # And save that fig-object is not created outside function
        externalfig = 'n'
    else:
        # There's an external fig, check if we want to plot at all
        externalfig = 'y'

        # If ax is not "no", then we plot the figure also
        # if ax is "no", then we don't plot anything and only the bins
        # are sent out
        if ax != 'no':

            ax.plot(radial_range,temperature_bins,'k')

            ax.fill_between(
                radial_range,
                temperature_min,
                temperature_max,
                color='b',
                alpha=0.2
            )

            ax.fill_between(
                radial_range,
                temperature_bins-temperature_std,
                temperature_bins+temperature_std,
                color='b',
                alpha=0.4
            )

            # Plot star's radius here
            ax.plot([Rstar,Rstar],[0,4000],'r:',linewidth=1)
            
            # Plot settings
            ax.set(
                ylabel=r'Gas and dust temperature (K)',
                xlabel=r'Distance (AU)',
                ylim=(0,4000),
                xlim=(0,radial_range.max()+0.5)
            )

    # if figure object is created outside, don't return any figs or ax
    # if no figure at all, just send data
    if externalfig == 'y':
        return temperature_bins,temperature_std,temperature_max,temperature_min,radial_range
    
    # Else return only fig and ax-objects
    if externalfig == 'n':
        return fig,ax
    


# Plot grain sizes of R3D cells against radius
def plot_grainsize_radius(
        gridpath:str='../r3dresults/st28gm06n052/grid_distances.csv',
        amrpath:str='../r3dresults/st28gm06n052/amr_grid.inp',
        grainsizepath:str='../grain_sizes_186.dat'
    ):
    """
    TODO info
    """

    # Load radius-array
    griddistances = load_griddistances(
        gridpath=gridpath,
        amrpath=amrpath,
    )
    radius_au = griddistances[:,0]/AUcm

    # Load grainsizes in um
    # Load grainsizes
    sizes,Nleafs = load_grainsizes(
        grainsize_path=grainsizepath
    )
    sizes_um = sizes*1e4

    # Plot
    fig, ax = plt.figure(), plt.axes()

    ax.plot(radius_au,sizes_um,
        '.', markersize = 1
    )
    ax.set(
        xlabel=r'Distance (AU)',
        ylabel=r'Grain size ($\mu$m)',
    )
    return fig, ax


# Plot absorption, scattering, and angles of the various species
def plot_onekappa(
        specie_name:str = '',
        specie_number:int = 0,
        path:str = '../'
    ):
    """
    Load and plots absorption/scattering/scattering angles of one dust specie
    
    INPUT
    specie_name: a string
    OR
    specie_number: an integer
    path: path to folder containing dustkappa and dustopac-files

    OUTPUT
    Figure.
    """

    # Load abs-scat-data
    specie_name,wavelengths,kappadata = load_onekappa(
        specie_name=specie_name,
        specie_number=specie_number,
        path=path
    )

    # Load and plots r3d density data for ONE dust specie
    fig, ax = plt.figure(), plt.axes()
    ax.plot(wavelengths,kappadata[0],'b')
    legendtext = ['Absorption']

    if len(kappadata) > 1:
        ax.plot(wavelengths,kappadata[1],'r')
        legendtext.append('Scattering')

    if len(kappadata) > 2:
        ax.plot(wavelengths,kappadata[2],'g')
        legendtext.append(r'Scattering angle: $\left< \cos \theta \right>$')

    ax.legend(legendtext)
    ax.set(
        ylabel=r'$\kappa_{\rm abs}$, $\kappa_{\rm scat}$ (cm$^2$/g)',
        xlabel=r'Wavelength ($\mu$m)',
        title=f'Kappa of {specie_name}',
        xscale='log',yscale='log'
    )
    fig.show();


# Plot absorption, scattering, and angles of the various species
def plot_allkappa(
        path:str='../'
    ):
    """
    TODO write more info
    Load and plots absorption/scattering/scattering angles of all species in dustopac 
    at folder of path
    
    ARGUMENTS
      path: path to folder containing dustkappa and dustopac-files

    RETURNS
      Figure and Axes-objects, plot with fig.show() 
    """

    # Automatically add / to end of path if it's missing
    if path[-1] != '/':
        path += '/'

    # Extract species names from dustopac-file
    counter = 5
    specie_names = []
    with open(f'{path}dustopac.inp', 'r') as f:
        for nn,line in enumerate(f.readlines()):
            
            # Specie-names are on line numbers 1+specie_number*4
            if nn == counter:
                specie_names.append(line.strip())
                counter +=4

    # Extract all opacities
    kappas = []
    Nkappa = []

    for specie_name in specie_names:
        specie_name,wavelengths,kappadata = load_onekappa(
            specie_name=specie_name,
            path=path
        )
        # Save number of data in each dataset
        Nkappa.append(len(kappadata))

        # Save all data of each dataset
        for kappa in kappadata:
            kappas.append(kappa)
    
    
    # Plot all data in subplots
    # One plot for absorption
    # one for scattering
    # and one final for scattering angles


    # List of subplot titles
    figuretitles = ['Absorption', 'Scattering', 'Extinction', 'Average scattering angle']

    # List of line markers
    from matplotlib.lines import Line2D
    markerlist = list(Line2D.markers)

    # List of colours
    import matplotlib.colors as mcolors
    colourlist = []
    for colour in list(mcolors.TABLEAU_COLORS):
        colourlist.append(colour.split(':')[1])


    counter = 0
    fig, ax = plt.subplots(
        max(Nkappa),1,
        num='abs_scat_angle.pdf',
        figsize=(6,13)
    )
    for nkappa in Nkappa:
        for nn in range(nkappa):

            linecounter = int(counter/nkappa)

            # Plot cyrves
            ax[nn].plot(
                wavelengths,
                kappas[counter + nn],
                color=colourlist[linecounter],
                linestyle='-',
                label='_nolegend_',
                zorder=1
            )

            # Plot symbols on each
            indexnumber = round(0.24*len(wavelengths))+30*linecounter
            ax[nn].plot(
                wavelengths[indexnumber],
                kappas[counter + nn][indexnumber],
                color=colourlist[linecounter],
                marker=markerlist[2*linecounter],
                zorder=5
            )

            ax[nn].set(
                title=figuretitles[nn],
                xscale='log',yscale='log'
            )
        counter += nkappa
    
    # Change scattering angles to linear scale
    if max(Nkappa) > 3:
        ax[3].set(yscale='linear')

    # Change species name strings:
    grainsizes_legend = []
    for grainsize in specie_names:
        # I use the exponent in my rounding to make sure that the decimals are correct
        grainsize = float(grainsize.split('_')[1])
        grainsizes_legend.append(rf'{grainsize:.3f} $\mu$m')
    ax[nkappa-1].legend(labels=grainsizes_legend)

    fig.tight_layout()

    return fig,ax



# Plot SED
def plot_sed(
        path:str='../spectrum.out',
        distance:float=1
    ):
    """
    Input:
    path: path to spectrum.out
    distance: distance to sources in pc

    Output:
    SED-figure and ax-objects
    fluxmax, wavemax: maximum flux in Jy and corresponding wavelength in micrometer    
    """

    # Load SED
    wavelengths,spectrum = load_spectrum(path=path,distance=distance)

    # Extract maximum flux and wavelength
    maxindex = np.argmax(spectrum)
    
    # plot SED
    fig, ax = plt.figure(), plt.axes()
    ax.plot(
        wavelengths,spectrum,'b'
    )
    ax.set(
        ylabel=f'Flux density (Jy at {distance} pc)',
        xlabel=r'Wavelength ($\mu$m)',
        title='Output SED',
        xscale='log',yscale='log'
    )

    return fig,ax,spectrum[maxindex],wavelengths[maxindex]


# Function to plot a list of SEDs ontop of eachother in the same plot
def plot_sedsmany(
        pathlist:list=['../spectrum.out'],
        legendlist:list=['spectrum1'],
        distance:float=1
    ):
    """
    Input:
    pathlist: list of paths to all spectrum.out files
    distance: distance to source in pc

    Output:
    SED-figure and ax-objects
    """

    fig, ax = plt.figure(), plt.axes()

    for path in pathlist:

        # Load SED
        wavelengths,spectrum = load_spectrum(path=path,distance=distance)

        # plot SED
        ax.plot(
            wavelengths,spectrum
        )

    ax.set(
        ylabel=f'Flux density (Jy at {distance} pc)',
        xlabel=r'Wavelength ($\mu$m)',
        title='Output SEDs',
        xscale='log',yscale='log'
    )
    ax.legend(legendlist)

    return fig,ax


# Plot images
def plot_images(
        path:str='../',
        images:list=['image.out'],
        distance:float=1
    ):
    """
    Plots list of images with Z-scale following automatic gamma function

    ARGUMENTS
      path: path to folder containing images
      images: list of image-file names inside your folder
      distance: distance to source in pc (default 1 pc)

    RETURNS
      Figure with 1 plot with automatically gamma corrected
      Total flux of images in Jy at chosen distance
    """

    # Initiate list for flux densities of images
    fluxtotal = []

    # Loop through all images
    for image in images:

        image2d,image2dlog,flux,axisplot = load_images(
            path=path,
            image=image,
            distance=distance
        )
        # Save all total flux densities
        fluxtotal.append(flux)

        # Compute and apply gamma function of each image
        scale = np.max(image2d)-np.min(image2d)
        gamma = 0.3*np.log(image2d.max())/np.log(image2d.mean())
        imageplot = ((image2d / scale)**gamma) * scale

        # Plot each image
        fig, ax = plt.figure(dpi=150), plt.axes()

        im0 = ax.imshow(
            imageplot, origin='lower', extent=axisplot, cmap=plt.get_cmap('hot')
        )
        ax.set(
            title=f"{image.replace('image_', '').replace('.out', '')}", 
            xlabel='Offset (AU)',
            ylabel='Offset (AU)'
        )
        fig.tight_layout()

        return fig, ax, fluxtotal



def plot_imageslinlog(
        path:str='../',
        images:list=['image.out'],
        distance:float=1
    ):
    """
    Plots list of images

    ARGUMENTS
      path: path to folder containing images
      images: list of image-file names inside your folder
      distance: distance to source in pc (default 1 pc)

    RETURNS
      Figure with 2 subplots, linear and logarithmic scales.
      Total flux of images in Jy at chosen distance
    """

    fluxtotal = []

    # Loop through all images
    for image in images:

        image2d,image2dlog,flux,axisplot = load_images(
            path=path,
            image=image,
            distance=distance
        )

        fluxtotal.append(flux)


        # Plot each image (one in linear and one log)
        fig, ax = plt.subplots(
            1,2, 
            dpi = 150, 
            num = path+image
        )
        
        im0 = ax[0].imshow(
            image2d, origin='lower', extent=axisplot, cmap=plt.get_cmap('hot')
        )
        ax[0].set(
            title=f"{path[-4:-1]}: {image.replace('image_', '').replace('.out', '')} (Lin)", 
            xlabel='Offset (AU)',
            ylabel='Offset (AU)'
        )
                
        ax[1].imshow(
            image2dlog, origin='lower', extent=axisplot, cmap=plt.get_cmap('hot')
        )
        ax[1].set(
            title=f"{path[-4:-1]}: {image.replace('image_', '').replace('.out', '')} (Log)", 
            xlabel='Offset (AU)',
        )

        #cb0 = plt.colorbar(im0, orientation = 'vertical',shrink=0.4,pad=0.15)
        #cb0.set_label(label = rf'Flux at {distance} pc (Jy/asec$^2$)',fontsize= 10)

        # Change figure size
    fig.tight_layout()

    return fig, ax, fluxtotal


# Plot a list of images in subplots (vertically)
def plot_imagesubplots(
        imagelist:list=['../image.out'],
        distance:float=1
    ):
    """
    Plots a list of images from R3d in subplots (max 3 in one row)

    ARGUMENTS
      imagelist: list of paths to .out files
      distance: distance to source in pc

    RETURNS
      gives matplotlib objects: fig and ax, plot with fig.show()
    """

    # Number of plots
    Nplots = len(imagelist)

    # Set number of rows and columns in subplots (max 3 cols)
    Nrows = int(np.ceil(Nplots/3))
    if Nplots <= 3:
        Ncols = Nplots
    if Nplots > 3:
        Ncols = 3


    # Set fig-and-axis settings for subplots
    fig, ax = plt.subplots(
        Nrows,Ncols, 
        dpi = 100, 
        figsize = (4*Ncols,4*Nrows),
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
        image2d,image2dlog,flux,axisplot = load_images(
            path=f'{path}{modelname}/{phase}',
            image=imagefilename,
            distance=distance
        )
        # Compute and apply gamma function of each image
        scale = np.max(image2d)-np.min(image2d)
        gamma = 0.3*np.log(image2d.max())/np.log(image2d.mean())
        imageplot = ((image2d / scale)**gamma) * scale

        # Plot image at spot nn, set title and axis labels
        im0 = ax.ravel()[nn].imshow(
            imageplot, 
            origin='lower', extent=axisplot, 
            cmap=plt.get_cmap('hot')
        )
        ax.ravel()[nn].set(
            title=f'i:{incl}, phi:{phi}, {wavelengthum} um, {int(np.round(flux/1e6))} MJy', 
            xlabel='Offset (au)',
            ylabel='Offset (au)',
        )

    # Return fig and ax
    return fig, ax






def plot_imagecrosssections(
        path:str='../',
        images:list=['image.out'],
        distance:float=1
    ):
    """
    Plots cross sections through center along x and y-axes of images

    INPUT
    path: path to folder containing images
    images: list of image-file names inside your folder
    distance: distance to source in pc (default 1 pc)

    OUTPUT
    Figure with 4 subplots, 2 linear and 2 logarithmic scales.
    """

    fluxtotal = []

    # Loop through all images
    for image in images:

        # Declare list
        image1d = []

        # Load images
        with open(path+image, 'r') as f:
            for nl,line in enumerate(f.readlines()):
                
                # row 1: pixels by pixels
                if nl == 1:
                    npixels = int(line.split()[0])
                
                # row 3: pixel size in cm, divide by AUcm for AU
                if nl == 3:
                    pixelsize_au = float(line.split()[0])/AUcm
                
                # row 4: wavelenght in um
                # NOTE might be useful also :)
                #if nl == 4:
                #    wavelength = float(line)
                
                # row 6 onward: pixels
                if nl > 5:
                    # Some rows are empty (and they contain space and \n, so strip them lines)
                    if len(line.strip()) > 0:
                        image1d.append(float(line.strip()))

        # Extract some useful quantities
        # pixel size in mas
        pixelsize_mas = pixelsize_au / distance
        
        # Size of whole image in AU
        size_au = pixelsize_au * npixels
        axisplot  = np.linspace(0.5*size_au,-0.5*size_au,npixels)

        # Total flux density of the image
        fluxtotal.append(sum(image1d) * 1.e23 * 2.35044305391e-11 * pixelsize_mas**2)

        # Create 2D arrays
        image2d = np.zeros((npixels,npixels))
        image2dlog = np.zeros((npixels,npixels))
        nx,ny = 0,0

        for flux in image1d:
            image2d[nx,ny] = flux * 1.e23 * 2.35044305391e-11 * pixelsize_mas**2
            image2dlog[nx,ny] = np.log10(flux * 1.e23 * 2.35044305391e-11 * pixelsize_mas**2)

            # Move nx and ny
            nx = nx + 1
            if nx == npixels:
                nx = 0
                ny = ny + 1

        half_npixels = int(np.round(0.5*npixels))

        # Plot each image (one in linear and one log)

        # Along X-axis
        fig, ax = plt.subplots(2,2, dpi=150)

        ax[0][0].plot(
            axisplot,image2d[:,half_npixels]
        )
        ax[0][0].set(
            title=f"X-cross at {image.replace('image_', '').replace('.out', '')} (Lin)", 
            xlabel='Offset (AU)',
            ylabel='Flux density (Jy at 1 pc)'
        )

        ax[0][1].plot(
            axisplot,image2d[:,half_npixels]
        )
        ax[0][1].set(
            title=f"X-cross at {image.replace('image_', '').replace('.out', '')} (Log)", 
            yscale='log',
            xlabel='Offset (AU)',
            ylabel='Flux density (Jy at 1 pc)'
        )

        # Along Y-axis
        ax[1][0].plot(
            axisplot,image2d[half_npixels,:]
        )
        ax[1][0].set(
            title=f"Y-cross at {image.replace('image_', '').replace('.out', '')} (Lin)", 
            xlabel='Offset (AU)',
            ylabel='Flux density (Jy at 1 pc)'
        )

        ax[1][1].plot(
            axisplot,image2d[half_npixels,:]
        )
        ax[1][1].set(
            title=f"Y-cross at {image.replace('image_', '').replace('.out', '')} (Log)", 
            yscale='log',
            xlabel='Offset (AU)',
            ylabel='Flux density (Jy at 1 pc)'
        )

        fig.tight_layout()
        fig.show()


def plot_opticalthick(
        path:str = '../',
    ):
    """
    Plots optical thickness along an average line-of-sight and returns an estimate of
    the average radius of the star.

    ARGUMENTS
      path:str = path to folder containing all r3d-model-data.

    RETURNS
      star_surface:float = average radial distance to stellar surface in cm
      fig,ax = matplotlibobjects for plotting
    """

    # Automatically add / to end of path if it's missing
    if path[-1] != '/':
        path += '/'

    # Load densities (first one and then the rest to get the Nspecies)
    print('Loading densities')
    Ncells, Nspecies, density = load_dustdensity(
        path=path+'dust_density.inp',
        numb_specie=1,
    )
    densities = np.zeros((Ncells,Nspecies))
    densities[:,0] = density
    if Nspecies > 1:
        for nspecie in range(1,Nspecies):
            Ncells, Nspecies, density = load_dustdensity(
                path=path+'dust_density.inp',
                numb_specie=nspecie+1,
            )
            # Add all densities to one 2D array
            densities[:,nspecie] = density

    # Load Absorptions
    print('Loading opacities')
    kappa = 0
    kappas = np.zeros(Ncells)
    specieindeces = []
    for nspecie in range(Nspecies):
        specie_name,wavelengths,kappadata = load_onekappa(
            specie_number=nspecie+1,
            path=path
        )
        kappadata = np.array(kappadata)

        # Save maximum opacity of each specie
        kappa = np.max(np.array(kappadata[0] + kappadata[1]))

        # Create an Ncells long array with the opacities of each species at correct cell index
        # For this I can use the densities (where they are not 0)
        specieindeces.append(np.where(densities[:,nspecie] > 0)[0])
        
        # Save kappa-values and multiply with corresponding densities
        kappas[specieindeces[nspecie]] = kappa * densities[[specieindeces[nspecie]],nspecie]


    # Load griddistances [radial,x,y,z]
    print('Loading grid-info')
    griddistances = load_griddistances(
        gridpath=path+'../grid_distances.csv',
        amrpath=path+'../amr_grid.inp'
    )


    # Set up arrays for final optical depth computations
    print('Computing optical depths')
    # size of each spherical shell needs to be ~0.05 AU to work
    Nradius = int(np.round(np.max(griddistances[:,0])/(AUcm*0.05)))
    dxarray = np.linspace(0,np.max(griddistances[:,0]),Nradius)
    dtauarray = np.zeros(Nradius)

    for nx in range(2,np.shape(dxarray)[0]+1):
        nx *= -1

        # Extract indeces for each bin
        binindeces = np.where(
            (dxarray[nx] < griddistances[:,0]) & (griddistances[:,0] <= dxarray[nx+1])
        )

        # Average optical depth difference for each bin is then
        if np.isnan(np.mean(kappas[binindeces])) == True:
            dtauarray[nx] = 0
        else:
            dtauarray[nx] = np.mean(kappas[binindeces]) * (dxarray[nx+1] - dxarray[nx])


    # Optical depth along line of sight is then the sum of all outer values
    #dtauarray[-1] = dtauarray[-2]
    dtauarray[-1] = 0
    for nx in range(1,np.shape(dxarray)[0]+1):
        nx *= -1
        dtauarray[nx] += dtauarray[nx+1]


    # Find surface of star (if possible)
    if len(np.where(dtauarray < 1)[0]) > 0:
        star_surface_index = np.where(dtauarray < 1)[0][0]
        star_surface = 0.5* (dxarray[star_surface_index-1] + dxarray[star_surface_index])
    else:
        star_surface_index = 0
        star_surface = 0

    # Plot average optical thickness along LOS
    print('Creating figure')
    fig, ax = plt.figure('Optical thickness along average LOS', figsize=(6, 4)), plt.axes()

    # Optical thickness
    ax.plot(
        dxarray/AUcm,dtauarray,'b.'
    )

    # Stellar radius
    if star_surface_index > 0:
        ax.plot(
            [star_surface/AUcm,star_surface/AUcm],[dtauarray[0],dtauarray[-1]],'r'
        )
        ax.plot(
            [dxarray[star_surface_index-1]/AUcm,dxarray[star_surface_index]/AUcm],[1,1],'r'
        )

    ax.set(
        ylabel=r'$\tau$',
        xlabel=r'Distance along LOS (AU)',
        yscale='log'
    )
    fig.tight_layout()

    print('Returning: star_surface, fig, ax')
    return star_surface,fig,ax


# Function to plot and compute Pearson-Chi2-comparison between two sets of data
def plot_chisquare(
        simulation:list,
        observation:list,
        obssigma:list,
        xdata:list
    ):
    """
    PLOTS
      chi2 = (simulation - observation)^2 / obssigma^2

    ARGUMENTS
      simulation: list/array with simulated data to compare with observation
      observations: list/array with observed data
      obssigma: List/array with error bar data for observations
      Length of both must be the same

    RETURNS
      fig,ax,chisq_array,chiaq_reduced
        fig-object and ax-object
        The Chi^2-array that is plotted
        Total chi^2-number
          1/N * sum( (simulation - observation)^2 / obssigma^2 )
          where N = length of input arrays.
    """

    chisq_array, chiaq_reduced =  compute_chisquare(
        simulation,
        observation,
        obssigma
    )

    fig, ax = plt.figure(figsize=(6, 4)), plt.axes()

    # Chi2-array:
    ax.plot(xdata,chisq_array)

    # Line to singify 1-sigma limit
    ax.plot([xdata[0],xdata[-1]],[1,1],'r--')

    ax.set(
        ylabel=r'$\chi ^2$',
        xlabel=r'X-data',
        yscale='log',
        xscale='log',
        ylim=[1e-3,chisq_array.max()]
    )
    fig.tight_layout()

    return fig,ax,chisq_array,chiaq_reduced


# ------------------------------------------------------------ #
# Compute different quantities

def compute_sed_luminosity(
        path:str = '../spectrum.out',
        distance:float = 1
    ):
    """
    Insert a spectrum.out from r3d and get the bolometric luminosity in Watt
    Cuts wavelength grid below 0.3 um because of some strange UV-bumps. Doesnt
    change the lum much however.

    ARGUMENTS
      path = path to spectrum-file including file name
      distance = distance the SED is observed from

    RETURNS
      luminosity in Watt
    """

    # Load spectrum
    wavelengths_temp,spectrum_temp = load_spectrum(path,1)
    wavelengths_temp = np.array(wavelengths_temp)
    spectrum_temp = np.array(spectrum_temp)
    # Restric wavelength range from 0.3um due to UV-bump in some
    # spectra
    wavelengths = wavelengths_temp[np.where(wavelengths_temp >= 0.3)[0]]
    spectrum = spectrum_temp[np.where(wavelengths_temp >= 0.3)[0]]
    nwave = len(wavelengths)

    # Integrate the SED (using trapezoidal method, and change units to SI units)
    sedintegral = 0
    for nn in range(nwave-1):
        # 1.499e-12 = 0.5 * 1e-26 * 1e6 * c which are the corrections for units and the trapezoid-half.
        # Wavelength is summed in reversed order because its a sum over frequency
        sedintegral += (spectrum[nn] + spectrum[nn+1]) * (1/wavelengths[nn] - 1/wavelengths[nn+1])*1.499e-12

    # Compute bolometric luminosity
    luminosity = 4.*np.pi*(distance*pc)**2. * sedintegral

    return luminosity


def compute_luminosity(
        wavelengths:list,
        spectrum:list,
        distance:float = 1
    ):
    """
    ARGUMENTS
      wavelengthum: List of wavelengths in um
      spectrum: List of flux densities in Jy
      distance: distance to source in pc

    RETURNS
      luminosity in Watt
    """

    nwave = len(wavelengths)

    if nwave == len(spectrum):
        # Integrate the SED (using trapezoidal method, and change units to SI units)
        sedintegral = 0

        for nn in range(nwave-1):
            # 1.499e-12 = 0.5 * 1e-26 * 1e6 * c which are the corrections for units and the trapezoid-half.
            # Wavelength is summed in reversed order because its a sum over frequency
            sedintegral += (spectrum[nn] + spectrum[nn+1]) * (1/wavelengths[nn] - 1/wavelengths[nn+1])*1.499e-12

        # Compute bolometric luminosity
        luminosity = 4.*np.pi*(distance*pc)**2. * sedintegral

        return luminosity

    else:
        raise ValueError('ERROR, wavelengths and spectrum have different lengths')


# Function to plot and compute Pearson-Chi2-comparison between two sets of data
def compute_chisquare(
        simulation:list,
        observation:list,
        obssigma:list,
    ):
    """
    Computes chi2-numbers for comparisons betweens arrays/models/data.

    ARGUMENTS
      simulation: list/array with simulated data to compare with observation
      observations: list/array with observed data
      obssigma: List/array with error bar data for observations
        Length of all input lists/arrays must be the same

    RETURNS
      chisq_array,chiaq_reduced
        chisq_array: chi2 for each point in xdata-array
        chiaq_reduced: Total chi^2-number (ie reduced sum):
            1/N * sum( (simulation - observation)^2 / obssigma^2 )
            where N = length of input arrays.
    """

    # Check so that length of both are correct
    if len(simulation) != len(observation):
        raise ValueError('ERROR: your input data arrays are not of equal length')

    # Put data in np.arrays if needed
    if type(simulation) == list:
        simulation = np.array(simulation)
    if type(observation) == list:
        observation = np.array(observation)
    if type(obssigma) == list:
        obssigma = np.array(obssigma)

    # Compute chi2-array and total number
    chisq_array = ((simulation - observation) / obssigma)**2
    chiaq_reduced = 1/chisq_array.size * chisq_array.sum()

    return chisq_array,chiaq_reduced


# Function to load and translate filter to r3d-grid
# and normalise it to integral=1
def translate_filter(
        filterpath:str,
        wavelengthr3d:str
    ):
    """
    Translate filterprofiles as downloaded from SVO-page
    http://svo2.cab.inta-csic.es/svo/theory/fps3/index.php
    for usage on Radmc3d-computed SEDs.
    Also see sect 8.5 of RADMC3D manual, pages 90-91 for v2.0.

    ARGUMENTS
      filterpath: string with path to filter file from current folder
      wavelengthr3d: string to r3d-wavelengthfile, wavelength_micron.inp

    RETURNS
      file with filter within r3d-grid and filter-integral normalised to 1
    """
    # Extract filterfolder and filename
    filterfile = filterpath.split('/')[-1]
    filterfolder = ''
    for fold in filterpath.split('/')[:-1]:
        filterfolder += f'{fold}/'

    # Load r3d-wavelengths
    r3dwaves = np.loadtxt(wavelengthr3d)[1:]

    # Load filter profile and wavelengths translated to um from Å
    filterwave = np.loadtxt(filterpath)[:,0]*1e-4
    filterorig = np.loadtxt(filterpath)[:,1]

    # Interpolate to r3d-grid
    r3dfilter = np.interp(
        x=r3dwaves,
        xp=filterwave,
        fp=filterorig
    )

    # Normalise filter to 1
    r3dfilter = r3dfilter/np.trapz(r3dfilter,r3dwaves)

    # Save in new files, r3d-um-wavelength and filter side by side
    with open(f'{filterfolder}r3d_{filterfile}', 'w') as fr3d:
        # Write header
        fr3d.writelines(f'# {filterfile} band profile translated to radmc3d-grid\n# and integral normalised to 1\n')
        for r3dwave,r3dfilt in zip(r3dwaves,r3dfilter):
            fr3d.writelines(f'{r3dwave}    {r3dfilt}\n')


# Function to remove those strange spikes we're getting by combining two or more
# SEDs
def remove_sedspikes(
        paths:list = ['../r3dsims/spectrum.out'],
        save:bool=False
    ):
    """
    Combine spectra to remove random spikes.
    Obviously the spectra should come with the same physical settings, but with
    variations in e.g. number of cores or photons.

    ARGUMENTS
        paths: list of paths to the SEDs you want to combine into one
        save: False or True if you want to save the combined SED into a new file
              This file is put in the subfolder of your script.

    RETURNS
        wavelengths: list or array with wavelength grid in um
        sed_final: list or array with combined SEDs flux densities in Jy
    """

    # Declare list of SEDs
    seds = []

    # Load SEDs, uses default distance of 1pc
    for path in paths:
        seds.append(load_spectrum(path=path)[1])
    
    # And load wavelenths
    wavelength = load_spectrum(path=paths[0])[0]

    # Loop through wavelengths and all SEDs and extract min-values at each wavelength
    sed_final = []

    for nn in range(len(seds[0])):
        all_fluxes = []

        for sed in seds:
            all_fluxes.append(sed[nn])

        sed_final.append(np.min(all_fluxes))

    # Save as new sed-out-file
    #
    # Format:
    #            1
    #         1000
    #  
    #   0.1000000000000E+00   0.1238187064847E-23
    #    lambda in um          flux in Jy *1e-23 * distance[pc]**2  (ie normalised to 1pc)
    #
    if save == True:
        with open('../spectrum.out','w') as f:

            # Write header
            f.write(f'           1\n        1000\n\n')

            # Write wavelength and flux
            for nn in range(len(seds[0])):
                f.write(f'  {wavelength[nn]}   {sed_final[nn]*1e-23}\n')

    # Return final SED and wavelengths
    return wavelength,sed_final

# Function to remove spikes from images
def remove_imagespikes(
        folders:list=['../'],
        imagefilename:str='image.out'
    ):
    """
    Combines images made from different random seeds to remove flux spikes

    ARGUMENTS
      folders: list of strings with paths to folder containing image files
      imagefilename: string with file name of image you want to clean

    RETURNS
      file with new corrected image
      image1d: 1D np.array with flux densities of image
    """

    # Number of images to combine
    Nimage = len(folders)

    # Delcare Lists
    images = []
    npix = []
    header = []

    # Load images and put in list
    for nfolder,folder in enumerate(folders):
        # Create image pixel flux list
        image1d = []

        # line 1 is: npix-x    npix-y
        # line >5 are pixel fluxes
        with open(folder+imagefilename,'r') as imagefile:
            for nline,line in enumerate(imagefile.readlines()):
                
                # Save header lines
                if nline < 6 and nfolder == 0:
                    header.append(line)

                # Extract resolution
                if nline == 1:
                    npix.append(int(re.findall('\d+', line)[0]))

                # Extract fluxes
                if nline > 5:
                    if len(line.strip()) > 0:
                        image1d.append(float(line.strip()))

        # error check
        if nfolder > 0:
            if npix[nfolder] != npix[nfolder-1]:
                raise ValueError('    ERROR: your images are not the same size')

        # Save each image as lilsts in one list
        images.append(image1d)


    # Loop through all images simulatneously and save the minimum flux density of each
    # in an array
    newimagefluxes = np.zeros(npix[0]**2)

    for nn in range(npix[0]**2):

        # Reset list for pixelwise flux density of each image
        comparefluxes = []

        for nimage in range(Nimage):
            comparefluxes.append(images[nimage][nn])
        
        # Save minimum flux density of each image pixel
        newimagefluxes[nn] = np.min(comparefluxes)


    # Save new image
    with open(f'../{imagefilename}', 'w') as fnewimage:
        
        # First write copy of header
        for line in header:
            fnewimage.writelines(line)

        # Then write fluxes
        for newflux in newimagefluxes:
            fnewimage.writelines(f'   {newflux}\n')

    # Return image1D, np-array with flux densities in cgs-units at 1pc
    return newimagefluxes


# Loads two SEDs and returns comparisons between them
def compare_seds(
        sedpaths:list=[
            '../spectrum0.out',
            '../spectrum1.out'
        ]
    ):
    """
    general info TODO

    ARGUMENTS
      - sedpaths:list containing paths to two separate spectrum.out files from radmc3d
    RETURNS
      - radmc3d-output style file containing (sed1-sed0)/sed0 for all wavelengths
      - Prints Pearson's chi2-number defined as (log(sed1) - log(sed0))^2/log(sed0)
        with log10 on the flux densities to get numbers that are more easily compared
        with eachother and we for example assume
            sed0 = expected (original model)
            sed1 = model to compare with
    """
    print('Running: a3d.compare_seds')

    # Check if number of included SEDs are correct
    if len(sedpaths) != 2:
        raise ValueError('ERROR: there are not 2 SEDs included; can only compare two.')

    else:
        # Extract file names
        spectrum0 = re.split('/', sedpaths[0])[-1]
        spectrum1 = re.split('/', sedpaths[1])[-1]

        # Load SEDs and change to np.arrays
        wavelengths, sed0 = load_spectrum(
            path = sedpaths[0]
        )
        sed0 = np.array(sed0)
        wavelengths, sed1 = load_spectrum(
            path = sedpaths[1]
        )
        sed1 = np.array(sed1)

        Nlambda = len(wavelengths)

        # Compute comparisons
        comparison = (sed1-sed0)/sed0
        pearsonchi = np.sum( (np.log10(sed1) - np.log10(sed0))**2/np.log10(sed0) )

        # Save new SED-file
        with open(f'../{spectrum0[:-4]}_wrt_{spectrum1[:-4]}.out','w') as f:

            # Write header
            f.write(f'           1\n        {Nlambda}\n\n')

            # Write wavelength and flux
            for nn in range(Nlambda):
                f.write(f'  {wavelengths[nn]}   {comparison[nn]}\n')

        print(f'    {spectrum0[:-4]}_wrt_{spectrum1[:-4]}.out\nDONE\n')
        return pearsonchi


# Compares two images
def compare_images(
        imagepaths:list=[
            '../image0.out',
            '../image1.out'
        ]
    ):
    """
    Subtracts fluxes of two images and saves residuals, either cut to a minimum
    of 0 per pixel or with absolute vlue taken.

    ARGUMENTS
      imagepaths: list with two strings containing paths to the two images you
                  want to compare
    
    RETURNS
      Two new image files in radmc3d-formats:
        ../residual_zerolimit_{imagename} : residual of path[0] minus path[1] image
                                            with fluxes cut at 0.
        ../residual_abs_{imagename} : residual of  path[0] minus path[1] image with
                                      absolute value taken of all fluxes
    """
    print('Running: a3d.compare_images()')

    # Extract image file names
    image0 = re.split('/', imagepaths[0])[-1]
    image1 = re.split('/', imagepaths[1])[-1]

    # Delcare Lists
    images = []
    npix = []
    header = []



    # Load images and put in list
    for nimage,imagepath in enumerate(imagepaths):
        # Create image pixel flux list
        image1d = []

        # line 1 is: npix-x    npix-y
        # line >5 are pixel fluxes
        with open(imagepath,'r') as fimage:
            for nline,line in enumerate(fimage.readlines()):
            
                # Save header lines of first image
                if nline < 6 and nimage == 0:
                    header.append(line)

                # Extract resolution
                if nline == 1:
                    npix.append(int(re.findall('\d+', line)[0]))

                # Extract pixel fluxes
                if nline > 5:
                    if len(line.strip()) > 0:
                        image1d.append(float(line.strip()))

        # error check: both images must have same resolution!
        if nimage > 0:
            if npix[nimage] != npix[nimage-1]:
                raise ValueError('    ERROR: your images are not the same size')

        # Save each image as lists in one list
        images.append(image1d)


    # Extract residuals of images, one cut at 0 and one absolute
    #
    # declare array for both residuals
    residualfluxes = np.zeros((2,npix[0]**2))

    # Loop through image pixels
    for nn in range(npix[0]**2):

        # Reset list for pixelwise flux density of each image
        comparefluxes = []

        for nimage in range(2):
            comparefluxes.append(images[nimage][nn])
        
        # Save image1 minus image2 at each pixel
        # Set to zero if below zero for one
        # Set to absolute for the other
        #
        residual = comparefluxes[0] - comparefluxes[1]
        if residual < 0:
            residualfluxes[0,nn] = 0
        else:
            residualfluxes[0,nn] = residual

        # 1: and absolute
        residual = np.abs(comparefluxes[0] - comparefluxes[1])
        residualfluxes[1,nn] = residual

    # Save new image
    with open(f'../image_resid_zerolim_{image0[:-4]}_minus_{image1[:-4]}.out', 'w') as fnewimage:
        
        # First write copy of header
        for line in header:
            fnewimage.writelines(line)

        # Then write fluxes
        for newflux in residualfluxes[0,:]:
            fnewimage.writelines(f'   {newflux}\n')
    # Save new image
    with open(f'../image_resid_abs_{image0[:-4]}_minus_{image1[:-4]}.out', 'w') as fnewimage:
        
        # First write copy of header
        for line in header:
            fnewimage.writelines(line)

        # Then write fluxes
        for newflux in residualfluxes[1,:]:
            fnewimage.writelines(f'   {newflux}\n')

    print(f'A3D compare images:\n    image_resid_zerolim_{image0[:-4]}_minus_{image1[:-4]}.out\n    image_resid_abs_{image0[:-4]}_minus_{image1[:-4]}.out\nDONE\n')


# Compute a black body function based on SED peak flux and wavelength
def compute_blackbody_freq(
        peak_flux:float = 1e6,
        peak_wavelength:float = 1e6,
        wavelengths:list = [0.1,1,10,100]
    ):
    """
    Compute a frequency based black body function (of wavelength)

    ARGUMENTS
      peak_flux: peak flux in any units, this is used to normalise the output BB to these
                 units.
      peak_wavelength: in microns! Corresponding wavelength of peak flux.
      wavelengths: in microns! Array with wavelength grid.

    RETURNS
      Array with black body function within given wavelength grid normalised to
      given peak flux.
    """

    # Check type of wavelength array
    if type(wavelengths) == list:
        wavelengths = np.array(wavelengths)

    # Translate wavelengths to frequency
    peak_freq = c/(peak_wavelength*1e-6)
    frequencies = c/(wavelengths*1e-6)

    # Compute temperature from Wiens displacement law
    BBtemperature = peak_freq / 5.879e10

    # Compute black body function
    BBfreq = 2*hplanck*frequencies**3/c**2 * (np.exp( hplanck*frequencies/(kboltz*BBtemperature)) -1 )**-1

    # Normalise BB to input units
    BBfit = peak_flux/np.max(BBfreq) * BBfreq

    return BBfit


# Load image and compute position of photocentre
def compute_photocentre(
        image_link:str='../r3dresults/image.out',
        beam_size:float=2.5
    ):
    """
    Loads one r3d-image and computes "photo-centre"-coordinates
    with a radius of beam_size. Everything is in units AU.

    ARGUMENTS
      image_link:str - Link to image.out-file
      beam_size:float - radius of circle (in AU) within which
                        photocentre is found. Set to 0 to ignore!
    RETURNS
      posX  x-coordinates in AU
      posY  y-coordinates in AU
      posR  Radial distances in AU
    """

    # Load image
    image_folder = '/'.join(image_link.split('/')[:-1])+'/'
    image_file = image_link.split('/')[-1]
    image2d,image2dlog,totalflux,axisplot = load_images(
        path=image_folder,
        image=image_file,
        distance=1
    )
    # Compute photocentreposition
    # Uses Eqs 6 and 7 from Beguin+2024
    (Nx,Ny) = np.shape(image2d)
    xrange = np.linspace(axisplot[0],axisplot[1],Nx)
    yrange = np.linspace(axisplot[2],axisplot[3],Ny)

    nominatorX = 0
    nominatorY = 0
    demoninator = 0

    for nx in range(Nx):
        for ny in range(Ny):
            if beam_size != 0:
                if np.sqrt(xrange[nx]**2 + yrange[ny]**2) < beam_size:
                    # Xcoord-sums
                    nominatorX += image2d[ny,nx] * xrange[nx]
                    # Ycoord-sums
                    nominatorY += image2d[ny,nx] * yrange[ny]
                    # Common demoninator
                    demoninator += image2d[ny,nx]
            else:
                # Xcoord-sums
                nominatorX += image2d[ny,nx] * xrange[nx]
                # Ycoord-sums
                nominatorY += image2d[ny,nx] * yrange[ny]
                # Common demoninator
                demoninator += image2d[ny,nx]
            if beam_size > np.max(xrange):
                print('  WARNING: beam size is larger than image.')

    posX = nominatorX/demoninator
    posY = nominatorY/demoninator
    posR = np.sqrt(posX*posX + posY*posY)

    return posX,posY,posR


# Compute and write time & angle dependent photocentre positions
def write_photocentre_files(
        model_path:str='../r3dresults/',
        wavelength:str='02',
        angles:list=['i000_phi000'],
        beam_size:float=2.5
    ):
    """
    Computes and saves in a file, photocentre positions for many snapshots and angles
    of images.    

    ARGUMENTS
      model_path:str - Path to folder with all snapshot-subfolders
      wavelength:str - Wavelength of images written as eg '01' or '10' for 1 or 10um
      angles:list of str - List of strings with angles, eg 
                           ['i000_phi000','i180_phi000']
      beam_size:float - Radius of circle within which photocentre is computed
                        set to 0 to ignore (take whole image)

    RETURNS
      File: photocentre_{wavelength}um.dat in model_path-folder
    """
    print('Running: write_photocentre_files')
    # Automatically add / to end of path if it's missing
    if model_path[-1] != '/':
        model_path += '/'

    # Load all snapshots
    snapshots = [int(filename) for filename in os.listdir(model_path) if os.path.isdir(model_path+filename)]
    snapshots.sort()
    Nsnapshots = len(snapshots)


    # Open file to write all coords
    with open(f'{model_path}photocentre_{wavelength}um.dat','w') as fphotoc:
        print('  Writing photocentre file')

        # Write header
        fphotoc.writelines(f'# Photocentre coordinates at lambda = {wavelength} um\n')
        fphotoc.writelines(f'# Within circle with radius {beam_size} au (=0 -> whole image included).\n')
        fphotoc.writelines(f'Nsnapshots={Nsnapshots}\n')
        fphotoc.writelines(f'# Angles : \n')
        for angle in angles:
            fphotoc.writelines(f'    {angle}')
        fphotoc.writelines('\n# Snapshot      X-Y-R for each angle in AU')

        # Loop over time
        for snapshot in snapshots:

            # Print snapshot number
            fphotoc.writelines(f'\n{snapshot}    ')

            # Loop over angles
            for angle in angles:

                # Extract each image's photocentre
                pcX,pcY,pcR = compute_photocentre(
                    image_link=f'{model_path}{snapshot}/image_{angle}_{wavelength}um.out',
                    beam_size=beam_size
                )
                # Write X-Y-R- coords for each angle, in AU
                fphotoc.writelines(f'{pcX:.5f}  {pcY:.5f}  {pcR:.5f}    ')

    print('Done!')
