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

# ------------------------------------------------------------ #
# List of functions
#
#
# Functions that load various r3d input data
# ------------------------------------------
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
# plot_temperaturebins_raduis(
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
#
# ------------------------------------------------------------ #
# Functions that load various r3d input data

# Load grid properties
def load_grid_properties(
        amrpath:str='../amr_grid.inp'
    ):
    """
    Loads basic proparties from amr_grid.inp

    INPUT
    -----
    amrpath: str, path to amr_grid.inp

    OUTPUT
    ------
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
        print(f'ERROR: load_griddistances can not find {gridpath}.')


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
    Load and extracts densities of one dust specie of dust_density.inp

    INPUT
    -----
    path: string with path and filename of density file
    numb_specie: integer with number of the specie you want to load (default=1)
    
    OUTPUT
    ------
    Ncells: Number of cells in grid (nleafs)
    Nspec: Number of species in dust_density file
    dust_densities: np.array containing densities of all cells for specie numb_specie

    """

    if numb_specie == 0:
        print(f'ERROR numb of specie cant be 0')
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

        # Check that the chosen number of specie exists
        if numb_specie > Nspecies:
            print('\nERROR\nThis dust specie doesnt exist.\n')

        else:
            # Reduce specie number by 1 (index starts at 0)
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
      nspecie: number of specie you want to load
    
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
    
    # Check that specie number is correct
    if numb_specie <= Nspecies and numb_specie > 0:
        # Create density np.array
        dust_temperatures = np.zeros(Ncells)

        # Extract dust densities (reduce specie number by 1, since python starts at index=0)
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
        path:str='../r3dsims/spectrum.out',
        distance:float=1
    ):
    """
    Loads and returns SED and wavelength.

    ARGUMENTS
      path: str with path and filename of spectrum.out file
      distance: float, distance to star in pc

    RETURNS
      wavelengths: array with wavelength grid in microns
      spectrum: array with SED in Jy normalised to indicated distance
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
                npixels = int(line.split()[0])
            
            # row 3: pixel size in cm, divide by AUcm for AU
            if nl == 3:
                pixelsize_au = float(line.split()[0])/AUcm
            
            # NOTE might be useful later, commented out for now
            # row 4: wavelenght in um
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
    axisplot  = [0.5*size_au,-0.5*size_au,-0.5*size_au,0.5*size_au]

    # Total flux density of the image in Jy
    totalflux = sum(image1d) * 1.e23 * 2.35044305391e-11 * pixelsize_mas**2

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


    return image2d,image2dlog,totalflux,axisplot


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

    Inputs
    ------
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

    INPUT
    density_path: path to density.inp-file
    grid_path: path to grid_distances.csv'
    amr_path: path to amr_grid.inp
    numb_specie: number of the dust specie you want to plot

    OUTPUT
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
            title=f'Dust specie {nn+1}'
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
            title=f'Dust specie {nn+1}'
        )
    fig2.tight_layout()
    return fig1, fig2


# Plots temperature of radial bins (for paper-usage)
def plot_temperaturebins_radius(
        temperature_path:str='../dust_temperature.dat',
        grid_path:str='../grid_distances.csv',
        amr_path:str='../amr_grid.inp',
        numb_specie:int=1
    ):
    """
    Plots average temperature of 100 spherical shells, 
    and max-min-values, and STD of each shell.
    Can take a few minutes when loading the larger data sets.

    ARGUMENTS
      temperature_path:str='../dust_temperature.dat',
      grid_path:str='../grid_distances.csv',
      amr_path:str='../amr_grid.inp',
      numb_specie:int=1

    RETURNS
      fig,ax, temperature_bins,temperature_std,minmax_bins,radial_range
        figure-object
        axes-object
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
    Mstar,Rstar,Lstar = a5d.load_star_information(
        savpath='../co5bold_data/dst28gm06n052/st28gm06n052_186.sav',
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

    # Create figure-ax-objects
    fig, ax = plt.figure('Average radial temperatures', figsize=(6, 4)), plt.axes()
    # ax.plot
    # ax.set

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


    # For output: 
    #   create a min-max-array with average
    #   difference between min-max
    minmax_bins = 0.5*(temperature_max - temperature_min)


    return fig,ax, temperature_bins,temperature_std,minmax_bins,radial_range


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
    Load and plots absorption/scattering/scattering angles of all species in dustopac at folder of path
    
    ARGUMENTS
      path: path to folder containing dustkappa and dustopac-files

    RETURNS
      Figure and Axes-objects, plot with fig.show() 
    """

    # Automatically add / to end of path if it's missing
    if path[-1] != '/':
        path += '/'

    # Extract specie names from dustopac-file
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
    
    # Change specie name strings:
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
    #fig.set_figheight(3)
    #fig.set_figwidth(7)
    fig.tight_layout()

    return fig, ax, fluxtotal


# Plot a list of images in subplots (vertically)
def plot_imagesubplots(
        imagelist:list=['../image.out'],
        distance:float=1,
        scale:str='lin'
    ):
    """
    Plots a list of images from R3d in vertical subplots

    ARGUMENTS
      imagelist: list of paths to .out files
      distance: distance to source in pc
      scale: choce 'lin' or 'log' for lineare or logarithmic

    RETURNS
      gives matplotlib objects: fig and ax, plot with fig.show()
    """

    # Number of plots
    Nplots = len(imagelist)

    # Set fig-and-axis settings for subplots
    fig, ax = plt.subplots(
        Nplots,1, 
        dpi = 150, 
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

        # Only use lin or log depending on choice
        if scale == 'lin':
            # Save linear data in list (list of arrays)
            imagedata = image2d
        if scale == 'log':
            imagedata = image2dlog
        
        # Plot image at spot nn, set title and axis labels
        im0 = ax[nn].imshow(
            imagedata, 
            origin='lower', extent=axisplot, 
            cmap=plt.get_cmap('hot')
        )
        ax[nn].set(
            title=f'{modelname}_{phase}: i:{incl}, phi:{phi}, {wavelengthum} um', 
            xlabel='Offset (AU)',
            ylabel='Offset (AU)',
        )

        # Set colour bar settings and label
        divider = make_axes_locatable(ax[nn])
        cax = divider.append_axes('right', size='5%', pad=0.05)
        cb0 = plt.colorbar(im0, cax=cax, orientation = 'vertical',shrink=0.6,pad=0.15)
        cb0.set_label(label = rf'Flux at {distance} pc (Jy/asec$^2$)',fontsize= 10)

    # Change figure size
    fig.set_figheight(10)
    fig.set_figwidth(10)
    fig.tight_layout()

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

        # Create an Ncells long array with the opacities of each specie at correct cell index
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
    Nradius = 100
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
    dtauarray[-1] = dtauarray[-2]
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


# Function to plot Pearson-Chi2-comparison between two sets of data
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

    # Check so that length of both are correct
    if len(simulation) != len(observation):
        return 'ERROR: your input data arrays are not of equal length'

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

    ARGUMENTS
      path = path to spectrum-file including file name
      distance = distance the SED is observed from

    RETURNS
      luminosity in Watt
    """

    # Load spectrum
    wavelengths,spectrum = load_spectrum(path,1)
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
    INPUT
    wavelengthum: List of wavelengths in um
    spectrum: List of flux densities in Jy
    distance: distance to source in pc

    OUTPUT
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
        print('ERROR, wavelengths and spectrum have different lengths')

