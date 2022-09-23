# functions and tools for loading and analyzing output
# from co5bold.

# Import various libraries
import cython
import numpy as np
from scipy.io.idl import readsav
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm, xlabel, xscale, yscale
import os
import re

# My own libraries
import create_r3d_functions as c3d
import analyze_r3d_functions as a3d

# Basic definitions (AU as cython due to slow funcs below)
AUcm = cython.declare(cython.float ,1.49598e13) # cm
Msol = 1.989e33 # g
Rsol = 6.955e10 # cm
Lsol = 3.828e26 # W

# Note
# Rstar = 1.651AU (355 Rsun)
#
# Cython info: might be needed later
# @cython.cfunc
# this decorator when declaring what's included in the functions
# @cython.locals(a=cython.int)
# AUcm = cython.declare(cython.float ,1.49598e13) # cm

# ------------------------------------------------------------ #
# List of functions
# 
#
# Load C5D-data and saves in arrays
# ---------------------------------
#
# load_grid_properties(
#    savpath:str='../co5bold_data/dst28gm06n056/st28gm06n056_140.sav'
# )
#
# function that outputs star's props, but instead lum, mass, temperature, radius
# load_star_information(
#    savpath:str='../co5bold_data/dst28gm06n056/st28gm06n056_140.sav',
# )
#
# load_star_properties(
#    savpath:str='../co5bold_data/dst28gm06n056/st28gm06n056_140.sav'
# )
#
# load_dustspecies_names(
#    savpath:str='../co5bold_data/dst28gm06n056/st28gm06n052_186.sav'
# )
#
# load_dust_densitytemperature()
#    savpath:str='../co5bold_data/dst28gm06n052/st28gm06n052_186.sav',
#    nspecies:int=0
# )
#
#
# Plot c5d-data in r3d-grid
# -------------------------
#
# plot_opakapparadius(
#    path:str='../'
# )
#
# plot_opakappatemperature(
#    path:str='../'
# )
#
# plot_densitytemperature(
#    Teff:float=2700,
#    path:str='../r3dresults/st28gm06n056/140/',
#    nspecie=1
# )
#
#
# Functions to create r3d-data from c5d-data
# ------------------------------------------
#
# create_star(
#    savpath:str='../co5bold_data/dst28gm06n052/st28gm06n052_186.sav',
#    amrpath:str='../amr_grid.inp',
#    gridpath:str='../grid_distances.csv',
#    sizepath:str='../grid_cellsizes.csv',
# )
#
# NOTE New function, combines star's density and Rosseland Opacity into R3D-density
#      Star's Opacity is then 1 for all Wavelengths
# TODO Change the opacity later to some normalized function that varies realistically
#      to reproduce the star's SED more than as a BB/Grey body
# create_staropadensity(
#    pathopacity:str='../star_opacities.dat',
#    pathstardensity:str='../dust_density_onestarstar.inp',
#    pathwavelength:str='../wavelength_micron.inp',
#    corrfac:float=1.0
# )
#
#
# Smoothing of stellar data
# -------------------------
#
# smooth_stellardata(
#    path = '../r3dresults/st28gm06n056/',
#    phases = [140,141,142]
# )
#
# smooth_opacity(
#    path:str='../star_opacities.dat',
#    smooth_out:int = 4,
#    smooth_in:int = 0,
#    smooth_tolerance_log:float = 0.1
# )
#
# smooth_temperature(
#    path:str = '../dust_temperature.dat',
#    smooth_out:int = 10,
#    smooth_in:int = 3,
#    smooth_tolerance:float=1.5
# )
#
# smooth_density(
#    path:str = '../dust_density.inp',
#    smooth_out:int = 9,
#    smooth_in:int = 3,
#    smooth_tolerance:float=1.5
# )
#
#
# Funcs to load and create dusty envelope
# ---------------------------------------
#
# create_dustfiles(
#    savpath:str='../co5bold_data/dst28gm06n052/st28gm06n052_186.sav',
#    amrpath:str='../amr_grid.inp',
#    gridpath:str='../grid_distances.csv',
#    sizepath:str='../grid_cellsizes.csv',
#    Nspecies:int=1,
#    monomermasses:list=[2.3362e-22]
# )
#
# TODO
# extract_grainsizes()
#   borrow stuff from create_star and create_dustfiles, or incorporate in create dustfiles?
#
# ============================================================
# Functions that load C5D-data and saves them in arrays

# Load c5d grid properties
@cython.cfunc
def load_grid_properties(
        savpath:str='../co5bold_data/dst28gm06n056/st28gm06n056_140.sav'
    ):
    """
    TODO INFO
    outputorder: c5dgrid,cellcourners,cellsize
    """
    
    # Load sav-file
    c5ddata = readsav(savpath)

    # Extract data - TODO MAY HAVE TO CHOSE THIS MANUALLY DEPENDING ON FILE
    c5ddata = c5ddata['ful']

    # Get number of gridcells from co5bold data
    nc5dedge = cython.declare(cython.int, np.size(c5ddata['Z'][0][0][16]))

    # Declare output arrays
    c5dgrid = np.zeros((nc5dedge,3))
    cellcourners = np.zeros((nc5dedge,3))
    cellsizesx = []
    cellsizesy = []
    cellsizesz = []

    # Declare variables
    nn = cython.declare(cython.int)

    for nn in range(nc5dedge):

        # Save x,y,z distances in cm
        c5dgrid[nn,0] = c5ddata['Z'][0][0][16][0][0][nn]
        c5dgrid[nn,1] = c5ddata['Z'][0][0][19][0][nn][0]
        c5dgrid[nn,2] = c5ddata['Z'][0][0][22][nn][0][0]

        # Save coordinates of cell courners
        cellcourners[nn,0] = c5ddata['Z'][0][0][25][0][0][nn]
        cellcourners[nn,1] = c5ddata['Z'][0][0][28][0][nn][0]
        cellcourners[nn,2] = c5ddata['Z'][0][0][31][nn][0][0]


        # Extract cell sizes
        if nn > 0:
            cellsizesx.append(
                (cellcourners[nn,0] - cellcourners[nn-1,0])
            )
            cellsizesy.append(
                (cellcourners[nn,1] - cellcourners[nn-1,1])
            )
            cellsizesz.append(
                (cellcourners[nn,2] - cellcourners[nn-1,2])
            )

    # Add final grid courner
    cellcourners[-1,0] = c5ddata['Z'][0][0][25][0][0][-1]
    cellcourners[-1,1] = c5ddata['Z'][0][0][28][0][-1][0]
    cellcourners[-1,2] = c5ddata['Z'][0][0][31][-1][0][0]

    cellsizesx.append(
        (cellcourners[-1,0] - cellcourners[-2,0])
    )
    cellsizesy.append(
        (cellcourners[-1,1] - cellcourners[-2,1])
    )
    cellsizesz.append(
        (cellcourners[-1,2] - cellcourners[-2,2])
    )

    # Extract minimum grid size
    cellsize = (min(cellsizesx) + min(cellsizesy) + min(cellsizesz))/3

    return c5dgrid,cellcourners,cellsize


# function that outputs star's props, but instead lum, mass, temperature, radius
def load_star_information(
        savpath:str='../co5bold_data/dst28gm06n056/st28gm06n056_140.sav',
        printoutput:str='y'
    ):
    """
    Returns and prints a string with the basic properties of the star in the c5d-sim.

    ARGUMENT
      savpath:str = path to sav-file you want to work with

    RETURNS
      Mstar in gram
      Rstar in cm
      Lstar in Watt
    """

    c5ddata = readsav(savpath)
    c5ddata = c5ddata['ful']
    c5ddata = str(c5ddata['PAR'][0][0][2][1])[1:]

    if printoutput == 'y':
        print(f'Stellar info: {c5ddata}')

    # Extract all numbers in string
    numbers = re.findall("[0-9]+", c5ddata)
    
    # Save and return mass, radius and luminosity in gram, cm, and Watt
    Mstar = float(numbers[1]) * Msol
    Rstar = float(numbers[2]) * Rsol
    Lstar = float(numbers[3]) * Lsol

    return Mstar,Rstar,Lstar


# Extract co5bold densities into a separate array 
# - note this is probably faster than loading it in the r3d-file-writing functions
@cython.cfunc
def load_star_properties(
        savpath:str='../co5bold_data/dst28gm06n056/st28gm06n056_140.sav'
    ):
    """
    Loads c5d-star's densities, temperatures, and opacities, and puts them into 3D-arrays

    INPUT
    savpath: path to sav-file
    nbins: number of bins to put the opacities in, ie number of duststar-species for the star

    OUTPUT
    c5dstar_densities: 3D array with star's densities
    c5dstar_temperatures: 3D array with star's temperatures 
    c5dstar_opacities: 3D array with star's opacities
    """
    
    # Load sav-file
    c5ddata = readsav(savpath)

    # Extract data
    c5ddata = c5ddata['ful']

    # Get number of gridcells from co5bold data
    nc5dedge = cython.declare(cython.int, np.size(c5ddata['Z'][0][0][16]))

    # Declare variables
    nx = cython.declare(cython.int)
    ny = cython.declare(cython.int)
    nz = cython.declare(cython.int)

    # Declare np.array
    c5dstar_densities = np.zeros((nc5dedge,nc5dedge,nc5dedge))
    c5dstar_temperatures = np.zeros((nc5dedge,nc5dedge,nc5dedge))
    c5dstar_opacities = np.zeros((nc5dedge,nc5dedge,nc5dedge))

    # Extract densities - This can take time, some 2min per property
    for nx in range(nc5dedge):
        for ny in range(nc5dedge):
            for nz in range(nc5dedge):
                c5dstar_densities[nx,ny,nz] = c5ddata['Z'][0][0][34][nx][ny][nz]
                c5dstar_temperatures[nx,ny,nz] = c5ddata['EOS'][0][0][1][nx][ny][nz]
                c5dstar_opacities[nx,ny,nz] = c5ddata['OPA'][0][0][0][nx][ny][nz]
    
    return c5dstar_densities, c5dstar_temperatures, c5dstar_opacities



# Function that just lists the number and names of the dust species available in the data
def load_dustspecies_names(
        savpath:str='../co5bold_data/dst28gm06n052/st28gm06n052_186.sav'
    ):
    """
    Extract number of dust species and corresponding list of specie-names from c5d-data

    INPUT
    savpath:str = path to c5d-sav-file

    OUTPUT
    number of dust species
    list of strings with dust specie names
    """
    # Extract number of dust species in data
    c5ddata = readsav(savpath)
    c5ddata = c5ddata['ful']
    Nc5dspecies = int((len(c5ddata['Z'][0][0]) - 40)/3)

    # Loop over number of species and extract specie names also
    speciesnames = []
    for nspecies in range(Nc5dspecies):
        speciesnames.append(str(c5ddata['Z'][0][0][42+3*nspecies])[4:-1])

    # Return number of species and list of specie names
    return Nc5dspecies, speciesnames



# Function for loading one dust specie from c5d-data
@cython.cfunc
@cython.locals(nspecies=cython.int)
def load_dust_densitytemperature(
        savpath:str = '../co5bold_data/dst28gm06n052/st28gm06n052_186.sav',
        nspecies:int = 0
    ):
    """
    Loads c5d-data and extracts number density of dust monomers and dust temperature

    ARGUMENTS
      savpath:str = path to sav-file
      nspecies:int = number of the specie to extract, start with 0!

    RETURNS
      c5ddust_densities: array with monomer number density in the c5d-grid
      c5ddust_temperatures: array with dust temperatures within c5d-grid
    """

    # Load sav-file
    c5ddata = readsav(savpath)

    # Extract data
    c5ddata = c5ddata['ful']

    # Get number of gridcells from co5bold data
    nc5dedge = cython.declare(cython.int, np.size(c5ddata['Z'][0][0][16]))

    # Declare variables
    nx = cython.declare(cython.int)
    ny = cython.declare(cython.int)
    nz = cython.declare(cython.int)

    # Declare np.arrays for number density of dust monomers and temperatures
    c5ddust_densities = np.zeros((nc5dedge,nc5dedge,nc5dedge))
    c5ddust_temperatures = np.zeros((nc5dedge,nc5dedge,nc5dedge))

    # Extract densities and temperatures (this is time-demanding!)
    for nx in range(nc5dedge):
        for ny in range(nc5dedge):
            for nz in range(nc5dedge):

                # Densities:
                c5ddust_densities[nx,ny,nz] = c5ddata['Z'][0][0][40+3*nspecies][nx][ny][nz]

                # Temperatures (only save those cells where there is dust!)
                if c5ddust_densities[nx,ny,nz] > 0:
                    c5ddust_temperatures[nx,ny,nz] = c5ddata['EOS'][0][0][1][nx][ny][nz]
                else: 
                    c5ddust_temperatures[nx,ny,nz] = 0

    # Return density-temperature arrays
    return c5ddust_densities, c5ddust_temperatures


# Function that loads and extracts gas densities and dust densities of a chosen dust specie
# Primarily to use for getting grain sizes
@cython.cfunc
@cython.locals(nspecies=cython.int)
def load_dustgas_densities(
        savpath:str = '../co5bold_data/dst28gm06n052/st28gm06n052_186.sav',
        nspecies:int = 0
    ):
    """
    TODO
    fill here ...

    ARGUMENTS
      ...
    RETURNS
      ...    
    """

    # Load sav-file
    c5ddata = readsav(savpath)

    # Extract data
    c5ddata = c5ddata['ful']

    # Get number of gridcells from co5bold data
    nc5dedge = cython.declare(cython.int, np.size(c5ddata['Z'][0][0][16]))

    # Declare variables
    nx = cython.declare(cython.int)
    ny = cython.declare(cython.int)
    nz = cython.declare(cython.int)
    monomer_density = cython.declare(cython.float)
    gas_density = cython.declare(cython.float)

    # Declare np.arrays for number density of dust monomers and temperatures
    c5ddust_densities = np.zeros((nc5dedge,nc5dedge,nc5dedge))
    c5dstar_densities = np.zeros((nc5dedge,nc5dedge,nc5dedge))

    # Extract densities and temperatures (this is time-demanding!)
    for nx in range(nc5dedge):
        for ny in range(nc5dedge):
            for nz in range(nc5dedge):

                monomer_density = c5ddata['Z'][0][0][40+3*nspecies][nx][ny][nz]

                # check already before if data is non-zero?
                if monomer_density > 0:

                    # Dust monomer densities:
                    c5ddust_densities[nx,ny,nz] = monomer_density

                    # Gas density
                    gas_density = c5ddata['Z'][0][0][34][nx][ny][nz]
                    c5dstar_densities[nx,ny,nz] = gas_density
    
    return c5dstar_densities, c5ddust_densities



# ==========================================================================
# Plot c5d-data (but in r3d-grid)

# Plot c5d opacities (in r3d-grid)
def plot_opakapparadius(
        path:str='../'
    ):
    """
    Plots star_opacity.dat
    I.e. the kappa per r3dgridcell as given by c5d.
    Input: 
    path to folder of your r3d-model where amr_grid.inp, grid_distcances.csv, 
    and star_opacities.dat are included.
    """

    # Automatically add / to end of path if it's missing
    if path[-1] != '/':
        path += '/'

    # load opacity.dat
    opacity = c3d.load_staropacities(
        path = path+'star_opacities.dat'
    )

    # load grid
    griddistances = a3d.load_griddistances(
        amrpath=path+'amr_grid.inp',
        gridpath=path+'../grid_distances.csv'
    )
    radiusau = griddistances[:,0]/AUcm

    # Load and plots r3d density data for ONE dust specie
    fig, ax = plt.subplots(2,1)

    
    ax[0].plot(
        radiusau,opacity,
        linestyle='',marker='.',markersize=1
    )
    ax[0].set(
        ylabel=r'$\kappa_{\rm abs}$ (cm$^2$/g)',
        xlabel=r'Distance (AU)',
        title='Grid cell OPA'
    )

    ax[1].plot(
        radiusau,opacity,
        linestyle='',marker='.',markersize=1
    )
    ax[1].set(
        ylabel=r'$\kappa_{\rm abs}$ (cm$^2$/g)',
        xlabel=r'Distance (AU)',
        yscale='log'
    )

    fig.show()


# Plot c5d-OPA vs r3d-temperature
def plot_opakappatemperature(
        path:str='../'
    ):

    # Automatically add / to end of path if it's missing
    if path[-1] != '/':
        path += '/'

    # Load opacity.dat
    opacity = []
    with open(path+'star_opacities.dat', 'r') as fopacity:
        for line in fopacity.readlines():
            if line[0] != '#':
                opacity.append(float(line))
    opacity = np.array(opacity)


    # Load all opacities as given for R3D

    # First extract specie names from dustopac-file
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

    for specie_name in specie_names:
        specie_name,wavelengths,kappadata = a3d.load_onekappa(
            specie_name=specie_name,
            path=path
        )

        # Save kappa from each specie
        kappas.append(np.mean(np.array(kappadata[0])))


    # Load dust_temperature.dat
    temperature_path = path+'dust_temperature.dat'
    # Load first dust_temperature
    Ncells,Nspec,temperature = a3d.load_temperature(
        path=temperature_path,
        numb_specie=1
    )
    temperatures = [temperature]

    # Load rest of temperatures
    if Nspec > 1:
        for numb_specie in range(2,Nspec+1):
            temperatures.append(
                a3d.load_temperature(
                    path=temperature_path,
                    numb_specie=numb_specie
                )[2]
            )
    
    # Control colours of each density distribution
    colour = cm.rainbow(np.linspace(0, 1, Nspec))

    # Set objects for plot with all in the same figure
    fig, ax = plt.figure(), plt.axes()
    
    for nn, c in enumerate(colour):

        temperature = np.array(temperatures[nn])

        ax.plot(
            temperature[np.where(temperature > 0)[0]],
            opacity[np.where(temperature > 0)[0]],
            markeredgecolor=c,linestyle='',marker='.',markersize=1
        )

        # Plot kappa of each bin
        ax.plot(
            [
                np.min(temperature[np.where(temperature > 0)[0]]),
                np.max(temperature[np.where(temperature > 0)[0]])
            ],
            [kappas[nn],kappas[nn]],
            'k-'
        )

    ax.set(
        xlabel=r'Cell temperature (K)',
        ylabel=r'$\kappa_{\rm abs}$ (cm$^2$/g)',
        yscale='log',
        title=f'Dust species 1 to {Nspec}'
    )
    fig.tight_layout()
    fig.show()

# Plot Density vs Temperature

# Plot c5d-OPA vs r3d-temperature
def plot_densitytemperature(
        modelname:str='st28gm06n056',
        phase:str='140',
        nspecie=1
    ):
    """
    Plots density vs temperature of a given model, given specie
    """

    # Set path to r3d-data
    path = f'../r3dresults/{modelname}/{phase}/'

    Mstar,Rstar,Lstar = load_star_information(
        savpath = f'../co5bold_data/d{modelname}/{modelname}_{phase}.sav'
    )

    # Get paths to necessary files
    density_path = path+'dust_density.inp'
    grid_path = path+'../grid_distances.csv'
    amr_path = path+'amr_grid.inp'
    

    # load dust_density of specie nspecie
    Ncells,Nspec,density = a3d.load_dustdensity(
        path=density_path,
        numb_specie=nspecie
    )

    # Load dust_temperature of nspecie
    temperature_path = path+'dust_temperature.dat'
    Ncells,Nspec,temperature = a3d.load_temperature(
        path=temperature_path,
        numb_specie=nspecie
    )

    # Load radial distance to middle for each cell
    griddistances = a3d.load_griddistances(
        gridpath= grid_path,
        amrpath= amr_path
    )[:,0]



    # Separate density and temperature to two lists, inside and outside the star
    # ie inside or outside Rstar
    density_insidestar = []
    temperature_insidestar = []

    density_outsidestar = []
    temperature_outsidestar = []

    for nn in range(Ncells):
        if griddistances[nn] < Rstar:
            density_insidestar.append(density[nn])
            temperature_insidestar.append(temperature[nn])
        
        else:
            density_outsidestar.append(density[nn])
            temperature_outsidestar.append(temperature[nn])


    # Plot density vs temperature
    fig, ax = plt.figure(), plt.axes()

    # NOTE
    # I've been testing different ways of getting transperent dots, doesn't make a
    # difference but might be nice to keep here
    # Outside star
    ax.plot(
        temperature_outsidestar,
        density_outsidestar,
        mec=(0,0,1,0.3),
        #markeredgecolor='b',
        linestyle='',marker='.',markersize=1#, alpha=0.7
    )

    # Plot red dots inside the star

    ax.plot(
        temperature_insidestar,
        density_insidestar,
        mec=(1,0,0,0.3),
        #markeredgecolor='r',
        linestyle='',marker='.',markersize=1#, alpha=0.7
    )

    # Set final settings
    ax.set(
        xlabel=r'Cell temperature (K)',
        ylabel=r'$\rho$ (g cm$^{-3}$)',
        yscale='log',
        title=f'Dust species {nspecie}'
    )
    fig.tight_layout()
    fig.show()





# ==========================================================================
# Functions to create r3d-data from c5d-data

# Extract data on star and create r3d-files from it
def create_star(
        savpath:str='../co5bold_data/dst28gm06n052/st28gm06n052_186.sav',
        amrpath:str='../amr_grid.inp',
        gridpath:str='../grid_distances.csv',
        sizepath:str='../grid_cellsizes.csv',
    ):
    """
    Extracts data from Co5bold sav-files and creates a dust_star for Radmc3d.

    INPUT
    -----
    amrpath: path to amr_grid.inp
    gridpath: path to grid_distances.csv
    sizepath: path to grid_cellsizes.csv
    savpath: path to co5bold sav-file

    OUTPUT
    ------
    R3D density file: dust_density_onestar.inp
    R3D temperature file: dust_temperature_onestar.dat
    Useful file with list of extracted opacities 'kappa': star_opacities.dat
    """

    # Extract phase-designation from savpath
    phase = savpath.split('_')[-1].split('.')[0]

    # Load R3D grid
    print('Loading R3D grid')
    nleafs = a3d.load_grid_properties(amrpath=amrpath)[2]
    r3ddistances = a3d.load_griddistances(amrpath=amrpath,gridpath=gridpath)
    r3dcellsizes = a3d.load_cellsizes(amrpath=amrpath,sizepath=sizepath)

    # Load C5D grid
    print('Loading C5D grid properties')
    c5dgrid, c5dcellcourners, c5dcellsize = load_grid_properties(savpath=savpath)

    # Compute distance to courner of c5d-grid
    c5d_gridcournerdist = np.sqrt(
        np.max(c5dcellcourners[:,0])**2 + \
        np.max(c5dcellcourners[:,1])**2 + \
        np.max(c5dcellcourners[:,2])**2
    )

    # Compute distance to nearest basecell
    r3d_nearestbasecell = np.min(np.abs(
        r3ddistances[np.where(r3dcellsizes == r3dcellsizes.max())[0],0]
    ))

    # Check if the R3D-grid's nearest basecell is inside the c5d-grid
    if r3d_nearestbasecell + r3dcellsizes.max() >= c5d_gridcournerdist:
        print('\nERROR')
        print('    R3D basecells are all outside the CO5BOLD-grid, this will not work')
        print('    when translating gas densities from CO5BOLD to R3D. Stopping.')

    # Check so that the smallest c5dcells are not larger than the r3d's smallest cells
    elif r3dcellsizes.min() <= c5dcellsize:
        print('\nERROR')
        print('    R3D grid resolution is higher than C5D grid, stopping')
        print('    No output is given. Change your R3D grid cells to something larger.\n')

    else:
        # Load C5D star densities, temperatures and opacities
        print('Loading C5D star properties (density, temperature, opacity)')
        c5dstar_densities,c5dstar_temperatures,c5dstar_opacities = load_star_properties(savpath=savpath)

        # Start working :)
        print('Translating C5D gas data to R3D data')

        # Declare stuff for the loops
        r3d_densities = 0
        r3d_temperatures = 0
        r3d_opacities = 0
        progbar = 0

        # Open r3d data files
        with open(f'../dust_density_onestar_{phase}.inp', 'w') as fdensity, \
             open(f'../dust_temperature_onestar_{phase}.dat', 'w') as ftemperature, \
             open(f'../star_opacities_{phase}.dat', 'w') as fopacity:

            # Write headers:
            # 1
            # nleafs
            # number dust species
            fdensity.write(f'1\n{int(nleafs)}\n1\n')
            ftemperature.write(f'1\n{int(nleafs)}\n1\n')
            fopacity.write('# List of c5d-opacities translated to r3d-spatial grid.\n# Use as input when separating one-specie-density_star-file into several species\n# and creating dust-star opacity files.\n')

            # Loop over r3d grid
            for nr3d in range(nleafs):

                # Extract size range for current r3dcell
                r3dxrange = [
                    r3ddistances[nr3d,1]-0.5*r3dcellsizes[nr3d],
                    r3ddistances[nr3d,1]+0.5*r3dcellsizes[nr3d]
                ]
                r3dyrange = [
                    r3ddistances[nr3d,2]-0.5*r3dcellsizes[nr3d],
                    r3ddistances[nr3d,2]+0.5*r3dcellsizes[nr3d]
                ]
                r3dzrange = [
                    r3ddistances[nr3d,3]-0.5*r3dcellsizes[nr3d],
                    r3ddistances[nr3d,3]+0.5*r3dcellsizes[nr3d]
                ]   

                # Extract indeces of all c5dcells within current r3dcell
                c5dxrange = np.argwhere(r3dxrange[0] <= c5dgrid[np.argwhere(c5dgrid[:,0] <= r3dxrange[1]),0])[:,0]
                c5dyrange = np.argwhere(r3dyrange[0] <= c5dgrid[np.argwhere(c5dgrid[:,1] <= r3dyrange[1]),1])[:,0]
                c5dzrange = np.argwhere(r3dzrange[0] <= c5dgrid[np.argwhere(c5dgrid[:,2] <= r3dzrange[1]),2])[:,0]

                # Number of c5dcells within r3dcell (which is faster?)
                #nchildcells = c5dxrange.size*c5dyrange.size*c5dzrange.size
                nchildcells = 0

                # Then loop through c5dcells within r3dcell
                for nnz in c5dzrange:
                    for nny in c5dyrange:
                        for nnx in c5dxrange:

                            # Save sum of all c5d-data for each r3d-cell
                            r3d_densities += c5dstar_densities[nnx,nny,nnz]
                            r3d_temperatures += c5dstar_temperatures[nnx,nny,nnz]
                            r3d_opacities += c5dstar_opacities[nnx,nny,nnz]

                            # Number of cells
                            nchildcells += 1

                # Check if there actually are any c5dcells within r3dcell
                # If not, then your r3dgrid is probably smaller than the c5dgrid
                # and then the density and temperature will be zero for some cells
                # This shouldn't happen since I have the if statement earlier that breaks
                # this function would that happen
                if nchildcells > 0:
                    # Otherwise save the average of the c5d-data of each r3dcell
                    r3d_densities /= nchildcells
                    r3d_temperatures /= nchildcells
                    r3d_opacities /= nchildcells

                else:
                    r3d_densities = c5dstar_densities[c5dxrange[0],c5dyrange[0],c5dzrange[0]]
                    r3d_temperatures = c5dstar_temperatures[c5dxrange[0],c5dyrange[0],c5dzrange[0]]
                    r3d_opacities = c5dstar_opacities[c5dxrange[0],c5dyrange[0],c5dzrange[0]]
                    print(f'nchildcells = {nchildcells}')

                # Then write data to r3d files
                fdensity.write(f'{r3d_densities}\n')
                ftemperature.write(f'{r3d_temperatures}\n')
                fopacity.write(f'{r3d_opacities}\n')

                # Reset data
                r3d_densities = 0
                r3d_temperatures = 0
                r3d_opacities = 0
                nchildcells = 0

                # Some progress bar info
                if int(nr3d/nleafs*100) == 25 and progbar == 0:
                    progbar += 1
                    print('Finished 25 per cent of the grid.')

                if int(nr3d/nleafs*100) == 50 and progbar == 1:
                    progbar += 1
                    print('Finished 50 per cent of the grid.')

                if int(nr3d/nleafs*100) == 75 and progbar == 2:
                    progbar += 1
                    print('Finished 75 per cent of the grid.')

    print(f'C5D Dust-star:\n    dust_density_onestar_{phase}.inp\n    dust_temperature_onestar_{phase}.dat\n    star_opacities_{phase}.dat\nDONE\n')


# For when creating several stars/phases
def create_stars(
        modelnames:list = ['st28gm06n056'],
        phases:list = [140,141,142]
    ):
    """
    Function for creating many star-models, from lists of them. 
    Also moves files to correct folders. Also useful when only doing one model and phase.
        NOTE this is not implemented correctly yet, only do one model (several phases are OK)!

    INPUT
    modelnames:list = list of modelnames
    phases:list = list of phases, listed in the same order as models are listed
        NOTE this is not implemented correctly yet, only do one model!
    """

    # TODO there should be a list of phases per model inputted in 

    for modelname in modelnames:
        for phase in phases:

            print(f'    Translating C5D star model {modelname} and phase {phase} to R3D')

            savpath = f'../co5bold_data/d{modelname}/{modelname}_{phase}.sav'
            path = f'../r3dresults/{modelname}/'

            create_star(
                savpath = savpath,
                amrpath = f'{path}amr_grid.inp',
                gridpath = f'{path}grid_distances.csv',
                sizepath = f'{path}grid_cellsizes.csv'
            )

            # Move files to specific folders
            print('Moving files\n')

            # Files from create star
            os.system(
                f'mv ../dust_density_onestar_{phase}.inp {path}{phase}/dust_density_onestar.inp'
            )
            os.system(
                f'mv ../dust_temperature_onestar_{phase}.dat {path}{phase}/dust_temperature_onestar.dat'
            )
            os.system(
                f'mv ../star_opacities_{phase}.dat {path}{phase}/star_opacities.dat'
            )

            # Copy files from create grid-part above to phase-folders
            os.system(
                f'cp {path}amr_grid.inp {path}{phase}/'
            )
            os.system(
                f'cp {path}wavelength_micron.inp {path}{phase}/'
            )


# Takes the stellar density data and opacity data from c5d
# creates a density_star-file where the "density" depends on each cell's mass density
# and the opacity-kappa of each cell as given by c5d.
# Creates ONE kappa_star that is only = 1 over all wavelengths.
def create_staropadensity(
        pathopacity:str='../star_opacities.dat',
        pathstardensity:str='../dust_density_onestarstar.inp',
        pathwavelength:str='../wavelength_micron.inp',
        phase:str=1,
        corrfac:float=1.0
    ):
    """
    INPUT
    pathopacity: path to star_opacities.dat',
    pathstardensity: path to dust_density_onestarstar.inp',
    pathwavelength: path to wavelength_micron.inp',

    OUTPUT
    dust_density_opastar_{phase}.inp
    dustopac_star_{phase}.inp
    dustkappa_opastar_{phase}.inp
    """

    print('Loading density, opacity, wavelengths')

    # Load star densities (in r3d-grid)
    Ncells,Nspec,star_densities = a3d.load_dustdensity(path=pathstardensity,numb_specie=1)

    # load star opacities (in r3d-grid)
    opacity = []
    with open(pathopacity, 'r') as fopacity:
        for line in fopacity.readlines():
            if line[0] != '#':
                opacity.append(float(line))
    opacity = np.array(opacity)

    # Load wavelengthgrid
    wavelengths,Nwave = a3d.load_wavelengthgrid(path=pathwavelength)
    wavelengths = np.array(wavelengths)

    # Create gas opacity
    #kappa = corrfact * np.median(wavelengths) / wavelengths

    # Assume
    # r3dopacity * r3ddensity = c5dopacity * c5ddensity
    #
    # Opacity file is set to 1 for all wavelengths
    # r3dopacity = 1
    # This means that
    # r3ddensity = c5dopacity * c5ddensity
    # for each cell
    #
    # Thus we obtain a simple grey body model for the opacity for all cells without having
    # to add hundreds and hundreds of opacity files.

    print(f'Change density to densityr3d = {corrfac} * kappac5d * densityc5d')
    #print(f'Change density to densityr3d = kappac5d^{kramer_exponent} * densityc5d')

    opacity_densities = np.zeros(Ncells)
    #meandensity = np.mean(star_densities)

    for nn in range(Ncells):
        opacity_densities[nn] = corrfac * opacity[nn] * star_densities[nn]

    # Print new star-opacity-density file
    print('Writing new radmc3d-files')

    with open(f'../dust_density_opastar_{phase}.inp', 'w') as fdensity:
        
        # Write headers
        fdensity.write(f'1\n{int(Ncells)}\n1\n')

        # Write densities and temperatures
        for nn,dens in enumerate(opacity_densities):
            fdensity.write(f'{dens}\n')

    # Print new dustopac_starbins.inp file
    print('Writing opacity files for the binned star.')
    with open(f'../dustopac_star_{phase}.inp', 'w') as fopac:

        # Print header
        fopac.write(f'2\n1\n-----------------------------\n')
        # Print specie name
        fopac.write(f'1\n0\nopastar\n-----------------------------\n')

    with open(f'../dustkappa_opastar_{phase}.inp', 'w') as fopac:

        # Write header (1 for no scattering in this and number of wavelengths)
        fopac.write(f'1\n{Nwave}\n')

        # Write wavelength, abscoeff, scattercoeff
        for nn,wave in enumerate(wavelengths):
            fopac.write(f'{wave}    1.0    0.0\n')
    
    print(f'C5D create star opacities densities:\n    dust_density_opastar_{phase}.inp\n    dustopac_star_{phase}.inp\n    dustkappa_opastar_{phase}.inp\nDONE\n')

# Smoothing of data-functions ----------------------------------------------------
#
# The c5d-data contains various cells in/near the surface of the star which gives
# very high luminosity and temperature. I smooth the r3d-data files to remove this.
# 
# There are functions for smoothin density, opacity and temperatures but the main
# function only does this for the opacity and density since the temperature 
# smoothing gives sketchy results and I don't trust them. Too many artifacts with
# very little smoothing in the r3d-data, and not much changes in luminosity.


# Main stellar-data-smoothing function, smooths opacity and density from negative spikes
# using the optimal settings I wound after a few weeks of exploring :P
def smooth_stellardata(
        path:str='../r3dresults/st28gm06n056/',
        phases:list=[140,141,142],
        starradii:list=[1,1,1],
        griddistances:list=[0],
        clean_data:str='y'
    ):
    """
    Removes spikes in stellar gas-data's temperature, opacity and mass density. 
    These spikes results in strange and high luminosities. Compares cells with 
    median values of cells around the spike in the r3d-lists and changes the spike
    to the median according to certain specifics based on general analyzes of
    all phases of st28gm06n056.

    ARGUMENTS
      path:str = path to model-folder with r3d-data.
      phases:list = list of phases within model-folder where data needs smoothing
      starradii:list = radii of the star in each phase in centimeters
      griddistances:list = list of distances to all R3D grid cells in centimeters
      clean_data:str = Default 'y', will REMOVE all intermediate files!

    RETURNS
      dust_temperature_onestar_smoothed.dat
      star_opacities_smoothed.dat
      dust_density_opastar.inp
      dustkappa_opastar.inp
      dustopac_opastar.inp
    """

    # Automatically add / to end of path if it's missing
    if path[-1] != '/':
        path += '/'

    for nphase,phase in enumerate(phases):

        print(f'    Smoothing star for phase {phase}')

        # Stellar radius is then (in centimeters)
        starradius = starradii[nphase]

        # Smooth a5d-Temperatures, remove pos spikes
        smooth_temperature(
            path = f'{path}{phase}/dust_temperature_onestar.dat',
            starradius = starradius,
            phase = phase,
            griddistances = griddistances,
            smooth_out = 10,
            smooth_in = 3,
            smooth_tolerance = 1.0
        )
        os.system(f'mv ../dust_temperature_smoothed_{phase}.dat {path}{phase}/dust_temperature_onestar_smoothed.dat')
    

        # Smooth a5d-Opacity, remove neg spikes
        smooth_opacity(
            path = f'{path}{phase}/star_opacities.dat',
            starradius=starradius,
            phase = phase,
            griddistances=griddistances,
            smooth_out = 7,
            smooth_in = 6,
            smooth_tolerance_log = 0
        )
        os.system(f'mv ../star_opacities_smoothed_{phase}.dat {path}{phase}/star_opacities_smoothed.dat')

        # Smooth a5d-densities, remove negative spikes
        smooth_density(
            path = f'{path}{phase}/dust_density_onestar.inp',
            starradius=starradius,
            phase = phase,
            griddistances=griddistances,
            smooth_out = 37,
            smooth_in = 36,
            smooth_tolerance = 0.4
        )
        os.system(f'mv ../dust_density_smoothed_{phase}.inp {path}{phase}/dust_density_onestar_smoothed.inp')

        # Combine c5d-density and opacity into one and save in r3d-density file
        create_staropadensity(
            pathopacity = f'{path}{phase}/star_opacities_smoothed.dat',
            pathstardensity = f'{path}{phase}/dust_density_onestar_smoothed.inp',
            pathwavelength = f'{path}wavelength_micron.inp',
            phase = phase,
            corrfac = 1.0
        )
    
        # Clean intermediate files?
        # if yes, then remove all intermediate star-files
        if clean_data == 'y':
            os.system(f'rm {path}{phase}/dust_density_onestar_smoothed.inp')
            os.system(f'rm {path}{phase}/dust_density_onestar.inp')
            os.system(f'rm {path}{phase}/dust_temperature_onestar.dat')
            os.system(f'rm {path}{phase}/star_opacities.dat')

        # Move/copy all files to correct places
        os.system(f'mv ../dust_density_opastar_{phase}.inp {path}{phase}/dust_density_opastar.inp')
        os.system(f'mv ../dustkappa_opastar_{phase}.inp {path}{phase}/dustkappa_opastar.inp')
        os.system(f'mv ../dustopac_star_{phase}.inp {path}{phase}/dustopac_opastar.inp')

    print(f'smooth_stellardata ({phase}):\n    dust_temperature_onestar_smoothed.dat\n    star_opacities_smoothed.dat\n    dust_density_opastar.inp\n    dustkappa_opastar.inp\n    dustopac_opastar.inp\nDONE\n')


# Smooths the opacity file, removes spikes
def smooth_opacity(
        path:str='../star_opacities.dat',
        starradius:float=1,
        phase:str = 1,
        griddistances:list=[0],
        smooth_out:int = 4,
        smooth_in:int = 0,
        smooth_tolerance_log:float = 0.1
    ):
    """
    Remove outlier cells with large negative change in opacity within the star's surface

    INPUT
    Path: path to star_opacities.dat
    smooth_tolerance_log: Limit in order of number of negative orders of magnitude from median
    Higher smooth tolerance > more tolerant for spikes

    OUTPUT
    New file: star_opacities_smoothed_{phase}.dat
    """
    print('Removing opacity spikes')

    # Load opacity-file
    opacity = c3d.load_staropacities(path = path)

    # Declarations
    Ncells = opacity.size
    counter = 0

    # Loop over all cells (except the outermost cells)
    for nn in range(smooth_out,Ncells-smooth_out):

        # Check if nn is inside the star (plus a small tolerance factor)
        if griddistances[nn] < 1.01*starradius:

            ## Using mean opacity
            #mean_opacity = 0.25 * (opacity[nn-2] + opacity[nn-1] + opacity[nn+1] + opacity[nn+2])
            #if opacity[nn] < 10**-smooth_tolerance_log * mean_opacity \
            #or opacity[nn] > 10**smooth_tolerance_log  * mean_opacity:
            #    opacity[nn] = mean_opacity
            #    counter += 1
        
            # Median of a range of numbers
            # Extract indeces, all within range except current index
            nindeces = [
                nmedian for nmedian in range(nn-smooth_out,nn+smooth_out+1) if nmedian < (nn-smooth_in) or nmedian > (nn+smooth_in)
            ]
            median_list = []

            for nmedian in nindeces:
                median_list.append(opacity[nmedian])

            median_opacity = np.median(np.array(median_list))

            if opacity[nn] < 10**-smooth_tolerance_log * median_opacity :
                opacity[nn] = median_opacity
                counter += 1

    print(f'{smooth_in}-{smooth_out}: Number of smoothed OPAcells: {counter}')

    # Write new file
    with open(f'../star_opacities_smoothed_{phase}.dat', 'w') as fopacity:

        # Write header
        fopacity.write('# List of c5d-opacities translated to r3d-spatial grid.\n# Use as input when separating one-specie-density_star-file into several species\n# and creating dust-star opacity files.\n')

        # Write new data
        for nn in range(Ncells):
            fopacity.write(f'{opacity[nn]}\n')

    print(f'C5D smooth opacities:\n    star_opacities_smoothed_{phase}.dat\nDONE\n')


# Smoothing, removing spikes in temperatures
def smooth_temperature(
        path:str = '../dust_temperature.dat',
        starradius:float=1,
        phase:str = 1,
        griddistances:list=[0],
        smooth_out:int = 10,
        smooth_in:int = 3,
        smooth_tolerance:float=1.5
    ):
    """
    Remove outlier cells with high temperatures

    Higher smooth_tolerance > more tolerant for spikes

    TODO write more info
    
    """
    print('Removing temperature spikes')

    # load temperature
    Ncells,Nspecies,temperatures = a3d.load_temperature(path=path)

    counter = 0

    # Loop over eventual species
    for nspecie in range(Nspecies):
        
        # Loop over grid cells of each specie (except the outermost cells)
        for ncell in range(smooth_out,Ncells-smooth_out):

            # Check if nn is inside the star (plus a small tolerance factor)
            if griddistances[ncell] < 1.01*starradius:

                # Index of cell in total list
                nn = ncell + Ncells*nspecie

                nindeces = [
                    nmedian for nmedian in range(nn-smooth_out,nn+smooth_out+1) if nmedian < (nn-smooth_in) or nmedian > (nn+smooth_in)
                ]
                median_list = []

                for nmedian in nindeces:
                    median_list.append(temperatures[nmedian])

                median_temperature = np.median(np.array(median_list))


                if temperatures[nn] > smooth_tolerance * median_temperature:
                    temperatures[nn] = median_temperature
                    counter += 1

    print(f'{smooth_in}-{smooth_out}: Number of smoothed Tcells: {counter}')

    # Print new temperature file
    with open(f'../dust_temperature_smoothed_{phase}.dat', 'w') as ftemperature:
        # Write headers:
        # 1
        # nleafs
        # number dust species
        ftemperature.write(f'1\n{int(Ncells)}\n{Nspecies}\n')

        for nn in range(Ncells*Nspecies):
            ftemperature.write(f'{temperatures[nn]}\n')

    print(f'C5D smooth temperatures:\n    dust_temperature_smoothed_{phase}.dat\nDONE\n')


# Smooth density
def smooth_density(
        path:str = '../dust_density_onestar.inp',
        starradius:float=1,
        phase:str=1,
        griddistances:list=[0],
        smooth_out:int = 9,
        smooth_in:int = 3,
        smooth_tolerance:float=1.5
    ):
    """
    Remove outlier cells with low densities, ie remove negative spikes

    Only use _onestar-density file for this!

    Smaller smooth_tolerance > more tolerant for spikes (is this correct?)

    TODO write more&better info
    """
    print('Removing density spikes')

    # load density
    Ncells,Nspecies,densities = a3d.load_dustdensity(path=path,numb_specie=1)

    counter = 0
        
    # Loop over grid cells (except outermost cells)
    for nn in range(smooth_out,Ncells-smooth_out):

        # Check if nn is inside the star (plus a small tolerance factor)
        if griddistances[nn] < 1.01*starradius:

            nindeces = [
                nmedian for nmedian in range(nn-smooth_out,nn+smooth_out+1) if nmedian < (nn-smooth_in) or nmedian > (nn+smooth_in)
            ]
            median_list = []

            for nmedian in nindeces:
                median_list.append(densities[nmedian])

            median_densities = np.median(np.array(median_list))


            if densities[nn] < smooth_tolerance * median_densities:
                densities[nn] = median_densities
                counter += 1


    print(f'{smooth_in}-{smooth_out}: Number of smoothed Density cells: {counter}')

    # Print new temperature file
    with open(f'../dust_density_smoothed_{phase}.inp', 'w') as fdensity:
        # Write headers:
        # 1
        # nleafs
        # number dust species
        fdensity.write(f'1\n{int(Ncells)}\n{Nspecies}\n')

        for nn in range(Ncells*Nspecies):
            fdensity.write(f'{densities[nn]}\n')

    print(f'C5D smooth densities:\n    dust_density_smoothed_{phase}.inp\nDONE\n')


# ====================================================================
# Funcs to load and create dusty envelope

# Extract and construct dust_density-files from C5D-data
# uses ['Z'][0][0][40][x][y][z] and ['Z'][0][0][43][x][y][z]
# or rather use
# ['Z'][0][0][40+3*nspecie][x][y][z]
# where nspecie is 0, 1, 2 ... up until (number of dust species)-1
# and the data are in number densities cm^-3 of monomers
# so the mass density is the number density * mass of dust forming molecule

# Also:
# str(teststar['Z'][0][0][40+3*nspecie+2])[4:-1]
# gives a string with the name of the specie!
# so this can also print the dustkappa-list-file for r3d!

def create_dustfiles(
        savpath:str='../co5bold_data/dst28gm06n052/st28gm06n052_186.sav',
        amrpath:str='../amr_grid.inp',
        gridpath:str='../grid_distances.csv',
        sizepath:str='../grid_cellsizes.csv',
        Nspecies:int=1,
        monomermasses:list=[2.3362e-22]
    ):
    """
    Creates R3D-files dust_density, dust_temperature and dustopac-list from C5D dust envelope-data.

    ARGUMENTS
      savpath:str = path to c5d-sav-file
      amrpath:str = path to r3d-amr-grid-file
      gridpath:str = path to corresponding file with list of r3d radial grid-distances
      sizepath:str = path to corresponding file with list of r3d cell sizes
      Nspecies:int = number of dust species to include (this runs through them in order... Might change later)
      monomermasses:list = list of masses of dust species monomers (check which species are included in c5d with a5d.load_dustspecies_names)

    RETURNS
      dust_density_dust_{phase}.inp
      dust_temperature_dust_{phase}.dat
      dustopac_dust_{phase}.inp
    """

    # Extract phase-designation from savpath
    phase = savpath.split('_')[-1].split('.')[0]

    # Load R3D grid
    print('Loading R3D grid')
    nleafs = a3d.load_grid_properties(amrpath=amrpath)[2]
    r3ddistances = a3d.load_griddistances(amrpath=amrpath,gridpath=gridpath)
    r3dcellsizes = a3d.load_cellsizes(amrpath=amrpath,sizepath=sizepath)

    # Load C5D grid
    print('Loading C5D grid properties')
    c5dgrid, c5dcellcourners, c5dcellsize = load_grid_properties(savpath=savpath)

    # Load number of dust species and specie-names
    Nc5dspecies, specienames =  load_dustspecies_names(savpath=savpath)

    # Check so that the smallest c5dcells are not larger than the r3d's smallest cells
    if r3dcellsizes.min() <= c5dcellsize:
        print('\nERROR')
        print('    R3D grid resolution is higher than C5D grid, stopping')
        print('    No output is given. Change your R3D grid cells to something larger.\n')
    
    # Check so that the number of dust species declared is not larger than available
    elif Nspecies > Nc5dspecies:
        print('\nERROR')
        print('    You asked for more dust species than available in C5D-data')
        print(f'    Number of dust species in C5D-data: {Nc5dspecies}')

    else:
        print(f'Translating C5D dust data to R3D dust data ({phase})')

        # Open r3d data files
        with open(f'../dust_density_dust_{phase}.inp', 'w') as fdensity, \
             open(f'../dust_temperature_dust_{phase}.dat', 'w') as ftemperature, \
             open(f'../dustopac_dust_{phase}.inp', 'w') as fopac:

            # Write headers:
            #
            # Density:
            # 1
            # nleafs
            # number dust species
            fdensity.write(f'1\n{int(nleafs)}\n{int(Nspecies)}\n')

            # Temperature:
            # 1
            # nleafs
            # number dust species
            ftemperature.write(f'1\n{int(nleafs)}\n{int(Nspecies)}\n')

            # dustopac:
            # 2
            # Number of species
            # -----------------------
            fopac.write(f'2\n{int(Nspecies)}\n-----------------------------\n')

            # Loop through the number of species you want to include
            # and write densities and temperatures in the files.
            # Each specie is listed in same files but after eachother.
            for nspecies in range(Nspecies):

                # Declare stuff for the loops
                r3d_density = 0
                r3d_temperature = 0
                progbar = 0
                speciesname = specienames[nspecies]

                # Load c5d-dust densities and temperatures
                c5ddensities, c5dtemperatures = load_dust_densitytemperature(savpath=savpath)

                # Some output
                print(f'Writing dust specie number {nspecies+1}:')
                print(f'    {speciesname}')
                print(f'    Monomer mass: {monomermasses} g')

                # Write the dustopac file
                # 1
                # 0
                # speciename
                # ---------------
                fopac.write(f"1\n0\n{speciesname}\n-----------------------------\n")

                # Loop over the r3d-grid
                for nr3d in range(nleafs):

                    # Extract size range for current r3dcell
                    r3dxrange = [
                        r3ddistances[nr3d,1]-0.5*r3dcellsizes[nr3d],
                        r3ddistances[nr3d,1]+0.5*r3dcellsizes[nr3d]
                    ]
                    r3dyrange = [
                        r3ddistances[nr3d,2]-0.5*r3dcellsizes[nr3d],
                        r3ddistances[nr3d,2]+0.5*r3dcellsizes[nr3d]
                    ]
                    r3dzrange = [
                        r3ddistances[nr3d,3]-0.5*r3dcellsizes[nr3d],
                        r3ddistances[nr3d,3]+0.5*r3dcellsizes[nr3d]
                    ]   

                    # Extract indeces of all c5dcells within current r3dcell
                    c5dxrange = np.argwhere(r3dxrange[0] <= c5dgrid[np.argwhere(c5dgrid[:,0] <= r3dxrange[1]),0])[:,0]
                    c5dyrange = np.argwhere(r3dyrange[0] <= c5dgrid[np.argwhere(c5dgrid[:,1] <= r3dyrange[1]),1])[:,0]
                    c5dzrange = np.argwhere(r3dzrange[0] <= c5dgrid[np.argwhere(c5dgrid[:,2] <= r3dzrange[1]),2])[:,0]

                    # Number of c5dcells within r3dcell (with data)
                    ndustcells = 0
                    ntempcells = 0

                    # Then loop through c5dcells within r3dcell
                    for nnz in c5dzrange:
                        for nny in c5dyrange:
                            for nnx in c5dxrange:
                                # Sum all densities and temperatures (only those with data)

                                if c5ddensities[nnx,nny,nnz] > 0:
                                    r3d_density += c5ddensities[nnx,nny,nnz]
                                    ndustcells += 1

                                if c5dtemperatures[nnx,nny,nnz] > 0:
                                    r3d_temperature += c5dtemperatures[nnx,nny,nnz]
                                    ntempcells += 1

                    # Take average of those cells with data (and average with respect
                    # to number of cells containing data only)
                    if ndustcells > 0:
                        # Recalculate number density of monomers to mass density
                        # eg
                        # Mg2SiO4: 2*24.305u + 28.085u + 4*15.999u = 140.69u = 2.3362e-22 gram
                        r3d_density *= monomermasses[nspecies] / ndustcells

                    # Also for temperature cells
                    if ntempcells > 0:
                        r3d_temperature /= ntempcells

                    # Write data to r3d files
                    fdensity.write(f'{r3d_density}\n')
                    ftemperature.write(f'{r3d_temperature}\n')

                    # Reset data
                    r3d_density = 0
                    r3d_temperature = 0
                    ndustcells = 0
                    ntempcells = 0

                    # Some progress bar info
                    if int(nr3d/nleafs*100) == 25 and progbar == 0:
                        progbar += 1
                        print('Finished 25 per cent of the grid.')

                    if int(nr3d/nleafs*100) == 50 and progbar == 1:
                        progbar += 1
                        print('Finished 50 per cent of the grid.')

                    if int(nr3d/nleafs*100) == 75 and progbar == 2:
                        progbar += 1
                        print('Finished 75 per cent of the grid.')

    # End functions with aknowledgements
    print(f'C5D Dust-data:\n    dust_density_dust_{phase}.inp\n    dust_temperature_dust_{phase}.dat\n    dustopac_dust_{phase}.inp\nDONE\n')


# TODO
# Function that extracts and saves grain sizes per R3d-Cell
# will also bin these by a number of grain sizes
# but first i will have to rwite something that extracts max and min grain sizes, 
# see if this is a logarithmic range or not, such things
#
def extract_grainsizes(
        amrpath:str='../r3dresults/st28gm06n052/amr_grid.inp',
        gridpath:str='../r3dresults/st28gm06n052/grid_distances.csv',
        sizepath:str='../r3dresults/st28gm06n052/grid_cellsizes.csv',
        savpath:str='../co5bold_data/dst28gm06n052/st28gm06n052_186.sav',
        Amon:float=2.3362e-22,
        rhomon:float=3.27,
        nHnd:float=3e-16,
        mH:float=1.6726e-27,
        epsilonHe:float=0.1
    ):
    """
    Info about the grain sizes, what is the equation?
    TODO
    
    ARGUMENTS
      PATHS
        amrpath:str=
        gridpath:str=
        sizepath:str=
        savpath:str=
      CONSTANTS
        Amon = 2.3362e-22 # g
        rhomon = 3.27 # g cm-3
        nHnd = 3e-16
        mH = 1.6726e-27 # g
        epsilonHe = 0.1

    RETURNS
      file: grain_sizes.dat
    """
    # Extract phase-designation from savpath
    phase = savpath.split('_')[-1].split('.')[0]

    # Compute constants
    grainsize_constants = 3/(4*np.pi) * Amon/rhomon * nHnd * mH * (1+epsilonHe)

    # Load R3D grid
    print('Loading R3D grid')
    nleafs = a3d.load_grid_properties(amrpath=amrpath)[2]
    r3ddistances = a3d.load_griddistances(amrpath=amrpath,gridpath=gridpath)
    r3dcellsizes = a3d.load_cellsizes(amrpath=amrpath,sizepath=sizepath)

    # Load C5D grid
    print('Loading C5D grid properties')
    c5dgrid, c5dcellcourners, c5dcellsize = load_grid_properties(savpath=savpath)

    # Check so that the smallest c5dcells are not larger than the r3d's smallest cells
    if r3dcellsizes.min() <= c5dcellsize:
        print('\nERROR')
        print('    R3D grid resolution is higher than C5D grid, stopping')
        print('    No output is given. Change your R3D grid cells to something larger.\n')
    
    else:
        print(f'Computing grain sizes and saving them in R3D-grid ({phase})')

        # Open data files
        with open(f'../grain_sizes_{phase}.dat', 'w') as fsizes:

            # Write header
            fsizes.write(f'# List of grain sizes of each cell.\n# Same order as in R3D density and temperature files.\n# As extracted from {savpath}\n')

            # Declare stuff for the loops
            monomer_density = 0
            gas_densities = 0
            progbar = 0

            # Some output
            print(f'Extracting gas and monomer densities from CO5BOLD-files ({phase})')

            # Extract gas and dust monomer densities
            gas_densities, monomer_densities = load_dustgas_densities()

            # Loop over the r3d-grid
            for nr3d in range(nleafs):

                # Extract size range for current r3dcell
                r3dxrange = [
                    r3ddistances[nr3d,1]-0.5*r3dcellsizes[nr3d],
                    r3ddistances[nr3d,1]+0.5*r3dcellsizes[nr3d]
                ]
                r3dyrange = [
                    r3ddistances[nr3d,2]-0.5*r3dcellsizes[nr3d],
                    r3ddistances[nr3d,2]+0.5*r3dcellsizes[nr3d]
                ]
                r3dzrange = [
                    r3ddistances[nr3d,3]-0.5*r3dcellsizes[nr3d],
                    r3ddistances[nr3d,3]+0.5*r3dcellsizes[nr3d]
                ]   

                # Extract indeces of all c5dcells within current r3dcell
                c5dxrange = np.argwhere(r3dxrange[0] <= c5dgrid[np.argwhere(c5dgrid[:,0] <= r3dxrange[1]),0])[:,0]
                c5dyrange = np.argwhere(r3dyrange[0] <= c5dgrid[np.argwhere(c5dgrid[:,1] <= r3dyrange[1]),1])[:,0]
                c5dzrange = np.argwhere(r3dzrange[0] <= c5dgrid[np.argwhere(c5dgrid[:,2] <= r3dzrange[1]),2])[:,0]

                # Number of c5dcells within r3dcell (with data)
                ndustcells = 0
                ngascells = 0

                # Then loop through c5dcells within r3dcell
                for nnz in c5dzrange:
                    for nny in c5dyrange:
                        for nnx in c5dxrange:
                            # Sum all densities and temperatures (only those with data)

                            if monomer_densities[nnx,nny,nnz] > 0:
                                monomer_density += monomer_densities[nnx,nny,nnz]
                                ndustcells += 1

                            if gas_densities[nnx,nny,nnz] > 0:
                                gas_densities += gas_densities[nnx,nny,nnz]
                                ngascells += 1

                # Take average of those cells with data (and average with respect
                # to number of cells containing data only)
                if ndustcells > 0:
                    monomer_density /= ndustcells
                if ngascells > 0:
                    gas_densities /= ngascells

                # Save ratio times grainsize_constants
                grain_sizes = (grainsize_constants * monomer_density / gas_densities)**(1/3)

                # Write data to r3d files
                fsizes.write(f'{grain_sizes}\n')

                # Reset data
                monomer_density = 0
                gas_densities = 0
                ndustcells = 0
                ngascells = 0

                # Some progress bar info
                if int(nr3d/nleafs*100) == 25 and progbar == 0:
                    progbar += 1
                    print('Finished 25 per cent of the grid.')

                if int(nr3d/nleafs*100) == 50 and progbar == 1:
                    progbar += 1
                    print('Finished 50 per cent of the grid.')

                if int(nr3d/nleafs*100) == 75 and progbar == 2:
                    progbar += 1
                    print('Finished 75 per cent of the grid.')

    # End functions with aknowledgements
    print(f'C5D grain sizes:\n    grain_sizes_{phase}.dat\nDONE\n')



    # Put these 




    #sizes = np.repeat(grainsize_constants,Ncells)

    #for nn in range(Ncells):
    #    sizes[nn] *= monomer_densities[nn]/gas_densities




    # TODO for now we return the array
    # late we will have to make a binned list
    # and round each number in this array to nearest of the bins
    #return sizes


# TODO
# function that bins grain sizes
#
# save these in a file
# grain_sizes_binned.dat

def bin_grainsizes(
        sizes=np.array([1,2,3]),
        nbins:int=10,
    ):

    # number of grain size bins
    nbins = 10
    min_size = sizes.min
    max_size = sizes.max

    size_bins = np.logspace(min_size,max_size,nbins)

    for nn,size in enumerate(sizes):
        # some intelligent way of rounding each size to nearest in the bins
        1+1
    
    return 'hej'



