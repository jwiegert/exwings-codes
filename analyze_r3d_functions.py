# Various functions useful for analyzing in and out-data of RADMC3D
# ------------------------------------------------------------ #
# Useful packages
import os
import csv
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm

# Useful numbers
c = 2.998e8 # speed of light in m/s
pc = 30.857e15 # 1 parsec in m
AUcm = 1.49598e13 # AU in cm

# ------------------------------------------------------------ #
# Functions that load various r3d input data

# TODO?
# Both these below need the number of species?

# Load grid properties
def load_gridprops(amrpath:str='../amr_grid.inp'):
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
        return f'ERROR: load_gridprops can not find {amrpath}.'


# Load griddistances
def load_griddistances(
        gridpath:str='../grid_distances.csv',
        amrpath:str='../amr_grid.inp',
    ):
    """
    Loads array of distances to cells of the current grid.
    Distances are from centrum of grid (radial,x,y, and z) in centimeters.

    Inputs
    ------
    gridpath: path to grid_distances.csv
    amrpath: path to amr_grid.inp

    OUTPUT
    ------
    griddistances: np.array of size nleafs*4. Each column is distances in cm
    """

    # Check if file exists
    if os.path.exists(gridpath) == True:

        # Extract necessary info from amr_grid
        nleafs = load_gridprops(amrpath)[2]

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
        return f'ERROR: load_griddistances can not find {gridpath}.'

# Load list of grid cell sizes
def load_cellsizes(
        sizepath:str='../grid_cellsizes.csv',
        amrpath:str='../amr_grid.inp',
    ):
    """
    Loads array of grid cell sizes in centimeters
    
    Inputs
    ------
    sizepath: path to grid_cellsizes.csv
    amrpath: path to amr_grid.inp

    OUTPUT
    ------
    gridsizes: np.array of size nleafs*1 with cell sizes in same order as in dust_density
    """

    # Check if file exists
    if os.path.exists(sizepath) == True:

        # Extract necessary info from amr_grid
        nleafs = load_gridprops(amrpath)[2]

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
def load_wavelengthgrid(path:str='../wavelength_micron.inp'):
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

def load_dustdensity(
        path:str='../dust_density.inp',
        numb_specie:int=1
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

    # Read header
    with open(path,'r') as f:
        for nn,line in enumerate(f.readlines()):

            # Number of cells
            if nn == 1:
                Ncells = int(line)

            # Total number of species
            if nn == 2:
                Nspec = int(line)

    # Check that the chosen number of specie exists
    if numb_specie > Nspec:
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

        return Ncells,Nspec,dust_densities

# ------------------------------------------------------------ #
# Load output data from R3D

# Load dust_temperature-file
def load_dusttemperature(
        path:str='../dust_temperature.dat',
    ):
    """
    Load and extracts temperatures in output file dust_temperature.dat

    INPUT
    -----
    path: string with path and filename of temperature file
    
    OUTPUT
    ------
    Ncells: Number of cells in grid (nleafs)
    dust_temperatures: np.array containing temperatures in grid
    """

    # Read header
    with open(path,'r') as f:
        for nn,line in enumerate(f.readlines()):

            # Number of cells
            if nn == 1:
                Ncells = int(line)



    # Create density np.array
    dust_temperatures = np.zeros(Ncells)

    # Extract dust densities
    with open(path,'r') as f:
        for nn,line in enumerate(f.readlines()):
            if nn > 2:
                dust_temperatures[nn-3] = float(line)

    return Ncells,dust_temperatures


# Load SED
def load_spectrum(
        path:str='../r3dsims/spectrum.out',
        distance:float=1
    ):
    """
    Loads and returns SED and wavelength.

    INPUT
    -----
    path: str with path and filename of spectrum.out file
    distance: float, distance to star in pc

    OUTPUT
    ------
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


# ------------------------------------------------------------ #
# Plot various details of input and output data of R3D

def plot_onedensity_radius(
        density_path:str='../dust_density.inp',
        grid_path:str='../grid_distances.csv',
        amr_path:str='../amr_grid.inp',
        numb_specie:int=1
    ):
    """
    Plots one figure with radian density distribution of one dust species.

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
        density_path:str='../dust_density.inp',
        grid_path:str='../grid_distances.csv',
        amr_path:str='../amr_grid.inp',
    ):
    """
    Plots one figure with radial density distribution of all dust species.

    INPUT
    density_path: path to density.inp-file
    grid_path: path to grid_distances.csv'
    amr_path: path to amr_grid.inp

    OUTPUT
    Shows figure
    """

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
        ax.plot(
            griddistances[:,0],densities[nn],
            markeredgecolor=c,
            linestyle='',marker='.',markersize=1
        )
    ax.set(
        ylabel=r'Density (g\,cm$^{-3}$)',
        xlabel=r'Distance (AU)',
        title=f'Dust species 1 to {Nspec}'
    )
    fig.show()

    # Set objects for subplots (two columns, increasing number of rows)
    # ax[rows,columns]
    fig,ax = plt.subplots((-(-Nspec//2)),2)

    for nn, c in enumerate(colour):
        ax.ravel()[nn].plot(
            griddistances[:,0],densities[nn],
            markeredgecolor=c,
            linestyle='',marker='.',markersize=1
        )
        ax.ravel()[nn].set(
            ylabel=r'Density (g\,cm$^{-3}$)',
            xlabel=r'Distance (AU)',
            title=f'Dust specie {nn+1}'
        )
    fig.show()


def plot_temperature_radius(
        temperature_path:str='../dust_temperature.dat',
        grid_path:str='../grid_distances.csv',
        amr_path:str='../amr_grid.inp',
    ):
    """
    Plots one figure with radial temperature distribution.

    INPUT
    temperature_path: path to dust_temperature.dat
    grid_path: path to grid_distances.csv'
    amr_path: path to amr_grid.inp

    OUTPUT
    Shows figure
    """

    # Load griddistances
    griddistances = load_griddistances(
        gridpath=grid_path,
        amrpath=amr_path
    )/AUcm

    # load dust_temperature
    # TODO
    #temperature = 
    #def load_dusttemperature(
    #    path:str='../dust_temperature.dat',
    #):

    # Load and plots r3d density data for ONE dust specie
    fig, ax = plt.figure(), plt.axes()
    ax.plot(
        griddistances[:,0],temperature,
        linestyle='',marker='.',markersize=1
    )
    ax.set(
        ylabel=r'Density (g\,cm$^{-3}$)',
        xlabel=r'Distance (AU)',
        title='Grid temperatures'
    )
    fig.show()





# ------------------------------------------------------------ #
# Compute different storheter

def compute_luminosity(path:str='../r3dsims/spectrum.out',distance:float=1):
    """
    INFO HERE
    """

    # Load spectrum
    wavelengths,spectrum = load_spectrum(path,1)
    nwave = len(wavelengths)

    # Integrate the SED (using trapezoidal method, and change units to SI units)
    sedintegral = 0
    for nn in range(nwave-1):
        sedintegral += 0.5*(spectrum[nn] + spectrum[nn+1])*1e-26 * (c/wavelengths[nn] - c/wavelengths[nn+1])*1e6

    # Compute bolometric luminosity(?)
    luminosity = 4.*np.pi*(distance*pc)**2. * sedintegral

    return luminosity
