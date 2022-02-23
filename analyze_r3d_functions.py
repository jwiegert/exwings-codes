# Various functions useful for analyzing in and out-data of RADMC3D
# ------------------------------------------------------------ #
# Useful packages
import os
import csv
import numpy as np

# Useful numbers
c = 2.998e8 # speed of light in m/s
pc = 30.857e15 # 1 parsec in m

# ------------------------------------------------------------ #
# Functions that load various stuff

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

                # Number of refinements, number of cells (nleafs), number of branches
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

# Load SED
def load_spectrum(path:str='../r3dsims/spectrum.out',distance:float=1):
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
