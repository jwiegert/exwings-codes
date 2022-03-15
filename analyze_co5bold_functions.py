# functions and tools for loading and analyzing output
# from co5bold.

# Import various libraries
import cython
import numpy as np
from scipy.io.idl import readsav

# My own libraries
import create_r3d_functions as c3d

# Define useful numbers
AUcm = cython.declare(cython.float ,1.49598e13) # cm

# ============================================================

# TODO
# do I need a function that extracts things like cellsizes also?


# Load co5bold grid cell coordinates
def load_grid_coordinates(
        savpath:str='../co5bold_data/dst28gm06n056/st28gm06n056_140.sav'
    ):
    
    # Load sav-file
    c5ddata = readsav(savpath)

    # Extract data - TODO MAY HAVE TO CHOSE THIS MANUALLY DEPENDING ON FILE
    c5ddata = c5ddata['ful']

    # Get number of gridcells from co5bold data
    nc5dedge = np.size(c5ddata['Z'][0][0][16])

    # Declare output array
    c5dgrid = np.zeros((nc5dedge,3))

    for nn in range(nc5dedge):

        # Save x,y,z distances in cm
        c5dgrid[nn,0] = c5ddata['Z'][0][0][16][0][0][nn]
        c5dgrid[nn,1] = c5ddata['Z'][0][0][19][0][nn][0]
        c5dgrid[nn,2] = c5ddata['Z'][0][0][22][nn][0][0]

    return c5dgrid


# Load co5bold grid cell courners and cell sizes
def load_grid_cellsizes(
        savpath:str='../co5bold_data/dst28gm06n056/st28gm06n056_140.sav'
    ):
    
    # Load sav-file
    c5ddata = readsav(savpath)

    # Extract data - TODO MAY HAVE TO CHOSE THIS MANUALLY DEPENDING ON FILE
    c5ddata = c5ddata['ful']

    # Get number of gridcells from co5bold data
    nc5dedge = np.size(c5ddata['Z'][0][0][16])

    # Declare array for cellcourners
    cellcourners = np.zeros((nc5dedge,3))
    
    # And lists for cellsizes
    cellsizesx = []
    cellsizesy = []
    cellsizesz = []

    for nn in range(nc5dedge):
        cellsizesx.append(
            (c5ddata['Z'][0][0][25][0][0][nn+1] - c5ddata['Z'][0][0][25][0][0][nn])
        )
        cellsizesy.append(
            (c5ddata['Z'][0][0][28][0][nn+1] - c5ddata['Z'][0][0][28][0][nn])[0]
        )
        cellsizesz.append(
            (c5ddata['Z'][0][0][31][nn+1] - c5ddata['Z'][0][0][31][nn])[0][0]
        )

        cellcourners[nn,0] = c5ddata['Z'][0][0][25][0][0][nn]
        cellcourners[nn,1] = c5ddata['Z'][0][0][28][0][nn][0]
        cellcourners[nn,2] = c5ddata['Z'][0][0][31][nn][0][0]
    
    # Add final grid courner
    cellcourners[-1,0] = c5ddata['Z'][0][0][25][0][0][-1]
    cellcourners[-1,1] = c5ddata['Z'][0][0][28][0][-1][0]
    cellcourners[-1,2] = c5ddata['Z'][0][0][31][-1][0][0]

    # Extract minimum grid size
    # TODO change this later to lists/arrays with cell sizes instead
    cellsize = (min(cellsizesx) + min(cellsizesy) + min(cellsizesz))/3

    return cellsize,cellcourners


# Extract co5bold densities into a separate array
@cython.cfunc
def load_star_densities(
        savpath:str='../co5bold_data/dst28gm06n056/st28gm06n056_140.sav'
    ):
    
    # Load sav-file
    c5ddata = readsav(savpath)

    # Extract data - TODO MAY HAVE TO CHOSE THIS MANUALLY DEPENDING ON FILE
    c5ddata = c5ddata['ful']

    # Get number of gridcells from co5bold data
    nc5dedge = cython.declare(cython.int, np.size(c5ddata['Z'][0][0][16]))

    # Declare variables
    nx = cython.declare(cython.int)
    ny = cython.declare(cython.int)
    nz = cython.declare(cython.int)

    # Declare np.array
    c5dstar_densities = np.zeros((nc5dedge,nc5dedge,nc5dedge))

    # Extract densities
    for nx in range(nc5dedge):
        for ny in range(nc5dedge):
            for nz in range(nc5dedge):
                c5dstar_densities[nx,ny,nz] = c5ddata['Z'][0][0][34][nx][ny][nz]
    
    return c5dstar_densities


# TODO
# Extract co5bold densities to r3d grid
@cython.cfunc
def create_star():

    # Extract and create input data with
    # duststar densities
    # duststar temperatures
    # duststar opacities

    return 'hej'




# this decorator when declaring what's included in the functions
# @cython.locals(a=cython.int)







