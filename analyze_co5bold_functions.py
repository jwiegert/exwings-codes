# functions and tools for loading and analyzing output
# from co5bold.

# Import various libraries
import cython
import numpy as np
from scipy.io.idl import readsav

# My own libraries
import create_r3d_functions as c3d
import analyze_r3d_functions as a3d

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
# create grid in this function also? NO - just use the pre-made function
# But also create a function that gives numbers from c5d-grid that can
# be used for r3d


# nc5dedge = np.size(c5dgrid[:,0])
# c5dstar_densities = a5d.load_star_densities()


def create_star(
        # All necessary inputs
        # path to 
    ):


    # Load R3D grid
    print('Loading R3D grid')
    # amrpath:str='../amr_grid.inp'
    nxyz,nrefines,nleafs,nbranch,gridedge = a3d.load_gridprops()
    # gridpath:str='../grid_distances.csv',
    # amrpath:str='../amr_grid.inp',
    r3ddistances = a3d.load_griddistances()
    # sizepath:str='../grid_cellsizes.csv',
    # amrpath:str='../amr_grid.inp',
    r3dcellsizes = a3d.load_cellsizes()

    # Load C5D grid
    print('Loading C5D grid properties')
    # savpath:str='../co5bold_data/dst28gm06n056/st28gm06n056_140.sav'
    # perhaps also the name of the dict-parameter
    c5dcellsize,c5dcourners = load_grid_cellsizes()
    # savpath:str='../co5bold_data/dst28gm06n056/st28gm06n056_140.sav'
    c5dgrid = load_grid_coordinates()
    nc5dedge = np.size(c5dgrid[:,0])

    # Load C5D densities
    print('Loading C5D star-densities')
    # savpath:str='../co5bold_data/dst28gm06n056/st28gm06n056_140.sav'
    c5dstar_densities = load_star_densities()


    # Start working :)
    print('Translating C5D data to R3D data')

    # Declare stuff for the loops
    stararray = np.zeros(nleafs)
    #startemperature = np.zeros(nleafs) # not used yet, next step
    progbar = 0

    # Adaptive range used for when cellsizes are similar or equal - fixes
    # bug where I got a lot of zero-cells
    adaptive_range = c5dcellsize/r3dcellsizes.min() * 1.2

    # Check so that the the c5dcells are not larger than the r3d's smallest cells
    if r3dcellsizes.min() < c5dcellsize:
        print('\nERROR')
        print('    R3D grid resolution is higher than C5D grid, stopping')
        print('    No output is given. Change your R3D grid cells to something larger.\n')
    
    else:
        # Otherwise loop over r3d grid
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

            # Number of c5dcells within r3dcell
            nchildcells = c5dxrange.size*c5dyrange.size*c5dzrange.size

        # TODO HÄR ÄR JAG








    # Extract and create input data with
    # TODO duststar densities
    # TODO duststar temperatures
    # TODO duststar opacities

    return 'hej'



# @cython.cfunc
# this decorator when declaring what's included in the functions
# @cython.locals(a=cython.int)
#







