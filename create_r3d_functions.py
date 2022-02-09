# Functions for creating input data for R3d simulations for Exwings project
import os
import numpy as np
import scipy as s
from scipy.integrate import quad


# ------------------------------------------------------------ #
def movecoordinates(nxyz,nx,ny,nz):
    """
    Recurring function to move through the coordinates of the grid.
    nxyz : integer, the number of base grid cells (largest cubes)
    nx, ny, nz: integers, the current grid coordinates
    """

    nx     += 1
    if nx  == nxyz:
        nx =  0
        ny += 1
    if ny  == nxyz:
        nx =  0
        ny =  0
        nz += 1
    return nx,ny,nz
# ------------------------------------------------------------ #

# Req inputs:
# base grid size (in centimeters!)
# number of refinements (nrefines)
# number of base cubes (on one side of the grid)
# position of refinements, radial distances to the center of the grid

def create_grid(basecubesize, nxyz, refinementlist):
    """
    Info here!
    """

    # Basic definitions
    AUcm = 1.49598e13
    nrefines = len(refinementlist)

    # Info text
    print("Creating amr_grid with octree refinement.")
    print(f"Number refinements: {nrefines}")
    print(f"Size of base cell: {basecubesize} AU")
    print(f"Number of base cells along one side of the grid: {nxyz}")

    # Change units to cm for basecubesize
    basecubesize *= AUcm

    # Make sure the nxyz is even, if not warn and change:
    if nxyz%2 != 0:
        nxyz += 1
        print(f"Warning, number of base cells is not even, it's now {nxyz}.\n")
    
    # Create basic parameters of the grid
    #   nbasecubes : total number of base cells
    #     gridedge : total size of the grid side
    # gridcourners : coordinates of base grid courners
    #     griddist : list of distances to center of grid (not for R3D)
    nbasecubes     = int(nxyz * nxyz * nxyz)
    gridedge       = nxyz * basecubesize
    gridcourners   = np.linspace(-gridedge*0.5,gridedge*0.5,nxyz+1)

    griddist       = np.zeros(nxyz)
    for nn in range(nxyz):
        griddist[nn]  = 0.5 * (gridcourners[nn] + gridcourners[nn+1])

    # Compute grid properties

    # Base cube size
    basecubesize   = gridcourners[1] - gridcourners[0]

    # Children cube sizes
    smallcubesize = []
    smallcubesize.append(0.5 * basecubesize)
    print(f"Child cell size 0: {smallcubesize[0]/AUcm} AU")

    for nn in range(nrefines-1):
        smallcubesize.append(0.5 * smallcubesize[nn])
        print(f"Child cell size {nn+1}: {smallcubesize[nn+1]/AUcm} AU")
        
    
    # Children grid coordinates (see R3D manual for the matrix for the order of
    # child cells inside a parent cell).
    gridcorrx = []
    gridcorry = []
    gridcorrz = []

    for nn in range(nrefines):
        gridcorrx.append(np.array([-1,1,-1,1,-1,1,-1,1]) * 0.5 * smallcubesize[nn])
        gridcorry.append(np.array([-1,-1,1,1,-1,-1,1,1]) * 0.5 * smallcubesize[nn])
        gridcorrz.append(np.array([-1,-1,-1,-1,1,1,1,1]) * 0.5 * smallcubesize[nn])

    # Define refinement arrays
    # I.e. add refinements within a list of radii
    refinematrix = np.zeros(nbasecubes,nrefines)

    # Initial grid coordinates
    nx,ny,nz = 0,0,0

    for nr in range(nrefines):
        for nn in range(nbasecubes):
            # TODO: change my refinement scheme to these loops
    
    
