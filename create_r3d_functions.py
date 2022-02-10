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
    print(f"Size of base cell: {basecubesize} AU")
    print(f"Number of base cells along one side of the grid: {nxyz}")
    print(f"Distances to refinement limits: {refinementlist} AU")
    print(f"Number refinements: {nrefines}")

    # Change units to cm
    basecubesize *= AUcm
    refinementlist = [dist*AUcm for dist in refinementlist]

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
    print(f"Length of total side of whole grid: {gridedge/AUcm:.2f} AU")
    print("")

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
    print("")

    # Children grid coordinate corrections (see R3D manual for the matrix 
    # for the order of child cells inside a parent cell).
    # Each sublist is a list of coord-corrections for each layer
    # gridcorrx[level of refinement][coordinate correction of smaller cube]
    # First level are zeros so that 
    gridcorrx = [np.array([0,0,0,0,0,0,0,0])]
    gridcorry = [np.array([0,0,0,0,0,0,0,0])]
    gridcorrz = [np.array([0,0,0,0,0,0,0,0])]

    for nr in range(nrefines):
        gridcorrx.append(np.array([-1,1,-1,1,-1,1,-1,1]) * 0.5 * smallcubesize[nr])
        gridcorry.append(np.array([-1,-1,1,1,-1,-1,1,1]) * 0.5 * smallcubesize[nr])
        gridcorrz.append(np.array([-1,-1,-1,-1,1,1,1,1]) * 0.5 * smallcubesize[nr])

    # Define refinement arrays
    # I.e. add refinements within a list of radii
    refinearrayprev = np.array([0])

    for nr in range(nrefines+1):

        # Initial refinement list
        refinearray = np.zeros(nbasecubes + 8*np.size(np.where(refinearrayprev > 0)))

        # Initial grid coordinates
        nx,ny,nz = 0,0,0

        for nn in range(nbasecubes):

            # Add refinement cell coordinate corrections
            gridx = 0
            gridy = 0
            gridz = 0

            # Compute distances to each refined cell for each level
            # Up to the current refinement level

            # TODO: doesn't enter this now
            # The main loop needs to run at least once
            # which creates the problem that refinementlist[nr] is out of index
            # why I have this if statement
            # so now the refinearray has the correct size but I don't
            # get any ones in it. There should be 8 ones
            # each folloed by 8 zeros
            if nr < nrefines:
                for nnr in range(nr+1):
                    for nsmall in range(8):

                        gridx += gridcorrx[nnr][nsmall]
                        gridy += gridcorry[nnr][nsmall]
                        gridz += gridcorrz[nnr][nsmall]
                        
                        # Add distances to base cell distances
                        refdistance = np.sqrt(
                            (griddist[nx] + gridx)**2 + 
                            (griddist[ny] + gridy)**2 + 
                            (griddist[nz] + gridz)**2)

                        # Add layers of refinements
                        if refdistance <= refinementlist[nr]:
                            refinearray[nn + nnr*nsmall] = nr+1

            # Move base cell coordinates
            nx,ny,nz = movecoordinates(nxyz,nx,ny,nz)
        
        refinearrayprev = refinearray
    
    print(refinearray)
    print(nxyz**3)
    print(len(refinearray))
        
    



            



