# Functions for creating input data for R3d simulations for Exwings project
# ------------------------------------------------------------ #
# Various useful packages
import os
import sys
import csv
from matplotlib.pyplot import grid
import numpy as np
from datetime import datetime

# Might be used later
#import os
#import scipy as s
#from scipy.integrate import quad

# ------------------------------------------------------------ #
# Shorter simpler functions

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


# Load grid proparties
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




# function that combines different dust-specie density files into one
    # headre of dust_density
    # 1
    # nleafs
    # number dust species


# ------------------------------------------------------------ #

def create_grid(basecubesize:float, nxyz:int, refinementlist:list, savegrid:str='y'):
    """
    Creates grid for Radmc3D simulations and related informaton 
    files used for analysis and creating dust envelopes. Grid is 
    octree cubic. It is refined (higher resolution) closer to 
    the central of the grid. A maximum of four (4) levels of 
    refinements is allowed in this version of the function.
    
    INPUTS
    ------
    basecubesize: length of side of base cells in AU (these are cubes) [int or float]
    
    nxyz: number of base cells along one side of the whole grid [even number, int]
    
    refinementlist: list of radial distances in AU to each level of refinement [float,float], no more than 4 numbers!
    
    savegrid: default set to 'y' [str]. If not 'y', then no grid_distances.csv or grid_cellsizes.csv will be saved. These are useful for analysing inputs and outputs of R3D!

    OUTPUTS
    -------
    amr_grid.inp : grid file for R3D-simulations

    Optional: 
    grid_distances.csv : array of radial, x, y, and z distances to each grid cell in cm
    grid_cellsizes.csv : array of sizes of each grid cell in cm
    (both have same order as in dust_density.inp and dust_temperature)
    """

    # Basic definitions
    AUcm = 1.49598e13
    nrefines = len(refinementlist)

    if nrefines > 4:
        sys.exit(f'ERROR: this is hard coded to allow a maximum of 4 grid refinements. You specified {nrefines} refinements. STOPPING')

    # Info text
    print('Creating amr_grid with octree refinement.')
    print(f'Size of base cell: {basecubesize} AU')
    print(f'Number of base cells along one side of the grid: {nxyz}')
    print(f'Distances to refinement limits: {refinementlist} AU')
    print(f'Number refinements: {nrefines}')

    # Change units to cm
    basecubesize *= AUcm
    refinementlist = [dist*AUcm for dist in refinementlist]

    # Make sure the nxyz is even, if not warn and change:
    if nxyz%2 != 0:
        nxyz += 1
        print(f'Warning, number of base cells is not even, it is now {nxyz}.\n')
    
    # Create basic parameters of the grid
    #   nbasecubes : total number of base cells
    #     gridedge : total size of the grid side
    # gridcourners : coordinates of base grid courners
    #     griddist : list of distances to center of grid (not for R3D)
    nbasecubes     = int(nxyz * nxyz * nxyz)
    gridedge       = nxyz * basecubesize
    gridcourners   = np.linspace(-gridedge*0.5,gridedge*0.5,nxyz+1)
    print(f'Length of total side of whole grid: {gridedge/AUcm:.2f} AU')
    print('')

    griddist       = np.zeros(nxyz)
    for nn in range(nxyz):
        griddist[nn]  = 0.5 * (gridcourners[nn] + gridcourners[nn+1])

    # Compute grid properties

    # Base cube size
    basecubesize   = gridcourners[1] - gridcourners[0]

    # Children cube sizes
    smallcubesize = []
    smallcubesize.append(0.5 * basecubesize)
    print(f'Child cell size 1: {smallcubesize[0]/AUcm} AU')

    for nn in range(nrefines-1):
        smallcubesize.append(0.5 * smallcubesize[nn])
        print(f'Child cell size {nn+2}: {smallcubesize[nn+1]/AUcm} AU')
    print('')

    # Children grid coordinate corrections (see R3D manual for the matrix 
    # for the order of child cells inside a parent cell).
    # Each sublist is a list of coord-corrections for each layer
    # gridcorrx[level of refinement][coordinate correction of smaller cube]
    gridcorrx = []
    gridcorry = []
    gridcorrz = []

    for nr in range(nrefines):
        gridcorrx.append(np.array([-1,1,-1,1,-1,1,-1,1]) * 0.5 * smallcubesize[nr])
        gridcorry.append(np.array([-1,-1,1,1,-1,-1,1,1]) * 0.5 * smallcubesize[nr])
        gridcorrz.append(np.array([-1,-1,-1,-1,1,1,1,1]) * 0.5 * smallcubesize[nr])

    # Define refinement arrays
    # I.e. add refinements within a list of radii
    # I need to fix this to something neater, ie recurring functions instead

    if nrefines > 0:

        # Add first refinement:
        # Define first arrays
        refinearray0 = np.zeros(nbasecubes)
        nx,ny,nz = 0,0,0

        # Loop over base cells
        for nn in range(nbasecubes):

            # Compute distances to each base cell
            refdistance = np.sqrt(
                griddist[nx]**2 + griddist[ny]**2 + griddist[nz]**2)

            # Add ones for each refined cell
            if refdistance < refinementlist[0]:
                refinearray0[nn] = 1

            # Move coordinates
            nx,ny,nz = movecoordinates(nxyz,nx,ny,nz)
        
        # Create new refinearray
        refinearray1 = np.zeros(nbasecubes + 8*np.size(np.where(refinearray0 > 0)))

        # Reset counters
        nx,ny,nz = 0,0,0
        counter = [0 for _ in range(nrefines)]

        for nn in range(nbasecubes):

            # Add 8 zeros after each 1 in refinearray
            if refinearray0[nn] == 1:
                refinearray1[nn + counter[0]*8] = 1

                # Add second level of refinements
                if nrefines > 1:

                    # Loop over these children cells
                    for child1 in range(8):

                        # Check distance to cells
                        refdistance = np.sqrt(
                            (griddist[nx] + gridcorrx[0][child1])**2 + 
                            (griddist[ny] + gridcorry[0][child1])**2 + 
                            (griddist[nz] + gridcorrz[0][child1])**2
                        )

                        if refdistance < refinementlist[1]:
                            refinearray1[nn + sum(counter)*8 + child1+1] = 2
                    
                # Update counter
                counter[0] += 1

            # Move coordinates
            nx,ny,nz = movecoordinates(nxyz,nx,ny,nz)

        if nrefines > 1:

            # Create new refinearray
            refinearray2 = np.zeros(nbasecubes + 8*np.size(np.where(refinearray1 > 0)))

            # Reset counters
            nx,ny,nz = 0,0,0
            counter = [0 for _ in range(nrefines)]


            for nn in range(nbasecubes):

                # Add 1s to correct positions in new refinearray
                if refinearray0[nn] == 1:
                    refinearray2[nn + sum(counter)*8] = 1

                    for child1 in range(8):
                        if refinearray1[nn + counter[0]*8 + child1+1] == 2:

                            refinearray2[nn + sum(counter)*8 + child1+1] = 2

                            # Search for second level of refinements
                            if nrefines > 2:

                                # Then loop over second level children
                                for child2 in range(8):
                                
                                    # Check distance to cells
                                    refdistance = np.sqrt(
                                        (griddist[nx] + gridcorrx[0][child1] + gridcorrx[1][child2])**2 + 
                                        (griddist[ny] + gridcorry[0][child1] + gridcorry[1][child2])**2 + 
                                        (griddist[nz] + gridcorrz[0][child1] + gridcorrz[1][child2])**2
                                    )

                                    # Add 3s to the refined cells
                                    if refdistance < refinementlist[2]:
                                        refinearray2[nn + sum(counter)*8 +child1+child2+2] = 3
                                
                            counter[1] += 1
                    counter[0] += 1
                # Move coordinates
                nx,ny,nz = movecoordinates(nxyz,nx,ny,nz)
            
            if nrefines > 2:

                # Create new refinearray
                refinearray3 = np.zeros(nbasecubes + 8*np.size(np.where(refinearray2 > 0)))

                # Reset counters
                nx,ny,nz = 0,0,0
                counter = [0 for _ in range(nrefines)]

                for nn in range(nbasecubes):

                    # Add 1s to correct positions in new refinearray
                    if refinearray0[nn] == 1:
                        refinearray3[nn + sum(counter)*8] = 1

                        for child1 in range(8):
                            
                            # Add 2s to the correct positions
                            if refinearray1[nn + 8*counter[0] + child1+1] == 2:
                                refinearray3[nn + sum(counter)*8 + child1+1] = 2

                                for child2 in range(8):

                                    # Add 3s to the correct positions
                                    if refinearray2[nn + 8*(counter[0]+counter[1]) + child1+child2+2] == 3:
                                        refinearray3[nn + sum(counter)*8 + child1+child2+2] = 3

                                        # Search for third refinements
                                        if nrefines > 3:
                                            for child3 in range(8):

                                                refdistance = np.sqrt(
                                                    (griddist[nx] + gridcorrx[0][child1] + gridcorrx[1][child2] + gridcorrx[2][child3])**2 + 
                                                    (griddist[ny] + gridcorry[0][child1] + gridcorry[1][child2] + gridcorry[2][child3])**2 + 
                                                    (griddist[nz] + gridcorrz[0][child1] + gridcorrz[1][child2] + gridcorrz[2][child3])**2
                                                )

                                                if refdistance < refinementlist[3]:
                                                    refinearray3[nn + sum(counter)*8 + child1+child2+child3+3] = 4
                                        
                                        # Update counters
                                        counter[2] += 1
                                counter[1] += 1
                        counter[0] += 1
                    # Move coordinates
                    nx,ny,nz = movecoordinates(nxyz,nx,ny,nz)

                # Create final refinementarray and save data
                # (Yes this is ugly, I will have to solve this some day with
                # nested functions instead but now it was just way too complex)
                refinearraysave = np.zeros(nbasecubes + 8*np.size(np.where(refinearray3 > 0)))
                counter = 0
                for nn in range(np.size(refinearray3)):
                    if refinearray3[nn] == 1:
                        refinearraysave[nn+counter] = 1
                    if refinearray3[nn] == 2:
                        refinearraysave[nn+counter] = 1
                    if refinearray3[nn] == 3:
                        refinearraysave[nn+counter] = 1
                    if refinearray3[nn] == 4:
                        refinearraysave[nn+counter] = 1
                        counter += 8

            else:
                refinearraysave = np.zeros(nbasecubes + 8*np.size(np.where(refinearray2 > 0)))
                counter = 0
                for nn in range(np.size(refinearray2)):
                    if refinearray2[nn] == 1:
                        refinearraysave[nn+counter] = 1
                    if refinearray2[nn] == 2:
                        refinearraysave[nn+counter] = 1
                    if refinearray2[nn] == 3:
                        refinearraysave[nn+counter] = 1
                        counter += 8

        else:
            refinearraysave = np.zeros(nbasecubes + 8*np.size(np.where(refinearray1 > 0)))
            counter = 0
            for nn in range(np.size(refinearray1)):
                if refinearray1[nn] == 1:
                    refinearraysave[nn+counter] = 1
                if refinearray1[nn] == 2:
                    refinearraysave[nn+counter] = 1
                    counter += 8

    # Print amr_grid
    print('Writing amr_grid.inp')

    nbranch = int(np.size(refinearraysave))
    nleafs = int(nbranch - (nbranch - nxyz*nxyz*nxyz) / 8)

    # Save time of creation of amr_grid.inp
    dt_string = datetime.now().strftime('%Y-%m-%d %H:%M')

    with open('../amr_grid.inp', 'w') as f:

        # Write header containing:
        # Number of Base cells in earch axis
        # Number of cells and children cells etc
        f.write(f'1\n\n1\n0\n0\n\n1 1 1\n{nxyz} {nxyz} {nxyz}\n\n{nrefines} {nleafs} {nbranch}\n\n')

        # Write coordinates of courners of base cells in cm
        for _ in range(3):
            for gridc in gridcourners:
                f.write(f'{gridc}\n')
            # Add empty line between the 3 dimensions
            f.write('\n')
        # Add additional empty line before list of refinements
        f.write('\n')

        # Write refinements
        for gridc in refinearraysave:
            f.write(f'{int(gridc)}\n')
    
    print('Finished amr_grid.inp')
    print('')
    
    # Print grid_distances.csv and grid_cellsizes.csv

    if savegrid == 'y' or savegrid == 'yes' or savegrid == 'Y':
        print('Writing grid_distances.csv and grid_cellsizes.csv (not necessary for Radmc3d, but useful for pre/portprocessing of your model files. They have the same order as dust_densities.inp)')

        # Declare an array for the distances to the centre of each cell
        # Radial and x,y,z distances
        griddistances = np.zeros((nleafs,4))
        # And an array for the sizes of each cell
        gridsizes = np.zeros(nleafs)

        # Declare a few counters
        # nbig is the index of all cells
        nbig,nx,ny,nz = 0,0,0,0
        counter = [0 for _ in range(nrefines)]

        for nn in range(np.size(refinearray0)):

            # Distances to base cells
            if refinearray0[nn] == 0:
                
                # Distances along axes
                griddistances[nbig,1] = griddist[nx]
                griddistances[nbig,2] = griddist[ny]
                griddistances[nbig,3] = griddist[nz]

                # Radial distance
                griddistances[nbig,0] = np.sqrt(
                    griddistances[nbig,1]**2 + griddistances[nbig,2]**2 + griddistances[nbig,3]**2)

                # Cell sizes
                gridsizes[nbig] = basecubesize

                # Move coordinates
                nx,ny,nz = movecoordinates(nxyz,nx,ny,nz)
                nbig += 1
            
            # Distances to first refinement
            if refinearray1[nn + counter[0]*8] == 1:
                for child1 in range(8):
                    if refinearray1[nn + counter[0]*8 + child1+1] == 0:

                        # Distances along axes
                        griddistances[nbig,1] = griddist[nx] + gridcorrx[0][child1]
                        griddistances[nbig,2] = griddist[ny] + gridcorry[0][child1]
                        griddistances[nbig,3] = griddist[nz] + gridcorrz[0][child1]

                        # Radial distance
                        griddistances[nbig,0] = np.sqrt(
                            griddistances[nbig,1]**2 + griddistances[nbig,2]**2 + griddistances[nbig,3]**2)

                        # Cell sizes
                        gridsizes[nbig] = smallcubesize[0]

                        # Move main index
                        nbig += 1
                    
                    # Check for second refinment
                    if nrefines > 1 and refinearray2[nn + (counter[0]+counter[1])*8 + child1+1] == 2:
                        for child2 in range(8):
                            if refinearray2[nn + (counter[0]+counter[1])*8 + child1+child2+2] == 0:

                                # Distances along axes
                                griddistances[nbig,1] = griddist[nx] + gridcorrx[0][child1] + gridcorrx[1][child2]
                                griddistances[nbig,2] = griddist[ny] + gridcorry[0][child1] + gridcorry[1][child2]
                                griddistances[nbig,3] = griddist[nz] + gridcorrz[0][child1] + gridcorrz[1][child2]

                                # Radial distance
                                griddistances[nbig,0] = np.sqrt(
                                    griddistances[nbig,1]**2 + griddistances[nbig,2]**2 + griddistances[nbig,3]**2)

                                # Cell sizes
                                gridsizes[nbig] = smallcubesize[1]

                                # Move main index
                                nbig += 1
                            
                            # Check for third refinement
                            if nrefines > 2 and refinearray3[nn + sum(counter)*8 + child1+child2+2] == 3:
                                for child3 in range(8):
                                    if refinearray3[nn + sum(counter)*8 + child1+child2+child3+3] == 0:

                                        # Distances along axes
                                        griddistances[nbig,1] = griddist[nx] + gridcorrx[0][child1] + gridcorrx[1][child2] + gridcorrx[2][child3]
                                        griddistances[nbig,2] = griddist[ny] + gridcorry[0][child1] + gridcorry[1][child2] + gridcorry[2][child3]
                                        griddistances[nbig,3] = griddist[nz] + gridcorrz[0][child1] + gridcorrz[1][child2] + gridcorrz[2][child3]

                                        # Radial distance
                                        griddistances[nbig,0] = np.sqrt(
                                            griddistances[nbig,1]**2 + griddistances[nbig,2]**2 + griddistances[nbig,3]**2)

                                        # Cell sizes
                                        gridsizes[nbig] = smallcubesize[2]

                                        # Move main index
                                        nbig += 1

                                    # Check for fourth and final refinement
                                    if nrefines > 3 and refinearray3[nn + sum(counter)*8 + child1+child2+child3+3] == 4:
                                        for child4 in range(8):
                        
                                            # Distances along axes
                                            griddistances[nbig,1] = griddist[nx] + gridcorrx[0][child1] + gridcorrx[1][child2] + gridcorrx[2][child3] + gridcorrx[3][child4]
                                            griddistances[nbig,2] = griddist[ny] + gridcorry[0][child1] + gridcorry[1][child2] + gridcorry[2][child3] + gridcorry[3][child4]
                                            griddistances[nbig,3] = griddist[nz] + gridcorrz[0][child1] + gridcorrz[1][child2] + gridcorrz[2][child3] + gridcorrz[3][child4]

                                            # Radial distance
                                            griddistances[nbig,0] = np.sqrt(
                                                griddistances[nbig,1]**2 + griddistances[nbig,2]**2 + griddistances[nbig,3]**2)

                                            # Cell sizes
                                            gridsizes[nbig] = smallcubesize[3]

                                            # Move main index
                                            nbig += 1

                                # Increase counters
                                counter[2] += 1
                        counter[1] += 1
                counter[0] += 1

                # Move base cell coordinates
                nx,ny,nz = movecoordinates(nxyz,nx,ny,nz)
        
        # Save distances to grid cells
        # This has the same order as the dust_densities later on
        with open('../grid_distances.csv', 'w') as f:
            f.write('# Distances to centrum of grid cells (same order as in dust_densities.inp) in cm from centrum of grid\n')
            f.write('# Radial, x-distance, y-distance, z-distance\n')
            f.write('# ----------\n')
            for gridc in griddistances:
                f.write(f'{gridc[0]},{gridc[1]},{gridc[2]},{gridc[3]}\n')
        print('Finished grid_distances.csv')

        # Save sizes of grid cells
        # This also has the same order as the dust_densities files
        with open('../grid_cellsizes.csv', 'w') as f:
            f.write('# Sizes of grid cells (same order as in dust_densities.inp) in cm\n')
            f.write('# ----------\n')
            for gridsize in gridsizes:
                f.write(f'{gridsize}\n')
        print('Finished grid_cellsizes.csv\n')
    
    # Print grid_info-file, text file with summary o grid info in readable form
    print('Writing grid_info.txt')
    with open('../grid_info.txt', 'w') as f:
        f.writelines([
            f'Information to amr_grid created at {dt_string}\n\n',
            f'                             Size of base cell: {basecubesize/AUcm} AU\n',
            f'Number of base cells (of one side of the grid): {nxyz}\n',
            f'                        Length of side of grid: {gridedge/AUcm} AU\n'
            f'                            Number refinements: {nrefines}\n',
            f'         Total number of cells (numb of leafs): {nleafs}\n\n'
        ])

        for nn,refdist in enumerate(refinementlist):
            f.write(f'Radial distance to refinement {nn+1}: {refdist/AUcm} AU\n')
        f.write('\n')

        for nn,cellsize in enumerate(smallcubesize):
            f.write(f'Child cell size {nn+1}: {smallcubesize[nn]/AUcm} AU\n')

    print('Finished grid_info.txt\n')

    # Finish function
    print('Done')

# ------------------------------------------------------------ #

# def create duststar

def create_duststar(
        Mstar:float = 1,
        Rstar:float = 100,
        Teff:float = 2700
    ):
    """
    INFO HERE

    Should output
    dust_density.inp
    dust_temperature.dat
    dust_kappa


    sfärisk klump i mitten, 
    hög abs, 0 scat. 
    Jämför med svartkroppsstjärna. - 
    1Msun, 7000Lsun, 2700K
    """

    # Useful units
    #AUcm = 1.49598e13 # cm
    Msol = 1.989e33 # g
    Rsol = 6.955e10 # cm

    # Stellar props
    Mstar *= Msol # star mass in g
    Rstar *= Rsol # star radius in cm

    # Load griddistances and size of grid
    griddistances = load_griddistances()
    nleafs = np.size(griddistances[:,0])
  
    # Compute average density of star (R3d uses g/cm3)
    stardensity = Mstar / (4/3 * np.pi * Rstar**3)

    # Create and fill density and temperature array
    # fill with average density (future: rho(R))
    # and with Teff (future: T(R))
    densities = np.zeros(nleafs)
    temperatures = np.zeros(nleafs)

    for nn in range(nleafs):
        if griddistances[nn,0] <= Rstar:
            densities[nn] = stardensity
            temperatures[nn] = Teff
    
    # Write dust_density
    print('Writing dust_density_01.inp')
    with open('../dust_density_01.inp', 'w') as f:

        # Header of dust_density (and dusttemperature?)
        # 1
        # nleafs
        # number dust species
        f.write(f'1\n{int(nleafs)}\n1\n')
        for density in densities:
            f.write(f'{density}\n')

    print('Finished dust_density_01.inp')





"""
def create_spheredensity():

    # Load grid (should perhaps be its own function)

    # Compute normalization density (ie zero density)
    # 


    return 'hej'
""";