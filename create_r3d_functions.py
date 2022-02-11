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

def create_grid(basecubesize, nxyz, refinementlist):
    """
    Info here!
    """
    # TODO: add info in infotext

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

        # Add first refinement

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
                                        refinearray2[
                                            nn + sum(counter)*8 +child1+child2+2
                                        ] = 3
                                
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
                                                    refinearray3[
                                                        nn + sum(counter)*8 + child1+child2+child3+3
                                                    ] = 4
                                        
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

    # Check results
    print(refinearraysave)


    # Print amr_grid
    print("Writing amr_grid.inp")

    nbranch = int(np.size(refinearraysave))
    nleafs = int(nbranch - (nbranch - nxyz*nxyz*nxyz) / 8)

    with open('../amr_grid.inp', 'w'):

        # TODO update to nicer code here!



    f = open('amr_grid.inp', 'w')
    f.writelines(['1\n\n',\
                '1\n',\
                '0\n',\
                '0\n\n',\
                '1 1 1\n',\
                str(nxyz),' ',str(nxyz),' ',str(nxyz),'\n\n',\
                str(nrefines),' ',str(nleafs),' ',str(nbranch),'\n\n'])
    #
    # 2. The courners of the grid
    #
    for N in range(3):
        for nn in range(nxyz+1):
            f.writelines([str(gridcourners[nn]).strip('[]'),'\n'])
        f.writelines(['\n'])
    f.writelines(['\n'])
    #
    # 3. The refinements
    #
    refinearraysave = refinearraysave.astype(int)
    for nn in range(nbranch):
        f.writelines([str(refinearraysave[nn]).strip('[]'),'\n'])
    f.close()


        
# One function for computing distance to cell and adding number nr to the refinearray
# functions to compute the different indeces for refinearray, where to add numbers



            



