# Functions for creating input data for R3d simulations for Exwings project
# ------------------------------------------------------------ #
# Various useful packages
from asyncio import LimitOverrunError
import sys
import wave
import numpy as np
import scipy as s
from scipy.integrate import quad
from datetime import datetime
import os

# Might be used later
#import scipy as s
#from scipy.integrate import quad

# PyR3D-packages
# See manual for RADMC3D for references on these codes.
from bhmie.makedustopac import *

# My own functions
import analyze_r3d_functions as a3d

# Basic definitions
AUcm = 1.49598e13 # cm
Msol = 1.989e33 # g
Rsol = 6.955e10 # cm

# ------------------------------------------------------------ #
# List of functions
#
# Simpler functions
# -----------------
#
# find_zeroelements(
#    inputlist
# )
#
# movecoordinates(nxyz,nx,ny,nz)
#
# write_r3d_runscripts(
#    path = '../r3dresults/st28gm06n056/',
#    phase_list = [140,141,142],
#    sed_inclination_list = [0],
#    image_wavelength_list = [1],
#    image_inclination_list = [0],
#    image_sizeau = 7.4,
#    image_npix = 128,
# )   
#    
# 
# Grid functions
# --------------
#
# create_grid(
#    gridedge:float, 
#    nxyz:int, 
#    refinementlist:list,
#    inrefine:str=0.2,
#    savegrid:str=True
# )
# 
# create_wavelength(
#    wavelengthstart:float=0.1,
#    wavelengthend:float=2000.0,
#    nwave:int=100,
#    logscale:str='y'
# )
# 
# 
# Opacity-functions
# -----------------
#
# create_optoolscript(
#    wavelength_path:str='../wavelength_micron.inp',
#    grainum_sizes:list=[0.1],
#    grainsize_type:str='lognormal',
#    grainsize_na:int=21,
#    grainsize_power:float=-3.5,
#    specie:str='mg2sio4',
#    grain_type:str='mie'
# )
#
# create_kappaabs(
#    wavelengthpath:str='../wavelength_micron.inp',
#    optconstlist:list=['mg2sio4'],
#    agrainlist:list=[0.1],
#    matdens:float=2
# ) 
#
# load_staropacities(
#    path:str = '../star_opacities.dat'
# )
#
#
# Create simple R3D-data-functions
# --------------------------------
#
# create_duststar(
#    Mstar:float = 1,
#    Rstar:float = 100,
#    Teff:float = 2700
# )
# 
# TODO Not finished!
# create_spheredensity(
#    optconstlist:list=['mg2sio4'],
#    agrainlist:list=[0.1],
#    totaldustmass:float=1.989e33,
#    densitypower:float=-2,
#    inradius:float=3,
#    outradius:float=100
# )
# 
#
# Merge data files-functions
# --------------------------
#
# merge_dustdensities(
#    filenames:list=['dust_density.inp'],
#    modelname:str='st28gm06n052',
#    phases:list=[186,190,198],
# )
#
# merge_dustopac(
#    filenames:list=['dustopac.inp'],
#    modelname:str='st28gm06n052',
#    phases:list=[186,190,198],
# )
#
# merge_dusttemperatures(
#    filenames:list=['dust_temperature.dat'],
#    modelname:str='st28gm06n052',
#    phases:list=[186,190,198],
# )
#
# ------------------------------------------------------------ #
# Shorter simpler functions

# Searches and extracts sequences of 0s from lists/arrays
# TODO change so that it can extract sequences of any number?
def find_zeroelements(
        inputlist:list,
    ):
    """
    Searches for sequences of zeros and returns a list of lists where
    each list contains elements of sequences of zeros.

    ARGUMENTS
      inputlist: list: list or array to search through
    """

    # Declare lists to fill
    hole_lists = []
    hole_list = []

    for nn,number in enumerate(inputlist):

        # Save zeros
        if number == 0:
            hole_list.append(nn)

        # Save list when 0s ends and reset list
        if number != 0 and inputlist[nn-1] == 0 and len(hole_list) > 0:
            hole_lists.append(hole_list)
            hole_list = []
        
        # Save list if 0s continue to end of array
        if nn == len(inputlist)-1 and number == 0 and len(hole_list) > 0:
            hole_lists.append(hole_list)

    return hole_lists



# Move base cell coordinates
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


# Write runcommand-files, one for each phase and one main file that runs them in paralell
def write_r3d_runscripts(
        path:str = '../r3dresults/st28gm06n056/',
        phase_list:list = [140,141,142],
        sed_inclination_list:list = [0],
        image_wavelength_list:list = [1],
        image_inclination_list:list = [0],
        image_sizeau:float = 7.4,
        image_npix:int = 128,
    ):
    """
    Creates runcommand-files, one for each phase and one main file that runs them in paralell
    TODO: add lists for second viewing angle

    ARGUMENTS
      path = path from code-folder to main data folder (ie below the phase-folders, see example)
      phase_list = list of all phases for this model-star
      sed_inclination_list = list of inclination angles for the SED-simulations
      image_wavelength_list = list of image wavelengths in micrometres
      image_inclination_list = list of inclination angles for images
      image_sizeau = size of image-side in AU
      image_npix = number of pixels per side of images
    """

    # Automatically add / to end of path if it's missing
    if path[-1] != '/':
        path += '/'


    # Automatically adds 0 as inclination angle if none is given
    if len(image_inclination_list) < 1:
        image_inclination_list = [0]

    # Reset main runcommand file
    with open(f'{path}runcommand_main.sh', 'w') as fmain:
        fmain.write(f'#! /bin/bash\n#\n\n')

    # Loop through phase-list
    for phase in phase_list:

        # Write run-scripts for each phase
        with open(f'{path}runcommand{phase}.sh', 'w') as fr3d:

            # Enter bash-thingie
            fr3d.write(f'#! /bin/bash\n\n')

            # Enter phase-folder
            fr3d.write(f'cd {phase}\n\n')

            # Write SED-simulation-lines
            if len(sed_inclination_list) > 0:
                for inclination in sed_inclination_list:
                    fr3d.write(f'radmc3d sed incl {inclination}\n')
                    fr3d.write(f'mv spectrum.out spectrum_i{inclination}.out\n')

            # Add empty line
            fr3d.write('\n')

            # Write image-simulation lines (if any are given)
            if len(image_wavelength_list) > 0:

                for wavelength in image_wavelength_list:
                    for inclination in image_inclination_list:
                    
                        fr3d.write(f'radmc3d image nostar incl {inclination} lambda {wavelength} npix {image_npix} sizeau {image_sizeau}\n')
                        # TODO add more variables in image file name
                        fr3d.write(f'mv image.out image_i{inclination}_{wavelength}um.out\n')

        # Make phase-script-files executables
        os.system(
            f'chmod +x {path}runcommand{phase}.sh'
        )

        # Write main script that runs all phases in paralell
        with open(f'{path}runcommand_main.sh', 'a') as fmain:
            fmain.write(f'rm r3doutput_{phase}.txt\ntouch r3doutput_{phase}.txt\n./runcommand{phase}.sh | cat > r3doutput_{phase}.txt &\n\n')
    
    # Add a working ending to main-script
    with open(f'{path}runcommand_main.sh', 'a') as fmain:
        fmain.write(f'wait\necho All done\n')

    # Make main script executable
    os.system(
        f'chmod +x {path}runcommand_main.sh'
    )

    print('Finished writing run-r3d-scripts:\n    runcommand[PHASE].sh\n    runcommand_main.sh\n')


# ------------------------------------------------------------ #
# Grid functions

def create_grid(
        gridedge:float, 
        nxyz:int, 
        refinementlist:list,
        inrefine:float=0.2,
        savegrid:str=True
    ):
    """
    Creates grid for Radmc3D simulations and related informaton 
    files used for analysis and creating dust envelopes. Grid is 
    octree cubic. It is refined (higher resolution) closer to 
    the central of the grid. A maximum of four (4) levels of 
    refinements is allowed in this version of the function.
    
    INPUTS
    ------
    gridedge: total size of the grid side (in AU)
    nxyz: number of base cells along one side of the whole grid [even number, int]
    refinementlist: list of radial distances in AU to each level of refinement [float,float], no more than 4 numbers!
    inrefine: distance in AU between refinements inside the star
    savegrid: default set to 'y' [str]. If not 'y', then no grid_distances.csv or grid_cellsizes.csv will be saved. These are useful for analysing inputs and outputs of R3D!

    OUTPUTS
    -------
    amr_grid.inp : grid file for R3D-simulations
    OPTIONAL
    grid_distances.csv : array of radial, x, y, and z distances to each grid cell in cm
    grid_cellsizes.csv : array of sizes of each grid cell in cm
    (two last files have same order as in dust_density.inp and dust_temperature)
    """

    # Basic definitions
    nrefines = len(refinementlist)

    if nrefines > 4:
        sys.exit(f'ERROR: this is hard coded to allow a maximum of 4 grid refinements. You specified {nrefines} refinements. STOPPING')

    # Save the input nxyz-value for later
    oldnxyz = nxyz
    
    # Make sure the nxyz is integer
    if nxyz - int(nxyz) != 0:
        nxyz = int(round(nxyz))
        
    # Make sure nxyz is even
    if nxyz%2 != 0:
        nxyz += 1

    # Calculate new gridedge based on ratio between new and input nxyz
    if nxyz != oldnxyz:
        gridedge *= nxyz/oldnxyz
    
    # Calculate inner refinements around centrum of grid
    # Dive up given inrefine-radius in equal portions and list distances
    innerrefinements = [(nref+1)*inrefine/nrefines for nref in range(nrefines)]

    # Info text
    print('\n    Creating amr_grid with octree refinement.')
    print(f'Final length of total side of whole grid: {gridedge} AU')
    print(f'Number of base cells along one side of the grid: {nxyz}')
    print(f'Distances to outer refinement limits from centrum:\n    {refinementlist} AU')
    print(f'Distances to inner refinement limits from centrum:\n    {innerrefinements} AU')
    print(f'Number refinements: {nrefines}\n')

    # Change units to cm
    gridedge *= AUcm
    refinementlist = [dist*AUcm for dist in refinementlist]
    innerrefinements = [dist*AUcm for dist in innerrefinements]

    # Create basic parameters of the grid
    #   nbasecubes : total number of base cells
    #     gridedge : total size of the grid side
    # gridcourners : coordinates of base grid courners
    #     griddist : list of distances to center of grid (not for R3D)
    nbasecubes     = int(nxyz**3)
    gridcourners   = np.linspace(-gridedge*0.5,gridedge*0.5,nxyz+1)
    griddist       = np.zeros(nxyz)
    for nn in range(nxyz):
        griddist[nn]  = 0.5 * (gridcourners[nn] + gridcourners[nn+1])

    # Compute grid properties

    # Base cube size
    basecubesize   = gridcourners[1] - gridcourners[0]
    print(f'Size of base cell: {basecubesize/AUcm} AU')

    # Children cube sizes
    smallcubesize = []
    smallcubesize.append(0.5 * basecubesize)
    print(f'Child cell size 1: {smallcubesize[0]/AUcm} AU')

    for nn in range(nrefines-1):
        smallcubesize.append(0.5 * smallcubesize[nn])
        print(f'Child cell size {nn+2}: {smallcubesize[nn+1]/AUcm} AU')

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
            if refdistance < refinementlist[0] and refdistance >= innerrefinements[0]:
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

                        if refdistance < refinementlist[1] and refdistance >= innerrefinements[1]:
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
                                    if refdistance < refinementlist[2]  and refdistance >= innerrefinements[2]:
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

                                                if refdistance < refinementlist[3]  and refdistance >= innerrefinements[3]:
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
    
    # Print grid_distances.csv and grid_cellsizes.csv
    if savegrid == True:
        print('Writing grid_distances.csv and grid_cellsizes.csv\n(Not necessary for Radmc3d, but useful for pre/portprocessing of your model. They have the same order as dust_densities.inp)')

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

        # Save sizes of grid cells
        # This also has the same order as the dust_densities files
        with open('../grid_cellsizes.csv', 'w') as f:
            f.write('# Sizes of grid cells (same order as in dust_densities.inp) in cm\n')
            f.write('# ----------\n')
            for gridsize in gridsizes:
                f.write(f'{gridsize}\n')
    
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
            f.write(f'Radial distances to refinement {nn+1}: {innerrefinements[nn]/AUcm} - {refdist/AUcm} AU\n')
        f.write('\n')

        for nn,cellsize in enumerate(smallcubesize):
            f.write(f'Child cell size {nn+1}: {smallcubesize[nn]/AUcm} AU\n')

    # Finish function
    print('Create grid:\n    amr_grid.inp\n    grid_info.txt')
    if savegrid == True:
        print('    grid_distances.csv\n    grid_cellsizes.csv')
    print('DONE\n')


# Write wavelengthgrid
def create_wavelength(
        wavelengthstart:float=0.1,
        wavelengthend:float=2000.0,
        nwave:int=100,
        logscale:str='y'
    ):
    """
    Creates wavelength grid file for R3D. Lists all wavelength points in microns.

    INPUT
    -----
    wavelengthstart: Minimum wavelength in microns
    wavelengthend: Maximum wavelenth in microns
    nwave: Number of wavelength grid points
    logscale: 'y' (Yes) for logarithmic wavelength grid (Default)
    """

    print('Creating wavelength grid')

    if logscale == 'y' or logscale == 'yes' or logscale == 'Y' or logscale == 'Yes':

        print('Logarithmic wavelength grid was chosen')
        wavelengths = np.logspace(
            np.log10(wavelengthstart),
            np.log10(wavelengthend),
            nwave
        )

    else:

        print('Linear wavelength grid was chosen')
        wavelengths = np.linspace(
            wavelengthstart,
            wavelengthend,
            nwave
        )

    print('Writing wavelength_micron.inp')
    with open('../wavelength_micron.inp', 'w') as f:
        f.write(f'{nwave}\n')

        for wavelength in wavelengths:
            f.write(f'{wavelength}\n')
    
    print('Wavelength grid:\n    wavelength_micron.inp\nDONE.\n')

    # Return the wavelength grid or useage also
    return wavelengths

# ------------------------------------------------------------ #
# Opacity-functions


# Create optool-kappa-file-script
def create_optoolscript(
        wavelength_path:str='../wavelength_micron.inp',
        phase=186,
        grainum_sizes:list=[0.1],
        grainsize_type:str='normal',
        grainsize_na:int=21,
        grainsize_power:float=-3.5,
        specie:str='mg2sio4',
        grain_type:str='mie',
        polarisation_matrix:str='n'
    ):
    """
    Writes script to optool to write opacity files for Radmc3d and dustopac.inp with all
    dust grain sizes.

    ARGUMENTS
      grainum_sizes: Can be a list or str!
              list : list of mean-grainsizes in micrometres
               str : path to grain_sizes_binned_PHASE.dat - file, as created by
                     a5d.bin_grainsizes()

      grainsize_na:int : number of grain sizes per sizes in grainum_sizes
                       : intermediate-size steps are used to reduce number resonances from
                       : singular grain sizes

         grainsize_type:str : available grain size distributions:
                   'normal' : normal Gaussian >> grainsize_power is not used
                'lognormal' : log-normal Gaussian >> grainsize_power is not used
                      'mrn' : power-law, MRN or Dohnanyi-style
                            : e.g. grainsize_power = -3.5 gives an MRN (or Dohnanyi) -powerlaw

      specie:str : available species:
       'mg2sio4' : Forsterite

      grain_type:str : available kind of grains:
               'mie' : Mie-theory based spherical compact grains
               'dhs' : Hollow spheres, good for polarisation
    
      polarisation_matrix:str : include computation of polarisation matrix or not?
                          'y' : -s is included in optool-script and dustopac.inp is adapted for this
                          'n' : default setting, no polarisation matrix, standard dustopac.inp
    """

    lnk_path = f'~/program/optool/my_lnk_data/{specie}.lnk'

    # Load wavelenth_grid from wavelength.inp-file and use for the wavelengt-settings!
    wavelengths,nwave = a3d.load_wavelengthgrid(
        path=wavelength_path
    )
    # Wavelength grid limits
    lmin = min(wavelengths)
    lmax = max(wavelengths)

    
    # Check type of grainum_sizes variable and act accordingly
    if type(grainum_sizes) == str:

        #Then load sizes
        if os.path.exists(grainum_sizes) == True:

            print(f'Extracting grain sizes from {grainum_sizes}')
            grainsizes,Nleafs = a3d.load_grainsizes(
                grainsize_path=grainum_sizes
            )

            # Extract an array with the grain sizes only
            grainum_sizes = np.unique(grainsizes[np.where(grainsizes > 0)[0]])
            # Change unit to um
            grainum_sizes *= 1e4

            print(f'    Grain sizes are (um):\n{grainum_sizes}')
    
    elif type(grainum_sizes) == list:
        print(f'    Grain sizes are (um):\n{grainum_sizes}')
    
    else:
        print('ERROR! grainum_sizes must be either list or string (even with only one grain size)')
        return 'STOPPING'


    # MRN-dist needs Amin, Amax, power and number of sizez
    #   amin, max: given from amean
    #   power and na needs to be given : 2 numbers
    #
    # Lognormal needs Amin, Amax, Amean, Asig (always the same?), and number of sizes
    #   amin, amax, asig, given from amean
    #   needs number of sizes : 1 number
    #
    # Normal needs the same as Lognormal


    # Define grain size limits
    Ngrainsizes = len(grainum_sizes)

    if Ngrainsizes > 1:

        # Smallest size's lower limit:
        #   if the half-distance to next grain size is smaller than the smallest size
        #   take smallest-size minus half-distance
        #   else take one tenth of the smallest size
        #
        midsize = 0.5 * (grainum_sizes[0]+grainum_sizes[1]) - grainum_sizes[0]
        if midsize < grainum_sizes[0]:
            amin_um = [grainum_sizes[0] - midsize]
            print('yes')
        else: 
            amin_um = [0.1*grainum_sizes[0]]
        # Old alternative
        # amin_um = [0.5*grainum_sizes[0]]

        amax_um = []

        for nn in range(Ngrainsizes - 1):

            # Size in between the current and next size
            midsize = 0.5 * (grainum_sizes[nn]+grainum_sizes[nn+1])

            # Minimum sizes
            amin_um.append(midsize)

            # Maximum sizes (previous max size is the next's min-size, ie, the same in loop)
            amax_um.append(midsize)
        
        # Final maximum size is then this ie add the previous half distance on the max grain size
        amax_um.append(grainum_sizes[-1] + grainum_sizes[-1] - amax_um[-1])

        # Check if we have lognormal distr
        if grainsize_type == 'normal' or 'lognormal':
            #
            # Define sigma:
            #
            # Sigma can be defined from amin-max and that FWHM = 2.355 sigma
            #
            # set FWHM = amax-amin
            # sigma = (amax-amin) / 2.355
            #
            # ie we want the HALF of the HEIGHT of the Gaussian to be "connected" 
            # to the nextgrain size Gaussian.
            # This can ofcourse be varied
            # With the max-min-sizes there won't be any overlaps
            #
            # Because of the possibility of having log-grain-sizes I should make 
            # a list of asig also.
            #
            asigma_um = []
            for nn in range(Ngrainsizes):
                asigma_um.append( (amax_um[nn] - amin_um[nn])/2.355 )
                
                # Test a wider spread
                #asigma_um.append( 2*(amax_um[nn] - amin_um[nn])/2.355 )

    else:
        # If there is only one size, plus minus half of it
        amin_um = [0.5*grainum_sizes[0]]
        amax_um = [grainum_sizes[0] + amin_um[0]]

        if grainsize_type == 'normal' or 'lognormal':
            # Then sigma is just
            asigma_um = [(amax_um[0] - amin_um[0])/2.355]


    # Some (temporary) output 
    #if grainsize_type == 'normal' or 'lognormal':
    #    print('    Grain size ranges are then')
    #    for nn in range(Ngrainsizes):
    #        print(f'{amin_um[nn]:.3e} - {grainum_sizes[nn]:.3e} - {amax_um[nn]:.3e}')


    with open(f'../optool_script_{phase}.sh','w') as f, \
         open(f'../dustopac_{specie}_{phase}.inp', 'w') as fopac:

        # Write header of optool script
        f.write(f'# Script to run optool to create kappa-files for R3D\n#\n')
        f.write('#    Grain sizes (um):')
        for size in grainum_sizes:
            f.write(f' {size:.3e} ')
        f.write('\n')
        f.write(f'#    Grain size distribution-style: {grainsize_type}\n')
        f.write(f'#    Grain type: {grain_type}\n')
        f.write(f'#    Dust specie: {specie}\n#\n')

        # Write header of dustopac_dust_phase.inp
        # dustopac:
        # 2
        # Number of species
        # -----------------------
        fopac.write(f'2\n{int(Ngrainsizes)}\n-----------------------------\n')

        # Loop over grain sizes
        for nn in range(Ngrainsizes):

            # Grain size limits
            amin = amin_um[nn]
            amax = amax_um[nn]
            amean = grainum_sizes[nn]

            # Write the script with optool commands
            #
            # TODO
            # For now I cant get the full scattering matrix to work in R3D:
            #
            # Warning: The scattering opacities kappa_s in the file dustkapscatmat_mg2sio4_8.966e-01.inp
            #          have relative differences of up to    8.9569185783100469E-003  with the 
            #          scattering matrix elements in the same file. I will correct kappa_s.
            #
            # And then:
            #    STOP 8800
            #
            # So I only write normal dustkappa-files for now
            # otherwise you add     -s      before      -radmc
            #
            if grainsize_type == 'normal' or 'lognormal':
                # For Gaussian or log-normal distributions
                asigma = asigma_um[nn]

                if grainsize_type == 'normal':
                    # normal distribution
                    # -a amin amax amean:-asig [na]
                    if polarisation_matrix == 'y':
                        f.write(f'optool -c {lnk_path} -{grain_type} -a {amin:.3e} {amax:.3e} {amean:.3e}:-{asigma:.3e} {grainsize_na} -lmin {lmin} -lmax {lmax} -nlam {nwave} -s -radmc mg2sio4_{amean:.3e}\n')
                    else:
                        f.write(f'optool -c {lnk_path} -{grain_type} -a {amin:.3e} {amax:.3e} {amean:.3e}:-{asigma:.3e} {grainsize_na} -lmin {lmin} -lmax {lmax} -nlam {nwave} -radmc mg2sio4_{amean:.3e}\n')

                if grainsize_type == 'lognormal':
                    # log-normal distribution
                    # -a amin amax amean:asig [na]
                    if polarisation_matrix == 'y':
                        f.write(f'optool -c {lnk_path} -{grain_type} -a {amin:.3e} {amax:.3e} {amean:.3e}:{asigma:.3e} {grainsize_na} -lmin {lmin} -lmax {lmax} -nlam {nwave} -s -radmc mg2sio4_{amean:.3e}\n')
                    else:
                        f.write(f'optool -c {lnk_path} -{grain_type} -a {amin:.3e} {amax:.3e} {amean:.3e}:{asigma:.3e} {grainsize_na} -lmin {lmin} -lmax {lmax} -nlam {nwave} -radmc mg2sio4_{amean:.3e}\n')
            
            if grainsize_type == 'mrn':
                # For Mie theory, 
                # MRN-distribution, 3.5 means a^-3.5, minus is dropped
                # -a AMIN AMAX APOW NA
                if polarisation_matrix == 'y':
                    f.write(f'optool -c {lnk_path} -{grain_type} -a {amin:.3e} {amax:.3e} {-grainsize_power} {grainsize_na} -lmin {lmin} -lmax {lmax} -nlam {nwave} -s -radmc mg2sio4_{grainum_sizes[nn]}\n')
                else:
                    f.write(f'optool -c {lnk_path} -{grain_type} -a {amin:.3e} {amax:.3e} {-grainsize_power} {grainsize_na} -lmin {lmin} -lmax {lmax} -nlam {nwave} -radmc mg2sio4_{grainum_sizes[nn]}\n')

            # Write dustopac-file for these opacities
            # dustopac.inp :
            # 1 for dustkappa-files, 10 for dustkapscatmat-files
            # 0
            # speciename
            # ---------------
            if polarisation_matrix == 'y':
                fopac.write(f"10\n0\n{specie}_{amean:.3e}\n-----------------------------\n")
            else:
                fopac.write(f"1\n0\n{specie}_{amean:.3e}\n-----------------------------\n")

    # Make into executable
    os.system(f'chmod +x ../optool_script_{phase}.sh')

    print(f'create_optoolscript:\n    ../optool_script_{phase}.sh\n    ../dustopac_{specie}_{phase}.inp\nDONE\n')



# Creates dust-kappa-files
def create_kappaabs(
        wavelengthpath:str='../wavelength_micron.inp',
        optconstlist:list=['mg2sio4'],
        agrainlist:list=[0.1],
        matdens:float=2
    ):

    # Number of "species" is
    nrspec = len(optconstlist)*len(agrainlist)

    # Load wavelength-grid
    wavelengths,nwave = a3d.load_wavelengthgrid(wavelengthpath)

    # Change units from um to cm and change list to np.array
    wavelengths = np.array([wavelength*1e-4 for wavelength in wavelengths])

    # Logarithmic spread parameter
    logawidth = max(agrainlist)/len(agrainlist)

    # Go through the different dust species
    for optconst in optconstlist:

        # Load dust species
        optconst_file = f'../bhmie/lnk/{optconst}.lnk'

        # Load grain size grid
        for agrainum in agrainlist:

            #Change units to cm
            agraincm = agrainum*1e-4

            # Create all opacity files
            opacity = compute_opac_mie(
                optconst_file=optconst_file,
                matdens=matdens,
                agraincm=agraincm,
                lamcm = wavelengths,
                logawidth=logawidth,
                wfact=3,
                na=21,
                extrapolate=True,
                verbose=True
            )

            # Example code and info from makedustopac and bhmie for reference here
            #compute_opac_mie(optconst_file,matdens,agraincm,lamcm,
            #                     theta=None,logawidth=None,wfact=3.0,na=20,
            #                     chopforward=0.0,errtol=0.01,verbose=False,
            #                     extrapolate=False):
            """
                matdens       = Material density in g/cm^3
                agraincm      = Grain radius in cm
                lamcm         = Wavelength grid in cm (a numpy array)
                theta         = Optional angular grid (a numpy array) between 0 and 180
                                which are the scattering angle sampling points at 
                                which the scattering phase function is computed.
                logawidth     = Optional: if set, the size agrain will instead be a 
                                sample of sizes around agrain. This helps to smooth out
                                the strong wiggles in the phase function and opacity
                                of spheres at an exact size. Since in Nature it rarely
                                happens that grains all have exactly the same size, this
                                is quite natural. The value of logawidth sets the width
                                of the Gauss in ln(agrain), so for logawidth<<1 this
                                give a real width of logawidth*agraincm. 
                wfact         = (default=3.0) Grid width of na sampling points in units
                                of logawidth. The Gauss distribution of grain sizes is 
                                cut off at agrain * exp(wfact*logawidth) and
                                agrain * exp(-wfact*logawidth).
                na            = (default=20) Number of size sampling points 
                                (if logawidth set).
                extrapolate   = If set to True, then if the wavelength grid lamcm goes
                                out of the range of the wavelength grid of the 
                                optical constants file, then it will make a suitable
                                extrapolation: keeping the optical constants constant
                                for lamcm < minimum, and extrapolating log-log for
                                lamcm > maximum.
            """

            # Write dustkappa_*.inp
            write_radmc3d_kappa_file(opacity,f'{optconst}_{agrainum}')

            # Move kappaabs files
            os.system(f'mv dustkappa_{optconst}_{agrainum}.inp ../')

    # Write dustopac.inp, list of all species and names
    with open('dustopac.inp','w') as f:

        # Write header
        f.write(f'2\n{nrspec}\n-----------------------------\n')

        # Write each species and grain size
        for optconst in optconstlist:
            for agrainum in agrainlist:
                f.write(f'1\n0\n{optconst}_{agrainum}\n-----------------------------\n')

    # Move opac-file to subfolder of simulation folder
    os.system(f'mv dustopac.inp ../')


# Load the star_opacity-file that is created with a5d.create_star
# These are c5d-data, Rosseland opacity in a R3D-grid.
def load_staropacities(
        path:str = '../star_opacities.dat'
    ):

    # load opacity.dat
    opacity = []
    with open(path, 'r') as fopacity:
        for line in fopacity.readlines():
            if line[0] != '#':
                opacity.append(float(line))

    # Change to np.array
    opacity = np.array(opacity)

    return opacity


# ------------------------------------------------------------ #
# Create simple R3D-data-functions

# Creates a spherical dust blob that imitates an AGB-star in the centrum of the grid
def create_duststar(
        Mstar:float = 1,
        Rstar:float = 100,
        Teff:float = 2700
    ):
    """
    Creates a clump of dust in the central parts of the grid to emulate an AGBstar.
    For now, it's a spherical symmetric clump with constant mass (bulk) density and
    temperature. Might be more adaptable in the future.

    Mstar is propto the luminosity of the clump but only up until the clump is
    optically thick. With Rstar=100Rsol, Teff=2700K, this happens at around 2e-6Msol,
    which gives around 450Lsol.

    INPUTS
    ------    
    Mstar: Mass of duststar clump
    Rstar: radius in solar radii
    Teff: temperature of dust clump (effective temperature) in Kelvin, default 2700K

    OUTPUTS
    -------
    dust_density_star.inp
    dust_temperature_star.dat
    dust_kappa_star.inp
    wavelength_micron.inp

    Jämför med svartkroppsstjärna. - 
    1Msun, 7000Lsun, 2700K

    Creates wavelength grid also
    IMPORTANT: all create stellar input-functions should simultaneously
    create the wavelength grid. The wavelengthfunc is more of a suplementary
    function.
    OR I could design this so that one has to run wavelengthfunc first before creating
    stellar models? And this just loads the wavelength grid and all other already made data?
    """

    # Change units of stellar props
    Mstar *= Msol # star mass in g
    Rstar *= Rsol # star radius in cm

    print(f'Radius of star: {Rstar/AUcm:.2} AU ({Rstar/Rsol} Rsun)')
    
    # Load griddistances and size of grid
    #           return nxyz,nrefines,nleafs,nbranch,gridedge
    nleafs = a3d.load_grid_properties()[2]
    griddistances = a3d.load_griddistances()
  
    # Compute average density of star (R3d uses g/cm3)
    stardensity = Mstar / (4/3 * np.pi * Rstar**3)
    print(f'Average star density: {stardensity:.2} g/cm3')

    # Create and fill density and temperature array
    # fill with average density (future: rho(R))
    # and with Teff (future: T(R))
    densities = np.zeros(nleafs)
    temperatures = np.zeros(nleafs)

    counter = 0
    for nn in range(nleafs):
        if griddistances[nn,0] <= Rstar:
            densities[nn] = stardensity
            temperatures[nn] = Teff
            counter += 1
    print(f'Cells inside star: {counter}')
    
    # Write dust_density
    print('Writing dust_density_star.inp')
    with open('../dust_density_star.inp', 'w') as f:

        # Header of dust_density (and dusttemperature)
        # 1
        # nleafs
        # number dust species
        f.write(f'1\n{int(nleafs)}\n1\n')
        for density in densities:
            f.write(f'{density}\n')

    # Write dust_temperature
    print('Writing dust_temperature_star.dat')
    with open('../dust_temperature_star.dat', 'w') as f:

        # Header of dust_temperature (and dusttemperature)
        # 1
        # nleafs
        # number dust species
        f.write(f'1\n{int(nleafs)}\n1\n')

        # Add temperatures
        for temperature in temperatures:
            f.write(f'{temperature}\n')

    print('Writing dustkappa_star.inp')
    # Constants in SI-units
    c = 2.998e8
    kb = 1.3806503e-23
    hplanck = 6.626068e-34

    # Peak of BB (in terms of frequency, but corresponding wavelength, in microns)
    wavelengthmax = hplanck*c / (2.821*kb*Teff) * 1e6

    # kappaabs_0 = 5000 cm2/g
    kappaabszero = 1
    #kappascatzero = 0

    # Create wavelength grid
    wavelengths = create_wavelength()
    nwave = len(wavelengths)

    # Declare arrays for "star"s kappa
    kappaabs_star = np.zeros(nwave)
    kappascat_star = np.zeros(nwave)
    # Note that I don't bother with any g[lambda] here, ie mean scattering angle
    # <cos(theta)> for the "star"'s kappafile.

    for nn in range(nwave):
        if wavelengths[nn] > wavelengthmax:
            kappaabs_star[nn] = kappaabszero * (wavelengthmax/wavelengths[nn])**1
            #kappascat_star[nn] = kappascatzero * (wavelengthmax/wavelengths[nn])**4
        else:
            kappaabs_star[nn] = kappaabszero
            #kappascat_star[nn] = kappascatzero
    
    # Chose scatter mode
    if kappascat_star[0] == 0:
        scattermode = 1
    else:
        scattermode = 2
    
    # Write dustkappa_star.inp
    with open('../dustkappa_star.inp', 'w') as f:

        # Write header
        # iformat: 1 - reads only lambda and abs
        #          2 - reads lambda, abs and scat
        #          3 - reads lambda, abs, scat and g
        # nlam: number of wavelength points
        f.write(f'{scattermode}\n{nwave}\n')

        # Write data
        for nn in range(nwave):
            f.write(f'{wavelengths[nn]}    {kappaabs_star[nn]}    {kappascat_star[nn]}\n')
    
    print('Duststar: done.\n')




# TODO
# Function to create proof-of-concept, sphere with rho ~ r^-2 and agrain ~ r^k (k>0, see Darwin-papers)

def create_spheredensity(

        optconstlist:list=['mg2sio4'],
        agrainlist:list=[0.1],
        totaldustmass:float=1.989e33,
        densitypower:float=-2,
        inradius:float=3,
        outradius:float=100
    ):
    """

    INPUTS
    ------
    total mass of dust cloud in gram: default: 1Msun
    radial density parameter: rho ~ r^k
    inner radius (in au)
    outer radius (in au)
    grainsizelist, agrainlist, list of all grainsizes in um
    optconstlist, list of optical constants-names
    """

    # Change units to cgs
    inradius *= AUcm
    outradius *= AUcm


    
    # Load grid distances (radial, x,y,z distances)
    griddistances = a3d.load_griddistances(
        gridpath='../r3dsims/grid_distances.csv',
        amrpath='../r3dsims/amr_grid.inp'
    )
    # check if outradius is smaller than larges radial griddistances, if not, set outradius to max griddistance
    if outradius > np.max(griddistances[:,0]):
        outradius = np.max(griddistances[:,0])
    
    # Load grid cell sizes
    cellsizes = a3d.load_cellsizes(
        sizepath='../r3dsims/grid_cellsizes.csv',
        amrpath='../r3dsims/amr_grid.inp'
    )
    # Note: both returns np.arrays
    nleafs = np.size(cellsizes)

    # Number of grainsizes
    Nagrain = len(agrainlist)

    # Number of dust species
    nrspec = len(optconstlist)*Nagrain

    # Minimum grain size
    agrainmin = min(agrainlist)



    # Compute normalization density (ie zero density)
    # TODO this is the same for all dust species, allow for different morphs? Ie different
    # radii and radial dependence later?
    # TODO this divides the zero density equally between dust species and grain sizes
    # The r^(2+inputvalue) is due to spherical coordinates!
    radiusintegral = s.integrate.quad(lambda x: x**(2+densitypower), inradius, outradius)
    zerodensity = totaldustmass / nrspec * inradius**densitypower / (4.*np.pi * radiusintegral[0])


    # Create density distribution
    # ns, species, choses grain sizes and chemical compositions

    print('Writing dust_density.inp')

    totaldustmass = 0
    densitymatrix = np.zeros(nrspec*nleafs)

    #setted_list = [2,9,6,20,15]
    #value_chosen = 17
    #min(setted_list, key=lambda x:abs(x-value_chosen)) > gives nearest value

    # settedlist: agrainlist
    # valuechosen: agrain(R) [something simple for now]: min(agrainlist) * R/inradius

    # Add densities
    for ns in range(len(optconstlist)):
        for nn,griddistance in enumerate(griddistances[:,0]):
            if griddistance >= inradius and griddistance <= outradius:
                
                # Find which grainsize is closest to radial grainsize relation
                grainsize = min(agrainlist, key=lambda x:abs(x-(agrainmin*griddistance/inradius)))
                ngrain = np.where(np.array(agrainlist) == grainsize)[0][0]

                print(agrainmin*griddistance/inradius,grainsize,ngrain)

                # Allocate densities to bins of grain sizes, and for each grain chemical specie
                densitymatrix[nn + ngrain*nleafs + ns*nleafs*Nagrain] = zerodensity * (griddistance/inradius)**densitypower


                # Compute real total dust mass
                totaldustmass += densitymatrix[nn + ngrain*nleafs + ns*nleafs*Nagrain]*cellsizes[nn]**3



    # Open density file
    with open('../dust_density.inp','w') as f:

        # Write header
        f.write(f'1\n{nleafs}\n{nrspec}\n')

        # Add densitymatrixloop here

    print('Finished dust_density.inp')
    print(f'Total dust mass is {totaldustmass} g ({totaldustmass/Msol} Msol)')

    


    # For each radial distance, round grainsize to nearest 1/10 of max grainsize.
    # No, round to the nearest value in the grainsizelist that should be inputed also


    # Output also the REAL total dust mass with the help of cellsizes array since
    # the sphere might well be cut in the outer parts of the grid
    # also because the sphere is in a cubic grid



    print('create_spheredensity: Done')

    # TODO create opacity files directly here?

# ------------------------------------------------------------ #
# Functions that combines different dust-specie files


# Function that merges density.inp-files into one
def merge_dustdensities(
        workpath:str = '../r3dresults/',
        filenames:list=['dust_density.inp'],
    ):
    """
    Merges listed dust_density-files in a workfolder into one inp-file

    ARGUMENTS
      workfolder:str = path to specific work/project folder
      filenames:list = list of dust_density-filenames
    
    RETURNS
      One dust_density.inp for each phase in phases-list moved to work folder
    
    INFO
      Structure of dust_density is as follows (and line numbers)

    0 1 (yes, the number 1, nothing more or less)
    1 number of leafs (or cells)
    2 number dust species
    3 onward-> list of densities in g/cm3 starting from cell 0,0,0
      and in order of dust specie. Must be the same order as in dustopac.inp
    """

    # Automatically add / to end of path if it's missing
    if workpath[-1] != '/':
        workpath += '/'

    # Number of density files to merge
    Nfiles = len(filenames)

    if Nfiles < 2:
        print(f'Number of density files is {Nfiles}. There is nothing to merge')
    else:
        print(f'Merging {Nfiles} dust_density_*.inp files in {workpath}')
        
        densitylist = []
        nleafslist = []
        Nspeclist = []

        for filename in filenames:
            with open(f'{workpath}{filename}', 'r') as f_read:

                # read lines:
                # skip line 0
                #  
                # lines 3 to nleafs+3 are densities for spec1
                # lines nspecie*(nleafs+3) to (nspecie+1)*(nleafs+3) are for species n

                for nl,line in enumerate(f_read.readlines()):
                    if nl == 1:
                        # Line 1 is nleafs, number of grid cells
                        nleafslist.append(int(line))
                    if nl == 2:
                        # Line 2 is Nspecies, number of species within density file
                        Nspeclist.append(int(line))
                    if nl > 2:
                        # The rest are densities, one specie after another
                        # IE, the order IS IMPORTANT
                        densitylist.append(float(line))

        # Total number of species is the sum of them then
        Nspecies = sum(Nspeclist)

        # Check that number of cells (nleafs) doesn't vary with density-file
        if nleafslist[0] != sum(nleafslist)/len(filenames):
            print(f'Number of cells (nleafs) varies:\n    {nleafslist}')
            nleafs = 0

        else:
            # Number of leafs are correct then continue here
            nleafs = nleafslist[0]

            # Write the new total dust_density-file:
            # Open new dust_density:
            with open(f'../dust_density_merged.inp','w') as f_dust:

                # Write header
                # 1
                # nleafs
                # number dust species
                f_dust.write(f'1\n{int(nleafs)}\n{int(Nspecies)}\n')

                for density in densitylist:
                    f_dust.write(f'{density}\n')

        # Mode the new dust_density to its folder
        os.system(f'mv ../dust_density_merged.inp {workpath}dust_density_merged.inp')

        # Print final message
        print(f'merge_dustdensities:\n    {workpath}dust_density_merged.inp\nDONE\n')


# Function for merging a list of dustopac.inp files
# TODO rewrite so that it works with a specific path instead! more adaptable
def merge_dustopac(
        filenames:list=['dustopac.inp'],
        modelname:str='st28gm06n052',
        phases:list=[186,190,198],
    ):
    """
    Given a model, phase, and list of dustopac*.inp-files this
    will merge them into one, species in the same order as in
    the list of filenames.

    ARGUMENTS
      filenames:list = list of dustopac-filenames
      modelname:str = co5bold-model-name (folder in r3dresults)
      phases:list = list of phase-designations to loop through (folders inside modelfolder)
    
    RETURNS
      one dustopac.inp for each phase in phases-list, moved to each phase-folder

    INFO:
      structure (and line numbers) of dustopac.inp is the following

    0 iformat (set to 2)
    1 nspecies
    2 -----------------------------
    3 inputstyle[1] (set to 1 since I use dustkappa_NAME.inp)
    4 0
    5 <name of dust species 1>
    6 -----------------------------
    7 inputstyle[2]
    8 0
    9 <name of dust species 2>
    10-----------------------------
    etc
    """
    # Path to r3dresults-folder
    path = '../r3dresults/'

    # Number of density files to merge
    Nfiles = len(filenames)

    for phase in phases:

        if Nfiles < 2:
            print(f'Number of dustopac files (phase: {phase}) is {Nfiles}. There is nothing to merge')

        else:
            print(f'Merging {Nfiles} dustopac_*.inp files for phase {phase}')

            specienames = []            
            Nspeclist = []

            # First extract number of species in the files
            for filename in filenames:
                with open(f'{path}{modelname}/{phase}/{filename}', 'r') as f_read:
                    for nl,line in enumerate(f_read.readlines()):
                        if nl == 1:
                            # Line 1 is nspecies, number of species included here
                            Nspeclist.append(int(line))

            # Total number of species
            Nspecies = sum(Nspeclist)

            # Then extract the specie-names
            # These are on lines 5, 9, 13
            # ie, 1+4*(nspecie+1)
            for nf,filename in enumerate(filenames):
                with open(f'{path}{modelname}/{phase}/{filename}', 'r') as f_read:
                    for nl,line in enumerate(f_read.readlines()):
                        for nspecie in range(Nspeclist[nf]):
                            if nl == 1+4*(nspecie+1):
                                specienames.append(line.rstrip())

            # Write new dust_opac:
            with open(f'../dustopac_{phase}.inp','w') as f_opac:
                # Header
                f_opac.write(f'2\n{Nspecies}\n-----------------------------\n')
                # Then each specie
                for speciename in specienames:
                    f_opac.write(f'1\n0\n{speciename}\n-----------------------------\n')

            # Move new file to final folder
            os.system(f'mv ../dustopac_{phase}.inp {path}{modelname}/{phase}/dustopac.inp')

            # Print final message
            print(f'merge_dustopac:\n    {path}{modelname}/{phase}/dustopac.inp\nDONE\n')


# Function that merges dust_temperature-files in a given folder
def merge_dusttemperatures(
        filenames:list=['dust_temperature.dat'],
        modelname:str='st28gm06n052',
        phases:list=[186,190,198],
    ):
    """
    Given a model, phase, and list of dust_temperature*.dat-files this
    will merge them into one, species in the same order as in
    the list of filenames.

    ARGUMENTS
      filenames:list = list of dust_temperature-filenames
      modelname:str = co5bold-model-name (folder in r3dresults)
      phases:list = list of phase-designations to loop through (folders inside modelfolder)
    
    RETURNS
      one dust_temperature.dat for each phase in phases-list, moved to each phase-folder

    INFO:
      structure (and line numbers) of dust_temperature.dat is the following (with line
      numbers)

    0 1 (yes the number 1, nothin less, nothing more)
    1 number of leafs (ie number of cells)
    2 number dust species
    3 and onwards, temperatures of cells in K, starting at 0,0,0, listed in same order
      as in filenames-list.
    """
    '''
            workpath:str = '../r3dresults/',
            filenames:list=['dust_density.inp'],
        ):
        """
        Merges listed dust_density-files in a workfolder into one inp-file

        ARGUMENTS
        workfolder:str = path to specific work/project folder
        filenames:list = list of dust_density-filenames
        
        RETURNS
        One dust_density.inp for each phase in phases-list moved to work folder
        
        INFO
        Structure of dust_density is as follows (and line numbers)

        0 1 (yes, the number 1, nothing more or less)
        1 number of leafs (or cells)
        2 number dust species
        3 onward-> list of densities in g/cm3 starting from cell 0,0,0
        and in order of dust specie. Must be the same order as in dustopac.inp
        """
    '''




    # Path to r3dresults-folder
    path = '../r3dresults/'

    # Number of density files to merge
    Nfiles = len(filenames)

    for phase in phases:

        if Nfiles < 2:
            print(f'Number of dustopac files (phase: {phase}) is {Nfiles}. There is nothing to merge')

        else:
            print(f'Merging {Nfiles} dust_temperature_*.dat files for phase {phase}')

            temperaturelist = []
            nleafslist = []
            Nspeclist = []

            # Open each file in filenames
            for filename in filenames:
                with open(f'{path}{modelname}/{phase}/{filename}', 'r') as f_read:
                    #
                    # read lines:
                    # skip line 0
                    #  
                    # lines 3 to nleafs+3 are temperatures for spec1
                    # lines nspecie*(nleafs+3) to (nspecie+1)*(nleafs+3) are for species n
                    #
                    for nl,line in enumerate(f_read.readlines()):
                        if nl == 1:
                            # Line 1 is nleafs, number of grid cells
                            nleafslist.append(int(line))
                        if nl == 2:
                            # Line 2 is Nspecies, number of species within temperature file
                            Nspeclist.append(int(line))
                        if nl > 2:
                            # The rest are densities, one specie after another
                            # IE, the order IS IMPORTANT
                            temperaturelist.append(float(line))

            # Total number of species is the sum of them then
            Nspecies = sum(Nspeclist)

            # Check that number of cells (nleafs) doesn't vary with density-file
            if nleafslist[0] != sum(nleafslist)/len(filenames):
                print(f'Number of cells (nleafs) varies:\n    {nleafslist}')
                nleafs = 0

            else:
                # Number of leafs are correct then continue here
                nleafs = nleafslist[0]

                # Write the new total dust_density-file:
                # Open new dust_temperature:
                with open(f'../dust_temperature_{phase}.dat','w') as f_temperature:

                    # Write header
                    # 1
                    # nleafs
                    # number dust species
                    f_temperature.write(f'1\n{int(nleafs)}\n{int(Nspecies)}\n')

                    for temperature in temperaturelist:
                        f_temperature.write(f'{temperature}\n')

        # Mode the new dust_density to its folder
        os.system(f'mv ../dust_temperature_{phase}.dat {path}{modelname}/{phase}/dust_temperature.dat')

        # Print final message
        print(f'merge_dustdensities:\n    {path}{modelname}/{phase}/dust_temperature.dat\nDONE\n')
