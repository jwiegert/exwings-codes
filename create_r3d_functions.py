# Functions for creating input data for R3d simulations for Exwings project
# ------------------------------------------------------------ #
# Various useful packages
import sys
import numpy as np
import scipy as s
from scipy.integrate import quad
from datetime import datetime

# Might be used later
#import os
#import scipy as s
#from scipy.integrate import quad

# My own functions
import analyze_r3d_functions as a3d

# Basic definitions
AUcm = 1.49598e13 # cm
Msol = 1.989e33 # g
Rsol = 6.955e10 # cm

# ------------------------------------------------------------ #
# Shorter simpler functions

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

# TODO
# function that combines different dust-specie density files into one
# def merge_dustdensity():
    # headre of dust_density
    # 1
    # nleafs
    # number dust species

# TODO
# function that combines different dust tmeperature files
# def merge_dusttemperature():


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
        print(f'Warning, number of base cells is not even, it is now {nxyz}.')
    
    # Create basic parameters of the grid
    #   nbasecubes : total number of base cells
    #     gridedge : total size of the grid side
    # gridcourners : coordinates of base grid courners
    #     griddist : list of distances to center of grid (not for R3D)
    nbasecubes     = int(nxyz * nxyz * nxyz)
    gridedge       = nxyz * basecubesize
    gridcourners   = np.linspace(-gridedge*0.5,gridedge*0.5,nxyz+1)
    print(f'Length of total side of whole grid: {gridedge/AUcm:.2f} AU')

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
    
    # Print grid_distances.csv and grid_cellsizes.csv

    if savegrid == 'y' or savegrid == 'yes' or savegrid == 'Y':
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
        print('Finished grid_distances.csv')

        # Save sizes of grid cells
        # This also has the same order as the dust_densities files
        with open('../grid_cellsizes.csv', 'w') as f:
            f.write('# Sizes of grid cells (same order as in dust_densities.inp) in cm\n')
            f.write('# ----------\n')
            for gridsize in gridsizes:
                f.write(f'{gridsize}\n')
        print('Finished grid_cellsizes.csv')
    
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

    print('Finished grid_info.txt')

    # Finish function
    print('Create grid: done.\n')

# ------------------------------------------------------------ #

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
    
    print('Wavelength grid: done.\n')

    # Return the wavelength grid or useage also
    return wavelengths

# ------------------------------------------------------------ #

# Creates a dust blog that imitates an AGB-star in the centrum of the grid
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
    nleafs = a3d.load_gridprops()[2]
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

    print('Finished dust_density_star.inp')

    # Write dust_temperature
    print('Writing dust_temperature_star.dat') # TODO: check the correct filenames!
    with open('../dust_temperature_star.dat', 'w') as f:

        # Header of dust_temperature (and dusttemperature)
        # 1
        # nleafs
        # number dust species
        f.write(f'1\n{int(nleafs)}\n1\n')

        # Add temperatures
        for temperature in temperatures:
            f.write(f'{temperature}\n')
    print('Finished dust_temperature_star.dat')

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
    # TODO check the math here, why r^(2+inputvalue)? Is it due to some derivative? No, due to spherical coordinates!
    # TODO this divides the zero density equally between dust species and grain sizes
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
