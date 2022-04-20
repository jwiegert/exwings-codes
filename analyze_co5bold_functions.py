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

# Cython info: might be needed later
# @cython.cfunc
# this decorator when declaring what's included in the functions
# @cython.locals(a=cython.int)
# AUcm = cython.declare(cython.float ,1.49598e13) # cm



# ============================================================


# Load c5d grid properties
@cython.cfunc
def load_grid_properties(
        savpath:str='../co5bold_data/dst28gm06n056/st28gm06n056_140.sav'
    ):
    """
    TODO INFO
    outputorder: c5dgrid,cellcourners,cellsize
    """
    
    # Load sav-file
    c5ddata = readsav(savpath)

    # Extract data - TODO MAY HAVE TO CHOSE THIS MANUALLY DEPENDING ON FILE
    c5ddata = c5ddata['ful']

    # Get number of gridcells from co5bold data
    nc5dedge = cython.declare(cython.int, np.size(c5ddata['Z'][0][0][16]))

    # Declare output arrays
    c5dgrid = np.zeros((nc5dedge,3))
    cellcourners = np.zeros((nc5dedge,3))
    cellsizesx = []
    cellsizesy = []
    cellsizesz = []

    # Declare variables
    nn = cython.declare(cython.int)

    for nn in range(nc5dedge):

        # Save x,y,z distances in cm
        c5dgrid[nn,0] = c5ddata['Z'][0][0][16][0][0][nn]
        c5dgrid[nn,1] = c5ddata['Z'][0][0][19][0][nn][0]
        c5dgrid[nn,2] = c5ddata['Z'][0][0][22][nn][0][0]

        # Save coordinates of cell courners
        cellcourners[nn,0] = c5ddata['Z'][0][0][25][0][0][nn]
        cellcourners[nn,1] = c5ddata['Z'][0][0][28][0][nn][0]
        cellcourners[nn,2] = c5ddata['Z'][0][0][31][nn][0][0]


        # Extract cell sizes
        if nn > 0:
            cellsizesx.append(
                (cellcourners[nn,0] - cellcourners[nn-1,0])
            )
            cellsizesy.append(
                (cellcourners[nn,1] - cellcourners[nn-1,1])
            )
            cellsizesz.append(
                (cellcourners[nn,2] - cellcourners[nn-1,2])
            )

    # Add final grid courner
    cellcourners[-1,0] = c5ddata['Z'][0][0][25][0][0][-1]
    cellcourners[-1,1] = c5ddata['Z'][0][0][28][0][-1][0]
    cellcourners[-1,2] = c5ddata['Z'][0][0][31][-1][0][0]

    cellsizesx.append(
        (cellcourners[-1,0] - cellcourners[-2,0])
    )
    cellsizesy.append(
        (cellcourners[-1,1] - cellcourners[-2,1])
    )
    cellsizesz.append(
        (cellcourners[-1,2] - cellcourners[-2,2])
    )

    # Extract minimum grid size
    cellsize = (min(cellsizesx) + min(cellsizesy) + min(cellsizesz))/3

    return c5dgrid,cellcourners,cellsize


# Extract co5bold densities into a separate array 
# - note this is probably faster than loading it in the r3d-file-writing functions
@cython.cfunc
def load_star_properties(
        savpath:str='../co5bold_data/dst28gm06n056/st28gm06n056_140.sav'
    ):
    """
    TODO info
    """
    
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
    c5dstar_temperatures = np.zeros((nc5dedge,nc5dedge,nc5dedge))

    # Extract densities - This can take time, some 2min per property
    for nx in range(nc5dedge):
        for ny in range(nc5dedge):
            for nz in range(nc5dedge):
                c5dstar_densities[nx,ny,nz] = c5ddata['Z'][0][0][34][nx][ny][nz]
                c5dstar_temperatures[nx,ny,nz] = c5ddata['EOS'][0][0][1][nx][ny][nz]
    
    return c5dstar_densities, c5dstar_temperatures


# Extract and create input data with
# TODO duststar opacities



# Function for loading one dust specie from c5d-data
@cython.cfunc
@cython.locals(nspecies=cython.int)
def load_dust_density(
        savpath:str='../co5bold_data/dst28gm06n052/st28gm06n052_186.sav',
        nspecies:int=0
    ):
    """
    TODO instructions
    
    """

    # Load sav-file
    c5ddata = readsav(savpath)

    # Extract data
    c5ddata = c5ddata['ful']

    # Get number of gridcells from co5bold data
    nc5dedge = cython.declare(cython.int, np.size(c5ddata['Z'][0][0][16]))

    # Declare variables
    nx = cython.declare(cython.int)
    ny = cython.declare(cython.int)
    nz = cython.declare(cython.int)

    # Declare np.arrays, tmeperature is a placeholder here TODO
    c5ddust_densities = np.zeros((nc5dedge,nc5dedge,nc5dedge))
    c5ddust_temperatures = 0# np.zeros((nc5dedge,nc5dedge,nc5dedge))

    # Extract densities - This can take some time
    for nx in range(nc5dedge):
        for ny in range(nc5dedge):
            for nz in range(nc5dedge):
                c5ddust_densities[nx,ny,nz] = c5ddata['Z'][0][0][40+3*nspecies][nx][ny][nz]

                # TODO: find the dust temperature in the c5d-data
                #c5ddust_temperatures[nx,ny,nz] = c5ddata['EOS'][0][0][1][nx][ny][nz]
    
    return c5ddust_densities, c5ddust_temperatures






def create_star(
        # All necessary inputs
        # paths!
        savpath:str='../co5bold_data/dst28gm06n056/st28gm06n056_140.sav',
        amrpath:str='../amr_grid.inp',
        gridpath:str='../grid_distances.csv',
        sizepath:str='../grid_cellsizes.csv'
    ):
    """
    Extracts data from Co5bold sav-files and creates a dust_star for Radmc3d.

    INPUT
    -----
    amrpath: path to amr_grid.inp
    gridpath: path to grid_distances.csv
    sizepath: path to grid_cellsizes.csv
    savpath: path to co5bold sav-file

    OUTPUT
    ------
    R3D density file: dust_density_star.inp
    R3D temperature file: dust_temperature_star.dat
    TODO: star's  opacities
    """

    # Load R3D grid
    print('Loading R3D grid')
    nleafs = a3d.load_gridprops(amrpath=amrpath)[2]
    r3ddistances = a3d.load_griddistances(amrpath=amrpath,gridpath=gridpath)
    r3dcellsizes = a3d.load_cellsizes(amrpath=amrpath,sizepath=sizepath)

    # Load C5D grid
    print('Loading C5D grid properties')
    c5dgrid, c5dcellcourners, c5dcellsize = load_grid_properties(savpath=savpath)

    # Check so that the smallest c5dcells are not larger than the r3d's smallest cells
    if r3dcellsizes.min() <= c5dcellsize:
        print('\nERROR')
        print('    R3D grid resolution is higher than C5D grid, stopping')
        print('    No output is given. Change your R3D grid cells to something larger.\n')

    else:
        # Load C5D densities

        # TODO: instead of loading data in arrays
        # just load the specific sav-file and then loop in that array
        # then I only need one loop instead of one per data plus
        # these that loads it into the data here below?
        # or that could be slower since now I work with nparrays...
        print('Loading C5D star properties (density, temperature, opacity)')
        c5dstar_densities,c5dstar_temperatures = load_star_properties(savpath=savpath)

        # Start working :)
        print('Translating C5D data to R3D data')

        # Declare stuff for the loops
        r3d_densities = 0
        r3d_temperatures = 0
        progbar = 0

        # Adaptive range used for when cellsizes are too similar - fixes bug 
        # where I got a lot of zero-cells because my loop misses the c5dcells inside
        # the current r3dcell
        adaptive_range = c5dcellsize/r3dcellsizes.min() * 1.1

        # Open r3d data files
        with open('../dust_density_star.inp', 'w') as fdensity, open('../dust_temperature_star.dat', 'w') as ftemperature:

            # Write headers:
            # 1
            # nleafs
            # number dust species
            fdensity.write(f'1\n{int(nleafs)}\n1\n')
            ftemperature.write(f'1\n{int(nleafs)}\n1\n')

            # Loop over r3d grid
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
                c5dxrange = np.argwhere(r3dxrange[0] <= c5dgrid[np.argwhere(c5dgrid[:,0] < r3dxrange[1]),0])[:,0]
                c5dyrange = np.argwhere(r3dyrange[0] <= c5dgrid[np.argwhere(c5dgrid[:,1] < r3dyrange[1]),1])[:,0]
                c5dzrange = np.argwhere(r3dzrange[0] <= c5dgrid[np.argwhere(c5dgrid[:,2] < r3dzrange[1]),2])[:,0]

                # Number of c5dcells within r3dcell
                nchildcells = c5dxrange.size*c5dyrange.size*c5dzrange.size

                # Check so that there are any c5dcells within r3dcell
                if nchildcells < 1:

                    # When r3dcellsize approx c5dcellsize the loop misses the c5dcells inside r3dcells
                    # due to numerical errors. So here I increase the spatial range to loop.
                    temprange = [r3dxrange[0]*adaptive_range,r3dxrange[1]*adaptive_range]
                    c5dxrange = np.argwhere(temprange[0] <= c5dgrid[np.argwhere(c5dgrid[:,0] < temprange[1]),0])[:,0]
                    temprange = [r3dyrange[0]*adaptive_range,r3dyrange[1]*adaptive_range]
                    c5dyrange = np.argwhere(temprange[0] <= c5dgrid[np.argwhere(c5dgrid[:,1] < temprange[1]),1])[:,0]
                    temprange = [r3dzrange[0]*adaptive_range,r3dzrange[1]*adaptive_range]
                    c5dzrange = np.argwhere(temprange[0] <= c5dgrid[np.argwhere(c5dgrid[:,2] < temprange[1]),2])[:,0]

                    # New number of c5dcells within r3dcell
                    nchildcells = c5dxrange.size*c5dyrange.size*c5dzrange.size

                # Then loop through c5dcells within r3dcell
                for nnz in c5dzrange:
                    for nny in c5dyrange:
                        for nnx in c5dxrange:

                            # Sum all densities
                            r3d_densities += c5dstar_densities[nnx,nny,nnz]

                            # Sum all temperatures
                            r3d_temperatures += c5dstar_temperatures[nnx,nny,nnz]

                # Average the data of each r3dcell by number of c5dcells
                r3d_densities /= nchildcells
                r3d_temperatures /= nchildcells

                # Write data to r3d files
                fdensity.write(f'{r3d_densities}\n')
                ftemperature.write(f'{r3d_temperatures}\n')

                # Reset data
                r3d_densities = 0
                r3d_temperatures = 0

                # Some progress bar info
                if int(nr3d/nleafs*100) == 25 and progbar == 0:
                    progbar += 1
                    print('Finished 25 per cent of the grid.')

                if int(nr3d/nleafs*100) == 50 and progbar == 1:
                    progbar += 1
                    print('Finished 50 per cent of the grid.')

                if int(nr3d/nleafs*100) == 75 and progbar == 2:
                    progbar += 1
                    print('Finished 75 per cent of the grid.')

    print('C5D Dust-star: done.\n')






# Extract and construct dust_density-files from C5D-data
# uses ['Z'][0][0][40][x][y][z] and ['Z'][0][0][43][x][y][z]
# or rather use
# ['Z'][0][0][40+3*nspecie][x][y][z]
# where nspecie is 0, 1, 2 ... up until (number of dust species)-1
# and the data are in number densities cm^-3
# Mass density in each cell = 4/3 pi agrain**3 * rhograin * numberdensity
# So another input based on what I used for dustkappa for each specie
# Also:
# str(teststar['Z'][0][0][40+3*nspecie+2])[4:-1]
# gives a string with the name of the specie!
# so this can also print the dustkappa-list-file for r3d!

def create_dust_files(
        savpath:str='../co5bold_data/dst28gm06n052/st28gm06n052_186.sav',
        amrpath:str='../amr_grid.inp',
        gridpath:str='../grid_distances.csv',
        sizepath:str='../grid_cellsizes.csv',
        Nspecies:int=1,
        grainsizeum:list=[0.1],
        graindensity:list=[2]
    ):
    """
    grainsizecm & graindensity should be lists, one for each specie
    """

    # Load R3D grid
    print('Loading R3D grid')
    nleafs = a3d.load_gridprops(amrpath=amrpath)[2]
    r3ddistances = a3d.load_griddistances(amrpath=amrpath,gridpath=gridpath)
    r3dcellsizes = a3d.load_cellsizes(amrpath=amrpath,sizepath=sizepath)

    # Load C5D grid
    print('Loading C5D grid properties')
    c5dgrid, c5dcellcourners, c5dcellsize = load_grid_properties(savpath=savpath)

    # Load C5D data
    c5ddata = readsav(savpath)
    c5ddata = c5ddata['ful']

    # Number of dust species in data:
    Nc5dspecies = int((len(c5ddata['Z'][0][0]) - 40)/3)

    # Check so that the smallest c5dcells are not larger than the r3d's smallest cells
    if r3dcellsizes.min() <= c5dcellsize:
        print('\nERROR')
        print('    R3D grid resolution is higher than C5D grid, stopping')
        print('    No output is given. Change your R3D grid cells to something larger.\n')
    
    # Check so that the number of dust species declared is not larger than available
    elif Nspecies > Nc5dspecies:
        print('\nERROR')
        print('    You asked for more dust species than available in C5D-data')
        print(f'    Number of dust species in C5D-data: {Nc5dspecies}')

    else:
        # Adaptive range used for when cellsizes are too similar - fixes bug 
        # where I got a lot of zero-cells because my loop misses the c5dcells inside
        # the current r3dcell
        adaptive_range = c5dcellsize/r3dcellsizes.min() * 1.1

        # Open r3d data files
        with open('../dust_density_dust.inp', 'w') as fdensity, open('../dustopac_dust.inp', 'w') as fopac:

            # Write headers:
            # Density:
            # 1
            # nleafs
            # number dust species
            fdensity.write(f'1\n{int(nleafs)}\n{int(Nspecies)}\n')

            # dustopac:
            # 2
            # Number of species
            # -----------------------
            fopac.write(f'2\n{int(Nspecies)}\n-----------------------------\n')

            # Loop through the number of species you want to include
            for nspecies in range(Nspecies):

                # Declare stuff for the loops
                r3d_density = 0
                progbar = 0

                # Output species-number, species name, grain size and grain mass density
                speciesname = str(c5ddata['Z'][0][0][42+3*nspecies])[4:-1]
                grainsizecm = grainsizeum[nspecies]*1e-4

                print(f'Writing dust specie number {nspecies+1}:')
                print(f'    {speciesname}')
                print(f'    Grain size: {grainsizecm*1e4} um')
                print(f'    Grain mass density: {graindensity[nspecies]} g cm-3')

                # Load c5d-dust densities and temperatures
                #c5ddensities, c5dtemperatures = load_dust_density(savpath=savpath, nspecies=nspecies)

                # Write the dustopac file
                # 1
                # 0
                # name
                # ---------------
                fopac.write(f"1\n0\n{speciesname}\n-----------------------------\n")

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
                    c5dxrange = np.argwhere(r3dxrange[0] <= c5dgrid[np.argwhere(c5dgrid[:,0] < r3dxrange[1]),0])[:,0]
                    c5dyrange = np.argwhere(r3dyrange[0] <= c5dgrid[np.argwhere(c5dgrid[:,1] < r3dyrange[1]),1])[:,0]
                    c5dzrange = np.argwhere(r3dzrange[0] <= c5dgrid[np.argwhere(c5dgrid[:,2] < r3dzrange[1]),2])[:,0]

                    # Number of c5dcells within r3dcell
                    nchildcells = c5dxrange.size*c5dyrange.size*c5dzrange.size

                    # Check so that there are any c5dcells within r3dcell
                    if nchildcells < 1:

                        # When r3dcellsize approx c5dcellsize the loop misses the c5dcells inside r3dcells
                        # due to numerical errors. So here I increase the spatial range to loop.
                        temprange = [r3dxrange[0]*adaptive_range,r3dxrange[1]*adaptive_range]
                        c5dxrange = np.argwhere(temprange[0] <= c5dgrid[np.argwhere(c5dgrid[:,0] < temprange[1]),0])[:,0]
                        temprange = [r3dyrange[0]*adaptive_range,r3dyrange[1]*adaptive_range]
                        c5dyrange = np.argwhere(temprange[0] <= c5dgrid[np.argwhere(c5dgrid[:,1] < temprange[1]),1])[:,0]
                        temprange = [r3dzrange[0]*adaptive_range,r3dzrange[1]*adaptive_range]
                        c5dzrange = np.argwhere(temprange[0] <= c5dgrid[np.argwhere(c5dgrid[:,2] < temprange[1]),2])[:,0]

                        # New number of c5dcells within r3dcell
                        nchildcells = c5dxrange.size*c5dyrange.size*c5dzrange.size

                    # Then loop through c5dcells within r3dcell
                    for nnz in c5dzrange:
                        for nny in c5dyrange:
                            for nnx in c5dxrange:
                                # Sum all densities
                                #r3d_density += c5ddensities[nnx,nny,nnz]
                                r3d_density += c5ddata['Z'][0][0][40+3*nspecies][nnx][nny][nnz]

                    # Average the density of each r3dcell by number of c5dcells
                    # and
                    # Recalculate number density to mass density
                    # 4/3 pi agrain**3 * rhograin * numberdensity
                    # Units are in cgs.
                    r3d_density *= \
                        4.1887902047863905 * grainsizecm**3 * graindensity[nspecies] / nchildcells
                    # TODO: solve why my ass is crazy high...

                    # Write data to r3d files
                    fdensity.write(f'{r3d_density}\n')

                    # Reset data
                    r3d_density = 0

                    # Some progress bar info
                    if int(nr3d/nleafs*100) == 25 and progbar == 0:
                        progbar += 1
                        print('Finished 25 per cent of the grid.')

                    if int(nr3d/nleafs*100) == 50 and progbar == 1:
                        progbar += 1
                        print('Finished 50 per cent of the grid.')

                    if int(nr3d/nleafs*100) == 75 and progbar == 2:
                        progbar += 1
                        print('Finished 75 per cent of the grid.')

                    if int(nr3d/nleafs*100) > 99 and progbar == 3:
                        print('Finished 99 per cent of the grid.')

    print('C5D Dust-data: done.\n')
