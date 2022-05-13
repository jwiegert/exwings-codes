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
    Loads c5d-star's densities, temperatures, and opacities, and puts them into 3D-arrays

    INPUT
    savpath: path to sav-file
    nbins: number of bins to put the opacities in, ie number of duststar-species for the star

    OUTPUT
    c5dstar_densities: 3D array with star's densities
    c5dstar_temperatures: 3D array with star's temperatures 
    c5dstar_opacities: 3D array with star's opacities
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

    # Declare np.array
    c5dstar_densities = np.zeros((nc5dedge,nc5dedge,nc5dedge))
    c5dstar_temperatures = np.zeros((nc5dedge,nc5dedge,nc5dedge))
    c5dstar_opacities = np.zeros((nc5dedge,nc5dedge,nc5dedge))

    # Extract densities - This can take time, some 2min per property
    for nx in range(nc5dedge):
        for ny in range(nc5dedge):
            for nz in range(nc5dedge):
                c5dstar_densities[nx,ny,nz] = c5ddata['Z'][0][0][34][nx][ny][nz]
                c5dstar_temperatures[nx,ny,nz] = c5ddata['EOS'][0][0][1][nx][ny][nz]
                c5dstar_opacities[nx,ny,nz] = c5ddata['OPA'][0][0][0][nx][ny][nz]
    
    return c5dstar_densities, c5dstar_temperatures, c5dstar_opacities


# Function for loading one dust specie from c5d-data
@cython.cfunc
@cython.locals(nspecies=cython.int)
def load_dustdensity(
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

# ==========================================================================
# Functions to create r3d-data from c5d-data
# Extract data on star and create r3d-files from it
def create_star(
        # All necessary inputs
        # paths!
        savpath:str='../co5bold_data/dst28gm06n052/st28gm06n052_186.sav',
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
        c5dstar_densities,c5dstar_temperatures,c5dstar_opacities = load_star_properties(savpath=savpath)

        # Start working :)
        print('Translating C5D data to R3D data')

        # Declare stuff for the loops
        r3d_densities = 0
        r3d_temperatures = 0
        r3d_opacities = 0
        progbar = 0

        # Adaptive range used for when cellsizes are too similar - fixes bug 
        # where I got a lot of zero-cells because my loop misses the c5dcells inside
        # the current r3dcell
        adaptive_range = c5dcellsize/r3dcellsizes.min() * 1.1

        # Open r3d data files
        with open('../dust_density_star.inp', 'w') as fdensity, open('../dust_temperature_star.dat', 'w') as ftemperature, open('../star_opacities.dat', 'w') as fopacity:

            # Write headers:
            # 1
            # nleafs
            # number dust species
            fdensity.write(f'1\n{int(nleafs)}\n1\n')
            ftemperature.write(f'1\n{int(nleafs)}\n1\n')
            fopacity.write('# List of c5d-opacities translated to r3d-spatial grid.\n# Use as input when separating one-specie-density_star-file into several species\n# and creating dust-star opacity files.\n')

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

                            # Sum all opacities
                            r3d_opacities += c5dstar_opacities[nnx,nny,nnz]

                # Check again if there are any c5dcells within r3dcell
                # If not, then your r3dgrid is probably larger than the c5dgrid
                # and then you can keep the density and temperature at zero
                if nchildcells > 0:
                    # Otherwise average the data of each r3dcell by number of c5dcells
                    r3d_densities /= nchildcells
                    r3d_temperatures /= nchildcells
                    r3d_opacities /= nchildcells

                # Then write data to r3d files
                fdensity.write(f'{r3d_densities}\n')
                ftemperature.write(f'{r3d_temperatures}\n')
                fopacity.write(f'{r3d_opacities}\n')

                # Reset data
                r3d_densities = 0
                r3d_temperatures = 0
                r3d_opacities = 0

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


# Takes the stellar data from c5d-data that has been translated to r3d but as one specie, and
# transforms it into separate species on a chosen number of bins of opacities as given by the
# c5d-data.
def create_staropacity(
        pathopacity:str='../star_opacities.dat',
        pathstardensity:str='../dust_density_star.inp',
        pathwavelength:str='../wavelength_micron.inp',
        pathtemperature:str='../dust_temperature.dat',
        nbins:int=5
    ):
    """
    nbins = number of species the star will consist of

    OUTPUT
    dust_density_star_opabins.inp
    dust_temperature_star_opabins.dat
    star_opacities_bins.dat
    dustopac_starbins.inp
    dustkappa_star{no+1}.inp
    """

    if nbins > 20:
        return 'WARNING too many files! Stopping'

    else:
        # Load denstiy_star.inp
        Ncells,Nspec,star_densities = a3d.load_dustdensity(path=pathstardensity,numb_specie=1)

        # Load temperature_star
        Ncells,temperatures = a3d.load_temperature(path=pathtemperature)

        # load opacity.dat
        opacity = []
        with open(pathopacity, 'r') as fopacity:
            for line in fopacity.readlines():
                if line[0] != '#':
                    opacity.append(float(line))

        # Load wavelengthgrid
        wavelengths,Nwave = a3d.load_wavelengthgrid(path=pathwavelength)

        # Bin the opacities
        print(f'Binning star-opacities to {nbins} species.')

        # Take all bins from histogram above the average bin-height
        # Only take one from those below (change later to some value dependant on the range?)
        # Adapt length (number bins) until the chosen number of bins-1 are achieved (done here with while)
        # The last bin is the average of the rest.
        valuestemp_length = 0
        bins = 0
        while valuestemp_length < nbins-1:
            bins += 1

            # Extract bins from histogram
            histvalues = np.histogram(opacity, range=(1,max(opacity)), bins=bins)

            # Take bins of those larger than average
            valuestemp = [
                histvalues[1][nn] for nn in np.where(histvalues[0] > np.mean(histvalues[0]))
            ]
            valuestemp_length = len(valuestemp[0])

        # Save the average value of each bin
        values = np.zeros(nbins)

        values[0] = np.array(opacity)[
            np.where(opacity < histvalues[1][1])[0]
        ].mean()

        for nn in range(1, nbins-1):
            values[nn] = np.array(opacity)[np.where(
                (opacity > histvalues[1][nn]) & (opacity < histvalues[1][nn+1])
            )[0]].mean()

        # For the last bin, take one average of all higher values
        values[-1] = np.array(opacity)[
            np.where(opacity > histvalues[1][nbins])[0]
        ].mean()

        # Bin the opacities to those that are nearest the bins
        opacitybins = np.zeros(len(opacity))
        for no,opac in enumerate(opacity):
            nn = (np.abs(values - opac)).argmin()
            opacitybins[no] = values[nn]

        print('Creating new density, temperature and binned opacity files.')
        # Create new density file with densities separated by species decided by the binning
        # of the opacities of the star
        new_densities = np.zeros(nbins*Ncells)
        new_temperatures = np.zeros(nbins*Ncells)
        for nn in range(Ncells):
            for no,opac in enumerate(values):
                if opacitybins[nn] == opac:
                    new_densities[nn + Ncells*no] = star_densities[nn]
                    new_temperatures[nn + Ncells*no] = temperatures[nn]

        # Print new star-density file
        with open('../dust_density_star_opabins.inp', 'w') as fdensity,\
             open('../dust_temperature_star_opabins.dat', 'w') as ftemperature,\
             open('../star_opacities_bins.dat', 'w') as fopacity:
            
            # Write headers
            fdensity.write(f'1\n{int(Ncells)}\n{int(nbins)}\n')
            ftemperature.write(f'1\n{int(Ncells)}\n{int(nbins)}\n')
            fopacity.write('# Binned opacities for the star(dust) distributions.\n# In same order as density file.\n# For reference and quality checks.\n')

            # Write densities and temperatures
            for nn,dens in enumerate(new_densities):
                fdensity.write(f'{dens}\n')
                ftemperature.write(f'{new_temperatures[nn]}\n')

                # And opacity bins, not necessary for simulations but useful for checking input data later
                if nn < Ncells:
                    fopacity.write(f'{opacitybins[nn]}\n')

        # Print new dustopac_starbins.inp file
        print('Writing opacity files for the binned star.')
        with open('../dustopac_starbins.inp', 'w') as fopac:

            # Print header
            fopac.write(f'2\n{nbins}\n-----------------------------\n')

            # Print star's opacity / species names
            for nn in range(nbins):
                fopac.write(f'1\n0\nstar{nn+1}\n-----------------------------\n')

        # Print opacity files, star-kappa-files
        # TODO perhaps add a factor to compensate since this kappa is just some constant number
        for no,opac in enumerate(values):
            with open(f'../dustkappa_star{no+1}.inp', 'w') as fopac:

                # Write header (1 for no scattering in this and number of wavelengths)
                fopac.write(f'1\n{Nwave}\n')

                # Write wavelength, abscoeff, scattercoeff
                for wave in wavelengths:
                    fopac.write(f'{wave}    {opac}    0.0\n')
        
        print('Done with: dust_density_star_opabins.inp, dust_temperature_star_opabins.inp, star_opacities_bins.dat, dustopac_starbins.inp, dustkappa_starX.inp')



# ====================================================================
# Funcs to load and create dusty envelope

# Extract and construct dust_density-files from C5D-data
# uses ['Z'][0][0][40][x][y][z] and ['Z'][0][0][43][x][y][z]
# or rather use
# ['Z'][0][0][40+3*nspecie][x][y][z]
# where nspecie is 0, 1, 2 ... up until (number of dust species)-1
# and the data are in number densities cm^-3 of monomers
# so the mass density is the number density * mass of dust forming molecule

# Also:
# str(teststar['Z'][0][0][40+3*nspecie+2])[4:-1]
# gives a string with the name of the specie!
# so this can also print the dustkappa-list-file for r3d!

def create_dustfiles(
        savpath:str='../co5bold_data/dst28gm06n052/st28gm06n052_186.sav',
        amrpath:str='../amr_grid.inp',
        gridpath:str='../grid_distances.csv',
        sizepath:str='../grid_cellsizes.csv',
        Nspecies:int=1,
        monomermasses:list=[2.3362e-22]
    ):
    """
    #graindensity:list=[2],
    #grainsizeum:list=[0.1],
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
                #grainsizecm = grainsizeum[nspecies]*1e-4

                print(f'Writing dust specie number {nspecies+1}:')
                print(f'    {speciesname}')
                print(f'    Monomer mass: {monomermasses} g')

                # Load c5d-dust densities and temperatures
                #c5ddensities, c5dtemperatures = load_dustdensity(savpath=savpath, nspecies=nspecies)

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
                    # Recalculate number density of monomers to mass density
                    # 4/3 pi agrain**3 * rhograin * numberdensity
                    # Units are in cgs.
                    #r3d_density *= \
                    #    4.1887902047863905 * grainsizecm**3 * graindensity[nspecies] / nchildcells
                    # TODO: add 
                    # Mg2SiO4: 2*24.305u + 28.085u + 4*15.999u = 140.69u = 2.3362e-22 gram
                    r3d_density *= monomermasses[nspecies] / nchildcells


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

                    if int(nr3d/nleafs*100) == 99 and progbar == 3:
                        print('Finished 99 per cent of the grid.')

    print('C5D Dust-data: done.\n')
