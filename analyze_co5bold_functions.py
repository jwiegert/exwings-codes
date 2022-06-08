# functions and tools for loading and analyzing output
# from co5bold.

# Import various libraries
import cython
import numpy as np
from scipy.io.idl import readsav
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm, xlabel, xscale, yscale

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

# Plot c5d opacities (in r3d-grid)
def plot_opakapparadius(
        path:str='../'
    ):
    """
    Plots star_opacity.dat
    I.e. the kappa per r3dgridcell as given by c5d.
    Input: 
    path to folder of your r3d-model where amr_grid.inp, grid_distcances.csv, 
    and star_opacities.dat are included.
    """

    # Automatically add / to end of path if it's missing
    if path[-1] != '/':
        path += '/'

    # load opacity.dat
    opacity = []
    with open(path+'star_opacities.dat', 'r') as fopacity:
        for line in fopacity.readlines():
            if line[0] != '#':
                opacity.append(float(line))

    # load grid
    griddistances = a3d.load_griddistances(
        amrpath=path+'amr_grid.inp',
        gridpath=path+'../grid_distances.csv'
    )
    radiusau = griddistances[:,0]/AUcm

    # Load and plots r3d density data for ONE dust specie
    fig, ax = plt.subplots(2,1)

    
    ax[0].plot(
        radiusau,opacity,
        linestyle='',marker='.',markersize=1
    )
    ax[0].set(
        ylabel=r'$\kappa_{\rm abs}$ (cm$^2$/g)',
        xlabel=r'Distance (AU)',
        title='Grid cell OPA'
    )

    ax[1].plot(
        radiusau,opacity,
        linestyle='',marker='.',markersize=1
    )
    ax[1].set(
        ylabel=r'$\kappa_{\rm abs}$ (cm$^2$/g)',
        xlabel=r'Distance (AU)',
        yscale='log'
    )

    fig.show()

# Plot c5d-OPA vs r3d-temperature
def plot_opakappatemperature(
        path:str='../'
    ):

    # Automatically add / to end of path if it's missing
    if path[-1] != '/':
        path += '/'

    # Load opacity.dat
    opacity = []
    with open(path+'star_opacities.dat', 'r') as fopacity:
        for line in fopacity.readlines():
            if line[0] != '#':
                opacity.append(float(line))
    opacity = np.array(opacity)


    # Load all opacities as given for R3D

    # First extract specie names from dustopac-file
    counter = 5
    specie_names = []
    with open(f'{path}dustopac.inp', 'r') as f:
        for nn,line in enumerate(f.readlines()):
            
            # Specie-names are on line numbers 1+specie_number*4
            if nn == counter:
                specie_names.append(line.strip())
                counter +=4

    # Extract all opacities
    kappas = []

    for specie_name in specie_names:
        specie_name,wavelengths,kappadata = a3d.load_onekappa(
            specie_name=specie_name,
            path=path
        )

        # Save kappa from each specie
        kappas.append(np.mean(np.array(kappadata[0])))


    # Load dust_temperature.dat
    temperature_path = path+'dust_temperature.dat'
    # Load first dust_temperature
    Ncells,Nspec,temperature = a3d.load_temperature(
        path=temperature_path,
        numb_specie=1
    )
    temperatures = [temperature]

    # Load rest of temperatures
    if Nspec > 1:
        for numb_specie in range(2,Nspec+1):
            temperatures.append(
                a3d.load_temperature(
                    path=temperature_path,
                    numb_specie=numb_specie
                )[2]
            )
    
    # Control colours of each density distribution
    colour = cm.rainbow(np.linspace(0, 1, Nspec))

    # Set objects for plot with all in the same figure
    fig, ax = plt.figure(), plt.axes()
    
    for nn, c in enumerate(colour):

        temperature = np.array(temperatures[nn])

        ax.plot(
            temperature[np.where(temperature > 0)[0]],
            opacity[np.where(temperature > 0)[0]],
            markeredgecolor=c,linestyle='',marker='.',markersize=1
        )

        # Plot kappa of each bin
        ax.plot(
            [
                np.min(temperature[np.where(temperature > 0)[0]]),
                np.max(temperature[np.where(temperature > 0)[0]])
            ],
            [kappas[nn],kappas[nn]],
            'k-'
        )

    ax.set(
        xlabel=r'Cell temperature (K)',
        ylabel=r'$\kappa_{\rm abs}$ (cm$^2$/g)',
        yscale='log',
        title=f'Dust species 1 to {Nspec}'
    )
    fig.tight_layout()
    fig.show()

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
    R3D density file: dust_density_onestar.inp
    R3D temperature file: dust_temperature_onestar.dat
    Useful file with list of extracted opacities 'kappa': star_opacities.dat
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
        with open('../dust_density_onestar.inp', 'w') as fdensity, open('../dust_temperature_onestar.dat', 'w') as ftemperature, open('../star_opacities.dat', 'w') as fopacity:

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

    print('C5D Dust-star:\n    dust_density_onestar.inp\n    dust_temperature_onestar.dat\n    star_opacities.dat\nDONE\n')


# Takes the stellar data from c5d-data that has been translated to r3d but as one specie, and
# transforms it into separate species on a chosen number of bins of opacities as given by the
# c5d-data.
def create_staropacity(
        pathopacity:str='../star_opacities.dat',
        pathstardensity:str='../dust_density_onestarstar.inp',
        pathwavelength:str='../wavelength_micron.inp',
        pathtemperature:str='../dust_temperature_onestar.dat',
        temperaturebins:list=[],
        N_opabins:int=2
    ):
    """
    INPUT
    pathopacity: path to star_opacities.dat',
    pathstardensity: path to dust_density_onestarstar.inp',
    pathwavelength: path to wavelength_micron.inp',
    pathtemperature: path to dust_temperature_onestar.dat',
    temperaturebins: list of temperature range bins (it adds lowest and higehst temperature automatically)
    N_opabins: number of kappa-bins spread logarithimically between min and max

    OUTPUT
    dust_density_starbins.inp
    dust_temperature_starbins.dat
    star_opacities_bins.dat
    dustopac_starbins.inp
    dustkappa_star{no+1}.inp
    """

    Tbins = len(temperaturebins)

    if Tbins > 20:
        return 'WARNING too many files! Stopping'
    
    elif Tbins == 0:
        print('Only one bin, will only create r3d-opacity files, use *_onestar.inp files for density etc')
        
        # load opacity.dat
        opacity = []
        with open(pathopacity, 'r') as fopacity:
            for line in fopacity.readlines():
                if line[0] != '#':
                    opacity.append(float(line))

        # Load wavelengthgrid
        wavelengths,Nwave = a3d.load_wavelengthgrid(path=pathwavelength)

        # Take a SINGLE average of ALL opacity in c5d opacity-data
        opacityvalue = np.mean(np.array(opacity))

        # Save as constant as function of wavelength/frequency

        # Print new dustopac_starbins.inp file
        print('Writing opacity files for the star.')
        with open('../dustopac_starbins.inp', 'w') as fopac:

            # Print header
            fopac.write(f'2\n1\n-----------------------------\n')

            # Print star's opacity / species names
            fopac.write(f'1\n0\nstar\n-----------------------------\n')

        # Print opacity file, star-kappa-file
        with open(f'../dustkappa_star.inp', 'w') as fopac:

            # Write header (1 for no scattering in this and number of wavelengths)
            fopac.write(f'1\n{Nwave}\n')

            # Write wavelength, abscoeff, scattercoeff
            for wave in wavelengths:
                fopac.write(f'{wave}    {opacityvalue}    0.0\n')

        print('C5D create star opacities:\n    dustopac_starbins.inp\n    dustkappa_star.inp\nDONE\n')

    else:
        # Load denstiy_star.inp
        Ncells,Nspec,star_densities = a3d.load_dustdensity(path=pathstardensity,numb_specie=1)

        # Load temperature_star
        Ncells,Nspec,temperatures = a3d.load_temperature(path=pathtemperature)

        # load opacity.dat
        opacity = []
        with open(pathopacity, 'r') as fopacity:
            for line in fopacity.readlines():
                if line[0] != '#':
                    opacity.append(float(line))
        opacity = np.array(opacity)

        # Load wavelengthgrid
        wavelengths,Nwave = a3d.load_wavelengthgrid(path=pathwavelength)

        # Add lowest temperature to temperature list, change to array and update Tbins
        # The lowest temperature is not 0, this if statements makes sure it's the smallest
        # real temperature that is saved.
        if np.min(temperatures) == 0:
            temperaturebins.append(np.sort(temperatures)[1])
        else: 
            temperaturebins.append(np.min(temperatures))

        temperature_bins = np.array(temperaturebins)
        Tbins = np.size(temperature_bins)

        print('Temperature bin-ranges are (K):')
        for temp in temperature_bins:
            print(f'        {temp}')

        # Define opacity bins, no more than one per order of magnitude
        # Minimum opacity bin is the minimum value of opacity
        opamin = np.log10(min(opacity))

        # Maximum opacity is defined from above the highest tmeperature limit
        opamax = np.log10(np.mean(
            opacity[ np.where(temperatures > temperature_bins[0])[0]]
        ))
        
        # Array with OPA-bin-limits (numb of limits are one more)
        if N_opabins < 1:
            N_opabins = 1
        Obins = N_opabins + 1
        opabins = np.logspace(opamax, opamin, Obins)

        print(f'Opacity bin-ranges are (cm2/g):')
        for opa in opabins:
            print(f'        {opa}')

        # Bin the opacities by temperature ranges (the stellar surface is never above eg 5000k)
        # and by opacity range

        # Figure out number of bins (at least one bin always)
        nbin = 0
        nbins = 1
        listbins = np.zeros(Tbins*Obins)
        listbins[0] = 1

        for nT in range(1,Tbins):
            for nO in range(1,Obins):
                
                nbin += 1
                sumopacity = 0

                for nn in range(Ncells):
                    if \
                        temperature_bins[nT-1] > temperatures[nn] > temperature_bins[nT] and \
                        opabins[nO-1] > opacity[nn] > opabins[nO]:

                        sumopacity += opacity[nn]
                
                if sumopacity != 0:
                    nbins += 1

                    # Listbins lists the corresponding bin-number for each combination of nT and nO
                    # 0 where there is none
                    listbins[nbin] = nbins
                   

        # Start binning
        print(f'Binning to {nbins} species.')

        new_temperatures = np.zeros(nbins*Ncells)
        new_densities = np.zeros(nbins*Ncells)
        new_opacities = np.zeros(nbins*Ncells)

        for nn in range(Ncells):

            # First bin is inside the star, above highest T or OPA
            if temperatures[nn] > temperature_bins[0] or opacity[nn] > opabins[0]:
                new_temperatures[nn] = temperatures[nn]
                new_densities[nn] = star_densities[nn]
                new_opacities[nn] = opacity[nn]

            bincounter = 0

            # Remaining bins are outside and within these ranges
            for nT in range(1,Tbins):
                for nO in range(1,Obins):

                    bincounter += 1

                    if temperature_bins[nT-1] > temperatures[nn] > temperature_bins[nT] and \
                        opabins[nO-1] > opacity[nn] > opabins[nO]:

                        nbin = int(listbins[bincounter] - 1)

                        new_temperatures[nn + Ncells*nbin] = temperatures[nn]
                        new_densities[nn + Ncells*nbin] = star_densities[nn]
                        new_opacities[nn + Ncells*nbin] = opacity[nn]


        # Create new opacities from these limits
        opacityvalues = np.zeros(nbins)

        for nbin in range(1,nbins+1):

            # If there are no numbers then just let the opacity be 0
            if np.where(new_opacities[Ncells*(nbin-1):Ncells*nbin] > 0)[0].size > 0:

                opacityvalues[nbin-1] = np.mean(
                    new_opacities[
                        Ncells*(nbin-1) + np.where(new_opacities[Ncells*(nbin-1):Ncells*nbin] > 0)[0]
                    ]
                )

            #print(f'{nbin}: {Ncells*(nbin-1)}:{Ncells*nbin} OPAC: {opacityvalues[nbin-1]}')
            
        # Print new star-density, temperature and opacitybin-files
        print('Writing new radmc3d-files')

        with open('../dust_density_starbins.inp', 'w') as fdensity,\
             open('../dust_temperature_starbins.dat', 'w') as ftemperature,\
             open('../star_opacities_bins.dat', 'w') as fopacity:
            
            # Write headers
            fdensity.write(f'1\n{int(Ncells)}\n{int(nbins)}\n')
            ftemperature.write(f'1\n{int(Ncells)}\n{int(nbins)}\n')
            fopacity.write('# Binned opacities for the star(dust) distributions.\n# In same order as density file.\n# For reference and quality checks.\n')

            # Write densities and temperatures
            for nn,dens in enumerate(new_densities):
                fdensity.write(f'{dens}\n')
                ftemperature.write(f'{new_temperatures[nn]}\n')
                fopacity.write(f'{new_opacities[nn]}\n')

        # Print new dustopac_starbins.inp file
        print('Writing opacity files for the binned star.')
        with open('../dustopac_starbins.inp', 'w') as fopac:

            # Print header
            fopac.write(f'2\n{nbins}\n-----------------------------\n')

            # Print star's opacity / species names
            for nn in range(nbins):
                fopac.write(f'1\n0\nstar{nn+1:02d}\n-----------------------------\n')

        # Print opacity files, star-kappa-files
        for nbin,opac in enumerate(opacityvalues):
            with open(f'../dustkappa_star{nbin+1:02d}.inp', 'w') as fopac:

                # Write header (1 for no scattering in this and number of wavelengths)
                fopac.write(f'1\n{Nwave}\n')

                # Write wavelength, abscoeff, scattercoeff
                for wave in wavelengths:
                    fopac.write(f'{wave}    {opac}    0.0\n')
        
        # TODO Write info-file with data on how the binning is done! stellarbinning_info.txt
        print('C5D create star opacities:\n    dust_density_starbins.inp\n    dust_temperature_starbins.inp\n    star_opacities_bins.dat\n    dustopac_starbins.inp\n    dustkappa_starN.inp\nDONE\n')


# Takes the stellar density data and opacity data from c5d
# creates a density_star-file where the "density" depends on each cell's mass density
# and the opacity-kappa of each cell as given by c5d.
# Creates ONE kappa_star that is only = 1 over all wavelengths.
def create_staropadensity(
        pathopacity:str='../star_opacities.dat',
        pathstardensity:str='../dust_density_onestarstar.inp',
        pathtemperature:str='../dust_temperature_onestar.dat',
        pathwavelength:str='../wavelength_micron.inp',
        Teff:float=2800,
        kramer_exponent:float=3.5
    ):
    """
    INPUT
    pathopacity: path to star_opacities.dat',
    pathstardensity: path to dust_density_onestarstar.inp',
    pathwavelength: path to wavelength_micron.inp',

    OUTPUT
    dust_density_opastar.inp
    dustopac_star.inp
    dustkappa_opastar.inp
    """

    print('Loading density, opacity, wavelengths')

    # Load star densities (in r3d-grid)
    Ncells,Nspec,star_densities = a3d.load_dustdensity(path=pathstardensity,numb_specie=1)

    # load star opacities (in r3d-grid)
    opacity = []
    with open(pathopacity, 'r') as fopacity:
        for line in fopacity.readlines():
            if line[0] != '#':
                opacity.append(float(line))
    opacity = np.array(opacity)

    # Load temperature_star
    Ncells,Nspec,temperatures = a3d.load_temperature(path=pathtemperature)

    # Load wavelengthgrid
    wavelengths,Nwave = a3d.load_wavelengthgrid(path=pathwavelength)
    wavelengths = np.array(wavelengths)

    # Create gas opacity
    #kappa = corrfact * np.median(wavelengths) / wavelengths


    print('Change density to kappa * density')

    # Assume
    # r3dopacity * r3ddensity = c5dopacity * c5ddensity
    #
    # Opacity file is set to 1 for all wavelengths
    # r3dopacity = 1
    # This means that
    # r3ddensity = c5dopacity * c5ddensity
    # for each cell
    #
    # Thus we obtain a simple grey body model for the opacity for all cells without having
    # to add hundreds and hundreds of opacity files.

    opacity_densities = np.zeros(Ncells)
    mintemperature = np.min(temperatures)
    meandensity = np.mean(star_densities)

    for nn in range(Ncells):
        #opacity_densities[nn] = opacity[nn] * star_densities[nn]

        if temperatures[nn] == 0:
            opacity_densities[nn] = opacity[nn] * star_densities[nn] * mintemperature**(-kramer_exponent)*Teff**(kramer_exponent)
        else:
            opacity_densities[nn] = opacity[nn] * star_densities[nn] * temperatures[nn]**(-kramer_exponent)*Teff**(kramer_exponent)

    # Print new star-opacity-density file
    print('Writing new radmc3d-files')

    with open('../dust_density_opastar.inp', 'w') as fdensity:
        
        # Write headers
        fdensity.write(f'1\n{int(Ncells)}\n1\n')

        # Write densities and temperatures
        for nn,dens in enumerate(opacity_densities):
            fdensity.write(f'{dens}\n')

    # Print new dustopac_starbins.inp file
    print('Writing opacity files for the binned star.')
    with open('../dustopac_star.inp', 'w') as fopac:

        # Print header
        fopac.write(f'2\n1\n-----------------------------\n')
        # Print specie name
        fopac.write(f'1\n0\nopastar\n-----------------------------\n')

    with open(f'../dustkappa_opastar.inp', 'w') as fopac:

        # Write header (1 for no scattering in this and number of wavelengths)
        fopac.write(f'1\n{Nwave}\n')

        # Write wavelength, abscoeff, scattercoeff
        for nn,wave in enumerate(wavelengths):
            fopac.write(f'{wave}    1.0    0.0\n')
    
    print('C5D create star opacities densities:\n    dust_density_opastar.inp\n    dustopac_star.inp\n    dustkappa_opastar.inp\nDONE\n')
















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
