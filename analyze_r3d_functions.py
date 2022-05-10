# Various functions useful for analyzing in and out-data of RADMC3D
# ------------------------------------------------------------ #
# Useful packages
import os
import csv
from turtle import distance
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm, xlabel, xscale, yscale

# Useful numbers
c = 2.998e8 # speed of light in m/s
pc = 30.857e15 # 1 parsec in m
AUcm = 1.49598e13 # AU in cm

# ------------------------------------------------------------ #
# Functions that load various r3d input data

# TODO?
# Both these below need the number of species?

# Load grid properties
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

                # Number of refinements,
                # number of cells (nleafs),
                # number of branches
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


# Load R3d wavelength grid
def load_wavelengthgrid(path:str='../wavelength_micron.inp'):
    """
    Ancillary function. Loads and extracts wavelengths from wavelength grid.
    
    Input
    -----
    path: str, path to wavelength_micron.inp

    Output
    ------
    wavelengths: list of wavelengths in micron
    nwave: number of wavelength grid points    
    """

    wavelengths = []

    with open(path,'r') as f:

        for nn,line in enumerate(f.readlines()):
            
            # Extract number of wavelengths
            if nn == 0:
                nwave = int(line)
            
            # Extract wavelengths
            else:
                wavelengths.append(float(line))
            
    return wavelengths,nwave

def load_dustdensity(
        path:str='../dust_density.inp',
        numb_specie:int=1
    ):
    """
    Load and extracts densities of one dust specie of dust_density.inp

    INPUT
    -----
    path: string with path and filename of density file
    numb_specie: integer with number of the specie you want to load (default=1)
    
    OUTPUT
    ------
    Ncells: Number of cells in grid (nleafs)
    Nspec: Number of species in dust_density file
    dust_densities: np.array containing densities of all cells for specie numb_specie

    """

    # Read header
    with open(path,'r') as f:
        for nn,line in enumerate(f.readlines()):

            # Number of cells
            if nn == 1:
                Ncells = int(line)

            # Total number of species
            if nn == 2:
                Nspec = int(line)

    # Check that the chosen number of specie exists
    if numb_specie > Nspec:
        print('\nERROR\nThis dust specie doesnt exist.\n')

    else:
        # Reduce specie number by 1 (index starts at 0)
        numb_specie -= 1

        # Create density np.array
        dust_densities = np.zeros(Ncells)

        # Extract dust densities
        with open(path,'r') as f:
            for nn,line in enumerate(f.readlines()):
                if nn > 2+numb_specie*Ncells and nn <= 2+(numb_specie+1)*Ncells:
                    dust_densities[nn-3-Ncells*numb_specie] = float(line)

        return Ncells,Nspec,dust_densities

# ------------------------------------------------------------ #
# Load output data from R3D

# Load dust_temperature-file
def load_temperature(
        path:str='../dust_temperature.dat',
        numb_specie:int=1
    ):
    """
    Load and extracts temperatures in output file dust_temperature.dat

    INPUT
    -----
    path: string with path and filename of temperature file
    nspecie: number of specie you want to load
    
    OUTPUT
    ------
    Ncells: Number of cells in grid (nleafs)
    dust_temperatures: np.array containing temperatures in grid
    """

    # Read header
    with open(path,'r') as f:
        for nn,line in enumerate(f.readlines()):

            # Number of cells
            if nn == 1:
                Ncells = int(line)
            
            # Number of species
            if nn == 2:
                Nspecies = int(line)
    
    # Check that specie number is correct
    if numb_specie <= Nspecies and numb_specie > 0:
        # Create density np.array
        dust_temperatures = np.zeros(Ncells)

        # Extract dust densities (reduce specie number by 1, since python starts at index=0)
        ncounter = 0
        numb_specie -= 1

        with open(path,'r') as f:
            for nn,line in enumerate(f.readlines()):
                if (2 + Ncells*numb_specie) < nn < (2+ Ncells*(numb_specie+1)):
                    dust_temperatures[ncounter] = float(line)
                    ncounter += 1

        return Ncells,Nspecies,dust_temperatures
    
    else:
        # Otherwise, return error message
        print(f'ERROR: number of species in data is {Nspecies}, your specie-number was {numb_specie}')



# Load SED
def load_spectrum(
        path:str='../r3dsims/spectrum.out',
        distance:float=1
    ):
    """
    Loads and returns SED and wavelength.

    INPUT
    -----
    path: str with path and filename of spectrum.out file
    distance: float, distance to star in pc

    OUTPUT
    ------
    wavelengths: array with wavelength grid in microns
    spectrum: array with SED in Jy normalised to indicated distance
    """

    # Declare lists
    wavelengths = []
    spectrum = []

    # Load data
    with open(path, 'r') as f:

        # Data starts at row 4, so first data at nn == 3, where each row has two numbers
        for nn,line in enumerate(f.readlines()):
            if nn > 2:

                # Extract wavelengths (in um)
                wavelengths.append(float(line.split('   ')[0]))

                # Extract and recompute flux density to Jy (at 1pc for now)
                spectrum.append(float(line.split('   ')[1]) * 1e23 / distance**2)
    
    # Return data
    return wavelengths,spectrum


# ------------------------------------------------------------ #
# Plot various details of input and output data of R3D

def plot_onedensity_radius(
        density_path:str='../dust_density.inp',
        grid_path:str='../grid_distances.csv',
        amr_path:str='../amr_grid.inp',
        numb_specie:int=1
    ):
    """
    Plots one figure with radian density distribution of one dust species.

    INPUT
    density_path: path to density.inp-file
    grid_path: path to grid_distances.csv'
    amr_path: path to amr_grid.inp
    numb_specie: number of the dust specie you want to plot

    OUTPUT
    Shows figure
    """

    # Load griddistances
    griddistances = load_griddistances(
        gridpath=grid_path,
        amrpath=amr_path
    )/AUcm

    # load dust_density
    density = load_dustdensity(
        path=density_path,
        numb_specie=numb_specie
    )[2]

    # Load and plots r3d density data for ONE dust specie
    fig, ax = plt.figure(), plt.axes()
    ax.plot(
        griddistances[:,0],density,
        linestyle='',marker='.',markersize=1
    )
    ax.set(
        ylabel=r'Density (g\,cm$^{-3}$)',
        xlabel=r'Distance (AU)',
        title=f'Dust species {numb_specie}'
    )
    fig.show()


def plot_alldensity_radius(
        density_path:str='../dust_density.inp',
        grid_path:str='../grid_distances.csv',
        amr_path:str='../amr_grid.inp',
    ):
    """
    Plots one figure with radial density distribution of all dust species.

    INPUT
    density_path: path to density.inp-file
    grid_path: path to grid_distances.csv'
    amr_path: path to amr_grid.inp

    OUTPUT
    Shows figure
    """

    # Get info on number of dust species in 

    # Load griddistances
    griddistances = load_griddistances(
        gridpath=grid_path,
        amrpath=amr_path
    )/AUcm

    # load first dust_density
    Ncells,Nspec,density = load_dustdensity(
        path=density_path,
        numb_specie=1
    )
    densities = [density]

    # Load the rest of the dust densities and put them in a list of np.arrays
    if Nspec > 1:
        for nn in range(2,Nspec+1):
            densities.append(
                load_dustdensity(
                    path=density_path,
                    numb_specie=nn
                )[2]
            )
    
    # Control colours of each density distribution
    colour = cm.rainbow(np.linspace(0, 1, Nspec))

    # Set objects for plot with all in the same figure
    fig, ax = plt.figure(), plt.axes()
    
    for nn, c in enumerate(colour):
        ax.plot(
            griddistances[:,0],densities[nn],
            markeredgecolor=c,
            linestyle='',marker='.',markersize=1
        )
    ax.set(
        ylabel=r'Density (g cm$^{-3}$)',
        xlabel=r'Distance (AU)',
        title=f'Dust species 1 to {Nspec}'
    )
    fig.show()

    # Set objects for subplots (two columns, increasing number of rows)
    # ax[rows,columns]
    fig,ax = plt.subplots((-(-Nspec//2)),2)

    for nn, c in enumerate(colour):
        ax.ravel()[nn].plot(
            griddistances[:,0],densities[nn],
            markeredgecolor=c,
            linestyle='',marker='.',markersize=1
        )
        ax.ravel()[nn].set(
            ylabel=r'Density (g cm$^{-3}$)',
            xlabel=r'Distance (AU)',
            title=f'Dust specie {nn+1}'
        )
    fig.tight_layout()
    fig.show()

def plot_onetemperature_radius(
        temperature_path:str='../dust_temperature.dat',
        grid_path:str='../grid_distances.csv',
        amr_path:str='../amr_grid.inp',
        numb_specie:int=1
    ):
    """
    Plots one figure with radial temperature distribution of one specie.

    INPUT
    temperature_path: path to dust_temperature.dat
    grid_path: path to grid_distances.csv'
    amr_path: path to amr_grid.inp

    OUTPUT
    Shows figure
    """

    # Load griddistances
    griddistances = load_griddistances(
        gridpath=grid_path,
        amrpath=amr_path
    )/AUcm

    # load dust_temperature
    temperatures = load_temperature(
        path=temperature_path,
        numb_specie=numb_specie
    )[2]

    # Load and plots r3d density data for ONE dust specie
    fig, ax = plt.figure(), plt.axes()
    ax.plot(
        griddistances[:,0],temperatures,
        linestyle='',marker='.',markersize=1
    )
    ax.set(
        ylabel=r'Temperature (K)',
        xlabel=r'Distance (AU)',
        title='Grid cell temperatures'
    )
    fig.show()

# Plot all temperatures
def plot_alltemperature_radius(
        temperature_path:str='../dust_temperature.dat',
        grid_path:str='../grid_distances.csv',
        amr_path:str='../amr_grid.inp'
    ):
    """
    Plots figures with radial temperature distribution of all species within temperature file.

    INPUT
    temperature_path: path to dust_temperature.dat
    grid_path: path to grid_distances.csv'
    amr_path: path to amr_grid.inp

    OUTPUT
    Shows figure
    """

    # Load griddistances
    griddistances = load_griddistances(
        gridpath=grid_path,
        amrpath=amr_path
    )/AUcm

    # Load first dust_temperature
    Ncells,Nspec,temperature = load_temperature(
        path=temperature_path,
        numb_specie=1
    )
    temperatures = [temperature]

    # Load rest of temperatures
    if Nspec > 1:
        for numb_specie in range(2,Nspec+1):
            temperatures.append(
                load_temperature(
                    path=temperature_path,
                    numb_specie=numb_specie
                )[2]
            )
    
    # Control colours of each density distribution
    colour = cm.rainbow(np.linspace(0, 1, Nspec))

    # Set objects for plot with all in the same figure
    fig, ax = plt.figure(), plt.axes()
    
    for nn, c in enumerate(colour):
        ax.plot(
            griddistances[:,0],temperatures[nn],
            markeredgecolor=c,
            linestyle='',marker='.',markersize=1
        )
    ax.set(
        ylabel=r'Cell temperature (K)',
        xlabel=r'Distance (AU)',
        title=f'Dust species 1 to {Nspec}'
    )
    fig.show()

    # Set objects for subplots (two columns, increasing number of rows)
    # ax[rows,columns]
    fig,ax = plt.subplots((-(-Nspec//2)),2)

    for nn, c in enumerate(colour):
        ax.ravel()[nn].plot(
            griddistances[:,0],temperatures[nn],
            markeredgecolor=c,
            linestyle='',marker='.',markersize=1
        )
        ax.ravel()[nn].set(
            ylabel=r'Cell temperature (K)',
            xlabel=r'Distance (AU)',
            title=f'Dust specie {nn+1}'
        )
    fig.tight_layout()
    fig.show()

# Plot absorption, scattering, and angles of the various species

# load by name or load by specie-number? or both?

def plot_onekappa(
        specie_name:str='',
        specie_number:int=0,
        path:str='../'
    ):

    if specie_number != 0:
        # Extract name from dustopac-file
        with open(f'{path}dustopac.inp', 'r') as f:

            # line number 1+specie_number*4




        # specie_name = blabla
        print(f'Extracting species number {specie_number}: {specie_name}')



    elif len(specie_name) != 0:
        # Or just extract species by given name
        print(f'Extracting species {specie_name}')
    
    else:
        # Nothing is given, print error
        print('ERROR: no species number nor name is given.')
    
    with open(f'{path}dustkappa_{specie_name}.inp', 'r') as f:
        print('blabla')






    return 'hej'





# Plot SED
def plot_sed(
        path:str='../spectrum.out',
        distance:float=1
    ):

    wavelengths,spectrum = load_spectrum(path=path,distance=distance)
    
    # Load and plots r3d density data for ONE dust specie
    fig, ax = plt.figure(), plt.axes()
    ax.plot(
        wavelengths,spectrum,'b'
    )
    ax.set(
        ylabel=f'Flux density (Jy at {distance} pc)',
        xlabel=r'Wavelength ($\mu$m)',
        title='Output SED',
        xscale='log',yscale='log'
    )
    fig.show();

# Plot images
def plot_images(
        path:str='../',
        images:list=['image.out'],
        distance:float=1
    ):
    """
    Plots list of images

    INPUT
    path: path to folder containing images
    images: list of image-file names inside your folder
    distance: distance to source in pc (default 1 pc)

    OUTPUT
    Figure with 2 subplots, linear and logarithmic scales.
    Total flux of images in Jy at chosen distance
    """

    fluxtotal = []

    # Loop through all images
    for image in images:

        # Declare list
        image1d = []

        # Load images
        with open(path+image, 'r') as f:
            for nl,line in enumerate(f.readlines()):
                
                # row 1: pixels by pixels
                if nl == 1:
                    npixels = int(line.split()[0])
                
                # row 3: pixel size in cm, divide by AUcm for AU
                if nl == 3:
                    pixelsize_au = float(line.split()[0])/AUcm
                
                # row 4: wavelenght in um
                # TODO might be useful also :)
                if nl == 4:
                    wavelength = float(line)
                
                # row 6 onward: pixels
                if nl > 5:
                    # Some rows are empty (and they contain space and \n, so strip them lines)
                    if len(line.strip()) > 0:
                        image1d.append(float(line.strip()))

        # Extract some useful quantities
        # pixel size in mas
        pixelsize_mas = pixelsize_au / distance
        
        # Size of whole image in AU
        size_au = pixelsize_au * npixels
        axisplot  = [0.5*size_au,-0.5*size_au,-0.5*size_au,0.5*size_au]

        # Total flux density of the image
        fluxtotal.append(sum(image1d) * 1.e23 * 2.35044305391e-11 * pixelsize_mas**2)

        # Create 2D arrays
        image2d = np.zeros((npixels,npixels))
        image2dlog = np.zeros((npixels,npixels))
        nx,ny = 0,0

        for flux in image1d:
            image2d[nx,ny] = flux * 1.e23 * 2.35044305391e-11 * pixelsize_mas**2
            image2dlog[nx,ny] = np.log10(flux * 1.e23 * 2.35044305391e-11 * pixelsize_mas**2)

            # Move nx and ny
            nx = nx + 1
            if nx == npixels:
                nx = 0
                ny = ny + 1

        # Plot each image (one in linear and one log)
        fig, ax = plt.subplots(1,2, dpi=150, num=image)
        
        im0 = ax[0].imshow(
            image2d, origin='lower', extent=axisplot
        )
        ax[0].set(
            title=f"{image.replace('image_', '').replace('.out', '')} (Lin)", 
            xlabel='Offset (AU)',
            ylabel='Offset (AU)'
        )
                
        ax[1].imshow(
            image2dlog, origin='lower', extent=axisplot
        )
        ax[1].set(
            title=f"{image.replace('image_', '').replace('.out', '')} (Log)", 
            xlabel='Offset (AU)',
        )

        cb0 = plt.colorbar(im0, orientation = 'vertical',shrink=0.6,pad=0.15)
        cb0.set_label(label = rf'Flux at {distance} pc (Jy/asec$^2$)',fontsize= 10)

        fig.tight_layout()
        fig.show()

    return fluxtotal


# ------------------------------------------------------------ #
# Compute different quantities

def compute_luminosity(path:str='../r3dsims/spectrum.out',distance:float=1):
    """
    Insert a spectrum.out from r3d and get the bolometric luminosity in Watt
    INPUT
    path = path to spectrum-file including file name
    OUTPUT
    luminosity in Watt
    """

    # Load spectrum
    wavelengths,spectrum = load_spectrum(path,1)
    nwave = len(wavelengths)

    # Integrate the SED (using trapezoidal method, and change units to SI units)
    sedintegral = 0
    for nn in range(nwave-1):
        sedintegral += 0.5*(spectrum[nn] + spectrum[nn+1])*1e-26 * (c/wavelengths[nn] - c/wavelengths[nn+1])*1e6

    # Compute bolometric luminosity(?)
    luminosity = 4.*np.pi*(distance*pc)**2. * sedintegral

    return luminosity
