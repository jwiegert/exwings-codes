# functions and tools for loading and analyzing output
# from co5bold.

# Import various libraries
import cython
import numpy as np
from scipy.io.idl import readsav

# My own libraries
import create_r3d_functions as c3d

# Define useful numbers
AUcm = cython.declare(cython.float ,1.49598e13) # cm

# ============================================================

# TODO
# dela upp i flera funktioner?
# en som bara returnerar en cell åt gången?

# Separate so that grid and densities are loaded separately
# Also, one number only so that I can loop outside the function

# this decorator when declaring what's included in the functions
# @cython.locals(a=cython.int)




# Load co5bold grid
@cython.cfunc
def load_co5boldgrid(
        savpath:str='../co5bold_data/dst28gm06n056/st28gm06n056_140.sav'
    ):
    
    # Load sav-file
    c5ddata = readsav(savpath)

    # Extract data - TODO MAY HAVE TO CHOSE THIS MANUALLY DEPENDING ON FILE
    c5ddata = c5ddata['ful']

    # Get number of gridcells from co5bold data
    nc5dedge = np.size(c5ddata['Z'][0][0][16])
    nc5d = nc5dedge**3

    # Declare array
    c5dgrid = np.zeros((nc5d,4))

    # Declare coordinates
    nx,ny,nz = 0,0,0

    for nn in range(nc5d):

        # Extract various data (yes I may not need these variables now)
        c5dx = c5ddata['Z'][0][0][16][0][0][nx]
        c5dy = c5ddata['Z'][0][0][19][0][ny][0]
        c5dz = c5ddata['Z'][0][0][22][nz][0][0]

        # Save radial distances in cm
        c5dgrid[nn,0] = np.sqrt(
            c5dx**2 + c5dy**2 + c5dz**2
        )

        # Save x,y,z distances in cm
        c5dgrid[nn,1] = c5dx
        c5dgrid[nn,2] = c5dy
        c5dgrid[nn,3] = c5dz
 
        # Move coordinates
        nx,ny,nz = c3d.movecoordinates(nc5dedge,nx,ny,nz)

    return c5dgrid



# Compare grids and return corresponding cell identity
@cython.cfunc
@cython.locals(
    c5dx = cython.double,
    c5dy = cython.double,
    c5dz = cython.double
)
def compare_grids(
    c5dx,c5dy,c5dz,
    r3dgrid
    ):

    nr = cython.declare(cython.int, 0)

    # Find c5d cells are close to r3d cells
    nr = np.argmin(
        (r3dgrid[:,1] - c5dx)**2 + \
        (r3dgrid[:,2] - c5dy)**2 + \
        (r3dgrid[:,3] - c5dz)**2
    )

    return nr



    
# Extract densities 
@cython.cfunc
def load_co5bolddensity2(
        c5dgrid,c5ddensities,
        r3dgrid
    ):

    # Declare arrays and settings on variables
    nc = cython.declare(cython.int)
    nr = cython.declare(cython.int)
    nr3d = cython.declare(cython.int, np.size(r3dgrid[:,0]))
    nc5d = cython.declare(cython.int, np.size(c5dgrid[:,0]))

    c5dx = cython.declare(cython.double)
    c5dy = cython.declare(cython.double)
    c5dz = cython.declare(cython.double)

    comparisonarray = np.zeros(nr3d)
    densityarray = np.zeros(nr3d)


    # Loop and extract densities
    for nc in range(nc5d):

        # Extract c5d coordinates
        c5dx = c5dgrid[nc,1]
        c5dy = c5dgrid[nc,2]
        c5dz = c5dgrid[nc,3]


        # Find c5d cells are close to r3d cells
        nr = np.argmin(
            (r3dgrid[:,1] - c5dx)**2 + \
            (r3dgrid[:,2] - c5dy)**2 + \
            (r3dgrid[:,3] - c5dz)**2
        )

        # Save number of c5dcells per r3dcell
        comparisonarray[nr] += 1

        # Save densities
        densityarray[nr] += c5ddensities[nc]

    # Average the densities by number of c5dcells per r3dcells
    densityarray = densityarray / comparisonarray

    return densityarray




"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.io.idl import readsav

a=readsav('st29gm04n002_seq.sav')
seq = a['seq']

# print(a.keys())  # get restored variables

# loaded sav file is as a dict file, with each restored variable as a numpy.recarray file
# (i.e. an ndarray that allows field access using attributes) so get the tag names
# (tag_names in IDL) via the dtype.names
tag_names = a['seq'].dtype.names
# print(tag_names)

# then can access tag data via
# B = a['seq']['tag_name'][0]

# misc
# can always check type(XX) and if array, check dims via .shape

# e.g. plot the radial velocity
B = a['seq']['v3'][0]
B = (1e-5 * seq.rhov3 / seq.rho)[0]
xx = a['seq']['time'][0] / (365*24*3600)
xx = xx - xx[0]
yy = seq.xc3_radial[0] / 6.98000e+10
hh = B

plt.pcolormesh(xx, yy, hh.transpose(), cmap=plt.cm.get_cmap('bwr'), shading='auto')
cbar=plt.colorbar()
cbar.set_label('Radial velocity ' + r'$[kms^{-1}]$', rotation=270, labelpad=25)
plt.xlabel(r'Time $[yr]$')
plt.ylabel(r'$R \: \: [R_{\odot}]$')
plt.tick_params(direction='in', top ='on', right = 'on')
plt.show()
"""


# TEMPORARY MAY REMOVE
# Load co5bold grid, returns radial, x, y, z distances, densities
@cython.cfunc
def load_co5boldgrid2(
        savpath:str='../co5bold_data/dst28gm06n056/st28gm06n056_140.sav'
    ):
    
    # Load sav-file
    teststar = readsav(savpath)

    # Extract data - TODO MAY HAVE TO CHOSE THIS MANUALLY DEPENDING ON FILE
    teststar = teststar['ful']

    # Get number of gridcells from co5bold data
    nc5dedge = np.size(teststar['Z'][0][0][16])
    nc5d = cython.declare(cython.int, nc5dedge**3)

    # Declare arrays   
    c5ddensities = np.zeros(nc5d)
    c5dgrid = np.zeros((nc5d,4))

    # Declare coordinates
    nx = cython.declare(cython.int, 0)
    ny = cython.declare(cython.int, 0)
    nz = cython.declare(cython.int, 0)

    c5dx = cython.declare(cython.double)
    c5dy = cython.declare(cython.double)
    c5dz = cython.declare(cython.double)

    for nn in range(nc5d):

        # Save densities
        c5ddensities[nn] = teststar['Z'][0][0][34][nx][ny][nz]

        c5dx = teststar['Z'][0][0][16][0][0][nx]
        c5dy = teststar['Z'][0][0][19][0][ny][0]
        c5dz = teststar['Z'][0][0][22][nz][0][0]

        # Save radial distances in cm
        c5dgrid[nn,0] = np.sqrt(
            c5dx**2 + c5dy**2 + c5dz**2
        )

        # Save x,y,z distances in cm
        c5dgrid[nn,1] = c5dx
        c5dgrid[nn,2] = c5dy
        c5dgrid[nn,3] = c5dz
 
        # Move coordinates
        nx,ny,nz = c3d.movecoordinates(nc5dedge,nx,ny,nz)

    return c5dgrid,c5ddensities