# For usage with Bash-scripting
# creates grid based on bash inputs
#
# > python3 scriptpy_creategrid.py $MODELNAME $PHASE1 ... $PHASE{N}
#
# And from data written in here, change codes manually if required!
#
# ------------------------------------------------------------------------
# Empty line in Bash to make some space
print('')
#
# Include inputs from bash

import sys
modelname = sys.argv[1]
phases = [phase for nn,phase in enumerate(sys.argv) if nn > 1]

# ------------------------------------------------------------------------
# Pythoncodes below

import create_r3d_functions as c3d
import analyze_co5bold_functions as a5d
import os

path = f'../r3dresults/{modelname}/'
AUcm = 1.49598e13 # cm

# Create folders
os.system(f'mkdir {path}')
for phase in phases:
    os.system(f'mkdir {path}{phase}')

# Extract star's radius for first phase (all of same model should have same settings
Mstar,Rstar,Lstar,Tstar = a5d.load_star_information(
    savpath = f'../../exwings_archivedata/co5bold_data/d{modelname}/{modelname}_{phases[0]}.sav',
    printoutput = 'y'
)
# Extract minimum cellsize for first phase (grid is the same for whole model)
c5dgrid,cellcourners,cellsize = a5d.load_grid_properties(
    savpath = f'../../exwings_archivedata/co5bold_data/d{modelname}/{modelname}_{phases[0]}.sav'
)


# Create spatial grid
#    Check numbers from headers of c5d-data-sets
#
#    use a5d.load_grid_properties() for this.
#    Gives minimum cellsize, cellcourner-coords, c5dgrid-xyz-distances to cells
#
# In this example:
#    Smallest c5d cells are 2*3.65/317 AU = 0.02302839116719243 AU
#    Radius of star: 1.651AU (355 Rsun)
#
#    light c5d-data: cube side ~  7AU ~  4Rstar
#    heavy c5d-data: cube side ~ 30AU ~ 18Rstar > max radius 9Rstar
#
# Settings for the grid

# Size of whole grid in AU
r3dedge = 2 * cellcourners.max() / AUcm 

# Size of base cells as based on smallest cells
basecell = 1.001 * cellsize * 2**4 / AUcm

# Number of base cells along one side of the cube
nxyz = r3dedge/basecell

# Radial distances to refinements based on stellar radius
refinementlist_au = [
    8*Rstar/AUcm,
    6*Rstar/AUcm,
    4*Rstar/AUcm,
    2*Rstar/AUcm
]

# Inner refinements up to 0.9 Rstar (divided equally radially by number of refinements)
inrefine_au = 0.6*Rstar / AUcm

# Create spatial grid
c3d.create_grid(
    gridedge=r3dedge, 
    nxyz=nxyz, 
    refinementlist=refinementlist_au, 
    inrefine=inrefine_au, 
    savegrid=True
)

# Move spatial grid files to correct folder
os.system(
    f'mv ../amr_grid.inp {path}'
)
os.system(
    f'mv ../grid_info.txt {path}'
)
os.system(
    f'mv ../grid_distances.csv {path}'
)
os.system(
    f'mv ../grid_cellsizes.csv {path}'
)


# Create wavelength grid
wavelengths = c3d.create_wavelength(
    wavelengthstart = 0.1,
    wavelengthend = 500.0,
    nwave = 1000,
    logscale = 'y'
)

# Move wavelengthgrid to data-folder
os.system(
    f'mv ../wavelength_micron.inp {path}'
)

# NOTE
# Resulting files are, for each ../r3dresults/[modelname]/
#
#   amr_grid.inp
#   grid_info.txt
#   grid_distances.csv
#   grid_cellsizes.csv
#   wavelength_micron.inp

