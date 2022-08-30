# For usage with Bash-scripting
# creates grid based on bash inputs
#
# > python3 scriptpy_creategrid.py $MODELNAME $PHASE1 ... $PHASE{N}
#
# And from data written in here, change codes manually if required!
#
# ------------------------------------------------------------------------
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
for phase in phases:
    os.system(f'mkdir {path}/{phase}')


# Extract star's radius for each phase
Mstar,Rstar,Lstar = a5d.load_star_information(
    savpath = f'../co5bold_data/d{modelname}/{modelname}_{phases[0]}.sav',
    printoutput = 'y'
)

# Extract minimu cellsize
c5dgrid,cellcourners,cellsize = a5d.load_grid_properties(
    savpath = f'../co5bold_data/d{modelname}/{modelname}_{phases[0]}.sav'
)

# Create spatial grid
#
# Check numbers from headers of c5d-data-sets
#
# use a5d.load_grid_properties() for this.
# Gives minimum cellsize, cellcourner-coords, c5dgrid-xyz-distances to cells
#
# In this example:
#
# Smallest c5d cells are 2*3.65/317 AU = 0.02302839116719243 AU
# Radius of star: 1.651AU (355 Rsun)
# Settings for the grid
r3dedge = 2.1 * Rstar # Size of whole grid
#basecell = 1.001 * cellsize * 2**4 # Size of base cells (based on smallest cells)
basecell = 2* 1.001 * cellsize # Size of base cells (based on smallest cells)
nxyz = r3dedge/basecell

# Refinements based on stellar radius
refinementlist_au = [
    3.0*Rstar/AUcm,
    2.5*Rstar/AUcm,
    2.0*Rstar/AUcm,
    1.5*Rstar/AUcm
]

# Inner refinements up to 0.9 Rstar
inrefine_au = 0.9*Rstar / (4*refinementlist_au[-1]*AUcm)

c3d.create_grid(
    gridedge=r3dedge/AUcm, 
    nxyz=nxyz, 
    refinementlist=refinementlist_au, 
    inrefine=inrefine_au, 
    savegrid=True
)


# Move grid files to correct folder
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
    wavelengthend = 1000.0,
    nwave = 100,
    logscale = 'y'
)

# Move wavelengthgrid to data-folder
os.system(
    f'mv ../wavelength_micron.inp {path}'
)

for phase in phases:
    # And copy wavelength-grid to each phase-folder
    os.system(
        f'cp {path}wavelength_micron.inp {path}{phase}/'
    )

