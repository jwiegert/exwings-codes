# For usage with Bash-scripting
# extracts dust and creates c3d-dust-files based on bash inputs
#
# > python3 scriptpy_createdust.py $MODELNAME $PHASE1 ... $PHASE{N}
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
phase = sys.argv[2]

# ------------------------------------------------------------------------
# Pythoncodes below

import analyze_co5bold_functions as a5d
import create_r3d_functions as c3d
import os


# First create dust files
a5d.create_dustfiles(
    savpath=f'../co5bold_data/d{modelname}/{modelname}_{phase}.sav',
    amrpath=f'../r3dresults/{modelname}/amr_grid.inp',
    gridpath=f'../r3dresults/{modelname}/grid_distances.csv',
    sizepath=f'../r3dresults/{modelname}/grid_cellsizes.csv',
    Nspecies=1,
    monomermasses=[2.3362e-22]
)
#    RETURNS
#      dust_density_dust_{phase}.inp
#      dust_temperature_dust_{phase}.dat
#      dustopac_dust_{phase}.inp


# Move them to correct folder
os.system(f'mv ../dust_density_dust_{phase}.inp ../r3dresults/{modelname}/{phase}/dust_density_dust.inp')
os.system(f'mv ../dust_temperature_dust_{phase}.dat ../r3dresults/{modelname}/{phase}/dust_temperature_dust.dat')
os.system(f'mv ../dustopac_dust_{phase}.inp ../r3dresults/{modelname}/{phase}/dustopac_dust.inp')

# NOTE
# Resulting files are, for each ../r3dresults/[modelname]/[phase]/
#
# dust_density_dust.inp
# dust_temperature_dust.dat
# dustopac_dust.inp
