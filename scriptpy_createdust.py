# For usage with Bash-scripting
# extracts dust and creates r3d-dust-files based on bash inputs
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

path = f'../r3dresults/{modelname}/'

# First extract grain sizes
a5d.extract_grainsizes(
    amrpath=f'{path}amr_grid.inp',
    gridpath=f'{path}grid_distances.csv',
    sizepath=f'{path}grid_cellsizes.csv',
    savpath=f'../co5bold_data/d{modelname}/{modelname}_{phase}.sav',
    Amon=2.3362e-22,
    rhomon=3.27,
    ndnH=3e-16,
    mH=1.6736e-24,
    epsilonHe=0.1
)
# Output is ../grain_sizes_{phase}.dat

# Bin grain sizes to 10 bins linearly
a5d.bin_grainsizes(
    grainsizepath=f'../grain_sizes_{phase}.dat',
    phase=phase,
    nbins=10,
    lin='y'
)
# Output is ../grain_sizes_binned_{phase}.dat


# Then create dust files 
a5d.create_dustfiles(
    savpath=f'../co5bold_data/d{modelname}/{modelname}_{phase}.sav',
    amrpath=f'../r3dresults/{modelname}/amr_grid.inp',
    gridpath=f'../r3dresults/{modelname}/grid_distances.csv',
    sizepath=f'../r3dresults/{modelname}/grid_cellsizes.csv',
    grainsizepath=f'../grain_sizes_binned_{phase}.dat',
    Nspecies=1,
    monomermasses=[2.3362e-22]
)
#    RETURNS
#      dust_density_dust_{phase}.inp
#      dust_temperature_dust_{phase}.dat


# Move them to correct folder
os.system(f'mv ../dust_density_dust_{phase}.inp ../r3dresults/{modelname}/{phase}/dust_density_dust.inp')
os.system(f'mv ../dust_temperature_dust_{phase}.dat ../r3dresults/{modelname}/{phase}/dust_temperature_dust.dat')


# Clean up
os.system(f'rm ../grain_sizes_{phase}.dat')
# Binned .dat is used for optool, save for that
#os.system(f'rm ../grain_sizes_binned_{phase}.dat')


# NOTE
# Resulting files are, for each ../r3dresults/[modelname]/[phase]/
#
# dust_density_dust.inp
# dust_temperature_dust.dat

