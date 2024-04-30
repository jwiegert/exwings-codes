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
modify_Tdust = sys.argv[3]

# ------------------------------------------------------------------------
# Pythoncodes below
print(f'Running scriptpy_createdust for {phase}.\n')

import analyze_co5bold_functions as a5d
import create_r3d_functions as c3d
import os

path = f'../r3dresults/{modelname}/'

# First extract grain sizes
a5d.extract_grainsizes(
    amrpath=f'{path}amr_grid.inp',
    gridpath=f'{path}grid_distances.csv',
    sizepath=f'{path}grid_cellsizes.csv',
    savpath=f'../../exwings_archivedata/co5bold_data/d{modelname}/{modelname}_{phase}.sav',
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
# The Nspecies is the number of chemical species of dust, not number of grain sizes
a5d.create_dustfiles(
    savpath=f'../../exwings_archivedata/co5bold_data/d{modelname}/{modelname}_{phase}.sav',
    amrpath=f'{path}amr_grid.inp',
    gridpath=f'{path}grid_distances.csv',
    sizepath=f'{path}grid_cellsizes.csv',
    grainsizepath=f'../grain_sizes_binned_{phase}.dat',
    Nspecies=1,
    monomermasses=[2.3362e-22]
)
#    RETURNS
#      dust_density_dust_{phase}.inp
#      dust_temperature_dust_{phase}.dat
# Move density file
os.system(f'mv ../dust_density_dust_{phase}.inp {path}{phase}/dust_density_dust.inp')


# Normalise radial profile of dust-temperatures with T(R)-Bladh2012
if modify_Tdust == 'yes':
    a5d.modify_dusttemperature(
        dusttemperature_path=f'../dust_temperature_dust_{phase}.dat',
        griddistance_path=f'{path}grid_distances.csv',
        sav_path=f'../../exwings_archivedata/co5bold_data/d{modelname}/{modelname}_{phase}.sav',
        gridinfo_path=f'{path}grid_info.txt',
        amr_path=f'{path}amr_grid.inp',
    )
    # Move modified temperature file
    os.system(f'mv ../dust_temperature_dust_modified_{phase}.dat {path}{phase}/dust_temperature_dust.dat')
    # And unmodified file with new name
    os.system(f'mv ../dust_temperature_dust_{phase}.dat {path}{phase}/dust_temperature_dust_unmodified.dat')

else:
    # Move unmodified temperature file
    os.system(f'mv ../dust_temperature_dust_{phase}.dat {path}{phase}/dust_temperature_dust.dat')


# Clean up
os.system(f'rm ../grain_sizes_{phase}.dat')
# Do not remove
# ../grain_sizes_binned_{phase}.dat 
# It is used for optool

# NOTE
# Resulting files are, for each ../r3dresults/[modelname]/[phase]/
#
# dust_density_dust.inp
# dust_temperature_dust.dat

