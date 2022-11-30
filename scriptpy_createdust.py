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

# Bin grain sizes to 10 bins linearly
a5d.bin_grainsizes(
    grainsizepath='../grain_sizes_186.dat',
    phase=phase,
    nbins=10,
    lin='y'
)
# Output is ../grain_sizes_binned_{phase}.dat

# Write optool-script to create all dust species
#
# NOTE
# The info in co5bold-data uses 'forsterite' as species
# while here and for optool I use the chemical composition
# automate this somehow? For now leave as is, I can check species
# for future co5bold-data manually and change this manually
#
c3d.create_optoolscript(
    wavelength_path=f'{path}wavelength_micron.inp',
    phase=phase,
    grainum_sizes=f'../grain_sizes_binned_{phase}.dat',
    grainsize_type='normal',
    grainsize_na=21,
    specie='mg2sio4',
    grain_type='mie'
)

# Move optool script to phase-folder
os.system(f'mv ../optool_script_{phase}.sh ../r3dresults/{modelname}/{phase}/optool_script.sh')


# Then create dust files 
a5d.create_dustfiles(
    savpath=f'../co5bold_data/d{modelname}/{modelname}_{phase}.sav',
    amrpath=f'../r3dresults/{modelname}/amr_grid.inp',
    gridpath=f'../r3dresults/{modelname}/grid_distances.csv',
    sizepath=f'../r3dresults/{modelname}/grid_cellsizes.csv',
    grainsizepath='../grain_sizes_binned_{phase}.dat',
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
