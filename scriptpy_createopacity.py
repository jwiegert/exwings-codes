# For usage with Bash-scripting
# extracts dust and creates r3d-dust-files based on bash inputs
#
# > python3 scriptpy_createopacity.py $MODELNAME $PHASE1 ... $PHASE{N}
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
print(f'Running scriptpy_createopacity for {phase}.\n')

import create_r3d_functions as c3d
import os

# Write optool-script to create all dust species
#
# NOTE
# The info in co5bold-data uses 'forsterite' as species
# while here and for optool I use the chemical composition
# automate this somehow? For now leave as is, I can check species
# for future co5bold-data manually and change this manually
#
# > this is why I write dustopac-file here instead of in a5d.create_dustfiles()
#
# Declare dust specie
specie='mg2sio4'

path=f'../r3dresults/{modelname}/'

c3d.create_optoolscript(
    wavelength_path=f'{path}wavelength_micron.inp',
    phase=phase,
    grainum_sizes=f'../grain_sizes_binned_{phase}.dat',
    grainsize_type='normal',
    grainsize_na=21,
    specie=specie,
    grain_type='dhs'
)

# Move dustopac-file
os.system(f'mv ../dustopac_{specie}_{phase}.inp {path}{phase}/dustopac_dust.inp')

# Merge with opastar dustopac-file
c3d.merge_dustopac(
    workpath = f'{path}{phase}/',
    filenames = ['dustopac_opastar.inp','dustopac_dust.inp']
)


# Move optool script to phase-folder
os.system(f'mv ../optool_script_{phase}.sh {path}{phase}/optool_script.sh')

# Run script
os.system(f'{path}{phase}/optool_script.sh')

# Check results
print('\nList of all opacity files:\n')
os.system(f'ls -1 *{specie}*')

# Move results
os.system(f'mv *mg2sio4* {path}{phase}/')

# Move grain size files to folders as well
os.system(f'mv ../grain_sizes_binned_{phase}.dat {path}{phase}/')

# Resulting files:
#    dustkapscatmat_{specie}_{grain size in um}.inp
#    dustopac.inp
#    optool_script.sh

