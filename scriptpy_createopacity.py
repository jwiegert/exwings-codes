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

c3d.create_optoolscript(
    wavelength_path=f'../r3dresults/opacities_{modelname}/wavelength_micron.inp',
    phase=phase,
    grainum_sizes=f'../grain_sizes_binned_{phase}.dat',
    grainsize_type='normal',
    grainsize_na=21,
    specie=specie,
    grain_type='mie'
)

# Move dustopac-file
os.system(f'mv ../dustopac_{specie}_{phase}.inp ../r3dresults/{modelname}/{phase}/dustopac_dust.inp')

# Merge with opastar dustopac-file
c3d.merge_dustopac(
    filenames=['dustopac_opastar.inp','dustopac_dust.inp'],
    modelname=modelname,
    phases=[phase],
)


# Move optool script to phase-folder
os.system(f'mv ../optool_script_{phase}.sh ../r3dresults/{modelname}/{phase}/optool_script.sh')

# Run script
os.system(f'../r3dresults/{modelname}/{phase}/optool_script.sh')

# Check results
print('\nList of all opacity files:\n')
os.system('ls -1 *mg2sio4*')

# Move results
os.system(f'mv *mg2sio4* ../r3dresults/{modelname}/{phase}/')

# Clean up?
#os.system(f'rm ../grain_sizes_binned_{phase}.dat')

# Resulting files:
#    dustkapscatmat_{specie}_{grain size in um}.inp
#    dustopac.inp
#    optool_script.sh
