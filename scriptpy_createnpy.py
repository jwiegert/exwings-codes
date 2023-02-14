# For usage with Bash-scripting
# Loads alld heavy co5bold-data and saves them in temporary *.npy-files for
# useage later. This makes the whole process a lot fast and less RAM-demaning.
#
# > python3 scriptpy_createnpy.py $MODELNAME $PHASE_N
#
# And from data written in here, change codes manually if required!
#
# ------------------------------------------------------------------------
# Empty line for some space in terminal
print('')
# Include inputs from bash

import sys
modelname = sys.argv[1]
phase = sys.argv[2]
include_dust = sys.argv[3]

# ------------------------------------------------------------------------
# Python code below

import analyze_co5bold_functions as a5d

a5d.load_c5dheavydata(
   savpath = f'../co5bold_data/d{modelname}/{modelname}_{phase}.sav',
   Nspecies = 1,
   gas_density = True,
   dust_density = include_dust,
   gas_opacity = True,
   temperature = True
)

# NOTE
# Resulting files are, for each [modelname] & [phase]
#
#    c5ddust_density_[phase].npy
#    c5dgas_density_[phase].npy
#    c5dgas_opacity_[phase].npy
#    c5d_temperature_[phase].npy
#

