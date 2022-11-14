# For usage with Bash-scripting
# Translates star data from c5d to r3d with help of bash-inputs
#
# > python3 scriptpy_creategrid.py $MODELNAME $PHASE_N
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

# ------------------------------------------------------------------------
# Python code below

import analyze_co5bold_functions as a5d
import analyze_r3d_functions as a3d

path = f'../r3dresults/{modelname}/'
AUcm = 1.49598e13 # cm

# Extract star's radius
print(f'\nExtracting {modelname}_{phase} information.')
Mstar,starradius,Lstar = a5d.load_star_information(
    savpath = f'../co5bold_data/d{modelname}/{modelname}_{phase}.sav',
    printoutput = 'n'
)

# Extract griddistances
griddistances = a3d.load_griddistances(
    gridpath= f'{path}grid_distances.csv',
    amrpath= f'{path}amr_grid.inp'
)[:,0]

# Create all r3d-data from c5d-data
# Must run a5d.heavydata first, ie scriptpy_createnpy.py [modelname] [phase]
a5d.create_stars(
    modelnames = [f'{modelname}'],
    phases = [phase]
)

# Remove negative spikes in opacity and density to correct final luminosity
a5d.smooth_stellardata(
    path = path,
    phases = [phase],
    starradii = [starradius],
    griddistances = griddistances,
    clean_data = 'y'
)

# NOTE
# Resulting files are, for each ../r3dresults/[modelname]/[phase]/
#
#   dust_temperature_onestar_smoothed.dat
#   star_opacities_smoothed.dat
#   dust_density_opastar.inp
#   dustkappa_opastar.inp
#   dustopac_opastar.inp

