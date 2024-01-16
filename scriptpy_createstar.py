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
savpath = f'../co5bold_data/d{modelname}/{modelname}_{phase}.sav'
AUcm = 1.49598e13 # cm

# Extract star's radius
print(f'\nExtracting {modelname}_{phase} information.')
Mstar,starradius,Lstar = a5d.load_star_information(
    savpath = savpath,
    printoutput = 'n'
)

# Extract griddistances
griddistances = a3d.load_griddistances(
    gridpath= f'{path}grid_distances.csv',
    amrpath= f'{path}amr_grid.inp'
)[:,0]


# Create all r3d-data from c5d-data
# Must run a5d.heavydata first, ie scriptpy_createnpy.py [modelname] [phase]
a5d.create_star(
    savpath = savpath,
    npypath = '../',
    amrpath = f'{path}amr_grid.inp',
    gridpath = f'{path}grid_distances.csv',
    sizepath = f'{path}grid_cellsizes.csv'
)
# OUTPUT FILES
#    R3D density file: dust_density_onestar_{phase}.inp
#    R3D temperature file: dust_temperature_onestar_{phase}.dat
#    Useful file with list of extracted opacities 'kappa': star_opacities_{phase}.dat
#
#
# Move output files to specific folders
print('  Moving files\n')
# Files from create star
os.system(
    f'mv ../dust_density_onestar_{phase}.inp {path}{phase}/dust_density_onestar.inp'
)
os.system(
    f'mv ../dust_temperature_onestar_{phase}.dat {path}{phase}/dust_temperature_onestar.dat'
)
os.system(
    f'mv ../star_opacities_{phase}.dat {path}{phase}/star_opacities.dat'
)
# Copy files from create grid-part above to phase-folders
os.system(
    f'cp {path}amr_grid.inp {path}{phase}/'
)
os.system(
    f'cp {path}wavelength_micron.inp {path}{phase}/'
)


# Reduce core temperature of the star to compensate for high liminosity when
# using monte carlo radiative transfer
a5d.reduce_coretemperature(
    savpath = savpath,
    temperaturepath = f'{path}{phase}/dust_temperature_onestar.dat,
    gridpath = f'{path}grid_distances.csv',
    amrpath = f'{path}amr_grid.inp',
    Rin = 0.9,
    Tin = 1000
)
# OUTPUT
#   dust_temperature_{phase}_{Rin}Rstar_{Tin}K.dat
#
# Move temperature-file
os.system(
    f'mv ../dust_temperature_{phase}_0.9Rstar_1000K.dat {path}{phase}/'
)


# Combine gas density and gas opacity into gas-optical-thickness, included as a
# dust_density-file in the R3D-forman.
a5d.create_staropadensity(
    pathopacity = f'{path}{phase}/star_opacities.dat',
    pathstardensity = f'{path}{phase}/dust_density_onestar.inp',
    pathwavelength = f'{path}wavelength_micron.inp',
    phase = phase,
    corrfac = 1.0
)
# OUTPUT
#
#  dust_density_opastar_{phase}.inp
#  dustopac_star_{phase}.inp
#  dustkappa_opastar_{phase}.inp
#
# Move final files to work/model folder
os.system(
    f'mv ../dust_density_opastar_{phase}.inp {path}{phase}/dust_density_opastar.inp'
)
os.system(
    f'mv ../dustopac_star_{phase}.inp {path}{phase}/dustopac_opastar.inp'
)
os.system(
    f'mv ../dustkappa_opastar_{phase}.inp {path}{phase}/dustkappa_opastar.inp'
)

# NOTE
# Resulting files are, for each ../r3dresults/[modelname]/[phase]/
#
#   dust_density_opastar.inp
#   dustkappa_opastar.inp
#   dustopac_opastar.inp
#   ./dust_temperature_{phase}_0.9Rstar_1000K.dat
#   star_opacities.dat


