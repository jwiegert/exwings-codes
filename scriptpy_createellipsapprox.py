import analyze_ellipsoidapprox as ael
import analyze_co5bold_functions as a5d
import sys
import os


# INSTRUCTIONS
#
# > python3 scriptpy_createellipsapprox.py $picklepath $workpath






# Path to pickle-file
picklepath = sys.argv[1]

# Path ro work-folder
workpath = sys.argv[2]
# Automatically add / to end of path if it's missing
if workpath[-1] != '/':
    workpath += '/'


# Paths to grid files
amrpath = f'{workpath}amr_grid.inp'
gridpath = f'{workpath}grid_distances.csv'
sizepath = f'{workpath}grid_cellsizes.csv'
# e.g.
#    amrpath = '../r3dresults/st28gm06n052_arief_tests/amr_grid.inp',
#    gridpath = '../r3dresults/st28gm06n052_arief_tests/grid_distances.csv',
#    sizepath = '../r3dresults/st28gm06n052_arief_tests/grid_cellsizes.csv',


# FIRST STEP
#
# Create initial R3D-files dust_density, dust_temperature and grain sizes per grid cell
ael.create_dustapproximation(
    picklepath = picklepath,
    amrpath = amrpath,
    gridpath = gridpath,
    sizepath = sizepath,
    monomermass = 2.3362e-22,
    Amon = 2.3362e-22,
    rhomon = 3.27,
    ndnH = 3e-16,
    mH = 1.6736e-24,
    epsilonHe = 0.1
)
# RETURNS
#   dust_density_approx.inp
#   dust_temperature_approx.dat
#   grainsizes_approx.dat
#
# Move files to work folder
os.system(f'mv ../dust_density_approx.inp {workpath}')
os.system(f'mv ../dust_temperature_approx.dat {workpath}')
os.system(f'mv ../grainsizes_approx.dat {workpath}')


# BIN GRAIN SIZES
#
# for simplicity I just put everything in the same folder, sometimes I have final data in a subfolder
Nbins = 10
a5d.bin_grainsizes(
    grainsizepath = f'{workpath}grain_sizes_approx.dat',
    phase = 'approx',
    nbins = Nbins,
    lin = 'y'
)
# RETURNS
#   grain_sizes_binned_{phase}.dat
#
# Move file to work folder
os.system(f'mv ../grain_sizes_binned_approx.dat {workpath}')


# Create all r3d-input data according to binned grain sizes
ael.bin_inpdata(
    grainsizes_path = f'{workpath}grain_sizes_binned_approx.dat',
    density_path = f'{workpath}dust_density_approx.inp',
    temperature_path = f'{workpath}dust_temperature_approx.dat',
    wavelength_path = f'{workpath}wavelength_micron.inp'
)
# RETURNS
#   dust_density_approx_{Nbins}bins.inp
#   dust_temperature_approx_{Nbins}bins.dat
#   optool_script_approx.sh
#   dustopac_mg2sio4_approx.inp
#
os.system(f'mv ../dust_density_approx_{Nbins}bins.inp {workpath}')
os.system(f'mv ../dust_temperature_approx_{Nbins}bins.dat {workpath}')
os.system(f'mv ../optool_script_approx.sh {workpath}')
os.system(f'mv ../dustopac_mg2sio4_approx.inp {workpath}')



# densities to merge:
# dust_density_approx_{Nbins}bins.inp
# and a c5d-star-file (opastar etc)

# temperatures to merge:
# dust_temperature_approx_{Nbins}bins.dat
# and a c5d-star-temperature

# TODO merge with a chosen stellar gas-model
#c3d.merge_dustdensities(
#    workpath:str = '../r3dresults/',
#    filenames=[
#        'dust_density_star.inp',
#        'dust_density_dust.inp'
#    ],
#)



#c3d.merge_dusttemperatures(
#    filenames=[
#        '../../../r3dsims/co5bold_tests/05-lowergridres/unmerged/dust_temperature_star.dat',
#        '../../../r3dsims/co5bold_tests/05-lowergridres/unmerged/dust_temperature_dust.dat'
#    ],
#    modelname='st28gm06n052_staranddust_1',
#    phases=[0],
#)

# TODO
# move the merged files too


# print some output, like "don, everything is in workpath, clean up and order your files"