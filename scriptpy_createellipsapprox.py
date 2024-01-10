import analyze_ellipsoidapprox as ael
import analyze_co5bold_functions as a5d
import create_r3d_functions as c3d
import sys
import os


# INSTRUCTIONS
#
# > python3 scriptpy_createellipsapprox.py $picklepath $workpath
#    picklepath: where to get ellisoidal dust data from
#    workpath: where everything should end up and where co5bold-star-dust-data are gathered
#
# Example:
# > python3 scriptpy_createellipsapprox.py '../../exwings_archivedata/co5bold_data/st28gm06n052_arief_tests/filled-radialdistribution-ellipsoids_st28gm06n052-032.pickle' '../r3dresults/st28gm06n052_arief_tests/032_186_radial_combine/' 


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


# MERGE DATA
#
# Merge ellips-dust-data with a stellar model from co5bold


# densities to merge:
# dust_density_approx_{Nbins}bins.inp
# and a c5d-star-file (opastar etc)

# temperatures to merge:
# dust_temperature_approx_{Nbins}bins.dat
# and a c5d-star-temperature


c3d.merge_dustdensities(
    workpath = workpath,
    filenames=[
        'dust_density_c5d.inp',
        'dust_density_approx.inp'
    ],
)
c3d.merge_dusttemperatures(
    workpath = workpath,
    filenames=[
        '../../../r3dsims/co5bold_tests/05-lowergridres/unmerged/dust_temperature_star.dat',
        '../../../r3dsims/co5bold_tests/05-lowergridres/unmerged/dust_temperature_dust.dat'
    ],
)
# merge dustopac also
c3d.merge_dustopac(
    workpath = workpath,
    filenames = [
        'dustopac_c5d.inp',
        'dustopac_mg2sio4_approx.inp'
    ]
)
# NOTE
# Moves the merged files back to work folder

print(f'DONE:\n  Files moved to {workpath}. Clean up and order your files for radmc3d-runs.\n  Run optool-script to create all opacity data')
