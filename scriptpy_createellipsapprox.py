import analyze_ellipsoidapprox as ael
import analyze_co5bold_functions as a5d
import os



# TODO
# obtain path-data to where data are to be saved
# (not movelnames and such, better with paths!)


# filled-mean-ellipsoids_st28gm06n052-032.pickle
# filled-radialdistribution-ellipsoids_st28gm06n052-032.pickle




ael.create_dustapproximation(
    picklepath = '../../exwings_archivedata/co5bold_data/st28gm06n052_arief_tests/filled-mean-ellipsoids_st28gm06n052-032.pickle'
)


# TODO move files


# Bin grain sizes

# 032_averages
# 032_radialdist

a5d.bin_grainsizes(
    grainsizepath = '../r3dresults/st28gm06n052_arief_tests/032_averages/grain_sizes_approx.dat',
    phase = 'approx'
)



# TODO move files


path = '../r3dresults/st28gm06n052_arief_tests/032_radialdist/'
ael.bin_inpdata(
    grainsizes_path = f'{path}grain_sizes_binned_approx.dat',
    density_path = f'{path}dust_density_approx.inp',
    temperature_path = f'{path}dust_temperature_approx.dat',
    wavelength_path = f'{path}wavelength_micron.inp'
)


# TODO move files


# TODO auto-run optool-script?



# TODO merge with a chosen stellar gas-model