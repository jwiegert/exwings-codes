# For usage with Bash-scripting
# merges all density, temperature and opac-files into one file for each phase
#
# > python3 scriptpy_mergedata.py $MODELNAME $PHASE1 ... $PHASE{N}
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

# Then merge them with stellar files
c3d.merge_dustdensities(
    filenames=['dust_density_opastar.inp','dust_density_dust.inp'],
    modelname=modelname,
    phases=[phase],
)

c3d.merge_dusttemperatures(
    filenames=['dust_temperature_onestar_smoothed.dat','dust_temperature_dust.dat'],
    modelname=modelname,
    phases=[phase],
)

c3d.merge_dustopac(
    filenames=['dustopac_opastar.inp','dustopac_dust.inp'],
    modelname=modelname,
    phases=[phase],
)

