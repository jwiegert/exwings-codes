# For usage with Bash-scripting
# creates various run-r3d-scripts
#
# > python3 scriptpy_write_r3d_runscripts.py $MODELNAME $PHASE1 ... $PHASE{N}
#
# And from data written in here, change codes manually if required!
#
# ------------------------------------------------------------------------
# Include inputs from bash

import create_r3d_functions as c3d
import sys

modelname = sys.argv[1]
phases = [phase for nn,phase in enumerate(sys.argv) if nn > 1]

c3d.write_r3d_runscripts(
    path = f'../r3dresults/{modelname}/',
    phase_list = phases,
    sed_inclination_list = [0,90,180,270],
    image_wavelength_list = [1,10,100],
    image_inclination_list = [0],
    image_sizeau = 30,
    image_npix = 256,
)

