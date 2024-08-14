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
import os

# Set up variables
modelname = sys.argv[1]
phases = [phase for nn,phase in enumerate(sys.argv) if nn > 1]

# Create folder to put files in
os.system(f'mkdir ../r3dsims/{modelname}/')

# Write scripts
c3d.write_r3d_runscripts(
    path = f'../r3dsims/{modelname}/',
    phase_list = phases,
    sed_angles_list = [
        [0,0],
        [90,0],
        [180,0],
        [270,0],
        [90,90],
        [90,270],
    ],
    image_wavelength_list = [1,10],
    image_angles_list = [
        [0,0],
        [90,0],
        [180,0],
        [270,0],
        [90,90],
        [90,270],
    ],
    image_sizeau = 30,
    image_npix = 512,
    Nscripts = 1
)






