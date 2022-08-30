#! /bin/bash

python3 scriptpy_creategrid.py st28gm06n056 140 141 142
wait

&
# create star-scripts here

python3 scriptpy_createstar_smoothing.py st28gm06n056 140

&

python3 scriptpy_createstar_smoothing.py st28gm06n056 141

&

python3 scriptpy_createstar_smoothing.py st28gm06n056 142

&
wait

# Create radmc3d-run scripts # Is this possible to do like this?
python3 c3d.write_r3d_runscripts(
    path = '../r3dresults/st28gm06n056/',
    phase_list = [140,141,142],
    sed_inclination_list = [0,90,180,270],
    image_wavelength_list = [10],
    image_inclination_list = [90],
    image_sizeau = 7.4,
    image_npix = 128,
)
