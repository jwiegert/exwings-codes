#! /bin/bash
#
# Use this script to create all r3d-data you need from c5d-data
# Change modelnames and phases and it goes automatically through all
# necessary python scripts.

python3 scriptpy_creategrid.py st28gm06n052 186 190 198
wait

# create star-scripts below
python3 scriptpy_createstar.py st28gm06n052 186 &
python3 scriptpy_createstar.py st28gm06n052 190 &
python3 scriptpy_createstar.py st28gm06n052 198 &

# Write r3d-runscripts afterwards
wait
python3 scriptpy_write_r3d_runscripts.py st28gm06n052 186 190 198

