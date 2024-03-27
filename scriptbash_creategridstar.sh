#! /bin/bash
#
# Use this script to create all r3d-data you need from c5d-data
# Change modelnames and phases and it goes automatically through all
# necessary python scripts.
#
# Define variables

modelname=st28gm06n052
#phase1=032
#phase0=186
phase1=190
#phase2=198

#modelname=st28gm06n056
#phase1=140

# Create R3D-grid
#python3 scriptpy_creategrid.py $modelname $phase0 $phase1 $phase2
python3 scriptpy_creategrid.py $modelname $phase1
wait

# Extract and create temporary npy-files for the remaining steps
#   gas densities, temperatures, opacity
include_dust=no
#python3 scriptpy_createnpy.py $modelname $phase0 $include_dust &
python3 scriptpy_createnpy.py $modelname $phase1 $include_dust &
#python3 scriptpy_createnpy.py $modelname $phase2 $include_dust &
wait

# Create stars
#python3 scriptpy_createstar.py $modelname $phase0 &
python3 scriptpy_createstar.py $modelname $phase1 &
#python3 scriptpy_createstar.py $modelname $phase2 &
wait

# Remove all npy-files
rm ../*.npy &

# Write r3d-runscripts
#python3 scriptpy_write_r3d_runscripts.py $modelname $phase0 $phase1 $phase2

wait
echo 'All done, press enter to finish.'

