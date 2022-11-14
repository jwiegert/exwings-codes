#! /bin/bash
#
# Use this script to create all r3d-data you need from c5d-data
# Change modelnames and phases and it goes automatically through all
# necessary python scripts.
#
# Define variables
modelname=st28gm06n056
phase0=140
phase1=141
phase2=142

# Create R3D-grid
python3 scriptpy_creategrid.py $modelname $phase0 $phase1 $phase2
wait

# Extract and create temporary npy-files for the remaining steps
#   gas & dust densities, temperatures, opacity
python3 scriptpy_createnpy.py $modelname $phase0 &
python3 scriptpy_createnpy.py $modelname $phase1 &
python3 scriptpy_createnpy.py $modelname $phase2 &
wait

# Create stars
python3 scriptpy_createstar.py $modelname $phase0 &
python3 scriptpy_createstar.py $modelname $phase1 &
python3 scriptpy_createstar.py $modelname $phase2 &
wait

# Remove all npy-files
rm ../*.npy &

# Write r3d-runscripts
python3 scriptpy_write_r3d_runscripts.py $modelname $phase0 $phase1 $phase2

