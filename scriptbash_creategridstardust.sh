#! /bin/bash
#
# Use this script to create all r3d-data you need from c5d-data
# Change modelnames and phases and it goes automatically through all
# necessary python scripts.
#
# Define variables
modelname=st28gm06n052
phase0=186
phase1=190
phase2=198

# Create R3D-grid
python3 scriptpy_creategrid.py $modelname $phase0 $phase1 $phase2
wait

# Create stars
python3 scriptpy_createstar.py $modelname $phase0 &
python3 scriptpy_createstar.py $modelname $phase1 &
python3 scriptpy_createstar.py $modelname $phase2 &

wait # remove this wait if I'm doing this on servers!
# Extract and create dust-files
python3 scriptpy_createdust.py $modelname $phase0 &
python3 scriptpy_createdust.py $modelname $phase1 &
python3 scriptpy_createdust.py $modelname $phase2 &

# Write r3d-runscripts
wait
python3 scriptpy_write_r3d_runscripts.py $modelname $phase0 $phase1 $phase2




