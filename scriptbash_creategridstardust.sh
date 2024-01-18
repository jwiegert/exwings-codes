#! /bin/bash
#
# Use this script to create all r3d-data you need from c5d-data
# Change modelnames and phases and it goes automatically through all
# necessary python scripts.
#
# Many users, nice this:
# $ nice -n 19 ./scriptbash_creategridstardust.sh
#
# Save in a logfile? Run this (create log.txt first)
# $ ./scriptbash_creategridstardust.sh | cat > logoutput.txt cat 2> logerror.txt
#
# All of it?
# $ nice -n 19 ./scriptbash_creategridstardust.sh | cat > logoutput.txt cat 2> logerror.txt
#
# Define variables
modelname=st28gm06n052
#phase0=186
phase0=190
#phase2=198


# Create R3D-grid
#python3 scriptpy_creategrid.py $modelname $phase0 $phase1 $phase2 &
python3 scriptpy_creategrid.py $modelname $phase0 &

# Extract and create temporary npy-files for the remaining steps
#   gas & dust densities, temperatures, opacity
include_dust=yes
python3 scriptpy_createnpy.py $modelname $phase0 $include_dust &
#python3 scriptpy_createnpy.py $modelname $phase1 $include_dust &
#python3 scriptpy_createnpy.py $modelname $phase2 $include_dust &


wait


# Extract and create star-files
python3 scriptpy_createstar.py $modelname $phase0 &
#python3 scriptpy_createstar.py $modelname $phase1 &
#python3 scriptpy_createstar.py $modelname $phase2 &
# Extract and create dust-files    TODO adapt dust-temperature to theory-model!
python3 scriptpy_createdust.py $modelname $phase0 &
#python3 scriptpy_createdust.py $modelname $phase1 &
#python3 scriptpy_createdust.py $modelname $phase2 &


wait

# Remove all npy-files
rm ../*.npy &

# and merge star and dust density and temperature data
python3 scriptpy_mergedata.py $modelname $phase0 &
#python3 scriptpy_mergedata.py $modelname $phase1 &
#python3 scriptpy_mergedata.py $modelname $phase2 &


wait


# Create dust opacity files and merge them with gas opacity
python3 scriptpy_createopacity.py $modelname $phase0 &
#python3 scriptpy_createopacity.py $modelname $phase1 &
#python3 scriptpy_createopacity.py $modelname $phase2 &

# Write r3d-runscripts
#python3 scriptpy_write_r3d_runscripts.py $modelname $phase0 $phase1 $phase2


