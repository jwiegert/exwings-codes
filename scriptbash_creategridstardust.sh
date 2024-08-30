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
# $ ./scriptbash_creategridstardust.sh 2> logerror.txt | cat > logoutput.txt
#

# Define variables: model and phases
modelname=st28gm06n052
#modelname=st28gm06n074
phase0=245

# Basic settings (for not including dust - used scriptpy_creategridstar.sh
include_dust=yes
modify_Tdust=yes

# Create R3D-grid
# One line for all phases!
python3 scriptpy_creategrid.py $modelname $phase0 &

# Extract and create temporary npy-files for the remaining steps
#   gas & dust densities, temperatures, opacity
# One line per phase
python3 scriptpy_createnpy.py $modelname $phase0 $include_dust &

wait

# Extract and create star-files
# One line per phase
python3 scriptpy_createstar.py $modelname $phase0 &

# Extract and create dust-files
# modify_Tdust -> yes if normalise dust temperature to T(R)-Bladh2012
# One line per phase
python3 scriptpy_createdust.py $modelname $phase0 $modify_Tdust &


wait


# Remove all npy-files
rm ../*.npy

# and merge star and dust density and temperature data
# One line per phase
python3 scriptpy_mergedata.py $modelname $phase0


# Create dust opacity files and merge them with gas opacity
# One at a time to not mix the opacity files!
# One line per phase
python3 scriptpy_createopacity.py $modelname $phase0

wait
echo 'All done, press enter to finish.'

