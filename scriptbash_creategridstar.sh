#! /bin/bash
python3 scriptpy_creategrid.py st28gm06n056 140 141 142
wait
# create star-scripts below

# TODO change create-star-functions so that everything is not written in the same
# files!!

python3 scriptpy_createstar.py st28gm06n056 140 &
python3 scriptpy_createstar.py st28gm06n056 141 &
python3 scriptpy_createstar.py st28gm06n056 142 &

# Write r3d-runscripts afterwards
wait
python3 scriptpy_write_r3d_runscripts.py st28gm06n056 140 141 142
