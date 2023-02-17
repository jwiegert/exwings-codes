# exwings-codes

Source codes for handling post processing radiative transfer simulation source codes for Exwings with Radmc-3D.  I.e. translate co5bold simulations into r3dmc3d-input data and simulate synthetic observations.

* TODO needs rewriting and updating from here!

# Running R3D

## Self checks

- Check your grid and density distribution using vtk-files

Even though you don't run any radiative transfer here you still need to have a radmc3d.inp file and kappa-files.

For grid only run

> radmc3d vtk_grid

To see dust density, where 1 is the dust specie you want to visualise.

> radmc3d vtk_dust_density 1

To see dust temperature run

> radmc3d vtk_dust_temperature 1

- Check vtk-data with paraview

> paraview FILENAME.vtk

Press apply. For just a grid, chose e.g. Representation Wireframe.

For a density file, press apply, and chose "threshold" (button, top left area). Chose some limit (increase minimum) and look around. Remember that R3D sets som 1e-80 ish g/cm3 i all grid cells since it can't be zero. So you have to increase the lower limit to see past all zero-density cells.

### Links

http://docs.cython.org/en/latest/src/tutorial/pure.html

https://cython.readthedocs.io/en/latest/src/userguide/language_basics.html

https://www.infoworld.com/article/3648539/

faster-python-made-easier-with-cythons-pure-python-mode.html

https://pandas.pydata.org/pandas-docs/stable/user_guide/enhancingperf.html

