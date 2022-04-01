# exwings-codes

Source codes for handling post processing radiative transfer simulation source codes for Exwings.

I.e. translate co5bold simulations into r3dmc3d-input data and simulate synthetic observations.

## create_r3d_functions

### create_grid()
```py
def create_grid(basecubesize:float, nxyz:int, refinementlist:list, savegrid:str='y'):
"""
    Creates grid for Radmc3D simulations and related informaton 
    files used for analysis and creating dust envelopes. Grid is 
    octree cubic. It is refined (higher resolution) closer to 
    the central of the grid. A maximum of four (4) levels of 
    refinements is allowed in this version of the function.
    
    INPUTS
    ------
    basecubesize: length of side of base cells in AU (these are cubes) [int or float]
    
    nxyz: number of base cells along one side of the whole grid [even number, int]
    
    refinementlist: list of radial distances in AU to each level of refinement [float,float], no more than 4 numbers!
    
    savegrid: default set to 'y' [str]. If not 'y', then no grid_distances.csv or grid_cellsizes.csv will be saved. These are useful for analysing inputs and outputs of R3D!

    OUTPUTS
    -------
    amr_grid.inp : grid file for R3D-simulations

    Optional: 
    grid_distances.csv : array of radial, x, y, and z distances to each grid cell in cm
    grid_cellsizes.csv : array of sizes of each grid cell in cm
    (both have same order as in dust_density.inp and dust_temperature)
    """
```


## plot_r3d_functions

```py
def plot_grid(
        gridpath:str='../grid_distances.csv',
        sizepath:str='../grid_cellsizes.csv',
        amrpath:str='../amr_grid.inp',
        nbins:int=10
    ):
    """
    Loads and plots the current grid.

    Inputs
    ------
    gridpath: path to grid_distances.csv
    amrpath: path to amr_grid.inp
    nbins: number of bins for radial histogram plot
    """
```

# Running R3D

## Self checks

- Check your grid and density distribution using vtk-files

Even though you don't run any radiative transfer here you still need to formulate a radmc3d.inp file.

For grid only run

> radmc3d vtk_grid

For density run, where 1 is the dust specie you want to visualise.

> radmc3d vtk_dust_density 1

For dust temperature run

> radmc3d vtk_dust_temperature 1

- Check vtk-data with paraview

> paraview FILENAME.vtk

Press apply. For just a grid, chose e.g. Representation Wireframe.

For a density file, press apply, and chose "threshold" (button, top left area). Chose some limit (increase minimum) and look around. Remember that R3D sets som 1e-80 ish g/cm3 i all grid cells since it can't be zero. So you have to increase the lower limit to see past all zero-density cells.

# Cython notes

### Pure Python with Cython

Examples:

Locals-decorator declares class of input variables

cython.declare to declare class of variables inside function

```py
import cython
@cython.cfunc
@cython.locals(
    c5dx = cython.double,
    c5dy = cython.double,
    c5dz = cython.double
)
def function_name(c5dx,c5dy,c5dz):
    nr = cython.declare(cython.int,0)
    foovariable1 = cython.declare(cython.double)
    foovariable2 = cython.declare(cython.float)

    # Here be dragons

    return 'something'
```

### Links

http://docs.cython.org/en/latest/src/tutorial/pure.html

https://cython.readthedocs.io/en/latest/src/userguide/language_basics.html

https://www.infoworld.com/article/3648539/

faster-python-made-easier-with-cythons-pure-python-mode.html

https://pandas.pydata.org/pandas-docs/stable/user_guide/enhancingperf.html

