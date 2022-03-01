# exwings-codes

SUBMIT WITH NEW COMPUTER?


Source codes for handling post processing radiative transfer simulation source codes for Exwings.

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

For density run

> radmc3d vtk_dust_density 1

Where 1 is the dust specie you want to visualise.

- Check vtk-data with paraview

> paraview FILENAME.vtk

Press apply. For just a grid, chose e.g. Representation Wireframe.

For a density file, press apply, and chose "threshold" (button, top left area). Chose some limit (increase minimum) and look around. Remember that R3D sets som 1e-80 ish g/cm3 i all grid cells since it can't be zero. So you have to increase the lower limit to see past all zero-density cells.