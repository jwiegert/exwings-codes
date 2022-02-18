# exwings-codes

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

