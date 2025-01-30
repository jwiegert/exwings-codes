import numpy as np


###########################################################################
# Ariefs function to create lists from the pickle-files
def data_unpacker(dim, coord_list, val_list):
    """
    Small function to unpack lists into grid;
    since saves a lot of memory to store clumps/ellipsoids as lists (rather than the eg. 679^3 grid).
    Available lists:
        coord_list - indices of points in the domain
        rho_gas_list - gas mass density in g/cm^3
        rho_dust_list - dust mass density in g/cm^3
        temperature_list - tenperature in K
        grainsizes_list - grain sizes in cm

    Example usage:

    # Load relevant pickle file
    # i.e. either as obtained from CO5BOLD (clumps); or derived ellipsoids
    filename='/home/ariah/dat/caucus/clumps-ellipsoids-comparison/'+'st28gm06n052-derivedclumps_200.pickle'  # ie on phy-exwings
    filename='/home/ariah/dat/caucus/clumps-ellipsoids-comparison/'+'st28gm06n052-derivedellipsoids_200.pickle'  # ie on phy-exwings
    with open(filename, 'rb') as handle:
        dustclumps = pickle.load(handle)

    grid_filled_rho_gas=data_unpacker(679, dustclumps['coord_list'], dustclumps['rho_gas_list'])
    grid_filled_rho_dust=data_unpacker(679, dustclumps['coord_list'], dustclumps['rho_dust_list'])
    grid_filled_temperature=data_unpacker(679, dustclumps['coord_list'], dustclumps['temperature_list'])
    grid_filled_grainsizes=data_unpacker(679, dustclumps['coord_list'], dustclumps['grainsizes_list'])

    The dimension supplied must match the dimension used in CO5BOLD (in this case 679^3).
    """
    array_shape = (dim, dim, dim)
    filled_arr = np.zeros(array_shape)

    for ii, _coords in enumerate(coord_list):
        _datain = np.hstack((_coords, val_list[ii].reshape(-1,1)))

        for coord in _datain:
            x, y, z, value = coord
            x=x.astype('int')
            y=y.astype('int')
            z=z.astype('int')
            filled_arr[x, y, z] += value
    return filled_arr

###########################################################################
# Function to extract and write dat-files with c5d-grid-coords
def write_c5dgrid_files(
        savpath:str='../../exwings_archivedata/co5bold_data/dst28gm06n052/st28gm06n052_150.sav',
        outpath:str='../',
    ):
    """
    Reads and prints dat-files with radial+xyz-coordinates of centre of CO5BOLD grid cells.
    Only one column, these coordinates repeat in all directions in the whole grid, so you
    can get all 3D-coords by repeated loops through these numbers.
    """
    import analyze_co5bold_functions as a5d

    c5dgrid,cellcourners,cellsize = a5d.load_grid_properties(
        savpath=savpath
    )
    Nside = np.shape(c5dgrid)[0]
    Ncells = Nside**3

    # Save and print c5dgrid

    with open(f'{outpath}gridc5d_distances.csv', 'w') as c5dfile:

        # Print header
        c5dfile.write('# Distances in cm to centre of CO5BOLD grid cells\n')
        c5dfile.write(f'# Total number of grid cells in CO5BOLD-grid: {Ncells}\n')
        c5dfile.write(f'# Number of cell along one axis: {Nside}\n')
        c5dfile.write('# --------\n')

        for ncell in range(Nside):
            c5dfile.write(f'{c5dgrid[ncell,0]}\n')

    # Print confirmation
    print(f'  write_c5dgrid_files: DONE: {outpath}gridc5d_distances.csv\n')

###########################################################################
# Func to merge final approx-dust-data with original star

def merge_final_data(
        workpath:str='../../r3dsims/st28gm06n052_arief/052_derivedclumps_150/',
    ):
    import create_r3d_functions as c3d
    import os

    c3d.merge_dustdensities(
        workpath = workpath,
        filenames=['dust_density_opastar.inp','dust_density_dust_binned.inp'],
    )
    c3d.merge_dusttemperatures(
        workpath = workpath,
        filenames=['dust_temperature_onestar_smoothed.dat','dust_temperature_dust_binned.dat'],
    )
    c3d.merge_dustopac(
        workpath = workpath,
        filenames = ['dustopac_opastar.inp','dustopac_dust.inp']
    )

