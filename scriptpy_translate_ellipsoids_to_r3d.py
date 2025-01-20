import pickle
import numpy as np
import data_unpacker

# Load relevant pickle file
# i.e. either as obtained from CO5BOLD (clumps); or derived ellipsoids
#filename='/home/ariah/dat/caucus/clumps-ellipsoids-comparison/'+'st28gm06n052-derivedclumps_200.pickle'
#filename='/home/ariah/dat/caucus/clumps-ellipsoids-comparison/'+'st28gm06n052-derivedellipsoids_200.pickle'

filename='../arief_data/st28gm06n052-derivedclumps_200.pickle'

with open(filename, 'rb') as handle:
    dustclumps = pickle.load(handle)

# Use Ariefs data unpacker to create new lists of data
grid_filled_rho_gas=data_unpacker.data_unpacker(679, dustclumps['coord_list'], dustclumps['rho_gas_list'])
#grid_filled_rho_dust=data_unpacker(679, dustclumps['coord_list'], dustclumps['rho_dust_list'])
#grid_filled_temperature=data_unpacker(679, dustclumps['coord_list'], dustclumps['temperature_list'])
#grid_filled_grainsizes=data_unpacker(679, dustclumps['coord_list'], dustclumps['grainsizes_list'])



# Check details
print(np.shape(grid_filled_rho_gas))
#print(np.shape(dustclumps['coord_list']))





# Load C5D-grid
#    needs to create grid-info-files first with a5d.load_grid_properties()
#    outputs 
#      c5dgrid       coordinates, xyz of each cell centre
#      cellcourners  similar coordinates for courners of cells
#      cellsize      just size of minimum cell size


# figure out a way to loop through ariefs lists and c5d-grid



# Load R3D-grid



# Map data to R3D-grid and save r3d-input-files

# grid_filled_rho_gas -> gas_density-file
# grid_filled_rho_dust -> dust density file
# grid_filled_temperature -> gas temperature file
# grid_filled_grainsizes -> grain size file



# Move data to appropiate folders
