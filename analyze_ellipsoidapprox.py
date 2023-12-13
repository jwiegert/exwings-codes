# Loop through coords of one of the ellipsoids and save in the r3d-grid
# In cm
import numpy as np
import pickle
import analyze_r3d_functions as a3d
import os

AUcm = 1.49598e13 # AU in cm cm
monomermass = 2.3362e-22 # Forsterite mass in g
testrun = 2



############################################################################################

if testrun == 2:

    print(f'Start time:')
    os.system('date')
    print('Loading data')

    # Load Ariefs pickles
    path = '../../exwings_archivedata/co5bold_data/st28gm06n052_arief_tests/filled-mean-ellipsoids_st28gm06n052-032.pickle'
    with open(path, 'rb') as f:
        ellipsoid_dict = pickle.load(f)


    # Load an r3d-grid
    r3ddistances = a3d.load_griddistances(
        gridpath = '../r3dresults/st28gm06n052_arief_tests/grid_distances.csv',
        amrpath = '../r3dresults/st28gm06n052_arief_tests/amr_grid.inp',
    )
    r3dcellsizes = a3d.load_cellsizes(
        amrpath = '../r3dresults/st28gm06n052_arief_tests/amr_grid.inp',
        sizepath = '../r3dresults/st28gm06n052_arief_tests/grid_cellsizes.csv',
    )
    nleafs = a3d.load_grid_properties(
        amrpath = '../r3dresults/st28gm06n052_arief_tests/amr_grid.inp'
    )[2]

    # Declare arrays
    # This was faster than looping through whole grid, writing zeros, again and again
    r3d_densities = np.zeros(nleafs)
    r3d_temperatures = np.zeros(nleafs)
    r3d_condenfrac = np.zeros(nleafs)
    r3d_mincellsize = r3dcellsizes.min()



    Nellipsoids = ellipsoid_dict['distance_centre_grid'].size
    print(f'  There are {Nellipsoids} dust clouds:')


    # Loop through each cloud
    for nellipsoid in range(Nellipsoids):

        Nc5d = ellipsoid_dict['coord_d_list'][nellipsoid][:,0].size
        print(f'  Cloud {nellipsoid+1} has {Nc5d} c5d-cells.')

        cloud_coordsX = np.array(ellipsoid_dict['coord_d_list'][nellipsoid][:,0])
        cloud_coordsY = np.array(ellipsoid_dict['coord_d_list'][nellipsoid][:,1])
        cloud_coordsZ = np.array(ellipsoid_dict['coord_d_list'][nellipsoid][:,2])

        cloud_densities = np.array(ellipsoid_dict['filled_nmonomer'][nellipsoid][:])

        # Extract cellindeces for r3d-cells around/in the cloud
        nr3drange = np.where(
            (r3ddistances[:,1] >= cloud_coordsX.min()) & \
            (r3ddistances[:,1] <= cloud_coordsX.max()) & \
            (r3ddistances[:,2] >= cloud_coordsY.min()) & \
            (r3ddistances[:,2] <= cloud_coordsY.max()) & \
            (r3ddistances[:,3] >= cloud_coordsZ.min()) & \
            (r3ddistances[:,3] <= cloud_coordsZ.max())
        )[0]

        print(f'      and {nr3drange.size} r3d-cells.')

        # Loop through the r3d-cells in cloud only
        # Extract c5d-cells with dust in them and
        # save at correct places
        for nr3d in nr3drange:

            # Check if nr3d-cell is in any c5d-cells:
            cloudcell = np.where(
                (cloud_coordsX >= r3ddistances[nr3d,1]-0.5*r3dcellsizes[nr3d]) & \
                (cloud_coordsX <= r3ddistances[nr3d,1]+0.5*r3dcellsizes[nr3d]) & \
                (cloud_coordsY >= r3ddistances[nr3d,2]-0.5*r3dcellsizes[nr3d]) & \
                (cloud_coordsY <= r3ddistances[nr3d,2]+0.5*r3dcellsizes[nr3d]) & \
                (cloud_coordsZ >= r3ddistances[nr3d,3]-0.5*r3dcellsizes[nr3d]) & \
                (cloud_coordsZ <= r3ddistances[nr3d,3]+0.5*r3dcellsizes[nr3d])
            )[0]

            if len(cloudcell) > 0:
                r3d_densities[nr3d] = cloud_densities[cloudcell].mean() * monomermass



    # Write density file
    print('  Prints densities')

    with open(f'../dust_density_arief032_long.inp', 'w') as fdensity:

        # Write headers:
        #
        # Density:
        # 1
        # nleafs
        # number dust species
        fdensity.write(f'1\n{int(nleafs)}\n1\n')

        for density in r3d_densities:
            fdensity.write(f'{density}\n')

    print('Done at:')
    os.system('date')
