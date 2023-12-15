# Loop through coords of one of the ellipsoids and save in the r3d-grid
# In cm
import numpy as np
import pickle
import analyze_r3d_functions as a3d
import os

AUcm = 1.49598e13 # AU in cm cm
monomermass = 2.3362e-22 # Forsterite mass in g





############################################################################################


def create_dustapproximation(
        picklepath:str = '../../exwings_archivedata/co5bold_data/st28gm06n052_arief_tests/filled-mean-ellipsoids_st28gm06n052-032.pickle',
        amrpath:str = '../r3dresults/st28gm06n052_arief_tests/amr_grid.inp',
        gridpath:str = '../r3dresults/st28gm06n052_arief_tests/grid_distances.csv',
        sizepath:str = '../r3dresults/st28gm06n052_arief_tests/grid_cellsizes.csv',
        monomermass:float = 2.3362e-22,
        Amon:float = 2.3362e-22,
        rhomon:float = 3.27,
        ndnH:float = 3e-16,
        mH:float = 1.6736e-24,
        epsilonHe:float = 0.1

    ):
    """
    Creates R3D-files dust_density, dust_temperature and grain sizes per grid cell
    from pickles-files with approximated (ellipsoidal) dust distributions

    ARGUMENTS
    Paths:
      picklepath:str = path to c5d pickle-file
      amrpath:str = path to r3d-amr-grid-file
      gridpath:str = path to corresponding file with list of r3d radial grid-distances
      sizepath:str = path to corresponding file with list of r3d cell sizes
    
    Dust grain properties:
      monomermass: float = mass of dust species monomers
      Amon = 2.3362e-22 # g
      rhomon = 3.27 # g cm-3
      ndnH = 3e-16
      mP = 1.6726e-24 # g (mH = 1.6736e-24 # g)
      epsilonHe = 0.1

    RETURNS
      dust_density_approx.inp
      dust_temperature_approx.dat
      grainsizes_approx.dat
    """

    # Some output (tests the time it takes also)
    print(f'Start time:')
    os.system('date')
    print('Loading data')


    # Load pickle-file
    path = picklepath
    with open(path, 'rb') as f:
        ellipsoid_dict = pickle.load(f)

    # Number of dust clouds in data
    Nellipsoids = ellipsoid_dict['distance_centre_grid'].size
    print(f'  There are {Nellipsoids} dust clouds:')


    # Load an r3d-grid
    r3ddistances = a3d.load_griddistances(
        gridpath = gridpath,
        amrpath = amrpath,
    )
    r3dcellsizes = a3d.load_cellsizes(
        amrpath = amrpath,
        sizepath = sizepath,
    )
    nleafs = r3ddistances.shape[0]

    # Declare arrays
    # This is faster than looping through whole grid, writing zeros, again and again
    # condenfrac is condensationfraction used to get grain size
    r3d_densities = np.zeros(nleafs)
    r3d_temperatures = np.zeros(nleafs)
    r3d_condenfrac = np.zeros(nleafs)


    # Loop through each cloud
    for nellipsoid in range(Nellipsoids):

        # Number of c5d-grid-cells per cloud
        Nc5d = ellipsoid_dict['coord_d_list'][nellipsoid][:,0].size
        print(f'  Cloud {nellipsoid+1} has {Nc5d} c5d-cells.')

        # Coordinates of all c5d-grid-cells of each cloud
        cloud_coordsX = np.array(ellipsoid_dict['coord_d_list'][nellipsoid][:,0])
        cloud_coordsY = np.array(ellipsoid_dict['coord_d_list'][nellipsoid][:,1])
        cloud_coordsZ = np.array(ellipsoid_dict['coord_d_list'][nellipsoid][:,2])

        # Extract lists of 
        cloud_densities = np.array(ellipsoid_dict['filled_nmonomer'][nellipsoid][:])
        cloud_temperatures = np.array(ellipsoid_dict['filled_temp'][nellipsoid][:])
        cloud_concenfrac = np.array(ellipsoid_dict['filled_quc'][nellipsoid][:])


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
        # save data at correct places in r3d-arrays
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
                r3d_temperatures[nr3d] = cloud_temperatures[cloudcell].mean()
                r3d_condenfrac[nr3d] = cloud_concenfrac[cloudcell].mean()

    # Compute grain sizes
    print('  Computes grain sizes')
    nHnd = 1./ndnH
    grainsize_constants = 3/(4*np.pi) * Amon/rhomon * nHnd * mH * (1+epsilonHe)
    r3d_grainsizes = (grainsize_constants * r3d_condenfrac)**(1/3)

    # Write all data files
    print('  Prints r3d-files')


    with open(f'../dust_density_approx.inp', 'w') as fdensity, \
         open(f'../dust_temperature_approx.dat', 'w') as ftemperature, \
         open(f'../grainsizes_approx.dat', 'w') as fgrainsizes:

        # Write headers:

        # Density:
        # 1
        # nleafs
        # number dust species
        fdensity.write(f'1\n{int(nleafs)}\n1\n')

        # Temperature:
        # 1
        # nleafs
        # number dust species
        ftemperature.write(f'1\n{int(nleafs)}\n1\n')

        # Grain sizes
        # Just some info
        fgrainsizes.write('# List of grain sizes in cm\n# Same order as standard R3D-input files')

        for nr3d in range(nleafs):
            fdensity.write(f'{r3d_densities[nr3d]}\n')
            ftemperature.write(f'{r3d_temperatures[nr3d]}\n')
            fgrainsizes.write(f'{r3d_grainsizes[nr3d]}\n')

    print('Done at:')
    os.system('date')
