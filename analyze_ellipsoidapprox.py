#I have now uploaded two pickle files on the phy-exwings server, located in: 
# /home/ariah/dat/caucus (sidenote: caucus stands for ChAracterising dUst CloUdS hahaha):
#
#    filled-mean-ellipsoids_st28gm06n052-032.pickle – ellipsoids filled in with 
#    mean values (i.e. uniform distribution)
#    filled-radialdistribution-ellipsoids_st28gm06n052-032.pickle – ellipsoids 
#    filled in with an inverse radial distribution (e.g. gas density higher 
#    further away from the star)
#
# These are for the snapshot 032 of the model st28gm06n052. I restructured the data 
# format following your suggestion yesterday, about separating each ellipsoids into
# different arrays rather than having one big [679^3] domain – which saved a lot of
# space! Now the data is mostly in a list of arrays. You can load the pickle files
# as follows:
#
###
#
#import pickle
#
#with open('[INSERTDIRECTORY]/filled-mean-ellipsoids_lists.pickle', 'rb') as handle:
#    b = pickle.load(handle)
#
###
#
# where you would have loaded a dictionary with the keys:
#
#    'centre_ellipsoid' (2D array): indices of the centre of the ellipsoids [X,Y,Z]
#    'distance_centre_grid' (1D array): distance from of the centre of the ellipsoids 
#       to the star (in cm)
#    'coord_list' (list of arrays): each element in the list is an array of the 
#       relevant indices for the eliipsoids
#    'coord_d_list' (list of arrays): each element in the list is an array of the 
#       distances from the star [X,Y,Z] (in cm)
#    'filled_rho' (list of arrays): each element in the list is an array for the 
#       values of the gas density in each coordinate
#    'filled_temp' (list of arrays): each element in the list is an array for the
#       values of the temperature in each coordinate
#    'filled_quc' (list of arrays): each element in the list is an array for the
#       values of the condensation fraction in each coordinate
#    'filled_nmonomer’ (list of arrays): each element in the list is an array for 
#       the values of the number of monomers in each coordinate
#
# There are a total of 22 ellipsoids in this snapshot. If you would like the data
# above in the format of the whole [679^3] domain, I have a small function which
# can do this -- please refer to function attached. If things are confusing then
# please let me know, we can look through the data together sometime today or
# whenever you get to this! Thank you :- )
#
# Important
#    `filled_nmonomer’ should be the total number of monomers, per grain; per cell
#    `filled_quc’ is actually just the monomer number density (cm-3); per cell


# Import standard packages
import numpy as np
import pickle
import os

import analyze_r3d_functions as a3d
import create_r3d_functions as c3d


# Define some constants
AUcm = 1.49598e13 # AU in cm
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
    with only one dust-bin.

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
      grain_sizes_approx.dat
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
    r3d_monomergasratio = np.zeros(nleafs)


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
        cloud_densities = np.array(ellipsoid_dict['filled_quc'][nellipsoid][:])
        cloud_temperatures = np.array(ellipsoid_dict['filled_temp'][nellipsoid][:])
        cloud_gasdensity = np.array(ellipsoid_dict['filled_rho'][nellipsoid][:])


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
                r3d_monomergasratio[nr3d] = cloud_densities[cloudcell].mean() / cloud_gasdensity[cloudcell].mean()

    # Compute grain sizes (in cm)
    print('  Computes grain sizes')
    nHnd = 1./ndnH
    grainsize_constants = 3/(4*np.pi) * Amon/rhomon * nHnd * mH * (1+epsilonHe)
    r3d_grainsizes = (grainsize_constants * r3d_monomergasratio)**(1/3)

    # Write all data files
    print('  Prints r3d-files')


    with open(f'../dust_density_approx.inp', 'w') as fdensity, \
         open(f'../dust_temperature_approx.dat', 'w') as ftemperature, \
         open(f'../grain_sizes_approx.dat', 'w') as fgrainsizes:

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
        fgrainsizes.write('# List of grain sizes in cm\n# Same order as standard R3D-input files\n')

        for nr3d in range(nleafs):
            fdensity.write(f'{r3d_densities[nr3d]}\n')
            ftemperature.write(f'{r3d_temperatures[nr3d]}\n')
            fgrainsizes.write(f'{r3d_grainsizes[nr3d]}\n')

    print('Done at:')
    os.system('date')

# -----------------------------------------------------------------------------
# Func to write new density-temperature files where dust is binned in species
# based on grain size bins
#
def bin_inpdata(
        grainsizes_path:str = '../grain_sizes_binned.dat',
        density_path:str = '../dust_density_approx.inp',
        temperature_path:str = '../dust_temperature_approx.dat',
        wavelength_path:str = '../wavelength_micron.inp'
    ):
    """
    INFO
    ARGUMENTS
      grainsizes_path: path to dat-file containing the list of binned grainsizes per grid cell
                       as outputted by a5d.bin_grainsizes()
      density_path: path the dust_density.inp file with all dust in one species as
                    outputted by create_dustapproximation() above
      temperature_path: path the dust_temperature.dat file with all dust in one species as
                        outputted by create_dustapproximation() above
      wavelength_path: path to your wavelength_micron.inp file
    RETURNS    
    """
    print('Start binning density and temperature')

    # Load binned grain sizes per grid cell in cm
    grainsizes_grid = []
    with open(grainsizes_path, 'r') as fgrainsizes:
        for line in fgrainsizes.readlines():
            if line[0] != '#':
                grainsizes_grid.append(float(line))
    grainsizes_grid = np.array(grainsizes_grid)

    # Extract unique sizes
    grainsizes = np.unique(grainsizes_grid)
    # Only save the grainsizes that aren't zero
    grainsizes = grainsizes[np.where(grainsizes > 0)[0]]
    Nbins = len(grainsizes)
    # Save separate list with grainsizes in um
    grainsizes_list_um = (grainsizes*1e4).tolist()
    

    # Load density and temperature, extract nleafs (number of cells)
    nleafs, Nspecies, densities_1bin = a3d.load_dustdensity(
        path = density_path,
        numb_specie = 1
    )
    nleafs, Nspecies, temperatures_1bin = a3d.load_temperature(
        path = temperature_path,
        numb_specie = 1
    )


    # Write optool-script and dustopac file
    print('  Writing optool script and dustopac-files')
    c3d.create_optoolscript(
        wavelength_path = wavelength_path,
        phase = 'approx',
        grainum_sizes = grainsizes_list_um,
        grainsize_type = 'normal',
        grainsize_na = 21,
        specie = 'mg2sio4',
        grain_type = 'dhs',
        polarisation_matrix = 'n',
    )


    # Loop through grid and create and write new density and temperature-files
    # All empty cells are with 0s, rest are separated by species
    print('  Creating new density and temperature arrays.')

    densities_bins = np.zeros(Nbins*nleafs)
    temperatures_bins = np.zeros(Nbins*nleafs)

    # Loop through each grain size
    for nbin,agrain in enumerate(grainsizes):

        # Grin cell counter for error checking
        graincell_counter = 0

        # Loop through whole grid (for each grain size)
        for nn in range(nleafs):

            # Save densities at positions of each species bin (in cm)
            if grainsizes_grid[nn] == agrain:
                graincell_counter += 1
                nntotal = nbin*nleafs + nn
                densities_bins[nntotal] = densities_1bin[nn]
                temperatures_bins[nntotal] = temperatures_1bin[nn]

        # For each grain size bin, check if any data were saved
        if graincell_counter == 0:
            print(f'\n  WARNING: no data saved for grain size number {nbin+1}: {agrain*1e4}um\n')

    # Print new dust_density and dust_temperature-files
    print('  Writing new density and temperature inp-files')
    with open(f'../dust_density_approx_{Nbins}bins.inp', 'w') as fdensity, \
         open(f'../dust_temperature_approx_{Nbins}bins.dat', 'w') as ftemperature:

        # Write headers:
        #
        # Density:
        # 1
        # nleafs
        # number dust species
        fdensity.write(f'1\n{int(nleafs)}\n{int(Nbins)}\n')

        # Temperature:
        # 1
        # nleafs
        # number dust species
        ftemperature.write(f'1\n{int(nleafs)}\n{int(Nbins)}\n')

        for nn in range(Nbins*nleafs):
            fdensity.write(f'{densities_bins[nn]}\n')
            ftemperature.write(f'{temperatures_bins[nn]}\n')

    print(f'  ../dust_density_approx_{Nbins}bins.inp\n  ../dust_temperature_approx_{Nbins}bins.dat\nDONE')







