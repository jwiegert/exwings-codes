import os
import pickle
import numpy as np
import data_unpacker
import analyze_r3d_functions as a3d
import create_r3d_functions as c3d
import analyze_co5bold_functions as a5d
#
# Script to load and translate A's picklefiles to R3D-format
print('Running: scriptpy to translate pickle-file data to r3d-format')

# Define model parameters and file names
modelname = 'st28gm06n052'
approxdesignations = [
    'derivedellipsoids'
]
snapshots = [
    '200_filled-outerellipsoid','300_filled-outerellipsoid'
]
snapshot_numbers = [
    200,300
]

# Loop over all approximations
for approxdesignation in approxdesignations:

    # Create model folder
    os.system(f'mkdir ../arief_data/{modelname}_{approxdesignation}')

    # Loop over all included snapshot numbers
    for snapshot,snapshot_number in zip(snapshots,snapshot_numbers):

        print(f'  Doing {modelname}_{approxdesignation}: {snapshot}')

        # Define and create output path
        outputpath = f'../arief_data/{modelname}_{approxdesignation}/{snapshot}'
        os.system(f'mkdir {outputpath}')

        # Load relevant pickle file
        # i.e. either as obtained from CO5BOLD (clumps); or derived ellipsoids
        print(f'  Loading pickle file: {modelname}-{approxdesignation}_{snapshot}.pickle')
        filename = f'../arief_data/{modelname}-{approxdesignation}_{snapshot}.pickle'

        with open(filename, 'rb') as handle:
            dustclumps = pickle.load(handle)

        # Use Ariefs data unpacker to create new lists of data
        # Skip rho_gas, use original star for comparisons
        grid_filled_rho_dust=data_unpacker.data_unpacker(
            679,
            dustclumps['coord_list'],
            dustclumps['rho_dust_list']
        )
        grid_filled_temperature=data_unpacker.data_unpacker(
            679,
            dustclumps['coord_list'],
            dustclumps['temperature_list']
        )
        grid_filled_grainsizes=data_unpacker.data_unpacker(
            679,
            dustclumps['coord_list'],
            dustclumps['grainsizes_list']
        )
        # Change to 3D numpy-array and rotate back to same angle-combo
        # that the original data give when translated directly from
        # co5bold approx-files

        c5ddust_densities = np.flip(np.rot90(np.rot90(np.rot90(
            np.array(grid_filled_rho_dust),
            k=-1, axes=(0,1)),
            k=-1, axes=(1,2)),
            k=-1, axes=(0,1)),
            axis=1
        )
        c5dstar_temperatures = np.flip(np.rot90(np.rot90(np.rot90(
            np.array(grid_filled_temperature),
            k=-1, axes=(0,1)),
            k=-1, axes=(1,2)),
            k=-1, axes=(0,1)),
            axis=1
        )
        c5ddust_grainsizes = np.flip(np.rot90(np.rot90(np.rot90(
            np.array(grid_filled_grainsizes),
            k=-1, axes=(0,1)),
            k=-1, axes=(1,2)),
            k=-1, axes=(0,1)),
            axis=1
        )

        # Load C5D-grid
        print('  Loading CO5BOLD grid')
        c5dgrid = np.loadtxt('../arief_data/gridc5d_distances.csv')
        Nc5dside = len(c5dgrid)

        # Load R3D-grid
        print('  Loading RADMC3D grid')
        nxyz,nrefines,nleafs,nbranch,gridedge = a3d.load_grid_properties(
            amrpath='../arief_data/amr_grid.inp'
        )
        r3ddistances = a3d.load_griddistances(
            gridpath='../arief_data/grid_distances.csv',
            amrpath='../arief_data/amr_grid.inp'
        )
        r3dcellsizes = a3d.load_cellsizes(
            amrpath='../arief_data/amr_grid.inp',
            sizepath='../arief_data/grid_cellsizes.csv'
        )

        # Start translating to R3D-grid to various inp-files
        # Start mapping data
        print('  Translating data to RADMC3D grid')
        r3d_dustdensities = 0
        r3d_temperatures = 0
        r3d_grainsizes = 0
        progbar = 0

        with open(f'../dust_density_dust.inp', 'w') as fdustdensity, \
            open(f'../dust_temperature_onestar.dat','w') as ftemperature, \
            open(f'../grain_sizes.dat','w') as fgrainsize:

            # Write headers:
            # 1
            # nleafs
            # number dust species (for density and temperature files)
            fdustdensity.write(f'1\n{int(nleafs)}\n1\n')
            ftemperature.write(f'1\n{int(nleafs)}\n1\n')
            # Info for grain size file
            fgrainsize.write('# List of grain sizes for each cell. Not-binned.\n')
            fgrainsize.write('# Same order as in R3D-density and temperature-files\n')

            # Loop over R3D-grid
            for nr3d in range(nleafs):

                # Extract size range for current r3dcell
                r3dxrange = [
                    r3ddistances[nr3d,1]-0.5*r3dcellsizes[nr3d],
                    r3ddistances[nr3d,1]+0.5*r3dcellsizes[nr3d]
                ]
                r3dyrange = [
                    r3ddistances[nr3d,2]-0.5*r3dcellsizes[nr3d],
                    r3ddistances[nr3d,2]+0.5*r3dcellsizes[nr3d]
                ]
                r3dzrange = [
                    r3ddistances[nr3d,3]-0.5*r3dcellsizes[nr3d],
                    r3ddistances[nr3d,3]+0.5*r3dcellsizes[nr3d]
                ]
                # Extract indeces of all c5dcells within current r3dcell
                c5dxrange = np.argwhere(
                    r3dxrange[0] <= c5dgrid[np.argwhere(c5dgrid <= r3dxrange[1])]
                )[:,0]
                c5dyrange = np.argwhere(
                    r3dyrange[0] <= c5dgrid[np.argwhere(c5dgrid <= r3dyrange[1])]
                )[:,0]
                c5dzrange = np.argwhere(
                    r3dzrange[0] <= c5dgrid[np.argwhere(c5dgrid <= r3dzrange[1])]
                )[:,0]
                # Counter for number of c5dcells within r3dcell
                nchildcells = 0

                # Then loop through c5dcells within r3dcell
                for nnz in c5dzrange:
                    for nny in c5dyrange:
                        for nnx in c5dxrange:

                            # Save sum of all c5d-data for each r3d-cell
                            r3d_dustdensities += c5ddust_densities[nnx,nny,nnz]
                            r3d_temperatures += c5dstar_temperatures[nnx,nny,nnz]
                            r3d_grainsizes += c5ddust_grainsizes[nnx,nny,nnz]

                            # Number of cells
                            nchildcells += 1

                # Check if there actually are any c5dcells within r3dcell
                # If not, then your r3dgrid is probably smaller than the c5dgrid
                # and then the density and temperature will be zero for some cells
                if nchildcells > 0:
                    # Otherwise save the average of the c5d-data of each r3dcell
                    r3d_dustdensities /= nchildcells
                    r3d_temperatures /= nchildcells
                    r3d_grainsizes /= nchildcells
                else:
                    r3d_dustdensities = c5ddust_densities[c5dxrange[0],c5dyrange[0],c5dzrange[0]]
                    r3d_temperatures = c5dstar_temperatures[c5dxrange[0],c5dyrange[0],c5dzrange[0]]
                    r3d_grainsizes = c5ddust_grainsizes[c5dxrange[0],c5dyrange[0],c5dzrange[0]]
                    print(f'  ERROR: nchildcells = {nchildcells}')

                # Then write data to r3d files
                fdustdensity.write(f'{r3d_dustdensities}\n')
                ftemperature.write(f'{r3d_temperatures}\n')
                fgrainsize.write(f'{r3d_grainsizes}\n')

                # Reset data
                r3d_dustdensities = 0
                r3d_temperatures = 0
                r3d_grainsizes = 0

                # Some progress bar info
                if int(nr3d/nleafs*100) == 25 and progbar == 0:
                    progbar += 1
                    print('    Finished 25 per cent of the grid.')

                if int(nr3d/nleafs*100) == 50 and progbar == 1:
                    progbar += 1
                    print('    Finished 50 per cent of the grid.')

                if int(nr3d/nleafs*100) == 75 and progbar == 2:
                    progbar += 1
                    print('    Finished 75 per cent of the grid.')

        print('    Finished writing inp/dat-files, moving to correct folders')
        os.system(f'mv ../dust_density_dust.inp {outputpath}/')
        os.system(f'mv ../dust_temperature_onestar.dat {outputpath}/')
        os.system(f'mv ../grain_sizes.dat {outputpath}/')


        print('  Running: binning of grain sizes to 10 bins and ')
        # Bin sizes
        a5d.bin_grainsizes(
            grainsizepath=f'{outputpath}/grain_sizes.dat',
            phase=snapshot,
            nbins=10,
            lin='y'
        )
        # Output is ../grain_sizes_binned_{phase}.dat
        # Move to correct folder
        os.system(f'mv ../grain_sizes_binned_{snapshot}.dat {outputpath}/grain_sizes_binned.dat')


        # Load r3d-dust density file, not-binned
        Ncells,Nspec,dust_density = a3d.load_dustdensity(
            path=f'{outputpath}/dust_density_dust.inp'
        )
        # Load r3d-temperature-file, not binned
        # NOTE! I don't care about re-normalising the temperature yet to bladh-approx!
        Ncells,Nspec,dust_temperature = a3d.load_temperature(
            path=f'{outputpath}/dust_temperature_onestar.dat'
        )

        # Load binned grain sizes (first check if they exist)
        print('    Loading grain sizes')
        if os.path.exists(f'{outputpath}/grain_sizes_binned.dat') == True:
            grainsizes,Nleafs = a3d.load_grainsizes(
                grainsize_path=f'{outputpath}/grain_sizes_binned.dat'
            )
            # Change unit to um
            grainsizes *= 1e4

            # Extract an array with the grain sizes only
            grainsizes_uniq = np.unique(grainsizes[np.where(grainsizes > 0)[0]])

            # List of grain sizes and number of bins (in micrometers!)
            Ngrainsizes = grainsizes_uniq.size
            print(f'    Available ({Ngrainsizes}) grain sizez are (um): {grainsizes_uniq}')
        else:
            # Create place-holder lists
            grainsizes = np.zeros(nleafs)
            grainsizes_uniq = [0]
            Ngrainsizes = 1


        # Open new r3d-files to bin data in
        with open(f'{outputpath}/dust_density_dust_binned.inp', 'w') as fdensity, \
            open(f'{outputpath}/dust_temperature_dust_binned.dat', 'w') as ftemperature:

            # Write headers:
            #
            # Density:
            # 1
            # nleafs
            # number dust species
            fdensity.write(f'1\n{int(nleafs)}\n{int(Ngrainsizes)}\n')

            # Temperature:
            # 1
            # nleafs
            # number dust species
            ftemperature.write(f'1\n{int(nleafs)}\n{int(Ngrainsizes)}\n')

            # Loop of grain sizes of each specie
            for nsize,size in enumerate(grainsizes_uniq):

                # Write densities and temperatures in files according to R3D-syntax.
                # Each grain size bin is listed in same files but after eachother.

                # Some output
                print(f'    Writing dust grain size numb {nsize}: {size:.3e}um')

                # Loop over the r3d-grid
                for nr3d in range(nleafs):

                    # kolla om nr3d-positionens grain size stämmer med size ovanför
                    # då kan jag skriva ned density och temperatur i rätt fil
                    # annars skriver den noll
                    # så delar jag upp  den i grain size bins
                    # TODO
                    # In case it's needed, adapt dust temperature here to
                    # Bladh-approximation

                    if grainsizes[nr3d] == size:
                        fdensity.write(f'{dust_density[nr3d]}\n')
                        ftemperature.write(f'{dust_temperature[nr3d]}\n')

                    else:
                        fdensity.write('0.0\n')
                        ftemperature.write('0.0\n')


        # Write opac-file and optool-script
        print('  Writing opac-files.')

        specie = 'mg2sio4'

        c3d.create_optoolscript(
            wavelength_path=f'../arief_data/wavelength_micron.inp',
            phase=snapshot,
            grainum_sizes=f'{outputpath}/grain_sizes_binned.dat',
            grainsize_type='normal',
            grainsize_na=21,
            specie=specie,
            grain_type='dhs'
        )
        # Move files
        os.system(f'mv ../optool_script_{snapshot}.sh {outputpath}/optool_script.sh')
        os.system(f'mv ../dustopac_{specie}_{snapshot}.inp {outputpath}/dustopac_dust.inp')
        # Run script
        os.system(f'{outputpath}/optool_script.sh')
        # Move results
        os.system(f'mv *mg2sio4* {outputpath}/')

        # Copy over star-only-files and merge with dust data
        print('  Copying star-only-files.')
        os.system(f'cp -v ../arief_data/st28gm06n052_nodust/{snapshot_number}/* {outputpath}')
        # Merge stellar and dust data
        data_unpacker.merge_final_data(
            workpath=outputpath
        )
        # Clean up
        os.system(f'rm {outputpath}/*opastar*')
        os.system(f'rm {outputpath}/*onestar*')

        print(f'  Done: {modelname}_{approxdesignation}: {snapshot}\n')
