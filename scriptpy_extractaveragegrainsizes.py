# extract average grain sizes
import analyze_co5bold_functions as a5d
import numpy as np
import pexpect
import os

# Total boxvolume is
AUcm = 1.49598e13 # cm
compbox_volume = (29.780836060395945 * AUcm)**3

# List of modelnames
models = [
    'st28gm06n052',
#    'st28gm06n074',
#    'st28gm06n075',
]
# Loop over models
for model in models:

    # Set paths
    workpath = f'../{model}/'
    modelpath = f'../../exwings/r3dresults/{model}/'

    # Extract all snapshots
    phases = [int(filename) for filename in os.listdir(modelpath) if os.path.isdir(modelpath+filename)]
    phases.sort()
    Nphases = len(phases)

    # For each model extract cellsizes
    cellsizes = np.loadtxt(
        f'{modelpath}grid_cellsizes.csv'
    )**3

    # Declare variable to total average grainsize
    total_average_grainsize = 0

    # Write file with all averages
    with open(f'{workpath}average_grainsizes.dat', 'w') as fgrainsize:
        # Write header
        fgrainsize.writelines(
            f'# Time evolution of volume-averaged-weighted grain sizes\n# In the whole computation cube\n# For model {model}\n'
        )
        # Loop over snapshots
        for snapshot in phases:
            print(f'===== Doing {model}_{snapshot:03d} =====')
            # copy save file
            cpsavcommand = pexpect.spawn(f'sudo cp -v /home/befre677/dat/rhd/d{model}/idl/{model}_{snapshot:03d}.sav {workpath}')
            cpsavcommand.expect('password')
            cpsavcommand.sendline('uppexwPlump123')
            cpsavcommand.interact()
            # Extract grainsizes from savfile
            a5d.load_c5dheavydata(
                savpath = f'{workpath}{model}_{snapshot:03d}.sav',
                Nspecies = 1,
                gas_density = True,
                dust_density = True,
                gas_opacity = False,
                temperature = False
            )
            a5d.extract_grainsizes(
                amrpath=f'{modelpath}amr_grid.inp',
                gridpath=f'{modelpath}grid_distances.csv',
                sizepath=f'{modelpath}grid_cellsizes.csv',
                savpath=f'{workpath}{model}_{snapshot:03d}.sav',
                Amon=2.3362e-22,
                rhomon=3.27,
                ndnH=3e-16,
                mH=1.6736e-24,
                epsilonHe=0.1
            )
            # Output is ../grain_sizes_SNAPSHOT.dat
            # move to work folder
            os.system(f'mv ../grain_sizes_{snapshot:03d}.dat {workpath}')
            #
            # Load all grainszies
            grainsizes_snapshot = np.loadtxt(
                f'{workpath}grain_sizes_{snapshot:03d}.dat'
            )
            # Find volume-averaged grain size for each snapshot
            # transform to um
            average_grainsize = np.sum(grainsizes_snapshot * cellsizes) / compbox_volume * 1e4
            #
            # Write to file
            fgrainsize.writelines(
                f'{snapshot:03d}    {average_grainsize:.5f}\n'
            )
            # Sum all snapshos's average grain sizes
            total_average_grainsize += average_grainsize
            # Remove npy-files
            os.system(f'rm -v ../*{snapshot:03d}.npy')
            # Remove sav-file
            rmsavcommand = pexpect.spawn(f'rm -v {workpath}{model}_{snapshot:03d}.sav')
            rmsavcommand.expect('remove write-protected regular file')
            rmsavcommand.sendline('y')
            rmsavcommand.interact()
        # Print total average grain size for the model
        fgrainsize.writelines(
            f'# -----\n# Total average value for {model} : {total_average_grainsize/Nphases} um\n#'
        )
