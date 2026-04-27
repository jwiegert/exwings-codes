# Test extract average grain size from r3d-input data
# Volume-average grain size (cell-size weighted)
import numpy as np
import os
#
# <agrain> = sum(agran(ncell) * volume(ncell))/Volume(box)
#
# Total boxvolume is
AUcm = 1.49598e13 # cm
compbox_volume = (29.780836060395945 * AUcm)**3
#
# List of modelnames
models = [
    'st28gm06n052',
    'st28gm06n074',
    'st28gm06n075',
]
for model in models:
    modelpath = f'../r3dresults/{model}/'
    # Extract all snapshots
    phases = [int(filename) for filename in os.listdir(modelpath) if os.path.isdir(modelpath+filename)]
    phases.sort()
    Nphases = len(phases)
    # For each model extract cellsizes
    cellsizes = np.loadtxt(f'{modelpath}grid_cellsizes.csv')**3
    # Declare variable to total average grainsize
    total_average_grainsize = 0
    # Write file with all averages
    with open(f'{modelpath}average_grainsizes.dat', 'w') as fgrainsize:
        # Write header
        fgrainsize.writelines(
            f'# Time evolution of volume-averaged-weighted grain sizes\n# In the whole computation cube\n# For model {model}\n'
        )
        # Loop over snapshots
        for snapshot in phases:
            # Load all grainszies
            grainsizes_snapshot = np.loadtxt(
                f'{modelpath}{snapshot:03d}/grain_sizes_binned_{snapshot:03d}.dat'
            )
            # Find volume-averaged grain size for each snapshot
            # transform to um
            average_grainsize = np.sum(grainsizes_snapshot * cellsizes) / compbox_volume * 1e4
            #
            # Write to file
            fgrainsize.writelines(
                f'{snapshot:03d}    {average_grainsize:.5f}\n'
            )
            total_average_grainsize += average_grainsize
    # Print total average grain size for the model
    print(f'    {model} : {total_average_grainsize/Nphases} um')

