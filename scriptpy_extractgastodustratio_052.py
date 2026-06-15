# Simple script to extract gas-to-dust ratios of all snapshots
# of each model and writing these numbers to dat-files for
# each model.
import analyze_r3d_functions as a3d
import numpy as np
import os

print('Started extracting gas-to-dust-ratios.')

# Set an inner radius to measure "all" gas masss
AUcm = 1.49598e13 # AU in cm
gasradius = 2*AUcm

# Set paths
path = '../r3dresults/'
#path = '../../exwings_archivedata/'
models = [
    'st28gm06n052',
    #'st28gm06n074',
    #'st28gm06n075',
]
for model in models:
    print(
        f'  Doing {model}'
    )
    # Extract all snapshots included for this model.
    modelpath = path+model+'/'
    #phases = [
    #    f'{int(filename):03d}' for filename in os.listdir(modelpath) if os.path.isdir(modelpath+filename)
    #]
    #phases.sort()
    # Instead of doing all, take subsamples
    phases = [
        f'{nphase:03d}' for nphase in range(1,100+1)
    ]
    # Load grid properties
    #
    # Load grid_distances (unit is cm)
    griddistances = a3d.load_griddistances(
        gridpath=f'{modelpath}grid_distances.csv',
        amrpath=f'{modelpath}amr_grid.inp'
    )
    gridradii = griddistances[:,0]
    #
    # Load grid-cellsizes
    gridsizes = a3d.load_cellsizes(
        sizepath=f'{modelpath}grid_cellsizes.csv',
        amrpath=f'{modelpath}amr_grid.inp',
    )
    # Open output file
    with open(f'../{model}_gastodust_ratio.dat', 'w') as fratios:
        # Print header
        fratios.writelines(f'# Gas to dust ratios for snapshots of model {model}.\n')
        fratios.writelines('# Also lists dust smallest dust formation radius i AU.\n')
        fratios.writelines('#\n')
        fratios.writelines('# Snapshot   Gas-in-dust-cells    Gas-outside-dust-form-R    Dust-form-R-(AU)\n')
        # Loop through all snapshots
        for phase in phases:
            print(f'  - Doing snapshot {phase}')
            #
            # Load gas densities*kappaross
            Ncells,Nspecies,gas_densityopacity = a3d.load_dustdensity(
                path=f'{modelpath}{phase}/dust_density_opastar.inp',
                numb_specie=1
            )
            # Load gas kappaross
            gas_opacity = np.loadtxt(
                f'{modelpath}{phase}/star_opacities_smoothed.dat'
            )
            # Load dust-densities at once
            dust_densities = np.zeros(Ncells)
            with open(f'{modelpath}{phase}/dust_density_dust.inp' , 'r') as fdustdens:
                # Set dust species counter
                numb_specie = 0
                for nn,line in enumerate(fdustdens.readlines()):
                    # Total number of dust species
                    if nn == 2:
                        Nspecies = int(line)
                    # Extract dust densities
                    if nn > 2:
                        # Update dust species counter when needed
                        if nn-3-Ncells*numb_specie >= Ncells:
                            numb_specie += 1
                        if numb_specie >= Nspecies:
                            raise ValueError('ERROR: too many dust species')
                        # Then fill dust density array
                        dust_densities[nn-3-Ncells*numb_specie] += float(line)
            #
            # Declare temporary variables
            dustmass = 0                 # total dust mass in grams
            gasmass_dustcells = 0        # Gas mass within dust cells only
            gasmass_dustradius = 0       # gas mass outside dust-form-radius
            dustform_radius = 30*AUcm    # Placeholder for dust radius outside box
            #
            # Extract dust and gas masses
            # Loop through cells
            for ncell in range(Ncells):
                #
                # Only look in cells with dust first
                if dust_densities[ncell] > 0:
                    #
                    # Grid cell volume in cm3
                    cellvolume = gridsizes[ncell]**3
                    #
                    # Add up dust masses
                    dustmass += dust_densities[ncell]*cellvolume
                    #
                    # Check for smallest altitude to centrum of grid
                    if dustform_radius > gridradii[ncell]:
                        dustform_radius = gridradii[ncell]
                    #
                    # Save gas masses within dust filled cells.
                    gasmass_dustcells += gas_densityopacity[ncell] / gas_opacity[ncell] * cellvolume
            #
            # Also add up gas masses of all cells outside dust forming radius
            for ncell in range(Ncells):
                if gridradii[ncell] >= dustform_radius:
                    # Grid cell volume in cm3
                    cellvolume = gridsizes[ncell]**3
                    # Dustmass of all cells outside Rdust
                    gasmass_dustradius += gas_densityopacity[ncell] / gas_opacity[ncell] * cellvolume
            #
            # Get gas-to-dust ratios and convert formation radius unit
            gastodust_dustcells = gasmass_dustcells / dustmass
            gastodust_dustradius = gasmass_dustradius / dustmass
            dustform_radius /= AUcm
            #
            # Write to file.
            fratios.writelines(
                f'  {phase}        {gastodust_dustcells:.4f}            {gastodust_dustradius:.4f}                  {dustform_radius:.4f}\n'
            )
print('DONE')
# Rename and move final files to each nospikes-model-folder under r3dresults
# in the terminal or shell-script.
#os.system('spd-say moo')
