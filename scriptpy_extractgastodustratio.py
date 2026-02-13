# Simple script to extract gas-to-dust ratios of all snapshots
# of each model and writing these numbers to dat-files for
# each model.
import analyze_r3d_functions as a3d
import numpy as np
import os

print('Started extracting gas-to-dust-ratios.')
AUcm = 1.49598e13 # AU in cm

#path = '../r3dresults/'
path = '../../exwings_archivedata/'
models = [
    'st28gm06n052',
    #'st28gm06n074',
    'st28gm06n075',
]
for model in models:
    print(
        f'  Doing {model}'
    )
    # Extract all snapshots included for this model.
    modelpath = path+model+'/'
    phases = [
        int(filename) for filename in os.listdir(modelpath) if os.path.isdir(modelpath+filename)
    ]
    phases.sort()
    Nphases = len(
        phases
    )
    # Load grid properties
    #
    # Load grid_distances (unit is cm)
    griddistances = a3d.load_griddistances(
        gridpath=f'{path}{model}/grid_distances.csv',
        amrpath=f'{path}{model}/amr_grid.inp'
    )
    gridradii = griddistances[:,0]
    #
    # Load grid-cellsizes
    gridsizes = a3d.load_cellsizes(
        sizepath=f'{path}{model}/grid_cellsizes.csv',
        amrpath=f'{path}{model}/amr_grid.inp',
    )
    # Open output file
    with open(f'../{model}_gastodust_ratio.dat', 'w') as fratios:
        # Print header
        fratios.writelines(f'# Gas to dust ratios for snapshots of model {model}.\n')
        fratios.writelines(f'# Also lists dust smallest dust formation radius i AU.\n')
        fratios.writelines(f'# All gas mass is all gas mass outside the dust form radius.\n')
        fratios.writelines(f'#\n')
        fratios.writelines(f'# Snapshot    All-gas    Gas-in-dust    Dust-formation-radius(AU)\n')
        # Loop through all snapshots
        for phase in phases:
            print(f'  Doing snapshot {phase:03d}')
            #
            # Load gas densities*kappaross
            Ncells,Nspecies,gas_densityopacity = a3d.load_dustdensity(
                path=f'{path}/{model}/{phase:03d}/dust_density_opastar.inp',
                numb_specie=1
            )
            # Load gas kappaross
            gas_opacity = np.loadtxt(
                f'{path}/{model}/{phase:03d}/star_opacities_smoothed.dat'
            )
            # Load dust-densities
            # First load number 1 and number of species
            Ncells,Nspecies,dust_densities = a3d.load_dustdensity(
                path=f'{path}/{model}/{phase:03d}/dust_density_dust.inp',
                numb_specie=1
            )
            # Then loop over rest of species, add all to same array
            for nspecies in range(2,Nspecies+1):
                Ncells,Nspecies,dust_density = a3d.load_dustdensity(
                    path=f'{path}/{model}/{phase:03d}/dust_density_dust.inp',
                    numb_specie=nspecies
                )
                dust_densities += dust_density
            # 
            # Declare temporary variables
            dustmass = 0       # total dust mass in grams
            gasmass_indust = 0 # gas mass in dustcontaining cells
            gasmass_all = 0    # "total" gas mass 
            dustform_radius = 30*AUcm
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
                    # Add up gas masses of dust-filled cells
                    gasmass_indust += gas_densityopacity[ncell]*gas_opacity[ncell]*cellvolume
                    #
                    # Check for smallest altitude to centrum of grid
                    if dustform_radius > gridradii[ncell]:
                        dustform_radius = gridradii[ncell]
            #
            # Then add up all gas masses outside dust formation radius
            for ncell in range(Ncells):
                if gridradii[ncell] >= dustform_radius:
                    gasmass_all += gas_densityopacity[ncell]*gas_opacity[ncell]*cellvolume
            #
            # Get gas-to-dust ratios and convert formation radius unit
            gastodust_indust = gasmass_indust/dustmass
            gastodust_allgas = gasmass_all/dustmass
            dustform_radius /= AUcm
            #
            # Write to file.
            fratios.writelines(
                f'  {phase:03d}         {gastodust_allgas:.3f}   {gastodust_indust:.5f}       {dustform_radius:.3f}'
            )
print('DONE')
# TODO
# rename and move final files to each nospikes-model-folder under r3dresults
os.system('spd-say moo')
