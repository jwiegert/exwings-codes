# Simple script to extract gas-to-dust ratios of all snapshots
# of each model and writing these numbers to dat-files for
# each model.
import analyze_r3d_functions as a3d
import numpy as np
import scipy
import os

print('Started extracting gas-to-dust-ratios.')

path = '../r3dresults/'
models = [
    'st28gm06n052',
    'st28gm06n074',
    'st28gm06n075',
]
for model in models:
    print(f'  Doing {model}')
    # Extract all snapshots included for this model.
    modelpath = path+model+'/'
    phases = [
        int(filename) for filename in os.listdir(modelpath) if os.path.isdir(modelpath+filename)
    ]
    phases.sort()
    Nphases = len(
        phases
    )
    # Open output file

    with open(f'../{model}_gastodust_ratio.dat', 'w') as fratios:
        # Print header
        fratios.writelines(f'# Gas to dust ratios for snapshots of moderl {model}.\n')
        fratios.writelines(f'# Based on {Nbins} bins of dust-filled cells.\n')
        fratios.writelines(f'#\n')
        fratios.writelines(f'# Snapshot  Mean     STD      Median   MAD\n')
        # Loop through all snapshots
        for phase in phases:
            # Load gas density-opacity-data
            Ncells,Nspecies,gas_densityopacity = a3d.load_dustdensity(
                path=f'{modelpath}{phase:03d}/dust_density_opastar.inp',
                numb_specie=1
            )
            # Load gas kappaross-data
            gas_opacity = np.loadtxt(
                f'{modelpath}{phase:03d}/star_opacities_smoothed.dat'
            )
            # Load dust densities and combine to one array
            # First load number 1 and number of species
            Ncells,Nspecies,dust_densities = a3d.load_dustdensity(
                path=f'{modelpath}{phase:03d}/dust_density_dust.inp',
                numb_specie=1
            )
            # Then loop over rest of speces, add all to same array
            for nspecies in range(2,Nspecies+1):
                Ncells,Nspecies,dust_density = a3d.load_dustdensity(
                    path=f'{modelpath}{phase:03d}/dust_density_dust.inp',
                    numb_specie=nspecies
                )
                dust_densities += dust_density
            #
            # Create a dust-to-gas list
            #
            dusttogas = []
            for ncell in range(Ncells):
                if dust_densities[ncell] > 0:
                    dusttogas.append(
                        dust_densities[ncell] / (gas_densityopacity[ncell]/gas_opacity[ncell])
                    )
            # Bin dust-to-gas cell-ratios into 100 bins
            Nbins = 100
            Ndustcells = len(
                dusttogas
            )
            binsize = int(
                np.round(Ndustcells/Nbins)
            )
            dusttogas_bins = np.zeros(
                Nbins
            )
            for nbin in range(Nbins):
                # Initial bin-index
                ind_init = int(
                    nbin*binsize
                )
                # Final bin-index (cut at end of list)
                if (nbin+1)*binsize < Ndustcells:
                    ind_final = int(
                        (nbin+1)*binsize
                    )
                else:
                    ind_final = int(
                        Ndustcells
                    )
                # Take the max dust-to-gas-ratio of each bin
                # to avoid noise with signularities
                dusttogas_bins[nbin] = np.max(
                    dusttogas[
                        ind_init:ind_final
                    ]
                )
            gastodust_bins = 1/dusttogas_bins
            # Extract mean, std, median, mad
            ratio_mean = np.mean(gastodust_bins)
            ratio_std = np.std(gastodust_bins)
            ratio_median = np.median(gastodust_bins)
            ratio_mad = scipy.stats.median_abs_deviation(gastodust_bins)
            # Write to file.
            fratios.writelines(
                f'  {phase:03d}       {ratio_mean:.1f}   {ratio_std:.1f}   {ratio_median:.1f}   {ratio_mad:.1f}'
            )
print('DONE')