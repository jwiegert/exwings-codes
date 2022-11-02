# A storage of old functions that might come handy later but clutter my main
# files otherwise

# Import various libraries
import cython
import numpy as np
from scipy.io.idl import readsav
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm, xlabel, xscale, yscale
import os
import re

# My own libraries
import create_r3d_functions as c3d
import analyze_r3d_functions as a3d

# Basic definitions (AU as cython due to slow funcs below)
AUcm = cython.declare(cython.float ,1.49598e13) # cm
Msol = 1.989e33 # g
Rsol = 6.955e10 # cm
Lsol = 3.828e26 # W

# ------------------------------------------------------------ #
# List of functions
#
# This func separates star into various bins and gives each bin a
# constant rosseland opacity which is the average opacity of that bin
# I instead incorprorate Rosseland Opacity in Density file
# create_staropacity(
#    pathopacity:str='../star_opacities.dat',
#    pathstardensity:str='../dust_density_onestarstar.inp',
#    pathwavelength:str='../wavelength_micron.inp',
#    pathtemperature:str='../dust_temperature_onestar.dat',
#    temperaturebins:list=[],
#    N_opabins:int=2
# )
#
#
# load_dust_densitytemperature()
#    savpath:str='../co5bold_data/dst28gm06n052/st28gm06n052_186.sav',
#    nspecies:int=0
# )
#
#
# load_star_properties(
#    savpath:str='../co5bold_data/dst28gm06n056/st28gm06n056_140.sav'
# )


# ============================================================


# OLD FUNCTION
# Not used any more since I instead incorporate the Rosseland opacity in the density file
# for the star.
#
# Takes the stellar data from c5d-data that has been translated to r3d but as one specie, and
# transforms it into separate species on a chosen number of bins of opacities as given by the
# c5d-data.
def create_staropacity(
        pathopacity:str='../star_opacities.dat',
        pathstardensity:str='../dust_density_onestarstar.inp',
        pathwavelength:str='../wavelength_micron.inp',
        pathtemperature:str='../dust_temperature_onestar.dat',
        temperaturebins:list=[],
        N_opabins:int=2
    ):
    """
    ARGUMENTS
      pathopacity: path to star_opacities.dat',
      pathstardensity: path to dust_density_onestarstar.inp',
      pathwavelength: path to wavelength_micron.inp',
      pathtemperature: path to dust_temperature_onestar.dat',
      temperaturebins: list of temperature range bins (it adds lowest and higehst temperature automatically)
      N_opabins: number of kappa-bins spread logarithimically between min and max

    RETURNS
      dust_density_starbins.inp
      dust_temperature_starbins.dat
      star_opacities_bins.dat
      dustopac_starbins.inp
      dustkappa_star{no+1}.inp
    """

    Tbins = len(temperaturebins)

    if Tbins > 20:
        return 'WARNING too many files! Stopping'
    
    elif Tbins == 0:
        print('Only one bin, will only create r3d-opacity files, use *_onestar.inp files for density etc')
        
        # load opacity.dat
        opacity = []
        with open(pathopacity, 'r') as fopacity:
            for line in fopacity.readlines():
                if line[0] != '#':
                    opacity.append(float(line))

        # Load wavelengthgrid
        wavelengths,Nwave = a3d.load_wavelengthgrid(path=pathwavelength)

        # Take a SINGLE average of ALL opacity in c5d opacity-data
        opacityvalue = np.mean(np.array(opacity))

        # Save as constant as function of wavelength/frequency

        # Print new dustopac_starbins.inp file
        print('Writing opacity files for the star.')
        with open('../dustopac_starbins.inp', 'w') as fopac:

            # Print header
            fopac.write(f'2\n1\n-----------------------------\n')

            # Print star's opacity / species names
            fopac.write(f'1\n0\nstar\n-----------------------------\n')

        # Print opacity file, star-kappa-file
        with open(f'../dustkappa_star.inp', 'w') as fopac:

            # Write header (1 for no scattering in this and number of wavelengths)
            fopac.write(f'1\n{Nwave}\n')

            # Write wavelength, abscoeff, scattercoeff
            for wave in wavelengths:
                fopac.write(f'{wave}    {opacityvalue}    0.0\n')

        print('C5D create star opacities:\n    dustopac_starbins.inp\n    dustkappa_star.inp\nDONE\n')

    else:
        # Load denstiy_star.inp
        Ncells,Nspec,star_densities = a3d.load_dustdensity(path=pathstardensity,numb_specie=1)

        # Load temperature_star
        Ncells,Nspec,temperatures = a3d.load_temperature(path=pathtemperature)

        # load opacity.dat
        opacity = []
        with open(pathopacity, 'r') as fopacity:
            for line in fopacity.readlines():
                if line[0] != '#':
                    opacity.append(float(line))
        opacity = np.array(opacity)

        # Load wavelengthgrid
        wavelengths,Nwave = a3d.load_wavelengthgrid(path=pathwavelength)

        # Add lowest temperature to temperature list, change to array and update Tbins
        # The lowest temperature is not 0, this if statements makes sure it's the smallest
        # real temperature that is saved.
        if np.min(temperatures) == 0:
            temperaturebins.append(np.sort(temperatures)[1])
        else: 
            temperaturebins.append(np.min(temperatures))

        temperature_bins = np.array(temperaturebins)
        Tbins = np.size(temperature_bins)

        print('Temperature bin-ranges are (K):')
        for temp in temperature_bins:
            print(f'        {temp}')

        # Define opacity bins, no more than one per order of magnitude
        # Minimum opacity bin is the minimum value of opacity
        opamin = np.log10(min(opacity))

        # Maximum opacity is defined from above the highest tmeperature limit
        opamax = np.log10(np.mean(
            opacity[ np.where(temperatures > temperature_bins[0])[0]]
        ))
        
        # Array with OPA-bin-limits (numb of limits are one more)
        if N_opabins < 1:
            N_opabins = 1
        Obins = N_opabins + 1
        opabins = np.logspace(opamax, opamin, Obins)

        print(f'Opacity bin-ranges are (cm2/g):')
        for opa in opabins:
            print(f'        {opa}')

        # Bin the opacities by temperature ranges (the stellar surface is never above eg 5000k)
        # and by opacity range

        # Figure out number of bins (at least one bin always)
        nbin = 0
        nbins = 1
        listbins = np.zeros(Tbins*Obins)
        listbins[0] = 1

        for nT in range(1,Tbins):
            for nO in range(1,Obins):
                
                nbin += 1
                sumopacity = 0

                for nn in range(Ncells):
                    if \
                        temperature_bins[nT-1] > temperatures[nn] > temperature_bins[nT] and \
                        opabins[nO-1] > opacity[nn] > opabins[nO]:

                        sumopacity += opacity[nn]
                
                if sumopacity != 0:
                    nbins += 1

                    # Listbins lists the corresponding bin-number for each combination of nT and nO
                    # 0 where there is none
                    listbins[nbin] = nbins
                   

        # Start binning
        print(f'Binning to {nbins} species.')

        new_temperatures = np.zeros(nbins*Ncells)
        new_densities = np.zeros(nbins*Ncells)
        new_opacities = np.zeros(nbins*Ncells)

        for nn in range(Ncells):

            # First bin is inside the star, above highest T or OPA
            if temperatures[nn] > temperature_bins[0] or opacity[nn] > opabins[0]:
                new_temperatures[nn] = temperatures[nn]
                new_densities[nn] = star_densities[nn]
                new_opacities[nn] = opacity[nn]

            bincounter = 0

            # Remaining bins are outside and within these ranges
            for nT in range(1,Tbins):
                for nO in range(1,Obins):

                    bincounter += 1

                    if temperature_bins[nT-1] > temperatures[nn] > temperature_bins[nT] and \
                        opabins[nO-1] > opacity[nn] > opabins[nO]:

                        nbin = int(listbins[bincounter] - 1)

                        new_temperatures[nn + Ncells*nbin] = temperatures[nn]
                        new_densities[nn + Ncells*nbin] = star_densities[nn]
                        new_opacities[nn + Ncells*nbin] = opacity[nn]


        # Create new opacities from these limits
        opacityvalues = np.zeros(nbins)

        for nbin in range(1,nbins+1):

            # If there are no numbers then just let the opacity be 0
            if np.where(new_opacities[Ncells*(nbin-1):Ncells*nbin] > 0)[0].size > 0:

                opacityvalues[nbin-1] = np.mean(
                    new_opacities[
                        Ncells*(nbin-1) + np.where(new_opacities[Ncells*(nbin-1):Ncells*nbin] > 0)[0]
                    ]
                )

            #print(f'{nbin}: {Ncells*(nbin-1)}:{Ncells*nbin} OPAC: {opacityvalues[nbin-1]}')
            
        # Print new star-density, temperature and opacitybin-files
        print('Writing new radmc3d-files')

        with open('../dust_density_starbins.inp', 'w') as fdensity,\
             open('../dust_temperature_starbins.dat', 'w') as ftemperature,\
             open('../star_opacities_bins.dat', 'w') as fopacity:
            
            # Write headers
            fdensity.write(f'1\n{int(Ncells)}\n{int(nbins)}\n')
            ftemperature.write(f'1\n{int(Ncells)}\n{int(nbins)}\n')
            fopacity.write('# Binned opacities for the star(dust) distributions.\n# In same order as density file.\n# For reference and quality checks.\n')

            # Write densities and temperatures
            for nn,dens in enumerate(new_densities):
                fdensity.write(f'{dens}\n')
                ftemperature.write(f'{new_temperatures[nn]}\n')
                fopacity.write(f'{new_opacities[nn]}\n')

        # Print new dustopac_starbins.inp file
        print('Writing opacity files for the binned star.')
        with open('../dustopac_starbins.inp', 'w') as fopac:

            # Print header
            fopac.write(f'2\n{nbins}\n-----------------------------\n')

            # Print star's opacity / species names
            for nn in range(nbins):
                fopac.write(f'1\n0\nstar{nn+1:02d}\n-----------------------------\n')

        # Print opacity files, star-kappa-files
        for nbin,opac in enumerate(opacityvalues):
            with open(f'../dustkappa_star{nbin+1:02d}.inp', 'w') as fopac:

                # Write header (1 for no scattering in this and number of wavelengths)
                fopac.write(f'1\n{Nwave}\n')

                # Write wavelength, abscoeff, scattercoeff
                for wave in wavelengths:
                    fopac.write(f'{wave}    {opac}    0.0\n')
        
        # TODO Write info-file with data on how the binning is done! stellarbinning_info.txt
        print('C5D create star opacities:\n    dust_density_starbins.inp\n    dust_temperature_starbins.inp\n    star_opacities_bins.dat\n    dustopac_starbins.inp\n    dustkappa_starN.inp\nDONE\n')



# Extract co5bold densities into a separate array 
# - note this is probably faster than loading it in the r3d-file-writing functions
def load_star_properties(
        savpath:str='../co5bold_data/dst28gm06n056/st28gm06n056_140.sav'
    ):
    """
    Loads c5d-star's densities, temperatures, and opacities, and puts them into 3D-arrays

    INPUT
    savpath: path to sav-file
    nbins: number of bins to put the opacities in, ie number of duststar-species for the star

    OUTPUT
    c5dstar_densities: 3D array with star's densities
    c5dstar_temperatures: 3D array with star's temperatures 
    c5dstar_opacities: 3D array with star's opacities
    """
    # NOTE
    # Takes 4 seconds
    
    # Load sav-file
    c5ddata = readsav(savpath)

    # Extract data
    c5ddata = c5ddata['ful']

    # Save np.arrays
    c5dstar_densities = c5ddata['Z'][0][0][34]
    c5dstar_temperatures = c5ddata['EOS'][0][0][1]
    c5dstar_opacities = c5ddata['OPA'][0][0][0]
    
    return c5dstar_densities, c5dstar_temperatures, c5dstar_opacities


# Function for loading one dust specie from c5d-data
def load_dust_densitytemperature(
        savpath:str = '../co5bold_data/dst28gm06n052/st28gm06n052_186.sav',
        nspecies:int = 0
    ):
    """
    Loads c5d-data and extracts number density of dust monomers and dust temperature

    ARGUMENTS
      savpath:str = path to sav-file
      nspecies:int = number of the specie to extract, start with 0!

    RETURNS
      c5ddust_densities: array with monomer number density in the c5d-grid
      c5ddust_temperatures: array with dust temperatures within c5d-grid
    """
    # NOTE
    # Takes 4 seconds

    # Load sav-file
    c5ddata = readsav(savpath)

    # Extract data
    c5ddata = c5ddata['ful']

    # Declare np.arrays for number density of dust monomers and temperatures
    c5ddust_densities = c5ddata['Z'][0][0][40+3*nspecies]
    c5ddust_temperatures = c5ddata['EOS'][0][0][1]

    # TODO add something so that only dust-filled cell-temperatures are saved?

    # Return density-temperature arrays
    return c5ddust_densities, c5ddust_temperatures




# Function that loads and extracts gas densities and dust densities of a chosen dust specie
# Primarily to use for getting grain sizes
def load_dustgas_densities(
        savpath:str = '../co5bold_data/dst28gm06n052/st28gm06n052_186.sav',
        nspecies:int = 0
    ):
    """
    TODO
    fill here ...

    ARGUMENTS
      ...
    RETURNS
      ...    
    """
    # NOTE
    # Takes 4.3 seconds

    # Load sav-file
    c5ddata = readsav(savpath)

    # Extract data
    c5ddata = c5ddata['ful']

    # Declare np.arrays for number density of dust monomers and temperatures
    c5ddust_densities = c5ddata['Z'][0][0][40+3*nspecies]
    c5dstar_densities = c5ddata['Z'][0][0][34]

    # TODO only save gas-densities where dust exists here?
    # Or just do that directly in the grain-size-extractor

    return c5dstar_densities, c5ddust_densities
