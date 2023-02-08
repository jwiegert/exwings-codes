# functions and tools for loading and analyzing DARWIN data, and to create 
# RADMC3D-input models from DARWIN data.

# Import various libraries
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm, xlabel, xscale, yscale

# My own libraries
import analyze_r3d_functions as a3d
import create_r3d_functions as c3d
#import analyze_co5bold_functions as a5d

# Basic definitions
AUcm = 1.49598e13 # cm
Msol = 1.989e33 # g
Rsol = 6.955e10 # cm
Lsol = 3.828e26 # W

# Note
# Co5bold-Rstar = 1.651AU (355 Rsun)

# ------------------------------------------------------------ #
# List of functions
# 
#
# load_darwindata(
#    modelname:str='M2n315u6'
# )
#
#




# ============================================================


# Load DARWIN-data and save in arrays
# -----------------------------------

def load_darwindata(
        modelname:str='M2n315u6',
        timestep:str='230791'
    ):
    """
    Loads and extracts densities and temperatures of one model and timestep

    ARGUMENTS
      modelname: str: Letter-number-code for the model.
      timestep: str: timestep number (ie model.XXXXXX.dat)

    RETURNS
      average_radius: np.array with radial range in cm
      average_density: np.array with time-averaged radial dependant gas densities in g/cm3
      average_temperature: np.array with time-averaged radial dependant gas temperatures in K
    """

    filename = f'../darwin_data/{modelname}_sel/model.{timestep}.dat'


    radius_cm = []
    gas_density = []
    gas_temperature = []
    # TODO gas_opacity

    with open(filename, 'r') as f:
        for nn,line in enumerate(f.readlines()):
            # 18 lines are the header
            if nn > 18:
                
                # Clean up each line
                line_list = line.strip().split('  ')

                # Save data
                radius_cm.append(float(line_list[0]))
                gas_density.append(float(line_list[1]))
                gas_temperature.append(float(line_list[2]))
                # TODO gas_opacity


    # Save a list of number of radial cells per time step
    Ncells = len(radius_cm)

    # Save in lists and reverse the radial direction
    radius_cm = np.array(radius_cm[::-1])
    gas_density = np.array(gas_density[::-1])
    gas_temperature = np.array(gas_temperature[::-1])
    # TODO all_opacities


    return radius_cm,gas_density,gas_temperature


# Take all time steps for model and make time-averaged data
def load_averagedarwindata(
        modelname:str='M2n315u6'
    ):
    """
    Loads all Darwin-data in a folder
    and takes time-average of all densities and temperatures

    ARGUMENTS
      modelname: string: Letter-number-code for the model.

    RETURNS
      average_radius: np.array with radial range in cm
      average_density: np.array with time-averaged radial dependant gas densities in g/cm3
      average_temperature: np.array with time-averaged radial dependant gas temperatures in K
    """

    # Given modelname, list all time-step-files
    filenames = os.listdir(f'../darwin_data/{modelname}_sel/')
    Ntimesteps = len(filenames)


    # Load all Darwin data
    Ncells = []
    all_radii = []
    all_density = []
    all_temperature = []
    # TODO all_opacities

    for filename in filenames:
        
        #Loop through files
        path = f'../darwin_data/{modelname}_sel/{filename}'

        # Reset lists
        radius_cm = []
        gas_density = []
        gas_temperature = []
        # TODO gas_opacity

        with open(path, 'r') as f:
            for nn,line in enumerate(f.readlines()):
                # 18 lines are the header
                if nn > 18:
                    
                    # Clean up each line
                    line_list = line.strip().split('  ')

                    # Save data
                    radius_cm.append(float(line_list[0]))
                    gas_density.append(float(line_list[1]))
                    gas_temperature.append(float(line_list[2]))
                    # TODO gas_opacity

        # Save a list of number of radial cells per time step
        Ncells.append(len(radius_cm))

        # Save in lists and reverse the radial direction
        all_radii.append(radius_cm[::-1])
        all_density.append(gas_density[::-1])
        all_temperature.append(gas_temperature[::-1])
        # TODO all_opacities
    


    # Average all data into one time-averaged set

    # Define grid to average into
    Naverage = max(Ncells)
    average_radius = np.linspace(
        min(min(all_radii)),
        max(max(all_radii)),
        Naverage
    )
    average_density = np.zeros(Naverage)
    average_temperature = np.zeros(Naverage)
    average_cellcounter = np.zeros(Naverage)
    # TODO average_opacity

    # Loop over all time steps
    for nstep in range(Ntimesteps):

        # Loop over all data in each time step
        for nn in range(Ncells[nstep]):

            # Check where in average_radius these data are
            for na in range(Naverage - 1):
                if all_radii[nstep][nn] >= average_radius[na] and all_radii[nstep][nn] < average_radius[na+1]:

                    # Sum all data at each radial point
                    average_density[na] += all_density[nstep][nn]
                    average_temperature[na] += all_temperature[nstep][nn]
                    # TODO average_opacity
                    average_cellcounter[na] += 1


    # Extract all zero-elements as sequence-based lists
    zero_lists = c3d.find_zeroelements(average_cellcounter)

    # Remove zeros in cell counter (to avoid infs, or NaN)
    average_cellcounter[np.where(average_cellcounter == 0)[0]] = 1

    # Average the data
    average_density /= average_cellcounter
    average_temperature /= average_cellcounter
    #TODO average_opacity /= average_cellcounter



    # Add average-data into the remaining "holes" in the averaged data sets
    for zero_list in zero_lists:

        # Make sure were not in the end of the array
        if zero_list[-1] != Naverage-1:

            # Take average of prev and next index (before and after "hole")
            average_density[zero_list] = 0.5*(average_density[zero_list[0]-1] + average_density[zero_list[-1]+1])
            average_temperature[zero_list] = 0.5*(average_temperature[zero_list[0]-1] + average_temperature[zero_list[-1]+1])
            # TODO average_opacity

        # For the end of the array, just take the final real value
        if zero_list[-1] == Naverage-1:
            average_density[zero_list] = average_density[zero_list[0]-1] 
            average_temperature[zero_list] = average_temperature[zero_list[0]-1]
            # TODO average_opacity


    return average_radius,average_density,average_temperature


# TODO
# func that computes rosseland opacities from density and temperatures
# and outputs array with density*opacity


# TODO
# Translate Darwindata to R3D-grid (and save files here?)
def darwin_to_radmc3d(
    darwin_radius:list,
    darwin_density:list,
    darwin_temperature:list,
    gridpath:str='../grid_distances.csv',
    amrpath:str='../amr_grid.inp'
):
    """
    Adapts (interpolates) 1D Darwin data to Radmc3d-grid and writes necessary 
    input files.

    ARGUMENTS
      darwin_radius:list-like: Distances to centre of star in cm (increasing order)
      darwin_density:list-like: Gas density * Rosseland opacities in g/cm2 * cm2/g
      darwin_temperature:list-like: Gas temperatures in Kelvin
      gridpath:str: Path to grid_distances.csv of your radmc3d-model
      amrpath:str: Path to amr_grid.inp of your radmc3d-model

    RETURNS
      dust_density_darwinstar.inp
      dust_temperature_darwinstar.dat
        Input files for Radmc3d, merge these with dust-files if necessary 
        (with c3d.merge_dustdensities and c3d.merge_dusttemperatures)
    """

    # Load r3d-grid
    print('Loading Radmc3d-grid.')
    r3d_griddistances = a3d.load_griddistances(
        gridpath=gridpath,
        amrpath=amrpath,
    )
    r3d_radius = r3d_griddistances[:,0]

    # Interpolate darwin-1d-data to r3d-radial grid
    # Default setting of interp is to keep first o last value outside the range
    r3d_darwindensity = np.interp(r3d_radius,darwin_radius,darwin_density)
    r3d_darwintemperature = np.interp(r3d_radius,darwin_radius,darwin_temperature)

    # Save in RADMC-3D inp-files
    nleafs = np.size(r3d_radius)
    progbar = 0

    print('Writing radmc3d-files:')
    with open('../dust_density_darwinstar.inp', 'w') as fdensity, \
        open('../dust_temperature_darwinstar.dat', 'w') as ftemperature:
        
        # Write headers:
        # 1
        # nleafs
        # number dust species
        fdensity.write(f'1\n{int(nleafs)}\n1\n')
        ftemperature.write(f'1\n{int(nleafs)}\n1\n')    

        # Write densities and temperatures

        # zip

        for nr3d in range(nleafs):
            fdensity.write(f'{r3d_darwindensity[nr3d]}\n')
            ftemperature.write(f'{r3d_darwintemperature[nr3d]}\n')

            # Some progress bar info
            if int(nr3d/nleafs*100) == 25 and progbar == 0:
                progbar += 1
                print('  Finished 25 per cent of the grid.')

            if int(nr3d/nleafs*100) == 50 and progbar == 1:
                progbar += 1
                print('  Finished 50 per cent of the grid.')

            if int(nr3d/nleafs*100) == 75 and progbar == 2:
                progbar += 1
                print('  Finished 75 per cent of the grid.')

    print(f'DARWIN Dust-star:\n    dust_density_darwinstar.inp\n    dust_temperature_darwinstar.dat\nDONE\n')

