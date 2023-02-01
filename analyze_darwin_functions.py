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
        modelname:str='M2n315u6'
    ):
    """
    Loads all Darwin-data in a folder
    and takes time-average of all data

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
# Translate Darwindata to R3D-grid (and save files here?)

def darwin_to_radmc3d():

    # Load r3d-grid
    # needs path to amr-file and radius-file

    # interpolate darwin-1d-data to r3d-radial grid

    # Multiply density with Rosseland-kappa!

    # save in r3d-style inp-files

    # Done!


    return 'hej'
