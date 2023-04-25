# functions and tools for loading and analyzing DARWIN data, and to create 
# RADMC3D-input models from DARWIN data.

# Import various libraries
import os
import math
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
Lsolcgs = 3.8280e33 # erg/s
sigma = 5.670374419e-8

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

    filename = f'../darwin_data/{modelname}_kappaross/model.{timestep}.rtk'

    radius_cm = []
    gas_density = []
    gas_temperature = []
    gas_opacity = []

    with open(filename, 'r') as f:
        for nn,line in enumerate(f.readlines()):
            # 18 lines are the header
            if nn > 8:
                
                # Clean up each line
                line_list = line.strip().split('  ')

                if len(line_list) > 1:
                    # Save data
                    radius_cm.append(float(line_list[0]))
                    gas_density.append(float(line_list[1]))
                    gas_temperature.append(float(line_list[2]))
                    gas_opacity.append(float(line_list[3]))

    # Save in arrays and reverse the radial direction to get growing radius
    radius_cm = np.array(radius_cm[::-1])
    gas_density = np.array(gas_density[::-1])
    gas_temperature = np.array(gas_temperature[::-1])
    gas_opacity = np.array(gas_opacity[::-1])

    return radius_cm,gas_density,gas_temperature,gas_opacity



# Translate Darwindata to R3D-grid (and save files here?)
def darwin_to_radmc3d(
        darwin_radius:list,
        darwin_density:list,
        darwin_temperature:list,
        darwin_opacity:list,
        input_radius:float,
        input_luminosity:float,
        gridpath:str='../grid_distances.csv',
        amrpath:str='../amr_grid.inp'
    ):
    """
    Adapts (interpolates) 1D Darwin data to Radmc3d-grid and writes necessary 
    input files.

    ARGUMENTS
      darwin_radius:list-like: Darwin-Distances to centre of star in cm (increasing order)
      darwin_density:list-like: Darwin-Gas density g/cm3
      darwin_temperature:list-like: Darwin-Gas temperatures in Kelvin
      darwin_opacity:list-like: Darwin-Rosseland opacities cm2/g

      input_radius:float: Rint, radius of Darwin-star in meters
      input_luminosity:float: Lext, (bol) luminosity of star in Watts

      gridpath:str: Path to grid_distances.csv of your radmc3d-model
      amrpath:str: Path to amr_grid.inp of your radmc3d-model

      
    RETURNS
      dust_density_darwinstar.inp
      dust_temperature_darwinstar.dat
        Input files for Radmc3d, merge these with dust-files if necessary 
        (with c3d.merge_dustdensities and c3d.merge_dusttemperatures)
    """

    # Compute effective temperature
    effective_temperature = (input_luminosity / (4*np.pi * input_radius**2 * sigma) )**0.25
    print(f'Effective temperature is: {effective_temperature} K')

    # Load r3d-grid
    print('Loading Radmc3d-grid.')
    r3d_griddistances = a3d.load_griddistances(
        gridpath=gridpath,
        amrpath=amrpath,
    )
    r3d_radius = r3d_griddistances[:,0]


    # Create density-opacity array for r3d-sims, since the gas-density
    # in r3d is density*rosselandopacity, and the gas kappa is just
    # kappaabs = 1 and kappascat = 0
    #
    # First check if density-opacity input data are np-arrays
    if type(darwin_density) == list:
        darwin_density = np.array(darwin_density)
    if type(darwin_opacity) == list:
        darwin_opacity = np.array(darwin_opacity)
    if type(darwin_temperature) == list:
        darwin_temperature = np.array(darwin_temperature)


    # Then write new array and divide both with number of r3d-cells per spherical shell

    # 0
    #print('normal and T-cut')
    #darwin_densityopacity = darwin_density * darwin_opacity
    #darwin_temperature[
    #    np.where(darwin_temperature >= effective_temperature)[0]
    #] = effective_temperature

    # 9
    print('rho-kappa divide by R and T-cut')
    darwin_densityopacity = darwin_density * darwin_opacity * (AUcm / darwin_radius)**2
    darwin_temperature[
        np.where(darwin_temperature >= effective_temperature)[0]
    ] = effective_temperature

    # 1
    #print('/r^2')
    #darwin_densityopacity = darwin_density * darwin_opacity * (AUcm / darwin_radius)**4
    #darwin_temperature *= (AUcm / darwin_radius)**2

    # 2
    #print('divide by r')
    #darwin_densityopacity = darwin_density * darwin_opacity * (AUcm / darwin_radius)**2
    #darwin_temperature *= (AUcm / darwin_radius)

    # 3
    # only divide rho-kappa by r4

    # 4
    #print('divide by Ncells')
    # Load number of cells per radial Darwin bin
    #Ncells = np.zeros(len(darwin_radius))
    #
    ## First number of cells are all inside smallest darwin cell
    #Ncells[0] = np.size(np.where(
    #    r3d_radius <= darwin_radius[0])[0]
    #)
    #for nn in range(1,len(darwin_radius)):
    #    Ncells[nn] = np.size(np.where(
    #        (r3d_radius >= darwin_radius[nn-1]) & (r3d_radius <= darwin_radius[nn]))[0]
    #    )
    #
    #    # Remove zeros and infs
    #    if Ncells[nn] == 0 or Ncells[nn] == math.inf:
    #        Ncells[nn] = 1
    #darwin_densityopacity = darwin_density * darwin_opacity / Ncells

    # 5
    #print('divide by delta R')
    #deltaRau = np.zeros(len(darwin_radius))
    #deltaRau[0] = darwin_radius[0]/AUcm
    #for nn in range(1,len(darwin_radius)):
    #    deltaRau[nn] = (darwin_radius[nn] - darwin_radius[nn-1])/AUcm
    #darwin_densityopacity = darwin_density * darwin_opacity / deltaRau
    #darwin_temperature /= deltaRau

    # 6
    #print('Divide T with R^2')
    #darwin_densityopacity = darwin_density * darwin_opacity
    #darwin_temperature *= (AUcm / darwin_radius)**2

    # 7
    #print('Divide T with R')
    #darwin_densityopacity = darwin_density * darwin_opacity
    #darwin_temperature *= (AUcm / darwin_radius)



    
    # Interpolate darwin-1d-data to r3d-radial grid
    # Default setting of interp is to keep first o last value outside the range
    r3d_darwindensityopa = np.interp(r3d_radius,darwin_radius,darwin_densityopacity)
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
            fdensity.write(f'{r3d_darwindensityopa[nr3d]}\n')
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

    # TEMP! Just to check!
    #return darwin_densityopacity