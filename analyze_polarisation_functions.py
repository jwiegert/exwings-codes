# Functions for analysing data from synthetic polarisation-images
# Import as apf for usage
# ------------------------------------------------------------ #
# Various useful packages
import matplotlib.pyplot as plt
import numpy as np
import scipy
import re
import os

# My own functions
import analyze_r3d_functions as a3d

# Basic definitions
c = 2.998e8                 # m/s
AUcm = 1.49598e13           # cm
Rsolcm = 6.955e10           # cm
Lsol = 3.828e26             # W
Rsolmeter = 6.955e8         # m
Msol = 1.989e33             # g
radian_mas = 206264800      # milliasec
baselineVLTI = 201.92       # metres
beamwidthVLTI_10um = 5.1e-3 # asec




#
# ------------------------------------------------------------ #
# List of functions
#
# apf.write_r3d_emptypolarisation(
#    wavelengthpath:str='../r3dsims/wavelength_micron.inp',
#    outputpath:str='../dustkapscatmat_opastar.inp'
# )
#
#
#
# ------------------------------------------------------------ #
# Functions be here

def write_r3d_emptypolarisation(
        wavelengthpath:str='../r3dsims/wavelength_micron.inp',
        outputpath:str='../dustkapscatmat_opastar.inp'
    ):
    """
    Write empty kappascatmat-file for usage for gas or other medium with no scattering
    when combined with dust species that uses the full mueller matrix for scattering in
    Radmc-3d.

    OUTPUT FILE STRUCTURE
    # header
    # 1, 1000, 181 (one species, 1000 wavelengths, 181 angles)
    #
    # 4 columns, 1000 long - wavelength(um), kappaabs(cm2/g), kappascat(cm2/g), angle
    #                          as is            ones             zeros           zeros
    #
    # 1 column, angles, 0 to 180, keep as is
    #
    # Z-matrix elements, one 4x4 matrix per wavelength AND angle, but only 6 elements from the matrix per line
    #   Z11  Z12   0    0
    #    0   Z22   0    0
    #    0    0   Z33  Z34
    #    0    0    0   Z44
    # but in this order:   Z11  Z12  Z22  Z33  Z34  Z44
    # so 1000x181 lines with six 0.0 per line.
        
    ARGUMENTS
      wavelengthpath:str path to r3d-wavelength-input file including filename; wavelength_micron.inp
      outputpath:str path to when outfile is put, including filename, eg dustkapscatmat_opastar.inp
    
    RETURNS
      dustkapscatmat_opastar.inp-file containing empty mueller matrices and no scattering and
      kappa_abs = 1 for all wavelengths. 
    """
    #
    # Load wavelength grid file and use for scatmat-file
    wavelengths,nwave = a3d.load_wavelengthgrid(
        path=wavelengthpath
    )
    # Number of matrices is then (one for each angle)
    Nmatrices = nwave*181
    #
    # Create and write file
    with open(outputpath, 'w') as fstar:
        #
        # Write Header
        fstar.writelines('#============================================================================\n')
        fstar.writelines('# Special opacity file for gas contents in co5bold data\n')
        fstar.writelines('# that includes only absorption (1 cm2/g) and no scattering\n')
        fstar.writelines('# and empty Mueller matrix to include no Stokes scattering but\n')
        fstar.writelines('# to be compatible with existing files with polarised scattering\n')
        fstar.writelines('# for the dust.\n')
        fstar.writelines('#============================================================================\n')
        fstar.writelines(f'       1\n    {nwave:0d}\n     181\n\n')
        #
        # Write wavelength grid, absorption=1, scattering=0, scat-angles = 0
        for wavelength in wavelengths:
            fstar.writelines(f'  {wavelength:.6E}  1.0  0.0  0.0\n')
        fstar.writelines('\n')
        #
        # Write angle grid
        for angle in range(181):
            fstar.writelines(f'  {angle}\n')
        fstar.writelines('\n')
        #
        # Write empty Mueller matrix
        for nmatrix in range(Nmatrices):
            fstar.writelines('  0.0  0.0  0.0  0.0  0.0  0.0\n')
    # Done!






def plot_onepolarisation(
        imagepath:str='../r3dresults_vltcompare/st28gm06n052_nospikes/166_optoolmgsio/',
        imagewave:float=0.65,
        polarization:str='l',
        distance:float=1,
    ):
    """
    Loads all angles of one snapshot and one wavelength and
    plots as 6 subplots

    ARGUMENTS
      imagepath:str path to folder that contains images you want to plot
      imagewave:float wavelength of image to plot in micrometres
      polarization:str Choice of polarization to plot intensity of, either 
        'i' : full intensity = Q2 + U2 + V2 (ie normal image, all incident light)
        'l' : all linear = sqrt(Q2 + U2)
        'v' : circular = sqrt(V2)
      distance:float normalisation to flux densities of image, in parsec

    RETURNS
      fig and axes-objects for plotting

    
    """

    # extract all images of certain wavelength
    # image_i000_phi000_0.65um.out

    images = []
    for file in os.listdir(imagepath):
        if re.findall(str(imagewave), file):
            images.append(file)

    Nimages = len(images)

    # Initiate subplot object
    fig,ax = plt.subplots(
        Nimages,1, 
        figsize=(6,6),
    )

    # TODO finish this

    return fig,ax



def plot_allpolarisations():
    print('hej')



