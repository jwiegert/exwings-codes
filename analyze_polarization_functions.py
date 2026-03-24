# Functions for analysing data from synthetic polarisation-images
# Import as apf for usage
# ------------------------------------------------------------ #
# Various useful packages
import matplotlib.pyplot as plt
import numpy as np
#import scipy
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
# apf.plot_onepolarization(
#    imagepath:str='../r3dresults_vltcompare/st28gm06n052_nospikes/166_optoolmgsio/',
#    imagewave:float=0.65,
#    polarization:str='l',
#    distance:float=1,
#)
#
# ------------------------------------------------------------ #
# Functions be here

def write_r3d_emptypolarization(
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


def plot_intensity_onepolarization(
        imagepath:str='../r3dresults_polarized/st28gm06n052_images/166/',
        imagewave:float=0.65,
        polarization:str='l',
        distance:float=1,
        vmaxcorr:float=1e-2
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
      vmaxcorr:float change max flux density limit for images

    RETURNS
      fix,ax,fluxtotal
        fig and axes-objects for plotting
        list of all subplots total flux densities in Jy set at distance.
    """
    # Automatically add / to end of path if it's missing
    if imagepath[-1] != '/':
        imagepath += '/'
        #
    # Check and set name of chosen polarization
    if polarization == 'i' or polarization == 'I':
        pol_name = 'Total intensity'
    elif polarization == 'l' or polarization == 'L':
        pol_name = 'Linear-pol intensity'
    elif polarization == 'v' or polarization == 'V':
        pol_name = 'Circular-pol intensity'
    else:
        raise ValueError('Incorrect polarisation, aborting!')
        #
    # Extract all images of certain wavelength
    images = []
    for file in os.listdir(imagepath):
        if re.findall(str(imagewave), file):
            images.append(file)
    Nimages = len(images)
    #
    # Initiate subplot object
    nrows = 2
    ncols = int(np.ceil(Nimages/2))
    fig,ax = plt.subplots(
        nrows,ncols,
        figsize=(6,6),
    )
    # Setup various lista and arrays
    angles=[]
    fluxtotal = []
    #
    for nimage,image in enumerate(images):
        # Set image frame position from nimage
        if nimage < ncols:
            nrow = 0
            ncol = nimage
        else:
            nrow = 1
            ncol = nimage - ncols
        intensity_1d = []        #
        # Extract angle names
        angles.append(
            image.split('_')[1]+' '+image.split('_')[2]
        )
        # Open image file
        with open(f'{imagepath}{image}', 'r') as fimage:
            # Loop through image
            for nl,line in enumerate(fimage.readlines()):
                #
                # Check if polarised image or not
                if nl == 0:
                    imagesetting = int(line.strip())
                    #
                # Then check polarisation setting
                if imagesetting == 3:
                    # 3 means full polarisation, continue with extracting data
                    #
                    # Row 1: image size, pixels by pixels
                    if nl == 1:
                        npixels_x = int(line.split()[0])
                        npixels_y = int(line.split()[1])
                        npixels = max([npixels_x,npixels_y])
                    # row 3: pixel size is in cm, divide by AUcm for AU
                    if nl == 3:
                        pixelsize_au = float(line.split()[0])/AUcm
                        #
                    # row 6 onward: pixel-values, four columns
                    # image.out's pixels has unit
                    # erg s-1 cm-2 Hz-1 ster-1
                    if nl > 5:
                        if len(line.split()) > 1:
                            # Loop over polarisation numbers
                            # Split line and extract each column, strip from \n and make floats
                            i_pixeldata = float(line.split()[0].strip())
                            q_pixeldata = float(line.split()[1].strip())
                            u_pixeldata = float(line.split()[2].strip())
                            v_pixeldata = float(line.split()[3].strip())
                            #
                            # Save chosen polarization data
                            if polarization == 'i' or polarization == 'I':
                                intensity_1d.append(
                                    i_pixeldata
                                )
                            elif polarization == 'l' or polarization == 'L':
                                intensity_1d.append(
                                    np.sqrt(q_pixeldata**2 + u_pixeldata**2)
                                )
                            elif polarization == 'v' or polarization == 'V':
                                intensity_1d.append(
                                    np.sqrt(v_pixeldata**2)
                                )
        # Continue with image plotting
        #
        # Abort if imagesetting is incorrect
        if imagesetting == 1:
            raise ValueError('Not polarised, load normal image with a3d.load_images()!')
        else:
            # Continue with plotting image
            #
            # Extract some useful quantities
            # pixel size in asec (pixelsize in au and distance in pc gives distance in asec)
            pixelsize_as = pixelsize_au / distance
            # Size of whole image in AU and image-axis-scales
            size_au = pixelsize_au * npixels
            axisplot  = [-0.5*size_au,0.5*size_au,-0.5*size_au,0.5*size_au]
            #
            # Total flux density of the image in Jy
            # Transform to Jy/pix and sum all
            # 1 Jy = 1e23 erg/(s cm2 Hz)
            # 1 asec = 1/(180/pi * 3600)^2 ster = 2.35044305391e-11 ster
            # 1 pixel = pixelsize_as^2  asec^2
            # So in effect unit is Jy/30x30AU2
            fluxtotal.append(
                sum(intensity_1d) * 1.e23 * 2.35044305391e-11 * pixelsize_as**2
            )
            # Create empty 2D array for the image and pixel counters
            image2d = np.zeros((npixels,npixels))
            nx,ny = 0,0
            #
            for pixelflux in intensity_1d:
                # Recompute unit to Jy/asec2
                image2d[nx,ny] = pixelflux * 1.e23 * 2.35044305391e-11 / distance**2
                #   
                # Move nx and ny
                nx = nx + 1
                if nx == npixels_x:
                    nx = 0
                    ny = ny + 1
            #
            # Plot each image
            #
            ax[nrow,ncol].imshow(
                image2d,
                origin='lower',
                extent=axisplot,
                cmap=plt.get_cmap('hot'),
                vmin=np.mean(image2d),
                vmax=np.max(image2d)*vmaxcorr
            )
            # Write out angle choice
            ax[nrow,ncol].set_title(
                angles[nimage], fontsize=14
            )
            # Set axis labels and tick sizes
            if nrow == 1:
                ax[nrow,ncol].set_xlabel(
                    'Offset (au)',fontsize=14
                )
            if ncol == 0:
                ax[nrow,ncol].set_ylabel(
                    'Offset (au)',fontsize=14
                )
            ax[nrow,ncol].tick_params(
                axis='both', which='major', labelsize=15
            )
            ax[nrow,ncol].tick_params(
                axis='both', which='major', labelsize=15
            )
    # Write out polarisation of the images
    fig.suptitle(
        pol_name, fontsize=14
    )
    return fig,ax,fluxtotal


# Plot a subplot with image intensity, linear intensity and circular polarization intensity
def plot_intensity_allpolarizations(
        imagepath:str='../r3dresults_polarized/st28gm06n052_images/166/',
        imagefile:str='image_i000_phi000_0.65um.out',
        distance:float=1,
        vmaxcorr:float=1e-1,
    ):
    """
    TODO

    vmaxcorr:float change max flux density limit for images
    """
    # Automatically add / to end of path if it's missing
    if imagepath[-1] != '/':
        imagepath += '/'
        #
    # Extract imagename for plot later
    imagename = f'{imagefile.split('_')[1]}{imagefile.split('_')[2]} {imagefile.split('_')[3][:-4]}'
    #
    # List polarisations
    polarization = [
        'Tot intensity',
        'Linear intensity',
        'Circular intensity',
    ]
    # Initiate subplot object
    fig,ax = plt.subplots(
        1,3,
        figsize=(9,3),
    )
    # Declare lists for image intensity, linear pol intensity and circ pol intensity
    # and for imageplotting
    I_intens_1d = []
    L_intens_1d = []
    V_intens_1d = []
    # Load data
    with open(imagepath+imagefile, 'r') as f:
        for nl,line in enumerate(f.readlines()):
            #
            # Check if polarised image or not
            if nl == 0:
                imagesetting = int(line.strip())
            if imagesetting == 3:
                # This is correct, otherwise raise error (see below)
                #            
                # Row 1: image size, pixels by pixels
                if nl == 1:
                    npixels_x = int(line.split()[0])
                    npixels_y = int(line.split()[1])
                    npixels = max([npixels_x,npixels_y])            
                # row 3: pixel size is in cm, divide by AUcm for AU
                if nl == 3:
                    pixelsize_au = float(line.split()[0])/AUcm
                #
                # row 6 onward: pixel-values, four columns
                # image.out's pixels has unit
                # erg s-1 cm-2 Hz-1 ster-1
                if nl > 5:
                    if len(line.split()) > 1:
                        # Loop over polarisation numbers
                        # Split line and extract each column, strip from \n and make floats
                        i_pixeldata = float(line.split()[0].strip())
                        q_pixeldata = float(line.split()[1].strip())
                        u_pixeldata = float(line.split()[2].strip())
                        v_pixeldata = float(line.split()[3].strip())
                        # Save each image-combo
                        I_intens_1d.append(
                            i_pixeldata
                        )
                        L_intens_1d.append(
                            np.sqrt(q_pixeldata**2 + u_pixeldata**2)
                        )
                        V_intens_1d.append(
                            np.sqrt(v_pixeldata**2)
                        )
    # Continue with plotting images
    if imagesetting == 1:
        # If not polarized, abort
        raise ValueError('Not polarised, load normal image with a3d.load_images()!')
    else:
        # Continue with plotting image
        #
        # Size of whole image in AU and image-axis-scales
        size_au = pixelsize_au * npixels
        axisplot  = [-0.5*size_au,0.5*size_au,-0.5*size_au,0.5*size_au]
        #
        # Create 2D arrays
        I_image2d = np.zeros((npixels,npixels))
        L_image2d = np.zeros((npixels,npixels))
        V_image2d = np.zeros((npixels,npixels))
        # Set pixel counter
        nx,ny = 0,0
        #
        for i_flux,l_flux,v_flux in zip(I_intens_1d,L_intens_1d,V_intens_1d):
            # Convert image1d to 2d and change unit to Jy/asec2, ie normalised to distance
            I_image2d[nx,ny] = i_flux * 1.e23 * 2.35044305391e-11 / distance**2
            L_image2d[nx,ny] = l_flux * 1.e23 * 2.35044305391e-11 / distance**2
            V_image2d[nx,ny] = v_flux * 1.e23 * 2.35044305391e-11 / distance**2
            # Move nx and ny
            nx = nx + 1
            if nx == npixels_x:
                nx = 0
                ny = ny + 1
        # Plot all three sub plots
        ax[0].imshow(
            I_image2d,
            origin='lower',
            extent=axisplot,
            cmap=plt.get_cmap('hot'),
            vmin=np.mean(I_image2d),
            vmax=np.max(I_image2d)*vmaxcorr
        )
        ax[0].set(
            title=polarization[0], 
            xlabel='Offset (AU)',
            ylabel='Offset (AU)'
        )
        ax[1].imshow(
            L_image2d,
            origin='lower',
            extent=axisplot,
            cmap=plt.get_cmap('hot'),
            vmin=np.mean(L_image2d),
            vmax=np.max(L_image2d)*vmaxcorr
        )
        ax[1].set(
            title=polarization[1],
            xlabel='Offset (AU)',
        )
        ax[2].imshow(
            V_image2d,
            origin='lower',
            extent=axisplot,
            cmap=plt.get_cmap('hot'),
            vmin=np.mean(V_image2d),
            vmax=np.max(V_image2d)*vmaxcorr
        )
        ax[2].set(
            title=polarization[2],
            xlabel='Offset (AU)',
        )
    fig.suptitle(
        imagename, fontsize=14
    )
    return fig,ax


