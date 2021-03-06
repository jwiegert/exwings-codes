import analyze_r3d_functions as a3d
import numpy as np
import os

# Import BH-codes
# See manual for RADMC3D for references on these codes.
from bhmie.makedustopac import *

def create_kappaabs(
        wavelengthpath:str='../wavelength_micron.inp',
        optconstlist:list=['mg2sio4'],
        agrainlist:list=[0.1],
        matdens:float=2
    ):

    # Number of "species" is
    nrspec = len(optconstlist)*len(agrainlist)

    # Load wavelength-grid
    wavelengths,nwave = a3d.load_wavelengthgrid(wavelengthpath)

    # Change units from um to cm and change list to np.array
    wavelengths = np.array([wavelength*1e-4 for wavelength in wavelengths])

    # Logarithmic spread parameter TODO not tested yet! works fine for 10 sizes
    logawidth = max(agrainlist)/len(agrainlist)

    # Go through the different dust species
    for optconst in optconstlist:

        # Load dust species
        optconst_file = f'../bhmie/lnk/{optconst}.lnk'

        # Load grain size grid
        for agrainum in agrainlist:

            #Change units to cm
            agraincm = agrainum*1e-4

            # Create all opacity files
            opacity = compute_opac_mie(
                optconst_file=optconst_file,
                matdens=matdens,
                agraincm=agraincm,
                lamcm = wavelengths,
                logawidth=logawidth,
                wfact=3,
                na=20,
                extrapolate=True,
                verbose=False
            )

            # Example code and info from makedustopac and bhmie for reference here
            #compute_opac_mie(optconst_file,matdens,agraincm,lamcm,
            #                     theta=None,logawidth=None,wfact=3.0,na=20,
            #                     chopforward=0.0,errtol=0.01,verbose=False,
            #                     extrapolate=False):
            """
                matdens       = Material density in g/cm^3
                agraincm      = Grain radius in cm
                lamcm         = Wavelength grid in cm (a numpy array)
                theta         = Optional angular grid (a numpy array) between 0 and 180
                                which are the scattering angle sampling points at 
                                which the scattering phase function is computed.
                logawidth     = Optional: if set, the size agrain will instead be a 
                                sample of sizes around agrain. This helps to smooth out
                                the strong wiggles in the phase function and opacity
                                of spheres at an exact size. Since in Nature it rarely
                                happens that grains all have exactly the same size, this
                                is quite natural. The value of logawidth sets the width
                                of the Gauss in ln(agrain), so for logawidth<<1 this
                                give a real width of logawidth*agraincm. 
                wfact         = (default=3.0) Grid width of na sampling points in units
                                of logawidth. The Gauss distribution of grain sizes is 
                                cut off at agrain * exp(wfact*logawidth) and
                                agrain * exp(-wfact*logawidth).
                na            = (default=20) Number of size sampling points 
                                (if logawidth set).
                extrapolate   = If set to True, then if the wavelength grid lamcm goes
                                out of the range of the wavelength grid of the 
                                optical constants file, then it will make a suitable
                                extrapolation: keeping the optical constants constant
                                for lamcm < minimum, and extrapolating log-log for
                                lamcm > maximum.
            """

            # Write dustkappa_*.inp
            write_radmc3d_kappa_file(opacity,f'{optconst}_{agrainum}')

            # Move kappaabs files
            os.system(f'mv dustkappa_{optconst}_{agrainum}.inp ../r3dsims/opacities/')

    # Write dustopac.inp, list of all species and names
    with open('dustopac.inp','w') as f:

        # Write header
        f.write(f'2\n{nrspec}\n-----------------------------\n')

        # Write each species and grain size
        for optconst in optconstlist:
            for agrainum in agrainlist:
                f.write(f'1\n0\n{optconst}_{agrainum}\n-----------------------------\n')

    # Move opac-file to subfolder of simulation folder
    os.system(f'mv dustopac.inp ../r3dsims/opacities')
