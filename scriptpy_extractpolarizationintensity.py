# Extract all integrated polarization intensities
# and write to file
import numpy as np
import os

models = [
    'st28gm06n052',
    'st28gm06n074',
    'st28gm06n075',
]
wavelengths = [
    '0.65',
    '0.82',
]
distance = 1
angles = [
    'i000_phi000',
    'i090_phi000',
    'i090_phi090',
    'i090_phi270',
    'i180_phi000',
    'i270_phi000',
]
for wavelength in wavelengths:
    for model in models:

        modelpath = f'../r3dresults_polarized/{model}_images'
        modelname = model

        # Automatically add / to end of path if it's missing
        if modelpath[-1] != '/':
            modelpath += '/'

        # Extract all snapshot numbers for this model
        snapshots = [int(filename) for filename in os.listdir(modelpath) if os.path.isdir(modelpath+filename)]
        snapshots.sort()
        Nsnapshots = len(snapshots)

        # Open file to write all coords
        with open(f'{modelpath}polarizationintensity_{wavelength}um.dat','w') as fpol:
            print(
                f'  Writing polarization file\n    {modelpath}polarizationintensity_{wavelength}um.dat'
            )
            # Write header
            fpol.writelines(f'# Polarization intensity at lambda = {wavelength} um for {modelname}\n')
            fpol.writelines(f'Nsnapshots={Nsnapshots}\n')
            fpol.writelines(f'# Angles : \n')
            for angle in angles:
                fpol.writelines(f'    {angle}')
            fpol.writelines('\n# Snapshot      Q U V L for each angle in Jy at 1pc - 2 spaces between intensity and 4 spaces between angles.')
            #
            # Loop over time
            for snapshot in snapshots:
                # Print some output
                if snapshot == 100 or snapshot == 200 or snapshot == 300 or snapshot == 400:
                    print(f'  At snapshot {snapshot}')
                #
                # Print snapshot number to dat file
                fpol.writelines(f'\n{snapshot:03d}    ')
                #
                # Define image folder path
                imagepath = f'{modelpath}{snapshot}/'
                #
                # Loop over angles
                for angle in angles:
                    # Define image file name
                    image = f'image_{angle}_{wavelength}um.out'
                    # Define lists for polarization intensities
                    q_intensity = []
                    u_intensity = []
                    v_intensity = []
                    l_intensity = []
                    # Read image
                    with open(f'{imagepath}{image}', 'r') as fimage:
                        # Loop through image
                        for nl,line in enumerate(fimage.readlines()):
                            # Check if polarized image or not
                            if nl == 0:
                                imagesetting = int(
                                    line.strip()
                                )
                            # Then check polarisation setting
                            if imagesetting == 3:
                                # 3 means full polarisation, continue with extracting data
                                #
                                # row 6 onward: pixel-values, four columns
                                # image.out's pixels has unit
                                # erg s-1 cm-2 Hz-1 ster-1
                                if nl > 5:
                                    if len(line.split()) > 1:
                                        # Loop over polarisation numbers
                                        # Split line and extract each column, strip from \n and make floats
                                        # and change unit to Jy at chosen distance
                                        q_pixeldata = float(line.split()[1].strip()) * 1.e23 * 2.35044305391e-11 / distance**2
                                        u_pixeldata = float(line.split()[2].strip()) * 1.e23 * 2.35044305391e-11 / distance**2
                                        v_pixeldata = float(line.split()[3].strip()) * 1.e23 * 2.35044305391e-11 / distance**2
                                        #
                                        # Save intensities in lists
                                        #
                                        q_intensity.append(
                                            np.sqrt(q_pixeldata*q_pixeldata)
                                        )
                                        u_intensity.append(
                                            np.sqrt(u_pixeldata*u_pixeldata)
                                        )
                                        v_intensity.append(
                                            np.sqrt(v_pixeldata*v_pixeldata)
                                        )
                                        l_intensity.append(
                                            np.sqrt(q_pixeldata*q_pixeldata + u_pixeldata*u_pixeldata)
                                        )
                            else:
                                raise ValueError(
                                    'Not polarized image!'
                                )
                    # Extract total polarized intensity
                    q_total = np.sum(q_intensity)
                    u_total = np.sum(u_intensity)
                    v_total = np.sum(v_intensity)
                    l_total = np.sum(l_intensity)
                    #
                    # Print polarization intensities to file
                    fpol.writelines(
                        f'{q_total:.3e}  {u_total:.3e}  {v_total:.3e}  {l_total:.3e}    '
                    )

        # 8 snapshots: 1min
        # 7.125sek per snapshot
        # 400 tar 50min
        # 450 tar 55min
        # så 3h för 3 modeller 1 våglängd
        # 6h för allt

