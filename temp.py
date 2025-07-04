# Load and check one figure
import analyze_r3d_functions as a3d
import numpy as np
import os

models = [
    'st28gm06n052_timedep',
    'st28gm06n074',
    'st28gm06n075'
]
angles = [
    'i000_phi000',
    'i090_phi000',
    'i090_phi090',
    'i090_phi270',
    'i180_phi000',
    'i270_phi000',
]

# Loop through models
for model in models:

    # Load stellar average flux at 10um
    starSEDs = np.loadtxt(f'../r3dresults/{model}_nodust/average_sed.dat')
    tenmicron = np.where(starSEDs[:,0] > 10)[0][0]-1
    startenflux = starSEDs[tenmicron,1]

    # Declare empty list
    fluxsums = []

    # Extract snapshot numbers and loop through them
    path = f'../r3dresults/{model}_nospikes/'
    snapshots = [int(filename) for filename in os.listdir(path) if os.path.isdir(path+filename)]

    for snapshot in snapshots:

        # Loop through angles
        for angle in angles:

            # Load image
            image2d,image2dlog,totalflux,axisplot = a3d.load_images(
                path=f'../r3dresults/{model}_nospikes/{snapshot:03d}',
                image=f'image_{angle}_10um.out'
            )
            # Size of image and pixels
            Npix = np.shape(image2d)[0]
            pixelsize = 2*axisplot[1] / Npix

            # Radius of "stellar circle" in pixels
            Rin_pixels = 2/pixelsize
            # And radius of outer part of important dust region
            Rout_pixels = 8/pixelsize

            # Loop through iamge and save total flux outside the star
            fluxsum = 0
            for nx in range(Npix):
                for ny in range(Npix):
                    # Check if pixel is within stellar circle
                    pixelcoord = np.sqrt((nx-0.5*Npix)**2 + (ny-0.5*Npix)**2)

                    # Sum all pixels outside this
                    # flux in Jy/pixel is image2d * pixelsize(asec2)
                    if pixelcoord > Rin_pixels and pixelcoord < Rout_pixels:
                        fluxsum += image2d[nx,ny] * pixelsize**2

            fluxsums.append(fluxsum)
    # Take average of all and compare with average flux density of 10um from the star
    print(f'{model}')
    print(f'    Mean contrast (dust/star):   {np.mean(fluxsums)/startenflux}')
    print(f'    Median constrat (dust/star): {np.median(fluxsums)/startenflux}')
    print(f'    Mean dust flux: {np.mean(fluxsums)}. Median dust flux: {np.median(fluxsums)}')

# OUTPUT notes (20min per model)
# st28gm06n052_timedep
#     Mean contrast (dust/star):   0.5951135548257981
#     Median constrat (dust/star): 0.5877750860782606
#     Mean dust flux: 12367792.552373314. Median dust flux: 12215282.735741451
# 
# st28gm06n074
#     Mean contrast (dust/star):   0.515662117720173
#     Median constrat (dust/star): 0.49675065351239295
#     Mean dust flux: 10775684.195016574. Median dust flux: 10380495.254496064
#
# st28gm06n075
#     Mean contrast (dust/star):   0.39002486640644257
#     Median constrat (dust/star): 0.387207184550776
#     Mean dust flux: 8121845.315783781. Median dust flux: 8063170.143629588

