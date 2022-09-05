# Plot images

import analyze_r3d_functions as a3d
import matplotlib.pyplot as plt

path = '../r3dresults/st28gm06n056/'
phases = [140,141,142]

# List of image file names:
# Yes, I could use list compr, but this is easier to read
images = []
for phase in phases:
    images.append(f'../r3dresults/st28gm06n056/{phase}/image_i0_10um.out')

# Plot images
fig,ax = testflux = a3d.plot_imagesubplots(
    imagelist = images,
    distance = 1,
    scale = 'lin'
)


plt.savefig(
    '20220905_st28gm06n056_allphases_10um.png', dpi=300, bbox_inches='tight'
)

