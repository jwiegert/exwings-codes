# Plots various figures for first co5bold-r3d-paper
import matplotlib.pyplot as plt
from scipy import ndimage
import numpy as np
import scipy.ndimage
import re
import os

from matplotlib import rc
from mpl_toolkits.axes_grid1 import make_axes_locatable

import analyze_co5bold_functions as a5d
import analyze_r3d_functions as a3d
import create_r3d_functions as c3d

# Figure settings
rc('font',**{'family':'serif','serif':['serif']})
rc('text', usetex=True)
rc('xtick.major',size=8)
rc('xtick.minor',size=4)
rc('ytick.major',size=8)
rc('ytick.minor',size=4)

# Constants
c = 2.998e8 # m/s
AUcm = 1.49598e13 # cm
Lsol = 3.828e26 # W
Rsolmeter = 6.955e8 # m
cubesize = 222757675648155.62/AUcm # Yes, hardcoded to large grid, change if needed for 
                                   # other grids, ie, from amr_grid, first coordinate 
                                   # is courner coordinate of first base cell
radian = 206264800 # milliasec
#baselineVLTI = 130.23 # metres
baselineVLTI = 201.92 # metres
diameterJWST = 6.5 # metres

# Dust model snapshots labels
phasetimes = [
    '29.95 yrs'
]


# Plot choices

# Processinginfo
plot_coboldgrid = 'n'
plot_opticalthickness = 'n'
list_smoothingchanges = 'n'
plot_2dslices = 'y'

# Grain properties
plot_grainsizehist = 'n'
plot_grainsizeradius = 'n'
plot_absscat = 'n'
plot_temperaturecompare = 'n'

# Plot co5bold-figs
plot_coboldsed = 'n'
plot_images_examples = 'n'

# Plot symmetric figs
plot_darwin_imagesed = 'n'
plot_darwin_comparesed = 'n'
plot_point_imagesed = 'n'

# Merge contour and images, only t2
plot_images_convolved_vlti = 'n'

# Observables
compute_luminosities = 'n'
compute_tenmicronfluxdensities = 'n'
compute_spectralindeces = 'n'
measuredustcloudflux = 'n'
plot_resolutiondistance = 'n'
check_smoothedimage_radius = 'n'
plot_smoothedimage_radius = 'n'



# ----------------------------------------------------------------
# FIG Cut through of CO5BOLD grid for st28gm06n052 with cell 
# sizes against distance from centre of the grid.

if plot_coboldgrid == 'y':

    c5dgrid,cellcourners,cellsize = a5d.load_grid_properties(
        savpath='../../exwings_archivedata/co5bold_data/dst28gm06n052/st28gm06n052_186.sav'
    #    savpath='../../exwings_archivedata/co5bold_data/dst28gm06n056/st28gm06n056_140.sav'
    )
    Ncourners = cellcourners[:,0].size - 1
    cellsizes = np.zeros(Ncourners)
    cellcoords = np.zeros(Ncourners)

    for nn in range(Ncourners):
        cellsizes[nn] = (cellcourners[nn+1,0] - cellcourners[nn,0])/AUcm
        cellcoords[nn] = (0.5*(cellcourners[nn+1,0] + cellcourners[nn,0]))/AUcm


    fig, ax = plt.figure('grid-cutthrough', figsize=(6, 4)), plt.axes()

    ax.plot(cellcoords,cellsizes,'.',markersize=2)

    ax.set(ylim=[0,0.14])
    ax.set_xlabel(r'Distance from centre along one axis (au)',fontsize=18)
    ax.set_ylabel(r'Cell size (au)',fontsize=18,)
    ax.tick_params(axis='both', which='major', labelsize=15)

    fig.tight_layout()
    fig.savefig('figs/grid-cutthrough.pdf', dpi=300, facecolor="white")
    fig.show()

# ----------------------------------------------------------------
#
# Plot Figure with average LOS optical thickness of gasmodel

if plot_opticalthickness == 'y':

    Reffective,fig,ax = a3d.plot_opticalthick(
        path='../r3dresults/st28gm06n056/140/'
    )
    print(f'{Reffective / 1.49598e13} AU')

    ax.set_xlabel(r'Distance along LOS (au)',fontsize=18)
    ax.set_ylabel(r'Optical thickness, $\tau$',fontsize=18)
    ax.tick_params(axis='both', which='major', labelsize=15)

    fig.tight_layout()
    fig.savefig('figs/optthick_los.pdf', dpi=300, facecolor="white")
    fig.show()






# List some numbers on max-changes to due to smoothing in the star
if list_smoothingchanges == 'y':
    # TODO
    # clean up this and make it work
    # Skip this all?

    # Temperature smoothing

    # Load c5d-temperature in r3d
    Ncells, Nspecies, temperatures = a3d.load_temperature(
        path='../dust_temperature_onestar_190.dat',
        numb_specie=1
    )

    # load starradius
    Mstar,starradius,Lstar,Tstar = a5d.load_star_information(
        savpath='../../exwings_archivedata/co5bold_data/dst28gm06n052/st28gm06n052_190.sav',
        printoutput = 'n'
    )
    AUcm = 1.49598e13 # cm

    # Load grid-distances
    gridcoords = a3d.load_griddistances(
        amrpath='../r3dresults/st28gm06n052_staranddust_1/amr_grid.inp',
        gridpath='../r3dresults/st28gm06n052_staranddust_1/grid_distances.csv',
    )
    griddistances = gridcoords[:,0]

    # Extract position of and strength of spikes (see func "smooth_stellardata")
    smooth_out = 10
    smooth_in = 2
    smooth_tolerance = 1.0

    # search and list "deviation-number" vs cell-number
    deviationlist = []
            
    # Loop over grid cells of each specie (except the outermost cells)
    for ncell in range(smooth_out,Ncells-smooth_out):

        # Check if nn is inside the star (plus a small tolerance factor)
        if griddistances[ncell] < 1.01*starradius:

            # Index of cell in total list (leftover from function)
            nn = ncell

            nindeces = [
                nmedian for nmedian in range(nn-smooth_out,nn+smooth_out+1) if nmedian < (nn-smooth_in) or nmedian > (nn+smooth_in)
            ]
            median_list = []

            for nmedian in nindeces:
                median_list.append(temperatures[nmedian])

            median_temperature = np.median(np.array(median_list))

            # Save deviating cell-number, temperature, deviation and corresponding median T
            if temperatures[nn] > smooth_tolerance * median_temperature:
                deviationlist.append([
                    nn,
                    temperatures[nn],
                    temperatures[nn]/median_temperature,
                    median_temperature
                ])

    # Change to array
    spikearray = np.zeros((len(deviationlist),4))
    for nn,line in enumerate(deviationlist):
        spikearray[nn,0] = int(line[0])
        spikearray[nn,1] = line[1]
        spikearray[nn,2] = line[2]
        spikearray[nn,3] = line[3]


    # Find maximum spike cellnumber, in terms of temperature and not deviation
    maxspike = (spikearray[:,1]).argmax()
    nn_maxspike = int(spikearray[
        maxspike,0
    ])


    # Print coords of max-spike
    print('Max-spike-coords:')
    print(f'    R: {gridcoords[nn_maxspike,0]/AUcm}')
    print(f'    X: {gridcoords[nn_maxspike,1]/AUcm}')
    print(f'    Y: {gridcoords[nn_maxspike,2]/AUcm}')
    print(f'    Z: {gridcoords[nn_maxspike,3]/AUcm}')
    print(f'    Tbefore: {spikearray[maxspike,1]}')
    print(f'    Tafter:  {spikearray[maxspike,3]}')
    print('')


    # Set up cell-index-numbers
    nn_range = 500
    nn_cells = np.linspace(
        nn_maxspike-nn_range+1,
        nn_maxspike+nn_range,
        nn_range*2
    )

    # Plot cell-temperatures in range
    plt.figure(1)
    plt.plot(nn_cells,temperatures[nn_maxspike-nn_range:nn_maxspike+nn_range],'b.',markersize=1)
    # and vs radius
    plt.figure(2)
    plt.plot(
        gridcoords[nn_cells.astype(int),0]/AUcm,
        temperatures[nn_maxspike-nn_range:nn_maxspike+nn_range],
        'b.',markersize=1
    )

    # Plot tmeperature of "spike-cells" before smoothing
    #print('R, X, Y, Z')
    for nn_deviation,nn_spikes in enumerate(spikearray[:,0]):

        nn_spikes = int(nn_spikes)

        for nn in nn_cells:
            if nn == nn_spikes:

                # Mark smoothed cells
                plt.figure(1)
                plt.plot(nn_spikes,temperatures[nn_spikes],'r.',markersize=2)
                plt.figure(2)
                plt.plot(
                    gridcoords[nn_spikes,0]/AUcm,
                    temperatures[nn_spikes],
                    'r.',markersize=2
                )
                # Plot new values to show change
                plt.plot(
                    gridcoords[nn_spikes,0]/AUcm,
                    spikearray[nn_deviation,3],
                    'g.',markersize=2
                )




    # Density-smoothing

    # Load c5d-density in r3d
    Ncells, Nspecies, densities = a3d.load_dustdensity(
        path='../dust_density_onestar_190.inp',
        numb_specie=1
    )

    # load starradius
    Mstar,starradius,Lstar,Tstar = a5d.load_star_information(
        savpath='../../exwings_archivedata/co5bold_data/dst28gm06n052/st28gm06n052_190.sav',
        printoutput = 'n'
    )
    AUcm = 1.49598e13 # cm

    # Load grid-distances
    gridcoords = a3d.load_griddistances(
        amrpath='../r3dresults/st28gm06n052_staranddust_1/amr_grid.inp',
        gridpath='../r3dresults/st28gm06n052_staranddust_1/grid_distances.csv',
    )
    griddistances = gridcoords[:,0]


    # Extract position of and strength of spikes (see func "smooth_stellardata")
    smooth_out = 37
    smooth_in = 36
    smooth_tolerance = 0.4

    # search and list "deviation-number" vs cell-number
    deviationlist = []
            
    # Loop over grid cells of each specie (except the outermost cells)
    for nn in range(smooth_out,Ncells-smooth_out):

        # Check if nn is inside the star (plus a small tolerance factor)
        if griddistances[nn] < 1.01*starradius:

            nindeces = [
                nmedian for nmedian in range(nn-smooth_out,nn+smooth_out+1) if nmedian < (nn-smooth_in) or nmedian > (nn+smooth_in)
            ]
            median_list = []

            for nmedian in nindeces:
                median_list.append(densities[nmedian])

            median_densities = np.median(np.array(median_list))

            # Save deviating cell-number, density, deviation and corresponding median T
            if densities[nn] < smooth_tolerance * median_densities:
                deviationlist.append([
                    nn,
                    densities[nn],
                    densities[nn]/median_densities,
                    median_densities
                ])


    # Change to array
    spikearray = np.zeros((len(deviationlist),4))
    for nn,line in enumerate(deviationlist):
        spikearray[nn,0] = int(line[0])
        spikearray[nn,1] = line[1]
        spikearray[nn,2] = line[2]
        spikearray[nn,3] = line[3]


    # Find maximum spike cellnumber, in terms of temperature and not deviation
    maxspike = (spikearray[:,1]).argmax()
    nn_maxspike = int(spikearray[
        maxspike,0
    ])


    # Print coords of max-spike
    print('Max-spike-coords:')
    print(f'    R: {gridcoords[nn_maxspike,0]/AUcm}')
    print(f'    X: {gridcoords[nn_maxspike,1]/AUcm}')
    print(f'    Y: {gridcoords[nn_maxspike,2]/AUcm}')
    print(f'    Z: {gridcoords[nn_maxspike,3]/AUcm}')
    print(f'    RHObefore: {spikearray[maxspike,1]}')
    print(f'    RHOafter:  {spikearray[maxspike,3]}')
    print('')



    # Set up cell-index-numbers
    nn_range = 1000
    nn_cells = np.linspace(
        nn_maxspike-nn_range+1,
        nn_maxspike+nn_range,
        nn_range*2
    )

    # Plot cell-density in range
    plt.figure(1)
    plt.plot(
        nn_cells,
        densities[nn_maxspike-nn_range:nn_maxspike+nn_range],
        'b.',markersize=1
    )
    # and vs radius
    plt.figure(2)
    plt.plot(
        gridcoords[nn_cells.astype(int),0]/AUcm,
        densities[nn_maxspike-nn_range:nn_maxspike+nn_range],
        'b.',markersize=1
    )

    # Plot tmeperature of "spike-cells" before smoothing
    #print('R, X, Y, Z')
    for nn_deviation,nn_spikes in enumerate(spikearray[:,0]):

        nn_spikes = int(nn_spikes)

        for nn in nn_cells:
            if nn == nn_spikes:

                # Mark smoothed cells
                plt.figure(1)
                plt.plot(
                    nn_spikes,
                    densities[nn_spikes],
                    'r.',markersize=2
                )
                plt.figure(2)
                plt.plot(
                    gridcoords[nn_spikes,0]/AUcm,
                    densities[nn_spikes],
                    'r.',markersize=2
                )
                # Plot new values to show change
                plt.plot(
                    gridcoords[nn_spikes,0]/AUcm,
                    spikearray[nn_deviation,3],
                    'g.',markersize=2
                )



    # Opacity-smoothing

    # Load c5d-opacity in r3d
    opacity = c3d.load_staropacities(
        path='../star_opacities_190.dat',
    )
    Ncells = opacity.size

    # load starradius
    Mstar,starradius,Lstar,Tstar = a5d.load_star_information(
        savpath='../../exwings_archivedata/co5bold_data/dst28gm06n052/st28gm06n052_190.sav',
        printoutput = 'n'
    )
    AUcm = 1.49598e13 # cm

    # Load grid-distances
    gridcoords = a3d.load_griddistances(
        amrpath='../r3dresults/st28gm06n052_staranddust_1/amr_grid.inp',
        gridpath='../r3dresults/st28gm06n052_staranddust_1/grid_distances.csv',
    )
    griddistances = gridcoords[:,0]


    # Extract position of and strength of spikes (see func "smooth_stellardata")
    smooth_out = 7
    smooth_in = 6
    smooth_tolerance_log = 0

    # search and list "deviation-number" vs cell-number
    deviationlist = []
            
    # Loop over grid cells of each specie (except the outermost cells)
    for nn in range(smooth_out,Ncells-smooth_out):

        # Check if nn is inside the star (plus a small tolerance factor)
        if griddistances[nn] < 1.01*starradius:

            nindeces = [
                nmedian for nmedian in range(nn-smooth_out,nn+smooth_out+1) if nmedian < (nn-smooth_in) or nmedian > (nn+smooth_in)
            ]
            median_list = []

            for nmedian in nindeces:
                median_list.append(opacity[nmedian])

            median_opacity = np.median(np.array(median_list))

            # Save deviating cell-number, opacity, deviation and corresponding median T
            if opacity[nn] < 10**-smooth_tolerance_log * median_opacity :
                deviationlist.append([
                    nn,
                    opacity[nn],
                    opacity[nn]/median_opacity,
                    median_opacity
                ])


    # Change to array
    spikearray = np.zeros((len(deviationlist),4))
    for nn,line in enumerate(deviationlist):
        spikearray[nn,0] = int(line[0])
        spikearray[nn,1] = line[1]
        spikearray[nn,2] = line[2]
        spikearray[nn,3] = line[3]

    # Find maximum spike cellnumber, in terms of temperature and not deviation
    maxspike = (spikearray[:,1]).argmax()
    nn_maxspike = int(spikearray[
        maxspike,0
    ])


    # Print coords of max-spike
    print('Max-spike-coords:')
    print(f'    R: {gridcoords[nn_maxspike,0]/AUcm}')
    print(f'    X: {gridcoords[nn_maxspike,1]/AUcm}')
    print(f'    Y: {gridcoords[nn_maxspike,2]/AUcm}')
    print(f'    Z: {gridcoords[nn_maxspike,3]/AUcm}')
    print(f'    KAPPAbefore: {spikearray[maxspike,1]}')
    print(f'    KAPPAafter:  {spikearray[maxspike,3]}')
    print('')


    # Set up cell-index-numbers
    nn_range = 1000
    nn_cells = np.linspace(
        nn_maxspike-nn_range+1,
        nn_maxspike+nn_range,
        nn_range*2
    )

    # Plot cell-opacity in range
    plt.figure(1)
    plt.plot(
        nn_cells,
        opacity[nn_maxspike-nn_range:nn_maxspike+nn_range],
        'b.',markersize=1
    )
    # and vs radius
    plt.figure(2)
    plt.plot(
        gridcoords[nn_cells.astype(int),0]/AUcm,
        opacity[nn_maxspike-nn_range:nn_maxspike+nn_range],
        'b.',markersize=1
    )

    # Plot tmeperature of "spike-cells" before smoothing
    #print('R, X, Y, Z')
    for nn_deviation,nn_spikes in enumerate(spikearray[:,0]):

        nn_spikes = int(nn_spikes)

        for nn in nn_cells:
            if nn == nn_spikes:

                # Mark smoothed cells
                plt.figure(1)
                plt.plot(
                    nn_spikes,
                    opacity[nn_spikes],
                    'r.',markersize=2
                )
                plt.figure(2)
                plt.plot(
                    gridcoords[nn_spikes,0]/AUcm,
                    opacity[nn_spikes],
                    'r.',markersize=2
                )
                # Plot new values to show change
                plt.plot(
                    gridcoords[nn_spikes,0]/AUcm,
                    spikearray[nn_deviation,3],
                    'g.',markersize=2
                )

    print('hej')


# Plot (in model-info part) unmodified pristine data
if plot_2dslices == 'y':
    # 2D slices of input data to r3d

    # Set up paths to data to plot
    gasdensitypath = '../r3dresults/st28gm06n052_staranddust_1/190/dust_density_onestar_190.inp'
    gastemperaturepath = '../r3dresults/st28gm06n052_staranddust_1/190/dust_temperature_onestar_190.dat'
    grainsizespath = '../r3dresults/st28gm06n052_staranddust_1/grain_sizes_190.dat'

    # And paths to grid information
    gridcoordspath = '../r3dresults/st28gm06n052_staranddust_1/grid_distances.csv'
    gridsizepath = '../r3dresults/st28gm06n052_staranddust_1/grid_cellsizes.csv'
    amrpath = '../r3dresults/st28gm06n052_staranddust_1/amr_grid.inp'

    # Load all cellsizes
    gridsizes = a3d.load_cellsizes(gridsizepath,amrpath)

    # Load coordinates of all cells
    gridcoords = a3d.load_griddistances(gridcoordspath,amrpath)
    # gridcoords[:,1] 2] and 3] : x,y, and z distances

    # Load densities
    gasdensity = a3d.load_dustdensity(
        path = gasdensitypath
    )[2]
    # Load temperature
    gastemperatures = a3d.load_temperature(
        path=gastemperaturepath,
    )[2]

    # Load grain sizes in cm
    grainsizes = a3d.load_grainsizes(
        grainsize_path=grainsizespath
    )[0]
    
    # Extract indeces of all close to Z = 0, or actually as close as possible to
    # base cell centre.
    #
    # NOTE: This is harcoded:
    # From amr_grid, middle plane is within
    # 0.5*-6187713212448.75 cm
    # and
    # 0.5*6187713212448.781 cm
    # and
    # N cells along outer edge: 72
    #
    # Or rather, the largest cells are in a plane within 0 and 6187713212448.781
    # Extract cells with centrum between 0 and -0.5*0.618...
    #
    # Chose observed axes
    # One along X-axis, one along Z-axis (eller?)
    los_axes = [1,3]
    # Below or above Z=0?
    basecell = 6187713212448.75
    #basecell = 6187713212448.781

    # Set some image settings based on grid settings
    # NOTE: hardcoded
    mincellsize = np.min(gridsizes)                  # Minimum cell size
    Nside = 72*16    # Number of finest cells along one side -> number of pixels along image side
    gridcourner = 222757675648155.62-0.5*mincellsize # Coords of grid courner in cm
    gridcournerAU = gridcourner/AUcm                 # Same but in AU
    gridsize = gridcourner*2                         # Length of image and grid size

    # Initiate fig-ax-objects
    fig, ax = plt.subplots(
        len(los_axes),3, 
        figsize = (12,9),
    )

    # Loop through the different LOSes
    for naxis,los_axis in enumerate(los_axes):

        # Positions of refinement 0 cells (-1/2)
        NNzeroZ_0 = np.argwhere(
            (gridcoords[:,los_axis] <= -0.99*(0.5)*basecell) & \
            (gridcoords[:,los_axis] >= -1.01*(0.5)*basecell)
        )
        # Positions of refinement 1 cells (-1/4)
        NNzeroZ_1 = np.argwhere(
            (gridcoords[:,los_axis] <= -0.99*(0.5+0.25)*basecell) & \
            (gridcoords[:,los_axis] >= -1.01*(0.5+0.25)*basecell)
        )
        # Positions of refinement 2 cells (-1/4 -1/8)
        NNzeroZ_2 = np.argwhere(
            (gridcoords[:,los_axis] <= -0.99*(0.5+0.25+0.125)*basecell) & \
            (gridcoords[:,los_axis] >= -1.01*(0.5+0.25+0.125)*basecell)
        )
        # Positions of refinement 3 cells (-1/4 -1/8 -1/16)
        NNzeroZ_3 = np.argwhere(
            (gridcoords[:,los_axis] <= -0.99*(0.5+0.25+0.125+0.0625)*basecell) & \
            (gridcoords[:,los_axis] >= -1.01*(0.5+0.25+0.125+0.0625)*basecell)
        )
        # Positions of refinement 4 cells (-1/4 -1/8 -1/16 -1/32)
        NNzeroZ_4 = np.argwhere(
            (gridcoords[:,los_axis] <= -0.99*(0.5+0.25+0.125+0.0625+0.03125)*basecell) & \
            (gridcoords[:,los_axis] >= -1.01*(0.5+0.25+0.125+0.0625+0.03125)*basecell)
        )
        # Merge arrays
        NNzeroZ = np.concatenate((
            NNzeroZ_0,
            NNzeroZ_1,
            NNzeroZ_2,
            NNzeroZ_3,
            NNzeroZ_4,
        ))

        # Create 2D arrays of X-Y-grid
        densities2D = np.zeros((Nside,Nside))
        temperatures2D = np.zeros((Nside,Nside))
        grainsizes2D = np.zeros((Nside,Nside))

        # Rotate "xy"-coordinates depending on los_axis
        #        r3d-z =  c5d-x
        #        r3d-y =  c5d-y
        #        r3d-x = -c5d-z
        if los_axis == 3:
            los_xaxis = 1
            los_yaxis = 2
        if los_axis == 2:
            los_xaxis = 1
            los_yaxis = 3
        if los_axis == 1:
            los_xaxis = 3
            los_yaxis = 2

        # Save Z=0-densities in 2D-arrays
        for nn in NNzeroZ:
            
            xindex = int(np.round((gridcoords[nn,los_xaxis]+gridcourner) / gridsize * Nside))-1
            yindex = int(np.round((gridcoords[nn,los_yaxis]+gridcourner) / gridsize * Nside))-1

            # Refinement levels
            # Radial distances to refinement 1: 0.2475659768178719 - 13.203518763619835 AU
            # Radial distances to refinement 2: 0.4951319536357438 - 9.902639072714877 AU
            # Radial distances to refinement 3: 0.7426979304536157 - 6.601759381809917 AU
            # Radial distances to refinement 4: 0.9902639072714876 - 3.3008796909049587 AU
            #
            # Number of pixels for each cell is
            #
            # 4th refinement
            if gridcoords[nn,0] >= 0.9902639072714876*AUcm and gridcoords[nn,0] <= 3.3008796909049587*AUcm:
                npixels = 1
            
            # 3rd refinement
            if gridcoords[nn,0] > 3.3008796909049587*AUcm and gridcoords[nn,0] <= 6.601759381809917*AUcm:
                npixels = 2
            if gridcoords[nn,0] >= 0.7426979304536157*AUcm and gridcoords[nn,0] < 0.9902639072714876*AUcm:
                npixels = 2

            # 2nd refinement
            if gridcoords[nn,0] > 6.601759381809917*AUcm and gridcoords[nn,0] <= 9.902639072714877*AUcm:
                npixels = 4
            if gridcoords[nn,0] >= 0.4951319536357438*AUcm and gridcoords[nn,0] < 0.7426979304536157*AUcm:
                npixels = 4

            # 1st refinement
            if gridcoords[nn,0] > 9.902639072714877*AUcm and gridcoords[nn,0] <= 13.203518763619835*AUcm:
                npixels = 8
            if gridcoords[nn,0] >= 0.2475659768178719*AUcm and gridcoords[nn,0] < 0.4951319536357438*AUcm:
                npixels = 8

            # Base cells
            if gridcoords[nn,0] < 0.2475659768178719*AUcm or gridcoords[nn,0] > 13.203518763619835*AUcm:
                npixels = 16

            # Fill cells
            # Log of densities
            # temperatures
            # grain sies in um
            densities2D[xindex-npixels:xindex+npixels,yindex-npixels:yindex+npixels] = np.log10(gasdensity[nn])
            temperatures2D[xindex-npixels:xindex+npixels,yindex-npixels:yindex+npixels] = np.log10(gastemperatures[nn])
            grainsizes2D[xindex-npixels:xindex+npixels,yindex-npixels:yindex+npixels] = grainsizes[nn]*1e4

        # Rotate image by 90 deg and set originparameter
        densities2D = ndimage.rotate(densities2D, 90)
        temperatures2D = ndimage.rotate(temperatures2D, 90)
        grainsizes2D = ndimage.rotate(grainsizes2D, 90)
        originparam = 'lower'

        # Switch xy-axis on second row
        if naxis == 1:
            # Rotate 90 deg again
            densities2D = ndimage.rotate(densities2D, -90)
            temperatures2D = ndimage.rotate(temperatures2D, -90)
            grainsizes2D = ndimage.rotate(grainsizes2D, -90)
            originparam = 'upper'


        # Set new axis ranges based on all the rotations and reversed axis.
        axisplot  = [
            -gridcournerAU,
            gridcournerAU,
            gridcournerAU,
            -gridcournerAU
        ]
        # Initialise labels and colour bars
        imbar = []
        colourbarlabels = [
            r'$\log ($ Gas density [g\,cm$^{-3}$] $)$',
            r'$\log ($ Gas temperature [K] $)$',
            r'Grain size ($\mu$m)'
        ]
        # Plot images and save colour bar info
        imbar.append(ax[naxis,0].imshow(
            densities2D,
            origin=originparam, extent=axisplot, 
            cmap=plt.get_cmap('pink'),
            vmin=-17,
            vmax=-6
        ))
        imbar.append(ax[naxis,1].imshow(
            temperatures2D,
            origin=originparam, extent=axisplot, 
            cmap=plt.get_cmap('hot'),
            vmin=2.5,
            vmax=4.7
        ))
        imbar.append(ax[naxis,2].imshow(
            grainsizes2D,
            origin=originparam, extent=axisplot, 
            cmap=plt.get_cmap('bone'),
        ))
        for nn in range(3):
            # Rotate axis to be consistent with freytag2023 (on first row)
            # and so that second row has postive in correct direction
            ax[naxis,nn].invert_xaxis()
            ax[naxis,nn].invert_yaxis()
            ax[naxis,nn].invert_xaxis()
            # and change x-lim and ylim to skip emtpy edges
            ax[naxis,nn].set_xlim(-gridcournerAU+basecell/AUcm,gridcournerAU)
            ax[naxis,nn].set_ylim(-gridcournerAU+basecell/AUcm,gridcournerAU)

            # Set ticksettings    
            ax[naxis,nn].tick_params(axis='both', which='major', labelsize=15)

    # Set colour bars, only top row
    for nn in range(3):
        divider = make_axes_locatable(ax[0][nn])
        cax = divider.append_axes(
            'top', 
            size='4%', 
            pad=0.7,
        )
        cb0 = plt.colorbar(
            imbar[nn], 
            cax=cax, 
            orientation = 'horizontal', 
        )
        cb0.set_label(
            label = colourbarlabels[nn], fontsize=15
        )
        cb0.ax.tick_params(labelsize=15)


    # Set ylabel texts, both rows
    ax[0,0].set_ylabel('Y (au)',fontsize=18)
    ax[1,0].set_ylabel('Y (au)',fontsize=18)
    # Set xlabel text and tick params for all
    for nn in range(3):
        ax[0,nn].set_xlabel('X (au)',fontsize=18) # only bottom row
        ax[1,nn].set_xlabel('Z (au)',fontsize=18) # only bottom row

    fig.tight_layout()
    #Save figure
    fig.savefig(f'figs/co5bold_2dslices.pdf', dpi=300, facecolor="white")

# ----------------------------------------------------------------
#
# Plot Figure with Mass-vs-grain size histogram
#
# Warning: slow because of loading lots dust densities


# Plots 
if plot_grainsizehist == 'y':

    # Initiate fig-axes-object
    fig,ax = plt.figure(
        figsize=(6, 4)
    ), plt.axes()


    # Load binned grain sizes per cell
    grain_sizes,Ncells = a3d.load_grainsizes(
        grainsize_path = '../r3dresults/st28gm06n052_staranddust_1/grain_sizes_binned_190.dat'
    )

    # And save grain sizes in um
    grainsize_bins = np.unique(grain_sizes[np.where(grain_sizes > 0)[0]])*1e4

    # And number of grain sizes, ie species
    Nspecies = np.size(grainsize_bins)

    # And create array for masses per grain size bin
    dustmass_bins = np.zeros(Nspecies)

    # Compute volume per cell - load cell sizes
    cell_sizes = a3d.load_cellsizes(
        sizepath = '../r3dresults/st28gm06n052_staranddust_1/grid_cellsizes.csv',
        amrpath = '../r3dresults/st28gm06n052_staranddust_1/amr_grid.inp'
    )

    # Cell volumes
    cell_volumes = cell_sizes**3

    # Skip header, star, and load all dust data
    # (much faster with a custom density loader since my other density loader
    # opens and closes the file too much.)
    with open(f'../r3dresults/st28gm06n052_staranddust_1/190/dust_density.inp') as f:
        for nn,line in enumerate(f.readlines()):
            if nn > 3+Ncells:

                # Only look at cells with dust:
                data = float(line)

                if data > 0:

                    # Then we are at this species number
                    nspecie = int((nn-3-Ncells)/Ncells)

                    # And this cell
                    ncell = (nn-3-Ncells) - nspecie*Ncells

                    # Add the total dust mass of this cell to the size bin
                    dustmass_bins[nspecie] += cell_volumes[ncell]*data
                    # This is then the dust mass in grams per grain size

    # Plot figures as step-plot with symbols in the moddel
    ax.step(grainsize_bins,dustmass_bins,'b-o',where='mid')

    print(f'Finished figure object')
    
    # Tight layout for nicer image file
    fig.tight_layout()

    ax.set_ylabel(r'Dust mass (g)',fontsize=18)    
    ax.set_xlabel(r'Grain size ($\mu$m)',fontsize=18)
    ax.tick_params(axis='both', which='major', labelsize=15)

    fig.tight_layout()
    fig.savefig('figs/grainsize_hist.pdf', dpi=300, facecolor="white")

    fig.show()


# Plot figure with grain sizes vs R
if plot_grainsizeradius == 'y':

    grainsizes = '../../exwings_archivedata/r3dresults/st28gm06n052_generaldata/grain_sizes_190.dat'
    savpath = '../../exwings_archivedata/co5bold_data/dst28gm06n052/st28gm06n052_190.sav'
    grid_path = '../r3dresults/st28gm06n052_staranddust_1/grid_distances.csv'
    amr_path = '../r3dresults/st28gm06n052_staranddust_1/amr_grid.inp'
    phase = '190'
    fig, ax = plt.figure(figsize=(6,4)), plt.axes()



    # Load grid info
    # Load radial coordinates of R3D-cells and change to AU
    # Load half-size of cube to show where shells are partially outside the cube
    cellcoords = a3d.load_griddistances(
        gridpath=grid_path,
        amrpath=amr_path
    )
    radii = cellcoords[:,0]/AUcm

    # Load star's radius here
    Mstar,Rstar,Lstar,Tstar = a5d.load_star_information(
        savpath = savpath,
        printoutput='n'
    )
    Rstar /= AUcm


    # Fill plots
    Nbins = 100


    # Load grainsizes
    grain_sizes,Ncells = a3d.load_grainsizes(
        grainsize_path = grainsizes
    )
    # Change to um
    grain_sizes *= 1e4


    # N radial bins, set up radial grid
    radial_bins = np.linspace(0,radii.max(),Nbins+1)
    radial_range = np.zeros(Nbins)
    grainsize_bins = np.zeros(Nbins)
    grainsize_max = np.zeros(Nbins)
    grainsize_min = np.zeros(Nbins)
    grainsize_std = np.zeros(Nbins)


    # Loop through bins and save binned data
    for nn in range(Nbins):
        ncells = np.where((radii >= radial_bins[nn]) & (radii < radial_bins[nn+1]))[0]
        grainsize_bins[nn] = grain_sizes[ncells].mean()
        grainsize_max[nn] = grain_sizes[ncells].max()
        grainsize_min[nn] = grain_sizes[ncells].min()
        grainsize_std[nn] = grain_sizes[ncells].std()
        radial_range[nn] = radial_bins[nn] + 0.5 * radial_bins[1]


    
    # Plot max-min-range of grain sizes
    ax.fill_between(
        radial_range,
        grainsize_min,
        grainsize_max,
        color='g',
        alpha=0.2
    )

    # Plot std-range of grain sizes
    ax.fill_between(
        radial_range,
        grainsize_bins-grainsize_std,
        grainsize_bins+grainsize_std,
        color='g',
        alpha=0.4
    )

    # Plot average value
    ax.plot(radial_range,grainsize_bins,'k')

    # Add line for stellar surface and edge of cube
    #ax.plot([Rstar,Rstar],[0,grainsize_max.max()+0.05],'r:',linewidth=1)
    ax.plot([cubesize,cubesize],[0,grainsize_max.max()+0.05],'k:',linewidth=1)

    
    # Plot settings
    ax.set_ylabel(r'Grain sizes ($\mu$m)',fontsize=18)
    ax.set_xlabel(r'Distance (au)',fontsize=18)
    ax.set(
        xlim=(Rstar,radial_range.max()+0.5),
        ylim=(0,grainsize_max.max()+0.05)
    )
    ax.tick_params(axis='both', which='major', labelsize=15)

    print('Grain size fig object is finished, show with fig.show()')

    fig.tight_layout()
    fig.savefig('figs/grainsizes_radius.pdf', dpi=300, facecolor="white")
    fig.show()


# ----------------------------------------------------------------
#
# Plot Figure absorption and scattering, and angles

if plot_absscat == 'y':

    fig,ax = a3d.plot_allkappa(
        path='../r3dresults/opacities_mie_st28gm06n052/'
    )

    ax[0].set_ylabel(r'$\kappa_{\rm abs}$ (cm$^2$/g)',fontsize=18)    
    ax[1].set_ylabel(r'$\kappa_{\rm scat}$ (cm$^2$/g)',fontsize=18)
    ax[2].set_ylabel(r'$\kappa_{\rm abs} + \kappa_{\rm scat}$ (cm$^2$/g)',fontsize=18)
    ax[3].set_ylabel(r'$\left< \cos \theta \right>$',fontsize=18)

    ax[3].set_xlabel(r'Wavelength ($\mu$m)',fontsize=18)

    ax[0].tick_params(axis='both', which='major', labelsize=15)
    ax[1].tick_params(axis='both', which='major', labelsize=15)
    ax[2].tick_params(axis='both', which='major', labelsize=15)
    ax[3].tick_params(axis='both', which='major', labelsize=15)

    fig.tight_layout()
    fig.savefig('figs/abs_scat_angle.pdf', dpi=300, facecolor="white")

    fig.show()

# ----------------------------------------------------------------
#
# Plot Figure with radial dependence of temperatures


# Plot comparisons of temperature between co5bold and radmc3d
if plot_temperaturecompare == 'y':

    # The lowest altitude for dust formation at the included time snapshot corresponds 
    # to 3.4\,au, or 2.1\,$R_\star$, from the grid centre.
    # SO I start these plots at 3au.

    # Chose phase-designation
    phase = 190

    # Set up subplots for figure with only temperature (ratio later in another fig)
    fig,ax = plt.subplots(
        2,1,
        figsize=(6,8)
    )
    # Initiate figobjects for ratio-plot
    figratio, axratio = plt.figure(
        figsize=(6,4)
    ), plt.axes()

    # Set legends
    #legendlist = [
    #    r'$T($CO5BOLD$)$',
    #    r'$T($point$)$',
    #    r'$T($theory$)$'
    #]
    legendlist = [
        r'$T_{\rm CO5BOLD}$',
        r'$T_{\rm point}$',
        r'$T_{\rm theory}$'
    ]


    # Load star's radius here
    Mstar,Rstar,Lstar,Tstar = a5d.load_star_information(
        savpath = f'../../exwings_archivedata/co5bold_data/dst28gm06n052/st28gm06n052_{phase}.sav',
        printoutput = 'n'
    )
    Rstar /= AUcm


    # Load cobold-T and create subplot
    T_c5d, Tstd_c5d, Tmax_c5d, Tmin_c5d, radial_range = a3d.plot_temperaturebins_radius(
        temperature_path=f'../r3dresults/st28gm06n052_staranddust_1/{phase}/dust_temperature.dat',
        grid_path='../r3dresults/st28gm06n052_staranddust_1/grid_distances.csv',
        amr_path='../r3dresults/st28gm06n052_staranddust_1/amr_grid.inp',
        numb_specie = 1,
        ax=ax[0]
    )
    ax[0].set_ylabel(r'$T_{\rm CO5BOLD}$, K',fontsize=18)
    ax[0].set_xlabel(r'')
    ax[0].set(xlim=(3,26), ylim=(0,3000))
    ax[0].tick_params(axis='both', which='major', labelsize=15)

    # Load R3d-temperatures and bin each grain size's tmeperature
    # and take averages, max-min, stds of these
    # Re-use radial range from c5d-plot since the grid is identical

    # Load coordintaes of R3D-cells and change to AU
    # extract size of grid cube (round up because this is the centre of cell coord)
    # and create binned radius-array
    cellcoords = a3d.load_griddistances(
        gridpath='../r3dresults/st28gm06n052_pointtemperature/grid_distances.csv',
        amrpath='../r3dresults/st28gm06n052_pointtemperature/amr_grid.inp'
    )
    radii = cellcoords[:,0]/AUcm
    Nbins = 100
    radial_bins = np.linspace(0,radii.max(),Nbins+1)


    # Create temporary temperature arrays
    Nspecies = 10
    temperatures_avr = np.zeros((Nbins,Nspecies))
    temperatures_max = np.zeros((Nbins,Nspecies))
    temperatures_min = np.zeros((Nbins,Nspecies))
    temperatures_std = np.zeros((Nbins,Nspecies))

    # Bin each grain size temperature
    for nspecie in range(Nspecies):

        # Load each specie's R3D-temperature
        Ncells, Ntemp, temperatures = a3d.load_temperature(
            path=f'../r3dresults/st28gm06n052_pointtemperature/{phase}/dust_temperature.dat',
            numb_specie=nspecie+1
        )

        # Bin the temperatures and save in a Nbins*Nspecies arrays
        for nn in range(Nbins):
            ncells = np.where((radii >= radial_bins[nn]) & (radii < radial_bins[nn+1]))[0]
            temperatures_avr[nn,nspecie] = temperatures[ncells].mean()
            temperatures_max[nn,nspecie] = temperatures[ncells].max()
            temperatures_min[nn,nspecie] = temperatures[ncells].min()
            temperatures_std[nn,nspecie] = temperatures[ncells].std()


    # Save average/max-min/std's of each specie 
    # This saves the maximum standard deviation of eahc specie at each radial bin
    Tr3d_avr = np.zeros(Nbins)
    Tr3d_max = np.zeros(Nbins)
    Tr3d_min = np.zeros(Nbins)
    Tr3d_std = np.zeros(Nbins)

    for nn in range(Nbins):
        Tr3d_avr[nn] = temperatures_avr[nn,:].mean()
        Tr3d_max[nn] = temperatures_max[nn,:].max()
        Tr3d_min[nn] = temperatures_min[nn,:].min()
        Tr3d_std[nn] = temperatures_std[nn,:].max()


    ax[1].plot(radial_range,Tr3d_avr,'k')

    ax[1].fill_between(
        radial_range,
        Tr3d_min,
        Tr3d_max,
        color='r',
        alpha=0.2
    )

    ax[1].fill_between(
        radial_range,
        Tr3d_avr-Tr3d_std,
        Tr3d_avr+Tr3d_std,
        color='r',
        alpha=0.4
    )

    ax[1].set_ylabel(r'$T_{\rm point}$, K',fontsize=18)
    ax[1].set_xlabel(r'Distance (au)',fontsize=18)
    ax[1].set(xlim=(3,26), ylim=(0,3000))
    ax[1].tick_params(axis='both', which='major', labelsize=15)


    # A plot with Tc5d or Tr3d divided by T_theory
    Teff = 2800   # Table1
    pindex = -0.9 # Bladh 2012
    temperature_theory = Teff * (Rstar/(2*radial_range))**(2/(4+pindex))

    # A vertical line at limit of grid for all subplots
    for nn in range(2):
        # Grid cube
        ax[nn].plot([cubesize,cubesize],[0,4100],'k:',linewidth=1)
        # And stellar radius
        ax[nn].plot([Rstar,Rstar],[0,4100],'r:',linewidth=1)



    # Plot separate figure with temperature-ratios
    axratio.plot(radial_range,T_c5d/temperature_theory,'b',linewidth=2,label = legendlist[0])
    axratio.plot(radial_range,Tr3d_avr/temperature_theory,'r',linewidth=2,label = legendlist[1])
    axratio.plot([Rstar,radial_range[-1]],[1,1],'g:', label = legendlist[2])
    axratio.legend(fontsize=13)

    axratio.fill_between(
        radial_range,
        (T_c5d-Tstd_c5d)/temperature_theory,
        (T_c5d+Tstd_c5d)/temperature_theory,
        color='b',
        alpha=0.4
    )

    axratio.fill_between(
        radial_range,
        (Tr3d_avr-Tr3d_std)/temperature_theory,
        (Tr3d_avr+Tr3d_std)/temperature_theory,
        color='r',
        alpha=0.4
    )

    axratio.set_ylabel(r'$T_{\rm simulated}$ / $T_{\rm theory}$',fontsize=18)
    axratio.set_xlabel(r'Distance (au)',fontsize=18)
    axratio.set(
        ylim=(0.5,1.7),
        xlim=(3,26)
    )
    axratio.tick_params(axis='both', which='major', labelsize=15)

    # A vertical line at limit of grid
    # Grid cube
    axratio.plot([cubesize,cubesize],[0,3],'k:',linewidth=1)
    # And stellar radius
    axratio.plot([Rstar,Rstar],[0,3],'r:',linewidth=1)


    # Show all plots
    fig.tight_layout()
    figratio.tight_layout()
    fig.show()
    figratio.show()

    #Save figure
    fig.savefig(f'figs/temperatures_{phase}.pdf', dpi=300, facecolor="white")
    figratio.savefig(f'figs/temperatureratio_{phase}.pdf', dpi=300, facecolor="white")

    # Print some temperatures and ratios for a table in the paper
    # at 3 au, 6au, 9 au

    three_au = np.where(radial_range > 3)[0][0] - 1
    six_au = np.where(radial_range > 6)[0][0] - 1
    nine_au = np.where(radial_range > 9)[0][0] - 1
    fifteen_au = np.where(radial_range > 15)[0][0] - 1

    print(f' 3 au. Tc5d: {T_c5d[three_au]}   Tr3d: {Tr3d_avr[three_au]}   Tc5d/Tr3d: {T_c5d[three_au]/Tr3d_avr[three_au]}   Tc5d/Ttheory: {T_c5d[three_au]/temperature_theory[three_au]}   Tr3d/Ttheory: {Tr3d_avr[three_au]/temperature_theory[three_au]}')
    print(f' 6 au. Tc5d: {T_c5d[six_au]}   Tr3d: {Tr3d_avr[six_au]}   Tc5d/Tr3d: {T_c5d[six_au]/Tr3d_avr[six_au]}   Tc5d/Ttheory: {T_c5d[six_au]/temperature_theory[six_au]}   Tr3d/Ttheory: {Tr3d_avr[six_au]/temperature_theory[six_au]}')
    print(f' 9 au. Tc5d: {T_c5d[nine_au]}   Tr3d: {Tr3d_avr[nine_au]}   Tc5d/Tr3d: {T_c5d[nine_au]/Tr3d_avr[nine_au]}   Tc5d/Ttheory: {T_c5d[nine_au]/temperature_theory[nine_au]}   Tr3d/Ttheory: {Tr3d_avr[nine_au]/temperature_theory[nine_au]}')
    print(f'15 au. Tc5d: {T_c5d[fifteen_au]}   Tr3d: {Tr3d_avr[fifteen_au]}   Tc5d/Tr3d: {T_c5d[fifteen_au]/Tr3d_avr[fifteen_au]}   Tc5d/Ttheory: {T_c5d[fifteen_au]/temperature_theory[fifteen_au]}   Tr3d/Ttheory: {Tr3d_avr[fifteen_au]/temperature_theory[fifteen_au]}')





# ----------------------------------------------------------------
#
# Plot figs for co5bold-results

if plot_coboldsed == 'y':


    # Initiate figure
    fig, ax = plt.figure(figsize=(6,5)), plt.axes()

    # Set path info for SEDs
    path = '../r3dresults/st28gm06n052_staranddust_nospikes/'
    phase = '190'
    spectra = [
        '/spectrum_i000_phi000.out',
        '/spectrum_i090_phi000.out',
        '/spectrum_i090_phi090.out',
        '/spectrum_i090_phi270.out',
        '/spectrum_i180_phi000.out',
        '/spectrum_i270_phi000.out'
    ]
    legendlist = [
        r'0-0; \hspace{5mm} 3759',
        r'90-0; \hspace{3mm} 4597',
        r'90-90; \hspace{1mm} 5471',
        r'90-270; 7821',
        r'180-0; \hspace{1mm} 8559',
        r'270-0; \hspace{1mm} 6631',
    ]

    for nangle,spectrum in enumerate(spectra):

        wavelength,sed = a3d.load_spectrum(
            path = path+phase+spectrum
        )
        
        # Print all sed-luminosities
        luminosity = a3d.compute_luminosity(
            wavelengths=wavelength,
            spectrum=sed
        )
        print(f'{path}{spectrum}: {luminosity/Lsol}')
            
        ax.plot(wavelength,sed,label = legendlist[nangle])

    ax.set(xscale='log',yscale='log')
    ax.set_xlim(5e-1,6e1)
    ax.set_ylim(1e6,1.3e8)
    ax.tick_params(axis='both', which='major', labelsize=15)
    ax.set_xlabel(r'Wavelength ($\mu$m)', fontsize=18)    
        
    # Ylabel, titles and legend with angles
    ax.set_ylabel(r'$F_\nu$, Jy at 1 pc', fontsize=18)
    ax.set_title(rf'CO5BOLD', fontsize=16)
    ax.legend(title=r'$i$-$\phi$;   Int. flux dens. ($L_\odot$)', fontsize=13)

    fig.tight_layout()
    #fig.show()
   
    #Save figure
    fig.savefig(f'figs/seds_c5d.pdf', dpi=300, facecolor="white")


# Plot co5bold-images
if plot_images_examples == 'y':
    # Plot 2 images at each walanvength, one phase, lin and log

    # Chose wavelengths here
    wavelengths = ['01','10']

    # Suggested flux density limits for plots
    # lin and log
    fluxlimitslin = [
        [0,4],
        [0,1.4]
    ]
    fluxlimitslog = [
        [-9,7],
        [1,7]
    ]
    # Old numbers
    #  1um: 0 -> 3   (lin), -10 -> 6 (log)
    # 10um: 0 -> 1.5 (lin),  1  -> 6 (log)

    # Initialise fig
    fig, ax = plt.subplots(
        2,len(wavelengths),
        figsize = (10,8),
    )

    # Chose distnace
    distance = 1

    # Loop through wavelengths and plot images
    for nwave,wavelengthum in enumerate(wavelengths):

        # Image file name:
        image = f'../r3dresults/st28gm06n052_staranddust_1/190/image_i000_phi000_{wavelengthum}um.out'

        # Extract path and imagename from image
        imagestrings = re.split('/', image)
        path = f'{imagestrings[0]}/{imagestrings[1]}/'
        modelname = imagestrings[2]
        phase = imagestrings[3]
        imagefilename = imagestrings[4]

        # extract inclination and wavelength
        imagestrings = re.split('_', imagefilename)
        incl = imagestrings[1][1:]
        phi = imagestrings[2][3:]

        # Load data
        image2d,image2dlog,flux,axisplot = a3d.load_images(
            path=f'{path}{modelname}/{phase}',
            image=imagefilename,
            distance=distance
        )
        # Change to MJy per asec2
        image2d = image2d*1e-6

        # Plot linear scale image
        imlin = ax[0][nwave].imshow(
            image2d, 
            origin='lower', extent=axisplot, 
            cmap=plt.get_cmap('hot'),
            vmin=fluxlimitslin[nwave][0],
            vmax=fluxlimitslin[nwave][1]
        )
        ax[0][nwave].set_title(rf'{int(wavelengthum)}\,$\mu$m', fontsize=15)
        ax[0][nwave].tick_params(axis='both', which='major', labelsize=15)

        # Plot log-scale image
        imlog = ax[1][nwave].imshow(
            image2dlog, 
            origin='lower', extent=axisplot, 
            cmap=plt.get_cmap('hot'),
            vmin=fluxlimitslog[nwave][0],
            vmax=fluxlimitslog[nwave][1]
        )
        ax[1][nwave].set_xlabel('Offset (au)',fontsize=18)
        ax[1][nwave].tick_params(axis='both', which='major', labelsize=15)



        # Set colour bar settings and label, first row (linear)
        divider = make_axes_locatable(ax[0][nwave])
        cax = divider.append_axes(
            'right', 
            size='5%', 
            pad=0.05
        )
        cb0 = plt.colorbar(
            imlin, 
            cax=cax, 
            orientation = 'vertical',
        )
        cb0.set_label(
            label = rf'$F_\nu$(MJy/asec$^2$) at {distance} pc',fontsize= 15
        )
        cb0.ax.tick_params(labelsize=15)

        # Set colour bar settings and label, second row (logarithmic)
        divider = make_axes_locatable(ax[1][nwave])
        cax = divider.append_axes(
            'right', 
            size='5%', 
            pad=0.05
        )
        cb0 = plt.colorbar(
            imlog, 
            cax=cax, 
            orientation = 'vertical'
        )
        cb0.set_label(
            label = rf'$\log F$(Jy/asec$^2$) at {distance} pc',fontsize= 15
        )
        cb0.ax.tick_params(labelsize=15)


    ax[0][0].set_ylabel('Offset (au)',fontsize=18)
    ax[1][0].set_ylabel('Offset (au)',fontsize=18)

    fig.tight_layout()
    #fig.show()

    #Save figure
    fig.savefig(f'figs/images_{wavelengthum}umexamples.pdf', dpi=300, facecolor="white")





# ------------------------------------------------------------------------------
# Plot symmetric figs

# Plot Darwin SED and image
if plot_darwin_imagesed == 'y':

    # Initiate figure
    fig,ax = plt.subplots(2,1, figsize = (6, 9))

    # Set path info for SEDs
    path = '../r3dresults/st28gm06n052_darwinsource/'
    phase = '190'
    spectra = [
        '/spectrum_i000_phi000.out',
        '/spectrum_i090_phi000.out',
        '/spectrum_i090_phi090.out',
        '/spectrum_i090_phi270.out',
        '/spectrum_i180_phi000.out',
        '/spectrum_i270_phi000.out'
    ]
    legendlist = [
        '0-0',
        '90-0',
        '90-90',
        '90-270',
        '180-0',
        '270-0',
    ]
    # Set path and info for 10um image
    image = '../r3dresults/st28gm06n052_darwinsource/190/image_i000_phi000_10um.out'
    distance = 1 # 1 pc


    # Plot SED
    for nangle,spectrum in enumerate(spectra):

        wavelength,sed = a3d.load_spectrum(
            path = path+phase+spectrum
        )
        
        # Print all sed-luminosities
        luminosity = a3d.compute_luminosity(
            wavelengths=wavelength,
            spectrum=sed
        )
        print(f'{path}{spectrum}: {luminosity/Lsol}')
            
        ax[1].plot(wavelength,sed,label = legendlist[nangle])

    # Set SED axis settings
    ax[1].set(xscale='log',yscale='log')
    ax[1].set_xlim(5e-1,6e1)
    ax[1].set_ylim(1e6,1.3e8)
    ax[1].tick_params(axis='both', which='major', labelsize=15)
    ax[1].set_xlabel(r'Wavelength ($\mu$m)', fontsize=18)    
        
    # Ylabel, titles and legend with angles
    ax[1].set_ylabel(r'$F_\nu$, Jy at 1 pc', fontsize=18)
    ax[1].legend(title=r'$i$-$\phi$', fontsize=13)


    # Plot image

    # Extract path and imagename from image
    imagestrings = re.split('/', image)
    path = f'{imagestrings[0]}/{imagestrings[1]}/'
    modelname = imagestrings[2]
    phase = imagestrings[3]
    imagefilename = imagestrings[4]

    # extract inclination and wavelength
    imagestrings = re.split('_', imagefilename)
    incl = imagestrings[1][1:]
    phi = imagestrings[2][3:]
    wavelengthum = imagestrings[3][:-6]

    # Load data
    image2d,image2dlog,flux,axisplot = a3d.load_images(
        path=f'{path}{modelname}/{phase}',
        image=imagefilename,
        distance=distance
    )
    # Change to MJy per asec2
    image2d = image2d*1e-6

    # Plot image
    imlin = ax[0].imshow(
        image2d, 
        origin='lower', extent=axisplot, 
        cmap=plt.get_cmap('hot'),
        vmin=0,vmax=1.4
    )

    # Set axis settings
    ax[0].tick_params(axis='both', which='major', labelsize=15)
    ax[0].set_xlabel('Offset (au)', fontsize=18)
    ax[0].set_ylabel('Offset (au)', fontsize=18)

    # Set title
    ax[0].set_title(rf'DARWIN', fontsize=16)


    # Set colour bar settings and label
    divider = make_axes_locatable(ax[0])
    cax = divider.append_axes(
        'right', 
        size='4%', 
        pad=0.05,
    )
    cb0 = plt.colorbar(
        imlin, 
        cax=cax, 
        orientation = 'vertical', 
    )
    cb0.set_label(
        label = rf'$F_\nu$(MJy/asec$^2$) at {wavelengthum} $\mu$m \& {distance} pc', fontsize=15
    )
    cb0.ax.tick_params(labelsize=15)
    fig.tight_layout()
    
    #Save figure
    fig.savefig(f'figs/sedimage_darwin.pdf', dpi=300, facecolor="white")



# To plot comparison-SED-plot between sphere and c5d-star
if plot_darwin_comparesed == 'y':

    # Initiate figure
    fig, ax = plt.figure(figsize=(6,5)), plt.axes()

    # Folders with SEDs
    pathc5d = '../r3dresults/st28gm06n052_staranddust_nospikes/190/'
    pathdarwin = '../r3dresults/st28gm06n052_darwinsource/190/'

    # List of angles and legend list
    angles = [
        ['000','000'],
        ['090','000'],
        ['090','090'],
        ['090','270'],
        ['180','000'],
        ['270','000'],
    ]
    legendlist = [
        r'0-0',
        r'90-0',
        r'90-90',
        r'90-270',
        r'180-0',
        r'270-0',
    ]
    # Load all and subtract SEDs
    for nangle,angle in enumerate(angles):

        wavelengths, sedc5d = a3d.load_spectrum(
            path = f'{pathc5d}spectrum_i{angle[0]}_phi{angle[1]}.out'
        )
        wavelengths, seddrw = a3d.load_spectrum(
            path = f'{pathdarwin}spectrum_i{angle[0]}_phi{angle[1]}.out'
        )
        # Check that SEDs are equal length before subtracting
        if len(sedc5d) == len(seddrw):
            seddiff = np.zeros(len(wavelengths))
            
            for nn in range(len(wavelengths)):
                seddiff[nn] = seddrw[nn] - sedc5d[nn]

        else:
            raise ValueError('ERROR: SEDs are not of equal length')

        # Plot
        ax.plot(wavelengths, seddiff, label = legendlist[nangle])

    # Plot zero-line
    ax.plot([wavelengths[0],wavelengths[-1]],[0,0],'k--')

    # Plot settings
    ax.set(xscale='log')
    ax.set_xlim(5e-1,6e1)
    ax.set_xlabel(r'Wavelength ($\mu$m)', fontsize=18)
    ax.set_ylabel(r'$F_{\rm DARWIN} - F_{\rm CO5BOLD}$ (Jy)', fontsize=18)
    ax.tick_params(axis='both', which='major', labelsize=15)
    ax.legend(title=r'$i$-$\phi$', fontsize=13)
    fig.tight_layout()
   
    #Save figure
    fig.savefig(f'figs/compareseds_c5ddarwin.pdf', dpi=300, facecolor="white")


# Plot point source SED and image
if plot_point_imagesed == 'y':

    # Initiate figure
    fig,ax = plt.subplots(2,1, figsize = (6, 9))

    # Set path info for SEDs
    path = '../r3dresults/st28gm06n052_pointtemperature/'
    phase = '190'
    spectra = [
        '/spectrum_i000_phi000.out',
        '/spectrum_i090_phi000.out',
        '/spectrum_i090_phi090.out',
        '/spectrum_i090_phi270.out',
        '/spectrum_i180_phi000.out',
        '/spectrum_i270_phi000.out'
    ]
    legendlist = [
        '0, 0',
        '90, 0',
        '90, 90',
        '90, 270',
        '180, 0',
        '270, 0',
    ]
    # Set path and info for 10um image
    image = '../r3dresults/st28gm06n052_pointtemperature/190/image_i000_phi000_10um.out'
    distance = 1 # 1 pc


    # Plot SED
    for nangle,spectrum in enumerate(spectra):

        wavelength,sed = a3d.load_spectrum(
            path = path+phase+spectrum
        )
        
        # Print all sed-luminosities
        luminosity = a3d.compute_luminosity(
            wavelengths=wavelength,
            spectrum=sed
        )
        print(f'{path}{spectrum}: {luminosity/Lsol}')
            
        ax[1].plot(wavelength,sed,label = legendlist[nangle])

    # Set SED axis settings
    ax[1].set(xscale='log',yscale='log')
    ax[1].set_xlim(5e-1,6e1)
    ax[1].set_ylim(1e6,1.3e8)
    ax[1].tick_params(axis='both', which='major', labelsize=15)
    ax[1].set_xlabel(r'Wavelength ($\mu$m)', fontsize=18)    
        
    # Ylabel, titles and legend with angles
    ax[1].set_ylabel(r'$F_\nu$, Jy at 1 pc', fontsize=18)
    ax[1].legend(title=r'$i$, $\phi$', fontsize=13)


    # Plot image

    # Extract path and imagename from image
    imagestrings = re.split('/', image)
    path = f'{imagestrings[0]}/{imagestrings[1]}/'
    modelname = imagestrings[2]
    phase = imagestrings[3]
    imagefilename = imagestrings[4]

    # extract inclination and wavelength
    imagestrings = re.split('_', imagefilename)
    incl = imagestrings[1][1:]
    phi = imagestrings[2][3:]
    wavelengthum = imagestrings[3][:-6]

    # Load data
    image2d,image2dlog,flux,axisplot = a3d.load_images(
        path=f'{path}{modelname}/{phase}',
        image=imagefilename,
        distance=distance
    )
    # Change to MJy per asec2
    image2d = image2d*1e-6

    # Plot image
    imlin = ax[0].imshow(
        image2d, 
        origin='lower', extent=axisplot, 
        cmap=plt.get_cmap('hot'),
        vmin=0,vmax=0.7
    )

    # Set axis settings
    ax[0].tick_params(axis='both', which='major', labelsize=15)
    ax[0].set_xlabel('Offset (au)', fontsize=18)
    ax[0].set_ylabel('Offset (au)', fontsize=18)

    # Set title
    ax[0].set_title(rf'Point source', fontsize=16)

    # Set colour bar settings and label
    divider = make_axes_locatable(ax[0])
    cax = divider.append_axes(
        'right', 
        size='4%', 
        pad=0.05,
    )
    cb0 = plt.colorbar(
        imlin, 
        cax=cax, 
        orientation = 'vertical', 
    )
    cb0.set_label(
        label = rf'$F_\nu$(MJy/asec$^2$) at {wavelengthum} $\mu$m \& {distance} pc', fontsize=15
    )
    cb0.ax.tick_params(labelsize=15)

    fig.tight_layout()
    fig.show()
    
    #Save figure
    fig.savefig(f'figs/sedimage_point.pdf', dpi=300, facecolor="white")




# -----------------------------------------------------------------------
# Plot things to the Applications parts


if plot_images_convolved_vlti == 'y':
    #
    # Plot figure with convolved images, ie as observed
    # This one with 6 images with VLTI, all phases, three wavelength
    # VLTI at 10um and 200pc
    # VLTI at 3.5um and 200pc
    # VLTI at  1um and 200pc
    #
    # 198 at 10um has a spike at default seed
    #
    distance = 200 # parsec
    wavelengths = [1.6,3.5,10]

    imagelist = [
        f'../r3dresults/st28gm06n052_staranddust_1/190/image_i000_phi000_01.625um.out',
        '../r3dresults/st28gm06n052_staranddust_1/190/image_i000_phi000_03.5um.out',
        '../r3dresults/st28gm06n052_staranddust_1/190/image_i000_phi000_10um.out',
    ]
    # Maximum vertical scales
    vmaxima = [
        4,4,1
    ]
    # Initialise fig-ax for normal images
    figVLTI, axVLTI = plt.subplots(
        2,3, 
        figsize = (12,7.3),
    )


    # Loop through wavelengths/images
    for nwave,image in enumerate(imagelist):

        # Full-width at half-maximum: lambda/(2*baseline)
        fwhmVLTI = wavelengths[nwave]*1e-6 / (2*baselineVLTI) * radian
        sigmaVLTI = fwhmVLTI/2.355

        # For reference, print FWHM:
        print(f'VLTI {wavelengths[nwave]}um FWHM: {fwhmVLTI} mas')


        # Extract path and imagename from image
        imagestrings = re.split('/', image)
        modelname = imagestrings[2]
        imagefilename = imagestrings[4]
        phase = imagestrings[3]
        path = f'{imagestrings[0]}/{imagestrings[1]}/{modelname}/{phase}'

        # Load image
        image2d,image2dlog,flux,axisplot = a3d.load_images(
            path=path,
            image=imagefilename,
            distance=distance
        )
        # Change to Jy/mas2
        image2d = image2d * 1e-6

        # Extract props and compute distance-dependant scales
        Npix = np.shape(image2d)[0] # Number of pixels along one side
        auperpixel = np.abs(2*np.abs(axisplot[0])/Npix)  # number of au per pixel
        masperpixel = auperpixel/distance * 1000  # number of mas per pixel
        size_milliasec = masperpixel * Npix # size of image in mas
        maxflux = np.max(image2d) # max pixel flux

        # Image axis-limits
        axisplotmilliasec = [0.5*size_milliasec,-0.5*size_milliasec,-0.5*size_milliasec,0.5*size_milliasec]

        # Change sigma to order in number of pixels
        sigmaVLTI_pixels = sigmaVLTI / masperpixel

        # Convolve with Gaussian filter
        image2d = scipy.ndimage.gaussian_filter(
            image2d,
            sigma = sigmaVLTI_pixels
        )

        # Plot image
        imVLTI = axVLTI[0][nwave].imshow(
            image2d, 
            origin='lower', extent=axisplotmilliasec, 
            cmap=plt.get_cmap('hot'),
            vmin=0, vmax=vmaxima[nwave]
        )
        axVLTI[0][nwave].tick_params(axis='both', which='major', labelsize=15)

        # FWHM-cirle to show beam
        axVLTI[0][nwave].add_patch(
            plt.Circle(
                (0.75*axisplotmilliasec[0],0.75*axisplotmilliasec[1]), 
                radius=0.5*fwhmVLTI, color='cyan', fill=True)
        )
        # Star-radius-circle to show size of stellar disc (two circles to make it more visible)
        axVLTI[0][nwave].add_patch(
            plt.Circle(
                (0,0), 
                radius=1.65/distance*1000, 
                color='lime', fill=False, linestyle=':', linewidth=2
            )
        )
        axVLTI[0][nwave].add_patch(
            plt.Circle(
                (0,0), 
                radius=1.65/distance*1000, 
                color='b', fill=False, linestyle=':', linewidth=1
            )
        )

        # Plot contours in next line 
        # First set contour levels based on max flux
        starlevels = np.linspace(
            0.1*maxflux,maxflux,6
        )
        backlevels = np.logspace(
            np.log10(1e-6 * 0.01*maxflux),np.log10(0.01*maxflux),5
        )
        # Plot with background
        axVLTI[1][nwave].contour(
            image2d,
            origin='lower', extent=axisplotmilliasec,
            colors='k',
            levels=backlevels,
            linewidths=1
        )
        # Then contours for the star
        axVLTI[1][nwave].contour(
            image2d,
            origin='lower', extent=axisplotmilliasec,
            colors='r',
            levels=starlevels,
            linewidths=1
        )
        # Flip x-axis
        axVLTI[1][nwave].invert_xaxis()
        # Force equal axis
        axVLTI[1][nwave].set_aspect('equal', 'box')

        # FWHM-cirle to show beam
        axVLTI[1][nwave].add_patch(
            plt.Circle(
                (0.75*axisplotmilliasec[0],0.75*axisplotmilliasec[1]), 
                radius=0.5*fwhmVLTI, color='green', fill=True, zorder=10)
        )
        # Star-radius-circle to show size of stellar disc
        axVLTI[1][nwave].add_patch(
            plt.Circle(
                (0,0), 
                radius=1.65/distance*1000, 
                color='b', fill=False, linestyle='--', linewidth=2, 
                zorder=10
            )
        )
        axVLTI[1][nwave].tick_params(axis='both', which='major', labelsize=15)


        # Set titles for first row for both figs
        axVLTI[0][nwave].set_title(rf'{wavelengths[nwave]}\,$\mu$m', fontsize=15)

        # Set xlabels for final row for both figs
        axVLTI[1][nwave].set_xlabel('Offset (mas)',fontsize=18)

        # Set colour bar settings and label
        divider = make_axes_locatable(axVLTI[0][nwave])
        cax = divider.append_axes(
            'right', 
            size='5%', 
            pad=0.05
        )
        cb0 = plt.colorbar(
            imVLTI, 
            cax=cax, 
            orientation = 'vertical'
        )
        if nwave == 2:
            cb0.set_label(
                label = rf'$F_\nu$(Jy/mas$^2$) at {distance}\,pc',fontsize= 15
            )
        cb0.ax.tick_params(labelsize=15)

    # Set ylabels
    axVLTI[0][0].set_ylabel('Offset (mas)',fontsize=18)
    axVLTI[1][0].set_ylabel('Offset (mas)',fontsize=18)


    # 190 at 10um
    # x: 1 till 4       -> mas  5 -> 20  (15)
    # y: -2.5 till 0.5  -> mas  -12.5 -> 2.5 (15)
    axVLTI[1][2].add_patch(
        plt.Rectangle((5,-12.5), 15, 15, color='orange', fill=False, zorder=10)
    )

    # Final settings for figures, save and show if you want to
    figVLTI.tight_layout()

    #Save figure
    figVLTI.savefig(f'figs/images_VLTI.pdf', dpi=300, facecolor="white")

    #figVLTI.show()
    #figcontour.show()


# -------------------------------------------------------------------------------
# Compute various observables
#

if compute_luminosities == 'y':
    # Compute luminosities of all SEDs without dust

    # Load st28gm06n052_nodust
    #   186, 190, 198
    # Load st28gm06n056
    #   140, 141, 142
    #
    paths = [
        '../r3dresults/st28gm06n052_nodust/',
        '../r3dresults/st28gm06n056/',
    ]
    phaseses = [
        ['186','190','198'],
        ['140','141','142']
    ]
    spectra = [
        '/spectrum_i000_phi000.out',
        '/spectrum_i090_phi000.out',
        '/spectrum_i090_phi090.out',
        '/spectrum_i090_phi270.out',
        '/spectrum_i180_phi000.out',
        '/spectrum_i270_phi000.out'
    ]

    for npath, path in enumerate(paths):
        print(path)
        phases = phaseses[npath]
        for phase in phases:

            lumsum = 0
            
            for spectrum in spectra:

                wavelength,sed = a3d.load_spectrum(
                    path = path+phase+spectrum
                )

                # Print all sed-luminosities
                luminosity = a3d.compute_luminosity(
                    wavelengths=wavelength,
                    spectrum=sed
                )
                luminosity = luminosity/Lsol
                print(f'{phase}{spectrum}: {luminosity}')

                lumsum += luminosity
            
            print(f'    Average: {lumsum/len(spectra)}')


if compute_tenmicronfluxdensities == 'y':

    # Computes total stellar flux density in various directions at only 10um
    # and totalt dust flux density in various directions (and wavelengths?)
    # lists these and Fdust/Fstar-ratios
    #
    # In debris discs we used Ldust/Lstar-ratios also...

    # Load stellar flux denisty:
    # r3dresults/st28gm06n052_nodust

    # Load dust flux density:
    # r3dresults/st28gm06n052_nostar

    models = [
        'st28gm06n052_staranddust_nospikes',
        'st28gm06n052_nodust',
        'st28gm06n052_nostar'
    ]
    contents = [
        'Total flux'
        'Stellar flux',
        'Dust flux'
    ]
    phases = [
        '186',
        '190',
        '198'
    ]
    spectra = [
        'spectrum_i000_phi000.out',
        'spectrum_i090_phi000.out',
        'spectrum_i090_phi090.out',
        'spectrum_i090_phi270.out',
        'spectrum_i180_phi000.out',
        'spectrum_i270_phi000.out'
    ]
    angles = [
        '0-0',
        '90-0',
        '90-90',
        '90-270',
        '180-0',
        '270-0'
    ]

    # Loop through models and print lsits for tables
    print('Total flux                   Stellar Flux                   Dust flux                   Contrast (Fd/F*)                Ratio (Ld/L*)')
    for phase in phases:
        print(f'  {phase}                            {phase}                     {phase}')

        for nangle,spectrum in enumerate(spectra):

            # Reset numbers
            totalflux = 0
            starflux = 0
            dustflux = 0

            # Load total fluxes, at 200pc
            # Load stellar fluxes, at 200pc, comparing to images I show
            wavelengths,totalfluxes = a3d.load_spectrum(
                path = f'../r3dresults/st28gm06n052_staranddust_nospikes/{phase}/{spectrum}',
                distance = 200
            )
            tenmicron = np.where(np.array(wavelengths) > 10)[0][0]
            totalflux = 0.5 * (
                totalfluxes[tenmicron-1] + totalfluxes[tenmicron]
            )



            # Load stellar fluxes, at 200pc, comparing to images I show
            wavelengths,starfluxes = a3d.load_spectrum(
                path = f'../r3dresults/st28gm06n052_nodust/{phase}/{spectrum}',
                distance = 200
            )
            tenmicron = np.where(np.array(wavelengths) > 10)[0][0]
            starflux = 0.5 * (
                starfluxes[tenmicron-1] + starfluxes[tenmicron]
            )
            # Also load stellar luminosity
            starlum = a3d.compute_sed_luminosity(
                path = f'../r3dresults/st28gm06n052_nodust/{phase}/{spectrum}',
                distance = 200
            )



            # Load dust fluxes
            wavelengths,dustfluxes = a3d.load_spectrum(
                path = f'../r3dresults/st28gm06n052_nostar/{phase}/{spectrum}',
                distance = 200
            )
            tenmicron = np.where(np.array(wavelengths) > 10)[0][0]
            dustflux = 0.5 * (
                dustfluxes[tenmicron-1] + dustfluxes[tenmicron]
            )
            # Also load dust luminosity
            dustlum = a3d.compute_sed_luminosity(
                path = f'../r3dresults/st28gm06n052_nostar/{phase}/{spectrum}',
                distance = 200
            )

            # Print everything
            print(f'    {angles[nangle]}:{totalflux}           {starflux}           {dustflux}           {dustflux/starflux}           {dustlum/starlum}')


if compute_spectralindeces == 'y':
    #
    # https://en.wikipedia.org/wiki/Spectral_index
    # Given frequency \nu in Hz and radiative flux density S_{\nu } in Jy,
    # the spectral index \alpha is given implicitly by 
    # a = d log(S) / d log(freq)
    #
    # so that
    # S propto freq^a
    #
    # a spectral index of -0.1 to 2 at radio frequencies often indicates thermal emission,
    # for example, a black body, i the RJ-regime is
    # B propto freq^2

    # Load SEDs:
    # Normal with dust
    # Normal, no dust
    # Darwin with dust
    modelpaths = [
        '../r3dresults/st28gm06n052_staranddust_nospikes/',
        '../r3dresults/st28gm06n052_nodust/',
        '../r3dresults/st28gm06n052_darwinsource/'
    ]

    for modelpath in modelpaths:

        paths = [
            f'{modelpath}190/spectrum_i000_phi000.out',
            f'{modelpath}190/spectrum_i090_phi000.out',
            f'{modelpath}190/spectrum_i090_phi090.out',
            f'{modelpath}190/spectrum_i090_phi270.out',
            f'{modelpath}190/spectrum_i180_phi000.out',
            f'{modelpath}190/spectrum_i270_phi000.out',
        ]


        print(modelpath)

        for path in paths:

            # Load wavelength(um) and spectrum(Jy)
            wavelengths,spectrum = a3d.load_spectrum(
                path = path,
                distance = 1
            )
            wavelengths = np.array(wavelengths)
            spectrum = np.array(spectrum)

            # Extract some wavelengths
            thirtymicron = np.where(wavelengths > 30)[0][0]-1
            fiftymicron = np.where(wavelengths > 50)[0][0]
            hundredmicron =  np.where(wavelengths > 100)[0][0]-1
            twohundredmicron =  np.where(wavelengths > 200)[0][0]

            # Extract frequency
            freq = c/(wavelengths*1e-6)

            # Compute spectral indeces
            deltafreq = np.log10(freq[thirtymicron] - freq[fiftymicron])
            deltaflux = np.log10(spectrum[thirtymicron] - spectrum[fiftymicron])
            index_fifty = deltaflux/deltafreq

            deltafreq = np.log10(freq[thirtymicron] - freq[hundredmicron])
            deltaflux = np.log10(spectrum[thirtymicron] - spectrum[hundredmicron])
            index_hundred = deltaflux/deltafreq

            deltafreq = np.log10(freq[thirtymicron] - freq[twohundredmicron])
            deltaflux = np.log10(spectrum[thirtymicron] - spectrum[twohundredmicron])
            index_twohundred = deltaflux/deltafreq

            print(f'  {path[-24:-4]}')
            print(f'    30- 50um: {index_fifty}')
            print(f'    30-100um: {index_hundred}')
            print(f'    30-200um: {index_twohundred}')

        print()


if measuredustcloudflux == 'y':
    # Measures dust flux densities in patches that I define
    # Also plots the patches for reference

    fig, ax = a3d.plot_imagesubplots(
        imagelist = [
            '../r3dresults/st28gm06n052_staranddust_1/186/image_i000_phi000_10um.out',
            '../r3dresults/st28gm06n052_staranddust_1/190/image_i000_phi000_10um.out',
            '../r3dresults/st28gm06n052_staranddust_2/198/image_i000_phi000_10um.out',
        ],
        distance = 1,
        scale = 'lin'
    )

    # 186:
    # x: 5.5 till 4
    # y: -4.5 till -3
    ax[0].add_patch(
        plt.Rectangle((5.5,-4.5), -1.5, 1.5, color='cyan', fill=False)
    )

    # 190:
    # x: -1 till -4
    # y: -2.5 till 0.5
    ax[1].add_patch(
        plt.Rectangle((-1,-2.5), -3, 3, color='cyan', fill=False)
    )

    # 198:
    # x: -1 till -6
    # y: -6.5 - -5
    ax[2].add_patch(
        plt.Rectangle((-1,-6.5), -5, 1.5, color='cyan', fill=False)
    )

    # Extract fluxes of these areas
    #phases = [186,190,198]
    phases = [186,190,198]

    xranges = [
        [4,5.5],
        [-4,-1],
        [-6,-1]
    ]
    yranges = [
        [-4.5,-3],
        [-2.5,0.5],
        [-6.5,-5]
    ]

    # Initilize a second figure to check that I actually measure the correct areas
    fig2, ax2 = plt.subplots(1,3)

    # Update distance (in pc)
    distance = 200

    # Loop over phases
    for nphase,phase in enumerate(phases):

        # Load images themselves
        image2d,image2dlog,totalflux,axisplot = a3d.load_images(
            path = f'../r3dresults/st28gm06n052_staranddust_1/{phase}/',
            image = 'image_i000_phi000_10um.out',
            distance = distance
        )
        # Change to Jy/mas2
        image2d = image2d * 1e-6

        # Create axis-array
        axscale = np.linspace(axisplot[0],axisplot[1],image2d.shape[0])

        # rectangle-ranges:
        xrange = np.where((axscale >= xranges[nphase][0]) & (axscale <= xranges[nphase][1]))[0]
        yrange = np.where((axscale >= yranges[nphase][0]) & (axscale <= yranges[nphase][1]))[0]

        # Fluxes within
        fluxes = image2d[yrange[0]:yrange[-1],xrange[0]:xrange[-1]]

        # Extract total flux density of patch
        # Which is the sum of all fluxes times the area of each pixel in mas2
        pixelsize_mas = (np.abs(axisplot[0]) + np.abs(axisplot[1]))/image2d.shape[0]/distance * 1000
        total_patchflux = fluxes.sum() * pixelsize_mas**2

        print(f'{phase}, x: {xranges[nphase][0]} to {xranges[nphase][1]}')
        print(f'     y: {yranges[nphase][0]} to {yranges[nphase][1]}')
        print(f'     max: {fluxes.max()}     mean: {fluxes.mean()} Jy mas-2      Total flux: {total_patchflux} Jy')

        # Plot patches
        ax2[nphase].imshow(fluxes, origin='lower', cmap=plt.get_cmap('hot'))


    fig.show()
    fig2.show()




# Check "radius" of source at various images
if check_smoothedimage_radius == 'y':

    # Standard images
    #
    #    '../r3dresults/st28gm06n052_staranddust_1/186/image_i000_phi000_1.625um.out'
    #    '../r3dresults/st28gm06n052_staranddust_1/190/image_i000_phi000_1.625um.out'
    #    '../r3dresults/st28gm06n052_staranddust_1/198/image_i000_phi000_1.625um.out'
    #
    #    '../r3dresults/st28gm06n052_staranddust_1/186/image_i000_phi000_3.5um.out'
    #    '../r3dresults/st28gm06n052_staranddust_1/190/image_i000_phi000_3.5um.out'
    #    '../r3dresults/st28gm06n052_staranddust_1/198/image_i000_phi000_3.5um.out'
    #
    #    '../r3dresults/st28gm06n052_staranddust_1/186/image_i000_phi000_10um.out'
    #    '../r3dresults/st28gm06n052_staranddust_1/190/image_i000_phi000_10um.out'
    #    '../r3dresults/st28gm06n052_staranddust_2/198/image_i000_phi000_10um.out'

    distance = 200 # parsec
    relativelimit = 0.1
    phase = 190
    wavelength = 10 # um

    # Load correct seed also
    if phase == 198 and wavelength == 10:
        seed = 2
    else:
        seed = 1
    image = f'../r3dresults/st28gm06n052_staranddust_{seed}/{phase}/image_i000_phi000_{wavelength}um.out'

    # Angular resolution of "beam" is then approx:
    fwhmVLTI = wavelength*1e-6 / (2*baselineVLTI) * radian
    sigmaVLTI = fwhmVLTI/2.355


    # Star's radius here (for comparison plot)
    Rstar = 1.65


    # Extract path and imagename from image
    imagestrings = re.split('/', image)
    modelname = imagestrings[2]
    imagefilename = imagestrings[4]
    phase = imagestrings[3]
    path = f'{imagestrings[0]}/{imagestrings[1]}/{modelname}/{phase}'

    # Load image
    image2d,image2dlog,flux,axisplot = a3d.load_images(
        path=path,
        image=imagefilename,
        distance=distance
    )
    # Change to Jy/mas2
    image2d = image2d * 1e-6


    # Extract props and compute distance-dependant scales
    Npix = np.shape(image2d)[0] # Number of pixels along one side
    auperpixel = np.abs(2*axisplot[0]/Npix)  # number of au per pixel
    masperpixel = auperpixel/distance * 1000  # number of mas per pixel
    size_au = auperpixel * Npix # Size of image in au


    # Change sigma to order in number of pixels
    sigmaVLTI_pixels = sigmaVLTI / masperpixel

    # Convolve with Gaussian filter
    image2d = scipy.ndimage.gaussian_filter(
        image2d,
        sigma = sigmaVLTI_pixels
    )


    # Number of circular annulii are a quarter to make 
    # sure there are enough pixels in each annulus
    # And range is half the size of the image

    radial_range = np.linspace(auperpixel,size_au*0.5,int(Npix*0.25))
    radial_fluxes = np.zeros(int(Npix*0.25))
    radial_npixels = np.zeros(int(Npix*0.25))


    # List of pixel numbers with 0 in centrum, radially and along each axis
    rpix = list(range(0,int(Npix*0.5),2))
    xpix = list(range(int(-Npix*0.5),int(Npix*0.5)))
    ypix = list(range(int(-Npix*0.5),int(Npix*0.5)))


    # Loop through annulii, add fluxes of each
    print('Looping through image')
    for nn,nr in enumerate(rpix):
        for nx in xpix:
            for ny in ypix:
                if np.sqrt(nx**2 + ny**2) > nr and \
                    np.sqrt(nx**2 + ny**2) <= nr+2:
                    
                    imagex = int(nx+Npix*0.5)
                    imagey = int(ny+Npix*0.5)

                    # Add fluxes of each annulus and number of pixels
                    radial_fluxes[nn] += image2d[imagey,imagex]
                    radial_npixels[nn] += 1
    print('    done')

    # Average flux per annulus
    radial_fluxes /= radial_npixels

    # At what radius is half maximum? Or 25%? Or 10%?
    fluxlimit = relativelimit*image2d.max()
    if fluxlimit < radial_fluxes.max():
        radius_fluxlimit = radial_range[np.where(radial_fluxes >= fluxlimit)[0].max()]
    else:
        radius_fluxlimit = 0
    print('')
    print(f'Max flux of smoothed image: {image2d.max()} Jy')
    print(f'Max-flux ratio (image2d/annulus): {image2d.max()/radial_fluxes.max()}')
    print(f'{modelname}_{phase}, {imagefilename}')
    print(f'Obs-radius: {radius_fluxlimit} AU')

    # Plot to check
    plt.plot(radial_range,radial_fluxes)
    plt.plot([Rstar,Rstar],[0,radial_fluxes.max()])
    plt.plot([radius_fluxlimit,radius_fluxlimit],[0,fluxlimit])
    plt.show()



# Manually plot radius vs wavelength plots
if plot_smoothedimage_radius == 'y':

    # To plot source-radius vs wavelength
    wavelengths = [1.6,3.5,10.0]

    t1_50perc = [0,2.0587,2.4116]
    t2_50perc = [np.nan,np.nan,1.7057]
    t3_50perc = [0,1.4704,2.0071]

    t1_25perc = [1.9410,2.6469,2.9999]
    t2_25perc = [1.2351,1.9410,2.4116]  # This
    t3_25perc = [0,2.1763,2.6469]

    t1_10perc = [2.5293,2.9999,3.3528]
    t2_10perc = [1.9410,2.4116,3.5881]  # This
    t3_10perc = [2.0587,2.6469,3.1175]

    # Manually change these to whatever you want to plot
    plot_these_lists = [t2_50perc,t2_25perc,t2_10perc]
    legendlist = [
        r'50\,\%',
        r'25\,\%',
        r'10\,\%'
    ]

    # Initialise plot (add more styles if needed)
    fig, ax = plt.figure(figsize=(6, 4)), plt.axes()
    linestyles = ['-','--','-.',':']
    markerstyles = ['v','o','*','^']

    for nn,sourceradius in enumerate(plot_these_lists):

        ax.plot(
            wavelengths,sourceradius,
            linestyle=linestyles[nn],marker=markerstyles[nn],
            markersize = 8,
            label = legendlist[nn]
        )

    ax.plot([1,11],[1.65,1.65],'k:')
    ax.legend(fontsize=13)

    ax.set_xlabel(r'Wavelength ($\mu$m)',fontsize=18)
    ax.set_ylabel(r'Source radius (au)',fontsize=18,)
    ax.set_xlim([1,11])
    ax.tick_params(axis='both', which='major', labelsize=15)

    fig.tight_layout()

    fig.savefig('figs/source_radius.pdf', dpi=300, facecolor="white")
    fig.show()

