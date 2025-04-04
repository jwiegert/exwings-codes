# A local API for more easily explore the data within r3dresults.
#
# Tool for fast visualisation of local data kept in ../r3dresults/*
# Run code in python and open in brwoser at 
# http://127.0.0.1:8050/
# or at whatever port you want to set in the end of this code.
#
# -------------------------------------------------------------------
#
# Standard packages
import os
import re

# Custom packages
import analyze_r3d_functions as a3d

# Dashboard packages
from dash.dependencies import Input, Output
from dash import dcc, html
import dash
import dash_bootstrap_components as dbc

# Image specific packages
#
# Matplotlib's normal usage of TK is not compatible with Dashboards.
# Change this to Agg instead. Info at:
# https://stackoverflow.com/questions/27147300/matplotlib-tcl-asyncdelete-async-handler-deleted-by-the-wrong-thread
#
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import plotly
import plotly.tools as tls
import base64

# -------------------------------------------------------------------
# Definitions
path='../r3dresults/'
Lsol = 3.828e26 # W


# Settings for Dashboard

# List of models (folders) inside ../r3dresults
modelnames = [
    modelname for modelname in os.listdir(path=path) if os.path.isdir(path+modelname) == True
]
modelnames.sort()

# Themes
stylesheets = [dbc.themes.MATERIA]

# Create Dash object, app
app = dash.Dash(
    __name__, 
    external_stylesheets=stylesheets,
    meta_tags=[dict(name='viewport', content='width=device-width, initial-scale=1.0')]
)


# App layout 
#
# NOTE
# html.H1 = header 1, main title
# html.H2 = header 2, sub title
# html.P = paragraph
# dropdown creates a dropdown manue, 
# id is it's name, 
# value is its default value. 
# Options is a list of what you can chose.
# placeholder is a default placeholder text on unchosen menus
# Etc see more online
app.layout = dbc.Container([

    html.H2('R3D results-plotter:'),

    # Menu of model names, ie folders inside ../r3dresults/
    dcc.Dropdown(
        id='model-dropdown', 
        className='',
        options=[
            {'label':modelname, 'value':modelname}  for modelname in modelnames
        ],
        placeholder='Chose a model',
        value=modelnames[0]
    ),

    # Menu of folders (phases) inside chosen model
    dcc.Dropdown(
        id='phase-dropdown', 
        className='',
        placeholder='Chose a phase',
        value=0
    ),

    # Horisontal line
    html.Hr(),

    # Plot singular SED below here
    html.P('SEDs:'),

    # Menu of Inclinations and phi angles
    dcc.Dropdown(
        id='incl-dropdown', 
        className='',
        placeholder='Chose a spectrum',
        value='N/A'
    ),

    # SED plot
    dcc.Graph(
        id='sed-plot-one'
    ),

    # Print Bol. luminosity of the SED
    html.P('SED-luminosity (Lsol):'),
    html.Div(id='sed_luminosity'),


    # Horisontal line
    html.Hr(),

    # Plot singular image here
    html.P('Images:'),

    # Menu of image*.out-files
    dcc.Dropdown(
        id='image-dropdown', 
        className='',
        placeholder='Chose an image',
        value=0
    ),

    # Image plot
    html.Img(
        id='image-plot-one',
    ),


    # Horisontal line
    html.Hr(),

    # Plot a choice of SEDs in same figure
    html.P('Plot all SEDs on top of each other'),

    # SED-chooser
    dcc.RadioItems(
        id='sed-picker', 
        options=[
            {'label':'All phases (constant inclination)', 'value':'allphases'},
            {'label':'All inclinations (constant phase)', 'value':'allincls'}
        ],
        labelStyle={'display': 'block'},
        value=0
    ),

    # Plot all chosen SEDs
    dcc.Graph(
        id='sed-plot-all'
    ),

    # Print average luminosity of all plotted SEDs
    html.P('Mean SED-luminosity (Lsol):'),
    html.Div(id='sed_luminosity_average'),


    # Horisontal line
    html.Hr(),

    # Plot all chosen images vertically
    html.P('Plot all images  Inclination/Phase/Wavelength'),

    # Chose which set of images to plot
    dcc.RadioItems(
        id='image-picker', 
        options=[
            {'label':'All phases (constant incl & wavelength)', 'value':'allphases'},
            {'label':'All incls (constant phase & wavelength)', 'value':'allincls'},
            {'label':'All waves (constant phase & incl)', 'value':'allwaves'},
        ],
        labelStyle={'display': 'block'},
        value=0
    ),

    # Empty space between image-radioitems
    html.Br(), 

    # Image plot
    html.Img(
        id='image-plot-all',
    ),

])
# End of layout

# -----------------------------------------------------------------------
# Functions and Callbacks

# Given chosen modelname, change options for phases-menu
@app.callback(
    Output('phase-dropdown', 'options'),
    Input('model-dropdown', 'value'),
)
def create_phase_dict(modelname):

    # If statement is an attempt to stop looking for models when none is chocen
    if modelname != 0 :

        # List all folders within ../r3dresults/{modelname}
        phases = [
            phase for phase in os.listdir(path=f'{path}{modelname}/') if os.path.isdir(f'{path}{modelname}/{phase}') == True
        ]

        # Sort folders, or give string "no phases found"
        if len(phases) > 0:
            phases.sort()
        else:
            phases = ['No phases found']

        # Return dict with menu-options
        return [{'label':phase, 'value':phase} for phase in phases]


# Given model and phase, add list of images and SEDs
@app.callback(
    Output('image-dropdown', 'options'),
    Output('incl-dropdown', 'options'),
    Input('model-dropdown', 'value'),
    Input('phase-dropdown', 'value'),
)
def create_image_sed_dicts(modelname,phase):

    # If statement to stop search for folders if non are chosen
    if modelname != 0 or phase != 0:

        # Extract file names in folder
        filenames = os.listdir(path=f'{path}{modelname}/{phase}')
        filenames.sort()

        # Extract image*.out file names
        images = [image for image in filenames if image[0:5] == 'image']

        # Check if there are any images, and sort them, or return no-image-string
        if len(images) > 0:
            images.sort()
        else:
            images = ['No image found']

        # Extract inclinations of SEDs, ie 
        incls = []
        for filename in filenames:
            if filename[:8] == 'spectrum':
                incls.append(re.findall('spectrum_i.*.', filename)[0][9:-4])
        
        # Check if there are any SEDs
        if len(incls) > 0:
            incls.sort()
        else:
            incls = ['N/A']

        # return two dicts, one with image options and one with sed-incl-options
        return [{'label':image, 'value':image} for image in images], [{'label':f'Angles: {incl}', 'value':incl} for incl in incls]


# Given choices of model, phase and sed-inclination, plot the SED
@app.callback(
    Output('sed-plot-one', 'figure'),
    Output('sed_luminosity', 'children'),
    Input('model-dropdown', 'value'),
    Input('phase-dropdown', 'value'),
    Input('incl-dropdown', 'value'),
)
def plot_sed_one(modelname,phase,incl):

    # if-statement to stop searching if nothing is chosen
    if modelname != 0 or phase != 0 or incl != 'N/A':

        # Extract luminosity of SED
        lum = a3d.compute_sed_luminosity(
            path=f'{path}{modelname}/{phase}/spectrum_{incl}.out'
        )/Lsol

        # Create matplotlib-figure-object of SED, and convert to plotly-figure
        fig,ax,maxflux,maxwave = a3d.plot_sed(
            path=f'{path}{modelname}/{phase}/spectrum_{incl}.out'
        )
        plotly_fig = tls.mpl_to_plotly(fig)

        # Change some settings for the plotly plot
        plotly_fig.update_layout(
            title_font_size=18,
            hoverlabel_font_size=18,
            plot_bgcolor='white',
            yaxis = dict(tickfont = dict(size=16)),
            xaxis = dict(tickfont = dict(size=16))
        )
        plotly_fig.update_xaxes(
            titlefont_size=16,
            showgrid=True, gridwidth=1, gridcolor='lightgrey',
            ticks="outside", tickwidth=2, ticklen=10,
            showline=True, linewidth=2, linecolor='black', mirror=False
        )
        plotly_fig.update_yaxes(
            titlefont_size=16,
            showgrid=True, gridwidth=1, gridcolor='lightgrey',
            ticks="outside", tickwidth=2, ticklen=10,
            showline=True, linewidth=2, linecolor='black', mirror=False
        )

        # Return plotly figure obect and luminosity of the SED
        return plotly_fig, lum


# Given model, phase and image-file, plot this
@app.callback(
    Output('image-plot-one', 'src'),
    Input('model-dropdown', 'value'),
    Input('phase-dropdown', 'value'),
    Input('image-dropdown', 'value'),
)
def plot_image_one(modelname,phase,image):

    # If-statement to stop searching if nothing is chosen
    if modelname != 0 or phase != 0 or image != 0:

        # If temporary image-file already exists, remove it so the image updates
        if os.path.exists('temp.png') == True:
            os.system('rm temp.png')

        # Put image-choice to a list (since plot_images loops through list)
        image = [image]

        # Create image-figure-object and ax-object, save as temporary png-file
        fig,ax,testflux = a3d.plot_images(
            path = f'{path}{modelname}/{phase}/',
            images = image,
            distance = 1
        )
        plt.savefig(
            'temp.png', dpi=120, bbox_inches='tight'
        )

        # Change to base64-encoding that html.Img can read the temporary png
        tempbase64 = base64.b64encode(open('temp.png', 'rb').read()).decode('ascii')
        impath = 'data:image/png;base64,{}'.format(tempbase64)

        # Return path to base64-temp-image
        return impath


# Given model, plot all seds of your choice
@app.callback(
    Output('sed-plot-all', 'figure'),
    Output('sed_luminosity_average','children'),
    Input('model-dropdown', 'value'),
    Input('phase-dropdown', 'value'),
    Input('incl-dropdown', 'value'),
    Input('sed-picker', 'value'),
)
def plot_sed_all(modelname,phase,incl,choice):

    # If statement to stop searching if nothing is chosen
    if modelname != 0 or phase != 0 or incl != 'N/A' or choice != 0:

        # Initiate empty pathlist and legendlist 
        # (legend is either inclinations or phases, also used for legend in figure)
        # these depend on choice.
        pathlist = []
        legendlist = []

        # Create list of SEDs of all phases (constant - pre-chosen - incl)
        if choice == 'allphases':

            # List of phase-folders
            legendlist = [
                phase for phase in os.listdir(path=f'{path}{modelname}/') if os.path.isdir(f'{path}{modelname}/{phase}') == True
            ]

            # List of spectrum*.out-files
            pathlist = [f'{path}{modelname}/{phase}/spectrum_{incl}.out' for phase in legendlist]

        # Else, fill list with paths to chosen phsae and all incls there instead
        if choice == 'allincls':

            # List of files in phase-folder
            filenames = os.listdir(f'{path}{modelname}/{phase}/')

            # Extract all spectrum*.out files
            for filename in filenames:
                if filename[:8] == 'spectrum':
                    # Save all available inclinations (numbers as str only)
                    legendlist.append(re.findall('spectrum_i.*.', filename)[0][9:-4])

            # Save a list of all spectrum-files
            pathlist = [f'{path}{modelname}/{phase}/spectrum_{incl}.out' for incl in legendlist]

        # Create SED figure object and convert to plotly-figure
        fig,ax = a3d.plot_sedsmany(
            pathlist=pathlist,
            legendlist=legendlist,
            distance=1
        )
        plotly_fig_all = tls.mpl_to_plotly(fig)


        # Compute average luminosity of all chosen SEDs
        lumtot = 0
        for pathstr in pathlist:
            lum = a3d.compute_sed_luminosity(
                path=pathstr
            )
            lumtot += lum
        lumtot /= len(pathlist)*Lsol


        # Change some settings for the plotly plot
        plotly_fig_all.update_layout(
            title_font_size=18,
            hoverlabel_font_size=18,
            plot_bgcolor='white',
            yaxis = dict(tickfont = dict(size=16)),
            xaxis = dict(tickfont = dict(size=16))
        )
        plotly_fig_all.update_xaxes(
            titlefont_size=16,
            showgrid=True, gridwidth=1, gridcolor='lightgrey',
            ticks="outside", tickwidth=2, ticklen=10,
            showline=True, linewidth=2, linecolor='black', mirror=False
        )
        plotly_fig_all.update_yaxes(
            titlefont_size=16,
            showgrid=True, gridwidth=1, gridcolor='lightgrey',
            ticks="outside", tickwidth=2, ticklen=10,
            showline=True, linewidth=2, linecolor='black', mirror=False
        )

        # Return plotly-figure and average luminosity
        return plotly_fig_all,lumtot




# TODO fix so this plots with gamma function plotter
# Given model, phase, one-image-plotter, choice and scale, plot several images
@app.callback(
    Output('image-plot-all', 'src'),
    Input('model-dropdown', 'value'),
    Input('phase-dropdown', 'value'),
    Input('image-dropdown', 'value'),
    Input('image-picker', 'value'),
)
def plot_image_all(modelname,phase,image,choice):
    #
    # Modelname is foldername in r3dresults
    # Phase is phase-choice
    # Image is imagename, ie choice of wavelength and inclination
    # Choice chooses which two attributes (of three) are to be constant
    # Scale is lin or log
    #
    # If temporary image file already exists, remove it so that the image updates
    if os.path.exists('tempsubplots.png') == True:
        os.system('rm tempsubplots.png')

    # Extract chosen image inclination and wavelength
    # Filename is formatted as image_i{incl}_phi{phi}_{wavelength}um.out
    # EG image_i000_phi000_100um.out
    imageincl = re.split('_',image)[1][1:]
    imagephi = re.split('_',image)[2][3:]  #
    imagewave = re.split('_',image)[3][:-6]

    
    # Initate path to model-folder
    path = f'../r3dresults/{modelname}/'

    # Initate image list
    imagelist = []

    # Check the choice, do we want all phases, all incls or all wavelengths?
    #
    # Choice 1: all phases, rest are constant
    if choice == 'allphases':
        # Extract all phases in models-folder
        # and use chosen wavelength and inclination from image-choice earlier

        # Loop over phases listed in model-folder
        for phases in [
            folder for folder in os.listdir(path=path) if os.path.isdir(path+folder) == True
        ]:
            # Save one image-file-path from each phase-folder
            imagelist.append(f'{path}{phases}/image_i{imageincl}_phi{imagephi}_{imagewave}um.out')
            # And sort them
            imagelist.sort()

    # Choice 2: all inclinations
    if choice == 'allincls':
        
        # Save a list of images with chosen phase and wavelength, and all available incls
        for file in os.listdir(path=f'{path}{phase}/'):
            if file[:5] == 'image':
                if re.split('_', file)[3] == f'{imagewave}um.out':
                    imagelist.append(f'{path}{phase}/{file}')
                    imagelist.sort()

    # Choice 3: all wavelengths
    if choice == 'allwaves':

        # Save a list of images with chosen phase and inclination, and all available wavelengths
        for file in os.listdir(path=f'{path}{phase}/'):
            if file[:5] == 'image':
                if re.split('_', file)[1] == f'i{imageincl}':
                    imagelist.append(f'{path}{phase}/{file}')
                    imagelist.sort()


    # Create fig and ax-objects
    fig,ax = a3d.plot_imagesubplots(
        imagelist=imagelist,
        distance=1,
    )
    fig.tight_layout()
    # Save figure as temporary png-file
    plt.savefig(
        'tempsubplots.png', dpi=150, bbox_inches='tight'
    )
    # Change to base64-encoding that html.Img can read
    tempbase64 = base64.b64encode(open('tempsubplots.png', 'rb').read()).decode('ascii')
    imallpath = 'data:image/png;base64,{}'.format(tempbase64)

    # Return base64-image-path
    return imallpath


# Final part initiates the server only if this is the main running-script
if __name__ == "__main__":
    # Chose debug mode or not

    #app.run_server(debug=True)
    app.run_server(port=8050)
