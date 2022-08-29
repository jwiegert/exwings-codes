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

    # Menu of Inclinations
    dcc.Dropdown(
        id='incl-dropdown', 
        className='',
        placeholder='Chose an Inclination',
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
    html.P('Plots all SEDs on top of each other'),

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
    html.P('Plots all images  Inclination/Phase/Wavelength'),

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

    # Chose scale, linear or logarithmic
    dcc.RadioItems(
        id='image-scale-picker', 
        options=[
            {'label':'Linear flux scale', 'value':'lin'},
            {'label':'Logarithmic flux scale', 'value':'log'},
        ],
        labelStyle={'display': 'block'},
        value='lin'
    ),

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
                incls.append(re.findall('spectrum_i.*.', filename)[0][10:-4])
        
        # Check if there are any SEDs
        if len(incls) > 0:
            incls.sort()
        else:
            incls = ['N/A']

        # return two dicts, one with image options and one with sed-incl-options
        return [{'label':image, 'value':image} for image in images], [{'label':f'Inclination: {incl} deg', 'value':incl} for incl in incls]


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

        # TODO continue cleaning from here

        # Extract luminosity of SED
        lum = a3d.compute_sed_luminosity(
            path=f'{path}{modelname}/{phase}/spectrum_i{incl}.out'
        )/Lsol

        # Extract image
        fig,ax,maxflux,maxwave = a3d.plot_sed(
            path=f'{path}{modelname}/{phase}/spectrum_i{incl}.out'
        )
        plotly_fig = tls.mpl_to_plotly(fig)

        # Change some settings on the plotly plot
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

        return plotly_fig, lum



# Given model, phase and image-file, plot this
@app.callback(
    Output('image-plot-one', 'src'),
    Input('model-dropdown', 'value'),
    Input('phase-dropdown', 'value'),
    Input('image-dropdown', 'value'),
)
def plot_image_one(modelname,phase,image):
    if modelname != 0 or phase != 0 or image != 0:

        # If image already exists, remove it so that the image updates
        if os.path.exists('temp.png') == True:
            os.system('rm temp.png')

        # Change image to list
        image = [image]

        # Load and save image as png
        fig,ax,testflux = a3d.plot_images(
            path = f'{path}{modelname}/{phase}/',
            images = image,
            distance = 1
        )
        plt.savefig(
            'temp.png', dpi=120, bbox_inches='tight'
        )

        # Change to base64-encoding that html.Img can read
        tempbase64 = base64.b64encode(open('temp.png', 'rb').read()).decode('ascii')
        impath = 'data:image/png;base64,{}'.format(tempbase64)

        return impath



# Given model, phase and all phases or all inclinations, plot all seds of choice
@app.callback(
    Output('sed-plot-all', 'figure'),
    Output('sed_luminosity_average','children'),
    Input('model-dropdown', 'value'),
    Input('phase-dropdown', 'value'),
    Input('incl-dropdown', 'value'),
    Input('sed-picker', 'value'),
)
def plot_sed_all(modelname,phase,incl,choice):
    if modelname != 0 or phase != 0 or incl != 'N/A' or choice != 0:

        # pathlist and legendlist depends on choice.
        pathlist = []
        legendlist = []

        # All phases, ie inclination is chosen earlier, take all phases with this
        if choice == 'allphases':

            legendlist = [
                phase for phase in os.listdir(path=f'{path}{modelname}/') if os.path.isdir(f'{path}{modelname}/{phase}') == True
            ]

            for phase in legendlist:
                pathlist.append(f'{path}{modelname}/{phase}/spectrum_i{incl}.out')

        # Else, fill list with paths to chosen phsae and all incls there instead
        if choice == 'allincls':
            
            filenames = os.listdir(f'{path}{modelname}/{phase}/')

            for filename in filenames:
                if filename[:8] == 'spectrum':
                    legendlist.append(re.findall('spectrum_i.*.', filename)[0][10:-4])

            for incl in legendlist:
                pathlist.append(f'{path}{modelname}/{phase}/spectrum_i{incl}.out')

        # Create SED figure object
        fig,ax = a3d.plot_sedsmany(
            pathlist=pathlist,
            legendlist=legendlist,
            distance=1
        )
        plotly_fig_all = tls.mpl_to_plotly(fig)


        # Compute average luminosity of chosen SEDs
        lumtot = 0
        for pathstr in pathlist:
            lum = a3d.compute_sed_luminosity(
                path=pathstr
            )
            lumtot += lum
        lumtot /= len(pathlist)*Lsol


        # Change some settings on the plotly plot
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

        return plotly_fig_all,lumtot


@app.callback(
    Output('image-plot-all', 'src'),
    Input('model-dropdown', 'value'),
    Input('phase-dropdown', 'value'),
    Input('image-dropdown', 'value'),
    Input('image-picker', 'value'),
    Input('image-scale-picker', 'value')
)
def plot_image_all(modelname,phase,image,choice,scale):

    # modelname is foldername in r3dresults
    # phase is phase-choice
    # image is imagename, ie choice of wavelength and inclination
    # choice chooses which two attributes are to be constant
    # scale is lin or log

    # If image already exists, remove it so that the image updates
    if os.path.exists('tempsubplots.png') == True:
        os.system('rm tempsubplots.png')


    # Extract chocen image inclination and wavelength
    imageincl = re.split('_',image)[1][1:]
    imagewave = re.split('_',image)[2][:-6]

    # Initate path
    path = f'../r3dresults/{modelname}/'

    # Initate image list
    imagelist = []


    # TODO NÅGOT FUNKAR INTE HÄR!!
    # Check the choice, all phases, all incls or all waves?
    if choice == 'allphases':
        # Extract all phases in models-folder
        # and use chosen wavelength and inclination from image
        for phases in [
            folder for folder in os.listdir(path=path) if os.path.isdir(path+folder) == True
        ]:
            imagelist.append(f'{path}{phases}/image_i{imageincl}_{imagewave}um.out')
            imagelist.sort()

    # TODO något funkar inte här!!!
    # fixa...
    if choice == 'allincls':
        # Save a list of images with chosen phase and wavelength but all incls
        for file in os.listdir(path=f'{path}{phase}/'):
            if file[:5] == 'image':
                if re.split('_', file)[2] == f'{imagewave}um.out':
                    imagelist.append(f'{path}{phase}/{file}')
                    imagelist.sort()

    if choice == 'allwaves':
        # Save a list of images with chosen phase and inclination, but all wavelengths
        for file in os.listdir(path=f'{path}{phase}/'):
            if file[:5] == 'image':
                if re.split('_', file)[1] == f'i{imageincl}':
                    imagelist.append(f'{path}{phase}/{file}')
                    imagelist.sort()


    # Create fig and ax-objects
    fig,ax = a3d.plot_imagesubplots(
        imagelist=imagelist,
        distance=1,
        scale=scale
    )
    
    # Convert matplotlib-fig to png in base64
    plt.savefig(
        'tempsubplots.png', dpi=150, bbox_inches='tight'
    )
    # Change to base64-encoding that html.Img can read
    tempbase64 = base64.b64encode(open('tempsubplots.png', 'rb').read()).decode('ascii')
    imallpath = 'data:image/png;base64,{}'.format(tempbase64)

    return imallpath


# Chose debug mode or not
if __name__ == "__main__":
    #app.run_server(debug=True)
    app.run_server(port=8050)

