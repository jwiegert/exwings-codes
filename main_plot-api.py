# A local API for more easily explore the data within r3dresults.
#
# 20220817: 
# As of now, my results are only from various tests. But in the future there
# will be a lot of data, a lot of simulated images and SEDs, all in various
# viewing angles, resolutions etc etc. For this reason I need a much more
# convenient way of looking at these. As such I designed this simple API that I
# only run locally, to much faster look at results and inter-compare results.
#



#    API't ska automatiskt lista alla kataloger i r3dresults, sen alla i den man väljer där, och sen alla sed-filer, spectrum*.out och alla bildfiler, image*.out
#    Patherna innehåller info om modell, fas, inclination och våglängder
#    Plotta alla bilder med samma inklination, plotta alla bilder med samma våglängd
#    Plotta alla SEDer med samma inklination i samma figur
#        os.listdir(path='../r3dresults/st28gm06n056/140/')    
#    eller
#        import glob
#        glob.glob('../r3dresults/st28gm06n056/140/image*.out')
#    och för kataloger
#        glob.glob('../r3dresults/st28gm06n056/*/')
#    Ger då listor som kan vara inputs till mina menyer

# Standard (important) packages
import os

# Custom packages
import analyze_r3d_functions as a3d

# Dashboard packages
from dash.dependencies import Input, Output
from dash import dcc, html
import dash
import dash_bootstrap_components as dbc

# Image specific packages
import matplotlib.pyplot as plt
import io
import plotly
import plotly.tools as tls
import base64


# -------------------------------------------------------------------
# Definitions
path='../r3dresults/'
imcounter = 0

# List of models (folders) inside ../r3dresults
modelnames = [
    modelname for modelname in os.listdir(path=path) if os.path.isdir(path+modelname) == True
]
modelnames.sort()


# Dash board settings
stylesheets = [dbc.themes.MATERIA]

# Create a dash object called app
app = dash.Dash(__name__, external_stylesheets=stylesheets,
    meta_tags=[dict(name='viewport', content='width=device-width, initial-scale=1.0')]
)


# APP LAYOUT

# html.H1 = header 1, main title
# html.H2 = header 2, sub title
# html.P = paragraph
# dropdown creates a dropdown manue, id is it's name, value is
# its default value. Options is a list of what you can chose.


app.layout = dbc.Container([

    html.H2('R3D results-plotter:'),

    # Menu of model names inside ../r3dresults/
    dcc.Dropdown(
        id='model-dropdown', 
        className='',
        options=[
            {'label':modelname, 'value':modelname}  for modelname in modelnames
        ],
        placeholder='Chose a model'
    ),

    # Menu of folders inside chosen model
    dcc.Dropdown(
        id='phase-dropdown', 
        className='',
        placeholder='Chose a phase'
    ),

    html.Hr(),
    # Plot singular SED below here
    html.P('SEDs:'),

    # Menu of SED-files
    dcc.Dropdown(
        id='sed-dropdown', 
        className='',
        placeholder='Chose an SED'
    ),

    # SED plot
    dcc.Graph(
        id='sed-plot-one'
    ),
    html.P('SED-luminosity (Lsol):'),
    html.Div(id='sed_luminosity'),


    html.Hr(),
    # Plot singular image here
    html.P('Images:'),

    # Menu of image-files
    dcc.Dropdown(
        id='image-dropdown', 
        className='',
        placeholder='Chose an image'
    ),

    # Image plot
    html.Img(
        id='image-plot-one',
    ),

    # Horisontal line, separator
    html.Hr(),
    html.P('Plots all SEDs on top of each other'),

    # TODO
    # ie plot all SED-inclinations
    # add a choice, plot or not? if yes, plot all SEDs
    # Or on top of each other?!
    # Needs a new functions in a3d

    # Horisontal line, separator
    html.Hr(),
    html.P('Plots all images  Inclination/Phase/Wavelength')

    # TODO
    # if nothing is chosen, show nothing
    # how to do this?


])


# FUNCTIONS AND CALLBACKS

# Given chosen modelname, change options for phases-menu
@app.callback(
    Output('phase-dropdown', 'options'),
    Input('model-dropdown', 'value'),
)
def create_phase_dict(modelname):

    phases = [
        phase for phase in os.listdir(path=f'{path}{modelname}/') if os.path.isdir(f'{path}{modelname}/{phase}') == True
    ]

    if len(phases) > 0:
        phases.sort()

    else:
        phases = ['No phases found']

    return [{'label':phase, 'value':phase} for phase in phases]


# Given model and phase, add list of images and SEDs
@app.callback(
    Output('image-dropdown', 'options'),
    Output('sed-dropdown', 'options'),
    Input('model-dropdown', 'value'),
    Input('phase-dropdown', 'value'),
)
def create_image_sed_dicts(modelname,phase):

    # Extract file names in folder
    filenames = os.listdir(path=f'{path}{modelname}/{phase}')
    filenames.sort()


    # Extract image file names
    images = [image for image in filenames if image[0:5] == 'image']

    # Check if there are any
    if len(images) > 0:
        images.sort()
    else:
        images = ['No image found']


    # Extract SED file names
    seds = [sed for sed in filenames if sed[0:8] == 'spectrum']

    # Check if there are any SEDs
    if len(seds) > 0:
        seds.sort()
    else:
        seds = ['No SED found']

    return [{'label':image, 'value':image} for image in images], [{'label':sed, 'value':sed} for sed in seds]


# Given model, phase and sed-file, plot this
@app.callback(
    Output('sed-plot-one', 'figure'),
    Output('sed_luminosity', 'children'),
    Input('model-dropdown', 'value'),
    Input('phase-dropdown', 'value'),
    Input('sed-dropdown', 'value'),
)
def plot_sed_one(modelname,phase,sed):

    # Extract luminosity of SED
    lum = a3d.compute_sed_luminosity(
        path=f'{path}{modelname}/{phase}/{sed}'
    )/3.828e26

    # Extract image
    fig,ax,maxflux,maxwave = a3d.plot_sed(
        path=f'{path}{modelname}/{phase}/{sed}'
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
    plt.savefig('temp.png', dpi=120, bbox_inches='tight')

    # Change to base64-encoding that html.Img can read
    tempbase64 = base64.b64encode(open('temp.png', 'rb').read()).decode('ascii')
    impath = 'data:image/png;base64,{}'.format(tempbase64)

    return impath





if __name__ == "__main__":
    app.run_server(debug=True)

