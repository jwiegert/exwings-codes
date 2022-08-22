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


# Dashboard packages
from dash.dependencies import Input, Output
from dash import dcc, html
import dash
import dash_bootstrap_components as dbc

# -------------------------------------------------------------------
# Definitions
path='../r3dresults/'





# -------------------------------------------------------------------
# Settings for menues

# Dict of models (folders) inside ../r3dresults
modelnames = [
    modelname for modelname in os.listdir(path=path) if os.path.isdir(path+modelname) == True
]
modelnames.sort()
modelnames_dict = {modelname:modelname for modelname in modelnames}

# Initial phases_dict
modelname = modelnames[0]
phases = [
    phase for phase in os.listdir(path=f'{path}{modelname}/') if os.path.isdir(f'{path}{modelname}/{phase}') == True
]
phases.sort()
phases_dict = {phase:phase for phase in phases}



# Dicts of images and SEDs within phases
# TODO
#filenames = os.listdir(path=path+modelname+phase)
#filenames.sort()

#images = [image for image in filenames if image[0:5] == 'image']
#seds = [sed for sed in filenames if sed[0:8] == 'spectrum']




# Interactive stuff settings





# Dash board settings
stylesheets = [dbc.themes.MATERIA]

# Create a dash object called app
app = dash.Dash(__name__, external_stylesheets=stylesheets,
    meta_tags=[dict(name="viewport", content="width=device-width, initial-scale=1.0")]
)




# Define the layout
# html.H1 = header 1, main title
# html.H2 = header 2, sub title
# html.P = paragraph
# dropdown creates a dropdown manue, id is it's name (for later), value is
# its default value. Options is a list of what you can chose.
# So value='' gives the default value for each dcc.
# The objects in app.layout ends up on the page in the order we have them written here


app.layout = dbc.Container([

    dcc.Dropdown(
        id='model-dropdown', 
        className='',
        options=modelnames_dict,
        value=modelnames[0]
    ),

    dcc.Dropdown(
        id='phase-dropdown', 
        className='',
        options=phases_dict,
        value=''
    ),




])



# Given chosen modelname, change path for phases-menu
@app.callback(
    Output("modelname", "children"),
    Input("model-dropdown", "value"),
)
def create_phase_dict(modelname):
    return modelname








if __name__ == "__main__":
    app.run_server(debug=True)




'''
symbol_dict = dict(AAPL="Apple", NVDA="Nvidia", TSLA="Tesla", IBM="IBM")

stock_options_dropdown = [{
    "label": name, 
    "value": symbol} 
    for symbol, name in symbol_dict.items()]

df_dict = {symbol: stock_data_object.stock_dataframe(symbol)
    for symbol in symbol_dict}

# OHLC options : ie open high low close-options
# This choses which value to plot, the close time value or open time value etc.
# The columns in the data files are:
# [timestamp,open,high,low,close,volume]
# In the plot, we change which column is plotted with y=ohlc
ohlc_options = [{"label": option.capitalize(), "value": option} 
                for option in ["open", "high", "low", "close"]]

# Settings for time slider, slider marks must be a dict
# value=2 sets default value, which here is 1month
slider_marks = {i: mark for i,mark in enumerate(
    ["1 day", "1 week", "1 month", "3 months", "1 year", "5 years", "Max"]
)}


# Dash board settings
stylesheets = [dbc.themes.MATERIA]

# Create a dash object called app
app = dash.Dash(__name__, external_stylesheets=stylesheets,
    meta_tags=[dict(name="viewport", content="width=device-width, initial-scale=1.0")]
)

# Needed for Heroku to connect to
server = app.server

# Define the layout for our "app"
# html.H1 = header 1, main title
# html.H2 = header 2, sub title
# html.P = paragraph
# dropdown creates a dropdown manue, id is it's name (for later), value is
# its default value. Options is a list of what you can chose.
# So value='' gives the default value for each dcc.
# The objects in app.layout ends up on the page in the order we have them written here
app.layout = dbc.Container([

    # dbc-cards are like boxes in the page
    # with classname we can change margins (mt: margintop, m:margin all around), 
    # and style (text-primary etc)
    # lg and xl is to change length and size of boxes

    # Title box
    dbc.Card([
        dbc.CardBody(html.H1("Stocks viewer", 
            className='text-primary m-3')
        )
    ], className='mt-3'),

    # Different rows and columns for the page
    dbc.Row([
        dbc.Col(
            html.P("Choose a stock"), 
            className='mt-1',
            xl={'size':1, 'offset':2}
        ),
        dbc.Col(
            dcc.Dropdown(
                id='stock-picker-dropdown', className='',
                options=stock_options_dropdown,
                value='AAPL'
            ),
            lg='4', xl='3'
        ),
        dbc.Col(
            dbc.Card(
                dcc.RadioItems(id='ohlc-radio', className='',
                options=ohlc_options,
                value='close')
            ),
            lg='4', xl='3'
        )
    ], className='mt-4'),

    dbc.Row([
        dbc.Col([

            dcc.Graph(
                id='stock-graph', className=''
            ),
    
            dcc.Slider(
                id='time-slider', className='',
                min=0, max=6,
                step=None,
                value=2,
                marks=slider_marks
            )], 
            lg={"size":"6", "offset":1}, xl='6'
        ),
        # Interactive text below dropdown menu
        # we use header2 (because of order), and change font size under classname
        # to h5
        # text-dange/success choses colour
        dbc.Col([
            dbc.Card([
                html.H2("Highest value", className='h5 mt-3 mx-3'),
                html.P(id = "highest-value", className='mx-3 h1 text-success')
            ], className='mt-5 w-50'),
            dbc.Card([
                html.H2("Lowest value", className='h5 mt-3 mx-3'),
                html.P(id = "lowest-value", className='mx-3 h1 text-danger'),
            ], className='w-50')
        ])
    ]),


    # Stores an intermediate value on client's browser for sharing between callbacks
    # Needs a json object, see filter_df() below
    dcc.Store(id='filtered-df')
], fluid=True)

# Use callback decorators for extra functions

# This stores our time filtered columns in a temporary json-file
# This way we can use these data for several callback.
# In the previous version we filtered the data inside the plot function
# then we could only use it directly for the plot function
# Now we use the filtered data both for plot and for giving max/min values
@app.callback(Output("filtered-df", "data"), 
    Input("stock-picker-dropdown", "value"),
    Input("time-slider", "value")
)
def filter_df(stock, time_index):
    """Filters the dataframe and stores in intermediary for callbacks
    Returns:
        json object of filtered dataframe
    """
    # Extract stocks (see df_dict above, it uses class to extract data from csvfiles)
    dff_daily, dff_intraday = df_dict[stock]

    # intraday or daily depends on time index input (value)
    dff = dff_intraday if time_index <= 2 else dff_daily

    # Define a dictionary with days, maps 0-6 to number of days, max is not necessary
    days = {i: day for i,day in enumerate([1,7,30,90,365,5*365])}

    # This gives max instead
    dff = dff if time_index == 6 else filter_time(dff, days[time_index])

    return dff.to_json()


# Input 1 is our dropdown menu
# Output is the graph
# Input 2 is the time slider
@app.callback(
    Output("stock-graph", "figure"),
    Input("filtered-df", "data"),
    Input("stock-picker-dropdown", "value"),
    Input("ohlc-radio", "value")
)
def update_graph(json_df, stock, ohlc):
    
    # Change our temporary json file to dataframe
    dff = pd.read_json(json_df)

    # Plot figure using dff-data
    fig = px.line(dff, 
        x=dff.index, y=ohlc, 
        labels={"close":"Stock value (us-$)"}, 
        title=symbol_dict[stock]
    )

    return fig # fig object goes into output property, ie figure property


# This gives the number of the maximum and lowest value shown in the figure
@app.callback(
    Output("highest-value", "children"),
    Output("lowest-value", "children"),
    Input("filtered-df", "data"),
    Input("ohlc-radio", "value")
)
def highest_lowest_value(json_df, ohlc):

    dff = pd.read_json(json_df)

    highest_value = f"{dff[ohlc].max():.1f} $"
    lowest_value = f"{dff[ohlc].min():.1f} $"
    
    return highest_value, lowest_value

if __name__ == "__main__":
    app.run_server(debug=True)
'''





