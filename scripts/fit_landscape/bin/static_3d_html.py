#!/usr/bin/env python
# coding: utf-8
# Generate static 3D models of user-defined XYZ axes for HTML output/interaction
# Author: Brian Mann
# Last Modified: 2024-05-01

### IMPORT PYTHON PACKAGES ###

import re
import argparse

import pandas as pd
import plotly.express as px


### ASSIGN GLOBAL VARIABLES ###

infile = None
html = None
delimiter = "tab"
xfactor = None
yfactor = None
zfactor = None


### INPUT VARIABLES ###

parser = argparse.ArgumentParser()
parser.add_argument(
    "-i",
    "--infile",
    required=True,
    help="Provide filename for input data matrix"
)
parser.add_argument(
    "-d",
    "--delimiter",
    required=False,
    help="Indicate delimiter (and actual character) for provided data file. [Options: comma, tab, and space]"
)
parser.add_argument(
    "-v",
    "--html",
    required=False,
    help="Provide filename for HTML-formatted static 3D landscape output model"
)
parser.add_argument(
    "-x",
    "--xfactor",
    required=False,
    help="Provide specific factor/variable on the 3D model X-coordinate axis [Example: hd or pcd]"
)
parser.add_argument(
    "-y",
    "--yfactor",
    required=False,
    help="Provide specific factor/variable on the 3D model Y-coordinate axis [Example: SD-related factor]"
)
parser.add_argument(
    "-z",
    "--zfactor",
    required=False,
    help="Provide specific factor/variable on the 3D model Y-coordinate axis [Example: SD-related factor]"
)

args = parser.parse_args()

infile = args.infile
html = args.html
delimiter = args.delimiter
xfactor = args.xfactor
yfactor = args.yfactor
zfactor = args.zfactor


### IMPORT DATA ###

if args.delimiter is None:
    infile = pd.read_csv(infile, sep='\t')
elif args.delimiter == "comma":
    infile = pd.read_csv(infile, sep=",")
elif args.delimiter == "tab":
    infile = pd.read_csv(infile, sep="\t")
elif args.delimiter == "space":
    infile = pd.read_csv(infile, sep=" ")

# Import "core" data model
core = pd.DataFrame()
if str(infile['subtype'].unique()[0]) == "B vic":
    core = pd.read_csv('tbl/B_iconsolidated.txt', sep='\t')
elif str(infile['subtype'].unique()[0]) == "H1 swl":
    core = pd.read_csv('tbl/H1_iconsolidated.txt', sep='\t')
elif str(infile['subtype'].unique()[0]) == "H3":
    core = pd.read_csv('tbl/H3_iconsolidated.txt', sep='\t')


### DEFINE X, Y, AND Z FACTORS ###

processed = infile.copy()

# Define all possible, unique X "factors"
if args.xfactor is None:
    x_unique = list(
        pd.unique( 
            list(
                processed['x_factor']
            )
        )
    )
    x_selection = x_unique[0]
else:
    x_selection = args.xfactor

# Define all possible, unique Y "factors"
if args.yfactor is None:
    y_unique = list(
        pd.unique( 
            list(
                processed['y_factor']
            )
        )
    )
    y_selection = y_unique[0]
else:
    y_selection = args.yfactor

# Define all possible, unique Z "factors"
if args.zfactor is None:
    z_unique = list(
        pd.unique( 
            list(
                processed['z_factor']
            )
        )
    )
    z_selection = z_unique[0]
else:
    z_selection = args.zfactor


### SETUP 3D LANDSCAPE ###

# Define labels for X, Y, and Z axes
if x_selection == "hd":
    xlabel = 'Hamming Distance (HD)'
elif x_selection == "pcd":
    xlabel = 'Physiochemical Distance (PCD)'

temporal_context = re.findall(
    r'.*(seasonal|reference|srel|rrel).*',
    y_selection
)

atomic_context = re.findall(
    r'.*(calpha|carbon|heavy).*',
    y_selection
)

if 'sd' in y_selection:
    atomic_context = re.sub( 'Calpha', 'Ca', str(atomic_context[0]).capitalize() )
    ylabel = 'Structural Distance [SD, ' + str(temporal_context[0]).capitalize() + " (" + atomic_context + ")]"
elif 'total_score' in y_selection:
    if str(temporal_context[0]) == 'srel':
        ylabel = 'Rosetta Total Score [Relative, Seasonal]'
    elif str(temporal_context[0]) == 'rrel':
        ylabel = 'Rosetta Total Score [Relative, Reference]'
    else:
        ylabel = 'Rosetta Total Score'

testing_context = str(
    re.findall(
        r'.*(logfold|log|total\_score).*',
        z_selection
    )[0]
)

if "cladei|meani" in z_selection:
    aggregation_label = "(Clade/Mean)"
elif "cladei" in z_selection or "meani" in z_selection:
    aggregation_context = [
        str(
            re.findall(
                r'.*(clade|mean).*',
                z_selection
            )[0]
        ).capitalize()
    ]
    aggregation_label = "(" + str(aggregation_context[0]) + ")"
else:
    aggregation_label = ''

reagent_context = re.findall(
    r'.*(\d{4}\-\d{3}).*',
    z_selection
)

if len( reagent_context ) > 1:
    reagent_list = ",".str_join(reagent_context)
elif len( reagent_context ) == 1:
    reagent_list = str(reagent_context[0])
else:
    reagent_list = str()

if "log" in testing_context:
    antiserum_context = str(
        re.findall(
            r'.*(ha\_protein|lot).*',
            z_selection
        )[0]
    ).upper().replace('_',' ')
    zlabel = antiserum_context + ": " + reagent_list + " [" + testing_context.upper() + "] " + aggregation_label
elif "total_score" in testing_context:
    relative_context = str(
        re.findall(
            r'.*(rrel|srel).*',
            z_selection
        )
    )
    if str(relative_context[0]) == 'srel':
        zlabel = 'Rosetta Total Score [Relative, Seasonal]'
    elif str(relative_context[0]) == 'rrel':
        zlabel = 'Rosetta Total Score [Relative, Reference]'
    else:
        zlabel = 'Rosetta Total Score'

# Sub-Select surface 3D mesh model from user-provided input file
model = processed[ processed['x_factor'] == x_selection ]
model = model[ model['y_factor'] == y_selection ]
model = model[ model['z_factor'] == z_selection ]

# Model "landscape" and define all axes ranges/labels
fig = px.scatter_3d(
    x = model['x_predicted'],
    y = model['y_predicted'],
    z = model['z_predicted'],
    opacity = 0.0,
    labels = {
        'x': xlabel,
        'y': ylabel,
        'z': zlabel
    }
)

# Define "absolute" data points based on non-NULL XYZ coordinates
column_subset = [ 'variant_hash', x_selection, y_selection, z_selection ]
absolute_model = core[ column_subset ].dropna()

fig.add_scatter3d(
    x = absolute_model[ x_selection ],
    y = absolute_model[ y_selection ],
    z = absolute_model[ z_selection ],
    mode = 'markers',
    marker = dict(
        size = 5,
        color = '#A9A9A9'
    ),
    showlegend = False,
    text = absolute_model['variant_hash'],
    hovertemplate = "<b>%{text}</b><br>"
        + "<b>Coordinates:</b> (%{x}, %{y}, %{z})"
        + "<extra></extra>"
)

# Define surface 3D mesh model from landscape data set
model_pivot = pd.pivot_table(
    model,
    values = 'z_predicted',
    index = ['y_predicted'],
    columns = ['x_predicted']
)

fig.add_surface(
    x = list(model_pivot.columns),
    y = model_pivot.index,
    z = model_pivot.values,
    colorscale = 'balance_r',
    cmid = -2,
    opacity = 0.8,
    showlegend = False,
    hovertemplate = "<b>X: </b>%{x}<br>"
        + "<b>Y: </b>%{y}<br>"
        + "<b>Z: </b>%{z}<br>"
        + "<extra></extra>"
)

fig.update_layout(
    scene = dict(
        xaxis = dict( nticks = 8 ),
        yaxis = dict( nticks = 8 )
    ),
    font_family = "Arial",
    font_color = 'black'
)

fig.show()
