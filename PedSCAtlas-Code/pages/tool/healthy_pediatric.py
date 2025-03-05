import dash
from dash import Dash, html, dcc, Input, Output, State, dash_table, callback
import dash_bootstrap_components as dbc
import plotly.express as px
import plotly.graph_objects as go
import pandas as pd
import numpy as np
from itertools import cycle
import pickle
from pathlib import Path
import matplotlib.pyplot as plt
import base64
from io import BytesIO
from config import MODE

from load_data import DataHealthy, DataHealthyPD

# load helper functions
from helper.get_prop import * 
from helper.scatter import scatter_healthyPD
from helper.sc_exp import load_scExp_hb, load_scExp_h_adult
from helper.plot_options import options_hb

# page prefix, same across ONLINE + LOCAL
page_pre = '/'

dash.register_page(__name__, path = page_pre + 'tool/healthy_pediatric', title='PedSCAtlas', image='PedSCAtlas_logo.png')

# input for coloring cells on UMAP
sc_scatter_inputs = html.Div([
    dbc.Row([
        dbc.Label("Color Cells By"),
        dcc.Dropdown(
            id='sc_ind_hb',
            options = options_hb,
            value = "Cell Type"
        ),
    ]),
    html.Br(),
    html.Div("Uniform Manifold Approximation and Projection (UMAP), where each point represents a cell in the dataset.")
])

# input for grouping of Cell Type Proportions plot
prop_inputs = html.Div([
    dbc.Row([
        dbc.Label("Choose variable to show cell type proportions"),
        dcc.RadioItems(
            id = 'prop_type',
            options = ['Sample Id'],
            value = 'Sample Id', inline = False, labelStyle={'display': 'block'}, inputStyle={"margin-right": "5px"}
        )
    ]),
    html.Br(),
    dbc.Row([
        dbc.Label('Choose whether to show proportions or actual counts'),
        dcc.RadioItems(
            id = 'prop_method',
            options = ['Proportion','Actual'],
            value = 'Proportion', inline = False, labelStyle={'display': 'block'}, inputStyle={"margin-right": "5px"}
        )
    ])
]),

# inputs for sc gene expression - pediatric healthy
sc_gene_inputs = html.Div([
    dbc.Row([
        dbc.Label("Select Gene:"),
        dcc.Dropdown(
            id = 'sc_gene',
            options = DataHealthyPD.gene_names,
            value = "CD3D"
        )
    ]),
    html.Br(),
    dbc.Row([
        dbc.Col([
            dbc.Label("Choose Plot Type:"),
            dcc.RadioItems(
                id = 'sc_gene_type',
                options = ["Violin","Box-Plot","Feature"],
                value = "Violin", inline = False, labelStyle={'display': 'block'}, inputStyle={"margin-right": "5px"}
            ),
        ])
    ]),
    html.Br(),
    dbc.Row([
        dbc.Label("Group Cells By:"),
        dcc.Dropdown(
            id = 'sc_group_by',
            options = ['Cell Type','Sample Id','Cluster'],
            value = 'Cell Type' 
        )
    ]),
    html.Br(),
    dbc.Row([
        html.Div("Select a gene name and plot type to show its expression in Single-Cell RNA Healthy Pediatric dataset.")
    ])
])

# inputs for sc gene expression - healthy
healthy_gene_inputs = html.Div([
    dbc.Row([
        dbc.Label("Select Gene:"),
        dcc.Dropdown(
            id = 'healthy_gene',
            options = DataHealthy.gene_names,
            value = "CD3D"
        )
    ]),
    html.Br(),
    dbc.Row([
        dbc.Col([
            dbc.Label("Choose Plot Type:"),
            dcc.RadioItems(
                id = 'healthy_gene_type',
                options = ["Violin","Box-Plot","Feature"],
                value = "Violin", inline = False, labelStyle={'display': 'block'}, inputStyle={"margin-right": "5px"}
            )
        ])
    ]),
    html.Br(),

    dbc.Row([
        dbc.Col([
            dbc.Label("Choose Level:"),
            dcc.RadioItems(
                id = 'healthy_celltype',
                options = ['CellType_Broad','CellType_Fine'],
                value = 'CellType_Broad', inline = False, labelStyle={'display': 'block'}, inputStyle={"margin-right": "5px"}
            )
        ])
    ]),

    html.Br(),

    dbc.Row([
        html.Div("Select a gene name and plot type to show its expression in Single-Cell RNA Healthy Adult dataset.")
    ])
])


layout = html.Div([

    html.Br(),

    # first row - Dataset Description and SC UMAP
    dbc.Row([
        dbc.Col(width=1),

        dbc.Col(
            dbc.Container([
                html.H2("Healthy Pediatric"),
                html.Hr(),
                html.Div(children=dcc.Markdown(DataHealthyPD.dst_desc))#,
            ]), width = 3, style={}
        ),

        dbc.Col(
            dbc.Container([
                html.H2("Metadata UMAP"),
                html.Hr(),
                dbc.Row([
                    dbc.Col(sc_scatter_inputs, width=3),
                    dbc.Col(dbc.Spinner([dcc.Graph(id='sc_scatter_hb_ind')]), width = 9,)
                ])
            ]), width = 7, style={}
        )
    ]),

    # second row - sample cell type distributions
    dbc.Row([
        dbc.Col(width = 1),

        dbc.Col(
            dbc.Container([
                html.H2("Cell Type Proportions"),
                html.Hr(),
                dbc.Row([
                    dbc.Col(prop_inputs, width=3),
                    dbc.Col(dcc.Graph(id='ct_prop_hb_ind'), width = 9)
                ])
            ], fluid = True)
        , width = 10)
    ]),

    # third row - single-cell gene expression - cancer and healthy
    dbc.Row([
        dbc.Col(width=1),

        dbc.Col(
            dbc.Container([
                html.H2("Single-Cell Gene Expression in Healthy Pediatric BM"),
                html.Hr(),
                dbc.Row([
                    dbc.Col(sc_gene_inputs, width = 2),
                    dbc.Col(dbc.Spinner([dcc.Graph(id='sc_exp_fig_hb_ind')]), width = 10)
                ])
            ], fluid = True), width = 10
        )
    ]),

    html.Br(),

    # fourth row - single-cell gene expression - healthy tissue
    dbc.Row([
        dbc.Col(width=1),

        dbc.Col(
            dbc.Container([
                html.H2("Single-Cell Gene Expression in Healthy Adult BM"),
                html.Hr(),
                dbc.Row([
                    dbc.Col(healthy_gene_inputs, width = 2),
                    dbc.Col(dbc.Spinner([dcc.Graph(id='healthy_exp_fig_hb')]), width = 10)
                ])
            ], fluid = True), width = 10
        )
    ]),

    html.Br(),

    html.Hr(),

    # Footer
    html.Center([html.Footer(
        dbc.Row([
            dbc.Col(width=1),
            dbc.Col([
                html.Div([
                    dcc.Markdown('''
                        Please cite our abstract published in 
                        [blood](https://ashpublications.org/blood/article/140/Supplement%201/2278/489360/A-Single-Cell-Atlas-and-Interactive-Web-Resource), 
                        presented during the ASH 2022 conference. Manuscript is under review.''')
                ])
            ]),
            dbc.Col([
                dcc.Markdown('''
                    Developed by [Bhasin Systems Biomedicine Lab](https://bhasinlab.org/) at Emory University
                    in collaboration with the [Aflac Cancer & Blood Disorders Center](https://www.choa.org/medical-services/cancer-and-blood-disorders) at 
                    Children's Healthcare of Atlanta.
                ''')
            ]),
            dbc.Col(width=1)
        ])
    )]),

    html.Hr()
    

])


# sc UMAP
@callback(
        Output('sc_scatter_hb_ind','figure'),
        Input('sc_ind_hb','value')
)
def update_sc_figures(sc_ind_hb):        
    return scatter_healthyPD(sc_ind_hb,'ind')

# cell type proportions
@callback(
        Output('ct_prop_hb_ind','figure'),
        Input('prop_type','value'),
        Input('prop_method','value')
)
def update_ct_prop(prop_type, prop_method):
    prop_fig = get_prop(prop_type, DataHealthyPD.df, prop_method, DataHealthyPD.sc_colors.Color)
    return prop_fig


# sc Gene Exp - cancer
@callback(
        Output('sc_exp_fig_hb_ind', 'figure'),
        Input('sc_gene','value'),
        Input('sc_gene_type','value'),
        Input('sc_group_by','value'),
)
def update_sc_exp(sc_gene, sc_gene_type, sc_group_by):
    return load_scExp_hb(sc_gene, sc_gene_type, sc_group_by)

# sc Gene Exp - healthy
@callback(
        Output('healthy_exp_fig_hb', 'figure'),
        Input('healthy_gene','value'),
        Input('healthy_gene_type','value'),
        Input('healthy_celltype','value')
)
def update_healthy_exp(healthy_gene, healthy_gene_type, healthy_celltype):
    return load_scExp_h_adult(healthy_gene, healthy_gene_type, healthy_celltype)