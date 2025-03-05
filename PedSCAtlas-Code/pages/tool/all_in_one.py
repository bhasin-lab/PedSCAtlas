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

#from load_data import DataHealthy, DataHealthyPD, DataAcuteLeukemia, DataLambo
from load_data import DataAcuteLeukemia, DataLambo, DataHealthyPD, DataHealthy

# load helper functions
from helper.survival import *
from helper.get_prop import * 
from helper.utils import id_factory

# import scatter plot functions from helper/scatter.py
from helper.scatter import scatter_lambo, scatter_healthyPD, scatter_pan_leukemia

# import expression plot functions from helper/sc_exp.py
from helper.sc_exp import *

from helper.plot_options import *

# page prefix, same across ONLINE + LOCAL
page_pre = '/'

# button prefix
if MODE=="ONLINE":
    button_prefix="/PediatricSCAtlas"
elif MODE=="LOCAL":
    button_prefix=""

dash.register_page(__name__, path = page_pre + 'tool/overview', title='PedSCAtlas', image='PedSCAtlas_logo.png')

pageid=id_factory("Overview")

# input for coloring cells on UMAP
#sc_scatter_inputs = html.Div([
#    dbc.Row([
#        dbc.Label("Color Cells By"),
#        dcc.Dropdown(
#            id='sc_pt_color_o',
#            options = ['Cell Type','Sample Id','Cluster'],
#           value = "Cell Type"
#        ),
#        html.Div("Uniform Manifold Approximation and Projection (UMAP), where each point represents a cell in the dataset.")
#    ]),
#])

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

# inputs for sc gene expression
sc_gene_inputs = html.Div([
    dbc.Row([
        dbc.Label("Select Gene:"),
        dcc.Dropdown(
            id = 'sc_gene',
            options = DataHealthyPD.gene_names,
            value = "PTEN"
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
        html.Div("Select a gene name and plot type to show its expression in these Single-Cell RNA datasets.")
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

def scatterplot_template(title, desc, graph_id, d_options, sc_pt_color_id, page_link):

    scatterlayout = dbc.Col([
        dbc.Row([ # row for description
            dbc.Container([
                #html.H3(title),
                dbc.Button(title, color='primary', href=page_link, size='lg'),
                html.Hr(),
                html.Div(children=dcc.Markdown(desc), style = {"height":"20ex"}), #dcc.Markdown(DataHealthyPD.dst_desc)
            ])
        ]),
        dbc.Row([ # row for scatter inputs
            dbc.Container([
                html.Hr(),
                html.H5("Choose a category to color cells by:"),
                dcc.Dropdown(
                    id=sc_pt_color_id,
                    options = d_options,
                    value='Cell Type'
                )
            ])
        ]),
        dbc.Row([ # row for scatter
            dbc.Spinner([dcc.Graph(id=pageid(graph_id), figure={'layout':{
                            'autosize':True,
                            'height':500}})])
        ])
    ], width = 4)
    return scatterlayout

layout = html.Div([

    html.Br(),
    dbc.Row([
        dbc.Col(width=1),
        dbc.Col([ #
            dbc.Row(dbc.Container([# row for header and selection
                html.H2("Metadata UMAPs"),
                html.Hr(),
                #dbc.Row([
                #    dbc.Col(sc_scatter_in, width=6)
                #]),
                html.Br(),
                dbc.Row([ # row for aligned figures 
                    scatterplot_template("Pediatric Leukemia", DataAcuteLeukemia.dst_desc, "sc_scatter_pan", options_pan, "sc_opt_pan", button_prefix + 'tool/pediatric_leukemia'),
                    scatterplot_template("AML", DataLambo.dst_desc, "sc_scatter_lambo", options_lambo, "sc_opt_lambo", button_prefix + 'tool/aml_lambo_2023'),
                    scatterplot_template("Healthy Pediatric", DataHealthyPD.dst_desc, "sc_scatter_hb", options_hb, "sc_opt_hb", button_prefix + 'tool/healthy_pediatric'),
                ])
            ], fluid = True))
        ], width = 10)
    ]),
    html.Br(),

    # second row - sample cell type distributions
    dbc.Row([
        dbc.Col(width = 1),

        dbc.Col(
            dbc.Container([
                html.H2("Cell Type Proportions"),
                html.Hr(),
                dbc.Row([
                    dbc.Col(prop_inputs, width=2),
                    dbc.Col([
                        dcc.Graph(id=pageid("ct_prop_pan")),
                        dcc.Graph(id=pageid("ct_prop_lambo")),
                        dcc.Graph(id=pageid("ct_prop_hb"))
                    ], width = 10)
                ])
            ], fluid = True)
        , width = 10)
    ]),

    # third row - single-cell gene expression - cancer and healthy
    dbc.Row([
        dbc.Col(width=1),

        dbc.Col(
            dbc.Container([
                html.H2("Single-Cell Gene Expression in Datasets"),
                html.Hr(),
                dbc.Row([
                    dbc.Col(sc_gene_inputs, width = 2),
                    dbc.Col([
                        dbc.Spinner([dcc.Graph(id=pageid("sc_exp_fig_pan"))]),
                        dbc.Spinner([dcc.Graph(id=pageid("sc_exp_fig_lambo"))]),
                        dbc.Spinner([dcc.Graph(id=pageid("sc_exp_fig_hb"))])
                    ], width = 10)
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
                    dbc.Col(dbc.Spinner([dcc.Graph(id=pageid("healthy_exp_fig"))]), width = 10)
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


### UMAPs ###
# sc UMAP - Pan-Leukemia
@callback(
    Output(pageid("sc_scatter_pan"), 'figure'),
    Input('sc_opt_pan','value')
)
def update_sc_scatter_pan(sc_opt_pan):        
    sc_scatter_subtype = DataAcuteLeukemia.df.Diagnosis.unique()
    return scatter_pan_leukemia(sc_scatter_subtype, sc_opt_pan, page='overview')

# sc UMAP - Lambo AML
@callback(
    Output(pageid("sc_scatter_lambo"), 'figure'),
    Input('sc_opt_lambo','value')
)
def update_sc_scatter_lambo(sc_opt_lambo):        
    return scatter_lambo(sc_opt_lambo, page='overview')

# sc UMAP - Healthy Ped BM
@callback(
        Output(pageid("sc_scatter_hb"),'figure'),
        Input('sc_opt_hb','value')
)
def update_sc_scatter_hb(sc_opt_hb):
    return scatter_healthyPD(sc_opt_hb, page='overview')


### cell type proportions ###
@callback(
        Output(pageid("ct_prop_pan"),'figure'),
        Output(pageid("ct_prop_lambo"),'figure'),
        Output(pageid("ct_prop_hb"),'figure'),
        Input('prop_type','value'),
        Input('prop_method','value')
)
def update_ct_prop(prop_type, prop_method):
    prop_fig_pan = get_prop(prop_type, DataAcuteLeukemia.df, prop_method, DataAcuteLeukemia.sc_colors.Color) # , x_title=False
    prop_fig_lambo = get_prop(prop_type, DataLambo.df, prop_method, DataLambo.sc_colors.Color, title=False) #, x_title=False
    prop_fig_hb = get_prop(prop_type, DataHealthyPD.df, prop_method, DataHealthyPD.sc_colors.Color, title=False)
    return prop_fig_pan, prop_fig_lambo, prop_fig_hb


### sc Gene Exp - datasets ###
@callback(
        Output(pageid("sc_exp_fig_pan"), 'figure'),
        Output(pageid("sc_exp_fig_lambo"), 'figure'),
        Output(pageid("sc_exp_fig_hb"), 'figure'),
        Input('sc_gene','value'),
        Input('sc_gene_type','value'),
        Input('sc_group_by','value'),
        #Input('sc_exp_subtype','value')
)
def update_sc_exp(sc_gene, sc_gene_type, sc_group_by):
    sc_exp_fig_pan = load_scExp_pan(sc_gene, sc_gene_type, sc_group_by, DataAcuteLeukemia.df.Diagnosis.unique())
    sc_exp_fig_lambo = load_scExp_lambo(sc_gene, sc_gene_type, sc_group_by)
    sc_exp_fig_hb = load_scExp_hb(sc_gene, sc_gene_type, sc_group_by)
    return sc_exp_fig_pan, sc_exp_fig_lambo, sc_exp_fig_hb


### sc Gene Exp - healthy adult ###
@callback(
        Output(pageid("healthy_exp_fig"), 'figure'),
        Input('healthy_gene','value'),
        Input('healthy_gene_type','value'),
        Input('healthy_celltype','value')
)
def update_healthy_exp(healthy_gene, healthy_gene_type, healthy_celltype):
    return load_scExp_h_adult(healthy_gene, healthy_gene_type, healthy_celltype)
