import dash
from dash import html, dcc, callback, Input, Output, State, dash_table
import numpy as np
import pandas as pd
import dash_bootstrap_components as dbc
from datetime import date
import plotly.express as px
import plotly.graph_objects as go
from config import MODE

# file header - establish location of "Data" directory
if MODE=="ONLINE":
    dir = '/srv/dash-apps/PedSCAtlas-Data/'
elif MODE=="LOCAL":
    dir = '/opt/localdata/hmumme/PedSCAtlas/Data/'

# page prefix, same across ONLINE + LOCAL
page_pre = '/'

# load de tables
sc_de = pd.read_csv(dir + "AcuteLeukemia/sc_de.csv", index_col=0)
bulk_de = pd.read_csv(dir + "AcuteLeukemia/bulk_de.csv", index_col=0)
info = pd.read_csv(dir + "AcuteLeukemia/degs_info.csv", index_col=0)

dash.register_page(__name__, path=page_pre + 'DE', title='PedSCAtlas', image='PedSCAtlas_logo.png')

#dash.register_page(__name__, path='/PediatricSCAtlas-Dev/DE', title='PedSCAtlas', image='PedSCAtlas_logo.png')


PAGE_SIZE = 12

layout = html.Div(children=[

    html.Br(),

    dbc.Row([
        dbc.Col(width=1),
        dbc.Col(
            dbc.Container([
                html.H2("Choose a Disease Type"),
                html.Hr(),
                dcc.Markdown('''
                    Differential Expression (DE) testing was performed comparing each disease type in the **Pan-Leukemia** study against a healthy (or normal) control
                    in the Single-Cell and Bulk RNA-seq datasets. The exact parameters used are detailed on our [Github](https://github.com/bhasin-lab/PedSCAtlas/). 
                    Pick a disease type to explore its DE results.
                '''),
                dcc.Dropdown(
                    id='type',
                    options=np.unique(info.Disease_Type.values),
                    value=info.Disease_Type.values[0],
                    ),
                    html.Br(),
                html.P("Select log2FC threshold:"),
                dcc.Slider(min=0.25, max=2, step=0.25, value=0.25,
                    id='th_FC',
                    ),
                html.P("Select negative log10 of adjusted p-value (NLP) threshold:"),
                dcc.Slider(min=1, max=5, step=0.5,
                    id='th_NLP', value=1.5,
                    ),
                html.Br(),
                html.Div([
                    dcc.Markdown('''
                        Number of DEGs, samples, and cells for chosen disease type and thresholds.''')
                ]),
                dash_table.DataTable(id="type_info", data = info.to_dict('records'), columns = [{"name": i, "id": i} for i in info.columns])
            ]), width = 5, style={}
        ),

        dbc.Col(
            dbc.Container([
                html.H2("DEGs"),
                html.Hr(),
                html.Div([
                    dcc.Markdown('''Select a RNA Type to view DEGs in a table, based on thresholds chosen to the left.''')
                ]),
                dcc.RadioItems(['SC', 'Bulk'], 'SC', id = 'type_deg', 
                               inline = False, labelStyle={'display': 'block'}, inputStyle={"margin-right": "5px"}),
                html.Br(),                               
                dash_table.DataTable(id='table_deg', 
                                     columns = [{"name": i, "id": i} for i in sc_de.columns],
                                     page_current=0,
                                     page_size=PAGE_SIZE,
                                     page_action='custom',
                                     sort_action='custom',
                                     sort_mode='single',
                                     sort_by=[])
            ])
        ),

        dbc.Col(width=1),
    ]),

    html.Br(),

    dbc.Row([
        dbc.Col(width=1),
        dbc.Col(dbc.Container([
            html.H2("Volcano Plot"),
            html.Hr(),
            html.Div([
                dcc.Markdown('''Select a RNA Type to view DEGs on a volcano plot, based on thresholds chosen above. Note: many adjusted p-values were reported as zero in both SC and Bulk DEGs. For these values, the smallest non-zero p-value reported was substituted when calculating the negative log10 p-value (NLP).'''),
                dcc.RadioItems(['SC', 'Bulk'], 'SC', id = 'type_vol', 
                               inline = False, labelStyle={'display': 'block'}, inputStyle={"margin-right": "5px"}),
            ]),
        ]), width=3),
        dbc.Col(dbc.Container([
            dcc.Graph('volcano')
        ])),
        dbc.Col(width=1),
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
            dbc.Col(width=1),
        ])
    )]),
    
])

# Calculate number of DEGs and format for "Disase Type Information" table
@callback(
    Output('type_info', 'data'),
    Input('type', 'value'),
    Input('th_FC', 'value'),
    Input('th_NLP', 'value'))
def update_info(type, th_FC, th_NLP):
    # filter to correct type
    temp = info[info['Disease_Type'] == type].copy()


    # assign UP, DOWN, NS classes
    sc_de["SIG"] = "NS"
    sc_de.loc[(sc_de.NLP >th_NLP) & (sc_de.log2FC > th_FC), "SIG"] = "UP"
    sc_de.loc[(sc_de.NLP >th_NLP) & (sc_de.log2FC < th_FC*-1), "SIG"] = "DOWN"

    bulk_de["SIG"] = "NS"
    bulk_de.loc[(bulk_de.NLP >th_NLP) & (bulk_de.log2FC > th_FC), "SIG"] = "UP"
    bulk_de.loc[(bulk_de.NLP >th_NLP) & (bulk_de.log2FC < th_FC*-1), "SIG"] = "DOWN"

    # filter to correct disease type
    if type == "B-ALL":
        bulk_de_sub = bulk_de[bulk_de["Type"] == "Pre B-ALL"]
    else:
        bulk_de_sub = bulk_de[bulk_de["type"] == type]
    
    sc_de_sub = sc_de[sc_de["type"] == type + " Blast"]


    # calculate DEGs based on thresholds
    temp.loc[temp['RNA_Type'] == "SC","n_DEGs"] = len(np.where(sc_de_sub.SIG != "NS")[0])
    temp.loc[temp['RNA_Type'] == "Bulk","n_DEGs"] = len(np.where(bulk_de_sub.SIG != "NS")[0])
    return temp.to_dict("records")

# Show DEGs table in "DEGs" section
@callback(
    Output('table_deg', 'data'),
    Input('table_deg', 'page_current'),
    Input('table_deg', 'page_size'),
    Input('table_deg', 'sort_by'),
    Input('type', 'value'),
    Input('th_FC', 'value'),
    Input('th_NLP', 'value'),
    Input('type_deg', 'value')
)
def degs_table(page_current, page_size, sort_by, type, th_FC, th_NLP, type_deg):
    # filter to correct disease type
    if type_deg == "SC":
        type_name = type + " Blast"
        temp = sc_de.loc[sc_de["type"] == type_name].copy()
    else:
        temp = bulk_de.loc[bulk_de["type"] == type].copy()

    # assign UP, DOWN, NS classes
    temp["SIG"] = "NS"
    temp.loc[(temp.NLP >th_NLP) & (temp.log2FC > th_FC), "SIG"] = "UP"
    temp.loc[(temp.NLP >th_NLP) & (temp.log2FC < th_FC*-1), "SIG"] = "DOWN"
    
    # filter to not include NS
    temp = temp.loc[temp["SIG"] != "NS"]

    # custom sorting
    if len(sort_by):
        dff = temp.sort_values(
            sort_by[0]['column_id'],
            ascending=sort_by[0]['direction'] == 'asc',
            inplace=False
        )
    else:
        dff = temp
    
    # change format
    dff['p_val_adj'] = dff['p_val_adj'].map('{:e}'.format)


    return dff.round(4).iloc[
        page_current*page_size:(page_current+ 1)*page_size
    ].to_dict('records')

# Volcano Plot
@callback(
    Output('volcano', 'figure'),
    Input('type', 'value'),
    Input('th_FC', 'value'),
    Input('th_NLP', 'value'),
    Input('type_vol', 'value'))
def volcano(type, th_FC, th_NLP, type_vol):
    # filter to correct disease type
    if type_vol == "SC":
        type_name = type + " Blast"
        temp = sc_de.loc[sc_de["type"] == type_name].copy()
    else:
        type_name = type
        temp = bulk_de.loc[bulk_de["type"] == type].copy()

    # assign UP, DOWN, NS classes
    temp["SIG"] = "NS"
    temp.loc[(temp.NLP >th_NLP) & (temp.log2FC > th_FC), "SIG"] = "UP"
    temp.loc[(temp.NLP >th_NLP) & (temp.log2FC < th_FC*-1), "SIG"] = "DOWN"

    # scientific notation
    temp['p_val_adj'] = temp['p_val_adj'].map('{:e}'.format)

    fig=px.scatter(temp, x="log2FC", y='NLP', color='SIG',
        color_discrete_map={"NS":'grey', "DOWN":'blue', "UP":'red'},
        hover_name="gene", hover_data = ['p_val_adj', 'log2FC', 'NLP'],
        template='plotly_white',
        title=type_name + " DEGs - " + type_vol)

    fig.add_shape(type='line', 
        x0=min(temp.log2FC)-0.1, 
        x1=max(temp.log2FC)+0.1,
        y0=th_NLP, 
        y1=th_NLP,
        line = dict(color="black", width=2, dash="dot"),
        opacity=0.5)
    fig.add_shape(type='line', 
        x0=th_FC, 
        x1=th_FC,
        y0=min(temp.NLP), 
        y1=max(temp.NLP)+10,
        line = dict(color="black", width=2, dash="dot"),
        opacity=0.5)
    fig.add_shape(type='line', 
        x0=th_FC*-1, 
        x1=th_FC*-1,
        y0=min(temp.NLP), 
        y1=max(temp.NLP)+10,
        line = dict(color="black", width=2, dash="dot"),
        opacity=0.5)
    
    return fig