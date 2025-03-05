import dash
from dash import html, dcc, callback, Input, Output, State
import pandas as pd
import dash_bootstrap_components as dbc
from datetime import date
from config import MODE

# page prefix
if MODE=="ONLINE":
    prefix="/PediatricSCAtlas"
elif MODE=="LOCAL":
    prefix=""

dash.register_page(__name__, path='/', title='PedSCAtlas', image='PedSCAtlas_logo.png')


layout = html.Div(children=[
    html.Br(),html.Br(),

    # intro
    dbc.Row([
        dbc.Col([], width=2),
        dbc.Col([
            html.H2("Welcome to the Pediatric Single Cell Cancer Atlas!"),
            html.Br(),
            html.H4("(1) Choose an analysis page:"),
            html.Br(),
            html.Div([
                dcc.Dropdown([
                    {
                        "label": html.Span(['All Datasets'], style={'font-size':17}),
                        'value':'overview'
                    },
                    {
                        "label": html.Span(['Pediatric Leukemia'], style={'font-size':17}),
                        "value": "pediatric_leukemia"
                    },
                    {
                     "label": html.Span(['AML (Lambo, 2023)'], style={'font-size':17}),
                     "value": "aml_lambo_2023"
                    },
                    {
                     "label": html.Span(['Healthy Pediatric'], style={'font-size':17}),
                     "value": "healthy_pediatric"                        
                    }
                ], value='overview',id='dataset_choice', multi=False)
            ], style={"width":"35%"})
        ]),
        dbc.Col([], width=2),
    ], align="center"),

    html.Br(),html.Br(),

    # email
    dbc.Row([
        dbc.Col([], width=2),
        dbc.Col([
            html.H4(
                "(2) Please enter some information below to continue. This information is used for tool analytical purposes only."
            ),
            dbc.FormFloating([
                dbc.Input(type="email", placeholder="example@inst.edu", id="user_email"),
                dbc.Label("Email"),
            ]),
        ]),
        dbc.Col([], width=2)
    ]),

    html.Br(),

    # affiliation and country
    dbc.Row([
        dbc.Col([], width=2),
        dbc.Col([
            dbc.FormFloating([
                dbc.Input(type="text", placeholder="University", id="user_affil"),
                dbc.Label("Affiliation")
            ])
        ]),
        dbc.Col([
            dbc.FormFloating([
                dbc.Input(type="text", placeholder="Country", id="user_country"),
                dbc.Label("Country")
            ])
        ]),
        dbc.Col([html.Div(id='dummy1'),], width=2)
    ]),

    html.Br(),html.Br(),

    # enter button
    dbc.Row([
        dbc.Col([], width=2),
        dbc.Col([
            # dcc.Link(dbc.Button("Start Analysis", color="primary", id="enter", n_clicks=0), href="/tool", refresh=True)
            dcc.Link(dbc.Button("Start Analysis", id='enter', color='primary', n_clicks=0, outline=False), href=prefix+'/')
            # dcc.Link(dbc.Button("Start Analysis", color="primary", id="enter", n_clicks=0), href="/PediatricSCAtlas-Dev/tool", refresh=True)
        ]),
        dbc.Col([], width=2)
    ]),

    html.Br(),html.Br(),

    # github, more information
    dbc.Row([
        dbc.Col([], width=2),
        dbc.Col([
            html.H5(
                children=dcc.Markdown(
                    '''The **PedSCAtlas** was developed by [Bhasin Systems Biomedicine Lab](https://bhasinlab.org/) at Emory University
                    in collaboration with the [Aflac Cancer & Blood Disorders Center](https://www.choa.org/medical-services/cancer-and-blood-disorders) at 
                    Children's Healthcare of Atlanta. If you use our tool, please cite our abstract published in 
                    [*blood*](https://ashpublications.org/blood/article/140/Supplement%201/2278/489360/A-Single-Cell-Atlas-and-Interactive-Web-Resource), 
                    presented during the ASH 2022 conference. Manuscript is in development.
                    '''
                )
            )
        ]),
        dbc.Col([], width=2)
    ]),

    html.Br(),

    dbc.Row([
        dbc.Col([], width=2),
        dbc.Col([
            html.H5(
                children=dcc.Markdown(''' For more information, including a detailed user guide, please visit our [github](https://github.com/bhasin-lab/PedSCAtlas). You can post questions or issues there.
                ''')
            )
        ]),
        dbc.Col([], width=2)
    ])

    
])

@callback(
    Output('dummy1','children'),
    Input('enter', 'n_clicks'),
    State('user_email', 'value'),
    State('user_affil','value'),
    State('user_country','value')
)
def store_info(n_clicks,user_email, user_affil, user_country):
    email_txt = user_email
    affil_txt = user_affil
    country_txt = user_country
    date_txt = str(date.today())

    if n_clicks > 0:
        if user_email is None:
            email_txt = "N/A"
        if user_affil is None:
            affil_txt = "N/A"
        if user_country is None:
            country_txt = "N/A"
        
        f = open('./data/user/user_information.txt','a')
        line = [email_txt + "\t" + affil_txt + "\t" + country_txt + "\t" + date_txt + "\n"]
        f.writelines(line)
        f.close()

    return None

# button takes user to a certain page
@callback(
    Output("enter","href"),
    Input("dataset_choice","value"),
    prevent_initial_call=False
)
def button_link(value):
    out = prefix + "/tool/" + value
    # print(out)
    return out