from dash import html
import dash_bootstrap_components as dbc

# Define the navbar structure
from config import MODE

if MODE=="ONLINE":
    prefix = "/PediatricSCAtlas"
elif MODE=="LOCAL":
    prefix = ""


def Navbar():

    layout = html.Div([
        dbc.NavbarSimple(
            children=[
                # pages
                dbc.DropdownMenu([
                    dbc.DropdownMenuItem("Analysis pages", header=True),
                    dbc.DropdownMenuItem("All Datasets", href=prefix+'/tool/overview'),
                    dbc.DropdownMenuItem("Pediatric Leukemia", href=prefix+'/tool/pediatric_leukemia'),
                    dbc.DropdownMenuItem("AML (Lambo, 2023)", href=prefix+'/tool/aml_lambo_2023'),
                    dbc.DropdownMenuItem("Healthy Pediatric BM", href=prefix+'/tool/healthy_pediatric')
                ], nav=True, in_navbar=True, label="Analysis"),
                dbc.NavItem(dbc.NavLink("DE", href=prefix+"/DE")),
                dbc.NavItem(dbc.NavLink("Marker Testing", href=prefix+"/testing")),

                # external links
                dbc.NavItem(dbc.Button([html.I(className="fa-brands fa-github")], className="me-1", outline=False, href="https://github.com/bhasin-lab/PedSCAtlas")),
                dbc.NavItem(dbc.Button([html.I(className="fa-solid fa-envelope")], className="me-1", outline=False, href="mailto:hmumme@emory.edu")),
            ], 
            brand=html.H2("PedSCAtlas"),
            brand_href=prefix+"/",
            color="dark",
            dark=True,
        ), 
    ])

    return layout