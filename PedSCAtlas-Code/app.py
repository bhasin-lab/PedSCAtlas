from dash import Dash, html, dcc
import dash
import dash_bootstrap_components as dbc

from components import navbar

# set version - ONLINE means published, LOCAL means development
from config import MODE

if MODE=="ONLINE":
    app = Dash(__name__, external_stylesheets=[dbc.themes.LUMEN, dbc.icons.FONT_AWESOME], use_pages=True, requests_pathname_prefix='/PediatricSCAtlas/')
elif MODE=="LOCAL":
	app = Dash(__name__, use_pages=True, external_stylesheets=[dbc.themes.LUMEN, dbc.icons.FONT_AWESOME])

nav = navbar.Navbar()

app.layout = html.Div([
	nav,

	dash.page_container
])

if __name__ == '__main__':
    if MODE=="ONLINE":
         app.run(debug=False, host="192.168.1.108", port=8051, proxy="bhasinlab.bmi.emory.edu/PediatricSCAtlas")
    elif MODE=="LOCAL":
    	app.run_server(debug=True)