# var - variable to group x-axis of bar plot
# coords - dataframe for disease subtype
# type - either 'Actual' for acutal cell counts or 'Proportion' for proportion of total cells per var grouping

import plotly.express as px
import pandas as pd
import numpy as np
import plotly.graph_objects as go

def get_prop (var, coords, type, colors, title=True, x_title=True, dataset=""):
    # get actual counts of cells
    prop = coords.groupby(var, as_index=False, observed=True)[['Cell Type']].value_counts()

    # initiate to prevent error as described here: https://github.com/plotly/plotly.py/issues/3441
    fig = go.Figure(layout=dict(template='plotly'))    

    # plot or calculate proportions and plot
    if type == 'Actual':
        fig = px.bar(prop, x=var, y="count", color="Cell Type", title="Number of Cells Per Cell Type in " + var, category_orders={"Cell Type":coords["Cell Type"].cat.categories.values.tolist()},
                            template = 'plotly_white', color_discrete_sequence=colors,
                            labels={"count": "Cell Count" })

    elif type == 'Proportion':
        # calculate total counts per var
        total = coords[var].value_counts()

        # make new proportion column and calculate proportions per cell type
        prop['proportion'] = 0.0

        # get variables in var
        var_groups = coords[var].dropna().unique().tolist()

        # calculate proportion for each group
        for i in var_groups:
            prop.loc[prop[var] == i,"proportion"] = prop.loc[prop[var] == i,"count"] / total[i]

        # plot
        fig = px.bar(prop, x=var, y="proportion", color="Cell Type", title="Proportion of Cells Per Cell Type in " + var, category_orders={"Cell Type":coords["Cell Type"].cat.categories.values.tolist()},
                            template = 'plotly_white', color_discrete_sequence=colors,
                            labels={"proportion": "Proportion of Total" })
    
    if not x_title:
        fig.update_layout(xaxis_title=None)
    if not title: 
        fig.update_layout(title=None)
    return fig