# scatter functions for one-page view

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

#from load_data import DataLambo, DataHealthyPD, DataAcuteLeukemia
from load_data import DataAcuteLeukemia, DataLambo, DataHealthyPD

# load helper functions
from helper.survival import *
from helper.get_prop import * 
from helper.utils import id_factory

# scatter_pan_leukemia
# Pan-Leukemia Dataset
def scatter_pan_leukemia(sc_scatter_subtype, sc_pt_color, page):
    print("Pan-Leuk:", sc_pt_color)

    # initiate to prevent error as described here: https://github.com/plotly/plotly.py/issues/3441
    sc_scatter = go.Figure(layout=dict(template='plotly'))

    # set subtypes to plot
    if len(sc_scatter_subtype)<1:
        sc_scatter = go.Figure()
    else:
        df_sub = DataAcuteLeukemia.df.loc[DataAcuteLeukemia.df["Diagnosis"].isin(sc_scatter_subtype)].copy()
        df_sub["Diagnosis"] = df_sub["Diagnosis"].cat.remove_unused_categories()
        #print(len(df_sub))
        
        # set color palette
        if (sc_pt_color=="Cell Type") & (len(sc_scatter_subtype)<len(DataAcuteLeukemia.df.Diagnosis.unique())):
            df_sub[sc_pt_color] = df_sub[sc_pt_color].cat.remove_unused_categories()
            plot_color = DataAcuteLeukemia.sc_colors.loc[DataAcuteLeukemia.sc_colors.CellType.isin(df_sub[sc_pt_color].cat.categories)].Color.tolist()
        elif sc_pt_color=="Cell Type":
            plot_color = DataAcuteLeukemia.sc_colors.Color.tolist()
        else:
            plot_color = px.colors.qualitative.Dark24
        #print(plot_color)
        
        # check if there is any information after removing nans
        df_sub = df_sub[-df_sub[sc_pt_color].isna()]
        #print(df_sub.shape)
        if df_sub.shape[0] == 0:
            sc_scatter = go.Figure()

        # generate plot depending on type of grouping
        elif df_sub[sc_pt_color].dtype.name=="category":
            df_sub[sc_pt_color] = df_sub[sc_pt_color].cat.remove_unused_categories()
            sc_scatter = px.scatter(df_sub, x="UMAP_1", y="UMAP_2", 
                                    color=sc_pt_color, 
                                    category_orders={sc_pt_color:df_sub[sc_pt_color].cat.categories.values.tolist()},
                                    template = 'plotly_white', 
                                    color_discrete_sequence=plot_color, 
                                    hover_data=["Sample Id","Diagnosis"])
        else:
            sc_scatter = px.scatter(df_sub, x="UMAP_1", y="UMAP_2", color=sc_pt_color,
                                    template = 'plotly_white', color_discrete_sequence=plot_color, hover_data=["Sample Id","Diagnosis"])
    
    # move leged to bottom of plot on overview page
    if page=='overview':
        sc_scatter.update_layout(
            legend=dict(
                orientation='h',
                yanchor='top',
                y=-0.25, x=-0.5,
                xanchor='left'
            )
        )
    
    return sc_scatter


# scatter_lambo
# AML Lambo 2023
def scatter_lambo(sc_pt_color, page):
    print("Lambo:",sc_pt_color)
    df_sub = DataLambo.df.copy()

    # initiate to prevent error as described here: https://github.com/plotly/plotly.py/issues/3441
    sc_scatter_lambo = go.Figure(layout=dict(template='plotly'))

    if sc_pt_color=="Cell Type":
        plot_color = DataLambo.sc_colors.Color.tolist()
    else:
        plot_color = px.colors.qualitative.Dark24

    # check if there is any information after removing nans
    df_sub = df_sub[-df_sub[sc_pt_color].isna()]
    if df_sub.shape[0] == 0:
        sc_scatter_lambo = go.Figure()

    # generate plot depending on type of grouping
    elif df_sub[sc_pt_color].dtype.name=="category":
        df_sub[sc_pt_color] = df_sub[sc_pt_color].cat.remove_unused_categories()
        sc_scatter_lambo = px.scatter(df_sub, x="UMAP_1", y="UMAP_2", 
                                color=sc_pt_color, 
                                category_orders={sc_pt_color:df_sub[sc_pt_color].cat.categories.values.tolist()},
                                template = 'plotly_white', 
                                color_discrete_sequence=plot_color, 
                                hover_data=["Sample Id","Diagnosis"])
    else:
        sc_scatter_lambo = px.scatter(df_sub, x="UMAP_1", y="UMAP_2", color=sc_pt_color,
                                template = 'plotly_white', color_discrete_sequence=plot_color, hover_data=["Sample Id","Diagnosis"])
    
    # move leged to bottom of plot on overview page
    if page=='overview':
        sc_scatter_lambo.update_layout(
            legend=dict(
                orientation='h',
                yanchor='top',
                y=-0.25, x=-0.5,
                xanchor='left'
            )
        )

    return sc_scatter_lambo

# scatter_healthyPD 
# healthy pediatric BM
def scatter_healthyPD(sc_pt_color, page):
    print("Healthy Ped:", sc_pt_color)
    df_sub = DataHealthyPD.df.copy()

    # initiate to prevent error as described here: https://github.com/plotly/plotly.py/issues/3441
    sc_scatter_hb = go.Figure(layout=dict(template='plotly'))

    if sc_pt_color=="Cell Type":
        plot_color = DataHealthyPD.sc_colors.Color.tolist()
    else:
        plot_color = px.colors.qualitative.Dark24

    # check if there is any information after removing nans
    df_sub = df_sub[-df_sub[sc_pt_color].isna()]
    if df_sub.shape[0] == 0:
        sc_scatter_hb = go.Figure()

    # generate plot depending on type of grouping
    elif df_sub[sc_pt_color].dtype.name=="category":
        df_sub[sc_pt_color] = df_sub[sc_pt_color].cat.remove_unused_categories()
        sc_scatter_hb = px.scatter(df_sub, x="UMAP_1", y="UMAP_2", 
                                color=sc_pt_color, 
                                category_orders={sc_pt_color:df_sub[sc_pt_color].cat.categories.values.tolist()},
                                template = 'plotly_white', 
                                color_discrete_sequence=plot_color, 
                                hover_data=["Sample Id"])
    else:
        sc_scatter_hb = px.scatter(df_sub, x="UMAP_1", y="UMAP_2", color=sc_pt_color,
                                template = 'plotly_white', color_discrete_sequence=plot_color, hover_data=["Sample Id"])
    
    # move leged to bottom of plot on overview page
    if page=='overview':
        sc_scatter_hb.update_layout(
            legend=dict(
                orientation='h',
                yanchor='top',
                y=-0.25, x=-0.5,
                xanchor='left'
            )
        )

    return sc_scatter_hb


