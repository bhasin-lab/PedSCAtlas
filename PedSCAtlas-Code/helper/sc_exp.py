# Generate single-cell expression plots for datasets

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

from load_data import DataLambo, DataHealthyPD, DataAcuteLeukemia, DataHealthy

# load helper functions
from helper.survival import *
from helper.get_prop import * 
from helper.utils import id_factory

# load_scExp_lambo
# single-cell expression for AML Lambo 2023
def load_scExp_lambo(sc_gene, sc_gene_type, sc_group_by):
    # extract gene expression, add to coords DataLambo.df
    DataLambo.df['Gene_Expression'] = DataLambo.exp[DataLambo.gene_names.index(sc_gene),:].toarray()[0]

    #df_sub = DataLambo.df.loc[DataLambo.df["Diagnosis"].isin(sc_exp_subtype)].copy()
    #df_sub[sc_group_by] = df_sub[sc_group_by].cat.remove_unused_categories()

    df_sub = DataLambo.df.copy()

    # initiate to prevent error as described here: https://github.com/plotly/plotly.py/issues/3441
    sc_exp_fig_lambo = go.Figure(layout=dict(template='plotly'))  

    # plot title
    TITLE = sc_gene + " Expression - AML (Lambo 2023)"

    # check if cell type group is chosen
    #if (sc_group_by=="Cell Type") & (len(sc_exp_subtype)<len(DataLambo.df.Diagnosis.unique())):
    #    plot_color = DataLambo.sc_colors.loc[DataLambo.sc_colors.CellType.isin(df_sub[sc_group_by].cat.categories)].Color.values    
    if sc_group_by=="Cell Type":
        plot_color = DataLambo.sc_colors.Color                
    else:
        plot_color = px.colors.qualitative.Dark24

    # make figure depending on plot type
    if sc_gene_type == "Violin":
        sc_exp_fig_lambo = go.Figure()
        palette = cycle(plot_color)
        for i in df_sub[sc_group_by].cat.categories:
            df_i = df_sub.loc[df_sub[sc_group_by] == i]
            sc_exp_fig_lambo.add_trace(go.Violin(
                x=df_i[sc_group_by], y=df_i['Gene_Expression'],
                box_visible=True,
                name=i,
                spanmode='hard',
                marker_color=next(palette)
            ))
        sc_exp_fig_lambo.update_layout(template='plotly_white', title=TITLE)
        sc_exp_fig_lambo.update_yaxes(title_text='Log-Normalized Counts')

    elif sc_gene_type == "Box-Plot":
        sc_exp_fig_lambo = px.box(df_sub, x=sc_group_by, y='Gene_Expression', color=sc_group_by,
                            title = TITLE,
                            category_orders={sc_group_by:df_sub[sc_group_by].cat.categories.values.tolist()},
                            template = 'plotly_white', color_discrete_sequence=plot_color,
                            labels={"Gene_Expression":"Log-Normalized Counts"})
    elif sc_gene_type == "Feature":
        sc_exp_fig_lambo = px.scatter(df_sub, x="UMAP_1", y="UMAP_2", color='Gene_Expression', color_continuous_scale="Purples",
                            title = sc_gene + " Expression - AML (Lambo 2023)",
                            template = 'plotly_white',
                            hover_data=[sc_group_by,"Diagnosis"],
                            labels={"Gene_Expression":"Log-Normalized Counts"})
    #elif sc_gene_type == "Jitter":
    #    sc_exp_fig_lambo = px.box(DataLambo.df, x=sc_group_by, y='Gene_Expression', color=sc_group_by,
    #                        title = TITLE,
    #                        category_orders={sc_group_by:DataLambo.df[sc_group_by].cat.categories.values.tolist()},
    #                        template = 'plotly_white',
    #                        points='all', color_discrete_sequence=plot_color,
    #                        labels={"Gene_Expression":"Log-Normalized Counts"})
    return sc_exp_fig_lambo

# load_scExp_hb
# single-cell expression from Ped Healthy BM
def load_scExp_hb(sc_gene, sc_gene_type, sc_group_by):
    # extract gene expression, add to coords DataHealthyPD.df
    DataHealthyPD.df['Gene_Expression'] = DataHealthyPD.exp[DataHealthyPD.gene_names.index(sc_gene),:].toarray()[0]

    df_sub = DataHealthyPD.df.copy()

    # initiate to prevent error as described here: https://github.com/plotly/plotly.py/issues/3441
    sc_exp_fig_hb = go.Figure(layout=dict(template='plotly'))  

    # plot title
    TITLE = sc_gene + " Expression - Healthy Pediatric BM"

    if sc_group_by=="Cell Type":
        plot_color = DataHealthyPD.sc_colors.Color                
    else:
        plot_color = px.colors.qualitative.Dark24

    # make figure depending on plot type
    if sc_gene_type == "Violin":
        sc_exp_fig_hb = go.Figure()
        palette = cycle(plot_color)
        for i in df_sub[sc_group_by].cat.categories:
            df_i = df_sub.loc[df_sub[sc_group_by] == i]
            sc_exp_fig_hb.add_trace(go.Violin(
                x=df_i[sc_group_by], y=df_i['Gene_Expression'],
                box_visible=True,
                name=i,
                spanmode='hard',
                marker_color=next(palette)
            ))
        sc_exp_fig_hb.update_layout(template='plotly_white', title=TITLE)
        sc_exp_fig_hb.update_yaxes(title_text='Log-Normalized Counts')

    elif sc_gene_type == "Box-Plot":
        sc_exp_fig_hb = px.box(df_sub, x=sc_group_by, y='Gene_Expression', color=sc_group_by,
                            title = TITLE,
                            category_orders={sc_group_by:df_sub[sc_group_by].cat.categories.values.tolist()},
                            template = 'plotly_white', color_discrete_sequence=plot_color,
                            labels={"Gene_Expression":"Log-Normalized Counts"})
    elif sc_gene_type == "Feature":
        sc_exp_fig_hb = px.scatter(df_sub, x="UMAP_1", y="UMAP_2", color='Gene_Expression', color_continuous_scale="Purples",
                            title = sc_gene + " Expression - Healthy Pediatric BM",
                            template = 'plotly_white',
                            hover_data=[sc_group_by],
                            labels={"Gene_Expression":"Log-Normalized Counts"})

    return sc_exp_fig_hb

# load_scExp_pan
# single-cell expression from Pan-Leukemia
def load_scExp_pan(sc_gene, sc_gene_type, sc_group_by, sc_exp_subtype):

    # extract gene expression, add to coords DataAcuteLeukemia.df
    DataAcuteLeukemia.df['Gene_Expression'] = DataAcuteLeukemia.exp[DataAcuteLeukemia.gene_names.index(sc_gene),:].toarray()[0]

    # initiate to prevent error as described here: https://github.com/plotly/plotly.py/issues/3441
    sc_exp_fig = go.Figure(layout=dict(template='plotly'))  

    # plot title
    TITLE = sc_gene + " Expression - Pediatric Leukemia Atlas"

    # subtypes
    if len(sc_exp_subtype)<1:
        sc_exp_fig = go.Figure()
    else:
        df_sub = DataAcuteLeukemia.df.loc[DataAcuteLeukemia.df["Diagnosis"].isin(sc_exp_subtype)].copy()
        df_sub[sc_group_by] = df_sub[sc_group_by].cat.remove_unused_categories()

    # check if cell type group is chosen
        if (sc_group_by=="Cell Type") & (len(sc_exp_subtype)<len(DataAcuteLeukemia.df.Diagnosis.unique())):
            plot_color = DataAcuteLeukemia.sc_colors.loc[DataAcuteLeukemia.sc_colors.CellType.isin(df_sub[sc_group_by].cat.categories)].Color.values    
        elif sc_group_by=="Cell Type":
            plot_color = DataAcuteLeukemia.sc_colors.Color                
        else:
            plot_color = px.colors.qualitative.Dark24

        # make figure depending on plot type
        if sc_gene_type == "Violin":
            sc_exp_fig = go.Figure()
            palette = cycle(plot_color)
            for i in df_sub[sc_group_by].cat.categories:
                df_i = df_sub.loc[df_sub[sc_group_by] == i]
                sc_exp_fig.add_trace(go.Violin(
                    x=df_i[sc_group_by], y=df_i['Gene_Expression'],
                    box_visible=True,
                    name=i,
                    spanmode='hard',
                    marker_color=next(palette)
                ))
            sc_exp_fig.update_layout(template='plotly_white', title=TITLE)
            sc_exp_fig.update_yaxes(title_text='Log-Normalized Counts')

        elif sc_gene_type == "Box-Plot":
            sc_exp_fig = px.box(df_sub, x=sc_group_by, y='Gene_Expression', color=sc_group_by,
                                title = TITLE,
                                category_orders={sc_group_by:df_sub[sc_group_by].cat.categories.values.tolist()},
                                template = 'plotly_white', color_discrete_sequence=plot_color,
                                labels={"Gene_Expression":"Log-Normalized Counts"})
        elif sc_gene_type == "Feature":
            sc_exp_fig = px.scatter(df_sub, x="UMAP_1", y="UMAP_2", color='Gene_Expression', color_continuous_scale="Purples",
                                title = TITLE,
                                template = 'plotly_white',
                                hover_data=[sc_group_by,"Diagnosis"],
                                labels={"Gene_Expression":"Log-Normalized Counts"})
        #elif sc_gene_type == "Jitter":
        #    sc_exp_fig = px.box(DataAcuteLeukemia.df, x=sc_group_by, y='Gene_Expression', color=sc_group_by,
        #                        title = TITLE,
        #                        category_orders={sc_group_by:DataAcuteLeukemia.df[sc_group_by].cat.categories.values.tolist()},
        #                        template = 'plotly_white',
        #                        points='all', color_discrete_sequence=plot_color,
        #                        labels={"Gene_Expression":"Log-Normalized Counts"})

    return sc_exp_fig

# load_scExp_h_adult
# healthy Adult BM
def load_scExp_h_adult(healthy_gene, healthy_gene_type, healthy_celltype):
    # extract gene expression, add to coords DataAcuteLeukemia.df
    DataHealthy.df["Gene_Expression"] = DataHealthy.exp[DataHealthy.gene_names.index(healthy_gene),:].toarray()[0]

    # initiate to prevent error as described here: https://github.com/plotly/plotly.py/issues/3441
    healthy_exp_fig = go.Figure(layout=dict(template='plotly'))  

    # make figure depending on plot type
    if healthy_gene_type == "Violin":
        healthy_exp_fig = go.Figure()
        palette = cycle(px.colors.qualitative.Alphabet)
        for i in DataHealthy.df[healthy_celltype].unique():
            df_i = DataHealthy.df.loc[DataHealthy.df[healthy_celltype] == i]
            healthy_exp_fig.add_trace(go.Violin(
                x=df_i[healthy_celltype], y=df_i['Gene_Expression'],
                box_visible=True,
                name=i,
                spanmode='hard',
                marker_color=next(palette)
            ))
        healthy_exp_fig.update_layout(template='plotly_white', title=healthy_gene + " Expression")
        healthy_exp_fig.update_yaxes(title_text='Log-Normalized Counts')

    elif healthy_gene_type == "Box-Plot":
        healthy_exp_fig = px.box(DataHealthy.df, x=healthy_celltype, y="Gene_Expression", color=healthy_celltype,
                            title = healthy_gene + " Expression",
                            template = 'plotly_white', color_discrete_sequence=px.colors.qualitative.Alphabet,
                            labels={"Gene_Expression":"Log-Normalized Counts"})
    elif healthy_gene_type == "Feature":
        healthy_exp_fig = px.scatter(DataHealthy.df, x="UMAP_1", y="UMAP_2", color="Gene_Expression", color_continuous_scale="Purples",
                            title = healthy_gene + " Expression",
                            template = 'plotly_white', color_discrete_sequence=px.colors.qualitative.Alphabet,
                            labels={"Gene_Expression":"Log-Normalized Counts"})

    return healthy_exp_fig



# ADD OTHER DATASETS HERE

