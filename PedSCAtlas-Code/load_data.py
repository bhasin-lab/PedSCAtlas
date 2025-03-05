import pandas as pd
import numpy as np
from itertools import cycle
import pickle
from pathlib import Path
import matplotlib.pyplot as plt
import base64
from io import BytesIO
from config import MODE


if MODE=="ONLINE":
    dir = 'PedSCAtlas-Data/'
elif MODE=="LOCAL":
    dir = 'PedSCAtlas/Data/'


class PDAtlasData:
    def __init__(self, df, exp, gene_names, sc_table=None, dst_desc=None, sc_colors=None):
        self.df = df
        self.exp = exp
        self.gene_names = gene_names
        self.sc_table = sc_table
        self.dst_desc = dst_desc
        self.sc_colors = sc_colors


# ---------------------- load dataset Pan/AcuteLeukemia----------------------
#df = pd.read_csv(dir + 'AcuteLeukemia/scMD.csv')
def load_AcuteLeukemia():
    print("=======\nload_AcuteLeukemia=======\n=======\n=======\n")
    df = pd.read_csv(dir + 'AcuteLeukemia/scMD_wSub.csv') # file with "Subdiagnosis column"
    df = df.astype({'UMAP_1': 'float64', 'UMAP_2': 'float64', 'Cell Id':'str', 'Sample Id':'category','Cell Type':'category',
                    'Diagnosis':'category', 'Cluster':'category', 'Gene_Expression': 'float64',
                    'Source':'category',
                    'Time Point':'category','EOI MRD Status':'category',
                    'Mito Per':'float64','Subdiagnosis':'category'}) # assert types of columns
    #bulk_data = pickle.load(open(dir + "TARGET/bulk_EXP.p", "rb"))
    #bulk_md = pd.read_csv(dir + "TARGET/bulk_MD.csv", index_col=0)
    sc_exp = pickle.load(open(dir + "AcuteLeukemia/geneExp_norm.p", "rb" ))
    sc_gene_names = pd.read_csv(dir + 'AcuteLeukemia/sc_gene_names.csv')
    sc_gene_names = sc_gene_names.Gene.to_list()

    sc_table = pd.read_table(dir + "AcuteLeukemia/desc_table.txt")
    dst_desc = Path("./data/description/acuteLeukemia.txt").read_text()
    sc_colors = pd.read_csv(dir + "AcuteLeukemia/sc_colors.csv")
    
    return(PDAtlasData(df = df,
        exp = sc_exp,
        gene_names = sc_gene_names,
        sc_table = sc_table,
        dst_desc = dst_desc,
        sc_colors = sc_colors))

DataAcuteLeukemia = load_AcuteLeukemia()


# ---------------------- load datasets AML_Lambo----------------------
# AML Lambo 2023 Single-Cell
def load_Lambo():
    print("=======\nDataAML_Lambo=======\n=======\n=======\n")
    df = pd.read_csv(dir + 'AML_Lambo/scMD.csv')
    df = df.astype({'UMAP_1': 'float64', 'UMAP_2': 'float64', 'Cell Id':'str', 'Sample Id':'category','Cell Type':'category',
                    'Diagnosis':'category', 'Cluster':'category', 'Gene_Expression': 'float64',
                    'Source':'category', 'Malignant':'category','Time Point':'category','Patient':'category',
                    'Subdiagnosis':'category'}) # assert types of columns
    sc_exp = pickle.load(open(dir + "AML_Lambo/geneExp_norm.p", "rb" ))
    sc_gene_names = pd.read_csv(dir + 'AML_Lambo/sc_gene_names.csv')
    sc_gene_names = sc_gene_names.Gene.to_list()
    # Tables
    dst_desc = Path("./data/description/AML_Lambo.txt").read_text()
    sc_colors = pd.read_csv(dir + "AML_Lambo/sc_colors.csv")
    return(PDAtlasData(df = df,
        exp = sc_exp,
        gene_names = sc_gene_names,
        dst_desc = dst_desc,
        sc_colors = sc_colors))

DataLambo = load_Lambo()


# ---------------------- load dataset Healthy Ped----------------------
# Healthy Ped i.e. HB2
def load_HealthyPed(): 
    print("=======\nload_HealthyPed=======\n=======\n=======\n")
    df = pd.read_csv(dir + 'Healthy_Ped/scMD.csv')
    df = df.astype({'UMAP_1': 'float64', 'UMAP_2': 'float64', 'Cell Id':'str', 'Sample Id':'category','Cell Type':'category',
                    'Cluster':'category', 'Gene_Expression': 'float64'}) # assert types of columns
    sc_exp = pickle.load(open(dir + "Healthy_Ped/geneExp_norm.p", "rb" ))
    sc_gene_names = pd.read_csv(dir + 'Healthy_Ped/sc_gene_names.csv')
    sc_gene_names = sc_gene_names.Gene.to_list()

    # Tables
    dst_desc = Path("./data/description/Healthy_Ped.txt").read_text()
    sc_colors = pd.read_csv(dir + "Healthy_Ped/sc_colors.csv")
    return(PDAtlasData(df = df,
        exp = sc_exp,
        gene_names = sc_gene_names,
        dst_desc = dst_desc,
        sc_colors = sc_colors))

DataHealthyPD = load_HealthyPed()


# ---------------------- load datasets HealthyHCA----------------------
# Healthy HCA
def load_healthy():
    df_healthy = pd.read_table(dir + "HCA/coords.txt")
    healthy_exp = pickle.load(open(dir + "HCA/geneExp_norm.p", "rb" ))
    healthy_gene_names = pd.read_csv(dir + 'HCA/gene_names.csv')
    healthy_gene_names = healthy_gene_names.Gene.to_list()
    print("=======\nDataHealthy=======\n=======\n=======\n")
    return(PDAtlasData(df=df_healthy,
        exp = healthy_exp,
        gene_names = healthy_gene_names))
DataHealthy = load_healthy()