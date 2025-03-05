import dash
from dash import html, dcc, callback, Input, Output, State, dash_table
import pandas as pd
import dash_bootstrap_components as dbc
from datetime import date
import numpy as np
import plotly.express as px
import pickle

from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn import metrics
from random import sample

from config import MODE

# file header - establish location of "Data" directory
if MODE=="ONLINE":
    dir = '/srv/dash-apps/PedSCAtlas-Data/'
elif MODE=="LOCAL":
    dir = '/opt/localdata/hmumme/PedSCAtlas/Data/'

# page prefix, same across ONLINE + LOCAL
page_pre = '/'

dash.register_page(__name__, path=page_pre + 'testing', title='PedSCAtlas', image='PedSCAtlas_logo.png')

# load single-cell dataset
df = pd.read_csv(dir + 'AcuteLeukemia/scMD.csv')
df = df.astype({'UMAP_1': 'float64', 'UMAP_2': 'float64', 'Cell Id':'str', 'Sample Id':'category','Cell Type':'category',
                'Diagnosis':'category', 'Cluster':'category', 'Gene_Expression': 'float64',
                'Source':'category',
                'Time Point':'category','EOI MRD Status':'category','Mito Per':'float64'}) # Change type of columns
sc_exp = pickle.load(open(dir + "AcuteLeukemia/geneExp_norm.p", "rb" ))
sc_gene_names = pd.read_csv(dir + 'AcuteLeukemia/sc_gene_names.csv')
sc_gene_names = sc_gene_names.Gene.to_list()

# load de tables
sc_de = pd.read_csv(dir + "AcuteLeukemia/sc_de.csv", index_col=0)
bulk_de = pd.read_csv(dir + "AcuteLeukemia/bulk_de.csv", index_col=0)
#res = pd.read_csv(dir + 'AcuteLeukemia/res_de.csv', index_col=0)

# sc and bulk tested genes
sc_gene_names = pd.read_csv(dir + 'AcuteLeukemia/sc_gene_names.csv')
sc_gene_names = sc_gene_names.Gene.to_list()

bulk_gene_names = pd.read_csv(dir + 'AcuteLeukemia/bulk_gene_names.csv')
bulk_gene_names = bulk_gene_names.gene.to_list()

# placeholder summary table
de_res = pd.DataFrame({"RNA_Type":[""],"Comparison":[""],"FC(log2)":[0],"Adj_P_Val":[0]})

# intersection
all_genes = list(set(sc_gene_names).intersection(set(bulk_gene_names)))

# disease types
all_types = bulk_de.type.unique()
all_types[np.where(all_types == "Pre B-ALL")] = "B-ALL"

# setup model stats dataframe
conf_stats = pd.DataFrame({
    'Metric':['PPV','NPV','Sensitivity','Specificity'],
    'Value':[0,0,0,0]
})


layout = html.Div(children=[

    html.Br(),

    dbc.Row([
        dbc.Col(width=1),
        dbc.Col(
            dbc.Container([
                html.H2("Input gene of interest"),
                html.Hr(),
                html.Div([
                    dcc.Markdown('''
                        Marker testing in single-cell (SC) and bulk RNA-seq diseased and healthy data, from the **Pediatric Acute Leukemia** study. Input a gene and disease type to show expresion (differential and average) in different datasets.
                    ''')
                ]),
                dbc.Label("Input parameters:"),
                html.Div([
                    dcc.Dropdown(
                        id = 'g_test',
                        options = all_genes,
                        value = 'PTEN',
                        placeholder="Gene"
                )],style={"width": "50%"},),
                html.Br(),
                html.Div([
                    dcc.Dropdown(
                    id = 't_test',
                    options = all_types,
                    value = 'T/My MPAL',
                    placeholder="Disease Type")
                ],style={"width": "50%"},)

            ])
        ),

        dbc.Col(
            dbc.Container([
                html.H2("Differential Expression Testing"),
                html.Hr(),
                html.Div([
                    dcc.Markdown('''
                        DE results in SC and Bulk RNA-seq data. Comparisons and resulting log2FC, adjusted p-values are shown. Exact parameter information is located on our [Github](https://github.com/bhasin-lab/PedSCAtlas/). If results are shown as "Not Found", this means the log2FC result for this gene in this subtype was less than 0.25, and was therefore not reported during DE.
                    '''),
                ]),
                #dcc.Dropdown(),
                #html.Br(),
                dash_table.DataTable(id="de_res", data = de_res.to_dict('records'), columns = [{"name": i, "id": i} for i in de_res.columns])
            ])
        ),

        dbc.Col(width=1),
    ]),

    html.Br(),html.Br(),

    dbc.Row([
        dbc.Col(width=1),
        dbc.Col(
            dbc.Container([
                html.H2("Random Forest Classifier"),
                html.Hr(),
                html.Div([
                    dcc.Markdown('''A binary random forest model has been trained to classify the chosen leukemia 
                                 type's blast cells and healthy BM as blast and healthy BM cells based on your 
                                 chosen gene's expression. A *confusion matrix* is shown below (as a heatmap)
                                  with the test results. In addition, the testing results are summarized in a table.
                                  A "True Positive" means a Leukemia Blast was predicted to be a Leukemia Blast 
                                 based on only the selected gene's expression. See our 
                                 [Github](https://github.com/bhasin-lab/PedSCAtlas/) for exact 
                                 training/testing parameters. *Note*: the training and testing sets are randomly assigned 
                                 during each model formation; therefore, the results will change each time the application 
                                 parameters are adjusted or the page is refreshed.''')
                ]),
                dbc.Row([
                html.H3("Confusion Matrix"),
                dbc.Col(
                    html.Div([
                        dbc.Spinner([dcc.Graph(id='conf_matrix')])
                    ])
                ),
                dbc.Col(
                    html.Div([
                        html.H3('Testing Summary'),
                         dbc.Spinner([dash_table.DataTable(id="conf_stats", data = conf_stats.to_dict('records'), columns = [{"name": i, "id": i} for i in conf_stats.columns])])
                    ])
                )
                ]),          
            ])
        ),
        dbc.Col(width=1),
    ]),

    html.Hr(),

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

@callback(
    Output('de_res','data'),
    Input('g_test','value'),
    Input('t_test','value')
)
def user_de (g_test, t_test):
    de_res = pd.DataFrame({'RNA_Type':['SC','Bulk'], 'Comparison':['Malignant Cells vs. Healthy BM','Diseased vs. Normal Samples'], 'FC(log2)':[0.0,0.0], 'Adj_P_Val':[0.0,0.0]})
    # get sc results
    sc_u = sc_de.loc[(sc_de.type == t_test + " Blast") & (sc_de.gene == g_test)]
    
    # if gene not found, it means FC<0.25, add to de_res
    if len(sc_u.index) == 0:
         de_res.iloc[0,2] = 'Not Found'
         de_res.iloc[0,3] = 'Not Found'
    else:
         de_res.iloc[0,2] = round(sc_u["log2FC"].values[0],3)
         de_res.iloc[0,3] = round(sc_u["p_val_adj"].values[0],3)

    # get bulk results
         
    if t_test == "B-ALL":
        bulk_u = bulk_de.loc[(bulk_de.type == "Pre B-ALL") & (bulk_de.gene == g_test)]
    else:
        bulk_u = bulk_de.loc[(bulk_de.type == t_test) & (bulk_de.gene == g_test)]

    # add to     
    de_res.iloc[1,2]  = round(bulk_u["log2FC"].values[0],3)
    de_res.iloc[1,3]  = round(bulk_u["p_val_adj"].values[0],3)

    return de_res.to_dict('records')

@callback(
    Output('conf_matrix','figure'),
    Output('conf_stats','data'),
    Input('g_test','value'),
    Input('t_test','value')
)
def run_model(g_test, t_test):
    df_mod = df.copy()
    df_mod["Gene_Expression"] = sc_exp[sc_gene_names.index(g_test),:].toarray()[0]
    df_mod_u = df_mod.loc[(df_mod["Cell Type"] == t_test + " Blast") | (df_mod.Diagnosis == "Healthy BM")].copy()

    # set up binary classification
    df_mod_u["Class"] = "Blast"
    df_mod_u.loc[df_mod_u["Cell Type"] != t_test + " Blast", "Class"] = "Healthy Cell"

    # if AML, randomly downsample to 20k blast cells so classes are more balanced
    if t_test == "AML":
        healthy_inds = df_mod_u[df_mod_u['Class'] == "Healthy Cell"].index.tolist()
        aml_inds = sample(df_mod_u[df_mod_u['Class'] == "Blast"].index.tolist(),20000)
        df_mod_u = df_mod_u.loc[aml_inds+healthy_inds]

    # y 
    classes = df_mod_u["Class"].astype("str")

    # X array
    exp_data = df_mod_u.Gene_Expression.to_numpy().reshape(-1,1)

    # Perform train/test split (70/30)
    X_train, X_test, y_train, y_test = train_test_split(
        exp_data, classes,
        test_size=0.3, shuffle=True
    )

    # generate RF model
    clf = RandomForestClassifier()

    # Train (Fit)
    clf.fit(X_train, y_train)

    # Predict
    preds = clf.predict(X_test)

    # Confusion Matrix
    c_mat = metrics.confusion_matrix(y_test, preds, labels = ["Blast","Healthy Cell"])

    # plot heatmap
    conf_matrix = px.imshow(c_mat, text_auto=True,
                labels = dict(x='Predicted Class', y='True Class'),
                x=['Blast','Healthy Cell'],
                y=['Blast','Healthy Cell'],
                color_continuous_scale='reds',
                template = 'plotly_white')
    
    # generate statistics table
    conf_stats = pd.DataFrame({
        'Metric':['PPV','NPV','Sensitivity','Specificity'],
        'Value':[0.0,0.0,0.0,0.0]
    })
    conf_stats.loc[conf_stats.Metric == "PPV",'Value'] = round(c_mat[0,0] / (c_mat[0,0]+c_mat[0,1]),3)
    conf_stats.loc[conf_stats.Metric == "NPV",'Value'] = round(c_mat[1,1] / (c_mat[1,1]+c_mat[1,0]),3)
    conf_stats.loc[conf_stats.Metric == "Sensitivity",'Value'] = round(c_mat[0,0] / (c_mat[0,0]+c_mat[1,0]),3)
    conf_stats.loc[conf_stats.Metric == "Specificity",'Value'] = round(c_mat[1,1] / (c_mat[1,1]+c_mat[0,1]),3)

    del clf # make sure model is deleted after each call

    return conf_matrix, conf_stats.to_dict('records')