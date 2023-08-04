## PedSCAtlas (November 04 2022) ##

# load required libraries
library(shiny)
library(shinydashboard)
library(shinydashboardPlus)
library(shinycssloaders)
library(spsComps)
library(dashboardthemes)
library(plotly)
library(ggpubr)
library(ggalluvial)
library(limma)
library(tidyverse)
library(GSVA)
source("functions.R")
source("theme.R")

shinyApp(
   ui = dashboardPage(
    options = list(sidebarExpandOnHover = FALSE),

    header = dashboardHeader(
      title = "PedSCAtlas", titleWidth = "300px"
    ),

    sidebar = dashboardSidebar(collapsed = FALSE, minified = FALSE, width = "300px", fixed = T,
      sidebarMenu(id = "sideMenu",
        menuItem("SingleCell", tabName = "pedCancer", startExpanded = FALSE, selected = FALSE,
        
        # set data show options for each dataset
          conditionalPanel(
              condition = "input.dataset == 'leukemia' & input.view != 'ribbon'",
              checkboxGroupInput('show', label = "Choose Disease Type to Show:",
                choices = c("AML" = "AML",
                        "B-ALL" = "BALL",
                        "T-ALL" = "TALL",
                        "B/Myeloid MPAL" = "BMPAL",
                        "T/Myeloid MPAL" = "TMPAL",
                       "Healthy Control" = "Con"), selected = c("AML","TALL","BALL","BMPAL","TMPAL","Con")
              )
          ),
          conditionalPanel(
            condition = "input.dataset == 'wilmsTumor' & input.view != 'ribbon'",
            checkboxGroupInput("show.wt", label = "Choose Disease Type to Show:", choices = c("Anaplastic" = "Anaplastic", "Favorable" = "Favorable"),
            selected = c("Anaplastic","Favorable"))
          ),
          conditionalPanel(
            condition = "input.dataset == 'rms' & input.view != 'ribbon'",
            checkboxGroupInput("show.rms", label = "Choose Disease Type to Show:",
            choices = c("Alveolar RMS" = "ARMS", "Embryonal RMS" = "ERMS"),
            selected = c("ARMS","ERMS"))
          ),
          conditionalPanel(
            condition = "input.dataset == 'epn' & input.view != 'ribbon'",
            checkboxGroupInput("show.epn", label = "Choose Disease Type to Show:",
              choices = c("Posterior Fossa (PF)" = "Posterior Fossa (PF)", "Spinal (SP)" = "Spinal (SP)", "Supratentorial (ST)" = "Supratentorial (ST)"), 
              selected = c("Posterior Fossa (PF)","Spinal (SP)","Supratentorial (ST)")
            )
          ),
          conditionalPanel(
            condition = "input.dataset == 'preALL' & input.view != 'ribbon'",
            checkboxGroupInput("show.preALL", label = "Choose Disease Type to Show:",
              choices = c("Pre-B ALL (HDD)" = "Pre-B HDD", "Pre-B ALL (ETV6-RUNX1)" = "Pre-B ETV6-RUNX1",
              "Pre-T ALL" = "Pre-T", "Healthy Control" = "Healthy"),
              selected = c("Pre-B HDD","Pre-B ETV6-RUNX1","Pre-T","Healthy")
            )
          ),

          # set data viewing options with different plots available
          selectInput('view', label = "Viewing Format:",
              choices = c("UMAP with Groups" = "UMAP", "Violin with Gene Expression" = "violin", "UMAP with Gene Expression" = "GeneExp_Show", "Violin with Pathway Enrichment" = "v_path", "UMAP with Pathway Enrichment" = "u_path", "Ribbon Plot with Multiple Pathways" = "ribbon")),
          
          # choose which data aspect to group by for the datasets
          conditionalPanel(condition = "input.view != 'ribbon' & input.dataset == 'leukemia'",
            selectInput("groupBy", label = "Data Aspect to Show:",
              c("Clusters" = "clusters",
                      "Cell Type" = "CellType",
                      "Disease Type" = "Group",
                      "Sample ID" = "sample",
                      "Rem or Rel" = "state")
            )
          ),
          conditionalPanel(condition = "input.view != 'ribbon' & input.dataset == 'wilmsTumor'",
            selectInput("groupBy_wt", label = "Data Aspect to Show:",
              c("Clusters" = "clusters",
                      "Cell Type" = "CellType",
                      "Disease Type" = "Group",
                      "Sample ID" = "sample")
            )
          ),
          conditionalPanel(
            condition = "input.view != 'ribbon' & input.dataset == 'rms'",
            selectInput("groupBy_rms", label = "Data Aspect to Show:",
              c("Clusters"="clusters","Disease Type"="Group","Sample ID"="sample","Cell Type"="CellType","Age"="age",
                "Gender"="gender","Disease Stage"="stage_group","Risk Stratum"="risk_stratum",
                "Time Point"="timePoint","Sample Site"="sampleSite","Cell Type Overall" = "CellType_Overall"
              )
            )
          ),
          conditionalPanel(
            condition = "input.view != 'ribbon' & input.dataset == 'epn'",
            selectInput("groupBy_epn", label = "Data Aspect to Show:",
              c("Clusters"="clusters","Disease Type"="Group","Sample ID"="sample","Cell Type"="CellType","Age"="age",
                "Gender"="gender","Therapy"="therapy","Patient Outcome"="outcome", "Broad Group" = "GroupBroad", "Prognosis" = "Prognosis",
                "Time Point"="timePoint","Sample Site"="sampleSite","Cell Type Overall" = "CellType_Overall", "Fine Group" = "GroupFine", 
                "Overall Survival (Years)" = "OverallSurvival"
              )
            )
          ),
          conditionalPanel(
            condition = "input.view != 'ribbon' & input.dataset == 'preALL'",
            selectInput("groupBy_preALL", label = "Data Aspect to Show:",
              c("Clusters"="clusters","Disease Type"="Group","Sample ID"="sample","Cell Type"="CellType",
                "Prognosis"="Prognosis","Cell Cycle Phase"="CellCyclePhase","Overall Cell Type"="CellType_Overall", "Specific Cell Type" = "CellType_Specific"
              )
            )
          ),

          # for ribbon plots, choose which groups to compare for each dataset (when more than one group is available)
          conditionalPanel(condition = "input.view == 'ribbon' & input.dataset == 'leukemia'",
            selectInput("g_comp", label = "Select two groups to compare",
              choices = names(readRDS("datasets/leukemia_combinations.Rds"))
            )
          ),
          conditionalPanel(
            condition = "input.view == 'ribbon' & input.dataset == 'epn'",
            selectInput("g_comp_epn", label = "Select two groups to compare",
              choices = names(readRDS("datasets/epn_combinations.Rds"))
            )
          ),
          conditionalPanel(
            condition = "input.view == 'ribbon' & input.dataset == 'preALL'",
            selectInput("g_comp_preALL", label = "Select two groups to compare",
              choices = names(readRDS("datasets/preALL_combinations.Rda"))
            )
          ),

          # input gene to view expression
          conditionalPanel(
              condition = "input.view == 'violin' | input.view == 'GeneExp_Show'",
              textInput("Gene", "Enter Gene Name to show Expression:", value = ""),
          ),
          conditionalPanel(
            condition = "input.view == 'GeneExp_Show'",
            sliderInput("Filter","Filter cells on UMAP with Percentage based on Max Gene Expression:", min = 0, max = 100, value = 0)
          ),

          # add statistical comparisons
          conditionalPanel(
            condition = "(input.view == 'violin' | input.view == 'v_path') & input.groupBy == 'CellType' & input.dataset == 'leukemia'",
            selectInput("singlecell_comp_ct_l", label = "Choose a Statistical Comparison:",
              choices = c("No Comparisons" = "none",
              "Compare Against B/My MPAL Blasts (with Wilcox T-test)" = "B MPAL blast",
              "Compare Against T/My MPAL Blasts (with Wilcox T-test)" = "T MPAL blast",
              "Compare Against AML Blasts (with Wilcox T-test)" = "AML blast",
              "Compare Against B-ALL Blasts (with Wilcox T-test)" = "BALL blast",
              "Compare Against T-ALL Blasts (with Wilcox T-test)" = "TALL blast",
              "Compare Against NK/T Cells (with Wilcox T-test)" = "NK/T",
              "Compare Against B Cells (with Wilcox T-test)" = "B",
              "Compare Against Erythroid Cells (with Wilcox T-test)" = "eryth",
              "Compare Against T Cells (with Wilcox T-test)" = "T",
              "Compare Against Progenitor Cells (with Wilcox T-test)" = "progenitor",
              "Compare Against Monocytes (with Wilcox T-test)" = "monocyte",
              "Compare All (with One-Way ANOVA)" = "one_anova")
            )
          ),
          conditionalPanel(
            condition = "(input.view == 'violin' | input.view == 'v_path') & input.groupBy == 'Group' & input.dataset == 'leukemia'",
            selectInput("singlecell_comp_g_l", label = "Choose a Statistical Comparison:",
              choices = c("No Comparisons" = "none",
              "Compare Against B/My MPAL (with Wilcox T-test)" = "BMPAL",
              "Compare Against T/My MPAL(with Wilcox T-test)" = "TMPAL",
              "Compare Against AML (with Wilcox T-test)" = "AML",
              "Compare Against B-ALL (with Wilcox T-test)" = "BALL",
              "Compare Against T-ALL (with Wilcox T-test)" = "TALL",
              "Compare Against Control (with Wilcox T-test)" = "Con",
              "Compare All (with One-Way ANOVA)" = "one_anova")
            )
          ),
          conditionalPanel(
            condition = "(input.view == 'violin' | input.view == 'v_path') & input.groupBy_wt == 'CellType' & input.dataset == 'wilmsTumor'",
            selectInput("singlecell_comp_ct_wt", label = "Choose a Statistical Comparison:",
              choices = c("No Comparisons" = "none",
              "Compare Against Tumor Cells (with Wilcox T-test)" = "Tumor",
              "Compare Against Endothelial Cells (with Wilcox T-test)" = "Endothelial",
              "Compare Against Immune Cells (with Wilcox T-test)" = "Immune Cells",
              "Compare Against SMC/Fibro Cells (with Wilcox T-test)" = "SMC/Fibro",
              "Compare All (with One-Way ANOVA)" = "one_anova")
            )
          ),
          conditionalPanel(
            condition = "(input.view == 'violin' | input.view == 'v_path') & input.groupBy_wt == 'Group' & input.dataset == 'wilmsTumor'",
            selectInput("singlecell_comp_g_wt", label = "Choose a Statistical Comparison:",
              choices = c("No Comparisons" = "none",
              "Compare Against Anaplastic (with Wilcox T-test)" = "Anaplastic",
              "Compare Against Favorable (with Wilcox T-test)" = "Favorable",
              "Compare All (with One-Way ANOVA)" = "one_anova")
            )
          ),
          conditionalPanel(
            condition = "(input.view == 'violin' | input.view == 'v_path') & input.groupBy_rms == 'Group' & input.dataset == 'rms'",
            selectInput("singlecell_comp_g_rms", label = "Choose a Statistical Comparison:",
              choices = c("No Comparisons" = "none","Compare with Wilcox T-test" = "ERMS")
            )
          ),
          conditionalPanel(
            condition = "(input.view == 'violin' | input.view == 'v_path') & input.groupBy_rms == 'CellType' & input.dataset == 'rms'",
            selectInput("singlecell_comp_ct_rms", label = "Choose a Statistical Comparison:",
              choices = c("No Comparisons" = "none", "Compare Against ERMS Tumor cells (with Wilcox T-test)" = "ERMS Tumor", "Compare Against ARMS Tumor cells (with Wilcox T-test)" = "ARMS Tumor",
                "Compare Against Endothelial cells (with Wilcox T-test)" = "Endothelial", "Compare Against Monocyte cells (with Wilcox T-test)" = "Monocyte",
                "Compare Against T-cells (with Wilcox T-test)" = "Lymphocyte(T-cell)", "Compare Against Epithelial cells (with Wilcox T-test)" = "Epithelial",
                "Compare All (with One-Way ANOVA)" = "one_anova"
              )
            )
          ),
          conditionalPanel(
            condition = "(input.view == 'violin' | input.view == 'v_path') & input.groupBy_epn == 'Group' & input.dataset == 'epn'",
            selectInput("singlecell_comp_g_epn", label = "Choose a Statistical Comparison:",
              choices = c("No Comparisons" = "none", "Compare Against PF-A" = "PF-A", "Compare Against PF-B" = "PF-B", "Compare Against PF-SE" = "PF-SE",
                "Compare Against SP-MPE" = "SP-MPE", "Compare Against ST-RELA" = "ST-RELA", "Compare Against ST-YAP1" = "ST-YAP1", "Compare All (with One-Way ANOVA)" = "one_anova"
              )
            )
          ),
          conditionalPanel(
            condition = "(input.view == 'violin' | input.view == 'v_path') & input.groupBy_epn == 'CellType' & input.dataset == 'epn'",
            selectInput("singlecell_comp_ct_epn", label = "Choose a Statistical Comparison:",
              choices = c(
                "No Comparisons" = "none",
                "Compare Against Immune cells/Macrophages" = "Immune cells/Macrophages", "Compare Against Oligodendrocyte precursor cells" = "Oligodendrocyte precursor cells",
                "Compare Against Oligodendrocytes" = "Oligodendrocytes", "Compare Against T-cells" = "T-cells", "Compare Against Tumro Cells" = "Tumor", "Compare All (with One-Way ANOVA)" = "one_anova"
              )
            )
          ),
          conditionalPanel(
            condition = "(input.view == 'violin' | input.view == 'v_path') & input.groupBy_preALL == 'Group' & input.dataset == 'preALL'",
            selectInput("singlecell_comp_g_preALL", label = "Choose a Statistical Comparison:",
              choices = c("No Comparisons" = "none", "Compare Against Pre-B ALL (HDD)" = "Pre-B HDD",
              "Compare Against Pre-B ALL (ETV6-RUNX1)"="Pre-B ETV6-RUNX1", "Compare Against Pre-T ALL"="Pre-T",
              "Compare Against Healthy"="Healthy", "Compare All (with One-Way ANOVA)" = "one_anova"
              )
            )
          ),
          conditionalPanel(
            condition = "(input.view == 'violin' | input.view == 'v_path') & input.groupBy_preALL == 'CellType' & input.dataset == 'preALL'",
            selectInput("singlecell_comp_ct_preALL", label = "Choose a Statistical Comparison:",
              choices = c(
                "No Comparisons" = "none",
                "Compare Against Blasts"="blast", "Compare Against T/NK"="T/NK", 
                "Compare Against Erythrocytes"="Erythrocytes", "Compare Against B Cells and Monocytes"="Bcells/Mono",
                "Compare All (with One-Way ANOVA)" = "one_anova"
              )
            )
          ),

          conditionalPanel(
            condition = "(input.view == 'violin' | input.view == 'v_path') & input.groupBy_preALL == 'CellType_Specific' & input.dataset == 'preALL'",
            selectInput("singlecell_comp_cts_preALL", label = "Choose a Statistical Comparison:",
              choices = c(
                "No Comparisons" = "none",
                "Compare Against Pre-B ETV6-RUNX1 Blasts"="Pre-B ETV6-RUNX1 blast",
                "Compare Against Pre-B HDD Blasts"="Pre-B HDD blast",
                "Compare Against Pre-T Blasts"="Pre-T blast",
                "Compare All (with One-Way ANOVA)" = "one_anova"
              )
            )
          ),

          conditionalPanel(
              condition = "input.view == 'v_path' | input.view == 'u_path'",
              selectizeInput("Path", label = "Select a pathway to show its enrichment:",
              choices = readRDS("datasets/path_names.rds"))
          ),
          conditionalPanel(
            condition = "input.view == 'ribbon'",
            selectizeInput("Path_ribbon", label = "Select pathway(s) to show enrichment:",
              choices = readRDS("datasets/path_names.rds"), multiple = TRUE
            )
          ),
          actionButton("updatePlot","Produce Plot with Selected Options")
        ),
        menuItem("ImmuneCell", tabName = "immuneCell", startExpanded = FALSE, selected = FALSE,
          textInput("Gene_hca","Enter Gene Name:", value = ""),
          radioButtons(
            "view_hca", "Select Viewing Format:",
            choices = c("UMAP" = "umap","Violin Plot with Cell Types" = "violin")
          ),
          actionButton('updatePlot_hca','Press to Update Plot')            
        ),
        menuItem("Biomarkers", tabName = "biomarkers", startExpanded = FALSE, selected = FALSE,
          radioButtons('bio_choose',label="View a pre-loaded biomarker set or test your own:", choices=c("Pre-Loaded Biomarkers" = "pre_bio")),
          radioButtons('view_bio', label = 'Viewing Format:', choices = c("UMAP" = "UMAP_bio","Violin" = "violin_bio")),
          conditionalPanel(condition = "input.dataset == 'leukemia' & input.bio_choose == 'pre_bio'",
            radioButtons('gs_bio', label = 'Biomarker Set:',
            choices = c("B-MPAL" = "bm_bio", "T-MPAL" = 'tm_bio', "AML" = 'a_bio', "T-ALL" = 't_bio', 'B-ALL' = 'b_bio'))
          ),
          conditionalPanel(condition = "input.dataset == 'wilmsTumor' & input.bio_choose == 'pre_bio'",
            radioButtons('gs_bio_w', label = "Biomarker Set:",
            choices = c("Favorable" = "Favorable", "Anaplastic" = "Anaplastic"))
          ),
          conditionalPanel(condition = "input.dataset == 'rms' & input.bio_choose == 'pre_bio'",
            radioButtons('gs_bio_rms', label = "Biomarker Set:",
            choices = c("Alveolar RMS" = "ARMS_bio", "Embryonal RMS" = "ERMS_bio"))
          ),
          conditionalPanel(condition = "input.dataset == 'epn' & input.bio_choose == 'pre_bio'",
            radioButtons('gs_bio_epn', label = "Biomarker Set:",
              choices = c("Aggressive" = "Aggressive_bio", "Favorable" = "Favorable_bio")
            )
          ),
          conditionalPanel(
            condition = "input.dataset == 'preALL' & input.bio_choose == 'pre_bio'",
            radioButtons(
              "gs_bio_preALL", label = "Biomarker Set:",
              choices = c("Pre-B ALL (HDD)" = "Pre.B.HDD_bio", "Pre-B ALL (ETV6-RUNX1)" = "Pre.B.ETV6.RUNX1_bio",
              "Pre-T ALL"="Pre.T_bio")
            )
          ),

          conditionalPanel(condition = "bio_choose == 'user_bio'",
            textAreaInput("user_bioSet", label = "Input biomarker set")
          ),
          conditionalPanel(condition = "input.dataset != 'rms' & input.dataset != 'epn' & input.dataset != 'preALL'",
            selectInput("groupBy_bio", label = 'Data Aspect to Show:',
            c("Clusters" = "clusters",
                    "Cell Type" = "CellType",
                    "Disease Type" = "Group")),
          ),
          conditionalPanel(condition = "input.dataset == 'rms'",
            selectInput("groupBy_bio_rms", label = "Data Aspect to Show:",
              c("Clusters" = "clusters",
                "Disease Type" = "Group",
                "Sample ID" = "sample",
                "Cell Type" = "CellType",
                "Age" = "age",
                "gender" = "gender",
                "Disease Stage" = "stage_group",
                "Risk Stratum" = "risk_stratum",
                "Time Point" = "timePoint",
                "Sample Site" = "sampleSite",
                "Cell Type Overall" = "CellType_Overall")
            )
          ),
          conditionalPanel(condition = "input.dataset == 'epn'",
            selectInput("groupBy_bio_epn", label = "Data Aspect to Show:",
              c("Clusters"="clusters","Disease Type"="Group","Sample ID"="sample","Cell Type"="CellType","Age"="age",
                "Gender"="gender","Therapy"="therapy","Patient Outcome"="outcome", "Broad Group" = "GroupBroad", "Prognosis" = "Prognosis",
                "Time Point"="timePoint","Sample Site"="sampleSite","Cell Type Overall" = "CellType_Overall", "Fine Group" = "GroupFine", 
                "Overall Survival (Years)" = "OverallSurvival"
              )
            )
          ),
          conditionalPanel(
            condition = "input.dataset == 'preALL'",
            selectInput("groupBy_bio_preALL", label = "Data Aspect to Show:",
              c("Disease Type" = "Group", "Cell Type (Overall)" = "CellType_Overall", "Prognosis" = "Prognosis",
              "Cell Cycle Phase"="CellCyclePhase")
            )
          ),

          actionButton('updatePlot_bio', 'Press to Update Biomarker Selection')
        ),
        menuItem("BulkExpression", tabName = "bulkExp", startExpanded = FALSE, selected = FALSE,
          conditionalPanel(
            condition = "input.dataset != 'rms' & input.dataset != 'epn' & input.dataset != 'preALL'",
            selectInput("gene_select", label = "Gene Input Options:", c("Single Gene" = "single_gene","Gene Set" = "gene_set")),
            selectInput("exp_surv", label = "View Expression or Survival in Bulk Samples", c("Expression" = "exp", "Survival" = "surv")),
            conditionalPanel(
              condition = "input.gene_select == 'single_gene'",
              textInput("Gene_bulkExp","Enter Gene Name:", value = "")
            ),
            conditionalPanel(
              condition = "input.gene_select == 'gene_set'",
              textAreaInput("GeneSet_bulkExp", "Enter Gene Names, One per Line:", value = ""),
              textInput("GeneSet_name", "Name of Gene Set:", value = "user_set")
            ),
            conditionalPanel(condition = "input.dataset == 'leukemia' & input.exp_surv == 'exp'",
              selectInput("bulk_stats_l", label = "Choose which statistical comparisons to perform:",
                c("No Comparisons" = "none",
                    "Compare Against B/My MPAL (with Wilcox T-test)" = "to_bm",
                    "Compare Against T/My MPAL (with Wilcox T-test)" = "to_tm",
                    "Compare Against AML (with Wilcox T-test)" = "to_a",
                    "Compare Against Pre B-ALL (with Wilcox T-test)" = "to_preb",
                    "Compare All (with One-Way ANOVA)" = "one_anova"))
            ),
            conditionalPanel(condition = "input.dataset == 'wilmsTumor' & input.exp_surv == 'exp'",
              selectInput("bulk_stats_w", label = "Choose which statistical comparisons to perform:",
                c("No Comparisons" = "none", "Compare with Wilcox T-test" = "t")
              )
            ),
            conditionalPanel(condition = "input.dataset == 'leukemia' & input.exp_surv == 'surv'",
              selectInput("surv_cut_l", label = "Choose Method to Define Groups", c("median" = "median","cutP" = "cutP")),
              selectInput("surv_dataset_l", label = "Choose leukemia subtype for survival analysis:",
                c("B/My MPAL" = "B/M MPAL", "T/My MPAL" = "T/M MPAL", "AML" = "AML", "Pre-B ALL" = "B-Precursor", "All Subytpes" = "All")
              )
            ),
            conditionalPanel(condition = "input.dataset == 'wilmsTumor' & input.exp_surv == 'surv'",
              selectInput("surv_cut_w", label = "Choose Method to Define Groups", c("median" = "median","cutP" = "cutP")),
              selectInput("surv_dataset_w", label = "Choose Wilms Tumor subtype for survival analysis:",
                c("Favorable" = "Favorable", "Anaplstic" = "Anaplastic", "All Subtypes" = "All")
              )
            ),
            actionButton("updatePlot_bulkExp","Press to Update Plot with User Input")
          )
        ),
        menuItem("Prediction", tabName = "predict", startExpanded = FALSE, selected = FALSE,
          conditionalPanel(
            condition = "input.dataset == 'leukemia'",
            fileInput('user_in', "Input Sample(s) Expression of the Disease Biomarkers as a CSV (see Github Wiki for Formatting Instructions):",
              multiple = FALSE, accept = ".csv"),
            actionButton("go_predict", "Predict Leukemia Subtype of Sample(s)")
          )
         )
        )
      ),

    
    body = dashboardBody(
      #shinyDashboardThemes(theme = "grey_light"),
      tags$head(tags$link(rel="icon", type="image/png", href="PedSCAtlas_logo.png")),
      theme_grey_light,
      fluidRow(
        column(width = 3,
          box(title = "Messages", width = NULL, solidHeader = FALSE,
            img(src="PedSCAtlas_logo.png", align = "center", style = "height:100px"),
            hr(),
            span(textOutput("dataset_chosen"), style="color:#4f5379"),
            br(),
            span(textOutput('error_message'), style="color:#d90429")
          ),
          box(title = "Information", width = NULL, solidHeader = FALSE,
          span("Welcome to the Pediatric Single Cell Cancer Atlas! Please visit our Github repository to post issues/questions and view the user guide."), hr(),
          socialButton("https://github.com/bhasin-lab/PedSCAtlas", icon("github"))
        )),
        column(width = 9,
          box(title = "Plots", width = NULL, solidHeader = TRUE,
            uiOutput("plot_atlas")
          ),  
          box(title = "Tables", width = NULL, solidHeader = TRUE,
            dataTableOutput("pred_table")
          )
        )
              
      )
    ),

    controlbar = dashboardControlbar(
      id = "version",
      skin = "light",
      br(),
      controlbarItem("Version", "  Updated as of 11.04.2022   ")
    ),

    title = "PedSCAtlas",

    footer = dashboardFooter(
      # HTML to add links and change color
      left = tags$a(href="http://www.bhasinlab.org/","Bhasin Systems Biomedicine Lab at Emory University"),

      right = tags$a(href="https://www.choa.org/medical-services/cancer-and-blood-disorders","Aflac Cancer & Blood Disorders Center at Children's Healthcare of Atlanta")
    )
   ),
   server = function(input, output) { 
    ## choose dataset
    choose_dataset <- modalDialog(
      title = p("Before starting your analysis, select a dataset below. Then, choose a module on the tabs at the left side of the page to analyze the selected dataset.
        This popup will disappear after the datasets are loaded. Depending on the size of the dataset chosen, this may take a few minutes.", style = "color:navy"),
      size = "m",
      selectInput("dataset", label = "Select a pediatric cancer dataset to analyze",
        choices = c("Acute Leukemias" = "leukemia", 
          "Wilms Tumor" = "wilmsTumor", 
          "Rhabdomyosarcoma (RMS)" = "rms", 
          "Ependymoma (EPN)" = "epn",
          "Pre-B and Pre-T ALL" = "preALL"
          )
      ),
      easyClose = F,
      footer = tagList(
        actionButton("start", "Start Analysis")
      )
    )

    ## show the modal on start up
    showModal(choose_dataset) 

    # reactiveVal to start the chosen dataset
    atlas.data = reactiveVal()


    loading = addLoader$new("start", type = "spinner", color = "#d90429")

    observeEvent(input$start, {
      loading$show()
    })

    observeEvent(input$start, {
      
      if (input$dataset == "leukemia") {
        data = readRDS("datasets/acute_leukemias.Rda")
      } else if (input$dataset == "wilmsTumor") {
        data = readRDS("datasets/wilms_tumor.Rda")
      } else if (input$dataset == "rms") {
        data = readRDS("datasets/rms.Rda")
      } else if (input$dataset == "epn") {
        data = readRDS("datasets/epn.Rda")
      } else if (input$dataset == "preALL") {
        data = readRDS("datasets/preALL.Rda")
      }

      load(file="counter.Rdata")
      counter = counter + 1
      save(counter, file="counter.Rdata")

      if (input$dataset != "test") {
        data$coords_hca <- readRDS(file = "datasets/HCA/hca_umap.Rda")
        data$geneExp_hca <- readRDS(file = "datasets/HCA/hca_geneExp.Rda")

        atlas.data(data)
      
        print("dataset loaded ...")
      }
 
      removeModal()
      
    })

    # generate plots based on user input
    output$plot_user <- renderPlot({  
      input$updatePlot

      loaded.data = atlas.data()
      coords = loaded.data$coords
      geneExp = loaded.data$geneExp
      path_data = loaded.data$path

      if (input$dataset == "rms") {
        group_by = input$groupBy_rms
      } else if (input$dataset == "wilmsTumor") {
        group_by = input$groupBy_wt
      } else if (input$dataset == "epn") {
        group_by = input$groupBy_epn
      } else if (input$dataset == "preALL") {
        group_by = input$groupBy_preALL
      } else {
        group_by = input$groupBy
      }

      isolate({
        coords["Gene_Expression"] <- 0
        gn <- "No Gene Chosen"
        gnExp_filter <- 0

        if ((input$view == "GeneExp_Show" | input$view == "violin") & nchar(input$Gene) > 0) {
          gnFound = toupper(input$Gene) %in% rownames(geneExp)
        } else {
          gnFound = FALSE
        }       

        # assign statistical comparison
        comp = "none"
        if (input$view == "violin" | input$view == "v_path") {
          if (input$dataset == "leukemia" & input$groupBy == "CellType") comp = input$singlecell_comp_ct_l
          else if (input$dataset == "leukemia" & input$groupBy == "Group") comp = input$singlecell_comp_g_l
          else if (input$dataset == "wilmsTumor" & input$groupBy_wt == "CellType") comp = input$singlecell_comp_ct_wt
          else if (input$dataset == "wilmsTumor" & input$groupBy_wt == "Group") comp = input$singlecell_comp_g_wt
          else if (input$dataset == "rms" & input$groupBy_rms == "CellType") comp = input$singlecell_comp_ct_rms
          else if (input$dataset == "rms" & input$groupBy_rms == "Group") comp = input$singlecell_comp_g_rms
          else if (input$dataset == "epn" & input$groupBy_epn == "Celltype") comp = input$singlecell_comp_ct_epn
          else if (input$dataset == "epn" & input$groupBy_epn == "Group") comp = input$singlecell_comp_g_epn
          else if (input$dataset == "preALL" & input$groupBy_preALL == "CellType") comp = input$singlecell_comp_ct_preALL
          else if (input$dataset == 'preALL' & input$groupBy_preALL == "CellType_Specific") comp = input$singlecell_comp_cts_preALL
          else if (input$dataset == "preALL" & input$groupBy_preALL == "Group") comp = input$singlecell_comp_g_preALL
        }

        if (input$view == "GeneExp_Show" & nchar(input$Gene) > 0 & gnFound) {
          gn <- toupper(input$Gene)
          gn.ind <- which(rownames(geneExp) == gn)
          gnExp <- geneExp[gn.ind,]
          gnExp_filter <- input$Filter*max(gnExp)/100
          coords["Gene_Expression"] <- gnExp
          if (input$dataset == "leukemia") {
            coords_filt <- coords[apply(coords, 1, function(arr) any(arr %in% input$show)),]
          } else if (input$dataset == "wilmsTumor") {
            coords_filt <- coords[apply(coords, 1, function(arr) any(arr %in% input$show.wt)),]
          } else if (input$dataset == "rms") {
            coords_filt <- coords %>% filter(Group.Overall == input$show.rms)
          } else if (input$dataset == "epn") {
            coords_filt <- coords %>% filter(GroupBroad %in% input$show.epn)
          } else if (input$dataset == "preALL") {
            coords_filt <- coords %>% filter(Group %in% input$show.preALL)
          }
          
          plot = ggplot(coords_filt, aes_string(x="UMAP_1", y="UMAP_2", color="Gene_Expression")) + geom_point() 
          plot = plot + ggtitle(paste(input$Gene, " Expression")) + theme(text = element_text(size=20))

        } else if (input$view == "violin" & nchar(input$Gene) > 0 & gnFound) {
            gn <- toupper(input$Gene)
            gn.ind <- which(rownames(geneExp) == gn)
            gnExp <- geneExp[gn.ind,]
            coords["Gene_Expression"] <- gnExp
            if (input$dataset == "leukemia") {
              coords_filt <- coords[apply(coords, 1, function(arr) any(arr %in% input$show)),]
            } else if (input$dataset == "wilmsTumor") {
              coords_filt <- coords[apply(coords, 1, function(arr) any(arr %in% input$show.wt)),]
            } else if (input$dataset == "rms") {
              coords_filt <- coords %>% filter(Group.Overall == input$show.rms)
            } else if (input$dataset == "epn") {
              coords_filt <- coords %>% filter(GroupBroad %in% input$show.epn)
            } else if (input$dataset == "preALL") {
              coords_filt <- coords %>% filter(Group %in% input$show.preALL)
            }
            plot = ggplot(data = coords_filt, aes_string(x = group_by, y = "Gene_Expression", fill = group_by)) + geom_violin(trim = FALSE, scale = "width") 
            plot = plot + xlab(group_by) + ylab(paste(gn, "Expression")) + theme(text = element_text(size=20))

            # add statistical comparison
            if (comp == "none") plot = plot
            else if (comp == "one_anova") plot = plot + stat_compare_means(method = "anova")
            else {
              comp_list = getComparisons(comp, unique(coords_filt[group_by]))
              plot = plot + stat_compare_means(comparison = comp_list)
            }

            # if x axis group has more than 10 groups, don't show x axis labels
            if (nrow(unique(coords_filt[group_by])) > 10) {
              plot = plot + theme(axis.text.x=element_blank())
            }
          }

          else if (input$view == "v_path") {
            path = input$Path
            coords[path] = path_data[path]
            if (input$dataset == "leukemia") {
              coords_filt <- coords[apply(coords, 1, function(arr) any(arr %in% input$show)),]
            } else if (input$dataset == "wilmsTumor") {
              coords_filt <- coords[apply(coords, 1, function(arr) any(arr %in% input$show.wt)),]
            } else if (input$dataset == "rms") {
              coords_filt <- coords %>% filter(Group.Overall == input$show.rms)
            } else if (input$dataset == "epn") {
              coords_filt <- coords %>% filter(GroupBroad %in% input$show.epn)
            } else if (input$dataset == "preALL") {
              coords_filt <- coords %>% filter(Group %in% input$show.preALL)
            }

            #plot <- plot_ly(data = coords_filt, x = ~get(input$groupBy), y = ~get(path_title), split = ~get(input$groupBy), type = 'violin', box = list(visible = T), meanline = list(visible = T))
            #plot <- plot %>% layout(title = path_title, xaxis = list(title = input$groupBy), yaxis = list(title = path_title))
            plot = ggplot(data = coords_filt, aes_string(x = group_by, y = path, fill = group_by)) + geom_violin() 
            plot = plot + xlab(group_by) + ylab(paste(path, "Enrichment")) + theme(text = element_text(size=20))

            # add statistical comparison
            if (comp == "none") plot = plot
            else if (comp == "one_anova") plot = plot + stat_compare_means(method = "anova")
            else {
              comp_list = getComparisons(comp, unique(coords_filt[group_by]))
              plot = plot + stat_compare_means(comparison = comp_list)
            }

            # if x axis group has more than 10 groups, don't show x axis labels
            if (nrow(unique(coords_filt[group_by])) > 10) {
              plot = plot + theme(axis.text.x=element_blank())
            }
          }

          else if (input$view == "u_path") {
            path = input$Path
            coords[path] = path_data[path]
            if (input$dataset == "leukemia") {
              coords_filt <- coords[apply(coords, 1, function(arr) any(arr %in% input$show)),]
            } else if (input$dataset == "wilmsTumor") {
              coords_filt <- coords[apply(coords, 1, function(arr) any(arr %in% input$show.wt)),]
            } else if (input$dataset == "rms") {
              coords_filt <- coords %>% filter(Group.Overall == input$show.rms)
            } else if (input$dataset == "epn") {
              coords_filt <- coords %>% filter(GroupBroad %in% input$show.epn)
            } else if (input$dataset == "preALL") {
              coords_filt <- coords %>% filter(Group %in% input$show.preALL)
            }

            plot = ggplot(coords_filt, aes_string(x="UMAP_1",y="UMAP_2",color=path)) + geom_point() #
            plot = plot + ggtitle(paste(path, "Enrichment")) + theme(text = element_text(size=20))

          } else if (input$view == "ribbon") {
            req(length(input$Path_ribbon) > 1)
            if (input$dataset == "wilmsTumor") {
              g1 = "Favorable"
              g2 = "Anaplastic"
              group_ribbon = "Group"
            } else if (input$dataset == "rms") {
              g1 = "ARMS"
              g2 = "ERMS"
              group_ribbon = "Group.Overall"
            } else if (input$dataset == "epn") {
              vs_ind = gregexpr("vs",input$g_comp_epn)[[1]][1]
              g1 = substr(input$g_comp_epn, 1, vs_ind-1)
              g2 = substr(input$g_comp_epn, vs_ind+2, nchar(input$g_comp_epn))
              group_ribbon = "Group"
            } else if (input$dataset == "preALL") {
              vs_ind = gregexpr("vs",input$g_comp_preALL)[[1]][1]
              g1 = substr(input$g_comp_preALL, 1, vs_ind-1)
              g2 = substr(input$g_comp_preALL, vs_ind+2, nchar(input$g_comp_preALL))
              group_ribbon = "Group"
            }
            else {
              vs_ind = gregexpr("vs",input$g_comp)[[1]][1]
              g1 = substr(input$g_comp, 1, vs_ind-1)
              g2 = substr(input$g_comp, vs_ind+2, nchar(input$g_comp))
              group_ribbon = "Group"
            }
            
            path_names = input$Path_ribbon

            plot = ribbonPaths(path_data, coords, group_ribbon, g1, g2, path_names)

          }
          return(plot)
      })
      

    })

    # generate plots for UMAP with groups
    output$plot_umap <- renderPlot({
      input$updatePlot

      loaded.data = atlas.data()
      coords = loaded.data$coords

      isolate({
        if (input$dataset == "leukemia") {
          coords_filt <- coords[apply(coords, 1, function(arr) any(arr %in% input$show)),]
          group_by = input$groupBy
        } else if (input$dataset == "wilmsTumor") {
          coords_filt <- coords[apply(coords, 1, function(arr) any(arr %in% input$show.wt)),]
          group_by = input$groupBy_wt
        } else if (input$dataset == "rms") {
          coords_filt <- coords[apply(coords, 1, function(arr) any(arr %in% input$show.rms)),]
          group_by = input$groupBy_rms
        } else if (input$dataset == "epn") {
          coords_filt <- coords %>% filter(GroupBroad %in% input$show.epn)
          group_by = input$groupBy_epn
        } else if (input$dataset == "preALL") {
          coords_filt <- coords %>% filter(Group %in% input$show.preALL)
          group_by = input$groupBy_preALL
        }
        plot = ggplot(coords_filt, aes_string(x="UMAP_1", y="UMAP_2", color=group_by)) + geom_point() # + scale_fill_brewer(palette="Spectral")
        plot = plot + ggtitle(paste("UMAP with", group_by)) + theme(text = element_text(size=20))
        if (group_by == 'Group') {       
          plot = plot + scale_color_brewer(palette = 'Dark2')
        }
      })
      return(plot)
    })

    # generate plots for ImmuneCell tab
    output$plot_hca <- renderPlot({
      input$updatePlot_hca

      loaded.data = atlas.data()

      coords_hca = loaded.data$coords_hca
      geneExp_hca = loaded.data$geneExp_hca

      isolate({
        coords_hca["Gene_Expression"] <- 0
        gn_hca <- "No Gene Chosen"
        if (nchar(input$Gene_hca) > 0) {
          gn_hca <- toupper(input$Gene_hca)
          gn.ind_hca <- which(rownames(geneExp_hca) == gn_hca)
          gnExp_hca <- geneExp_hca[gn.ind_hca,]
          coords_hca["Gene_Expression"] <- gnExp_hca
          
          if (input$view_hca == "umap") {
            plot = ggplot(coords_hca, aes_string(x="UMAP_1", y="UMAP_2", color="Gene_Expression")) + geom_point()
            plot + ggtitle(paste(gn_hca, "Expression")) + scale_color_gradient(low = "grey", high = "purple")
          } else {
            plot = ggplot(coords_hca, aes_string(x="CellType", y="Gene_Expression", fill="CellType")) +
              geom_violin(scale = "width", trim = FALSE) + ggtitle(paste(gn_hca, "Expression"))
            plot + xlab("HCA Labeled Cell Type Annotation") + theme(axis.text.x = element_blank()) + theme(text = element_text(size=20))
          }
        
        }
      })
    })

    # generate pltos for Biomarkers tab
    output$plot_bio <- renderPlot({
      input$updatePlot_bio
      loaded.data = atlas.data()
      coords_bio = loaded.data$coords_bio

      if (input$dataset == "rms") {
        group_bio = input$groupBy_bio_rms
      } else if (input$dataset == "epn") {
        group_bio = input$groupBy_bio_epn
      } else if (input$dataset == "preALL") {
        group_bio = input$groupBy_bio_preALL
      } else {
        group_bio = input$groupBy_bio
      }

      isolate({
        if (input$dataset == "leukemia") {gs_bio = input$gs_bio}
          else if (input$dataset == "wilmsTumor") {
            gs_bio = input$gs_bio_w
            coords_bio$Anaplastic = unlist(coords_bio$Anaplastic)
            coords_bio$Favorable = unlist(coords_bio$Favorable)
          } else if (input$dataset == "rms") {
            gs_bio = input$gs_bio_rms
          } else if (input$dataset == "epn") {
            gs_bio = input$gs_bio_epn
          } else if (input$dataset == "preALL") {
            gs_bio = input$gs_bio_preALL
          }

        if (gs_bio == "tm_bio") {
            title = "T-MPAL"
          } else if (gs_bio == "bm_bio") {
            title = "B-MPAL"
          } else if (gs_bio == "a_bio") {
            title = "AML"
          } else if (gs_bio == "b_bio") {
            title = "BALL"
          } else if (gs_bio == "t_bio") {
            title = "TALL"
          } else if (gs_bio == "Favorable") {
            title = "Favorable"
          } else if (gs_bio == "Anaplastic") {
            title = "Anaplastic"
          } else if (gs_bio == "ARMS_bio") {
            title = "ARMS"
          } else if (gs_bio == "ERMS_bio") {
            title = "ERMS"
          } else if (gs_bio == "Aggressive_bio") {
            title = "Aggressive"
          } else if (gs_bio == "Favorable_bio") {
            title = "Favorable"
          } else if (gs_bio == "Pre.B.ETV6.RUNX1_bio") {
            title = "Pre-B ALL (ETV6-RUNX1)"
          } else if (gs_bio == "Pre.B.HDD_bio") {
            title = "Pre-B ALL (HDD)"
          } else if (gs_bio == "Pre.T_bio") {
            title = "Pre-T ALL"
          }

        if(input$view_bio == "violin_bio") {
          plot = ggplot(coords_bio, aes_string(x=group_bio, y=gs_bio, fill=group_bio)) + geom_violin() 
          plot = plot + ggtitle(paste(title,"Biomarker Score and",group_bio)) + theme(text = element_text(size=20))
        } else {
          plot = ggplot(coords_bio, aes_string(x="UMAP_1", y="UMAP_2", color=gs_bio)) + geom_point() 
          plot = plot + ggtitle(paste(title,"Biomarker Score")) + theme(text = element_text(size=20))
          plot = plot + scale_color_gradient2(low = "red", mid='grey', high = "blue")
        }
        return(plot)
      })
    })

    # generate plots for BulkExpression tab
    output$plot_bulk <- renderPlot({
      input$updatePlot_bulkExp

      loaded.data = atlas.data()
      bulkExp = loaded.data$bulk_data
      rownames(bulkExp) = toupper(rownames(bulkExp))
      bulkMD = loaded.data$bulk_md

      # assign choices
      if (input$dataset == "leukemia") {
        bulk_stats = input$bulk_stats_l
        surv_cut = input$surv_cut_l
        surv_dataset = input$surv_dataset_l
      } else if (input$dataset == "wilmsTumor") {
        bulk_stats = input$bulk_stats_w
        surv_cut = input$surv_cut_w
        surv_dataset = input$surv_dataset_w
      } 
      isolate({

        # single gene expression
        if (input$gene_select == "single_gene" & input$exp_surv == 'exp') {
          if (toupper(input$Gene_bulkExp) %in% toupper(rownames(bulkExp))) {
            #browser()
            bulk_gene = bulkExp[toupper(rownames(bulkExp)) == toupper(input$Gene_bulkExp),]
            bulkMD$geneExp = bulk_gene
            p = ggplot(bulkMD, aes(x = class, y = geneExp, fill = class)) + geom_boxplot() 
            p = p + xlab("Disease Subtype") + ylab(paste("Normalized", toupper(input$Gene_bulkExp), "Expression"))
            p = p + theme(text = element_text(size=20))

            if (bulk_stats == "one_anova") {
              p = p + stat_compare_means(method = "anova")
            } else if (bulk_stats == "to_bm") {
              comparisons = list(c("B/M MPAL","T/M MPAL"),c("B/M MPAL","AML"),c("B/M MPAL","B-Precursor"))
              p = p + stat_compare_means(comparisons = comparisons)
            } else if (bulk_stats == "to_tm") {
              comparisons = list(c("T/M MPAL","B/M MPAL"),c("T/M MPAL","AML"),c("T/M MPAL","B-Precursor"))
              p = p + stat_compare_means(comparisons = comparisons)            
            } else if (bulk_stats == "to_a") {
              comparisons = list(c("AML","T/M MPAL"),c("AML","B/M MPAL"),c("AML","B-Precursor"))
              p = p + stat_compare_means(comparisons = comparisons)            
            } else if (bulk_stats == "to_preb") {
              comparisons = list(c("B-Precursor","T/M MPAL"),c("B-Precursor","AML"),c("B-Precursor","B/M MPAL"))
              p = p + stat_compare_means(comparisons = comparisons)                       
            } else if (bulk_stats == "t") {
              p = p + stat_compare_means()
            }
            # (to-do) add different dataset options here
            p + ggtitle(paste("Bulk Normalized mRNA-seq Expression of", toupper(input$Gene_bulkExp)))
          }
        }       
        # gene set enrichment
        else if (input$gene_select == "gene_set" & input$exp_surv == 'exp') {
          gs = getGenes(input$GeneSet_bulkExp)
          dont_plot = FALSE # checking for error conditions below ...

          for (i in 1:length(gs)) {
            # check that gene is in our bulk dataset
            if (!(toupper(gs[i]) %in% toupper(rownames(bulkExp)))) {
              dont_plot = TRUE
            }
          }

          # make sure there are no commas
          if (any(grep(",", gs))) {
            dont_plot = TRUE
          }

          # make sure there is more than one gene
          if (length(gs) == 1) {
            dont_plot = TRUE
          }

          if (!dont_plot) {
            bulk_gs = bulkMD
            gs = list(gs)
            names(gs) = input$GeneSet_name

            # perform gene set enrichment analysis
            en_gs = gsva(bulkExp, gs, method = "ssgsea")
            bulk_gs$enrichment = t(en_gs)

            # generate boxplot
            p = ggplot(bulk_gs, aes(x=class, y=enrichment, fill = class)) + geom_boxplot() 
            p = p + xlab("Leukemia Subtype") + ylab("Enrichment")
            p = p + ggtitle(paste("Gene Set Enrichment of", input$GeneSet_name))
            p = p + theme(text = element_text(size=20))

            # add statistical comparisons
            if (bulk_stats == "one_anova") {
              p = p + stat_compare_means(method = "anova")
            } else if (bulk_stats == "to_bm") {
              comparisons = list(c("B/M MPAL","T/M MPAL"),c("B/M MPAL","AML"),c("B/M MPAL","B-Precursor"))
              p = p + stat_compare_means(comparisons = comparisons)
            } else if (bulk_stats == "to_tm") {
              comparisons = list(c("T/M MPAL","B/M MPAL"),c("T/M MPAL","AML"),c("T/M MPAL","B-Precursor"))
              p = p + stat_compare_means(comparisons = comparisons)            
            } else if (bulk_stats == "to_a") {
              comparisons = list(c("AML","T/M MPAL"),c("AML","B/M MPAL"),c("AML","B-Precursor"))
              p = p + stat_compare_means(comparisons = comparisons)            
            } else if (bulk_stats == "to_preb") {
              comparisons = list(c("B-Precursor","T/M MPAL"),c("B-Precursor","AML"),c("B-Precursor","B/M MPAL"))
              p = p + stat_compare_means(comparisons = comparisons)            
            } else if (bulk_stats == "t") {
              p = p + compare_means()
            }
            p
          }
        } 
        # single gene with survival
        else if (input$gene_select == "single_gene" & input$exp_surv == 'surv') {
          type = surv_dataset
          if (type == "All") {
            md_sub = bulkMD
          } else {
            md_sub = bulkMD %>% filter(class == type)
          }
        
          out = mCut(toupper(input$Gene_bulkExp), bulkExp, bulkMD, cut = surv_cut, type = type)
          md_sub["group"] = cutGroups(md_sub, out)
          plotSurv(md_sub, name = input$Gene_bulkExp, type = type)
        } 
        # gene set with survival
        else if (input$gene_select == "gene_set" & input$exp_surv == 'surv') {
          gs = getGenes(input$GeneSet_bulkExp)
          dont_plot = FALSE
          for (i in 1:length(gs)) {
            if (!(gs[i] %in% rownames(bulkExp))) {
              dont_plot = TRUE
            }
          }

          if (any(grep(",", gs))) {
            dont_plot = TRUE
          }

          if (length(gs) == 1) {
            dont_plot = TRUE
          }

          if (!dont_plot) {
            type = surv_dataset
            if (type == "All") {
              md_sub = bulkMD
            } else {
              md_sub = bulkMD %>% filter(class == type)
            }
        
            out = mSet(gs, bulkExp, md_sub, cut = surv_cut, type = type)
            md_sub["group"] = cutGroups(md_sub, out)
            plotSurv(md_sub, name = input$GeneSet_name, type = type)
          }
          
        }
      })
    })

    # Generate Cancer Subtype Predictions
    output$pred_table = renderDataTable({
      req(input$sidebarItemExpanded == "Prediction")

      loaded.data = atlas.data()
      genes = loaded.data$genes
      model = loaded.data$model

      file_in <- input$user_in # user input 
      ext <- tools::file_ext(file_in$datapath)

      req(file_in)
      user.data = read.csv(file_in$datapath, header = T, row.names = 1)
      colnames(user.data)[which(colnames(user.data) == "NME1.NME2")] = "NME1-NME2"

      # remove excess genes
      user.data = user.data[colnames(user.data) %in% genes]

      # generate predictions    
      preds <- cancerPred(model, user.data, genes)

      # put probabilities in a table to be shown
      #preds.table = data.frame("Sample" = rownames(user.data), "Prediction" = preds)
      preds.table = data.frame("Sample" = rownames(user.data))
      preds.table = cbind(preds.table, preds)
      rownames(preds.table) = NULL

      #colnames(preds.table) = c("Sample","Prob. of AML", "Prob. of B-ALL", "Prob. of B/M MPAL", "Prob. of T-ALL", "Prob. of T/M MPAL")
      colnames(preds.table) = c("Sample","Predicted Subtype")
      return(preds.table)

    })

    ## generate plot to show on dashboardBody
    output$plot_atlas <- renderUI({
      # retrieve chosen dataset
      loaded.data = atlas.data()

      req(input$sidebarItemExpanded)
      req(input$sidebarItemExpanded != 'Prediction')

      # for PediatricCancer Tab:
      if (input$sidebarItemExpanded == 'SingleCell') {
        req(input$updatePlot)
        if (input$dataset == "wilmsTumor" & input$groupBy == "state") {
          p = ggplot() + theme_void()
        }

        # pre-made html outputs
        else if (input$view == "UMAP") {
          if(length(loaded.data) > 0) {
            isolate({
              #p = getPage()
              #if (input$downloadPlot)
              p = shinycssloaders::withSpinner(plotOutput("plot_umap", height="500px"))
            })
          }
        }
        # generate plots based on user input
        else {
          p = shinycssloaders::withSpinner(plotOutput("plot_user", height="500px"))
        }
      }

      # for ImmuneCell Tab:
      else if (input$sidebarItemExpanded == "ImmuneCell") {
        req(input$updatePlot_hca)
        p = shinycssloaders::withSpinner(plotOutput("plot_hca", height="500px"))
      }

      # for Biomarkers Tab:
      else if (input$sidebarItemExpanded == "Biomarkers") {
        req(input$updatePlot_bio)
        p = shinycssloaders::withSpinner(plotOutput("plot_bio", height="500px"))
      }

      # for BulkExpression Tab:
      else if (input$sidebarItemExpanded == "BulkExpression" & input$dataset != "rms") {
        req(input$updatePlot_bulkExp)
        p = shinycssloaders::withSpinner(plotOutput("plot_bulk", height="500px"))
      }

      else  {
        p = ggplot() + theme_void()
      }

      return(p)
    })
    
   ## generate informative error messages
    output$error_message <- renderText({
      input$updatePlot

      req(input$sidebarItemExpanded)

      loaded.data = atlas.data()
      geneExp = loaded.data$geneExp

      coords_hca = loaded.data$coords_hca
      geneExp_hca = loaded.data$geneExp_hca

      genes = loaded.data$genes

      if (input$sidebarItemExpanded == "SingleCell") {
        isolate({
          req(input$updatePlot)
          gene_found <- length(which(rownames(geneExp) == toupper(input$Gene)))

          if (input$Gene == "" & input$view == "GeneExp_Show") {
            "Please enter a gene to view expression on UMAP."
          } else if(input$Gene != "" & gene_found == 0 & input$view != 'UMAP') {
            "The gene entered was not found in the expression data."
          } else if(input$view == "violin" & nchar(input$Gene) == 0) {
            "Please enter gene to view violin plot."
          } else if (input$dataset == "wilmsTumor" & input$groupBy == "state") {
            "This aspect is not available for this dataset. Please choose something else to group the plots by."
          } else if (input$view == "ribbon" & length(input$Path_ribbon) < 2) {
            "Please select at least two pathways to view a ribbon plot."
          } else {
            ""
          }
        })
      } else if (input$sidebarItemExpanded == 'ImmuneCell') {
          isolate({
            req(input$updatePlot_hca)
            gene_found_hca <- length(which(rownames(geneExp_hca) == toupper(input$Gene_hca)))
            if(input$Gene_hca == "") {
              "Please enter a gene to view expression in Healthy Immune Cells."
            } else if(input$Gene_hca != "" & gene_found_hca == 0) {
              "The gene entered was not found in the expression data."
            }
          })
      } else if (input$sidebarItemExpanded == "Biomarkers") {
        req(input$updatePlot_bio)
        if (input$dataset == "leukemia") {
          gs_bio = input$gs_bio
        } else if (input$dataset == "wilmsTumor") {
          gs_bio = input$gs_bio_w
        } else if (input$dataset == "rms") {
          gs_bio = input$gs_bio_rms
        } else if (input$dataset == "epn") {
          gs_bio = input$gs_bio_epn
        } else if (input$dataset == "preALL") {
          gs_bio = input$gs_bio_preALL
        }

        if(gs_bio == "bm_bio") {
          "B-MPAL Blast Biomarker Gene Set: S100A16, ATF3, GALNT2, KCNK12, MTRN2L12, HBEGF, CDKN1A, IFRD1, PLIN2, CD81, NR4A1, DDIT3, RFLNB, CBX4, UBE2S"
        } else if(gs_bio == "tm_bio") {
          "T-MPAL Blast Biomarker Gene Set: KDM5B, SMYD3, RAMP1, ZNF385D, MGLL, PTGER4, CD81, TUT1, NFE2, UBE2S, NME1-NME2"
        } else if(gs_bio == "a_bio") {
          "AML Blast Biomarker Gene Set: CDM5B, RNASEK, NM2, TUT1"
        } else if(gs_bio == "t_bio") {
          "T-ALL Blast Biomarker Gene Set: HES4, RCC1, CHI3L2, GALNT2, MATR3.1, LAT, RNASEK, NME2, TASP1"
        } else if (gs_bio == "b_bio") {
          "B-ALL Blast Biomarker Gene Set: S100A16, SPRY1, HBEGF, CDKN1A, NPY, IFRD1, PLIN2, TLE1, NR4A1, NEIL1, RASD1, CBX4"
        } else if (gs_bio == "Favorable") {
          "Favorable Tumor Biomarker Gene Set: LRFN5, HMCN1, CNTNAP2, PDGFC, ROBO2, FAT3, CACHD1, ERBB4, MECOM, UTRN"
        } else if (gs_bio == "Anaplastic") {
          "Anaplastic Tumor Biomarker Gene Set: PAX8, SLIT2, ITPR1, GPC5, DLG2, WT1, ITGA8, UNC5C, TRABD2B, TSIX, SLIT3, MEIS2, DPP6, MEIS1"
        } else if (gs_bio == "ARMS_bio") {
          "Alveolar Biomarker Gene Set: DIO2, DIO2-AS1, NOS1, KCNIP4, SORCS1, DLGAP1, RYR3, ANKS1B, NCKAP5, CSMD1"
        } else if (gs_bio == "ERMS_bio") {
          "Embryonal RMS Biomarker Gene Set: NRG1, HMGA2, RBMS3, ROBO1, LHFP, FAT3, CDH11, PSD3, SGCD, PRKG1"
        } else if (gs_bio == "Aggressive_bio") {
          "Aggressive Tumor Biomarker Gene Set: TNC, HNRNPA1, CCND1, STMN1, ELN, MEIS1, TUBB2B, LMO3, NREP, NTRK2"
        } else if (gs_bio == "Favorable_bio") {
          "Favorable Tumor Biomarker Gene Set: C9orf117, SPAG8, FBLN5, C9orf24, FRMPD2, WLS, BSCL2, CAPS, DNAAF1, AK1"
        } else if (gs_bio == "Pre.B.ETV6.RUNX1_bio") {
          "Pre-B ALL (ETV6-RUNX1) Blast Biomarker Gene Set: TERF2, CD24, SOCS2, UHRF1, CD74, HLA-DRA, HLA-DPB1, RCSD1, VPREB3, LAPTM5"
        } else if (gs_bio == "Pre.B.HDD_bio") {
          "Pre-B ALL (HDD) Blast Biomarker Gene Set: VPREB1, CD9, S100A16, TCF4, STIM2, LCN6, TCL1A, SOX4,ALOX5AP, KLF6"
        } else if (gs_bio == "Pre.T_bio") {
          "Pre-T ALL Blast Biomarker Gene Set: CHI3L2, FXYD2, TRBC2, TRBC1, CD3D, CD1E, MAL, TCBTB2, ALDH1A2, ITM2A"
        }

        # (to-do) add different dataset options here
      } else if(input$sidebarItemExpanded == "BulkExpression" & input$dataset != 'rms' & input$dataset != "epn" & input$dataset != "preALL") {
        req(input$updatePlot_bulkExp)
        bulkExp = loaded.data$bulk_data

        isolate({
          if (input$gene_select == "single_gene") {
            gene_found_bulk <- length(which(toupper(rownames(bulkExp)) == toupper(input$Gene_bulkExp)))
            if (input$Gene_bulkExp == "") {
              "Please enter gene or gene list to view bulk expression."
            } else if(input$Gene_bulkExp != "" & gene_found_bulk == 0) {
              "The gene entered was not found in the expression data."
            }
          }
          else if (input$gene_select == "gene_set") {
            gs = getGenes(input$GeneSet_bulkExp)

            not_found = "The following genes were not found in the expression data, please remove from input:"
            check_found = FALSE
            for (i in 1:length(gs)) {
              if (!(toupper(gs[i]) %in% toupper(rownames(bulkExp)))) {
                not_found = paste(not_found, gs[i])
                check_found = TRUE
              }
            }

            format_error = FALSE
            if (any(grep(",", gs))) {
              format_error = TRUE
              check_found = FALSE
            }

            single_error = FALSE
            print(gs)
            if (length(gs) == 1) {
              single_error = TRUE
              check_found = FALSE
            }

            if (check_found) {
              not_found
            } else if (format_error) {
              "Please check formatting. There should be one gene per line, no extra spaces or commas."
            } else if (single_error) {
              "If you would like to view the expression of one gene, please use the 'Single Gene' selection from the 'Select Gene Expression Viewing' dropdown."
            }
          }
        })
      } else if (input$sidebarItemExpanded == "Prediction" & input$dataset == "leukemia") {
        req(input$go_predict)

        file_in <- input$user_in # user input 
        ext <- tools::file_ext(file_in$datapath)
        genes = loaded.data$genes

        req(file_in)
        validate(need(ext == "csv", "Please upload a csv file"))

        user.data = read.csv(file_in$datapath, header = T, row.names = 1)
        colnames(user.data)[which(colnames(user.data) == "NME1.NME2")] = "NME1-NME2"

        check_genes = length(which(colnames(user.data) %in% genes))
        if (check_genes < 10) {
          "Please include at least 10 of the required biomarkers in your uploaded dataset."
        }
      } else {""}
    })

   ## show dataset chosen
   output$dataset_chosen <- renderText({
    if (length(input$dataset) == 0) {
      ""
    } else if (input$dataset == "leukemia") {
      "Analyzing Acute Leukemia Dataset. Refresh the page to choose a different dataset."
    } else if (input$dataset == "wilmsTumor") {
      "Analyzing Wilms Tumor Dataset. Refresh the page to choose a different dataset. The Prediction module is not available for this dataset."
    } else if (input$dataset == "rms") {
      "Analyzing Rhabdomyosarcoma (RMS) dataset. Refresh the page to choose a different dataset. The Prediction and BulkExpression modules are not available for this dataset."
    } else if (input$dataset == "epn") {
      "Analyzing Ependymoma (EPN) dataset. Refresh the page to choose a different dataset. The Prediction and BulkExpression modules are not available for this dataset."
    } else if (input$dataset == "preALL") {
      "Analyzing Pre-B and Pre-T ALL dataset. Refresh the page to choose a different dataset. The Prediction and BulkExpression modules are not available for this dataset."
    }
   })


  }
)