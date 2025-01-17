library(shiny)
library(bslib)
library(shinybusy)
library(plotly)

page_navbar(
  id = "navbar",
  title = readRDS("title.rds"),
  sidebar = sidebar(
    selectInput("model", "Model:", choices = NULL),
    selectInput("contrast", "Contrast:", choices = NULL),
    selectInput("geneset", "Geneset:", choices = NULL),
    open = F
  ),
  nav_panel("NB", htmlOutput("notebook")),
  nav_panel("QC", htmlOutput("qc")),
  nav_spacer(),
  nav_panel("MDS", selectInput("mds_var", "MDS Variable:", choices = NULL),
            plotlyOutput("mds")),
  nav_panel("Volcano Plot", textOutput("model_spec"),
            plotOutput("volcano_plot", height = "70vh")),
  nav_panel("GSEA Plot",
    textOutput("model_spec"),
    fixedRow(
      column(2, textInput("gsea_qfilt", "Q-value Filter:", value = "0.001"))),
    plotOutput("gsea_plot", height = "70vh")),
  nav_panel("DGE Tab", textOutput("model_spec"),
            DT::dataTableOutput("dg_tab")),
  nav_panel("GSEA Tab", textOutput("model_spec"),
            DT::dataTableOutput("gsea_tab")),
  nav_panel("Contrasts", textOutput("model_spec"), plotOutput("contrast_plot")),
  nav_panel("Compare",
    fixedRow(
      column(4, checkboxInput("contrast_overlap",
                              "Contrasts overlap", value = TRUE)),
      column(4, selectInput("model2", "Comparison model:", choices = NULL))),
    fixedRow(
      column(4, selectInput("contrast2", "Comparison contrast:",
                            choices = NULL)),
      column(4, selectInput("comp_type", "Comparison type:", choices = c(
        "Differential Expression", "Geneset Enrichment"
      ))),
      column(2, textInput("comp_pfilt", "P-value Filter:", value = "1"))),
    plotOutput("cor_plot"))
)
