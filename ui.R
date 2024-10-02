library(shiny)
library(shinyFiles)
library(bslib)
library(vcfR)
library(dplyr)
library(tibble)
library(tidyr)
library(pheatmap)
library(RColorBrewer)

ui <- fluidPage(
  theme = bslib::bs_theme(bootswatch = "sandstone"),#bg = "#FFFFFF",
                          #fg = "#000000",
                          #primary = "#0199F8",
                          #secondary = "#FF374B"),
  headerPanel(title = "P. falciparum heatmaps", windowTitle = "HeatMap"),
  #title = div(img(src="Cover.png", style="max-width: 100%; height: auto;"),window_title="ALeRT - SMART"),
  sidebarLayout(
    sidebarPanel(
       shinyDirButton("inputdata", "VCF directory", "VCF Directory"),
        uiOutput("file"),
        selectInput(inputId = "typemol",
                   label = "NT or AA",
                   choices = c("NT","AA")),
        sliderInput("af","Allele Frequency Average",value = 0.10, min = 0, max=1.00,step=0.01),
        style = "height: 100%",
        verbatimTextOutput("comando", placeholder = TRUE)
    ),
    mainPanel(
      #fluidRow(div(img(src="Cover.png", style="max-width: 100%; height: auto;"))),
      fluidRow(
        #wellPanel(
          #div(
            #tableOutput("finalR"),
            plotOutput("plot", height = "600px")#,
            #style = "overflow-x: auto; width: 100%; high: 100%"
          #)
          #downloadButton("downloadres", label = "Download Results", class = "butt"), uiOutput("id"),
          #style = "width: 97.5%; margin-bottom: 1%;"
        #)
      )
    )
  )  
)
