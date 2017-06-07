
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
# 
# http://www.rstudio.com/shiny/
#

library(shiny)
library(shinyBS)



shinyUI(pageWithSidebar(
  
  # Application title
  headerPanel(globalOptions[["title"]]),
  
  # Sidebar with a slider input for number of bins
  sidebarPanel(
    conditionalPanel(
      'input.tab == "Slice Series" || input.tab == "Single Slice"',
      sliderInput("range",
                  "Threshold",
                  min = 0,
                  max = 15,
                  step=0.1,
                  value = c(2, 15)),
      uiOutput("statsDescription")),
               
    #conditionalPanel(
    #  'input.tab == "Slice Series" || input.tab == "Single Slice"',
    #  sliderInput("high",
    #              "Upper Threshold",
    #              min = 0,
    #              max = 15,
    #              step=0.1,
    #              value = 20)),
    selectInput("statistic", "Statistic to display:", choices=c()),
                #choices=c("Neonatal" = "Neonatal", 
                #          "Sex" = "mouse.gender")),
    conditionalPanel(
      'input.tab == "Slice Series"',
      selectInput("dimension", "Dimension to display:",
                  choices=c("coronal" = 2, 
                            "sagittal" = 1, 
                            "axial" = 3))),
    conditionalPanel(
      'input.tab == "Slice Series"',
      sliderInput("sliceRange",
                  "Slice Range",
                  min=1, max=340, value=c(55, 280))),
    #conditionalPanel(
    #  'input.tab == "Slice Series"',
    #  sliderInput("end",
    #              "Last slice",
    #              min=-340, max=-1, value=-70)),
    conditionalPanel(
      'input.tab == "Slice Series"',
      sliderInput("rows", "Number of rows",
                  min=1, max=10, value=5)),
    conditionalPanel(
      'input.tab == "Slice Series"',
      sliderInput("columns", "Number of columns",
                  min=1, max=10, value=4)),
    conditionalPanel(
      'input.tab == "Single Slice"',
      checkboxInput("updatePlot", "Update plot", TRUE)
    ),
    conditionalPanel(
      'input.tab == "Single Slice" || input.tab == "Volumes"',
      selectInput("xvar", "X axis on plot", choices=globalOptions[["plotChoices"]])
                    #c("Group", "Treatment1","Treatment2"))

    ),
    conditionalPanel(
      'input.tab == "Single Slice" || input.tab == "Volumes"',
      selectInput("colour", "Colour on plot", choices=c("None", globalOptions[["plotChoices"]]))
    ),
    conditionalPanel(
      'input.tab == "Single Slice" || input.tab == "Volumes"',
      selectInput("fill", "fill on plot", choices=c("None", globalOptions[["plotChoices"]]))
    ),

    conditionalPanel(
      'input.tab == "Single Slice" || input.tab == "Volumes"',
      selectInput("graphType", "Graph type", 
                choices=c("boxplot", "point")))
    
  ),
  
  # Show a plot of the generated distribution
  mainPanel(
    tabsetPanel(
      id="tab",
      tabPanel("Slice Series", plotOutput("seriesPlot", height="800px", click="sliceSeriesClick"),
               bsModal("modalSliceSeries", "Export Slice Series", "exportButton", size="small",
                       selectInput("ssFileType", "Image Type", choices=c("PDF", "PNG", "TIFF")),
                       numericInput("ssHeight", "Height (in inches)", value=4, step=0.1),
                       numericInput("ssWidth", "Width (in inches)", value=4, step=0.1),
                       conditionalPanel('input.ssFileType != "PDF"',
                                        numericInput('ssRes', "Resolution (PPI)", value=150, step=1)),
                       textInput("ssTitle", "Plot Title", value=""),
                       downloadButton("ssDownload", "Download")),
               actionButton("exportButton", "Export")
      ),
      tabPanel("Single Slice", 
               style="background-color: black",
               fluidRow(
                 column(7,
                        plotOutput("coronalPlot", height="500px", click="plot_click")),
                 column(5,
                        plotOutput("axialPlot", height="500px", click="click_axial"))),
               fluidRow(
                 column(8,
                        plotOutput("sagittalPlot", click="click_sagittal")),
                 column(4,
                        plotOutput("graphPlot"))
               ),
      fluidRow(column(12, 
                      bsModal("modalSingleSlice", "Export Slices and Plot", "spexportButton", size="small",
                              selectInput("spOutputType", "Output", choices=c("Single Image" = "single",
                                                                              "Multiple Images in ZIP file" = "zip")),
                              selectInput("spFileType", "Image Type", choices=c("PDF", "PNG", "TIFF")),
                              numericInput("spHeight", "Height (in inches)", value=4, step=0.1),
                              numericInput("spWidth", "Width (in inches)", value=8, step=0.1),
                              conditionalPanel('input.spFileType != "PDF"',
                                               numericInput('spRes', "Resolution (PPI)", value=150, step=1)),
                              downloadButton("spDownload", "Download")),
                      actionButton("spexportButton", "Export"))),        
      fluidRow(

                 column(12,
                        verbatimTextOutput("summaryText")))),
      tabPanel("Volumes", 
               fluidRow(
                 DT::dataTableOutput(outputId="volumesTable")),
               fluidRow(
                 plotOutput("volumesPlot")
                 #column(6, 
                  #      plotOutput("volumesPlot2"))
               )
               )
               #),
      
    )
  )
))

