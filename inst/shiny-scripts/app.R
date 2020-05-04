library(shiny)

ui <- fluidPage(
  tags$h1(tags$b("MPLNClust:"),"Mixtures of Multivariate Poisson-Log Normal Model for Clustering Count Data"),
  tags$hr(),
  tabsetPanel(

    # Use of the uiOutput function to create drop-down menus that update as new
    # datasets are imported or created is based on code from the following RStudio
    # Shiny tutorial: https://shiny.rstudio.com/articles/dynamic-ui.html

    # A tab that permits the user to import a dataset
    tabPanel("Step I: Import Dataset",
             sidebarLayout(
               sidebarPanel(
                 tags$p("Import a dataset of counts, such that in the dataset rows
                         should correspond to observations and columns correspond to
                         variables in .csv format. File can be previewed after
                         importing."),
                 tags$p("An example maybe created using mplnDataGenerator() function
                        of MPLNClust R package and saved into .csv format using
                        write.csv() function."),
                 fileInput("file", "Select a dataset to import:"),
                 tags$p("Enter or select values required for clustering. Default
                        values are shown."),
                 textInput("Numgmin", "Enter gmin:", "1"),
                 textInput("Numgmax", "Enter gmax:", "2"),
                 textInput("NumnChains", "Enter nChains:", "3"),
                 selectInput('TypeinitMethod',
                             'Select initMethod:',
                             c("'kmeans' ",
                               "'random' ",
                               "'medoids' ",
                               "'clara' ",
                               "'fanny' ")),
                 textInput("NumnChains", "Enter nInitIterations:", "3"),
                 selectInput('Typenormalize',
                             'Select normalize:',
                             c("'Yes' ",
                               "'No' ")),
                 textInput("numNodes", "Enter numNodes:", "NA"),
               ),
               mainPanel(tableOutput("previewTable"))
             )

    ),

    # This tab allows the user to map imported peaks onto imported features
    # using mapPeaks
    tabPanel("Step II: Clustering",
             fluidRow(
               column(4, wellPanel(
                 tags$p("Datasets imported in the previous step will be clustered.)"),
                 tags$p("Results can be displayed in the next tab."),
                 uiOutput("selectPeaks"),
                 uiOutput("selectFeatures"),
                 actionButton(inputId = "mappingButton", label =
                                "Clustering  features"),
                 textOutput("mappingMessage"),
                 tags$p(),
                 uiOutput("selectResult"),
                 actionButton(inputId = "mappingDisplayButton",
                              label = "Display clustering progress"),
                 downloadButton("exportButton",
                                "Export clustering results as CSV file"),
                 tags$p(),
                 tags$p("(For file exporting, the filename should end with '.csv')")
               ))),
             fluidRow(
               column(12, tableOutput("mappingTable"))
             )
    ),

    # This tab allows the user to visualize the results of mapping peaks onto
    # closest features
    tabPanel("Step III: Visualize",
             sidebarLayout(
               sidebarPanel(
                 tags$p("There are  3 types of plots. Plots generated
                 plots can be saved as PNGs by clicking the 'Save' button."),
                 uiOutput("selectResultsForPlotting"),
                 uiOutput("selectPeaksForPlotting"),
                 selectInput('selectPlotType',
                             'Select the type of plot you want to make:',
                             c("1. Heatmap",
                               "2. Line plots",
                               "3. Bar plots")),

                 tags$p("The following parameters can be used to customize the
                 appearance of the generated plot. Selecting a new plot type
                 from the above menu will update the boxes for main title, x-axis
                 title, and y-axis title to the default titles for that plot
                 type."),
                 textInput("plotColor", "Enter plot color:", "Red"),
                 textInput("mainTitle", "Enter main plot title:"),
                 textInput("xTitle", "Enter plot x-axis title:"),
                 textInput("yTitle", "Enter plot y-axis title:"),
                 selectInput("backgroundStyle", "Select style of plot background:",
                             c("blackAndWhite", "grey", "minimal")),

                 tags$p("The following three parameters "),

                 numericInput("upstreamBins", "Enter the number of bins (plotted points)
                       to be used upstream of each feature start position:", 250),
                 numericInput("downstreamBins", "Enter the number of bins (plotted points)
                       to be used downstream of each feature start position:", 250),
                 numericInput("basesPerBin", "Enter the width of each bin, in base pairs:", 10),
                 actionButton(inputId = "plottingButton", label = "Generate plot"),
                 downloadButton("savePlot", "Save current plot as PNG"),
                 tags$p(),
                 tags$p("(For file exporting, the filename should end with '.png')")
               ),
               mainPanel(plotOutput("mainPlot", height = "750px"))
             )
    )
  )
)

server <- function(input, output, clientData, session) {

  output$mytable = DT::renderDataTable({
    output$file
  })
}

shinyApp(ui = ui, server = server)
