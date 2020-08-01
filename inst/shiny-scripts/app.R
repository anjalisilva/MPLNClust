library(shiny)
library(vroom)

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
                 tags$p("Import a dataset of counts, such that rows should
                         correspond to observations and columns correspond to
                         variables in .csv format."),
                 tags$p("A simulation dataset maybe created using mplnDataGenerator()
                        function of MPLNClust R package and saved into .csv format using
                        write.csv() function."),

                 # input
                 fileInput(inputId = "file1",
                           label = "Select a dataset to import:",
                           accept = c(".csv")),
                 tags$p("Enter or select values required for clustering. Default
                        values are shown."),
                 textInput(inputId = 'membership',
                           label = 'Enter memebrship vector seperated by commas,
                           e.g., 1,2,3...; If not available, leave as none:', "none"),
                 textInput(inputId = "ngmin",
                           label = "Enter gmin:", "1"),
                 textInput(inputId = "ngmax",
                           label = "Enter gmax:", "2"),
                 textInput(inputId = "nChains",
                           label = "Enter nChains:", "3"),
                 textInput(inputId = "nIterations",
                           label = "Enter nIterations:", "500"),
                 selectInput(inputId = 'typeinitMethod',
                             label = 'Select initMethod:',
                             choices = c("'kmeans' ",
                                         "'random' ",
                                         "'medoids' ",
                                         "'clara' ",
                                         "'fanny' ")),
                 textInput(inputId = "nInitIterations",
                           label = "Enter nInitIterations:", "3"),
                 selectInput(inputId = 'typenormalize',
                             label = 'Select normalize:',
                             choices = c("'Yes' ",
                                         "'No' ")),
                 textInput(inputId = "nNodes",
                           label = "Enter numNodes:", "NA"),



                 # output
                 plotOutput(outputId = "pairsplot"),

                 # actionButton
                 actionButton(inputId = "button2",
                              label = "Start Clustering"),

                 textOutput(outputId = "clustering")




               ),
               mainPanel(tableOutput(outputId = "previewTable"))
             )

    ),

    # This tab permits to track clustering
    # using mapPeaks
    tabPanel("Step II: Clustering",

    ),

    # This tab allows the user to visualize the results of mapping peaks onto
    # closest features
    tabPanel("Step III: Visualize",
             sidebarLayout(
               sidebarPanel(
                 tags$p("There are  3 types of plots. Plots generated
                 plots can be saved as PNGs by clicking the 'Save Plot'."),
                 selectInput(inputId = 'selectPlotType',
                             label = 'Select the type of plot:',
                             choices = c("1. Heatmap",
                                         "2. Line plots",
                                         "3. Bar plots")),

                 tags$p("The following parameters can be used to customize the
                 appearance of the generated plot. Selecting a new plot type
                 from the above menu will update the boxes for main title, x-axis
                 title, and y-axis title to the default titles for that plot
                 type."),
                 textInput(inputId = "mainTitle",
                           label = "Enter main plot title:"),
                 textInput(inputId = "xTitle",
                           label = "Enter plot x-axis title:"),
                 textInput(inputId = "yTitle",
                           label = "Enter plot y-axis title:"),
                 selectInput(inputId = "backgroundStyle",
                             label = "Select style of plot background:",
                             choices = c("blackAndWhite", "grey", "minimal")),

                 tags$p("The following three parameters "),


                 actionButton(inputId = "plottingButton",
                              label = "Generate Plot"),
                 downloadButton(outputId = "savePlot",
                                label = "Save Plot"),
                 tags$p()
               ),
               mainPanel(plotOutput(outputId = "threePlots"),
                         plotOutput(outputId = "threePlotsecond"))
             )
    )
  )
)

server <- function(input, output, clientData, session) {

  # code
  # https://community.rstudio.com/t/shiny-app-fileinput-getting-variables-from-within-that-file/2338/5

  # save input csv as a reactive
  matrixInput <- reactive({
    if (! is.null(input$file1))
      as.matrix(read.table(input$file1$datapath, sep = ",", header = TRUE))
  })


  # Pairsplot
  output$pairsplot <- renderPlot({
    if (! is.null(input$file1))
      pairs(matrixInput())
  })

  # Running mplnParallel
  # output$clustering <- renderText({
  #   if (! is.null(input$file1))
  #     anotheroutput <- MPLNClust::mplnParallel(
  #       dataset = matrixInput(),
  #       membership = "none",
  #       gmin = as.numeric(input$ngmin),
  #       gmax = as.numeric(input$ngmax),
  #       nChains = as.numeric(input$nChains),
  #       nIterations = as.numeric(input$nIterations),
  #       initMethod = "kmeans",
  #       nInitIterations = as.numeric(input$nInitIterations),
  #       normalize = "Yes",
  #       numNodes = NA)$normalization_factors
  #
  #   if (!is.null(input$file1))
  #     cat("\n Normalization factors", anotheroutput)
  # })


  startclustering <- shiny::eventReactive(eventExpr = input$button2, {
    if (class(input$membership) == "character") {
      MPLNClust::mplnVariational(
                dataset = matrixInput(),
                membership = "none",
                gmin = as.numeric(input$ngmin),
                gmax = as.numeric(input$ngmax),
                initMethod = "kmeans",
                nInitIterations = as.numeric(input$nInitIterations),
                normalize = "Yes")
    } else {
      MPLNClust::mplnVariational(
                dataset = matrixInput(),
                membership = as.numeric(input$membership),
                gmin = as.numeric(input$ngmin),
                gmax = as.numeric(input$ngmax),
                initMethod = "kmeans",
                nInitIterations = as.numeric(input$nInitIterations),
                normalize = "Yes")
    }})

  # Step II - printing  mplnParallel output
  output$clustering <- shiny::renderText({
    if (!is.null(input$file1))

    aa <- paste("BIC model selected is:", startclustering()[[12]][[2]], ";")
    bb <- paste("ICL model selected is:", startclustering()[[11]][[2]], ";")
    cc <- paste("AIC model selected is:", startclustering()[[13]][[2]], ";")
    dd <- paste("AIC3 model selected is:", startclustering()[[14]][[2]])
    paste(aa, bb, cc, dd, sep = "\n")
  })


  startplotting <- shiny::eventReactive(eventExpr = input$plottingButton, {
    MPLNClust::mplnVisualize(
      dataset = matrixInput(),
      plots = "heatmaps",
      probabilities = as.matrix(startclustering()[[7]][[1]][[4]][[9]]),
      clusterMembershipVector = as.numeric(startclustering()[[12]][[3]]),
      fileName = paste0("Plot_", date()),
      LinePlotColours = "black",
      format = "pdf")
  })

  # MPLNClust plots
  output$threePlots <- renderPlot({
    startplotting()[1]
  })

  output$threePlotsecond <- renderPlot({
    startplotting()[2]
  })

}

shinyApp(ui = ui, server = server)
