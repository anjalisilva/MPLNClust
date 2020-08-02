library(shiny)

# Define UI for random distribution app ----
ui <- fluidPage(

  # App title ----
  titlePanel(tags$h1(tags$b("MPLNClust:"),"Mixtures of MPLN for Clustering via Variational-EM")),

  # Sidebar layout with input and output definitions ----
  sidebarLayout(

    # Sidebar panel for inputs ----
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
      textInput(inputId = "ngmin",
                label = "Enter gmin:", "1"),
      textInput(inputId = "ngmax",
                label = "Enter gmax:", "2"),
      selectInput(inputId = 'typeinitMethod',
                  label = 'Select initMethod:',
                  choices = c("'kmeans' ",
                              "'random' ",
                              "'medoids' ",
                              "'clara' ",
                              "'fanny' ")),
      textInput(inputId = "nInitIterations",
                label = "Enter nInitIterations:", "1"),
      selectInput(inputId = 'typenormalize',
                  label = 'Select normalize:',
                  choices = c("'Yes' ",
                              "'No' ")),

      # actionButton
      actionButton(inputId = "button2",
                   label = "Start Clustering"),

      # br() element to introduce extra vertical spacing -
      br(),

    ), # End of side pannel


    # Main panel for displaying outputs
    mainPanel(

      # Output: Tabet
      tabsetPanel(type = "tabs",
                  tabPanel("Input", plotOutput("pairsplot")),
                  tabPanel("Cluster Results",
                           fluidRow(column(6, verbatimTextOutput('clustering')),
                                    column(6, plotOutput('pairsplot2')))),
                  tabPanel("Heatmap", plotOutput("heatmap")),
                  tabPanel("Barplot", plotOutput("barPlot")),
                  tabPanel("Lineplot", plotOutput("linePlot"))

      )
    )
  )
)

# Define server logic for random distribution app ----
server <- function(input, output) {

  # Reactive expression to generate the requested distribution ----
  # This is called whenever the inputs change. The output functions
  # defined below then use the value computed from this expression


  # Step I: save input csv as a reactive
  matrixInput <- reactive({
    if (! is.null(input$file1))
      as.matrix(read.table(input$file1$datapath, sep = ",", header = TRUE))
  })


  # Pairsplot
  output$pairsplot2 <- output$pairsplot <- shiny::renderPlot({
    if (! is.null(input$file1))
      pairs(matrixInput())
  })

  startclustering <- shiny::eventReactive(eventExpr = input$button2, {
    MPLNClust::mplnVariational(
      dataset = matrixInput(),
      membership = "none",
      gmin = as.numeric(input$ngmin),
      gmax = as.numeric(input$ngmax),
      initMethod = "kmeans",
      nInitIterations = as.numeric(input$nInitIterations),
      normalize = "Yes")
  })


  # Step II: clustering
  output$clustering <- shiny::renderText({
    if (!is.null(startclustering))

      aa <- paste("BIC model selected is:", startclustering()[[12]][[2]], "\n")

    bb <- paste("ICL model selected is:", startclustering()[[11]][[2]], "\n")

    cc <- paste("AIC model selected is:", startclustering()[[13]][[2]], "\n")

    dd <- paste("AIC3 model selected is:", startclustering()[[14]][[2]], "\n")
    paste(aa, bb, cc, dd, sep = "\n")
  })

  # Step III: visualize

  startplotting <- shiny::eventReactive(eventExpr = input$button2, {
    if (!is.null(startclustering))
      mplnVisualize(
        dataset = matrixInput(),
        plots = "heatmaps",
        probabilities = as.matrix(startclustering()$allResults$`G=2`$probaPost),
        clusterMembershipVector = as.numeric(startclustering()$allResults$`G=2`$clusterlabels),
        printPlot = FALSE)
  })

  # plot heatmap
  output$heatmap <- renderPlot({
    startplotting()
  })


  startplotting2 <- shiny::eventReactive(eventExpr = input$button2, {
    if (!is.null(startclustering))
      mplnVisualize(
        dataset = matrixInput(),
        plots = "bar",
        probabilities = as.matrix(startclustering()$allResults$`G=2`$probaPost),
        clusterMembershipVector = as.numeric(startclustering()$allResults$`G=2`$clusterlabels),
        printPlot = FALSE)
  })

  # plot bar
  output$barPlot <- renderPlot({
    startplotting2()
  })



  startplotting3 <- shiny::eventReactive(eventExpr = input$button2, {
    if (! is.null(startclustering))
      mplnVisualize(
        dataset = matrixInput(),
        plots = "lines",
        LinePlotColours = "multicolour",
        probabilities = as.matrix(startclustering()$allResults$`G=2`$probaPost),
        clusterMembershipVector = as.numeric(startclustering()$allResults$`G=2`$clusterlabels),
        printPlot = FALSE)
  })

  # plot bar
  output$linePlot <- renderPlot({
    startplotting3()[1]
  })

}

# Create Shiny app ----
shinyApp(ui, server)
