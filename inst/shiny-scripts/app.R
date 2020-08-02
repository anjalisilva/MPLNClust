library(shiny)

# Define UI for random distribution app ----
ui <- fluidPage(

  # App title ----
  titlePanel(tags$h1(tags$b("MPLNClust:"),"Mixtures of MPLN for Clustering via Variational-EM")),

  # Sidebar layout with input and output definitions ----
  sidebarLayout(

    # Sidebar panel for inputs ----
    sidebarPanel(

      tags$p("Provided a dataset of counts, perform clustering and show
              plots for model selected by BIC."),
      tags$p("Import a dataset of counts, such that rows should
                         correspond to observations and columns correspond to
                         variables in .csv format. Dataset must contain row
                         names (in first column) and column names."),
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
                  tabPanel("Input PairsPlot", plotOutput("pairsplot")),
                  tabPanel("Input Summary", verbatimTextOutput("textOut")),
                  tabPanel("Cluster Results",
                           fluidRow(column(6, verbatimTextOutput('clustering')),
                                    column(6, plotOutput('logL')))),
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


  startclustering <- eventReactive(eventExpr = input$button2, {
    withProgress(message = 'Clustering', value = 0, {
      # Number of times we'll go through the loop

      MPLNClust::mplnVariational(
        dataset = matrixInput()[ , 2:ncol(matrixInput())],
        membership = "none",
        gmin = as.numeric(input$ngmin),
        gmax = as.numeric(input$ngmax),
        initMethod = "kmeans",
        nInitIterations = as.numeric(input$nInitIterations),
        normalize = "Yes")

    })
  })

  # Textoutput
  output$textOut <- renderPrint({
    if (! is.null(startclustering))
      summary(startclustering()$dataset)
  })

  # Pairsplot
  output$pairsplot <- renderPlot({
    if (! is.null(startclustering))
      pairs(startclustering()$dataset)
  })


  # Step II: clustering
  output$clustering <- renderText({
    if (! is.null(startclustering))

      aa <- paste("BIC model selected is:", startclustering()$BICresults$BICmodelselected, "\n")

    bb <- paste("ICL model selected is:", startclustering()$ICLresults$ICLmodelselected, "\n")

    cc <- paste("AIC model selected is:", startclustering()$AICresults$AICmodelselected, "\n")

    dd <- paste("AIC3 model selected is:", startclustering()$AIC3results$AIC3modelselected, "\n")
    paste(aa, bb, cc, dd, sep = "\n")
  })

  # Step III: visualize

  output$logL <- renderPlot({
    if (! is.null(startclustering))

      plot(startclustering()$logLikelihood, type = "l",
           lty = 2, xlab = "G", ylab = "logL",
           main = paste("G vs log-likelihood"))
  })


  startplotting <- eventReactive(eventExpr = input$button2, {
    if (!is.null(startclustering))
      mplnVisualize(
        dataset = matrixInput()[ , 2:ncol(matrixInput())],
        plots = "heatmaps",
        probabilities = as.matrix(startclustering()$allResults$`G=2`$probaPost),
        clusterMembershipVector = as.numeric(startclustering()$allResults$`G=2`$clusterlabels),
        printPlot = FALSE)
  })

  # plot heatmap
  output$heatmap <- renderPlot({
    startplotting()
  })


  startplotting2 <- eventReactive(eventExpr = input$button2, {
    if (!is.null(startclustering))
      mplnVisualize(
        dataset = matrixInput()[ , 2:ncol(matrixInput())],
        plots = "bar",
        probabilities = as.matrix(startclustering()$allResults$`G=2`$probaPost),
        clusterMembershipVector = as.numeric(startclustering()$allResults$`G=2`$clusterlabels),
        printPlot = FALSE)
  })

  # plot bar
  output$barPlot <- renderPlot({
    startplotting2()
  })



  startplotting3 <- eventReactive(eventExpr = input$button2, {
    if (! is.null(startclustering))
      mplnVisualize(
        dataset = matrixInput()[ , 2:ncol(matrixInput())],
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
