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
              plots for models selected by BIC, ICL, AIC3 and AIC."),
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
                  choices = c("kmeans",
                              "random",
                              "medoids",
                              "clara",
                              "fanny")),
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
                  tabPanel("Cluster Results", verbatimTextOutput('clustering')),
                  tabPanel("Model Selection",
                           fluidRow(
                             splitLayout(cellWidths = c("50%", "50%"), plotOutput('BICvalues'), plotOutput('ICLvalues')),
                             splitLayout(cellWidths = c("50%", "50%"), plotOutput('AIC3values'), plotOutput('AICvalues')),
                           )),
                  tabPanel("Heatmap",
                           fluidRow(
                             splitLayout(cellWidths = c("50%", "50%"), plotOutput("heatmapBIC"), plotOutput('heatmapICL')),
                             splitLayout(cellWidths = c("50%", "50%"), plotOutput("heatmapAIC3"), plotOutput('heatmapAIC')),
                           )),
                  tabPanel("Barplot",
                           fluidRow(
                             splitLayout(cellWidths = c("50%", "50%"), plotOutput("barPlotBIC"), plotOutput('barPlotICL')),
                             splitLayout(cellWidths = c("50%", "50%"), plotOutput("barPlotAIC3"), plotOutput('barPlotAIC'))
                           ))

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
      as.matrix(read.csv(input$file1$datapath,
                         sep = ",",
                         header = TRUE,
                         row.names = 1))
  })


  startclustering <- eventReactive(eventExpr = input$button2, {
    withProgress(message = 'Clustering', value = 0, {
      # Number of times we'll go through the loop

      MPLNClust::mplnVariational(
        dataset = matrixInput(),
        membership = "none",
        gmin = as.numeric(input$ngmin),
        gmax = as.numeric(input$ngmax),
        initMethod = as.character(input$typeinitMethod),
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

  # plot logL
  output$logL <- renderPlot({
    if (! is.null(startclustering))

      if (length(startclustering()$logLikelihood) == 1) { # check if only one value
        if(as.numeric(input$ngmax) == 1) { # check if only one value is because gmax = 1
          plot(c(startclustering()$logLikelihood), type = "p",
               xlab = "G", ylab = "logL",
               main = paste("G vs log-likelihood"))
        } else { # check if only one value is because only one model is tested e.g., gmin = 4, gmax = 4
          plot(c(rep(NA, as.numeric(input$ngmax) - 1), startclustering()$logLikelihood),
               type = "p", xlab = "G", ylab = "logL",
               main = paste("G vs log-likelihood"))
        }
      } else { # ff more than one value
        plot(x = c(as.numeric(input$ngmin):as.numeric(input$ngmax)),
             y = startclustering()$logLikelihood, type = "l",
             lty = 2, xlab = "G", ylab = "logL",
             main = paste("G vs log-likelihood"), xaxt="n")
        axis(1, at = seq(as.numeric(input$ngmin), as.numeric(input$ngmax), by = 1))
      }
  })

  # plot ICL value
  output$ICLvalues <- renderPlot({
    if (! is.null(startclustering))
      if (length(startclustering()$logLikelihood) == 1) { # check if only one value
        if(as.numeric(input$ngmax) == 1) { # check if only one value is because gmax = 1
          plot(c(startclustering()$ICLresults$allICLvalues), type = "p",
               xlab = "G", ylab = "ICL value",
               main = paste("G vs ICL value"))
        } else { # check if only one value is because only one model is tested e.g., gmin = 4, gmax = 4
          plot(c(rep(NA, as.numeric(input$ngmax) - 1), startclustering()$ICLresults$allICLvalues),
               type = "p", xlab = "G", ylab = "ICL value",
               main = paste("G vs ICL value"))
        }
      } else { # ff more than one value
        plot(x = c(as.numeric(input$ngmin):as.numeric(input$ngmax)),
             y = startclustering()$ICLresults$allICLvalues, type = "l",
             lty = 2, xlab = "G", ylab = "ICL value",
             main = paste("G vs ICL value"), xaxt="n")
        axis(1, at = seq(as.numeric(input$ngmin), as.numeric(input$ngmax), by = 1))
      }
  })


  # plot BIC value
  output$BICvalues <- renderPlot({
    if (! is.null(startclustering))
      if (length(startclustering()$logLikelihood) == 1) { # check if only one value
        if(as.numeric(input$ngmax) == 1) { # check if only one value is because gmax = 1
          plot(c(startclustering()$BICresults$allBICvalues), type = "p",
               xlab = "G", ylab = "BIC value",
               main = paste("G vs BIC value"))
        } else { # check if only one value is because only one model is tested e.g., gmin = 4, gmax = 4
          plot(c(rep(NA, as.numeric(input$ngmax) - 1), startclustering()$BICresults$allBICvalues),
               type = "p", xlab = "G", ylab = "BIC value",
               main = paste("G vs BIC value"))
        }
      } else { # ff more than one value
        plot(x = c(as.numeric(input$ngmin):as.numeric(input$ngmax)),
             y = startclustering()$BICresults$allBICvalues, type = "l",
             lty = 2, xlab = "G", ylab = "BIC value",
             main = paste("G vs BIC value"), xaxt="n")
        axis(1, at = seq(as.numeric(input$ngmin), as.numeric(input$ngmax), by = 1))
      }
  })

  # plot AIC value
  output$AICvalues <- renderPlot({
    if (! is.null(startclustering))
      if (length(startclustering()$logLikelihood) == 1) { # check if only one value
        if(as.numeric(input$ngmax) == 1) { # check if only one value is because gmax = 1
          plot(c(startclustering()$AICresults$allAICvalues), type = "p",
               xlab = "G", ylab = "AIC value",
               main = paste("G vs AIC value"))
        } else { # check if only one value is because only one model is tested e.g., gmin = 4, gmax = 4
          plot(c(rep(NA, as.numeric(input$ngmax) - 1), startclustering()$AICresults$allAICvalues),
               type = "p", xlab = "G", ylab = "AIC value",
               main = paste("G vs AIC value"))
        }
      } else { # ff more than one value
        plot(x = c(as.numeric(input$ngmin):as.numeric(input$ngmax)),
             y = startclustering()$AICresults$allAICvalues, type = "l",
             lty = 2, xlab = "G", ylab = "AIC value",
             main = paste("G vs AIC value"), xaxt="n")
        axis(1, at = seq(as.numeric(input$ngmin), as.numeric(input$ngmax), by = 1))
      }
  })

  # plot AIC3 value
  output$AIC3values <- renderPlot({
    if (! is.null(startclustering))
      if (length(startclustering()$logLikelihood) == 1) { # check if only one value
        if(as.numeric(input$ngmax) == 1) { # check if only one value is because gmax = 1
          plot(c(startclustering()$AIC3results$allAIC3values), type = "p",
               xlab = "G", ylab = "AIC3 value",
               main = paste("G vs AIC3 value"))
        } else { # check if only one value is because only one model is tested e.g., gmin = 4, gmax = 4
          plot(c(rep(NA, as.numeric(input$ngmax) - 1), startclustering()$AIC3results$allAIC3values),
               type = "p", xlab = "G", ylab = "AIC3 value",
               main = paste("G vs AIC3 value"))
        }
      } else { # ff more than one value
        plot(x = c(as.numeric(input$ngmin):as.numeric(input$ngmax)),
             y = startclustering()$AIC3results$allAIC3values, type = "l",
             lty = 2, xlab = "G", ylab = "AIC3 value",
             main = paste("G vs AIC3 value"), xaxt="n")
        axis(1, at = seq(as.numeric(input$ngmin), as.numeric(input$ngmax), by = 1))
      }
  })



  # plot heatmap - BIC
  heatmapPlottingBIC <- eventReactive(eventExpr = input$button2, {
    if (!is.null(startclustering))
      mplnVisualize(
        dataset = matrixInput(),
        plots = "heatmaps",
        clusterMembershipVector = as.numeric(startclustering()$BICresults$BICmodelSelectedLabels),
        printPlot = FALSE)
  })

  # plot heatmap - BIC
  output$heatmapBIC <- renderPlot({
    heatmapPlottingBIC()
  })



  # plot bar - BIC
  barPlottingBIC <- eventReactive(eventExpr = input$button2, {
    if (!is.null(startclustering))
      if ((as.numeric(input$ngmax) - as.numeric(input$ngmin) + 1) == 1) {
        mplnVisualize(
          dataset = matrixInput(),
          plots = "bar",
          probabilities = as.matrix(startclustering()$allResults[[1]]$probaPost),
          clusterMembershipVector = as.numeric(startclustering()$BICresults$BICmodelSelectedLabels),
          printPlot = FALSE)
      } else {
        modelSelect <- which(seq(as.numeric(input$ngmin), as.numeric(input$ngmax), 1) == startclustering()$BICresults$BICmodelselected)
        mplnVisualize(
          dataset = matrixInput(),
          plots = "bar",
          probabilities = as.matrix(startclustering()$allResults[[as.numeric(modelSelect)]]$probaPost),
          clusterMembershipVector = as.numeric(startclustering()$BICresults$BICmodelSelectedLabels),
          printPlot = FALSE)
      }
  })

  # plot bar - BIC
  output$barPlotBIC <- renderPlot({
    barPlottingBIC()
  })






  # plot heatmap - ICL
  heatmapPlottingICL <- eventReactive(eventExpr = input$button2, {
    if (!is.null(startclustering))
      mplnVisualize(
        dataset = matrixInput(),
        plots = "heatmaps",
        clusterMembershipVector = as.numeric(startclustering()$ICLresults$ICLmodelSelectedLabels),
        printPlot = FALSE)
  })

  # plot heatmap - ICL
  output$heatmapICL <- renderPlot({
    heatmapPlottingICL()
  })



  # plot bar - ICL
  barPlottingICL <- eventReactive(eventExpr = input$button2, {
    if (!is.null(startclustering))
      if ((as.numeric(input$ngmax) - as.numeric(input$ngmin) + 1) == 1) {
        mplnVisualize(
          dataset = matrixInput(),
          plots = "bar",
          probabilities = as.matrix(startclustering()$allResults[[1]]$probaPost),
          clusterMembershipVector = as.numeric(startclustering()$ICLresults$ICLmodelSelectedLabels),
          printPlot = FALSE)
      } else {
        modelSelect <- which(seq(as.numeric(input$ngmin), as.numeric(input$ngmax), 1) == startclustering()$ICLresults$ICLmodelselected)
        mplnVisualize(
          dataset = matrixInput(),
          plots = "bar",
          probabilities = as.matrix(startclustering()$allResults[[as.numeric(modelSelect)]]$probaPost),
          clusterMembershipVector = as.numeric(startclustering()$ICLresults$ICLmodelSelectedLabels),
          printPlot = FALSE)
      }
  })

  # plot bar - ICL
  output$barPlotICL <- renderPlot({
    barPlottingICL()
  })










  # plot heatmap - AIC
  heatmapPlottingAIC <- eventReactive(eventExpr = input$button2, {
    if (!is.null(startclustering))
      mplnVisualize(
        dataset = matrixInput(),
        plots = "heatmaps",
        clusterMembershipVector = as.numeric(startclustering()$AICresults$AICmodelSelectedLabels),
        printPlot = FALSE)
  })

  # plot heatmap - AIC
  output$heatmapAIC <- renderPlot({
    heatmapPlottingAIC()
  })



  # plot bar - AIC
  barPlottingAIC <- eventReactive(eventExpr = input$button2, {
    if (!is.null(startclustering))
      if ((as.numeric(input$ngmax) - as.numeric(input$ngmin) + 1) == 1) {
        mplnVisualize(
          dataset = matrixInput(),
          plots = "bar",
          probabilities = as.matrix(startclustering()$allResults[[1]]$probaPost),
          clusterMembershipVector = as.numeric(startclustering()$AICresults$AICmodelSelectedLabels),
          printPlot = FALSE)
      } else {
        modelSelect <- which(seq(as.numeric(input$ngmin), as.numeric(input$ngmax), 1) == startclustering()$AICresults$AICmodelselected)
        mplnVisualize(
          dataset = matrixInput(),
          plots = "bar",
          probabilities = as.matrix(startclustering()$allResults[[as.numeric(modelSelect)]]$probaPost),
          clusterMembershipVector = as.numeric(startclustering()$AICresults$AICmodelSelectedLabels),
          printPlot = FALSE)
      }
  })

  # plot bar - AIC
  output$barPlotAIC <- renderPlot({
    barPlottingAIC()
  })









  # plot heatmap - AIC3
  heatmapPlottingAIC3 <- eventReactive(eventExpr = input$button2, {
    if (!is.null(startclustering))
      mplnVisualize(
        dataset = matrixInput(),
        plots = "heatmaps",
        clusterMembershipVector = as.numeric(startclustering()$AIC3results$AIC3modelSelectedLabels),
        printPlot = FALSE)
  })

  # plot heatmap - AIC3
  output$heatmapAIC3 <- renderPlot({
    heatmapPlottingAIC3()
  })



  # plot bar - AIC3
  barPlottingAIC3 <- eventReactive(eventExpr = input$button2, {
    if (!is.null(startclustering))
      if ((as.numeric(input$ngmax) - as.numeric(input$ngmin) + 1) == 1) {
        mplnVisualize(
          dataset = matrixInput(),
          plots = "bar",
          probabilities = as.matrix(startclustering()$allResults[[1]]$probaPost),
          clusterMembershipVector = as.numeric(startclustering()$AIC3results$AIC3modelSelectedLabels),
          printPlot = FALSE)
      } else {
        modelSelect <- which(seq(as.numeric(input$ngmin), as.numeric(input$ngmax), 1) == startclustering()$AIC3results$AIC3modelselected)
        mplnVisualize(
          dataset = matrixInput(),
          plots = "bar",
          probabilities = as.matrix(startclustering()$allResults[[as.numeric(modelSelect)]]$probaPost),
          clusterMembershipVector = as.numeric(startclustering()$AIC3results$AIC3modelSelectedLabels),
          printPlot = FALSE)
      }
  })

  # plot bar - AIC3
  output$barPlotAIC3 <- renderPlot({
    barPlottingAIC3()
  })

}

# Create Shiny app ----
shinyApp(ui, server)
