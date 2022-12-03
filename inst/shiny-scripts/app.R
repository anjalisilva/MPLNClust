library(shiny)
library(shinyalert)

# Define UI for random distribution app ----
ui <- fluidPage(

  # App title ----
  titlePanel(tags$h1(tags$b("MPLNClust:"),"Mixtures of MPLN for Clustering via Variational-EM")),

  # Sidebar layout with input and output definitions ----
  sidebarLayout(

    # Sidebar panel for inputs ----
    sidebarPanel(

      tags$p("Description: This is a Shiny App that is part of the MPLNClust R package
             (Silva et al., 2019). Most of the functions available via the package are made
             available with Shiny App. The MPLNClust is an R package for performing
             clustering using mixtures of multivariate Poisson-log normal
             (MPLN) distribution provided a count dataset. The observations of the
             dataset will be clustered into subgroups. The app permits to calculate
             Bayesian information criterion (BIC), Integrated Complete Likelihood (ICL)
             criterion, and Akaike Information Criterion (AIC) values, given log-likelihood,
             number of clusters, dimension of dataset, number of observations, and the
             probability. Provided the original dataset of counts, the dataset could
             be visualized."),

      # br() element to introduce extra vertical spacing ----
      br(),
      br(),

      # input
      tags$b("Instructions: Below, upload a count dataset, and enter or select
              values to perform cluster analysis. Default values are
              shown. Then press 'Run'. Navigate through the different tabs
              to the right to explore the results."),

      # br() element to introduce extra vertical spacing ----
      br(),
      br(),
      br(),
      # input
      shinyalert::useShinyalert(),  # Set up shinyalert
      uiOutput("tab2"),
      actionButton(inputId = "data1",
                   label = "Dataset 1 Details"),
      uiOutput("tab1"),
      actionButton(inputId = "data2",
                   label = "Dataset 2 Details"),
      fileInput(inputId = "file1",
                label = "Dataset: Select a count dataset to analyze. File should be
                in comma-separated value (.csv) format with rows corresponding
                to observations (e.g., genes) and columns to variables (e.g.,
                samples). You may download an example dataset above and explore first.",
                accept = c(".csv")),
      textInput(inputId = "ngmin",
                label = "ngmin: Enter the minimum value of components or clusters
                for clustering. This should be a positive integer.", "1"),
      textInput(inputId = "ngmax",
                label = "ngmax: Enter the maximum value of components or clusters.
                for clustering. This should be a positive integer, bigger
                than ngmin and less than or equal to number of total observations
                in the dataset.", "2"),
      selectInput(inputId = 'typeinitMethod',
                  label = 'initMethod: Select the initialization method.',
                  choices = c("kmeans",
                              "random",
                              "medoids",
                              "clara",
                              "fanny")),
      textInput(inputId = "nInitIterations",
                label = "nInitIterations: Enter the number of initial iterations.
                This should be a positive integer.", "1"),
      selectInput(inputId = 'typenormalize',
                  label = 'Normalization: Select whether to perform normalization
                  or not. Currently, normalization is performed using
                  calculting normalization factors via TMM method of edgeR
                   package (Robinson, et al., 2010). The option Yes is recommended
                   for raw RNA sequencing count data.',
                  choices = c("'Yes' ",
                              "'No' ")),

      # br() element to introduce extra vertical spacing ----
      br(),

      # actionButton
      actionButton(inputId = "button2",
                   label = "Run"),

      # br() element to introduce extra vertical spacing -
      br(),

    ), # End of side pannel


    # Main panel for displaying outputs
    mainPanel(

      # Output: Tabet
      tabsetPanel(type = "tabs",
                  tabPanel("Pairs Plot of Dataset",
                           h3("Instructions: Enter values and click 'Run' at the bottom left side."),
                           h3("Pairs Plot of Log-transformed Count Dataset:"),
                           br(),
                           plotOutput("pairsplot")),
                  tabPanel("Input Summary",
                           h3("Instructions: Enter values and click 'Run' at the bottom left side."),
                           h3("Summary of Count Dataset:"),
                           br(),
                           verbatimTextOutput("textOut")),
                  tabPanel("Cluster Results",
                           h3("Instructions: Enter values and click 'Run' at the bottom left side."),
                           h3("Summary of Clustering Results:"),
                           br(),
                           verbatimTextOutput('clustering')),
                  tabPanel("Information Criteria Plot",
                           h3("Instructions: Enter values and click 'Run' at the bottom left side."),
                           h3("Model Selection Results:"),
                           br(),
                           fluidRow(
                             splitLayout(cellWidths = c("50%", "50%"), plotOutput('BICvalues'), plotOutput('ICLvalues')),
                             splitLayout(cellWidths = c("50%", "50%"), plotOutput('AIC3values'), plotOutput('AICvalues')),
                           )),
                  tabPanel("Heatmap",
                           h3("Instructions: Enter values and click 'Run' at the bottom left side."),
                           h3("Heatmap of Input Dataset with Clusters:"),
                           h5("Note, the plots are in the order of models selected by: BIC (top, left), ICL (top, right) and AIC (bottom, left), AIC3 (bottom, right).
                              The cluster membership is indicated by the color legend to the left."),
                           br(),
                           fluidRow(
                             splitLayout(cellWidths = c("50%", "50%"), plotOutput("heatmapBIC"), plotOutput('heatmapICL')),
                             splitLayout(cellWidths = c("50%", "50%"), plotOutput("heatmapAIC3"), plotOutput('heatmapAIC')),
                           )),
                  tabPanel("Alluvial Plot",
                           h3("Instructions: Enter values and click 'Run' at the bottom left side."),
                           h3("Alluvial Plot of Input Dataset:"),
                           h5("Note, the x-axis values are in the order of BIC, ICL, AIC, AIC3.
                              Colors are assigned based on cluster membership of model selected via BIC."),
                           br(),
                           fluidRow(
                             splitLayout(cellWidths = c("100%"), plotOutput("alluvialPlot")),
                             h5("Note, the x-axis values are in the order of ICL, BIC, AIC, AIC3.
                              Colors are assigned based on cluster membership of model selected via ICL."),
                             splitLayout(cellWidths = c("100%"), plotOutput("alluvialPlot2")),
                             h5("Note, the x-axis values are in the order of AIC3, ICL, BIC, AIC
                              Colors are assigned based on cluster membership of model selected via AIC3."),
                             splitLayout(cellWidths = c("100%"), plotOutput("alluvialPlot3")),
                             h5("Note, the x-axis values are in the order of AIC, AIC3, ICL, BIC
                              Colors are assigned based on cluster membership of model selected via AIC."),
                             splitLayout(cellWidths = c("100%"), plotOutput("alluvialPlot4")),
                           )),
                  tabPanel("Barplot",
                           h3("Instructions: Enter values and click 'Run' at the bottom left side."),
                           h3("Barplot of Posterior Probabilities with Clusters:"),
                           h5("Note, the plots are in the order of models selected by: BIC (top, left), ICL (top, right) and AIC (bottom, left), AIC3 (bottom, right)."),
                           br(),
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
      mplnVisualizeHeatmap(dataset = matrixInput(),
                           clusterMembershipVector =
                           as.numeric(startclustering()$BICresults$BICmodelSelectedLabels),
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
        mplnVisualizeBar(
          dataset = matrixInput(),
          probabilities = as.matrix(startclustering()$allResults[[1]]$probaPost),
          clusterMembershipVector = as.numeric(startclustering()$BICresults$BICmodelSelectedLabels),
          printPlot = FALSE)
      } else {
        modelSelect <- which(seq(as.numeric(input$ngmin), as.numeric(input$ngmax), 1) == startclustering()$BICresults$BICmodelselected)
        mplnVisualizeBar(
          dataset = matrixInput(),
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
     mplnVisualizeHeatmap(dataset = matrixInput(),
                         clusterMembershipVector =
                         as.numeric(startclustering()$ICLresults$ICLmodelSelectedLabels),
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
        mplnVisualizeBar(
          dataset = matrixInput(),
          probabilities = as.matrix(startclustering()$allResults[[1]]$probaPost),
          clusterMembershipVector = as.numeric(startclustering()$ICLresults$ICLmodelSelectedLabels),
          printPlot = FALSE)
      } else {
        modelSelect <- which(seq(as.numeric(input$ngmin), as.numeric(input$ngmax), 1) == startclustering()$ICLresults$ICLmodelselected)
        mplnVisualizeBar(
          dataset = matrixInput(),
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
      mplnVisualizeHeatmap(dataset = matrixInput(),
                           clusterMembershipVector =
                             as.numeric(startclustering()$AICresults$AICmodelSelectedLabels),
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
        mplnVisualizeBar(
          dataset = matrixInput(),
          probabilities = as.matrix(startclustering()$allResults[[1]]$probaPost),
          clusterMembershipVector = as.numeric(startclustering()$AICresults$AICmodelSelectedLabels),
          printPlot = FALSE)
      } else {
        modelSelect <- which(seq(as.numeric(input$ngmin), as.numeric(input$ngmax), 1) == startclustering()$AICresults$AICmodelselected)
        mplnVisualizeBar(
          dataset = matrixInput(),
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
    mplnVisualizeHeatmap(dataset = matrixInput(),
                         clusterMembershipVector =
                         as.numeric(startclustering()$AIC3results$AIC3modelSelectedLabels),
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
        mplnVisualizeBar(
          dataset = matrixInput(),
          probabilities = as.matrix(startclustering()$allResults[[1]]$probaPost),
          clusterMembershipVector = as.numeric(startclustering()$AIC3results$AIC3modelSelectedLabels),
          printPlot = FALSE)
      } else {
        modelSelect <- which(seq(as.numeric(input$ngmin), as.numeric(input$ngmax), 1) == startclustering()$AIC3results$AIC3modelselected)
        mplnVisualizeBar(
          dataset = matrixInput(),
          probabilities = as.matrix(startclustering()$allResults[[as.numeric(modelSelect)]]$probaPost),
          clusterMembershipVector = as.numeric(startclustering()$AIC3results$AIC3modelSelectedLabels),
          printPlot = FALSE)
      }
  })

  # plot bar - AIC3
  output$barPlotAIC3 <- renderPlot({
    barPlottingAIC3()
  })



  # Alluvial plot
  alluvialPlotting <- eventReactive(eventExpr = input$button2, {
    if (!is.null(startclustering))
      mplnVisualizeAlluvial(nObservations = nrow(matrixInput()),
                            firstGrouping =
                            as.numeric(startclustering()$BICresults$BICmodelSelectedLabels),
                            secondGrouping =
                            as.numeric(startclustering()$ICLresults$ICLmodelSelectedLabels),
                            thirdGrouping =
                            as.numeric(startclustering()$AICresults$AICmodelSelectedLabels),
                            fourthGrouping =
                            as.numeric(startclustering()$AIC3results$AIC3modelSelectedLabels),
                            fileName = 'alluvial',
                            printPlot = FALSE)
  })

  alluvialPlotting2 <- eventReactive(eventExpr = input$button2, {
    if (!is.null(startclustering))
      mplnVisualizeAlluvial(nObservations = nrow(matrixInput()),
                            firstGrouping =
                              as.numeric(startclustering()$ICLresults$ICLmodelSelectedLabels),
                            secondGrouping =
                              as.numeric(startclustering()$BICresults$BICmodelSelectedLabels),
                            thirdGrouping =
                              as.numeric(startclustering()$AICresults$AICmodelSelectedLabels),
                            fourthGrouping =
                              as.numeric(startclustering()$AIC3results$AIC3modelSelectedLabels),
                            fileName = 'alluvial',
                            printPlot = FALSE)
  })

  alluvialPlotting3 <- eventReactive(eventExpr = input$button2, {
    if (!is.null(startclustering))
      mplnVisualizeAlluvial(nObservations = nrow(matrixInput()),
                            firstGrouping =
                              as.numeric(startclustering()$AIC3results$AIC3modelSelectedLabels),
                            secondGrouping =
                              as.numeric(startclustering()$ICLresults$ICLmodelSelectedLabels),
                            thirdGrouping =
                              as.numeric(startclustering()$BICresults$BICmodelSelectedLabels),
                            fourthGrouping =
                              as.numeric(startclustering()$AICresults$AICmodelSelectedLabels),
                            fileName = 'alluvial',
                            printPlot = FALSE)
  })


  alluvialPlotting4 <- eventReactive(eventExpr = input$button2, {
    if (!is.null(startclustering))
      mplnVisualizeAlluvial(nObservations = nrow(matrixInput()),
                            firstGrouping =
                              as.numeric(startclustering()$AICresults$AICmodelSelectedLabels),
                            secondGrouping =
                              as.numeric(startclustering()$AIC3results$AIC3modelSelectedLabels),
                            thirdGrouping =
                              as.numeric(startclustering()$ICLresults$ICLmodelSelectedLabels),
                            fourthGrouping =
                              as.numeric(startclustering()$BICresults$BICmodelSelectedLabels),
                            fileName = 'alluvial',
                            printPlot = FALSE)
  })


  # Alluvial Plot
  output$alluvialPlot <- renderPlot({
    alluvialPlotting()
  })

  output$alluvialPlot2 <- renderPlot({
    alluvialPlotting2()
  })

  output$alluvialPlot3 <- renderPlot({
    alluvialPlotting3()
  })

  output$alluvialPlot4 <- renderPlot({
    alluvialPlotting4()
  })




  # URLs for downloading data
  url1 <- a("Example Dataset 2", href="https://raw.githubusercontent.com/anjalisilva/TestingPackage/master/inst/extdata/GeneCountsData2.csv")
  output$tab1 <- renderUI({
    tagList("Download:", url1)
  })

  observeEvent(input$data2, {
    # Show a modal when the button is pressed
    shinyalert(title = "Example Dataset 2",
               text = "An RNAseq experiment conductd using bean plants from 2016 in Canada. This dataset has n = 30 genes along rows and d = 3 conditions or samples along columns. Data was generated at the University of Guelph, Canada in 2016. To save the file (from Chrome), click on link, then right click, select 'Save As...' and then save as a .csv file.
               Citation: Silva, A. (2020) TestingPackage: An Example R Package For BCB410H. Unpublished. URL https://github.com/anjalisilva/TestingPackage",
               type = "info")
  })

  url2 <- a("Example Dataset 1", href="https://drive.google.com/file/d/1jMBTPpsBwaigjR3mO49AMYDxzjVnNiAv/view?usp=sharing")
  output$tab2 <- renderUI({
    tagList("Download:", url2)
  })

  observeEvent(input$data1, {
    # Show a modal when the button is pressed
    shinyalert(title = "Example Dataset 1",
               text = "This is a simulated dataset generated from mixtures of multivariate Poisson log-normal
               distributions with G = 2 components. It has a size of n = 1000 observations along rows and d = 6
               samples along columns. Data was generated January, 2022. To save the file, click on link, then click 'Download' from the top right side.
               Citation: Silva, A., S. J. Rothstein, P. D. McNicholas, and S. Subedi (2019). A multivariate Poisson-log normal
               mixture model for clustering transcriptome sequencing data. BMC Bioinformatics. 2019;20(1):394. URL https://pubmed.ncbi.nlm.nih.gov/31311497/",
               type = "info")
  })


}

# Create Shiny app ----
shinyApp(ui, server)
