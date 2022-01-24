#' Alluvial Plot of Multiple Clustering Results
#'
#' A function to visualize clustering results via alluvial plots,
#' using the alluvial::alluvial() function. The function produces
#' an alluvial plot provided multiple clustering results for
#' the same group of observations. Up to four varying results
#' could be visualized. Minimum, one clustering result for
#' visualization is required. Maximum 10 colors (clusters) are
#' supported. Colors are assigned based on cluster membership
#' assigned for argument 'firstGrouping'.
#'
#' @param nObservations An integer specifying the total number of
#'    observations, N, in the dataset. Default value is 50L.
#' @param firstGrouping A vector of length nObservations, specifying the
#'     cluster membership of observations. This must be provided. Colors
#'     will be assigned based on cluster membership provided in this
#'     vector. Default value is a vector of length 50.
#' @param secondGrouping A vector of length N, specifying the cluster
#'     membership of N observations. This could be obtained via another
#'     clustering run or from a different model selection criteria.
#'     Default value is an empty vector.
#' @param thirdGrouping A vector of length N, specifying the cluster
#'     membership of N observations. This could be obtained via another
#'     clustering run or from a different model selection criteria.
#'     Default value is an empty vector.
#' @param fourthGrouping A vector of length N, specifying the cluster
#'     membership of N observations. This could be obtained via another
#'     clustering run or from a different model selection criteria.
#'     Default value is an empty vector.
#' @param printPlot Logical indicating if plot(s) should be saved in local
#'    directory. Default TRUE. Options TRUE or FALSE.
#' @param fileName Unique character string indicating the name for the plot
#'    being generated. Default is Plot_date, where date is obtained from
#'    date().
#' @param format Character string indicating the format of the image to
#'    be produced. Default 'pdf'. Options 'pdf' or 'png'.
#'
#' @return An alluvial plot should be returned.
#'
#' @examples
#' # Example 1
#' # Assign values for parameters
#' trueMu1 <- c(6.5, 6, 6, 6, 6, 6)
#' trueMu2 <- c(2, 2.5, 2, 2, 2, 2)
#'
#' trueSigma1 <- diag(6) * 2
#' trueSigma2 <- diag(6)
#'
#' # Generate simulated data for 500 x 6 dataset
#' simulatedCounts <- MPLNClust::mplnDataGenerator(nObservations = 500,
#'                                       dimensionality = 6,
#'                                       mixingProportions = c(0.79, 0.21),
#'                                       mu = rbind(trueMu1, trueMu2),
#'                                       sigma = rbind(trueSigma1, trueSigma2),
#'                                       produceImage = "No")
#'
#'  # Clustering data for G = 1:2
#'  MPLNClustResults <- MPLNClust::mplnVariational(
#'                               dataset = as.matrix(simulatedCounts$dataset),
#'                               membership = "none",
#'                               gmin = 1,
#'                               gmax = 2,
#'                               initMethod = "kmeans",
#'                               nInitIterations = 1,
#'                               normalize = "Yes")
#'
#'  # Visualize clustering results using alluvial plot
#'  # Access results using models selected via model selection criteria
#'  alluvialPlot <- mplnVisualizeAlluvial(nObservations = nrow(simulatedCounts$dataset),
#'                            firstGrouping = MPLNClustResults$BICresults$BICmodelSelectedLabels,
#'                            secondGrouping = MPLNClustResults$ICLresults$ICLmodelSelectedLabels,
#'                            thirdGrouping = MPLNClustResults$AIC3results$AIC3modelSelectedLabels,
#'                            fourthGrouping = MPLNClustResults$AICresults$AICmodelSelectedLabels,
#'                            fileName = paste0('Plot_',date()),
#'                            printPlot = FALSE,
#'                            format = 'pdf')
#'
#'  # Example 2
#'  # Perform clustering via K-means with centers = 2
#'  # Visualize clustering results using alluvial plot for
#'  # K-means and above MPLNClust results. Note, coloring
#'  # is set with respect to firstGrouping which is assinged
#'  # MPLNClust results.
#'
#'  set.seed(1234)
#'  alluvialPlotMPLNClust <- mplnVisualizeAlluvial(nObservations = nrow(simulatedCounts$dataset),
#'                                firstGrouping = MPLNClustResults$BICresults$BICmodelSelectedLabels,
#'                                secondGrouping = kmeans(simulatedCounts$dataset, 2)$cluster,
#'                                fileName = paste0('Plot_',date()),
#'                                printPlot = FALSE,
#'                                format = 'pdf')
#'
#'  # Note, coloring is set with respect to firstGrouping which is
#'  # assinged K-means results.
#'  set.seed(1234)
#'  alluvialPlotKmeans <- mplnVisualizeAlluvial(nObservations = nrow(simulatedCounts$dataset),
#'                            firstGrouping = kmeans(simulatedCounts$dataset, 2)$cluster,
#'                            secondGrouping = MPLNClustResults$BICresults$BICmodelSelectedLabels,
#'                            fileName = paste0('Plot_',date()),
#'                            printPlot = FALSE,
#'                            format = 'pdf')
#'
#' @author Anjali Silva, \email{a.silva@utoronto.ca}
#'
#' @references
#' Bojanowski,  M., R. Edwards (2016). alluvial: R Package for
#' Creating Alluvial Diagrams. R package version: 0.1-2,
#' \href{https://github.com/mbojan/alluvial}{Link}
#'
#' @export
#' @import graphics
#' @import alluvial
#' @importFrom dplyr case_when
#' @importFrom grDevices png
#' @importFrom grDevices pdf
#' @importFrom grDevices dev.off
mplnVisualizeAlluvial <- function(nObservations = 50L,
                                  firstGrouping = floor(runif(50, min = 1, max = 8)),
                                  secondGrouping = vector(mode = "integer", length = 0),
                                  thirdGrouping = vector(mode = "integer", length = 0),
                                  fourthGrouping = vector(mode = "integer", length = 0),
                                  fileName = paste0('Plot_',date()),
                                  printPlot = TRUE,
                                  format = 'pdf') {

  # Checking user input
  if (typeof(nObservations) != "integer") {
    stop("\n nObservations should be an integer")
  }

  if (is.vector(firstGrouping) != TRUE) {
    stop("\n firstGrouping should be a vector")
  }

  if (is.vector(secondGrouping) != TRUE) {
    stop("\n secondGrouping should be a vector")
  }

  if (is.vector(thirdGrouping) != TRUE) {
    stop("\n thirdGrouping should be a vector")
  }

  if (is.vector(fourthGrouping) != TRUE) {
    stop("\n fourthGrouping should be a vector")
  }

  # Obtaining path to save images
  pathNow <- getwd()

  # Setting colors based on first grouping
  coloursBarPlot <- c('#4363d8', '#f58231', '#911eb4',
                      '#46f0f0', '#f032e6',
                      '#bcf60c', '#fabebe', '#008080',
                      '#e6beff', '#9a6324',
                      '#fffac8', '#800000', '#aaffc3',
                      '#808000', '#ffd8b1',
                      '#000075', '#808080')
  setVectorColor <- firstGrouping
  if(max(firstGrouping) == 2) {
    colSetting <- dplyr::case_when(
                    setVectorColor == "1" ~ coloursBarPlot[1],
                    setVectorColor == "2" ~ coloursBarPlot[2],
                    TRUE ~ "orange")
  } else if(max(firstGrouping) == 3) {
    colSetting <- dplyr::case_when(
                      setVectorColor == "1" ~ coloursBarPlot[1],
                      setVectorColor == "2" ~ coloursBarPlot[2],
                      setVectorColor == "3" ~ coloursBarPlot[3],
                      TRUE ~ "orange")
  } else if(max(firstGrouping) == 4) {
    colSetting <- dplyr::case_when(
      setVectorColor == "1" ~ coloursBarPlot[1],
      setVectorColor == "2" ~ coloursBarPlot[2],
      setVectorColor == "2" ~ coloursBarPlot[3],
      setVectorColor == "4" ~ coloursBarPlot[4],
      TRUE ~ "orange")
  } else if(max(firstGrouping) == 5) {
    colSetting <- dplyr::case_when(
      setVectorColor == "1" ~ coloursBarPlot[1],
      setVectorColor == "2" ~ coloursBarPlot[2],
      setVectorColor == "2" ~ coloursBarPlot[3],
      setVectorColor == "4" ~ coloursBarPlot[4],
      setVectorColor == "5" ~ coloursBarPlot[5],
      TRUE ~ "orange")
  } else if(max(firstGrouping) == 6) {
    colSetting <- dplyr::case_when(
      setVectorColor == "1" ~ coloursBarPlot[1],
      setVectorColor == "2" ~ coloursBarPlot[2],
      setVectorColor == "2" ~ coloursBarPlot[3],
      setVectorColor == "4" ~ coloursBarPlot[4],
      setVectorColor == "5" ~ coloursBarPlot[5],
      setVectorColor == "6" ~ coloursBarPlot[6],
      TRUE ~ "orange")
  } else if(max(firstGrouping) == 7) {
    colSetting <- dplyr::case_when(
      setVectorColor == "1" ~ coloursBarPlot[1],
      setVectorColor == "2" ~ coloursBarPlot[2],
      setVectorColor == "2" ~ coloursBarPlot[3],
      setVectorColor == "4" ~ coloursBarPlot[4],
      setVectorColor == "5" ~ coloursBarPlot[5],
      setVectorColor == "6" ~ coloursBarPlot[6],
      setVectorColor == "7" ~ coloursBarPlot[7],
      TRUE ~ "orange")
  } else if(max(firstGrouping) == 8) {
    colSetting <- dplyr::case_when(
      setVectorColor == "1" ~ coloursBarPlot[1],
      setVectorColor == "2" ~ coloursBarPlot[2],
      setVectorColor == "2" ~ coloursBarPlot[3],
      setVectorColor == "4" ~ coloursBarPlot[4],
      setVectorColor == "5" ~ coloursBarPlot[5],
      setVectorColor == "6" ~ coloursBarPlot[6],
      setVectorColor == "7" ~ coloursBarPlot[7],
      setVectorColor == "8" ~ coloursBarPlot[8],
      TRUE ~ "orange")
  } else if(max(firstGrouping) == 9) {
    colSetting <- dplyr::case_when(
      setVectorColor == "1" ~ coloursBarPlot[1],
      setVectorColor == "2" ~ coloursBarPlot[2],
      setVectorColor == "3" ~ coloursBarPlot[3],
      setVectorColor == "4" ~ coloursBarPlot[4],
      setVectorColor == "5" ~ coloursBarPlot[5],
      setVectorColor == "6" ~ coloursBarPlot[6],
      setVectorColor == "7" ~ coloursBarPlot[7],
      setVectorColor == "8" ~ coloursBarPlot[8],
      setVectorColor == "9" ~ coloursBarPlot[9],
      TRUE ~ "orange")
  } else if(max(firstGrouping) == 10) {
    colSetting <- dplyr::case_when(
      setVectorColor == "1" ~ coloursBarPlot[1],
      setVectorColor == "2" ~ coloursBarPlot[2],
      setVectorColor == "3" ~ coloursBarPlot[3],
      setVectorColor == "4" ~ coloursBarPlot[4],
      setVectorColor == "5" ~ coloursBarPlot[5],
      setVectorColor == "6" ~ coloursBarPlot[6],
      setVectorColor == "7" ~ coloursBarPlot[7],
      setVectorColor == "8" ~ coloursBarPlot[8],
      setVectorColor == "9" ~ coloursBarPlot[9],
      setVectorColor == "10" ~ coloursBarPlot[10],
      TRUE ~ "orange")
  }


  # Saving cluster membership for each observation
  if(length(secondGrouping) == 0) {
    toVisualize <- data.frame(Observation = c(1:nObservations),
                              Method1 = firstGrouping,
                              Freq = rep(1, nObservations))

    plotAlluvial <- alluvial::alluvial(toVisualize[, c(1:2)],
                            freq = toVisualize$Freq,
                            border = NA, alpha = 0.5,
                            col = colSetting,
                            cex = 0.75,
                            axis_labels = c("Observation", paste0("G=", max(firstGrouping))))

  } else if(length(thirdGrouping) == 0) {
    toVisualize <- data.frame(Observation = c(1:nObservations),
                              Method1 = firstGrouping,
                              Method2 = secondGrouping,
                              Freq = rep(1, nObservations))

    plotAlluvial <- alluvial::alluvial(toVisualize[, c(1:3)],
                                       freq = toVisualize$Freq,
                                       border = NA, alpha = 0.5,
                                       col = colSetting,
                                       cex = 0.75,
                                       axis_labels = c("Observation",
                                                       paste0("G=", max(firstGrouping)),
                                                       paste0("G=", max(secondGrouping))))

  } else if(length(fourthGrouping) == 0) {
    toVisualize <- data.frame(Observation = c(1:nObservations),
                              Method1 = firstGrouping,
                              Method2 = secondGrouping,
                              Method3 = thirdGrouping,
                              Freq = rep(1, nObservations))

    plotAlluvial <- alluvial::alluvial(toVisualize[, c(1:4)],
                                       freq = toVisualize$Freq,
                                       border = NA, alpha = 0.5,
                                       col = colSetting,
                                       cex = 0.75,
                                       axis_labels = c("Observation",
                                                       paste0("G=", max(firstGrouping)),
                                                       paste0("G=", max(secondGrouping)),
                                                       paste0("G=", max(thirdGrouping))))

  } else {
    toVisualize <- data.frame(Observation = c(1:nObservations),
                              Method1 = firstGrouping,
                              Method2 = secondGrouping,
                              Method3 = thirdGrouping,
                              Method4 = fourthGrouping,
                              Freq = rep(1, nObservations))

    plotAlluvial <- alluvial::alluvial(toVisualize[, c(1:5)],
                                       freq = toVisualize$Freq,
                                       border = NA, alpha = 0.5,
                                       col = colSetting,
                                       cex = 0.75,
                                       axis_labels = c("Observation",
                                                       paste0("G=", max(firstGrouping)),
                                                       paste0("G=", max(secondGrouping)),
                                                       paste0("G=", max(thirdGrouping)),
                                                       paste0("G=", max(fourthGrouping))))

  }

  # printing plots
  if (printPlot == TRUE) {
    if (format == 'png') {
        grDevices::png(paste0(pathNow, "/AlluvialPlot_", fileName, ".png"))
    } else {
        grDevices::pdf(paste0(pathNow, "/AlluvialPlot_", fileName, ".pdf"))
    }
      plotAlluvial
      grDevices::dev.off()
    }

  class(plotAlluvial) <- "mplnAlluvialVisual"
  return(plotAlluvial)
}
