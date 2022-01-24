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
#' @return An alluvial plot should be returned. The x-axis values
#'    are in the order of vectors assigned (if any) to firstGrouping,
#'    secondGrouping, thirdGrouping and fourthGrouping, respectively.
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
#'  # is set with respect to argument firstGrouping, which is
#'  # assinged MPLNClust results.
#'
#'  set.seed(1234)
#'  alluvialPlotMPLNClust <- mplnVisualizeAlluvial(nObservations = nrow(simulatedCounts$dataset),
#'                                firstGrouping = MPLNClustResults$BICresults$BICmodelSelectedLabels,
#'                                secondGrouping = kmeans(simulatedCounts$dataset, 2)$cluster,
#'                                fileName = paste0('Plot_',date()),
#'                                printPlot = FALSE,
#'                                format = 'pdf')
#'
#'  # Note, coloring is set with respect to argument firstGrouping,
#'  # which is assinged K-means results.
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
#' Creating Alluvial Diagrams. R package version 0.1-2.
#' \href{https://github.com/mbojan/alluvial}{Link}
#'
#' Wickham, H., R. François, L. Henry and K. Müller (2021).
#' dplyr: A Grammar of Data Manipulation. R package version
#' 1.0.7. \href{https://CRAN.R-project.org/package=dplyr}{Link}
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


#' Visualize Clustered Results Via MPLN
#'
#' A function to visualize data and clustering results obtained
#' from a mixtures of multivariate Poisson-log normal (MPLN) model.
#' The function produces heatmaps and line plots of input data.
#' If cluster membership of observations are provided, this
#' information will be indicated in the figures. If a matrix of
#' probabilities for the observations belonging to each cluster
#' is provided, the option to produce a barplot of probabilities
#' is also available.
#'
#'
#' @param dataset A dataset of class matrix and type integer such that
#'    rows correspond to observations and columns correspond to variables.
#' @param clusterMembershipVector A numeric vector of length nrow(dataset)
#'    containing the cluster membership of each observation as generated by
#'    mpln(). Default is NA.
#' @param LinePlotColours Character string indicating if the line plots
#'    should be multicoloured or monotone, in black. Options are
#'    'multicolour' or 'black'. Default is 'black'.
#' @param printPlot Logical indicating if plot(s) should be saved in local
#'    directory. Default TRUE. Options TRUE or FALSE.
#' @param fileName Unique character string indicating the name for the plot
#'    being generated. Default is Plot_date, where date is obtained from
#'    date().
#' @param format Character string indicating the format of the image to
#'    be produced. Default 'pdf'. Options 'pdf' or 'png'.
#'
#' @return Plotting function provides the possibility for line and heatmap
#'    plots.
#'
#' @examples
#' # Example 1
#' # Setting parameters
#' trueMu1 <- c(6.5, 6, 6, 6, 6, 6)
#' trueMu2 <- c(2, 2.5, 2, 2, 2, 2)
#'
#' trueSigma1 <- diag(6) * 2
#' trueSigma2 <- diag(6)
#'
#' # Generate simulated data
#' simulatedCounts <- MPLNClust::mplnDataGenerator(nObservations = 100,
#'                                       dimensionality = 6,
#'                                       mixingProportions = c(0.79, 0.21),
#'                                       mu = rbind(trueMu1, trueMu2),
#'                                       sigma = rbind(trueSigma1, trueSigma2),
#'                                       produceImage = "No")
#'
#'  # Clustering data
#'  MPLNClustResults <- MPLNClust::mplnVariational(
#'                               dataset = as.matrix(simulatedCounts$dataset),
#'                               membership = "none",
#'                               gmin = 1,
#'                               gmax = 2,
#'                               initMethod = "kmeans",
#'                               nInitIterations = 1,
#'                               normalize = "Yes")
#'
#'  # Visualize data using line plot
#'  MPLNVisuals <- MPLNClust::mplnVisualizeLine(dataset = simulatedCounts$dataset,
#'                                          clusterMembershipVector =
#'                                          MPLNClustResults$allResults$`G=2`$clusterlabels,
#'                                          fileName = 'TwoClusterModel',
#'                                          printPlot = FALSE,
#'                                          format = 'png')
#'
#'  # Visualize data using line plot with multicolours
#'  # Use navigation buttons to see previous plots
#'  MPLNVisuals <- MPLNClust::mplnVisualizeLine(dataset = simulatedCounts$dataset,
#'                                          clusterMembershipVector =
#'                                          MPLNClustResults$allResults$`G=2`$clusterlabels,
#'                                          fileName = 'TwoClusterModel',
#'                                          LinePlotColours = "multicolour",
#'                                          printPlot = FALSE,
#'                                          format = 'png')
#'
#'  # Example 2
#'  # Carry out K-means clustering for same dataset as above
#'  # Use navigation buttons to see previous plots
#'  set.seed(1234)
#'  MPLNVisuals <- MPLNClust::mplnVisualizeLine(dataset = simulatedCounts$dataset,
#'                                          clusterMembershipVector = kmeans(simulatedCounts$dataset, 3)$cluster,
#'                                          fileName = 'ThreeClusterKmeansModel',
#'                                          LinePlotColours = "multicolour",
#'                                          printPlot = FALSE,
#'                                          format = 'png')
#'
#' @author Anjali Silva, \email{a.silva@utoronto.ca}
#'
#' @export
#' @import graphics
#' @importFrom grDevices png
#' @importFrom grDevices pdf
#' @importFrom grDevices dev.off
#' @importFrom RColorBrewer brewer.pal.info
#' @importFrom RColorBrewer brewer.pal
mplnVisualizeLine <- function(dataset,
                          clusterMembershipVector = NA,
                          fileName = paste0('Plot_',date()),
                          LinePlotColours = "black",
                          printPlot = TRUE,
                          format = 'pdf') {

  # Checking user input
  if (typeof(dataset) != "double" & typeof(dataset) != "integer") {
    stop("\n Dataset type needs to be integer")
  }

  if (is.matrix(dataset) != TRUE) {
    stop("\n Dataset needs to be a matrix")
  }

  if (is.logical(clusterMembershipVector) == TRUE) {
    cat("\n clusterMembershipVector is not provided.")
    clusterMembershipVector <- rep(1, nrow(dataset))

  } else if (is.numeric(clusterMembershipVector) == TRUE) {
    if (nrow(dataset) != length(clusterMembershipVector)) {
      stop("\n length(clusterMembershipVector) should match
          nrow(dataset)")
    }
  }

  # Obtaining path to save images
  pathNow <- getwd()

  # Saving cluster membership for each observation
  DataPlusLabs <- cbind(dataset, clusterMembershipVector)
  ordervector <- anothervector <- list()

  # Divide observations into each cluster based on membership
  for (i in 1:max(clusterMembershipVector)) {
    ordervector[[i]] <- which(DataPlusLabs[,
                                           ncol(dataset) + 1] == i)
    # divide observations as an integer based on cluster membership
    anothervector[[i]] <- rep(i,
                              length(which(DataPlusLabs[,
                                                        ncol(dataset) + 1] == i)))
  }

  # Setting the colours
  if(max(clusterMembershipVector) > 17) {
    qualColPals <- RColorBrewer::brewer.pal.info[brewer.pal.info$category == 'qual', ]
    coloursBarPlot <- unlist(mapply(RColorBrewer::brewer.pal,
                                    qualColPals$maxcolors,
                                    rownames(qualColPals)))
  } else {
    coloursBarPlot <- c('#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6',
                        '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324',
                        '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1',
                        '#000075', '#808080')
  }


  # empty plots
  linePlots <- NULL

  # Line Plots
  if (LinePlotColours == "multicolour") {
      linePlots <- list()
      for(cluster in unique(clusterMembershipVector)) {

        # Save how many observations below to each cluster size,
        # given by 'cluster'
        if (length(which(DataPlusLabs[, ncol(dataset) + 1] == cluster)) == 1) {
          toPlot2 <- as.matrix(DataPlusLabs[which(DataPlusLabs[,
                                                               ncol(dataset) + 1] == cluster), c(1:ncol(dataset))],
                               ncol = ncol(dataset))
          rownames(toPlot2) <- names(which(DataPlusLabs[, ncol(dataset) + 1] == cluster))
        } else if (length(which(DataPlusLabs[, ncol(dataset) + 1] == cluster)) > 1) {
          toPlot2 <- as.matrix(DataPlusLabs[which(DataPlusLabs[,
                                                               ncol(dataset) + 1] == cluster), c(1:ncol(dataset))],
                               ncol = ncol(dataset))
        }

        # Save column mean in last row
        toplot1 <- rbind(log(toPlot2 + 1), colMeans(log(toPlot2 + 1)))
        # If discontinunity is needed between samples (e.g. for 6 samples)
        # toplot1_space=cbind(toplot1[,c(1:3)],rep(NA,nrow(toPlot2)+1),
        # toplot1[,c(4:6)])


        if (printPlot == TRUE) {
          if (format == 'png') {
            grDevices::png(paste0(pathNow, "/LinePlot_Cluster", cluster,
                                  "_", fileName, ".png"))
          } else {
            grDevices::pdf(paste0(pathNow, "/LinePlot_Cluster", cluster,
                                  "_", fileName, ".pdf"))
          }

          linePlotMultiCol(dataset = dataset,
                           toplot1 = toplot1,
                           toPlot2 = toPlot2,
                           coloursBarPlot = coloursBarPlot,
                           cluster = cluster)
          grDevices::dev.off()
        }

        linePlots[[cluster]] <- linePlotMultiCol(dataset = dataset,
                                                 toplot1 = toplot1,
                                                 toPlot2 = toPlot2,
                                                 coloursBarPlot = coloursBarPlot,
                                                 cluster = cluster)
      }
    } else if (LinePlotColours == "black") {
      linePlots <- list()
      for(cluster in unique(clusterMembershipVector)) {

        # Save how many observations below to each cluster size,
        # given by 'cluster'
        if (length(which(DataPlusLabs[, ncol(dataset) + 1] == cluster)) == 1) {
          toPlot2 <- t(as.matrix(DataPlusLabs[which(DataPlusLabs[,
                                                                 ncol(dataset) + 1] == cluster), c(1:ncol(dataset))],
                                 ncol = ncol(dataset)))
          rownames(toPlot2) <- names(which(DataPlusLabs[, ncol(dataset) + 1] == cluster))
        } else if (length(which(DataPlusLabs[, ncol(dataset) + 1] == cluster)) > 1) {
          toPlot2 <- as.matrix(DataPlusLabs[which(DataPlusLabs[,
                                                               ncol(dataset) + 1] == cluster), c(1:ncol(dataset))],
                               ncol = ncol(dataset))
        }

        # Save column mean in last row
        toplot1 <- rbind(log(toPlot2 + 1), colMeans(log(toPlot2 + 1)))
        # If discontinunity is needed between samples (e.g. for 6 samples)
        # toplot1_space=cbind(toplot1[,c(1:3)],rep(NA,nrow(toPlot2)+1),
        # toplot1[,c(4:6)])

        if (printPlot == TRUE) {
          if (format == 'png') {
            grDevices::png(paste0(pathNow, "/LinePlot_Cluster", cluster,
                                  "_", fileName, ".png"))
          } else {
            grDevices::pdf(paste0(pathNow, "/LinePlot_Cluster", cluster,
                                  "_", fileName, ".pdf"))
          }
          linePlotMonoCol(dataset = dataset,
                          toplot1 = toplot1,
                          toPlot2 = toPlot2,
                          cluster = cluster)
          grDevices::dev.off()
        }
        linePlots[[cluster]] <- linePlotMonoCol(dataset = dataset,
                                                toplot1 = toplot1,
                                                toPlot2 = toPlot2,
                                                cluster = cluster)
      }
    }
  return(linePlots)
}

# Helper function
linePlotMultiCol <- function(dataset,
                             toplot1,
                             toPlot2,
                             coloursBarPlot,
                             cluster) {
  linePlotMultiCol <- graphics::matplot(t(toplot1), type = "l", pch = 1,
                                        col = c(rep(coloursBarPlot[cluster], nrow(toPlot2)), 7),
                                        xlab = "Samples", ylab = "Expression (log counts)", cex = 1,
                                        lty = c(rep(2, nrow(toPlot2)), 1),
                                        lwd = c(rep(3, nrow(toPlot2)), 4),
                                        xaxt = "n", xlim = c(1, ncol(toplot1)),
                                        main = paste("Cluster ", cluster))
  linePlotMultiCol <- linePlotMultiCol + axis(1, at = c(1:ncol(dataset)), labels = colnames(dataset))
  return(linePlotMultiCol)
}

# Helper function
linePlotMonoCol <- function(dataset,
                            toplot1,
                            toPlot2,
                            cluster) {
  linePlotMonoCol <- graphics::matplot(t(toplot1), type = "l", pch = 1,
                                       col = c(rep(1, nrow(toPlot2)), 7),
                                       xlab = "Samples", ylab = "Expression (log counts)", cex = 1,
                                       lty = c(rep(2, nrow(toPlot2)), 1),
                                       lwd = c(rep(3, nrow(toPlot2)), 4),
                                       xaxt = "n", xlim = c(1, ncol(toplot1)),
                                       main = paste("Cluster ", cluster))
  linePlotMonoCol <- linePlotMonoCol + axis(1, at = c(1:ncol(dataset)), labels = colnames(dataset))
  return(linePlotMonoCol)
}

# [END]
