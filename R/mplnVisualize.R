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


#' Visualize Clustered Results Via Line Plots
#'
#' A function to visualize clustering results via line plots.
#' Each cluster will have its own plot with mean of the
#' data in each cluster shown via a yellow line.
#'
#' @param dataset A dataset of class matrix and type integer such that
#'    rows correspond to observations and columns correspond to variables.
#' @param clusterMembershipVector A numeric vector of length nrow(dataset)
#'    containing the cluster membership of each observation. If not provided,
#'    all observations will be treated as belonging to one cluster. Default is NA.
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
#'  MPLNLineBlack <- MPLNClust::mplnVisualizeLine(dataset = simulatedCounts$dataset,
#'                                          clusterMembershipVector =
#'                                          MPLNClustResults$allResults$`G=2`$clusterlabels,
#'                                          fileName = 'TwoClusterModel',
#'                                          printPlot = FALSE,
#'                                          format = 'png')
#'
#'  # Visualize data using line plot with multicolours
#'  # Use navigation buttons to see previous plots
#'  MPLNLineColor <- MPLNClust::mplnVisualizeLine(dataset = simulatedCounts$dataset,
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
#'  KmeansLineColor <- MPLNClust::mplnVisualizeLine(dataset = simulatedCounts$dataset,
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

# Helper function for line plot
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

# Helper function for line plot
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





#' Visualize Posterior Probabilities via Bar Plot
#'
#' A function to produce a barplot of posterior probabilities
#' for each observation, after clustering via mixtures of
#' multivariate Poisson-log normal (MPLN) model.
#'
#' @param dataset A dataset of class matrix and type integer such that
#'    rows correspond to observations (N) and columns (C) correspond to variables.
#' @param probabilities A matrix of size N x C, such that rows correspond
#'    to N observations and columns correspond to C clusters. Each row
#'    should sum to 1. Default is NA.
#' @param clusterMembershipVector A numeric vector of length nrow(dataset)
#'    containing the cluster membership of each observation. Default is NA.
#' @param printPlot Logical indicating if plot(s) should be saved in local
#'    directory. Default TRUE. Options TRUE or FALSE.
#' @param fileName Unique character string indicating the name for the plot
#'    being generated. Default is Plot_date, where date is obtained from
#'    date().
#' @param format Character string indicating the format of the image to
#'    be produced. Default 'pdf'. Options 'pdf' or 'png'.
#'
#' @return A bar plot of posterior probabilities for each observation
#'    in the dataset.
#'
#' @examples
#' # Example 1
#' trueMu1 <- c(6.5, 6, 6, 6, 6, 6)
#' trueMu2 <- c(2, 2.5, 2, 2, 2, 2)
#'
#' trueSigma1 <- diag(6) * 2
#' trueSigma2 <- diag(6)
#'
#' # Generating simulated data
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
#'  # Visualize posterior probabilities via a bar plot
#'  MPLNVisuals <- MPLNClust::mplnVisualizeBar(dataset = simulatedCounts$dataset,
#'                                             probabilities =
#'                                             MPLNClustResults$allResults$`G=2`$probaPost,
#'                                             clusterMembershipVector =
#'                                             MPLNClustResults$allResults$`G=2`$clusterlabels,
#'                                             fileName = 'PosteriorProbMPLN',
#'                                             printPlot = FALSE,
#'                                             format = 'png')
#'
#' @author Anjali Silva, \email{a.silva@utoronto.ca}
#'
#' @export
#' @import graphics
#' @import ggplot2
#' @importFrom grDevices png
#' @importFrom grDevices pdf
#' @importFrom grDevices dev.off
#' @importFrom RColorBrewer brewer.pal.info
#' @importFrom RColorBrewer brewer.pal
#' @importFrom randomcoloR distinctColorPalette
#' @importFrom reshape melt
mplnVisualizeBar <- function(dataset,
                          probabilities = NA,
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

  if (is.logical(probabilities) == TRUE) {
    cat("\n Probabilities are not provided. Barplot of probabilities will not be produced.")
  } else if (is.matrix(probabilities) == TRUE) {
    if (nrow(probabilities) != length(clusterMembershipVector)) {
      stop("\n length(probabilities) should match nrow(dataset)")
    }
    if (any(rowSums(probabilities) >= 1.01)) {
      stop("\n rowSums(probabilities) reveals at least
          one observation has probability != 1.")
    }
    if (any(rowSums(probabilities) <= 0.99)) {
      stop("\n rowSums(probabilities) reveals at least
          one observation has probability != 1.")
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

  vec <- unlist(ordervector) # put observations in order of cluster membership
  colorsvector <- unlist(anothervector) # put all details together as integers

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

  # empty plot
  barPlot <- NULL

  # Bar plot
  tableProbabilities <- as.data.frame(cbind(Sample = c(1:nrow(probabilities)),
                                            Cluster = mclust::map(probabilities),
                                            probabilities))

  names(tableProbabilities) <- c("Sample", "Cluster",
                                 paste0("P", rep(1:(ncol(tableProbabilities) - 2))))

  tableProbabilitiesMelt <- reshape::melt(tableProbabilities,
                                          id.vars = c("Sample", "Cluster"))

  if (printPlot == TRUE) {
      barPlot <- barPlotFunction(tableProbabilitiesMelt = tableProbabilitiesMelt,
                                 coloursBarPlot = coloursBarPlot,
                                 probabilities = probabilities)
      ggplot2::ggsave(paste0(pathNow, "/barplot_", fileName, ".", format))
    }

  barPlot <- barPlotFunction(tableProbabilitiesMelt = tableProbabilitiesMelt,
                               coloursBarPlot = coloursBarPlot,
                               probabilities = probabilities)

  return(barPlot)
}

# Helper function for bar plot
barPlotFunction <- function(tableProbabilitiesMelt,
                            coloursBarPlot,
                            probabilities) {

  variable <- value <- Sample <- NULL

  if(is.data.frame(tableProbabilitiesMelt) != TRUE) {
    stop("tableProbabilitiesMelt should be a data frame")
  }

  if(is.character(coloursBarPlot) != TRUE) {
    stop("coloursBarPlot should be character")
  }

  if(is.matrix(probabilities) != TRUE) {
    stop("probabilities should be a matrix")
  }

  barPlot <- ggplot2::ggplot(data = tableProbabilitiesMelt,
                             ggplot2::aes(fill = variable, y = value, x = Sample))

  barPlot <-
    barPlot +
    ggplot2::geom_bar(position = "fill", stat = "identity") +
    scale_fill_manual(values = coloursBarPlot,
                      name = "Cluster") + theme_bw() +
    theme(text = element_text(size = 10),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(face = "bold"),
          axis.text.y = element_text(face = "bold")) +
    coord_cartesian(ylim = c(0, 1), xlim = c(1, nrow(probabilities))) +
    labs(x = "Observation") +
    scale_y_continuous(name = "Posterior probability", limits = c(0: 1))

  return(barPlot)
}







#' Visualize Clustered Results Via Heatmaps
#'
#' A function to produce heatmaps of data with clustering results.
#'
#' @param dataset A dataset of class matrix and type integer such that
#'    rows correspond to observations (N) and columns (C) correspond to variables.
#' @param probabilities A matrix of size N x C, such that rows correspond
#'    to N observations and columns correspond to C clusters. Each row
#'    should sum to 1. Default is NA.
#' @param clusterMembershipVector A numeric vector of length nrow(dataset)
#'    containing the cluster membership of each observation as generated by
#'    mpln(). Default is NA.
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
#' trueMu1 <- c(6.5, 6, 6, 6, 6, 6)
#' trueMu2 <- c(2, 2.5, 2, 2, 2, 2)
#'
#' trueSigma1 <- diag(6) * 2
#' trueSigma2 <- diag(6)
#'
#' # Generating simulated data
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
#'  # Visualize data via a Heatmap
#'  MPLNVisuals <- MPLNClust::mplnVisualizeHeatmap(dataset = simulatedCounts$dataset,
#'                                          MPLNClustResults$allResults$`G=2`$probaPost,
#'                                          clusterMembershipVector =
#'                                          MPLNClustResults$allResults$`G=2`$clusterlabels,
#'                                          fileName = 'TwoClusterModel',
#'                                          printPlot = FALSE,
#'                                          format = 'png')
#'
#' @author Anjali Silva, \email{a.silva@utoronto.ca}
#'
#' @export
#' @import graphics
#' @import ggplot2
#' @importFrom grDevices png
#' @importFrom grDevices pdf
#' @importFrom grDevices dev.off
#' @importFrom RColorBrewer brewer.pal.info
#' @importFrom RColorBrewer brewer.pal
#' @importFrom randomcoloR distinctColorPalette
#' @importFrom pheatmap pheatmap
#' @importFrom gplots heatmap.2
#' @importFrom gplots redgreen
#' @importFrom reshape melt
mplnVisualizeHeatmap <- function(dataset,
                          clusterMembershipVector = NA,
                          fileName = paste0('Plot_',date()),
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

  vec <- unlist(ordervector) # put observations in order of cluster membership
  colorsvector <- unlist(anothervector) # put all details together as integers

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
  heatmapOne <- heatmapTwo <- NULL

  # Heatmap 1
  if (printPlot == TRUE) {

      if (format == 'png') {
        grDevices::png(paste0(pathNow, "/heatmap1_", fileName, ".png"))
      } else {
        grDevices::pdf(paste0(pathNow, "/heatmap1_", fileName, ".pdf"))
      }

      gplots::heatmap.2(as.matrix(dataset[vec, ]),
                        dendrogram = "column",
                        trace = "none",
                        scale = "row",
                        Rowv = FALSE, col = rev(gplots::redgreen(75)),
                        RowSideColor = coloursBarPlot[colorsvector],
                        labRow = FALSE,
                        main = paste("Clustering results, G =",
                                     max(clusterMembershipVector)))
      graphics::par(xpd = TRUE)
      graphics::legend(xpd = TRUE, x = -0.1, y = 0,
                       legend = paste0("Cluster ", unique(colorsvector)),
                       col = unique(coloursBarPlot[colorsvector]),
                       lty = 1,
                       lwd = 5,
                       cex =.8, horiz = FALSE)

      grDevices::dev.off()
    }


    # Heatmap 2
    # Only produced if less than 17 clusters
    if(max(clusterMembershipVector) < 18) {
      # Defining annotation row
      annotation_row = data.frame(Cluster = factor(clusterMembershipVector[vec]))
      if(is.null(rownames(dataset)) == TRUE) {
        rownames(dataset)  = paste("Gene", c(1:nrow(dataset[vec, ])))
        rownames(annotation_row) = rownames(dataset[vec, ])
      } else {
        rownames(annotation_row) = rownames(dataset[vec, ])
      }

      # Define row annotation colours
      heatMap2RowAnnotation <- c("1" = coloursBarPlot[1], "2" = coloursBarPlot[2],
                                 "3" = coloursBarPlot[3], "4" = coloursBarPlot[4],
                                 "5" = coloursBarPlot[5], "6" = coloursBarPlot[6],
                                 "7" = coloursBarPlot[7], "8" = coloursBarPlot[8],
                                 "9" = coloursBarPlot[9], "10" = coloursBarPlot[10],
                                 "11" = coloursBarPlot[11], "12" = coloursBarPlot[12],
                                 "13" = coloursBarPlot[13], "14" = coloursBarPlot[14],
                                 "15" = coloursBarPlot[15], "16" = coloursBarPlot[16],
                                 "17" = coloursBarPlot[17])

      # Show row names or not based on dataset size
      if(nrow(dataset) < 50){
        showLabels = TRUE
      } else {
        showLabels = FALSE
      }


      if (printPlot == TRUE) {
        if (format == 'png') {
          grDevices::png(paste0(pathNow, "/heatmap2_", fileName, ".png"))
        } else {
          grDevices::pdf(paste0(pathNow, "/heatmap2_", fileName, ".pdf"))
        }

        heatmapFunctionTwo(dataset = dataset,
                           vec = vec,
                           showLabels = showLabels,
                           heatMap2RowAnnotation = heatMap2RowAnnotation,
                           annotation_row = annotation_row,
                           clusterMembershipVector = clusterMembershipVector)
        grDevices::dev.off()
      }


      heatmapTwo <- heatmapFunctionTwo(dataset = dataset,
                                       vec = vec,
                                       showLabels = showLabels,
                                       heatMap2RowAnnotation = heatMap2RowAnnotation,
                                       annotation_row = annotation_row,
                                       clusterMembershipVector = clusterMembershipVector)
    }

  return(list(heatmapOne,
              heatmapTwo))
}

# [END]
