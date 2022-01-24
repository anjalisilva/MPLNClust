#' Alluvial Plot of Multiple Clustering Results
#'
#' A function to visualize clustering results via alluvial plots.
#' The function produces an alluvial plot provided multiple
#' clustering results for the same group of observations. Up to
#' four varying results could be visualized. Minimum, one
#' clustering result for visualization is required.
#'
#' @param nObservations An integer specifying the total number of
#'    observations, N, in the dataset.
#' @param firstGrouping A vector of length N, specifying the cluster
#'     membership of N observations.
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
#' trueMu1 <- c(6.5, 6, 6, 6, 6, 6)
#' trueMu2 <- c(2, 2.5, 2, 2, 2, 2)
#'
#' trueSigma1 <- diag(6) * 2
#' trueSigma2 <- diag(6)
#'
#' # Generating simulated data
#' simulatedCounts <- MPLNClust::mplnDataGenerator(nObservations = 70,
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
#'  # Visualize data
#'  MPLNVisuals <- MPLNClust::mplnVisualize(dataset = simulatedCounts$dataset,
#'                                          plots = 'all',
#'                                          probabilities =
#'                                          MPLNClustResults$allResults$`G=2`$probaPost,
#'                                          clusterMembershipVector =
#'                                          MPLNClustResults$allResults$`G=2`$clusterlabels,
#'                                          fileName = 'TwoClusterModel',
#'                                          printPlot = FALSE,
#'                                          format = 'png')
#'
#' @author Anjali Silva, \email{anjali.silva@uhnresearch.ca}
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
mplnVisualizeAlluvial <- function(nObservations,
                                  firstGrouping,
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

  # Saving cluster membership for each observation

  if(length(secondGrouping) == 0){

  }
  toVisualize <- data.frame(Observation = c(1:nObservations),
                            Method1 = floor(runif(50, min=1, max=8)),
                            Method2 = floor(runif(50, min=1, max=5)),
                            Method3 = floor(runif(50, min=1, max=4)),
                            Method3 = floor(runif(50, min=1, max=4)),
                            Freq = rep(1, 50))
  testing$Method1[c(1:3)] <- 8
  testing$Method2[c(1:3)] <- 5
  testing$Method3[c(1:3)] <- 4

  alluvial(testing[, c(1:4)],
           freq=testing$Freq, border=NA, alpha = 0.5,
           col=case_when(testing$Method1 == "1" ~ "red",
                         testing$Method1 == "2" ~ "blue",
                         testing$Method1 == "3" ~ "purple",
                         testing$Method1 == "4" ~ "black",
                         testing$Method1 == "5" ~ "green",
                         testing$Method1 == "6" ~ "yellow",
                         testing$Method1 == "7" ~ "magenta",
                         testing$Method1 == "8" ~ "darkblue",
                         TRUE ~ "orange"),
           cex=0.75,
           axis_labels = c("Observation", "G=8", "G=5", "G=4"))







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
  heatmapOne <- heatmapTwo <- linePlots <- barPlot <- NULL

  if (plots == 'all' || plots == 'heatmaps') {

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
  }

  if (plots == 'all' || plots == 'lines') {
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
  }

  if (plots == 'all' || plots == 'bar') {

    if(is.logical(probabilities) == TRUE){
      stop("\n probabilities should be provided to make bar plot.")
    }

    # Bar plot
    tableProbabilities <- as.data.frame(cbind(Sample = c(1:nrow(probabilities)),
                                              Cluster = mclust::map(probabilities),
                                              probabilities))

    names(tableProbabilities) <- c("Sample", "Cluster",
                                   paste0("P", rep(1:(ncol(tableProbabilities)-2))))

    tableProbabilitiesMelt <- reshape::melt(tableProbabilities,
                                            id.vars = c("Sample","Cluster"))

    if (printPlot == TRUE) {
      barPlot <- barPlotFunction(tableProbabilitiesMelt = tableProbabilitiesMelt,
                                 coloursBarPlot = coloursBarPlot,
                                 probabilities = probabilities)
      ggplot2::ggsave(paste0(pathNow,"/barplot_", fileName,".",format))
    }

    barPlot <- barPlotFunction(tableProbabilitiesMelt = tableProbabilitiesMelt,
                               coloursBarPlot = coloursBarPlot,
                               probabilities = probabilities)
  }

  return(list(heatmapOne,
              heatmapTwo,
              linePlots,
              barPlot))
}
