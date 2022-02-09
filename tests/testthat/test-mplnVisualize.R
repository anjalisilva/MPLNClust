context("Checking for mplnVisualize performance")
library(MPLNClust)

test_that("Checking visualization via alluvial plot", {

  # Generating simulated data
  trueMu1 <- c(6.5, 6, 6, 6, 6, 6)
  trueMu2 <- c(1, 1.5, 1, 1, 1, 1)

  trueSigma1 <- diag(6) * 2
  trueSigma2 <- diag(6)

  set.seed(1234)
  simulatedCounts <- MPLNClust::mplnDataGenerator(nObservations = 100,
                                      dimensionality = 6,
                                      mixingProportions = c(0.6, 0.4),
                                      mu = rbind(trueMu1, trueMu2),
                                      sigma = rbind(trueSigma1, trueSigma2),
                                      produceImage = "No")
  set.seed(1234)
  mplnVariationalResults <- MPLNClust::mplnVariational(dataset = simulatedCounts$dataset,
                                                       membership = simulatedCounts$trueMembership,
                                                       gmin = 1,
                                                       gmax = 2,
                                                       initMethod = "kmeans",
                                                       normalize = "Yes")

  MPLNVisuals <- MPLNClust::mplnVisualizeBar(dataset = mplnVariationalResults$dataset,
                                             probabilities = mplnVariationalResults$allResults[[2]]$probaPost,
                                             clusterMembershipVector =
                                               mplnVariationalResults$allResults[[2]]$clusterlabels,
                                             fileName = 'PlotsWithProbability',
                                             printPlot = FALSE)

  expect_type(MPLNVisuals, "list")
  expect_length(MPLNVisuals, 9)
})

context("Checking for invalid user input for all plots")
test_that("Data clustering error upon invalid user input", {

  # Generating simulated data
  trueMu1 <- c(6.5, 6, 6, 6, 6, 6)
  trueMu2 <- c(2, 2.5, 2, 2, 2, 2)

  trueSigma1 <- diag(6) * 2
  trueSigma2 <- diag(6)

  set.seed(1234)
  simulatedCounts <- mplnDataGenerator(nObservations = 500,
                                       dimensionality = 6,
                                       mixingProportions = c(0.79, 0.21),
                                       mu = rbind(trueMu1, trueMu2),
                                       sigma = rbind(trueSigma1, trueSigma2),
                                       produceImage = "No")

  set.seed(1234)
  mplnVariationalResults <- MPLNClust::mplnVariational(dataset = simulatedCounts$dataset,
                                                       membership = simulatedCounts$trueMembership,
                                                       gmin = 1,
                                                       gmax = 2,
                                                       initMethod = "kmeans",
                                                       normalize = "Yes")

  # Dataset provided as character
  expect_error(MPLNClust::mplnVisualizeAlluvial(nObservations = "nrow(mplnVariationalResults$dataset)",
                                                firstGrouping = mplnVariationalResults$BICresults$BICmodelSelectedLabels,
                                                secondGrouping = mplnVariationalResults$ICLresults$ICLmodelSelectedLabels,
                                                thirdGrouping = mplnVariationalResults$AIC3results$AIC3modelSelectedLabels,
                                                fourthGrouping = mplnVariationalResults$AICresults$AICmodelSelectedLabels,
                                                fileName = paste0('Plot_',date()),
                                                printPlot = FALSE))


  # firstGrouping argument provided as logical
  expect_error(MPLNClust::mplnVisualizeAlluvial(nObservations = nrow(mplnVariationalResults$dataset),
                                                firstGrouping = TRUE,
                                                secondGrouping = mplnVariationalResults$ICLresults$ICLmodelSelectedLabels,
                                                thirdGrouping = mplnVariationalResults$AIC3results$AIC3modelSelectedLabels,
                                                fourthGrouping = mplnVariationalResults$AICresults$AICmodelSelectedLabels,
                                                fileName = paste0('Plot_',date()),
                                                printPlot = FALSE))


  # Dataset provided as wrong character
  expect_error(MPLNClust::mplnVisualizeHeatmap(dataset = "nrow(mplnVariationalResults$dataset)",
                                               clusterMembershipVector =
                                                 mplnVariationalResults$BICresults$BICmodelSelectedLabels,
                                               fileName = 'BICModel',
                                               printPlot = FALSE))


  # clusterMembershipVector length is larger than number of observations
  expect_error(MPLNClust::mplnVisualizeHeatmap(dataset = mplnVariationalResults$dataset,
                                               clusterMembershipVector =
                                                 c(mplnVariationalResults$BICresults$BICmodelSelectedLabels, 1),
                                               fileName = 'BICModel',
                                               printPlot = FALSE))

  # Dataset provided as wrong character
  expect_error(MPLNClust::mplnVisualizeLine(dataset = "mplnVariationalResults$dataset",
                                            clusterMembershipVector =
                                              mplnResults$allResults[[2]]$allResults$clusterlabels,
                                            LinePlotColours = "multicolour",
                                            fileName = 'LinePlot',
                                            printPlot = FALSE))

  # Dataset provided as wrong character
  expect_error(MPLNClust::mplnVisualizeBar(dataset = "mplnVariationalResults$dataset",
                                           probabilities = mplnVariationalResults$allResults[[2]]$probaPost,
                                           clusterMembershipVector =
                                             mplnVariationalResults$allResults[[2]]$clusterlabels,
                                           fileName = 'PlotsWithProbability',
                                           printPlot = FALSE))


  # clusterMembershipVector less than nObervations
  expect_error(MPLNClust::mplnVisualizeBar(dataset = mplnVariationalResults$dataset,
                                           probabilities = mplnVariationalResults$allResults[[2]]$probaPost,
                                           clusterMembershipVector =
                                             mplnVariationalResults$allResults[[2]]$clusterlabels[-1],
                                           fileName = 'PlotsWithProbability',
                                           printPlot = FALSE))


  # probabilities less than nObervations
  expect_error(MPLNClust::mplnVisualizeBar(dataset = mplnVariationalResults$dataset,
                                           probabilities = mplnVariationalResults$allResults[[2]]$probaPost[-1, ],
                                           clusterMembershipVector =
                                             mplnVariationalResults$allResults[[2]]$clusterlabels,
                                           fileName = 'PlotsWithProbability',
                                           printPlot = FALSE))

})

# [END]
