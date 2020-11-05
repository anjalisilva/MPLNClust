context("Checking for mplnVisualize performance")
library(MPLNClust)

test_that("Checking visualization", {

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

  MPLNVisuals <- MPLNClust::mplnVisualize(dataset = simulatedCounts$dataset,
                                          plots = 'all',
                                          probabilities = mplnVariationalResults$allResults$`G=2`$probaPost,
                                          clusterMembershipVector =
                                            mplnVariationalResults$allResults$`G=2`$clusterlabels,
                                          fileName = 'TwoClusterModel',
                                          printPlot = FALSE,
                                          format = 'png')

  expect_type(MPLNVisuals, "list")
  expect_length(MPLNVisuals, 4)
})

context("Checking for invalid user input")
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
  expect_error(mplnVisualize(dataset = "mplnVariationalResults$dataset",
                             plots = 'all',
                             probabilities = mplnVariationalResults$allResults$`G=2`$probaPost,
                             clusterMembershipVector =
                               mplnVariationalResults$allResults$`G=2`$clusterlabels,
                             fileName = 'TwoClusterModel',
                             printPlot = FALSE,
                             format = 'png'))

  # Plots provided as logical
  expect_error(mplnVisualize(dataset = mplnVariationalResults$dataset,
                             plots = TRUE,
                             probabilities = mplnVariationalResults$allResults$`G=2`$probaPost,
                             clusterMembershipVector =
                               mplnVariationalResults$allResults$`G=2`$clusterlabels,
                             fileName = 'TwoClusterModel',
                             printPlot = FALSE,
                             format = 'png'))

  # Plots provided as wrong character
  expect_error(mplnVisualize(dataset = mplnVariationalResults$dataset,
                             plots = "allplots",
                             probabilities = mplnVariationalResults$allResults$`G=2`$probaPost,
                             clusterMembershipVector =
                               mplnVariationalResults$allResults$`G=2`$clusterlabels,
                             fileName = 'TwoClusterModel',
                             printPlot = FALSE,
                             format = 'png'))

  # clusterMembershipVector less than nObervations
  expect_error(mplnVisualize(dataset = mplnVariationalResults$dataset,
                             plots = "heatmaps",
                             probabilities = mplnVariationalResults$allResults$`G=2`$probaPost,
                             clusterMembershipVector =
                               mplnVariationalResults$allResults$`G=2`$clusterlabels[-1],
                             fileName = 'TwoClusterModel',
                             printPlot = FALSE,
                             format = 'png'))

  # probabilities less than nObervations
  expect_error(mplnVisualize(dataset = mplnVariationalResults$dataset,
                             plots = "heatmaps",
                             probabilities = mplnVariationalResults$allResults$`G=2`$probaPost[-1, ],
                             clusterMembershipVector =
                               mplnVariationalResults$allResults$`G=2`$clusterlabels,
                             fileName = 'TwoClusterModel',
                             printPlot = FALSE,
                             format = 'png'))

})
# [END]
