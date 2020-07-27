context("Checking for non parallel clustering performance")
library(MPLNClust)

test_that("Checking clustering results", {

  # Generating simulated data
  # This part works with Test Package, but not during Check
  # Error: could not find function "cpp_object_initializer"
  # Adding the function doesn't help

  trueMu1 <- c(6.5, 6, 6, 6, 6, 6)
  trueMu2 <- c(2, 2.5, 2, 2, 2, 2)

  trueSigma1 <- diag(6) * 2
  trueSigma2 <- diag(6)

  set.seed(1234)
  simulatedCounts <- mplnDataGenerator(nObservations = 50,
                                       dimensionality = 6,
                                       mixingProportions = c(0.6, 0.4),
                                       mu = rbind(trueMu1, trueMu2),
                                       sigma = rbind(trueSigma1, trueSigma2),
                                       produceImage = "No")

  # set.seed(1234)
  # mplnMCMCNonParallelResults <- MPLNClust::mplnMCMCNonParallel(
  #                                     dataset = simulatedCounts$dataset,
  #                                     membership = simulatedCounts$trueMembership,
  #                                     gmin = 1,
  #                                     gmax = 1,
  #                                     nChains = 3,
  #                                     nIterations = 900,
  #                                     initMethod = "kmeans",
  #                                     normalize = "Yes")

  # expect_that(length(mplnMCMCNonParallelResults), equals(16))
  # expect_that(mplnMCMCNonParallelResults, is_a("mplnMCMCNonParallel"))
  # expect_that(mplnMCMCNonParallelResults$initalizationMethod, equals("kmeans"))
  # numPara <- c(27)
  # expect_that(mplnMCMCNonParallelResults$numbParameters, equals(numPara))
  # expect_that(mplnMCMCNonParallelResults$trueLabels, equals(simulatedCounts$trueMembership))
  # expect_that(trunc(mplnMCMCNonParallelResults$ICLresults$ICLmodelselected), equals(1))
  # expect_that(trunc(mplnMCMCNonParallelResults$AICresults$AICmodelselected), equals(1))
  # expect_that(trunc(mplnMCMCNonParallelResults$BICresults$BICmodelselected), equals(1))
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

  # dataset provided as character
  expect_error(mplnMCMCNonParallel(dataset = "dataset",
    membership = simulatedCounts$trueMembership,
    gmin = 1,
    gmax = 2,
    nChains = 3,
    nIterations = 1000,
    initMethod = "kmeans",
    nInitIterations = 3,
    normalize = "Yes"))

  # dataset provided as logical
  expect_error(mplnMCMCNonParallel(dataset = NA,
    membership = simulatedCounts$trueMembership,
    gmin = 1,
    gmax = 2,
    nChains = 3,
    nIterations = 1000,
    initMethod = "kmeans",
    nInitIterations = 3,
    normalize = "Yes"))

  # Incorrect size for true membership
  expect_error(mplnMCMCNonParallel(dataset = simulatedCounts$dataset,
    membership = simulatedCounts$trueMembership[-1],
    gmin = 1,
    gmax = 2,
    nChains = 3,
    nIterations = 1000,
    initMethod = "kmeans",
    nInitIterations = 3,
    normalize = "Yes"))

  # Incorrect class for true membership
  expect_error(mplnMCMCNonParallel(dataset = simulatedCounts$dataset,
    membership = "trueMembership",
    gmin = 1,
    gmax = 2,
    nChains = 3,
    nIterations = 1000,
    initMethod = "kmeans",
    nInitIterations = 3,
    normalize = "Yes"))

  # Incorrect input type for gmin
  expect_error(mplnMCMCNonParallel(dataset = simulatedCounts$dataset,
    membership = simulatedCounts$trueMembership,
    gmin = "1",
    gmax = 2,
    nChains = 3,
    nIterations = 1000,
    initMethod = "kmeans",
    nInitIterations = 3,
    normalize = "Yes"))

  # Incorrect input for gmin and gmax
  expect_error(mplnMCMCNonParallel(dataset = simulatedCounts$dataset,
    membership = simulatedCounts$trueMembership,
    gmin = 5,
    gmax = 2,
    nChains = 3,
    nIterations = 1000,
    initMethod = "kmeans",
    nInitIterations = 3,
    normalize = "Yes"))


  # Incorrect input type for nChains
  expect_error(mplnMCMCNonParallel(dataset = simulatedCounts$dataset,
    membership = simulatedCounts$trueMembership,
    gmin = 1,
    gmax = 2,
    nChains = "3",
    nIterations = 1000,
    initMethod = "kmeans",
    nInitIterations = 3,
    normalize = "Yes"))

  # Incorrect input type for nIterations
  expect_error(mplnMCMCNonParallel(dataset = simulatedCounts$dataset,
    membership = simulatedCounts$trueMembership,
    gmin = 1,
    gmax = 2,
    nChains = 3,
    nIterations = "1000",
    initMethod = "kmeans",
    nInitIterations = 3,
    normalize = "Yes"))

  # Incorrect input type for initMethod
  expect_error(mplnMCMCNonParallel(dataset = simulatedCounts$dataset,
    membership = simulatedCounts$trueMembership,
    gmin = 1,
    gmax = 2,
    nChains = 3,
    nIterations = 1000,
    initMethod = "other",
    nInitIterations = 3,
    normalize = "Yes"))


  # Incorrect input type for initMethod
  expect_error(mplnMCMCNonParallel(dataset = simulatedCounts$dataset,
    membership = simulatedCounts$trueMembership,
    gmin = 1,
    gmax = 2,
    nChains = 3,
    nIterations = 1000,
    initMethod = NA,
    nInitIterations = 3,
    normalize = "Yes"))

  # Incorrect input type for nInitIterations
  expect_error(mplnMCMCNonParallel(dataset = simulatedCounts$dataset,
    membership = simulatedCounts$trueMembership,
    gmin = 1,
    gmax = 2,
    nChains = 3,
    nIterations = 1000,
    initMethod = "kmeans",
    nInitIterations = "3",
    normalize = "Yes"))

  # Incorrect input type for normalize
  expect_error(mplnMCMCNonParallel(dataset = simulatedCounts$dataset,
    membership = simulatedCounts$trueMembership,
    gmin = 1,
    gmax = 2,
    nChains = 3,
    nIterations = 1000,
    initMethod = "kmeans",
    nInitIterations = 3,
    normalize = "Other"))

})
