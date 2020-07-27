context("Checking for parallel clustering performance")
library(MPLNClust)

test_that("Checking clustering results", {

  # Generating simulated data
  trueMu1 <- c(6.5, 6, 6, 6, 6, 6)
  trueMu2 <- c(1, 1.5, 1, 1, 1, 1)

  trueSigma1 <- diag(6) * 2
  trueSigma2 <- diag(6)

  set.seed(1234)
  simulatedCounts <- mplnDataGenerator(nObservations = 40,
                                      dimensionality = 6,
                                      mixingProportions = c(0.6, 0.4),
                                      mu = rbind(trueMu1, trueMu2),
                                      sigma = rbind(trueSigma1, trueSigma2),
                                      produceImage = "No")
  set.seed(1234)
  mplnMCMCResults <- MPLNClust::mplnMCMCParallel(dataset = simulatedCounts$dataset,
                                                 membership = simulatedCounts$trueMembership,
                                                 gmin = 2,
                                                 gmax = 2,
                                                 nChains = 3,
                                                 nIterations = 400,
                                                 initMethod = "kmeans",
                                                 normalize = "Yes",
                                                 numNodes = 2)
  # Setting numNodes = 2 based on the following entry, otherwise error.
  # "NB: you can’t use unexported functions and you shouldn’t open new graphics
  # devices or use more than two cores. Individual examples shouldn’t
  # take more than 5s."
  # https://stackoverflow.com/questions/41307178/error-processing-vignette-failed-with-diagnostics-4-simultaneous-processes-spa

  expect_that(length(mplnMCMCResults), equals(16))
  expect_that(mplnMCMCResults, is_a("mplnMCMCParallel"))
  expect_that(mplnMCMCResults$initalizationMethod, equals("kmeans"))
  numPara <- c(27)
  expect_that(mplnMCMCResults$numbParameters, equals(numPara))
  expect_that(mplnMCMCResults$trueLabels, equals(simulatedCounts$trueMembership))
  expect_that(trunc(mplnMCMCResults$ICLresults$ICLmodelselected), equals(2))
  expect_that(trunc(mplnMCMCResults$AICresults$AICmodelselected), equals(2))
  expect_that(trunc(mplnMCMCResults$AIC3results$AIC3modelselected), equals(2))
  expect_that(trunc(mplnMCMCResults$BICresults$BICmodelselected), equals(2))
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
  expect_error(mplnMCMCParallel(dataset = "dataset",
    membership = simulatedCounts$trueMembership,
    gmin = 1,
    gmax = 2,
    nChains = 3,
    nIterations = 1000,
    initMethod = "kmeans",
    nInitIterations = 3,
    normalize = "Yes"))

  # dataset provided as logical
  expect_error(mplnMCMCParallel(dataset = NA,
    membership = simulatedCounts$trueMembership,
    gmin = 1,
    gmax = 2,
    nChains = 3,
    nIterations = 1000,
    initMethod = "kmeans",
    nInitIterations = 3,
    normalize = "Yes"))

  # Incorrect size for true membership
  expect_error(mplnMCMCParallel(dataset = simulatedCounts$dataset,
    membership = simulatedCounts$trueMembership[-1],
    gmin = 1,
    gmax = 2,
    nChains = 3,
    nIterations = 1000,
    initMethod = "kmeans",
    nInitIterations = 3,
    normalize = "Yes"))

  # Incorrect class for true membership
  expect_error(mplnMCMCParallel(dataset = simulatedCounts$dataset,
    membership = "trueMembership",
    gmin = 1,
    gmax = 2,
    nChains = 3,
    nIterations = 1000,
    initMethod = "kmeans",
    nInitIterations = 3,
    normalize = "Yes"))

  # Incorrect input type for gmin
  expect_error(mplnMCMCParallel(dataset = simulatedCounts$dataset,
    membership = simulatedCounts$trueMembership,
    gmin = "1",
    gmax = 2,
    nChains = 3,
    nIterations = 1000,
    initMethod = "kmeans",
    nInitIterations = 3,
    normalize = "Yes"))

  # Incorrect input for gmin and gmax
  expect_error(mplnMCMCParallel(dataset = simulatedCounts$dataset,
    membership = simulatedCounts$trueMembership,
    gmin = 5,
    gmax = 2,
    nChains = 3,
    nIterations = 1000,
    initMethod = "kmeans",
    nInitIterations = 3,
    normalize = "Yes"))


  # Incorrect input type for nChains
  expect_error(mplnMCMCParallel(dataset = simulatedCounts$dataset,
    membership = simulatedCounts$trueMembership,
    gmin = 1,
    gmax = 2,
    nChains = "3",
    nIterations = 1000,
    initMethod = "kmeans",
    nInitIterations = 3,
    normalize = "Yes"))

  # Incorrect input type for nIterations
  expect_error(mplnMCMCParallel(dataset = simulatedCounts$dataset,
    membership = simulatedCounts$trueMembership,
    gmin = 1,
    gmax = 2,
    nChains = 3,
    nIterations = "1000",
    initMethod = "kmeans",
    nInitIterations = 3,
    normalize = "Yes"))

  # Incorrect input type for initMethod
  expect_error(mplnMCMCParallel(dataset = simulatedCounts$dataset,
    membership = simulatedCounts$trueMembership,
    gmin = 1,
    gmax = 2,
    nChains = 3,
    nIterations = 1000,
    initMethod = "other",
    nInitIterations = 3,
    normalize = "Yes"))


  # Incorrect input type for initMethod
  expect_error(mplnMCMCParallel(dataset = simulatedCounts$dataset,
    membership = simulatedCounts$trueMembership,
    gmin = 1,
    gmax = 2,
    nChains = 3,
    nIterations = 1000,
    initMethod = NA,
    nInitIterations = 3,
    normalize = "Yes"))

  # Incorrect input type for nInitIterations
  expect_error(mplnMCMCParallel(dataset = simulatedCounts$dataset,
    membership = simulatedCounts$trueMembership,
    gmin = 1,
    gmax = 2,
    nChains = 3,
    nIterations = 1000,
    initMethod = "kmeans",
    nInitIterations = "3",
    normalize = "Yes"))

  # Incorrect input type for normalize
  expect_error(mplnMCMCParallel(dataset = simulatedCounts$dataset,
    membership = simulatedCounts$trueMembership,
    gmin = 1,
    gmax = 2,
    nChains = 3,
    nIterations = 1000,
    initMethod = "kmeans",
    nInitIterations = 3,
    normalize = "Other"))

})

