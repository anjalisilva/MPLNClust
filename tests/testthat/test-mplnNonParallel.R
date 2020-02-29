context("Checking for non parallel clustering performance")
library(MPLNClust)

test_that("Checking clustering results", {

  # Generating simulated data
  trueMu1 <- c(6.5, 6, 6, 6, 6, 6)
  trueMu2 <- c(2, 2.5, 2, 2, 2, 2)

  trueSigma1 <- diag(6) * 2
  trueSigma2 <- diag(6)

  set.seed(1234)
  simulated_counts <- mplnDataGenerator(nObservations = 70,
                                        dimensionality = 6,
                                        mixingProportions = c(0.79, 0.21),
                                        mu = rbind(trueMu1, trueMu2),
                                        sigma = rbind(trueSigma1, trueSigma2),
                                        produceImage = "No")
  set.seed(1234)
  mplnNonParallelResults <- MPLNClust::mplnNonParallel(
                                dataset = simulated_counts$dataset,
                                membership = simulated_counts$trueMembership,
                                gmin = 1,
                                gmax = 2,
                                nChains = 3,
                                nIterations = 500,
                                initMethod = "kmeans",
                                nInitIterations = 1,
                                normalize = "Yes")

  expect_that(length(mplnNonParallelResults), equals(16))
  expect_that(mplnNonParallelResults, is_a("mplnNonParallel"))
  expect_that(mplnNonParallelResults$initalization_method, equals("kmeans"))
  normFactors <- c(-0.62996925, -0.15999746,  0.12852236,
    0.34015676,  0.28591670,  0.03537089)
  expect_that(mplnNonParallelResults$normalization_factors, equals(normFactors))
  numPara <- c(27, 55)
  expect_that(mplnNonParallelResults$numb_of_parameters, equals(numPara))
  expect_that(mplnNonParallelResults$true_labels, equals(simulated_counts$trueMembership))
  expect_that(trunc(mplnNonParallelResults$ICL_all$ICLmodelselected), equals(2))
  expect_that(trunc(mplnNonParallelResults$AIC_all$AICmodelselected), equals(2))
  expect_that(trunc(mplnNonParallelResults$BIC_all$BICmodelselected), equals(2))
})

context("Checking for invalid user input")
test_that("Data clustering error upon invalid user input", {

  # Generating simulated data
  trueMu1 <- c(6.5, 6, 6, 6, 6, 6)
  trueMu2 <- c(2, 2.5, 2, 2, 2, 2)

  trueSigma1 <- diag(6) * 2
  trueSigma2 <- diag(6)

  set.seed(1234)
  simulated_counts <- mplnDataGenerator(nObservations = 500,
    dimensionality = 6,
    mixingProportions = c(0.79, 0.21),
    mu = rbind(trueMu1, trueMu2),
    sigma = rbind(trueSigma1, trueSigma2),
    produceImage = "No")

  # dataset provided as character
  expect_error(mplnNonParallel(dataset = "dataset",
    membership = simulated_counts$trueMembership,
    gmin = 1,
    gmax = 2,
    nChains = 3,
    nIterations = 1000,
    initMethod = "kmeans",
    nInitIterations = 3,
    normalize = "Yes"))

  # dataset provided as logical
  expect_error(mplnNonParallel(dataset = NA,
    membership = simulated_counts$trueMembership,
    gmin = 1,
    gmax = 2,
    nChains = 3,
    nIterations = 1000,
    initMethod = "kmeans",
    nInitIterations = 3,
    normalize = "Yes"))

  # Incorrect size for true membership
  expect_error(mplnNonParallel(dataset = simulated_counts$dataset,
    membership = simulated_counts$trueMembership[-1],
    gmin = 1,
    gmax = 2,
    nChains = 3,
    nIterations = 1000,
    initMethod = "kmeans",
    nInitIterations = 3,
    normalize = "Yes"))

  # Incorrect class for true membership
  expect_error(mplnNonParallel(dataset = simulated_counts$dataset,
    membership = "trueMembership",
    gmin = 1,
    gmax = 2,
    nChains = 3,
    nIterations = 1000,
    initMethod = "kmeans",
    nInitIterations = 3,
    normalize = "Yes"))

  # Incorrect input type for gmin
  expect_error(mplnNonParallel(dataset = simulated_counts$dataset,
    membership = simulated_counts$trueMembership,
    gmin = "1",
    gmax = 2,
    nChains = 3,
    nIterations = 1000,
    initMethod = "kmeans",
    nInitIterations = 3,
    normalize = "Yes"))

  # Incorrect input for gmin and gmax
  expect_error(mplnNonParallel(dataset = simulated_counts$dataset,
    membership = simulated_counts$trueMembership,
    gmin = 5,
    gmax = 2,
    nChains = 3,
    nIterations = 1000,
    initMethod = "kmeans",
    nInitIterations = 3,
    normalize = "Yes"))


  # Incorrect input type for nChains
  expect_error(mplnNonParallel(dataset = simulated_counts$dataset,
    membership = simulated_counts$trueMembership,
    gmin = 1,
    gmax = 2,
    nChains = "3",
    nIterations = 1000,
    initMethod = "kmeans",
    nInitIterations = 3,
    normalize = "Yes"))

  # Incorrect input type for nIterations
  expect_error(mplnNonParallel(dataset = simulated_counts$dataset,
    membership = simulated_counts$trueMembership,
    gmin = 1,
    gmax = 2,
    nChains = 3,
    nIterations = "1000",
    initMethod = "kmeans",
    nInitIterations = 3,
    normalize = "Yes"))

  # Incorrect input type for initMethod
  expect_error(mplnNonParallel(dataset = simulated_counts$dataset,
    membership = simulated_counts$trueMembership,
    gmin = 1,
    gmax = 2,
    nChains = 3,
    nIterations = 1000,
    initMethod = "other",
    nInitIterations = 3,
    normalize = "Yes"))


  # Incorrect input type for initMethod
  expect_error(mplnNonParallel(dataset = simulated_counts$dataset,
    membership = simulated_counts$trueMembership,
    gmin = 1,
    gmax = 2,
    nChains = 3,
    nIterations = 1000,
    initMethod = NA,
    nInitIterations = 3,
    normalize = "Yes"))

  # Incorrect input type for nInitIterations
  expect_error(mplnNonParallel(dataset = simulated_counts$dataset,
    membership = simulated_counts$trueMembership,
    gmin = 1,
    gmax = 2,
    nChains = 3,
    nIterations = 1000,
    initMethod = "kmeans",
    nInitIterations = "3",
    normalize = "Yes"))

  # Incorrect input type for normalize
  expect_error(mplnNonParallel(dataset = simulated_counts$dataset,
    membership = simulated_counts$trueMembership,
    gmin = 1,
    gmax = 2,
    nChains = 3,
    nIterations = 1000,
    initMethod = "kmeans",
    nInitIterations = 3,
    normalize = "Other"))

})




