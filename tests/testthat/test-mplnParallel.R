context("Checking for parallel clustering performance")
library(MPLNClust)

test_that("Checking clustering results", {

  # Generating simulated data
  trueMu1 <- c(6.5, 6, 6, 6, 6, 6)
  trueMu2 <- c(2, 2.5, 2, 2, 2, 2)

  trueSigma1 <- diag(6) * 2
  trueSigma2 <- diag(6)

  set.seed(1234)
  simulated_counts <- mplnDataGenerator(nObservations = 50,
                                        dimensionality = 6,
                                        mixingProportions = c(0.79, 0.21),
                                        mu = rbind(trueMu1, trueMu2),
                                        sigma = rbind(trueSigma1, trueSigma2),
                                        produceImage = "No")

   mplnResults <- MPLNClust::mplnParallel(dataset = simulated_counts$dataset,
                                membership = simulated_counts$trueMembership,
                                gmin = 1,
                                gmax = 2,
                                nChains = 3,
                                nIterations = 500,
                                initMethod = "kmeans",
                                nInitIterations = 2,
                                normalize = "Yes")

  expect_that(length(mplnResults), equals(16))
  expect_that(mplnResults, is_a("mplnParallel"))
  expect_that(mplnResults$initalization_method, equals("kmeans"))
  normFactors <- c(-0.38090382,  0.04043787,  0.21095606,
    0.53566041, -0.44720551,  0.04105498)
  expect_that(mplnResults$normalization_factors, equals(normFactors))
  numPara <- c(27, 55)
  expect_that(mplnResults$numb_of_parameters, equals(numPara))
  expect_that(mplnResults$true_labels, equals(simulated_counts$trueMembership))
  expect_that(trunc(mplnResults$ICL_all$ICLmodelselected), equals(2))
  expect_that(trunc(mplnResults$AIC_all$AICmodelselected), equals(2))
  expect_that(trunc(mplnResults$BIC_all$BICmodelselected), equals(2))
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
  expect_error(mplnParallel(dataset = "dataset",
    membership = simulated_counts$trueMembership,
    gmin = 1,
    gmax = 2,
    nChains = 3,
    nIterations = 1000,
    initMethod = "kmeans",
    nInitIterations = 3,
    normalize = "Yes"))

  # dataset provided as logical
  expect_error(mplnParallel(dataset = NA,
    membership = simulated_counts$trueMembership,
    gmin = 1,
    gmax = 2,
    nChains = 3,
    nIterations = 1000,
    initMethod = "kmeans",
    nInitIterations = 3,
    normalize = "Yes"))

  # Incorrect size for true membership
  expect_error(mplnParallel(dataset = simulated_counts$dataset,
    membership = simulated_counts$trueMembership[-1],
    gmin = 1,
    gmax = 2,
    nChains = 3,
    nIterations = 1000,
    initMethod = "kmeans",
    nInitIterations = 3,
    normalize = "Yes"))

  # Incorrect class for true membership
  expect_error(mplnParallel(dataset = simulated_counts$dataset,
    membership = "trueMembership",
    gmin = 1,
    gmax = 2,
    nChains = 3,
    nIterations = 1000,
    initMethod = "kmeans",
    nInitIterations = 3,
    normalize = "Yes"))

  # Incorrect input type for gmin
  expect_error(mplnParallel(dataset = simulated_counts$dataset,
    membership = simulated_counts$trueMembership,
    gmin = "1",
    gmax = 2,
    nChains = 3,
    nIterations = 1000,
    initMethod = "kmeans",
    nInitIterations = 3,
    normalize = "Yes"))

  # Incorrect input for gmin and gmax
  expect_error(mplnParallel(dataset = simulated_counts$dataset,
    membership = simulated_counts$trueMembership,
    gmin = 5,
    gmax = 2,
    nChains = 3,
    nIterations = 1000,
    initMethod = "kmeans",
    nInitIterations = 3,
    normalize = "Yes"))


  # Incorrect input type for nChains
  expect_error(mplnParallel(dataset = simulated_counts$dataset,
    membership = simulated_counts$trueMembership,
    gmin = 1,
    gmax = 2,
    nChains = "3",
    nIterations = 1000,
    initMethod = "kmeans",
    nInitIterations = 3,
    normalize = "Yes"))

  # Incorrect input type for nIterations
  expect_error(mplnParallel(dataset = simulated_counts$dataset,
    membership = simulated_counts$trueMembership,
    gmin = 1,
    gmax = 2,
    nChains = 3,
    nIterations = "1000",
    initMethod = "kmeans",
    nInitIterations = 3,
    normalize = "Yes"))

  # Incorrect input type for initMethod
  expect_error(mplnParallel(dataset = simulated_counts$dataset,
    membership = simulated_counts$trueMembership,
    gmin = 1,
    gmax = 2,
    nChains = 3,
    nIterations = 1000,
    initMethod = "other",
    nInitIterations = 3,
    normalize = "Yes"))


  # Incorrect input type for initMethod
  expect_error(mplnParallel(dataset = simulated_counts$dataset,
    membership = simulated_counts$trueMembership,
    gmin = 1,
    gmax = 2,
    nChains = 3,
    nIterations = 1000,
    initMethod = NA,
    nInitIterations = 3,
    normalize = "Yes"))

  # Incorrect input type for nInitIterations
  expect_error(mplnParallel(dataset = simulated_counts$dataset,
    membership = simulated_counts$trueMembership,
    gmin = 1,
    gmax = 2,
    nChains = 3,
    nIterations = 1000,
    initMethod = "kmeans",
    nInitIterations = "3",
    normalize = "Yes"))

  # Incorrect input type for normalize
  expect_error(mplnParallel(dataset = simulated_counts$dataset,
    membership = simulated_counts$trueMembership,
    gmin = 1,
    gmax = 2,
    nChains = 3,
    nIterations = 1000,
    initMethod = "kmeans",
    nInitIterations = 3,
    normalize = "Other"))

})
