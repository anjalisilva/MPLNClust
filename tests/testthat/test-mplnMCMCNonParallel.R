context("Checking for non parallel clustering performance")
library(MPLNClust)

test_that("Checking clustering results", {

  # Generating simulated data
  trueMu1 <- c(6.5, 6, 6, 6, 6, 6)
  trueMu2 <- c(2, 2.5, 2, 2, 2, 2)

  trueSigma1 <- diag(6) * 2
  trueSigma2 <- diag(6)

  set.seed(1234)
  simulated_counts <- mplnDataGenerator(nObservations = 40,
                                        dimensionality = 6,
                                        mixingProportions = c(0.79, 0.21),
                                        mu = rbind(trueMu1, trueMu2),
                                        sigma = rbind(trueSigma1, trueSigma2),
                                        produceImage = "No")
  set.seed(1234)
  mplnMCMCNonParallelResults <- MPLNClust::mplnMCMCNonParallel(
                                dataset = simulated_counts$dataset,
                                membership = simulated_counts$trueMembership,
                                gmin = 2,
                                gmax = 2,
                                nChains = 3,
                                nIterations = 400,
                                initMethod = "kmeans",
                                normalize = "Yes")

  expect_that(length(mplnMCMCNonParallelResults), equals(16))
  expect_that(mplnMCMCNonParallelResults, is_a("mplnMCMCNonParallel"))
  expect_that(mplnMCMCNonParallelResults$initalization_method, equals("kmeans"))
  numPara <- c(27)
  expect_that(mplnMCMCNonParallelResults$numb_of_parameters, equals(numPara))
  expect_that(mplnMCMCNonParallelResults$true_labels, equals(simulated_counts$trueMembership))
  expect_that(trunc(mplnMCMCNonParallelResults$ICL_all$ICLmodelselected), equals(2))
  expect_that(trunc(mplnMCMCNonParallelResults$AIC_all$AICmodelselected), equals(2))
  expect_that(trunc(mplnMCMCNonParallelResults$BIC_all$BICmodelselected), equals(2))
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
  expect_error(mplnMCMCNonParallel(dataset = "dataset",
                               membership = simulated_counts$trueMembership,
                               gmin = 1,
                               gmax = 2,
                               nChains = 3,
                               nIterations = 1000,
                               initMethod = "kmeans",
                               nInitIterations = 3,
                               normalize = "Yes"))

  # dataset provided as logical
  expect_error(mplnMCMCNonParallel(dataset = NA,
                               membership = simulated_counts$trueMembership,
                               gmin = 1,
                               gmax = 2,
                               nChains = 3,
                               nIterations = 1000,
                               initMethod = "kmeans",
                               nInitIterations = 3,
                               normalize = "Yes"))

  # Incorrect size for true membership
  expect_error(mplnMCMCNonParallel(dataset = simulated_counts$dataset,
                               membership = simulated_counts$trueMembership[-1],
                               gmin = 1,
                               gmax = 2,
                               nChains = 3,
                               nIterations = 1000,
                               initMethod = "kmeans",
                               nInitIterations = 3,
                               normalize = "Yes"))

  # Incorrect class for true membership
  expect_error(mplnMCMCNonParallel(dataset = simulated_counts$dataset,
                               membership = "trueMembership",
                               gmin = 1,
                               gmax = 2,
                               nChains = 3,
                               nIterations = 1000,
                               initMethod = "kmeans",
                               nInitIterations = 3,
                               normalize = "Yes"))

  # Incorrect input type for gmin
  expect_error(mplnMCMCNonParallel(dataset = simulated_counts$dataset,
                               membership = simulated_counts$trueMembership,
                               gmin = "1",
                               gmax = 2,
                               nChains = 3,
                               nIterations = 1000,
                               initMethod = "kmeans",
                               nInitIterations = 3,
                               normalize = "Yes"))

  # Incorrect input for gmin and gmax
  expect_error(mplnMCMCNonParallel(dataset = simulated_counts$dataset,
                               membership = simulated_counts$trueMembership,
                               gmin = 5,
                               gmax = 2,
                               nChains = 3,
                               nIterations = 1000,
                               initMethod = "kmeans",
                               nInitIterations = 3,
                               normalize = "Yes"))


  # Incorrect input type for nChains
  expect_error(mplnMCMCNonParallel(dataset = simulated_counts$dataset,
                               membership = simulated_counts$trueMembership,
                               gmin = 1,
                               gmax = 2,
                               nChains = "3",
                               nIterations = 1000,
                               initMethod = "kmeans",
                               nInitIterations = 3,
                               normalize = "Yes"))

  # Incorrect input type for nIterations
  expect_error(mplnMCMCNonParallel(dataset = simulated_counts$dataset,
                               membership = simulated_counts$trueMembership,
                               gmin = 1,
                               gmax = 2,
                               nChains = 3,
                               nIterations = "1000",
                               initMethod = "kmeans",
                               nInitIterations = 3,
                               normalize = "Yes"))

  # Incorrect input type for initMethod
  expect_error(mplnMCMCNonParallel(dataset = simulated_counts$dataset,
                               membership = simulated_counts$trueMembership,
                               gmin = 1,
                               gmax = 2,
                               nChains = 3,
                               nIterations = 1000,
                               initMethod = "other",
                               nInitIterations = 3,
                               normalize = "Yes"))


  # Incorrect input type for initMethod
  expect_error(mplnMCMCNonParallel(dataset = simulated_counts$dataset,
                               membership = simulated_counts$trueMembership,
                               gmin = 1,
                               gmax = 2,
                               nChains = 3,
                               nIterations = 1000,
                               initMethod = NA,
                               nInitIterations = 3,
                               normalize = "Yes"))

  # Incorrect input type for nInitIterations
  expect_error(mplnMCMCNonParallel(dataset = simulated_counts$dataset,
                                membership = simulated_counts$trueMembership,
                                gmin = 1,
                                gmax = 2,
                                nChains = 3,
                                nIterations = 1000,
                                initMethod = "kmeans",
                                nInitIterations = "3",
                                normalize = "Yes"))

  # Incorrect input type for normalize
  expect_error(mplnMCMCNonParallel(dataset = simulated_counts$dataset,
                                membership = simulated_counts$trueMembership,
                                gmin = 1,
                                gmax = 2,
                                nChains = 3,
                                nIterations = 1000,
                                initMethod = "kmeans",
                                nInitIterations = 3,
                                normalize = "Other"))

})




