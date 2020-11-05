context("Checking for mplnVariational performance")
library(MPLNClust)

test_that("Checking clustering results", {

  # Generating simulated data
  trueMu1 <- c(6.5, 6, 6, 6, 6, 6)
  trueMu2 <- c(1, 1.5, 1, 1, 1, 1)

  trueSigma1 <- diag(6) * 2
  trueSigma2 <- diag(6)

  set.seed(1234)
  simulatedCounts <- mplnDataGenerator(nObservations = 1000,
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

  expect_type(mplnVariationalResults, "list")
  expect_s3_class(mplnVariationalResults, "mplnVariational")
  expect_length(mplnVariationalResults, 16)
  expect_identical(mplnVariationalResults$initalizationMethod, "kmeans")
  numPara <- c(27, 55)
  expect_identical(mplnVariationalResults$numbParameters, numPara)
  expect_identical(mplnVariationalResults$trueLabels, simulatedCounts$trueMembership)
  expect_identical(trunc(mplnVariationalResults$ICLresults$ICLmodelselected), 2)
  expect_identical(trunc(mplnVariationalResults$AICresults$AICmodelselected), 2)
  expect_identical(trunc(mplnVariationalResults$AIC3results$AIC3modelselected), 2)
  expect_identical(trunc(mplnVariationalResults$BICresults$BICmodelselected), 2)
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
  expect_error(mplnVariational(dataset = "dataset",
                              membership = simulatedCounts$trueMembership,
                              gmin = 1,
                              gmax = 2,
                              initMethod = "kmeans",
                              nInitIterations = 3,
                              normalize = "Yes"))

  # dataset provided as logical
  expect_error(mplnVariational(dataset = NA,
                                membership = simulatedCounts$trueMembership,
                                gmin = 1,
                                gmax = 2,
                                initMethod = "kmeans",
                                nInitIterations = 3,
                                normalize = "Yes"))

  # Incorrect size for true membership
  expect_error(mplnVariational(dataset = simulatedCounts$dataset,
                              membership = simulatedCounts$trueMembership[-1],
                              gmin = 1,
                              gmax = 2,
                              initMethod = "kmeans",
                              nInitIterations = 3,
                              normalize = "Yes"))

  # Incorrect class for true membership
  expect_error(mplnVariational(dataset = simulatedCounts$dataset,
                              membership = "trueMembership",
                              gmin = 1,
                              gmax = 2,
                              initMethod = "kmeans",
                              nInitIterations = 3,
                              normalize = "Yes"))

  # Incorrect input type for gmin
  expect_error(mplnVariational(dataset = simulatedCounts$dataset,
                              membership = simulatedCounts$trueMembership,
                              gmin = "1",
                              gmax = 2,
                              initMethod = "kmeans",
                              nInitIterations = 3,
                              normalize = "Yes"))

  # Incorrect input for gmin and gmax
  expect_error(mplnVariational(dataset = simulatedCounts$dataset,
                              membership = simulatedCounts$trueMembership,
                              gmin = 5,
                              gmax = 2,
                              initMethod = "kmeans",
                              nInitIterations = 3,
                              normalize = "Yes"))



  # Incorrect input type for initMethod
  expect_error(mplnVariational(dataset = simulatedCounts$dataset,
              membership = simulatedCounts$trueMembership,
              gmin = 1,
              gmax = 2,
              initMethod = "other",
              nInitIterations = 3,
              normalize = "Yes"))


  # Incorrect input type for initMethod
  expect_error(mplnVariational(dataset = simulatedCounts$dataset,
              membership = simulatedCounts$trueMembership,
              gmin = 1,
              gmax = 2,
              initMethod = NA,
              nInitIterations = 3,
              normalize = "Yes"))

  # Incorrect input type for nInitIterations
  expect_error(mplnVariational(dataset = simulatedCounts$dataset,
              membership = simulatedCounts$trueMembership,
              gmin = 1,
              gmax = 2,
              initMethod = "kmeans",
              nInitIterations = "3",
              normalize = "Yes"))

  # Incorrect input type for normalize
  expect_error(mplnVariational(dataset = simulatedCounts$dataset,
              membership = simulatedCounts$trueMembership,
              gmin = 1,
              gmax = 2,
              initMethod = "kmeans",
              nInitIterations = 3,
              normalize = "Other"))
})
# [END]
