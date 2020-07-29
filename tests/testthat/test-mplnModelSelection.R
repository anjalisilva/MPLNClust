context("Checking for model selection performance")
library(MPLNClust)

test_that("Checking AIC model selection", {

 # Generating simulated data
 trueMu1 <- c(6.5, 6, 6, 6, 6, 6)
 trueMu2 <- c(2, 2.5, 2, 2, 2, 2)

 trueSigma1 <- diag(6) * 2
 trueSigma2 <- diag(6)

 sampleData <- mplnDataGenerator(nObservations = 100,
                                 dimensionality = 6,
                                 mixingProportions = c(0.79, 0.21),
                                 mu = rbind(trueMu1, trueMu2),
                                 sigma = rbind(trueSigma1, trueSigma2),
                                 produceImage = "No")

 # Clustering
 mplnResults <- mplnVariational(dataset = sampleData$dataset,
                                membership = sampleData$trueMembership,
                                gmin = 1,
                                gmax = 2,
                                initMethod = "kmeans",
                                nInitIterations = 2,
                                normalize = "Yes")

 # Model selection
  AICmodel <- AICFunction(logLikelihood = mplnResults$logLikelihood,
                          nParameters = mplnResults$numbParameters,
                          clusterRunOutput = mplnResults$allResults,
                          gmin = mplnResults$gmin,
                          gmax = mplnResults$gmax,
                          parallel = FALSE)

  expect_that(length(AICmodel), equals(4))
  expect_that(AICmodel, is_a("AIC"))
  expect_that(length(unique(AICmodel$AICmodelSelectedLabels)), equals(2))
  expect_that(AICmodel$allAICvalues, is_a("numeric"))
  expect_that(trunc(AICmodel$AICmodelselected), equals(2))
 })



test_that("Checking AIC3 model selection", {

 # Generating simulated data
 trueMu1 <- c(6.5, 6, 6, 6, 6, 6)
 trueMu2 <- c(2, 2.5, 2, 2, 2, 2)

 trueSigma1 <- diag(6) * 2
 trueSigma2 <- diag(6)

 sampleData <- mplnDataGenerator(nObservations = 100,
                                 dimensionality = 6,
                                 mixingProportions = c(0.79, 0.21),
                                 mu = rbind(trueMu1, trueMu2),
                                 sigma = rbind(trueSigma1, trueSigma2),
                                 produceImage = "No")


 # Clustering
 mplnResults <- mplnVariational(dataset = sampleData$dataset,
                                membership = sampleData$trueMembership,
                                gmin = 1,
                                gmax = 2,
                                initMethod = "kmeans",
                                nInitIterations = 2,
                                normalize = "Yes")

 # Model selection
 AIC3model <- AIC3Function(logLikelihood = mplnResults$logLikelihood,
                           nParameters = mplnResults$numbParameters,
                           clusterRunOutput = mplnResults$allResults,
                           gmin = mplnResults$gmin,
                           gmax = mplnResults$gmax,
                           parallel = FALSE)

 expect_that(length(AIC3model), equals(4))
 expect_that(AIC3model, is_a("AIC3"))
 expect_that(length(unique(AIC3model$AIC3modelSelectedLabels)), equals(2))
 expect_that(AIC3model$allAIC3values, is_a("numeric"))
 expect_that(trunc(AIC3model$AIC3modelselected), equals(2))
})




test_that("Checking BIC model selection", {

  # Generating simulated data
  trueMu1 <- c(6.5, 6, 6, 6, 6, 6)
  trueMu2 <- c(2, 2.5, 2, 2, 2, 2)

  trueSigma1 <- diag(6) * 2
  trueSigma2 <- diag(6)

  sampleData <- mplnDataGenerator(nObservations = 100,
                                  dimensionality = 6,
                                  mixingProportions = c(0.79, 0.21),
                                  mu = rbind(trueMu1, trueMu2),
                                  sigma = rbind(trueSigma1, trueSigma2),
                                  produceImage = "No")

  # Clustering
  mplnResults <- mplnVariational(dataset = sampleData$dataset,
                                membership = sampleData$trueMembership,
                                gmin = 1,
                                gmax = 2,
                                initMethod = "kmeans",
                                nInitIterations = 2,
                                normalize = "Yes")

  # Model selection
  BICmodel <- BICFunction(logLikelihood = mplnResults$logLikelihood,
                          nParameters = mplnResults$numbParameters,
                          nObservations = nrow(mplnResults$dataset),
                          clusterRunOutput = mplnResults$allResults,
                          gmin = mplnResults$gmin,
                          gmax = mplnResults$gmax,
                          parallel = FALSE)

  expect_that(length(BICmodel), equals(4))
  expect_that(BICmodel, is_a("BIC"))
  expect_that(length(unique(BICmodel$BICmodelSelectedLabels)), equals(2))
  expect_that(BICmodel$allBICvalues, is_a("numeric"))
  expect_that(trunc(BICmodel$BICmodelselected), equals(2))
})


test_that("Checking ICL model selection", {

  # Generating simulated data
  trueMu1 <- c(6.5, 6, 6, 6, 6, 6)
  trueMu2 <- c(2, 2.5, 2, 2, 2, 2)

  trueSigma1 <- diag(6) * 2
  trueSigma2 <- diag(6)

  sampleData <- mplnDataGenerator(nObservations = 100,
                                  dimensionality = 6,
                                  mixingProportions = c(0.79, 0.21),
                                  mu = rbind(trueMu1, trueMu2),
                                  sigma = rbind(trueSigma1, trueSigma2),
                                  produceImage = "No")

  # Clustering
  mplnResults <- mplnVariational(dataset = sampleData$dataset,
                                 membership = sampleData$trueMembership,
                                 gmin = 1,
                                 gmax = 2,
                                 initMethod = "kmeans",
                                 nInitIterations = 2,
                                 normalize = "Yes")

  BICmodel <- BICFunction(logLikelihood = mplnResults$logLikelihood,
                          nParameters = mplnResults$numbParameters,
                          nObservations = nrow(mplnResults$dataset),
                          clusterRunOutput = mplnResults$allResults,
                          gmin = mplnResults$gmin,
                          gmax = mplnResults$gmax,
                          parallel = FALSE)

 # Model selection
 ICLmodel <- ICLFunction(resultsBIC = BICmodel,
                         clusterRunOutput = mplnResults$allResults,
                         gmin = mplnResults$gmin,
                         gmax = mplnResults$gmax,
                         parallel = FALSE)

 expect_that(length(ICLmodel), equals(4))
 expect_that(ICLmodel, is_a("ICL"))
 expect_that(length(unique(ICLmodel$ICLmodelSelectedLabels)), equals(2))
 expect_that(ICLmodel$allICLvalues, is_a("numeric"))
 expect_that(trunc(ICLmodel$ICLmodelselected), equals(2))
})

