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

  expect_type(AICmodel, "list")
  expect_s3_class(AICmodel, "AIC")
  expect_length(AICmodel, 4)
  expect_length(unique(AICmodel$AICmodelSelectedLabels), 2)
  expect_identical(class(AICmodel$allAICvalues), "numeric")
  expect_identical(trunc(AICmodel$AICmodelselected), 2)
 })
context("Checking for invalid user input for AIC")
test_that("AIC model selection error upon invalid user input", {

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

  # logLikelihood provided as character
  expect_error(AICFunction(logLikelihood = "mplnResults$logLikelihood",
                          nParameters = mplnResults$numbParameters,
                          clusterRunOutput = mplnResults$allResults,
                          gmin = mplnResults$gmin,
                          gmax = mplnResults$gmax,
                          parallel = FALSE))

  # nParameters provided as character
  expect_error(AICFunction(logLikelihood = mplnResults$logLikelihood,
                           nParameters = "mplnResults$numbParameters",
                           clusterRunOutput = mplnResults$allResults,
                           gmin = mplnResults$gmin,
                           gmax = mplnResults$gmax,
                           parallel = FALSE))

  # gmin provided as character
  expect_error(AICFunction(logLikelihood = mplnResults$logLikelihood,
                          nParameters = mplnResults$numbParameters,
                          clusterRunOutput = mplnResults$allResults,
                          gmin = "mplnResults$gmin",
                          gmax = mplnResults$gmax,
                          parallel = FALSE))

  # gmax provided as character
  expect_error(AICFunction(logLikelihood = mplnResults$logLikelihood,
                            nParameters = mplnResults$numbParameters,
                            clusterRunOutput = mplnResults$allResults,
                            gmin = mplnResults$gmin,
                            gmax = "mplnResults$gmax",
                            parallel = FALSE))


  # parallel provided as character
  expect_error(AICFunction(logLikelihood = mplnResults$logLikelihood,
                          nParameters = mplnResults$numbParameters,
                          clusterRunOutput = mplnResults$allResults,
                          gmin = mplnResults$gmin,
                          gmax = mplnResults$gmax,
                          parallel = "TRUE"))

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


 expect_type(AIC3model, "list")
 expect_s3_class(AIC3model, "AIC3")
 expect_length(AIC3model, 4)
 expect_length(unique(AIC3model$AIC3modelSelectedLabels), 2)
 expect_identical(class(AIC3model$allAIC3values), "numeric")
 expect_identical(trunc(AIC3model$AIC3modelselected), 2)
})


context("Checking for invalid user input for AIC3")
test_that("AIC3 model selection error upon invalid user input", {

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

  # logLikelihood provided as character
  expect_error(AIC3Function(logLikelihood = "mplnResults$logLikelihood",
                            nParameters = mplnResults$numbParameters,
                            clusterRunOutput = mplnResults$allResults,
                            gmin = mplnResults$gmin,
                            gmax = mplnResults$gmax,
                            parallel = FALSE))

  # nParameters provided as character
  expect_error(AIC3Function(logLikelihood = mplnResults$logLikelihood,
                            nParameters = "mplnResults$numbParameters",
                            clusterRunOutput = mplnResults$allResults,
                            gmin = mplnResults$gmin,
                            gmax = mplnResults$gmax,
                            parallel = FALSE))

  # gmin provided as character
  expect_error(AIC3Function(logLikelihood = mplnResults$logLikelihood,
                            nParameters = mplnResults$numbParameters,
                            clusterRunOutput = mplnResults$allResults,
                            gmin = "mplnResults$gmin",
                            gmax = mplnResults$gmax,
                            parallel = FALSE))

  # gmax provided as character
  expect_error(AIC3Function(logLikelihood = mplnResults$logLikelihood,
                            nParameters = mplnResults$numbParameters,
                            clusterRunOutput = mplnResults$allResults,
                            gmin = mplnResults$gmin,
                            gmax = "mplnResults$gmax",
                            parallel = FALSE))


  # parallel provided as character
  expect_error(AIC3Function(logLikelihood = mplnResults$logLikelihood,
                            nParameters = mplnResults$numbParameters,
                            clusterRunOutput = mplnResults$allResults,
                            gmin = mplnResults$gmin,
                            gmax = mplnResults$gmax,
                            parallel = "TRUE"))

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


  expect_type(BICmodel, "list")
  expect_s3_class(BICmodel, "BIC")
  expect_length(BICmodel, 4)
  expect_length(unique(BICmodel$BICmodelSelectedLabels), 2)
  expect_identical(class(BICmodel$allBICvalues), "numeric")
  expect_identical(trunc(BICmodel$BICmodelselected), 2)
})

context("Checking for invalid user input for BIC")
test_that("BIC model selection error upon invalid user input", {

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

  # logLikelihood provided as character
  expect_error(BICFunction(logLikelihood = "mplnResults$logLikelihood",
                           nParameters = mplnResults$numbParameters,
                           nObservations = nrow(mplnResults$dataset),
                           clusterRunOutput = mplnResults$allResults,
                           gmin = mplnResults$gmin,
                           gmax = mplnResults$gmax,
                           parallel = FALSE))

  # nParameters provided as character
  expect_error(BICFunction(logLikelihood = mplnResults$logLikelihood,
                          nParameters = "mplnResults$numbParameters",
                          nObservations = nrow(mplnResults$dataset),
                          clusterRunOutput = mplnResults$allResults,
                          gmin = mplnResults$gmin,
                          gmax = mplnResults$gmax,
                          parallel = FALSE))

  # gmin provided as character
  expect_error(BICFunction(logLikelihood = mplnResults$logLikelihood,
                          nParameters = mplnResults$numbParameters,
                          nObservations = nrow(mplnResults$dataset),
                          clusterRunOutput = mplnResults$allResults,
                          gmin = "mplnResults$gmin",
                          gmax = mplnResults$gmax,
                          parallel = FALSE))

  # gmax provided as character
  expect_error(BICFunction(logLikelihood = mplnResults$logLikelihood,
                          nParameters = mplnResults$numbParameters,
                          nObservations = nrow(mplnResults$dataset),
                          clusterRunOutput = mplnResults$allResults,
                          gmin = mplnResults$gmin,
                          gmax = "mplnResults$gmax",
                          parallel = FALSE))


  # parallel provided as character
  expect_error(BICFunction(logLikelihood = mplnResults$logLikelihood,
                          nParameters = mplnResults$numbParameters,
                          nObservations = nrow(mplnResults$dataset),
                          clusterRunOutput = mplnResults$allResults,
                          gmin = mplnResults$gmin,
                          gmax = mplnResults$gmax,
                          parallel = "TRUE"))

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

 # Model selection
  ICLmodel <- ICLFunction(logLikelihood = mplnResults$logLikelihood,
                          nParameters = mplnResults$numbParameters,
                          nObservations = nrow(mplnResults$dataset),
                          clusterRunOutput = mplnResults$allResults,
                          gmin = mplnResults$gmin,
                          gmax = mplnResults$gmax,
                          parallel = FALSE)


  expect_type(ICLmodel, "list")
  expect_s3_class(ICLmodel, "ICL")
  expect_length(ICLmodel, 4)
  expect_length(unique(ICLmodel$ICLmodelSelectedLabels), 2)
  expect_identical(class(ICLmodel$allICLvalues), "numeric")
  expect_identical(trunc(ICLmodel$ICLmodelselected), 2)
})

context("Checking for invalid user input for ICL")
test_that("ICL model selection error upon invalid user input", {

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

  # logLikelihood provided as character
  expect_error(ICLFunction(logLikelihood = "mplnResults$logLikelihood",
                            nParameters = mplnResults$numbParameters,
                            nObservations = nrow(mplnResults$dataset),
                            clusterRunOutput = mplnResults$allResults,
                            gmin = mplnResults$gmin,
                            gmax = mplnResults$gmax,
                            parallel = FALSE))

  # nParameters provided as character
  expect_error(ICLFunction(logLikelihood = mplnResults$logLikelihood,
                          nParameters = "mplnResults$numbParameters",
                          nObservations = nrow(mplnResults$dataset),
                          clusterRunOutput = mplnResults$allResults,
                          gmin = mplnResults$gmin,
                          gmax = mplnResults$gmax,
                          parallel = FALSE))

  # gmin provided as character
  expect_error(ICLFunction(logLikelihood = mplnResults$logLikelihood,
                            nParameters = mplnResults$numbParameters,
                            nObservations = nrow(mplnResults$dataset),
                            clusterRunOutput = mplnResults$allResults,
                            gmin = "mplnResults$gmin",
                            gmax = mplnResults$gmax,
                            parallel = FALSE))

  # gmax provided as character
  expect_error(ICLFunction(logLikelihood = mplnResults$logLikelihood,
                            nParameters = mplnResults$numbParameters,
                            nObservations = nrow(mplnResults$dataset),
                            clusterRunOutput = mplnResults$allResults,
                            gmin = mplnResults$gmin,
                            gmax = "mplnResults$gmax",
                            parallel = FALSE))


  # parallel provided as character
  expect_error(ICLFunction(logLikelihood = mplnResults$logLikelihood,
                          nParameters = mplnResults$numbParameters,
                          nObservations = nrow(mplnResults$dataset),
                          clusterRunOutput = mplnResults$allResults,
                          gmin = mplnResults$gmin,
                          gmax = mplnResults$gmax,
                          parallel = "TRUE"))

})
# [END]
