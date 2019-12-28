test_that("clustering", {
  library(MPLNClust)

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

  MPLNClustResults <- mpln(dataset = simulated_counts$dataset,
                           membership = simulated_counts$trueMembership,
                           gmin = 1,
                           gmax = 2,
                           nChains = 3,
                           nIterations = 600,
                           initMethod = "kmeans",
                           nInitIterations = 0,
                           normalize = "Yes")

  expect_that(length(MPLNClustResults), equals(16))
  expect_that(MPLNClustResults, is_a("MPLN"))
  expect_that(trunc(MPLNClustResults$ICL_all$ICLmodelselected), equals(2))
  expect_that(trunc(MPLNClustResults$AIC_all$AICmodelselected), equals(2))
  expect_that(trunc(MPLNClustResults$BIC_all$BICmodelselected), equals(2))
})
