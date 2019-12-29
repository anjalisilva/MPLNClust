test_that("data generation", {
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

  expect_that(length(simulated_counts), equals(9))
  expect_that(simulated_counts, is_a("mplnDataGenerator"))
  expect_that(trunc(simulated_counts$observations), equals(50))
  expect_that(trunc(simulated_counts$dimensionality), equals(6))
})
