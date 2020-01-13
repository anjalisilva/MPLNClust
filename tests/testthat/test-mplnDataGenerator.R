context("Checking for data simulation")
library(MPLNClust)

test_that("Data generation is as expected", {

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


context("Checking for invalid user input")
test_that("Data generate error upon invalid user input", {

  # nObservations provided as character
  expect_error(mplnDataGenerator(nObservations = "50",
                                        dimensionality = 6,
                                        mixingProportions = c(0.79, 0.21),
                                        mu = rbind(trueMu1, trueMu2),
                                        sigma = rbind(trueSigma1, trueSigma2),
                                        produceImage = "No"))


  # Generating simulated data - mu has incorrect dimension
  trueMu1 <- c(6.5, 6, 6, 6, 6)
  trueMu2 <- c(2, 2.5, 2, 2, 2)

  expect_error(mplnDataGenerator(nObservations = 50,
    dimensionality = 6,
    mixingProportions = c(0.79, 0.21),
    mu = rbind(trueMu1, trueMu2),
    sigma = rbind(trueSigma1, trueSigma2),
    produceImage = "No"))


  # Generating simulated data - sigma has incorrect dimension
  trueMu1 <- c(6.5, 6, 6, 6, 6, 6)
  trueMu2 <- c(2, 2.5, 2, 2, 2, 2)
  trueSigma1 <- diag(5) * 2
  trueSigma2 <- diag(5)

  expect_error(mplnDataGenerator(nObservations = 50,
    dimensionality = 6,
    mixingProportions = c(0.79, 0.21),
    mu = rbind(trueMu1, trueMu2),
    sigma = rbind(trueSigma1, trueSigma2),
    produceImage = "No"))


  # mixingProportions does not sum to 1
  expect_error(mplnDataGenerator(nObservations = 50,
    dimensionality = 6,
    mixingProportions = c(0.79, 0.2),
    mu = rbind(trueMu1, trueMu2),
    sigma = rbind(trueSigma1, trueSigma2),
    produceImage = "No"))

})
