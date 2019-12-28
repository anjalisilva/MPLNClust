# Generating simulated data

trueMu1 <- c(6.5, 6, 6, 6, 6, 6)
trueMu2 <- c(2, 2.5, 2, 2, 2, 2)

trueSigma1 <- diag(6) * 2
trueSigma2 <- diag(6)

simulated_counts <- mplnDataGenerator(nObservations = 50,
                                      dimensionality = 6,
                                      mixingProportions = c(0.79, 0.21),
                                      mu = rbind(trueMu1, trueMu2),
                                      sigma = rbind(trueSigma1, trueSigma2),
                                      produceImage = "No")

MPLNClustResults <- mpln(dataset = simulated_counts$dataset,
                         membership = simulated_counts$truemembership,
                         gmin = 1,
                         gmax = 1,
                         nChains = 3,
                         nIterations = 300,
                         initMethod = "kmeans",
                         nInitIterations = 0,
                         normalize = "Yes")
