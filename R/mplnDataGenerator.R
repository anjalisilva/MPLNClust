# Generating data using mixtures of MPLN
mplnDataGenerator <- function(nObservations,
                              dimensionality,
                              mixingProportions,
                              mu, sigma,
                              produceImage) {
  # Generate
  z <- t(stats::rmultinom(nObservations, size = 1, mixingProportions))

  y <- theta <- n_g <- vector("list", length = length(mixingProportions))
  # for visualization only
  theta2 <- matrix(NA, ncol = dimensionality, nrow = nObservations)

  for (i in seq_along(1:length(mixingProportions))) {
    n_g[[i]] <- which(z[, i] == 1)
    theta[[i]] <- mvtnorm::rmvnorm(n = length(n_g[[i]]),
                                   mean = mu[i, ],
                                   sigma = sigma[((i - 1) *
                                   dimensionality + 1):(i * dimensionality), ])
    theta2[n_g[[i]], ] <- mvtnorm::rmvnorm(n = length(n_g[[i]]),
                                  mean = mu[i, ],
                                  sigma = sigma[((i  -1) *
                                  dimensionality + 1):(i * dimensionality), ])
  }

  y <- matrix(NA, ncol = dimensionality, nrow = nObservations)
  for (i in seq_along(1:nObservations)) {
    for (j in seq_along(1:dimensionality)) {
      y[i, j] <- stats::rpois(1, exp(theta2[i, j]))
    }
  }

  norms <- log(edgeR::calcNormFactors(y))

  #generating counts with norm factors
  y2 <- matrix(NA, ncol = dimensionality, nrow = nObservations)
  for (i in seq_along(1:nObservations)) {
    for (j in seq_along(1:dimensionality)) {
      y2[i, j] <- stats::rpois(n = 1, lambda = exp(theta2[i, j] + norms[j]))
    }
  }

  if (produceImage == "Yes") {
    # Obtaining path to save images
    pathNow <- getwd()
    png(paste0(pathNow, "/PairsPlot.png"))
    pairs(log(y2), col = map(z) + 1,
          main = "Pairs plot of log-transformed data")
    dev.off()
  }

  results<-list(dataset = y2,
                truemembership = map(z),
                truenormfactors = norms,
                observations = nObservations,
                dimensionality = dimensionality,
                mixingProportions = mixingProportions,
                mu = mu,
                sigma = sigma)

  class(results) <- "mplnDataGenerator"
  return(results)
  # Developed by Anjali Silva
}



