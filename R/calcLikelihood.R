# Log likelihood calculation
calcLikelihood <- function(z, PI, dataset, 
                           mu_g, G, Sig_g, 
                           thetaStan, normFactors) { 
  nObservations <- nrow(dataset)
  like <- matrix(NA, nrow = nObservations, ncol = G)
  for (g in seq_along(1:G)) {
    for (i in seq_along(1:nObservations)) {
      x <- thetaStan[[g]][i, ]
      dimensionality <- ncol(dataset)
      like[i, g] <- (z[i, g] * (log(PI[g]) +
                    t(dataset[i, ]) %*% (x + normFactors) - 
                    sum(exp(x + normFactors)) - 
                    sum(lfactorial(dataset[i, ])) - 
                    dimensionality / 2 * log(2 * pi) - 
                    1 / 2 * log(det(Sig_g[((g - 1) * 
                    dimensionality + 1):(g * dimensionality), ])) - 
                    0.5 * t(x-mu_g[g, ]) %*% 
                    solve(Sig_g[((g - 1) * 
                    dimensionality + 1):(g * dimensionality), ]) %*% 
                    (x - mu_g[g, ])))
    }
  }
  loglike <- sum(rowSums(like))
  return(loglike)
  # Developed by Anjali Silva
}
