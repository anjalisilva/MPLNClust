# Log likelihood calculation
calcLikelihood <- function(z, PI, dataset, 
                           mu_g, G, Sig_g, 
                           thetaStan, normFactors) { 
  n <- nrow(dataset)
  like <- matrix(NA, nrow = n, ncol = G)
  for (g in 1:G) {
    for (i in 1:n) {
      x <- thetaStan[[g]][i, ]
      d <- ncol(dataset)
      like[i, g] <- (z[i,g] * (log(PI[g]) +
                    t(dataset[i, ]) %*% (x + normFactors) - 
                    sum(exp(x + normFactors)) - 
                    sum(lfactorial(dataset[i, ])) - d / 2 * log(2 * pi) - 
                    1 / 2 * log(det(Sig_g[((g - 1) * d + 1):(g * d), ])) - 
                    0.5 * t(x-mu_g[g, ]) %*% 
                    solve(Sig_g[((g - 1) * d + 1):(g * d), ]) %*% 
                    (x - mu_g[g, ])))
    }
  }
  loglike <- sum(rowSums(like))
  return(loglike)
  # Developed by Anjali Silva
}
