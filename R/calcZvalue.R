# Zvalue calculation
calcZvalue <- function(theta_Stan, dataset, G, 
                       mu_g, Sig_g, PI, 
                       normalizefactors) {
  d <- ncol(dataset)
  n <- nrow(dataset)
  forz <- matrix(NA,ncol=G,nrow=n)
  
  for (g in 1:G) {
    for (i in 1:n) {
      x <- theta_Stan[[g]][i,]
      # for zig calculation (the numerator part)
      forz[i, g] <- PI[g] * exp(t(dataset[i, ]) %*% (x + normalizefactors) -
                    sum(exp(x + normalizefactors)) - 
                    sum(lfactorial(dataset[i, ])) -
                    d / 2 * log(2 * pi) - 1 / 2 * 
                    log(det(Sig_g[((g - 1) * d + 1):(g*d), ])) - 
                    0.5 * t(x-mu_g[g, ]) %*% 
                    solve(Sig_g[((g - 1) * d + 1):(g * d), ]) %*% 
                    (x-mu_g[g, ]))
    }
    
    # check which forz == 0 and rowSums(forz)==0 and which of these
    # have both equalling to 0 (because 0/0 =NaN)
    if (G == 1) {
      errorpossible <- Reduce(intersect, list(which(forz == 0), 
                       which(rowSums(forz) == 0)))
      zvalue <- forz / rowSums(forz)
      zvalue[errorpossible, ] <- 1
    }else {
      zvalue <- forz / rowSums(forz)
    }
  }
  
  # check which forz == 0 and rowSums(forz)==0 and which of these
  # have both equalling to 0 (because 0/0 =NaN)
  if (G == 1) {
    errorpossible <- Reduce(intersect, 
                            list(which(forz == 0), which(rowSums(forz) == 0)))
    zvalue <- forz / rowSums(forz)
    zvalue[errorpossible, ] <- 1
  }else {
    zvalue <- forz / rowSums(forz)
  }
  return(zvalue)
  # Developed by Anjali Silva
}
