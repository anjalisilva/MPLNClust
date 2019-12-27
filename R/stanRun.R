# Stan sampling
stanRun <- function(model, gmin, gmax, dataset,
                    mu_all_outer, it_outer,
                    sigma_all_outer, numb_iterations,
                    n_chain = n_chain, normalizefacs) {

  fitrstan <- list()
  dimensionality <- ncol(dataset)

  for (g in seq_along(gmin:gmax)) {
    data1 = list(d = ncol(dataset),
    N = nrow(dataset),
    y = dataset,
    mu = mu_all_outer[[it_outer-1]][g, ],
    Sigma = sigma_all_outer[[it_outer-1]][((g - 1) *
            dimensionality + 1):(g*dimensionality), ],
    normfactors = as.vector(normalizefacs))
    stanproceed <- 0
    try = 1

    while (! stanproceed) {

      # cat("\nRstan generating sample at outer iteration",
      # it_outer, "for g: ",g , "try: ", try)
      # cat("\nNumber of iterations is", numb_iterations, "\n")
      fitrstan[[g]] <- rstan::sampling(object = model,
                                       data = data1,
                                       iter = numb_iterations,
                                       chains = n_chain,
                                       verbose = FALSE,
                                       refresh = -1)

      if (all(summary(fitrstan[[g]])$summary[, "Rhat"] < 1.1) == TRUE &&
          all(summary(fitrstan[[g]])$summary[, "n_eff"] > 100) == TRUE) {
          stanproceed <- 1
      } else if(all(summary(fitrstan[[g]])$summary[,"Rhat"] < 1.1) != TRUE ||           all(summary(fitrstan[[g]])$summary[, "n_eff"] > 100) != TRUE) {
          if(try == 10) { # stop after 10 attempts
            stanproceed = 1
          }
        numb_iterations = numb_iterations + 100
        try = try + 1
      }
    }
  }

  results <- list(fitrstan = fitrstan,
                  numb_iterations = numb_iterations)
  class(results) <- "RStan"
  return(results)
  # Developed by Anjali Silva
}
