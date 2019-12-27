# Initialization
initializationRun <- function(gmodel, dataset, init_method,
                              init_iterations, n_chain,
                              numb_iterations,
                              initialization = NA,
                              normalizefactors, mod) {

  z <- init_runs <- list()
  logL_init <- vector()
  nObservations <- nrow(dataset)
  dimensionality <- ncol(dataset)

  for(iterations in seq_along(1:init_iterations)) {
    # setting seed, to ensure if multiple iterations are selected by
    # user, then each run will give a different result.
    set.seed(iterations)
    if (init_method == "kmeans" | is.na(init_method)) {
      z[[iterations]] <- unmap(stats::kmeans(log(dataset + 1/3),
                               gmodel)$cluster)
    } else if (init_method == "random") {
      if(gmodel == 1) { # generating z if g=1
        z[[iterations]] <- as.matrix(rep.int(1, times = nObservations),
                                     ncol = gmodel,
                                     nrow = nObservations)
      } else { # generating z if g>1
        z_conv = 0
        while(! z_conv) {
          # ensure that dimension of z is same as G (i.e.,
          # if one column contains all 0s, then generate z again)
          z[[iterations]] <- t(stats::rmultinom(nObservations, size = 1,
                               prob = rep(1 / gmodel, gmodel)))
          if(length(which(colSums(z[[iterations]]) > 0)) == gmodel) {
            z_conv = 1
          }
        }
      }
    }else if (init_method == "medoids") {
      z[[iterations]] <- unmap(cluster::pam(log(dataset + 1/3),
                               k = gmodel)$cluster)
    }else if (init_method == "clara") {
      z[[iterations]] <- unmap(cluster::clara(log(dataset + 1/3),
                               k = gmodel)$cluster)
    }else if (init_method == "fanny") {
      z[[iterations]] <- unmap(cluster::fanny(log(dataset + 1/3),
                               k = gmodel)$cluster)
    }

    init_runs[[iterations]] <- mplnCluster(dataset = dataset,
                                            z = z[[iterations]],
                                            G = gmodel,
                                            nChains = n_chain,
                                            nIterations = numb_iterations,
                                            initialization = "init",
                                            normalizefac = normalizefactors,
                                            mod = mod)
    logL_init[iterations] <-
      unlist(tail((init_runs[[iterations]]$loglikelihood), n = 1))
  }

  initialization <- init_runs[[which(logL_init ==
                               max(logL_init, na.rm = TRUE))[1]]]
  return(initialization)
  # Developed by Anjali Silva
}
