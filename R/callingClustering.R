# Calling clustering

callingClustering <- function(data, gmin, gmax, nChains, 
                              nIterations = NA, initMethod=NA, 
                              nInitIterations = NA, 
                              normFactors, model) {
  ptm_inner = proc.time() 
  
  for (gmodel in seq_along(1:(gmax - gmin + 1))) {
    
    if(length(1:(gmax - gmin + 1)) == gmax) {
      clustersize <- gmodel
    } else if(length(1:(gmax - gmin + 1)) < gmax) {
      clustersize <- seq(gmin, gmax, 1)[gmodel]
    }
    
    if(nInitIterations != 0) {
      # cat("\nRunning initialization for G =", clustersize)
      initializeruns <- initializationRun(gmodel = clustersize, 
                                          dataset = data, 
                                          init_method = initMethod, 
                                          init_iterations = nInitIterations,
                                          n_chain = nChains, 
                                          numb_iterations = nIterations, 
                                          initialization = NA, 
                                          normalizefactors = normFactors, 
                                          mod = model)
      # cat("\nInitialization done for G =", clustersize)
      # cat("\nRunning clustering for G =", clustersize)
      allruns <- mplnCluster(dataset = data,
                              z = NA,
                              G = clustersize,
                              nChains = nChains,
                              nIterations = nIterations, 
                              initialization = initializeruns,
                              normalizefac = normFactors, 
                              mod = model)
      # cat("\nClustering done for G =", clustersize)
    } else if(nInitIterations == 0) {
      # cat("\nNo initialization done for G =", clustersize)
      # cat("\nRunning clustering for G =", clustersize)
      allruns <- mplnCluster(dataset = data, 
                              z = unmap(kmeans(log(data + 1/3),
                                               clustersize)$cluster), 
                              G = clustersize, 
                              nChains = nChains,
                              nIterations = nIterations, 
                              initialization = NA, 
                              normalizefac = normFactors, 
                              mod = model)
      # cat("\nClustering done for G =", clustersize)
    }
  }
  
  final_inner <- proc.time() - ptm_inner
  
  RESULTS <- list(gmin = gmin,
                  gmax = gmax,
                  initalization_method = initMethod,
                  all_results = allruns,
                  total_time = final_inner)
  
  class(RESULTS) <- "MPLN"
  return(RESULTS)
  # Developed by Anjali Silva
}
