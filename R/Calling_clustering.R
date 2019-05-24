# Calling clustering
calling_clustering <- function(y, Gmin, Gmax, n_chain, numb_iterations=NA, init_method=NA, init_iterations=NA, norm_factors, mod){
  
  ptm_inner = proc.time() 
  
  for (gmodel in 1:(Gmax-Gmin+1)){
    
    if(length(1:(Gmax-Gmin+1)) == Gmax){
      clustersize = gmodel
    }else if(length(1:(Gmax-Gmin+1)) < Gmax){
      clustersize = seq(Gmin, Gmax, 1)[gmodel]
    }
    
    if(init_iterations!=0){
      #cat("\nRunning initialization for G =", clustersize)
      initializeruns=initializationrun(gmodel=clustersize, y=y, init_method=init_method, init_iterations=init_iterations, n_chain=n_chain, numb_iterations=numb_iterations, initialization=NA, normalizefactors=norm_factors, mod=mod)
      #cat("\nInitialization done for G =", clustersize)
      #cat("\nRunning clustering for G =", clustersize)
      allruns=cluster_mpln(y=y,z=NA,G=clustersize,n_chain=n_chain,numb_iterations=numb_iterations, initialization=initializeruns,normalizefac=norm_factors, mod=mod)
      #cat("\nClustering done for G =", clustersize)
    }else if(init_iterations == 0){
      #cat("\nNo initialization done for G =", clustersize)
      #cat("\nRunning clustering for G =", clustersize)
      allruns=cluster_mpln(y=y, z=unmap(kmeans(log(y+1/3),clustersize)$cluster), G=clustersize, n_chain=n_chain, numb_iterations=numb_iterations, initialization=NA, normalizefac=norm_factors, mod=mod)
      #cat("\nClustering done for G =", clustersize)
    }
  }
  
  final_inner<-proc.time()-ptm_inner
  RESULTS <- list(  gmin = Gmin,
                    gmax = Gmax,
                    initalization_method = init_method,
                    allresults = allruns,
                    totaltime = final_inner)
  
  
  class(RESULTS) <- "MPLN"
  return(RESULTS)
  # Developed by Anjali Silva
}
