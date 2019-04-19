# calling the clustering function
calling_clustering = function(dataset, Gmin, Gmax, n_chains, n_iterations=NA, init_method=NA, n_init_iterations=NA, norm_factors, mod){
  
  ptm_inner = proc.time() 
  
  for (gmodel in 1:(Gmax-Gmin+1)){
    
    if(length(1:(Gmax-Gmin+1)) == Gmax){
      clustersize = gmodel
    }else if(length(1:(Gmax-Gmin+1)) < Gmax){
      clustersize = seq(Gmin, Gmax, 1)[gmodel]
    }
    
    if(n_init_iterations!=0){
      initializeruns=initializationrun(gmodel=clustersize, dataset=dataset, init_method=init_method, n_init_iterations=n_init_iterations, n_chains=n_chains, n_iterations=n_iterations, initialization=NA, normalizefactors=norm_factors, mod=mod)
      allruns=cluster_mpln(dataset=dataset,z=NA,G=clustersize,n_chains=n_chains,n_iterations=n_iterations, initialization=initializeruns, normalizefac=norm_factors, mod=mod)
    }else if(n_init_iterations == 0){
      allruns=cluster_mpln(dataset=dataset, z=unmap(kmeans(log(dataset+1/3),clustersize)$cluster), G=clustersize, n_chains=n_chains, n_iterations=n_iterations, initialization=NA, normalizefac=norm_factors, mod=mod)
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
}

