# stan sampling function
stanrun<-function(model, gmin, gmax, dataset, mu_all_outer, it_outer, sigma_all_outer, n_iterations, n_chains=n_chains, normalizefacs){
  
  fitrstan<-list()
  d<-ncol(dataset)
  
  for (g in gmin:gmax){
    data1=list(d=ncol(dataset),N=nrow(dataset),y=dataset,mu=mu_all_outer[[it_outer-1]][g,],Sigma=sigma_all_outer[[it_outer-1]][((g-1)*d+1):(g*d),], normfactors=as.vector(normalizefacs))
    stanproceed<-0
    try=1
    
    while (!stanproceed){
      cat("\nRstan generating sample at outer iteration", it_outer, "for g: ", g)
      cat("\nNumber of iterations is", n_iterations, "\n")
      fitrstan[[g]]<-sampling(object=model,
        data=data1,
        iter=n_iterations, chains = n_chains, verbose=FALSE, refresh=-1)
      
      if (all(summary(fitrstan[[g]])$summary[,"Rhat"] < 1.1) == TRUE && all(summary(fitrstan[[g]])$summary[,"n_eff"]>100) == TRUE){
        stanproceed<-1
      } else if(all(summary(fitrstan[[g]])$summary[,"Rhat"] < 1.1) != TRUE || all(summary(fitrstan[[g]])$summary[,"n_eff"]>100) != TRUE){
        if(try == 10){ # stop after 10 tries
          stanproceed = 1
        }
        n_iterations = n_iterations+100
        try=try+1
      }
    }
  } # close g loop
  
  
  results <- list(fitrstan = fitrstan,
    n_iterations = n_iterations)
  class(results) <- "RStan"
  return(results)
  
  return(results)
}
