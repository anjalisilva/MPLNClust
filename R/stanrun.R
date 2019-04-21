stanrun<-function(model, gmin, gmax, y, mu_all_outer, it_outer, sigma_all_outer, numb_iterations, n_chain=n_chain, normalizefacs){
  fitrstan<-list()
  d<-ncol(y)
  
  
  for (g in gmin:gmax){
    data1=list(d=ncol(y),N=nrow(y),y=y,mu=mu_all_outer[[it_outer-1]][g,],Sigma=sigma_all_outer[[it_outer-1]][((g-1)*d+1):(g*d),], normfactors=as.vector(normalizefacs))
    stanproceed<-0
    try=1
    
    while (!stanproceed){
      
      cat("\nRstan generating sample at outer iteration", it_outer, "for g: ",g , "try: ", try)
      cat("\nNumber of iterations is", numb_iterations, "\n")
      fitrstan[[g]]<-sampling(object=model,
        data=data1,
        iter=numb_iterations, chains = n_chain, verbose=FALSE, refresh=-1)
      
      if (all(summary(fitrstan[[g]])$summary[,"Rhat"] < 1.1) == TRUE && all(summary(fitrstan[[g]])$summary[,"n_eff"]>100) == TRUE){
        stanproceed<-1
      } else if(all(summary(fitrstan[[g]])$summary[,"Rhat"] < 1.1) != TRUE || all(summary(fitrstan[[g]])$summary[,"n_eff"]>100) != TRUE){
        if(try == 10){ # stop after 10 tries
          stanproceed = 1
        }
        numb_iterations = numb_iterations+100
        try=try+1
      }
    }
  } # close g loop
  
  
  results <- list(fitrstan = fitrstan,
    numb_iterations = numb_iterations)
  class(results) <- "RStan"
  return(results)
  
  return(results)
}
