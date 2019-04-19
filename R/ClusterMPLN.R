# clustering function
cluster_mpln<-function(dataset,z,G,n_chains,n_iterations, initialization, normalizefac, mod){
  
  d<-ncol(dataset)
  n<-nrow(dataset)
  
  norm_mu_outer<-norm_sigma_outer<-vector() # for convergence calculation
  median_mu_outer<-median_sigma_outer<-list()
  mu_all_outer<-sigma_all_outer<-list() # for saving mu and sigma values
  obs<-PI<-logL<-vector()
  it_outer<-2 # the starting value of interation for outer loop
  conv_outer<-0
  
  
  if (all(is.na(initialization))==TRUE || all(initialization =="init")){
    mu_all_outer[[1]]<-mu_g <- matrix(log(mean(dataset)), ncol=d, nrow=G) # mean for both t and normal distribution
    sigma_all_outer[[1]]<-Sig_g <- do.call("rbind", rep(list(cov(log(dataset+1))*d), G)) # sig for sigma of t distribtuion
  }else{
    mu_all_outer[[1]]<-mu_g <- initialization$finalmu
    sigma_all_outer[[1]]<-Sig_g <- initialization$finalsigma
    z=initialization$probaPost
  }
  
  while(!conv_outer){ 
    cat("************** Running for G =",G,"and Iteration =",it_outer,"******************")
    obs=apply(z, 2, sum) # number of observations in each group
    PI=sapply(obs, function(x) x/n)  # obtain probability of each group
    
    theta_Stan<-E_theta2<-list()
    rstan_results<-stanrun(model=mod, gmin=1,gmax=G,dataset=dataset,mu_all_outer=mu_all_outer, it_outer=it_outer, sigma_all_outer=sigma_all_outer, n_iterations=n_iterations, n_chains=n_chains, normalizefacs=normalizefac)
    
    # turn results into a matrix
    tt=lapply(as.list(c(1:G)), function(x) as.matrix(rstan_results$fitrstan[[x]]) ) 
    
    for (g in 1:G){  
      # expected value of theta
      theta_Stan[[g]]=t(sapply(c(1:n), function(i) colMeans(tt[[g]][,c(i,(t(sapply(c(1:n), function(x) c(1:(d-1))*n+x)))[i,])]) )) 
      # expected value of theta theta-transform
      E_theta2[[g]]=lapply(as.list(c(1:n)), function(i) z[i,g]*t(tt[[g]][,c(i,(t(sapply(c(1:n), function(x) c(1:(d-1))*n+x)))[i,])])%*%tt[[g]][,c(i,(t(sapply(c(1:n), function(x) c(1:(d-1))*n+x)))[i,])]/((0.5*rstan_results$n_iterations)*n_chains))
    }
    
    # updating value of mu
    mu_all_outer[[it_outer]]= mu_g=t(sapply(c(1:G), function(g) colSums(z[,g]*theta_Stan[[g]])/sum(z[,g])))
    
    # updating value of sigma
    sigma_all_outer[[it_outer]]=Sig_g=do.call(rbind, lapply(as.list(c(1:G)), function(g) Reduce("+",E_theta2[[g]])/sum(z[,g])-mu_g[g,]%*%t(mu_g[g,])))
    
    # updating log-likelihood
    logL[it_outer]<-calc_likelihood(z=z, PI=PI, dataset=dataset, mu_g=mu_all_outer[[it_outer]], G=G, Sig_g=sigma_all_outer[[it_outer]], theta_Stan=theta_Stan, normalizefactors=normalizefac)
    
    # convergence of outer loop
    threshold_outer<-10
    if(it_outer>(threshold_outer)){
      
      if (all(heidel.diag(logL[-1], eps=0.1, pvalue=0.05)[,c(1,4)]==1) || it_outer>100){
        
        programclust<-vector()
        programclust<-map(z)
        
        # checking for empty clusters
        J <- 1:ncol(z)
        K <- as.logical(match(J, sort(unique(programclust)), nomatch = 0))
        if(length(J[!K])>0){ # J[!K] tells which are empty clusters
          z<-z[,-J[!K]]
          programclust<-map(z)
        }
        
        conv_outer<-1
      } 
    }
    
    # if running for initialization, need to stop after 1 iteration
    if(it_outer==2 && all(is.na(initialization) !=TRUE)){
      if(all(initialization == "init")){
        programclust<-vector()
        programclust<-map(z)
        conv_outer<-1
      }
    }
    
    # only update until convergence, not after
    if(conv_outer!=1){
      z<-zvalue_calculation(theta_Stan=theta_Stan,dataset=dataset,G=G,mu_g=mu_g,Sig_g=Sig_g,PI=PI, normalizefactors=normalizefac)
      it_outer<-it_outer+1 # updating outer loop iteration
      n_iterations = n_iterations+10
    }
    
  } # end of outer loop
  
  
  results <- list(finalmu=mu_all_outer[[it_outer]]+ matrix(rep(normalizefac,nrow(mu_all_outer[[it_outer]])),byrow=TRUE,ncol=ncol(mu_all_outer[[it_outer]])), 
    finalsigma=sigma_all_outer[[it_outer]],
    allmu = lapply(mu_all_outer, function(x) (x+matrix(rep(normalizefac,nrow(mu_all_outer[[it_outer]])),byrow=TRUE,ncol=ncol(mu_all_outer[[it_outer]])))),      
    allsigma = sigma_all_outer, 
    clusterlabels = programclust,
    iterations = it_outer, 
    proportion = PI, 
    loglikelihood = logL,
    probaPost = z)
  
  class(results) <- "MPLNcluster"
  return(results)
}
