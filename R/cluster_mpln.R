cluster_mpln<-function(y,z,G,n_chain,numb_iterations, initialization, normalizefac, mod){
  
  d<-ncol(y)
  n<-nrow(y)
  
  norm_mu_outer<-norm_sigma_outer<-vector() # for convergence calculation
  median_mu_outer<-median_sigma_outer<-list()
  mu_all_outer<-sigma_all_outer<-list() # for saving mu and sigma values
  obs<-PI<-logL<-vector()
  it_outer<-2 # the starting value of interation for outer loop
  conv_outer<-0
  
  
  if (all(is.na(initialization))==TRUE || all(initialization =="init")){
    mu_all_outer[[1]]<-mu_g <- matrix(log(mean(y)), ncol=d, nrow=G) # mean for both t and normal distribution
    sigma_all_outer[[1]]<-Sig_g <- do.call("rbind", rep(list(cov(log(y+1))*d), G)) # sig for sigma of t distribtuion
  }else{
    mu_all_outer[[1]]<-mu_g <- initialization$finalmu
    sigma_all_outer[[1]]<-Sig_g <- initialization$finalsigma
    z=initialization$probaPost
  }
  
  while(!conv_outer){ 
    for(g in 1:G){
      obs[g]<-sum(z[,g]) # number of observations in each group
      PI[g]<-obs[g]/n  # obtain probability of each group
    }
    
    theta_Stan<-E_theta2<-list()
    rstan_results<-stanrun(model=mod, gmin=1, gmax=G, y=y, mu_all_outer=mu_all_outer, it_outer=it_outer, sigma_all_outer=sigma_all_outer, numb_iterations=numb_iterations, n_chain=n_chain, normalizefacs=normalizefac)
    
    fit = rstan_results$fitrstan
    numb_iterations = rstan_results$numb_iterations
    
    for (g in 1:G){
      tt<-as.matrix(fit[[g]])
      theta_Stan[[g]]<-matrix(NA,nrow=n,ncol=d)
      E_theta2[[g]]<-list()
      
      for (i in 1:n){
        zz<-c(1:(d-1))*n+i
        theta_mat<-tt[,c(i,zz)]
        theta_Stan[[g]][i,]<-colMeans(theta_mat)
        E_theta2[[g]][[i]]<-z[i,g]*t(tt[,c(i,zz)])%*%tt[,c(i,zz)]/((0.5*numb_iterations)*n_chain)
      }
      
      mu_g[g,]<-colSums(z[,g]*theta_Stan[[g]])/sum(z[,g])
      Sig_g[((g-1)*d+1):(g*d),]<-Reduce("+",E_theta2[[g]])/sum(z[,g])-mu_g[g,]%*%t(mu_g[g,])
    }
    
    
    mu_all_outer[[it_outer]]<-mu_g
    sigma_all_outer[[it_outer]]<-Sig_g
    
    logL[it_outer]<-calc_likelihood(z=z, PI=PI, y=y, mu_g=mu_all_outer[[it_outer]], G=G, Sig_g=sigma_all_outer[[it_outer]], theta_Stan=theta_Stan, normalizefactors=normalizefac)
    
    # convergence of outer loop
    norm_mu_outer[it_outer]<-norm((mu_all_outer[[it_outer]]-mu_all_outer[[it_outer-1]]),type="F")
    norm_sigma_outer[it_outer]<-norm(sigma_all_outer[[it_outer]]-sigma_all_outer[[it_outer-1]],type="F")
    median_mu_outer[[it_outer]]<-median(norm_mu_outer, na.rm = TRUE)
    median_sigma_outer[[it_outer]]<-median(norm_sigma_outer, na.rm = TRUE)
    #par(mfrow=c(1,2))
    #plot(norm_mu_outer, main=paste0("Norm outer mean, G=", G), type="l", ylab="median(norm_mu_outer)", xlab="iterations")
    #plot(norm_sigma_outer, main=paste0("Norm outer sigma, G=", G), type="l", ylab="median(norm_sigma_outer)", xlab="iterations")
    
    
    threshold_outer<-2
    if(it_outer>(threshold_outer+1)){
      
      cat("\nMedian difference of mean and sigma in outer loop respectively ", c(abs(median_mu_outer[[it_outer-threshold_outer]]-median_mu_outer[[it_outer]]))) 
      if( ( (abs(median_mu_outer[[it_outer-threshold_outer]]-median_mu_outer[[it_outer]])<5) && (abs(median_sigma_outer[[it_outer-threshold_outer]]-median_sigma_outer[[it_outer]])<5) ) || it_outer>100){
        cat("\nConvergence of mu and sigma at outer loop iteration ", it_outer) # take out absolute value
        programclust<-vector()
        programclust<-map(z)
        
        # checking for spurious clusters and getting rid of them
        #keep<-as.numeric(names(which(table(programclust)>5))) 
        #if ( (length(keep) !=length(unique(programclust))) && (length(keep) !=0) ){
        #  z<-as.matrix(z[,keep])
        #  z<-z/rowSums(z)
        #  programclust<-map(z)
        #}
        
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
    
    # if running for initialization, need to stop after 10 iterations
    if(it_outer==10 && all(is.na(initialization) !=TRUE)){
      if(all(initialization == "init")){
        programclust<-vector()
        programclust<-map(z)
        conv_outer<-1
      }
    }
    
    if(conv_outer!=1){ # only update until convergence, not after
      z<-zvalue_calculation(theta_Stan=theta_Stan,y=y,G=G,mu_g=mu_g,Sig_g=Sig_g,PI=PI, normalizefactors=normalizefac)
      it_outer<-it_outer+1 # updating outer loop iteration
      numb_iterations = numb_iterations+10
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
    probaPost = z,
    stanresults = fit)
  
  class(results) <- "MPLNcluster"
  return(results)
} 