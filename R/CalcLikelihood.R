# likelihood calculation function
calc_likelihood<-function(z, PI, dataset, mu_g, G, Sig_g, theta_Stan, normalizefactors){ 
  
  n<-nrow(dataset)
  like<-matrix(NA, nrow=n, ncol=G)
  d<-ncol(dataset)
  like=sapply(c(1:G), function(g) sapply(c(1:n), function(i) z[i,g] *(log(PI[g]) +
      t(dataset[i,])%*%(theta_Stan[[g]][i,]+normalizefactors)-sum(exp(theta_Stan[[g]][i,]+normalizefactors))-sum(lfactorial(dataset[i,]))-
      d/2*log(2*pi)-1/2*log(det(Sig_g[((g-1)*d+1):(g*d),]))-0.5*t(theta_Stan[[g]][i,]-mu_g[g,])%*%solve(Sig_g[((g-1)*d+1):(g*d),])%*%(theta_Stan[[g]][i,]-mu_g[g,])) ) )
  
  loglike<-sum(rowSums(like))
  return(loglike)
}
