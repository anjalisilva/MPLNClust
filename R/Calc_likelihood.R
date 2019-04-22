# Log likelihood calculation
calc_likelihood <- function(z, PI, y, mu_g, G, Sig_g, theta_Stan, normalizefactors){ 
  n<-nrow(y)
  like<-matrix(NA, nrow=n, ncol=G)
  for (g in 1:G){
    for (i in 1:n){
      x<-theta_Stan[[g]][i,]
      d<-ncol(y)
      like[i,g]<-(z[i,g] *(log(PI[g]) +
                             t(y[i,])%*%(x+normalizefactors)-sum(exp(x+normalizefactors))-sum(lfactorial(y[i,]))-
                             d/2*log(2*pi)-1/2*log(det(Sig_g[((g-1)*d+1):(g*d),]))-0.5*t(x-mu_g[g,])%*%solve(Sig_g[((g-1)*d+1):(g*d),])%*%(x-mu_g[g,])))
    }
  }
  loglike<-sum(rowSums(like))
  return(loglike)
}
