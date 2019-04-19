# z value calculation function
zvalue_calculation<-function(theta_Stan,dataset,G,mu_g,Sig_g,PI, normalizefactors){
  
  d<-ncol(dataset)
  n<-nrow(dataset)
  forz=sapply(c(1:G), function(g) sapply(c(1:n), function(i) PI[g]*exp(t(dataset[i,])%*%(theta_Stan[[g]][i,]+normalizefactors)-sum(exp(theta_Stan[[g]][i,]+normalizefactors))-sum(lfactorial(dataset[i,]))-
      d/2*log(2*pi)-1/2*log(det(Sig_g[((g-1)*d+1):(g*d),]))-0.5*t(theta_Stan[[g]][i,]-mu_g[g,])%*%solve(Sig_g[((g-1)*d+1):(g*d),])%*%(theta_Stan[[g]][i,]-mu_g[g,])) ) )
  if (G==1){
    errorpossible<-Reduce(intersect, list(which(forz==0),which(rowSums(forz)==0)))
    zvalue<-forz/rowSums(forz)
    zvalue[errorpossible,]<-1
  }else {zvalue<-forz/rowSums(forz)}
  
  return(zvalue)
}

