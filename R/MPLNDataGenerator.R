
#### function ####
Datagenerator<-function(i, N, d, pi_g, means, sigmas, ProduceImage){
  
  source("PackageCheck.R")
  
  set.seed(i)
  z<-t(rmultinom(N,size=1,pi_g))
  
  # loading needed packages
  LoadCheckPkg(pckgs=c("mvtnorm","clusterGeneration","edgeR","mclust"))
  
  y<-theta<-n_g <- vector("list", length = length(pi_g)) 
  theta2<-matrix(NA,ncol=d,nrow=N) # for visualization only
  
  for (i in 1:length(pi_g)){
    n_g[[i]]<-which(z[,i]==1)
    theta[[i]]<-rmvnorm(length(n_g[[i]]), mean=means[i,], sigma=sigmas[((i-1)*d+1):(i*d),])
    theta2[n_g[[i]],]<-rmvnorm(length(n_g[[i]]), mean=means[i,], sigma=sigmas[((i-1)*d+1):(i*d),])
  }
  
  y<-matrix(NA,ncol=d,nrow=N)
  for (i in 1:N){
    for (j in 1:d){
      y[i,j]<-rpois(1,exp(theta2[i,j])) 
    }
  }
  
  norms <- log(calcNormFactors(y))

  #generating counts with norm factors
  y2<-matrix(NA,ncol=d,nrow=N)
  for (i in 1:N){
    for (j in 1:d){
      y2[i,j]<-rpois(1,exp(theta2[i,j]+norms[j])) 
    }
  }

  if (ProduceImage=="Yes"){
    # Obtaining path to save images
    pathNow<-getwd()
    
    png(paste0(pathNow,"/PairsPlot.png"))
    pairs(log(y2), col=map(z)+1, main="Pairs plot of log-transformed data")
    dev.off() 
  }
  
  results<-list(dataset=y2,
                truemembership=map(z),
                truenormfactors =norms,
                observations = N,
                dimensionality = d,
                pi_g = pi_g,
                means = means,
                sigmas = sigmas)
  
  class(results) <- "MPLNDataGenerator"
  return(results)
}  



