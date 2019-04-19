# initialization function
initializationrun<-function(gmodel, dataset, init_method, n_init_iterations, n_chains, n_iterations, initialization=NA, normalizefactors, mod){
  z<-init_runs<-list()
  logL_init<-vector()
  n<-nrow(dataset)
  d<-ncol(dataset)
  
  for(iterations in 1:n_init_iterations){
    if (init_method=="kmeans" | is.na(init_method)){
      if (!require(mclust)) suppressWarnings(install.packages('mclust')) # loading needed packages
      suppressWarnings(library(mclust))
      z[[iterations]]<-unmap(kmeans(log(dataset+1/3),gmodel)$cluster)
    }else if (init_method=="random"){
      if(gmodel==1){ # generating z if g=1
        z[[iterations]] <- as.matrix(rep.int(1, times=n), ncol=gmodel, nrow=n)
      } else { # generating z if g>1
        z_conv=0
        while(!z_conv){ # ensure that dimension of z is same as G (i.e. 
          # if one column contains all 0s, then generate z again)
          z[[iterations]] <- t(rmultinom(n, size = 1, prob=rep(1/gmodel,gmodel))) 
          if(length(which(colSums(z[[iterations]])>0)) ==gmodel){
            z_conv=1
          }
        }
      }
    }else if (init_method=="medoids"){
      if (!require(cluster)) install.packages('cluster') 
      library(cluster)
      
      if (!require(mclust)) suppressWarnings(install.packages('mclust')) # loading needed packages
      suppressWarnings(library(mclust))
      
      z[[iterations]]<-unmap(pam(log(dataset+1/3),k=gmodel)$cluster)
    }else if (init_method=="clara"){
      if (!require(cluster)) install.packages('cluster') 
      library(cluster)
      
      z[[iterations]]<-unmap(clara(log(dataset+1/3),k=gmodel)$cluster)
    }else if (init_method=="fanny"){
      if (!require(cluster)) install.packages('cluster') 
      library(cluster)
      
      z[[iterations]]<-unmap(fanny(log(dataset+1/3),k=gmodel)$cluster)
    }
    
    init_runs[[iterations]]=cluster_mpln(dataset=dataset,z=z[[iterations]],G=gmodel,n_chains=n_chains,n_iterations=n_iterations, initialization="init", normalizefac=normalizefactors, mod=mod)
    logL_init[iterations] <- unlist(tail((init_runs[[iterations]]$loglikelihood), n=1)) 
  }
  
  initialization<-init_runs[[which(logL_init==max(logL_init, na.rm = TRUE))[1]]]
  return(initialization)
}  
