
# AIC calculation
AIC_function <- function(ll, k, run, gmin, gmax){
  AIC <- -2*ll+ 2*k
  AICmodel<-seq(gmin, gmax, 1)[grep(min(AIC,na.rm = TRUE), AIC)]
  AICmodel_labels<-run[[grep(min(AIC,na.rm = TRUE), AIC)]]$allresults$clusterlabels
  AICMessage<-NA
  
  if (max(AICmodel_labels)!=AICmodel){
    AICmodel<-max(AICmodel_labels)
    AICMessage<-"Spurious or empty cluster resulted."
  }
  
  AICresults<-list(allAICvalues=AIC,
                   AICmodelselected=AICmodel,
                   AICmodelselected_labels=AICmodel_labels,
                   AICMessage=AICMessage)
  class(AICresults) <- "AIC"
  return(AICresults)
}  


# AIC3 calculation
AIC3_function <- function(ll, k, run, gmin, gmax){
  AIC3 <- -2*ll+ 3*k
  AIC3model<-seq(gmin, gmax, 1)[grep(min(AIC3,na.rm = TRUE), AIC3)]
  AIC3model_labels<-run[[grep(min(AIC3,na.rm = TRUE), AIC3)]]$allresults$clusterlabels
  AIC3Message<-NA
  
  if (max(AIC3model_labels)!=AIC3model){
    AIC3model<-max(AIC3model_labels)
    AIC3Message<-"Spurious or empty cluster resulted."
  }
  AIC3results<-list(allAIC3values=AIC3,
                    AIC3modelselected=AIC3model,
                    AIC3modelselected_labels=AIC3model_labels,
                    AIC3Message=AIC3Message)
  class(AIC3results) <- "AIC3"
  return(AIC3results)
}


# BIC calculation 
BIC_function <- function(ll, k, n, run, gmin, gmax){
  BIC <- -2*ll+ (k* log(n))
  BICmodel<-seq(gmin, gmax, 1)[grep(min(BIC,na.rm = TRUE), BIC)]
  BICmodel_labels<-run[[grep(min(BIC,na.rm = TRUE), BIC)]]$allresults$clusterlabels
  BICMessage<-NA
  
  if (max(BICmodel_labels)!=BICmodel){
    BICmodel<-max(BICmodel_labels)
    BICMessage<-"Spurious or empty cluster resulted."
  }
  
  BICresults<-list(allBICvalues=BIC,
                   BICmodelselected=BICmodel,
                   BICmodelselected_labels=BICmodel_labels,
                   BICMessage=BICMessage)
  class(BICresults) <- "BIC"
  return(BICresults)
}

# ICL calculation
ICL_function <- function(bIc, gmax, gmin, run){
  ICL<-vector()
  for (g in 1:(gmax-gmin+1)){
    z<-run[[g]]$allresults$probaPost
    mapz<-mclust::unmap(run[[g]]$allresults$clusterlabels)
    forICL<-function(g){sum(log(z[which(mapz[,g]==1),g]))}
    ICL[g] <- bIc$allBICvalues[g] + sum(sapply(1:ncol(mapz),forICL))
  }
  ICLmodel<-seq(gmin, gmax, 1)[grep(min(ICL, na.rm = TRUE), ICL)]
  ICLmodel_labels<-run[[grep(min(ICL, na.rm = TRUE), ICL)]]$allresults$clusterlabels
  ICLMessage<-NA
  
  if (max(ICLmodel_labels)!=ICLmodel){
    ICLmodel<-max(ICLmodel_labels)
    ICLMessage<-"Spurious or empty cluster resulted."
  }
  
  
  ICLresults<-list(allICLvalues=ICL,
                   ICLmodelselected=ICLmodel,
                   ICLmodelselected_labels=ICLmodel_labels,
                   ICLMessage=ICLMessage)
  class(ICLresults) <- "ICL"
  return(ICLresults)
}


