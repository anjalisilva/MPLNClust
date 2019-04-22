
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

