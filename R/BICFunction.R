
# BIC function 
BIC_function = function(ll, k, n, run, gmin, gmax){
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
