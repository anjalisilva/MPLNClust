# AIC3 function
AIC3_function = function(ll, k, run, gmin, gmax, dataset){
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
