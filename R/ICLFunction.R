
# ICL function
ICL_function = function(bIc, gmax, gmin, run){
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
