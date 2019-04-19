# Generate parallel coordinates plots
data_visualization<-function(theta_Stan,dataset,G,mu_g,Sig_g,PI, normalizefactors){
  
  # loading needed packages
  LoadCheckPkg(pckgs=c("MASS","RColorBrewer"))
  
  #library(RColorBrewer)
  #display.brewer.pal(n = 12, name = 'Paired')
  #abc<-brewer.pal(12, "Paired")
  #abc<-brewer.pal(12, "Set3")
  
  #opar<-par(no.readonly=T)
  #par(mfrow = c(1,2))
  #par(oma=c(1.2,2,2,0))
  #par(mar=c(2,2,0.1,0.1))
  ## For each color (cluster) in the random network
  for (i in 1:2){
    input = log(nointernal_sim_2clusters[[1]]$dataset+1)[which(nointernal_sim_2clusters[[1]]$BIC.all$BICmodelselected_labels==i),]
    color = rep(abc[i], ncol(input))
    colnames(input)<-c(1:6)
    parcoord(input, lty = 1, var.label = FALSE, col = color)
    axis(2,at=seq(0,1,length.out=5),labels=round(seq(min(input),max(input), length.out=5),0))
  }
  
}

