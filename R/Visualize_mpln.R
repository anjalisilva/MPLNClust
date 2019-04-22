# Visualize clustered results
visualize_mpln<-function(dataset, ClusterMembershipVector){
  
  # Checking/ Loading needed packages
  LoadCheckPkg(pckgs=c("pheatmap","gplots","RColorBrewer","MASS"))
  
  # Obtaining path to save images
  pathNow<-getwd()
  
  # Saving cluster membership for each observation
  DataPlusLabs=cbind(dataset,ClusterMembershipVector)
  ordervector = list()
  anothervector = list()
  
  for (i in 1:max(ClusterMembershipVector)){
    ordervector[[i]]=which(DataPlusLabs[,ncol(dataset)+1]==i)
    anothervector[[i]]=rep(i,length(which(DataPlusLabs[,ncol(dataset)+1]==i)))
  }
  
  vec<-unlist(ordervector)
  colorsvector<-unlist(anothervector)
  
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  
  # Heatmap 1
  png(paste0(pathNow,"/Clustering_Heatmap1.png"))
  heatmap.2(as.matrix(dataset[vec,]),dendrogram="column",trace="none",scale="row",
    Rowv=FALSE,  col = rev(redgreen(75)), RowSideColors=col_vector[colorsvector+1])
  par(xpd=TRUE)
  legend("topleft",      
    legend = paste0("Cluster ", unique(colorsvector)),
    col = unique(col_vector[colorsvector+1]), 
    lty= 1,             
    lwd = 5,           
    cex=.5,  xpd = TRUE, horiz = FALSE
  )
  dev.off()
  
  # Heatmap 2
  annotation_row = data.frame(
    Cluster = factor(ClusterMembershipVector[vec]))
  if(is.null(rownames(dataset)) == TRUE){
    rownames(dataset)  = paste("Gene",c(1:nrow(dataset[vec,])))
    rownames(annotation_row) = rownames(dataset[vec,])
  }else{
    rownames(annotation_row) = rownames(dataset[vec,])
  }
  
  png(paste0(pathNow,"/Clustering_Heatmap2.png"))
  pheatmap(as.matrix(dataset[vec,]), show_colnames = T, labels_col=colnames(dataset), annotation_row =annotation_row , fontface="italic", legend = T, scale ="row",border_color = "black", cluster_row = FALSE, cluster_col = FALSE, color =  rev(redgreen(1000)) )
  dev.off()
  
  # Line Plots
  png(paste0(pathNow,"/Clustering_LinePlots.png"))
  par(mfrow=c(2,max(ClusterMembershipVector)))
  for(cluster in 1:max(ClusterMembershipVector)){
      # Save how many observations below to each cluster size, given by 'cluster'
      toplot_1=DataPlusLabs[which(DataPlusLabs[,ncol(dataset)+1]==cluster),c(1:ncol(dataset))]
      # Save column mean in last row
      toplot1=rbind(log(toplot_1+1), colMeans(log(toplot_1+1)))
      # If discontinunity is needed between samples (e.g. for 6 samples)
      # toplot1_space=cbind(toplot1[,c(1:3)],rep(NA,nrow(toplot_1)+1),toplot1[,c(4:6)])
      matplot(t(toplot1), type="l", pch=1, col=c(rep(1,nrow(toplot_1)),7), xlab="Samples", ylab="Expression (log counts)", cex=1, lty=c(rep(2,nrow(toplot_1)),1),lwd=c(rep(1,nrow(toplot_1)),3), xaxt="n", xlim=c(1,ncol(toplot1)), main=paste("Cluster ",cluster))
      axis(1,at = c(1:ncol(dataset)), labels=colnames(dataset))
    }
  dev.off()
  
}
