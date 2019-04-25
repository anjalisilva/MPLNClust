# Calculate and remove rows with zeros
remove_zero_counts<-function(dataset){
  
  if(length(which(rowSums(dataset)==0))>0){
    dataset <- dataset[-which(rowSums(dataset)==0),]
  }
  
  RESULTS <- list(dataset= dataset)
  class(RESULTS) <- "MPLN_ZerosRemoved"
  return(RESULTS)
                  
}
