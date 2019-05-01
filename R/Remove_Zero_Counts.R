# Calculate and remove rows with zeros
remove_zero_counts<-function(dataset){

  zeroSUMrows <- which(rowSums(dataset)==0)
  nbrZeroSUMrows <- length(zeroSUMrows)
  
  if (nbrZeroSUMrows > 0) {
      dataset <- dataset[-zeroSUMrows,]
      cat(paste(nbrZeroSUMrows,"row(s) removed from the dataset because sum=0",'\n'))
  }
  
  RESULTS <- list(dataset= dataset)
  class(RESULTS) <- "MPLN_ZerosRemoved"
  return(RESULTS)
                  
}
