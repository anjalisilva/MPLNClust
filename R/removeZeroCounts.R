# Calculate and remove rows with zeros
removeZeroCounts <- function(dataset, output = FALSE) {

  zeroSUMrows <- which(rowSums(dataset) == 0)
  nbrZeroSUMrows <- length(zeroSUMrows)

  if (output) {
    write.csv(rownames(dataset[zeroSUMrows, ]), 
              file=paste0("zeroSUMrows_", 
              format(Sys.time(), "%d%b%Y"),".csv") )
  }

  if (nbrZeroSUMrows > 0) {
      dataset <- dataset[- zeroSUMrows, ]
      cat(paste(nbrZeroSUMrows,
        "row(s) removed from the dataset because sum = 0",'\n'))
  }

 
  RESULTS <- list(dataset = dataset)
  class(RESULTS) <- "MPLN_ZerosRemoved"
  return(RESULTS)
  # Developed by Anjali Silva
}
