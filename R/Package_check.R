# For CheckPackageOnly Function

# # Function: check whether a package is available in the system, if not prints an error.

# Input:
# RegularPckgs: list of packages to be downloaded
# BioManagerPckgs: list of packages that needs to be downloaded using biocLite

# Output: Prints if the package is loaded or not. 


RegularPckgs <- c("coda",
                 "capushe",
                 "cluster",
                 "clusterGeneration",
                 "mclust",
                 "mvtnorm",
                 "Rcpp",
                 "rstan",
                 "parallel")

BioManagerPckgs<- c("edgeR")

CheckPackageOnly <- function(RegularPckgs=NA, BioManagerPckgs=NA){
  fail = FALSE
  # Code developed by A.Silva based on code provided by Dr. Marcelo Ponce, April 2019
  if (!is.na((RegularPckgs)[1])){
    for (pckg in RegularPckgs) {
      cat("Checking package", pckg, "\n")
      # check whether the package is NOT loaded
      if (! paste('package:',pckg,sep="") %in% search()) {
        # check whether the package is available in the system
        if (pckg %in%  .packages(all.available = TRUE)) {
          # load the package
          cat("Loading library",pckg,"... \n")
          library(pckg, verbose=FALSE, character.only=TRUE)
        } else {
          cat("**************** Package: ",pckg, "is not found. Install using",paste("install.packages('",pckg, "') **************** \n",sep="")) 
          fail =TRUE
        }
      }
    }
  }
  
  if (!is.na((BioManagerPckgs)[1])){
    for (pckg in BioManagerPckgs) {
      cat("Checking package", pckg, "\n")
      # check whether the package is NOT loaded
      if (! paste('package:',pckg,sep="") %in% search()) {
        # check whether the package is available in the system
        if (pckg %in%  .packages(all.available = TRUE)) {
          # load the package
          cat("Loading library",pckg,"... \n")
          library(pckg, character.only=TRUE, verbose=FALSE)
        } else {
          cat("**************** Package: ",pckg, "is not found. Install using",paste("BiocManager::install('",pckg,"') ****************\n",sep="")) 
          fail =TRUE
        }
      }
    }
  }
  
  if (fail) stop()
}

# Function to check the versions of a given set of packages
CheckVersion <- function(pckges) {
  
  for (pck in pckges){
    cat(pck, as.character(packageVersion(pck)), '\n')
  }
  
}

CheckPackageOnly(RegularPckgs = RegularPckgs, BioManagerPckgs=BioManagerPckgs)

CheckVersion(c(RegularPckgs,BioManagerPckgs))
