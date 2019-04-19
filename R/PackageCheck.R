
# A function to check whether a package is available in the system
LoadCheckPkg <- function(pckgs){
  
  fail = FALSE
  for (pckg in pckgs) {
    # check whether the package is NOT loaded
    if (! paste('package:',pckg,sep="") %in% search()) {
      # check whether the package is available in the system
      if (pckg %in%  .packages(all.available = TRUE)) {
        # load the package
        cat("Loading library",pckg,"... \n")
        library(pckg, character.only=TRUE)
      } else {
        msg <- paste("Package:",pckg, "not found! You need to install this package using",paste("install.package('",pckg,"') \n",sep=""))
        cat(msg)
        
        fail =TRUE
      }
    }
  }
  
  if (fail) stop()
}