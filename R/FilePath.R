# A function to establish the location of the script
# Code developed by Dr. Marcelo Ponce, received in April 2019

ALLargs <- function() {
  cmdArgs <- commandArgs(trailingOnly = FALSE)
  scrPATHkw <- "--file="
  match <- grep(scrPATHkw, cmdArgs)
  if (length(match) > 0) {
    # Rscript
    scriptLOC <- dirname(normalizePath(sub(scrPATHkw, "", cmdArgs[match])))
  } else {
    # 'source'd via R console
    scriptLOC <- normalizePath(sys.frames()[[1]]$ofile)
  }
  
  cmdArgs <- commandArgs(trailingOnly =TRUE)
  return(list(scriptLOC,cmdArgs))
}

allArgs <- ALLargs()

print(allArgs)


# Printing #
cat(" MPLNClust v2.0 (2019)",'\n')
cat("----------------------------",'\n')

# Load utilities file with fns. defns.

utilitiesFile <- paste(allArgs[[1]],"/utils_RACS-IGR.R",sep='')
print(utilitiesFile)
#print(dirname(sys.frame(1)$ofile)
#sourceDir <- getSrcDirectory(function(dummy) {dummy})
#print(sourceDir)


if (file.exists(utilitiesFile)) {
  source(utilitiesFile)
} else {
  stop("Critical ERRROR: '",utilitiesFile,"' NOT found!")
}
