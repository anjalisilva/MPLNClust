# Checking/Loading needed packages
LoadCheckPkg(pckgs=c("parallel","rstan","Rcpp","mclust","mvtnorm","edgeR","capushe","clusterGeneration","coda"))

# function to load auxiliary functions
source("Calc_likelihood.R")
source("Calculate_parameters.R")
source("Calling_clustering.R")
source("Cluster_mpln.R")
source("Initialization_run.R")
source("Main_mpln.R")
source("Model_selection.R")
source("MPLNdata_generator.R")
source("Remove_Zero_Counts.R")
source("Stan_run.R")
source("Visualize_mpln.R")
source("Zvalue_calculation.R")
