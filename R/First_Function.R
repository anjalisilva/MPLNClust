
testing_dataset <- simulated_counts$dataset # Assign test dataset using the variable name 

# Making RStan model #
LoadCheckPkg(pckgs=c("parallel","rstan"))
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Making RStan model #
mod = stan_model("MPLN.stan")

# Calculate the number of cores
no_cores = detectCores()-1

# Initiate cluster
cl = makeCluster(no_cores) 

print("Doing clusterExport")
clusterExport(cl,c("mod", "testing_dataset","zvalue_calculation", "calc_likelihood", "stanrun", "initializationrun", "BIC_function","ICL_function","AIC_function","AIC3_function", "calculate_parameters", "cluster_mpln", "calling_clustering"))


print("Doing clusterEvalQ")
#other packages need to be downloaded using clusterEvalQ
clusterEvalQ(cl, library(rstan))
clusterEvalQ(cl, library(Rcpp))
clusterEvalQ(cl, library(mclust))
clusterEvalQ(cl, library(mvtnorm))
clusterEvalQ(cl, library(edgeR))
clusterEvalQ(cl, library(capushe))
clusterEvalQ(cl, library(clusterGeneration))
clusterEvalQ(cl, library(coda))
