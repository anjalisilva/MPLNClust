# **Should have all the function R files and MPLN.stan in the same working directory as this file

# Reading needed code
source("AIC_function.R")
source("AIC3_function.R")
source("BIC_function.R")
source("Calc_likelihood.R")
source("Calculate_parameters.R")
source("Calling_clustering.R")
source("Cluster_mpln.R")
source("ICL_function.R")
source("Initialization_run.R")
source("Main_mpln.R")
source("MPLNdata_generator.R")
source("Package_check.R")
source("Stan_run.R")
source("Visualize_mpln.R")
source("Zvalue_calculation.R")


# Values for data simulation
true_mu1 <- c(6.5,6,6,6,6,6)  
true_mu2 <- c(2,2.5,2,2,2,2) 

true_sigma1 <- diag(6) * 2
true_sigma2 <- diag(6)

# Data simulated is saved as 'simulated_counts'
simulated_counts <- Datagenerator_mpln(N = 50, d = 6, pi_g = c(0.79,0.21), means = rbind(true_mu1,true_mu2), sigmas = rbind(true_sigma1,true_sigma2), ProduceImage="Yes")

# Checking/Loading needed packages
LoadCheckPkg(pckgs=c("parallel","rstan","Rcpp","mclust","mvtnorm","edgeR","capushe","clusterGeneration","coda"))

# Making RStan model 
mod = stan_model("MPLN.stan")

# Calculate the number of cores
no_cores = detectCores()-1

# Initiate cluster
cl = makeCluster(no_cores) 

# Doing clusterExport
clusterExport(cl,c("mod", "simulated_counts","AIC_function","AIC3_function","BIC_function","calc_likelihood","calculate_parameters","calling_clustering","cluster_mpln","ICL_function",
  "initializationrun","main_mpln","LoadCheckPkg","stanrun","zvalue_calculation"))

# Doing clusterEvalQ
clusterEvalQ(cl, library(parallel))
clusterEvalQ(cl, library(rstan))
clusterEvalQ(cl, library(Rcpp))
clusterEvalQ(cl, library(mclust))
clusterEvalQ(cl, library(mvtnorm))
clusterEvalQ(cl, library(edgeR))
clusterEvalQ(cl, library(capushe))
clusterEvalQ(cl, library(clusterGeneration))
clusterEvalQ(cl, library(coda))

# Running clustering
MPLNClust_results <- main_mpln(dataset=simulated_counts$dataset, 
                               membership=simulated_counts$truemembership, 
                               Gmin=1, 
                               Gmax=1, 
                               n_chains=3, 
                               n_iterations=100, 
                               init_method="kmeans", 
                               n_init_iterations=0, 
                               normalize=NA)

# To visualize clustered data
visualize_mpln(dataset=simulated_counts$dataset, ClusterMembershipVector=MPLNClust_results$BIC.all$BICmodelselected_labels)

#Saving results with date as .RData file
#save.image(paste0("MPLNClust_results,"_",format(Sys.time(), "%d%b%Y"),".RData"))