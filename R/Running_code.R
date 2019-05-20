# **Should have all the function R files and MPLN.stan in the same working directory as this file

source("Setup.R")

#####################################  DATA GENERATION/LOADING  #####################################
# Values for data simulation
true_mu1 <- c(6.5,6,6,6,6,6)  
true_mu2 <- c(2,2.5,2,2,2,2) 

true_sigma1 <- diag(6) * 2
true_sigma2 <- diag(6)

# Data simulated is saved as 'simulated_counts'
simulated_counts <- Datagenerator_mpln(N = 50, d = 6, pi_g = c(0.79,0.21), means = rbind(true_mu1,true_mu2), sigmas = rbind(true_sigma1,true_sigma2), ProduceImage="Yes")
#####################################################################################################

# Making RStan model 
mod = stan_model("MPLN.stan")

# Calculate the number of cores
no_cores = detectCores()-1

# Initiate cluster
cl = makeCluster(no_cores) 

# Doing clusterExport
clusterExport(cl,c("calc_likelihood","calculate_parameters","calling_clustering","cluster_mpln","initializationrun","main_mpln","mod","AIC_function","AIC3_function","BIC_function","ICL_function","remove_zero_counts","stanrun","zvalue_calculation"))

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
                               n_iterations=200, 
                               init_method="kmeans", 
                               n_init_iterations=5, 
                               normalize="TMM")

# To visualize clustered data
visualize_mpln(dataset=simulated_counts$dataset, ClusterMembershipVector=MPLNClust_results$BIC.all$BICmodelselected_labels)

#Saving results with date as .RData file
save.image(paste0("MPLNClust_results_",format(Sys.time(), "%d%b%Y"),".RData"))
