
source("AIC3Function.R")
source("AICFunction.R")
source("BICFunction.R")
source("CalcLikelihood.R")
source("CalculateParameters.R")
source("CallingClustering.R")
source("ClusterMPLN.R")
source("ICLFunction.R")
source("InitializationRun.R")
source("MPLNClustering.R")
source("MPLNDataGenerator.R")
source("PackageCheck.R")
source("StanRun.R")
source("ZValueCalculation")

# Generating data

true_mu1 <- c(6.5,6,6,6,6,6)  
true_mu2 <- c(2,2.5,2,2,2,2) 

true_sigma1 <- diag(6) * 2
true_sigma2 <- diag(6)

simulated_counts <- Datagenerator(i = 1, N = 50, d = 6, pi_g = c(0.79,0.21), means = rbind(true_mu1,true_mu2), sigmas = rbind(true_sigma1,true_sigma2), ProduceImage="Yes")



testing_dataset <- simulated_counts$dataset # Assign test dataset using the variable name 'testing_dataset'
clus_results <- MPLNClustering(dataset = testing_dataset, Gmin = 1, Gmax = 5, n_chains = 3, n_iterations=1000, membership = NA, init_method = "kmeans", n_init_iterations = 5, normalize = "TMM")



