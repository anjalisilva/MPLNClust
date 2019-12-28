# Generating simulated data

trueMu1 <- c(6.5, 6, 6, 6, 6, 6)  
trueMu2 <- c(2, 2.5, 2, 2, 2, 2) 

trueSigma1 <- diag(6) * 2
trueSigma2 <- diag(6)

simulated_counts <- mplnDataGenerator(nObservations = 70, 
                                      dimensionality = 6, 
                                      mixingProportions = c(0.79, 0.21), 
                                      mu = rbind(trueMu1, trueMu2), 
                                      sigma = rbind(trueSigma1, trueSigma2), 
                                      produceImage = "No")

MPLNClustResults <- mpln(dataset = simulated_counts$dataset, 
                         membership = "none", 
                         gmin = 1, 
                         gmax = 1, 
                         nChains = 3, 
                         nIterations = 300,
                         initMethod = "kmeans", 
                         nInitIterations = 0,
                         normalize = "Yes")

MPLNVisuals <- mplnVisualize(dataset = MPLNClustResults$dataset,
  clusterMembershipVector =
    MPLNClustResults$all_results[[1]]$all_results$clusterlabels,
  fileName = 'RegularDataSet', plots = 'line',
  format = 'png')
