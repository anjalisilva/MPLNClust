#' Model Selection Via Akaike Information Criterion
#'
#' Performs model selection using Akaike Information Criterion (AIC).
#' Formula: - 2 * logLikelihood + 2 * nParameters.
#'
#' @param logLikelihood A vector with value of final log-likelihoods for
#'      each cluster size.
#' @param nParameters A vector with number of parameters for each
#'      cluster size.
#' @param clusterRunOutput Output from mplnVariational, mplnMCMCParallel, or
#'    mplnMCMCNonParallel, if available. Default value is NA. If provided,
#'    the vector of cluster labels obtained by mclust::map() for best model
#'    will be provided in the output.
#' @param gmin A positive integer specifying the minimum number of components
#'    to be considered in the clustering run.
#' @param gmax A positive integer, >gmin, specifying the maximum number of
#'    components to be considered in the clustering run.
#' @param parallel TRUE or FALSE indicating if MPLNClust::mplnMCMCParallel
#'    has been used.
#'
#' @return Returns an S3 object of class MPLN with results.
#' \itemize{
#'   \item allAICvalues - A vector of AIC values for each cluster size.
#'   \item AICmodelselected - An integer specifying model selected by AIC.
#'   \item AICmodelSelectedLabels - A vector of integers specifying cluster labels
#'     for the model selected. Only provided if user input clusterRunOutput.
#'   \item AICMessage - A character vector indicating if spurious clusters are
#'     detected. Otherwise, NA.
#' }
#'
#' @examples
#' # Generating simulated data
#'
#' trueMu1 <- c(6.5, 6, 6, 6, 6, 6)
#' trueMu2 <- c(2, 2.5, 2, 2, 2, 2)
#'
#' trueSigma1 <- diag(6) * 2
#' trueSigma2 <- diag(6)
#'
#' sampleData <- MPLNClust::mplnDataGenerator(nObservations = 100,
#'                                  dimensionality = 6,
#'                                  mixingProportions = c(0.79, 0.21),
#'                                  mu = rbind(trueMu1, trueMu2),
#'                                  sigma = rbind(trueSigma1, trueSigma2),
#'                                  produceImage = "No")
#'
#' # Clustering
#' mplnResults <- MPLNClust::mplnVariational(dataset = sampleData$dataset,
#'                                 membership = sampleData$trueMembership,
#'                                 gmin = 1,
#'                                 gmax = 2,
#'                                 initMethod = "kmeans",
#'                                 nInitIterations = 2,
#'                                 normalize = "Yes")
#'
#' # Model selection
#'  AICmodel <- MPLNClust::AICFunction(logLikelihood = mplnResults$logLikelihood,
#'                          nParameters = mplnResults$numbParameters,
#'                          clusterRunOutput = mplnResults$allResults,
#'                          gmin = mplnResults$gmin,
#'                          gmax = mplnResults$gmax,
#'                          parallel = FALSE)
#'
#' @author {Anjali Silva, \email{anjali.silva@uhnresearch.ca}}
#'
#' @references
#'
#' Akaike, H. (1973). Information theory and an extension of the maximum likelihood
#' principle. In \emph{Second International Symposium on Information Theory}, New York, NY,
#' USA, pp. 267–281. Springer Verlag.
#'
#' @export
#'
AICFunction <- function(logLikelihood,
                        nParameters,
                        clusterRunOutput = NA,
                        gmin,
                        gmax,
                        parallel = FALSE) {

  # Performing checks
  if(is.numeric(gmin) != TRUE || is.numeric(gmax) != TRUE) {
    stop("Class of gmin and gmin should be numeric.")
  }

  if (gmax < gmin) {
    stop("gmax cannot be less than gmin.")
  }

  if(is.numeric(nParameters) != TRUE) {
    stop("nParameters should be a vector of integers indicating
      number of parameters for each cluster size.")
  }

  if(is.numeric(logLikelihood) != TRUE) {
    stop("logLikelihood should be a vector of numeric values.")
  }

  if(length(logLikelihood) != (gmax - gmin + 1)) {
    stop("logLikelihood should be a vector of length (gmax - gmin + 1).")
  }

  if(is.logical(parallel) != TRUE) {
    stop("Should be logical, TRUE or FALSE indicating if
          MPLNClust::mplnMCMCParallel has been used.")
  }



  AIC <- - 2 * logLikelihood + 2 * nParameters
  AICmodel <- seq(gmin, gmax, 1)[grep(min(AIC, na.rm = TRUE), AIC)]
  AICMessage <- NA # For spurious clusters


  # Obtain cluster labels for best model if clusterRunOutput is provided
  if(all(is.na(clusterRunOutput)) == FALSE) {
    if(isTRUE(parallel) == "FALSE") {
      # If non parallel MCMC-EM or Variational EM run
      AICmodelLabels <- clusterRunOutput[[grep(min(AIC, na.rm = TRUE), AIC)]]$clusterlabels
    } else {
      # If parallel MCMC-EM run
      AICmodelLabels <- clusterRunOutput[[grep(min(AIC, na.rm = TRUE), AIC)]]$allResults$clusterlabels
    }

    # Check for spurious clusters, only possible if cluster labels provided
    if (max(AICmodelLabels) != AICmodel) {
      AICmodel <- max(AICmodelLabels)
      AICMessage <- "Spurious or empty cluster resulted."
    }
  } else {
    AICmodelLabels <- "clusterRunOutput not provided"
  }


  AICresults<-list(allAICvalues = AIC,
                   AICmodelselected = AICmodel,
                   AICmodelSelectedLabels = AICmodelLabels,
                   AICMessage = AICMessage)
  class(AICresults) <- "AIC"
  return(AICresults)
}


#' Model Selection Via Akaike Information Criterion by Bozdogan (1994)
#'
#' Performs model selection using Akaike Information Criterion by
#' Bozdogan (1994), called AIC3. Formula: - 2 * logLikelihood + 3 * nParameters.
#'
#' @param logLikelihood A vector with value of final log-likelihoods for
#'      each cluster size.
#' @param nParameters A vector with number of parameters for each
#'      cluster size.
#' @param clusterRunOutput Output from mplnVariational, mplnMCMCParallel, or
#'    mplnMCMCNonParallel, if available. Default value is NA. If provided,
#'    the vector of cluster labels obtained by mclust::map() for best model
#'    will be provided in the output.
#' @param gmin A positive integer specifying the minimum number of components
#'    to be considered in the clustering run.
#' @param gmax A positive integer, >gmin, specifying the maximum number of
#'    components to be considered in the clustering run.
#' @param parallel TRUE or FALSE indicating if MPLNClust::mplnMCMCParallel
#'    has been used.
#'
#' @return Returns an S3 object of class MPLN with results.
#' \itemize{
#'   \item allAIC3values - A vector of AIC3 values for each cluster size.
#'   \item AIC3modelselected - An integer specifying model selected by AIC3.
#'   \item AIC3modelSelectedLabels - A vector of integers specifying cluster labels
#'     for the model selected. Only provided if user input clusterRunOutput.
#'   \item AIC3Message - A character vector indicating if spurious clusters are
#'     detected. Otherwise, NA.
#' }
#'
#' @examples
#' # Generating simulated data
#'
#' trueMu1 <- c(6.5, 6, 6, 6, 6, 6)
#' trueMu2 <- c(2, 2.5, 2, 2, 2, 2)

#' trueSigma1 <- diag(6) * 2
#' trueSigma2 <- diag(6)

#' sampleData <- MPLNClust::mplnDataGenerator(nObservations = 100,
#'                                  dimensionality = 6,
#'                                  mixingProportions = c(0.79, 0.21),
#'                                  mu = rbind(trueMu1, trueMu2),
#'                                  sigma = rbind(trueSigma1, trueSigma2),
#'                                  produceImage = "No")

#' # Clustering
#' mplnResults <- MPLNClust::mplnVariational(dataset = sampleData$dataset,
#'                                 membership = sampleData$trueMembership,
#'                                 gmin = 1,
#'                                 gmax = 2,
#'                                 initMethod = "kmeans",
#'                                 nInitIterations = 2,
#'                                 normalize = "Yes")
#'
#' # Model selection
#' AIC3model <- MPLNClust::AIC3Function(logLikelihood = mplnResults$logLikelihood,
#'                            nParameters = mplnResults$numbParameters,
#'                            clusterRunOutput = mplnResults$allResults,
#'                            gmin = mplnResults$gmin,
#'                            gmax = mplnResults$gmax,
#'                            parallel = FALSE)
#'
#' @author {Anjali Silva, \email{anjali.silva@uhnresearch.ca}}
#'
#' @references
#'
#' Akaike, H. (1973). Information theory and an extension of the maximum likelihood
#' principle. In \emph{Second International Symposium on Information Theory}, New York, NY,
#' USA, pp. 267–281. Springer Verlag.
#'
#' #' Bozdogan, H. (1994). Mixture-model cluster analysis using model selection criteria
#' and a new informational measure of complexity. In \emph{Proceedings of the First US/Japan
#' Conference on the Frontiers of Statistical Modeling: An Informational Approach:
#' Volume 2 Multivariate Statistical Modeling}, pp. 69–113. Dordrecht: Springer Netherlands.
#'
#' @export
#'
AIC3Function <- function(logLikelihood,
                         nParameters,
                         clusterRunOutput = NA,
                         gmin,
                         gmax,
                         parallel = FALSE) {

  # Performing checks
  if(is.numeric(gmin) != TRUE || is.numeric(gmax) != TRUE) {
    stop("Class of gmin and gmin should be numeric.")
  }

  if (gmax < gmin) {
    stop("gmax cannot be less than gmin.")
  }

  if(is.numeric(nParameters) != TRUE) {
    stop("nParameters should be a vector of integers indicating
      number of parameters for each cluster size.")
  }

  if(is.numeric(logLikelihood) != TRUE) {
    stop("logLikelihood should be a vector of numeric values.")
  }

  if(length(logLikelihood) != (gmax - gmin + 1)) {
    stop("logLikelihood should be a vector of length (gmax - gmin + 1).")
  }

  if(is.logical(parallel) != TRUE) {
    stop("Should be logical, TRUE or FALSE indicating if
      MPLNClust::mplnMCMCParallel has been used.")
  }


  AIC3 <- - 2 * logLikelihood + 3 * nParameters
  AIC3model <- seq(gmin, gmax, 1)[grep(min(AIC3,na.rm = TRUE), AIC3)]
  AIC3Message <- NA # For spurious clusters

  # Obtain cluster labels for best model if clusterRunOutput is provided
  if(all(is.na(clusterRunOutput)) == FALSE) {
    if(isTRUE(parallel) == "FALSE") {
      # If non parallel MCMC-EM or Variational EM run
    AIC3modelLabels <- clusterRunOutput[[grep(min(AIC3,na.rm = TRUE), AIC3)]]$clusterlabels
  } else {
    # If parallel MCMC-EM run
    AIC3modelLabels <- clusterRunOutput[[grep(min(AIC3,na.rm = TRUE), AIC3)]]$allResults$clusterlabels
  }
    # Check for spurious clusters, only possible if cluster labels provided
    if (max(AIC3modelLabels) != AIC3model) {
      AIC3model <- max(AIC3modelLabels)
      AIC3Message <- "Spurious or empty cluster resulted."
    }
  } else {
    AIC3modelLabels <- "clusterRunOutput not provided"
  }


  AIC3results <- list(allAIC3values = AIC3,
                      AIC3modelselected = AIC3model,
                      AIC3modelSelectedLabels = AIC3modelLabels,
                      AIC3Message = AIC3Message)
  class(AIC3results) <- "AIC3"
  return(AIC3results)
}



#' Model Selection Via Bayesian Information Criterion
#'
#' Performs model selection using Bayesian Information Criterion (BIC) by
#' Schwarz (1978). Formula: - 2 * logLikelihood + (nParameters * log(nObservations)).
#'
#' @param logLikelihood A vector with value of final log-likelihoods for
#'      each cluster size.
#' @param nParameters A vector with number of parameters for each
#'      cluster size.
#' @param nObservations A positive integer specifying the number of observations
#'      in the dataset analyzed.
#' @param clusterRunOutput Output from mplnVariational, mplnMCMCParallel, or
#'    mplnMCMCNonParallel, if available. Default value is NA. If provided,
#'    the vector of cluster labels obtained by mclust::map() for best model
#'    will be provided in the output.
#' @param gmin A positive integer specifying the minimum number of components
#'    to be considered in the clustering run.
#' @param gmax A positive integer, >gmin, specifying the maximum number of
#'    components to be considered in the clustering run.
#' @param parallel TRUE or FALSE indicating if MPLNClust::mplnMCMCParallel
#'    has been used.
#'
#' @return Returns an S3 object of class MPLN with results.
#' \itemize{
#'   \item allBICvalues - A vector of BIC values for each cluster size.
#'   \item BICmodelselected - An integer specifying model selected by BIC
#'   \item BICmodelSelectedLabels - A vector of integers specifying cluster labels
#'     for the model selected. Only provided if user input clusterRunOutput.
#'   \item BICMessage - A character vector indicating if spurious clusters are
#'     detected. Otherwise, NA.
#' }
#'
#' @examples
#' # Generating simulated data
#'
#' trueMu1 <- c(6.5, 6, 6, 6, 6, 6)
#' trueMu2 <- c(2, 2.5, 2, 2, 2, 2)

#' trueSigma1 <- diag(6) * 2
#' trueSigma2 <- diag(6)

#' sampleData <- MPLNClust::mplnDataGenerator(nObservations = 100,
#'                                  dimensionality = 6,
#'                                  mixingProportions = c(0.79, 0.21),
#'                                  mu = rbind(trueMu1, trueMu2),
#'                                  sigma = rbind(trueSigma1, trueSigma2),
#'                                  produceImage = "No")

#' # Clustering
#' mplnResults <- MPLNClust::mplnVariational(dataset = sampleData$dataset,
#'                                 membership = sampleData$trueMembership,
#'                                 gmin = 1,
#'                                 gmax = 2,
#'                                 initMethod = "kmeans",
#'                                 nInitIterations = 2,
#'                                 normalize = "Yes")
#'
#' # Model selection
#' BICmodel <- MPLNClust::BICFunction(logLikelihood = mplnResults$logLikelihood,
#'                          nParameters = mplnResults$numbParameters,
#'                          nObservations = nrow(mplnResults$dataset),
#'                          clusterRunOutput = mplnResults$allResults,
#'                          gmin = mplnResults$gmin,
#'                          gmax = mplnResults$gmax,
#'                          parallel = FALSE)
#'
#' @author {Anjali Silva, \email{anjali.silva@uhnresearch.ca}}
#'
#' @references
#'
#' Schwarz, G. (1978). Estimating the dimension of a model. \emph{The Annals of Statistics}
#' 6.
#'
#' @export
#'
BICFunction <- function(logLikelihood,
                        nParameters,
                        nObservations,
                        clusterRunOutput = NA,
                        gmin,
                        gmax,
                        parallel = FALSE) {

  # Performing checks
  if(is.numeric(gmin) != TRUE || is.numeric(gmax) != TRUE) {
    stop("Class of gmin and gmin should be numeric.")
  }

  if (gmax < gmin) {
    stop("gmax cannot be less than gmin.")
  }

  if(is.numeric(nParameters) != TRUE) {
    stop("nParameters should be a vector of integers indicating
      number of parameters for each cluster size.")
  }

  if(is.numeric(logLikelihood) != TRUE) {
    stop("logLikelihood should be a vector of numeric values.")
  }

  if(length(logLikelihood) != (gmax - gmin + 1)) {
    stop("logLikelihood should be a vector of length (gmax - gmin + 1).")
  }

  if(is.logical(parallel) != TRUE) {
    stop("Should be logical, TRUE or FALSE indicating if
      MPLNClust::mplnMCMCParallel has been used.")
  }



  BIC <- - 2 * logLikelihood + (nParameters * log(nObservations))
  BICmodel <- seq(gmin, gmax, 1)[grep(min(BIC, na.rm = TRUE), BIC)]
  BICMessage <- NA # For spurious clusters

  # Obtain cluster labels for best model if clusterRunOutput is provided
  if(all(is.na(clusterRunOutput)) == FALSE) {
    if(isTRUE(parallel) == "FALSE") {
      # If non parallel MCMC-EM or Variational EM run
      BICmodelLabels <- clusterRunOutput[[grep(min(BIC, na.rm = TRUE),
      BIC)]]$clusterlabels
  } else {
    # If parallel MCMC-EM run
    BICmodelLabels <- clusterRunOutput[[grep(min(BIC, na.rm = TRUE),
      BIC)]]$allResults$clusterlabels
  }
  # Check for spurious clusters, only possible if cluster labels provided
  if (max(BICmodelLabels) != BICmodel) {
      BICmodel <- max(BICmodelLabels)
      BICMessage <- "Spurious or empty cluster resulted."
    }
  } else {
    BICmodelLabels <- "clusterRunOutput not provided"
  }

  BICresults <- list(allBICvalues = BIC,
                     BICmodelselected = BICmodel,
                     BICmodelSelectedLabels = BICmodelLabels,
                     BICMessage = BICMessage)
  class(BICresults) <- "BIC"
  return(BICresults)
}



#' Model Selection Via Integrated Completed Likelihood
#'
#' Performs model selection using integrated completed likelihood (ICL) by
#' Biernacki et al., (2000).
#'
#' @param logLikelihood A vector with value of final log-likelihoods for
#'      each cluster size.
#' @param nParameters A vector with number of parameters for each
#'      cluster size.
#' @param nObservations A positive integer specifying the number of observations
#'      in the dataset analyzed.
#' @param clusterRunOutput Output from MPLNClust::mplnVariational,
#'    MPLNClust::mplnMCMCParallel, or MPLNClust::mplnMCMCNonParallel functions.
#'    Either clusterRunOutput or probaPost must be provided.
#' @param probaPost A list that is length (gmax - gmin + 1) containing posterior
#'    probability at each g, for g = gmin:gmax. This argument is useful if
#'    clustering output have been generated non-serially, e.g., g = 1:5 and
#'    g = 6:10 rather than g = 1:10. Either clusterRunOutput or probaPost
#'    must be provided.
#' @param gmin A positive integer specifying the minimum number of components
#'    to be considered in the clustering run.
#' @param gmax A positive integer, > gmin, specifying the maximum number of
#'    components to be considered in the clustering run.
#' @param parallel TRUE or FALSE indicating if MPLNClust::mplnMCMCParallel
#'    has been used.
#'
#' @return Returns an S3 object of class MPLN with results.
#' \itemize{
#'   \item allICLvalues - A vector of ICL values for each cluster size.
#'   \item ICLmodelselected - An integer specifying model selected by ICL.
#'   \item ICLmodelSelectedLabels - A vector of integers specifying cluster labels
#'     for the model selected. Only provided if user input clusterRunOutput.
#'   \item ICLMessage - A character vector indicating if spurious clusters are
#'     detected. Otherwise, NA.
#' }
#'
#' @examples
#' # Generating simulated data
#'
#' trueMu1 <- c(6.5, 6, 6, 6, 6, 6)
#' trueMu2 <- c(2, 2.5, 2, 2, 2, 2)
#'
#' trueSigma1 <- diag(6) * 2
#' trueSigma2 <- diag(6)
#'
#' sampleData <- MPLNClust::mplnDataGenerator(nObservations = 100,
#'                                            dimensionality = 6,
#'                                            mixingProportions = c(0.79, 0.21),
#'                                            mu = rbind(trueMu1, trueMu2),
#'                                            sigma = rbind(trueSigma1, trueSigma2),
#'                                            produceImage = "No")
#'
#' # Clustering
#' mplnResults <- MPLNClust::mplnVariational(dataset = sampleData$dataset,
#'                                           membership = sampleData$trueMembership,
#'                                           gmin = 1,
#'                                           gmax = 2,
#'                                           initMethod = "kmeans",
#'                                           nInitIterations = 2,
#'                                           normalize = "Yes")
#'
#' # Model selection
#' ICLmodel <- MPLNClust::ICLFunction(logLikelihood = mplnResults$logLikelihood,
#'                                    nParameters = mplnResults$numbParameters,
#'                                    nObservations = nrow(mplnResults$dataset),
#'                                    clusterRunOutput = mplnResults$allResults,
#'                                    gmin = mplnResults$gmin,
#'                                    gmax = mplnResults$gmax,
#'                                    parallel = FALSE)
#'
#' @author {Anjali Silva, \email{anjali.silva@uhnresearch.ca}}
#'
#' @references
#'
#' Biernacki, C., G. Celeux, and G. Govaert (2000). Assessing a mixture model for
#' clustering with the integrated classification likelihood. \emph{IEEE Transactions
#' on Pattern Analysis and Machine Intelligence} 22.
#'
#' @export
#' @importFrom mclust unmap
#'
ICLFunction <- function(logLikelihood,
                        nParameters,
                        nObservations,
                        clusterRunOutput = NA,
                        probaPost = NA,
                        gmax,
                        gmin,
                        parallel = FALSE) {

  # Performing checks
  if (gmax < gmin) {
    stop("gmax cannot be less than gmin.")
  }

  if(is.numeric(nParameters) != TRUE) {
    stop("nParameters should be a vector of integers indicating
      number of parameters for each cluster size.")
  }

  if(is.numeric(logLikelihood) != TRUE) {
    stop("logLikelihood should be a vector of numeric values.")
  }

  if(length(logLikelihood) != (gmax - gmin + 1)) {
    stop("logLikelihood should be a vector of length (gmax - gmin + 1).")
  }

  if(is.logical(parallel) != TRUE) {
    stop("Should be logical, TRUE or FALSE indicating if
      MPLNClust::mplnMCMCParallel has been used.")
  }

  if(all(is.na(probaPost)) != TRUE) {
    if(length(probaPost) != (gmax - gmin + 1)) {
      stop("probaPost must be a list of length (gmax - gmin + 1)
      containing posterior probability at each g.")
    }
  }

  if(all(is.na(clusterRunOutput)) == TRUE && all(is.na(probaPost)) == TRUE) {
    stop("Either clusterRunOutput or probaPost must be provided.")
  }

  BIC <- - 2 * logLikelihood + (nParameters * log(nObservations))

  ICL <- vector()

  # if clusterRunOutput is provided by user
  if(all(is.na(clusterRunOutput)) != TRUE) {
    for (g in 1:(gmax - gmin + 1)) {
      if(isTRUE(parallel) == "FALSE") {
        # If non parallel run
        z <- clusterRunOutput[[g]]$probaPost
        mapz <- mclust::unmap(clusterRunOutput[[g]]$clusterlabels)
      } else {
          # If parallel run
          z <- clusterRunOutput[[g]]$allResults$probaPost
          mapz <- mclust::unmap(clusterRunOutput[[g]]$allResults$clusterlabels)
      }
      forICL <- function(g){sum(log(z[which(mapz[, g] == 1), g]))}
      ICL[g] <- BIC[g] + sum(sapply(1:ncol(mapz), forICL))
    }
    ICLmodel <- seq(gmin, gmax, 1)[grep(min(ICL, na.rm = TRUE), ICL)]

    if(isTRUE(parallel) == "FALSE") {
      # If non parallel MCMC-EM or Variational EM run
      ICLmodelLabels <- clusterRunOutput[[grep(min(ICL, na.rm = TRUE),
        ICL)]]$clusterlabels
    } else {
      # If parallel MCMC-EM
      ICLmodelLabels <- clusterRunOutput[[grep(min(ICL, na.rm = TRUE),
        ICL)]]$allResults$clusterlabels
    }
  }

  # if probaPost is provided by user
  if(all(is.na(probaPost)) != TRUE) {

    for (g in 1:(gmax - gmin + 1)) {
      z <- probaPost[[g]]
      mapz <- mclust::unmap(mclust::map(probaPost[[g]]))
      forICL <- function(g){sum(log(z[which(mapz[, g] == 1), g]))}
      ICL[g] <- BIC[g] + sum(sapply(1:ncol(mapz), forICL))
    }
    ICLmodel <- seq(gmin, gmax, 1)[grep(min(ICL, na.rm = TRUE), ICL)]
    ICLmodelLabels <- mclust::map(probaPost[[grep(min(ICL, na.rm = TRUE), ICL)]])
  }

  # Check for spurious clusters
  ICLMessage <- NA
  if (max(ICLmodelLabels) != ICLmodel) {
    ICLmodel <- max(ICLmodelLabels)
    ICLMessage <- "Spurious or empty cluster resulted."
  }

  ICLresults <- list(allICLvalues = ICL,
                     ICLmodelselected = ICLmodel,
                     ICLmodelSelectedLabels = ICLmodelLabels,
                     ICLMessage = ICLMessage)
  class(ICLresults) <- "ICL"
  return(ICLresults)
}
# [END]
