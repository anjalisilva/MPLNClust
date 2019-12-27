#' Model-Based Clustering Using MPLN
#'
#' Performs clustering using mixtures of multivariate Poisson-log
#' normal (MPLN) distribution and model selection using AIC, AIC3,
#' BIC and ICL.
#'
#' @param dataset A dataset of class matrix and type integer such that
#'    rows correspond to observations and columns correspond to variables.
#' @param membership A numeric vector of length nrow(dataset) containing the
#'    cluster membership of each observation. If not available, leave as NA.
#' @param gmin A positive integer specifying the minimum number of components
#'    to be considered in the clustering run.
#' @param gmax A positive integer, >gmin, specifying the maximum number of
#'    components to be considered in the clustering run.
#' @param nChains A positive integer specifying the number of Markov chains.
#'    Default is 3, the recommended minimum number.
#' @param nIterations A positive integer specifying the number of iterations
#'    for each chain (including warmup). The value should be greater than 40.
#'    The upper limit will depend on size of dataset.
#' @param initMethod An algorithm for initialization. Current options are
#'    "kmeans", "random", "medoids", "clara", or "fanny". Default is "kmeans".
#' @param nInitIterations A positive integer specifying the number of
#'     initialization runs to be considered.
#' @param normalize A string with options "Yes" or "No" specifying
#'     if normalization should be performed. Currently, normalization factors
#'     are calculated using TMM method of edgeR package. Default is "Yes".
#'
#' @return Returns an S3 object of class MPLN with results.
#' \itemize{
#'   \item dataset - The input dataset on which clustering is performed.
#'   \item dimensionality - Dimensionality of the input dataset.
#'   \item normalization_factors - A vector of normalization factors used
#'      for input dataset.
#'   \item gmin - Minimum number of components considered in the clustering
#'      run
#'   \item gmax - Maximum number of components considered in the clustering
#'      run
#'   \item initalization_method - Method used for initialization.
#'   \item all_results - A list with all results.
#'   \item loglikelihood - A vector with value of final log-likelihoods for
#'      each cluster size.
#'   \item numb_of_parameters - A vector with number of parameters for each
#'      cluster size.
#'   \item true_labels - The vector of true labels, if provided by user.
#'   \item ICL_all - A list with all ICL model selection results.
#'   \item BIC_all - A list with all BIC model selection results.
#'   \item AIC_all - A list with all AIC model selection results.
#'   \item AIC3_all - A list with all AIC3 model selection results.
#'   \item slope_heuristics - If more than 10 models are considered, slope heuristic
#'      results.
#'   \item Djumpmodel_selected - If more than 10 models are considered, slope heuristic
#'      results.
#'   \item DDSEmodel_selected - If more than 10 models are considered, slope heuristic
#'      results.
#'   \item total_time - Total time used for clustering and model selection.
#' }
#'
#' @examples
#' # Generating simulated data
#' # Not run
#' # trueMu1 <- c(6.5, 6, 6, 6, 6, 6)
#' # trueMu2 <- c(2, 2.5, 2, 2, 2, 2)
#'
#' # trueSigma1 <- diag(6) * 2
#' # trueSigma2 <- diag(6)
#'
#' # sampleData <- mplnDataGenerator(nObservations = 100,
#' #                                 dimensionality = 6,
#' #                                 mixingProportions = c(0.79, 0.21),
#' #                                 mu = rbind(trueMu1, trueMu2),
#' #                                 sigma = rbind(trueSigma1, trueSigma2),
#' #                                 produceImage = "No")
#'
#' # Running clustering
#' # mplnResults <- mpln(dataset = sampleData$dataset,
#' #                     membership = sampleData$truemembership,
#' #                     gmin = 1,
#' #                     gmax = 2,
#' #                     nChains = 3,
#' #                     nIterations = 1000,
#' #                     initMethod = "kmeans",
#' #                     nInitIterations = 2,
#' #                     normalize = "Yes")
#'
#' @references
#' Aitchison, J. and C. H. Ho (1989). The multivariate Poisson-log normal distribution.
#' \emph{Biometrika} 76.
#'
#' Akaike, H. (1973). Information theory and an extension of the maximum likelihood
#' principle. In \emph{Second International Symposium on Information Theory}, New York, NY,
#' USA, pp. 267–281. Springer Verlag.
#'
#' Biernacki, C., G. Celeux, and G. Govaert (2000). Assessing a mixture model for
#' clustering with the integrated classification likelihood. \emph{IEEE Transactions
#' on Pattern Analysis and Machine Intelligence} 22.
#'
#' Bozdogan, H. (1994). Mixture-model cluster analysis using model selection criteria
#' and a new informational measure of complexity. In \emph{Proceedings of the First US/Japan
#' Conference on the Frontiers of Statistical Modeling: An Informational Approach:
#' Volume 2 Multivariate Statistical Modeling}, pp. 69–113. Dordrecht: Springer Netherlands.
#'
#' Schwarz, G. (1978). Estimating the dimension of a model. \emph{The Annals of Statistics}
#' 6.
#'
#' Silva, A. et al. (2019). A multivariate Poisson-log normal mixture model
#' for clustering transcriptome sequencing data. \emph{BMC Bioinformatics} 20.
#' \href{https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-2916-0}{Link}
#'
#' @export
#' @import coda
#' @importFrom capushe capushe
#' @import cluster
#' @import clusterGeneration
#' @importFrom edgeR calcNormFactors
#' @import MASS
#' @import mclust
#' @import mvtnorm
#' @import Rcpp
#' @importFrom rstan sampling
#' @importFrom rstan stan_model
#' @import parallel
#' @import stats
#' @importFrom utils tail
#' @importFrom utils write.csv
#'
mpln <- function(dataset, membership = NA, gmin = 1, gmax = 2,
                 nChains = 3, nIterations = NA,
                 initMethod = "kmeans", nInitIterations = 0,
                 normalize = "Yes") {

  ptm <- proc.time()

  # Performing checks
  if (typeof(dataset) != "double" & typeof(dataset) != "integer") {
    stop("Dataset type needs to be integer")
  }

  if (class(dataset) != "matrix") {
    stop("Dataset needs to be a matrix")
  }

  dataset <- removeZeroCounts(dataset = dataset)$dataset

  if (gmax < gmin) {
    stop("gmax cannot be less than gmin")
  }

  if(is.na(nIterations)) {
    nIterations <- 1000
  }

  # if(is.na(nChains) || nChains<3) {
  #  nChains <- 3
  #  cat("Recommended number of chains is minimum 3. nChains' set to 3")}

  if(nIterations < 40) {
    stop("RStan nIterations argument should be greater than 40")}


  dimensionality <- ncol(dataset)
  nObservations <- nrow(dataset)

  if(all(is.na(membership) != TRUE) && length(membership) != nObservations) {
    stop("Length of membership character vector and
      sample size of dataset should match")
  }

  if(all(is.na(membership) != TRUE) &&
      all((diff(sort(unique(membership))) == 1) != TRUE) ) {
    stop("Cluster memberships in the membership vector
      are missing a cluster, e.g. 1,3,4,5,6 is missing cluster 2")
  }

  if(length(which(apply(dataset, 1, function(x) all(x == 0)) == TRUE)) != 0) {
    cat("\nDataset row(s)",
      c(which(apply(dataset, 1, function(x) all(x == 0)) == TRUE)),
      "will be removed as this/these contain(s) all zeros")

    if(all(is.na(membership) == FALSE)) {
      membership <- membership[- c(which(apply(dataset, 1, function(x)
        all(x == 0)) == TRUE))]
    }
    dataset <- dataset[- c(which(apply(dataset, 1, function(x)
      all(x == 0)) == TRUE)), ]
    nObservations <- nrow(dataset)
  }

  if(all(is.na(membership) == TRUE)) {
    membership <- "Not provided"
  }

  if (gmax > nObservations) {
    stop("gmax cannot be larger than nrow(dataset)")
  }

  # Calculating normalization factors
  if(normalize == "Yes") {
    normFactors <- log(as.vector(edgeR::calcNormFactors(as.matrix(dataset),
      method = "TMM")))
  } else if(normalize == "No") {
    normFactors <- rep(0, dimensionality)
  } else{
    stop("normFactors should be 'Yes' or 'No' ")
  }

  # Construct a Stan model
  stancode <- 'data{int<lower=1> d; // Dimension of theta
                    int<lower=0> N; //Sample size
                    int y[N,d]; //Array of Y
                    vector[d] mu;
                    matrix[d,d] Sigma;
                    vector[d] normfactors;
                }
                parameters{
                    matrix[N,d] theta;
                }
                model{ //Change values for the priors as appropriate
                    for (n in 1:N){
                      theta[n,]~multi_normal(mu,Sigma);
                    }
                    for (n in 1:N){
                      for (k in 1:d){
                        real z;
                        z=exp(normfactors[k]+theta[n,k]);
                        y[n,k]~poisson(z);
                      }
                    }
                }'

  mod <- rstan::stan_model(model_code = stancode, verbose = FALSE)

  # Constructing parallel code
  mplnParallel <- function(g) {
    ## ** Never use set.seed(), use clusterSetRNGStream() instead,
    # to set the cluster seed if you want reproducible results
    # clusterSetRNGStream(cl = cl, iseed = g)
    mplnParallelRun <- callingClustering(data = dataset,
                                        gmin = g,
                                        gmax = g,
                                        nChains = nChains,
                                        initMethod = initMethod,
                                        nInitIterations = nInitIterations,
                                        nIterations = nIterations,
                                        normFactors = normFactors,
                                        model = mod)
    return(mplnParallelRun)
  }

  # Calculate the number of cores and initiate cluster
  cl <- parallel::makeCluster(parallel::detectCores() - 1)

  # Doing clusterExport
  parallel::clusterExport(cl = cl,
                          varlist = c("nInitIterations",
                                      "dataset",
                                      "initMethod",
                                      "normFactors",
                                      "nChains",
                                      "nIterations",
                                      "calcLikelihood",
                                      "calcParameters",
                                      "callingClustering",
                                      "mplnCluster",
                                      "initializationRun",
                                      "mpln",
                                      "mod",
                                      "AICFunction",
                                      "AIC3Function",
                                      "BICFunction",
                                      "ICLFunction",
                                      "removeZeroCounts",
                                      "stanRun",
                                      "calcZvalue"),
                                      envir = environment())

  # Doing clusterEvalQ
  parallel::clusterEvalQ(cl, library(parallel))
  parallel::clusterEvalQ(cl, library(rstan))
  parallel::clusterEvalQ(cl, library(Rcpp))
  parallel::clusterEvalQ(cl, library(mclust))
  parallel::clusterEvalQ(cl, library(mvtnorm))
  parallel::clusterEvalQ(cl, library(edgeR))
  parallel::clusterEvalQ(cl, library(capushe))
  parallel::clusterEvalQ(cl, library(clusterGeneration))
  parallel::clusterEvalQ(cl, library(coda))

  parallelRun = list()
  cat("\nRunning parallel code now.")
  parallelRun = parallel::clusterMap(cl = cl,
                                     fun = mplnParallel,
                                     g = gmin:gmax)
  cat("\nDone parallel code.")
  #parallel::stopCluster(cl)

  BIC <- ICL <- AIC <- AIC3 <- Djump <- DDSE <- k <- ll <- vector()

  for(g in seq_along(1:(gmax - gmin + 1))) {
    # save the final log-likelihood
    ll[g] <- unlist(tail(parallelRun[[g]]$all_results$loglikelihood,
             n = 1))

    k[g] <- calcParameters(numberG = g,
            dimensionality = dimensionality)

    if (g == max(1:(gmax - gmin + 1))) { # starting model selection
      bic <- BICFunction(ll = ll,
                         k = k,
                         n = nObservations,
                         run = parallelRun,
                         gmin = gmin,
                         gmax = gmax)

      icl <- ICLFunction(bIc = bic,
                         gmin = gmin,
                         gmax = gmax,
                         run = parallelRun)

      aic <- AICFunction(ll = ll,
                         k = k,
                         run = parallelRun,
                         gmin = gmin,
                         gmax = gmax )

      aic3 <- AIC3Function(ll = ll,
                           k = k,
                           run = parallelRun,
                           gmin = gmin,
                           gmax = gmax)
    }
  }

  # for Djump and DDSE
  if((gmax - gmin + 1) > 10 ) {
    # adapted based on HTSCluster package 2.0.8 (25 Oct 2016)
    PMM <- parallelRun
    runs <- gmin:gmax
    gmax <- gmax
    logLike.final <- suppressWarnings(do.call("cbind",
      lapply(PMM, function(x) x$loglikelihood)))
    # gives log-likelihood for each cluster at each run
    logLike.val <- apply(logLike.final, 1, max)


    message("Note: diagnostic plots for results corresponding
      to model selection via slope heuristics (Djump and DDSE)
      should be examined to ensure that sufficiently complex
      models have been considered.")
    Kchoice <- gmin:gmax
    k <- k # number of parameters
    mat <- cbind(Kchoice, k/nObservations, k/nObservations,
                 - logLike.val)
    ResCapushe <- capushe::capushe(mat, nObservations)
    DDSEmodel <- ResCapushe@DDSE@model
    Djumpmodel <- ResCapushe@Djump@model
    final <- proc.time() - ptm

    RESULTS <- list(dataset = dataset,
                    dimensionality = dimensionality,
                    normalization_factors = normFactors,
                    gmin = gmin,
                    gmax = gmax,
                    initalization_method = initMethod,
                    all_results = parallelRun,
                    loglikelihood = ll,
                    numb_of_parameters = k,
                    true_labels = membership,
                    ICL_all = icl,
                    BIC_all = bic,
                    AIC_all = aic,
                    AIC3_all = aic3,
                    slope_heuristics = ResCapushe,
                    Djumpmodel_selected = ResCapushe@Djump@model,
                    DDSEmodel_selected = ResCapushe@DDSE@model,
                    total_time = final)

  } else {# end of Djump and DDSE
    final <- proc.time() - ptm
    RESULTS <- list(dataset = dataset,
                    dimensionality = dimensionality,
                    normalization_factors=normFactors,
                    gmin = gmin,
                    gmax = gmax,
                    initalization_method = initMethod,
                    all_results = parallelRun,
                    loglikelihood = ll,
                    numb_of_parameters = k,
                    true_labels = membership,
                    ICL_all = icl,
                    BIC_all = bic,
                    AIC_all = aic,
                    AIC3_all = aic3,
                    slope_heuristics = "Not used",
                    total_time = final)
  }

  class(RESULTS) <- "MPLN"
  return(RESULTS)
  # Developed by Anjali Silva
}
