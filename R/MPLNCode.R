#' Model-Based Clustering Using MPLN via Parallel Performance
#'
#' Performs clustering using mixtures of multivariate Poisson-log
#' normal (MPLN) distribution and model selection using AIC, AIC3,
#' BIC and ICL. Coarse grain parallelization is employed, such that
#' when a range of components/clusters (g = 1,...,G) are considered, each
#' G is run on a different processor. This can be performed because
#' each component/cluster size is independent from another. All
#' components/clusters in the range to be tested have been parallelized
#' to run on a seperate core using the *parallel* R package. The number of
#' nodes to be used for clustering can be specified or calculated using
#' *parallel::detectCores() - 1*.
#'
#' @param dataset A dataset of class matrix and type integer such that
#'    rows correspond to observations and columns correspond to variables.
#' @param membership A numeric vector of length nrow(dataset) containing the
#'    cluster membership of each observation. If not available,
#'    leave as "none".
#' @param gmin A positive integer specifying the minimum number of components
#'    to be considered in the clustering run.
#' @param gmax A positive integer, >gmin, specifying the maximum number of
#'    components to be considered in the clustering run.
#' @param nChains A positive integer specifying the number of Markov chains.
#'    Default is 3, the recommended minimum number.
#' @param nIterations A positive integer specifying the number of iterations
#'    for each MCMC chain (including warmup). The value should be greater than
#'    40. The upper limit will depend on size of dataset.
#' @param initMethod An algorithm for initialization. Current options are
#'    "kmeans", "random", "medoids", "clara", or "fanny". Default is "kmeans".
#' @param nInitIterations A positive integer or zero, specifying the number
#'    of initialization runs to be considered. Default is 0.
#' @param normalize A string with options "Yes" or "No" specifying
#'     if normalization should be performed. Currently, normalization factors
#'     are calculated using TMM method of edgeR package. Default is "Yes".
#' @param numNodes A positive integer indicating the number of nodes to be
#'     used from the local machine to run the clustering algorithm. Else
#'     leave as NA, so default will be detected as
#'     parallel::makeCluster(parallel::detectCores() - 1).
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
#' # Clustering
#' # mplnResults <- mplnParallel(dataset = sampleData$dataset,
#' #                     membership = sampleData$trueMembership,
#' #                     gmin = 1,
#' #                     gmax = 2,
#' #                     nChains = 3,
#' #                     nIterations = 1000,
#' #                     initMethod = "kmeans",
#' #                     nInitIterations = 2,
#' #                     normalize = "Yes")
#'
#' @author {Anjali Silva, \email{anjali.silva@uhnresearch.ca}, Sanjeena Dang,
#'          \email{sdang@math.binghamton.edu}. }
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
#' @importFrom mclust unmap
#' @importFrom mclust map
#' @importFrom mvtnorm rmvnorm
#' @import Rcpp
#' @importFrom rstan sampling
#' @importFrom rstan stan_model
#' @import parallel
#' @import stats
#' @importFrom utils tail
#' @importFrom utils write.csv
#'
mplnParallel <- function(dataset, membership = "none", gmin = 1, gmax = 2,
  nChains = 3, nIterations = 1000,
  initMethod = "kmeans", nInitIterations = 0,
  normalize = "Yes", numNodes = NA) {

  ptm <- proc.time()

  # Performing checks
  if (typeof(dataset) != "double" & typeof(dataset) != "integer") {
    stop("Dataset type needs to be integer.")
  }

  if (class(dataset) != "matrix") {
    stop("Dataset needs to be a matrix.")
  }

  dataset <- removeZeroCounts(dataset = dataset)$dataset

  dimensionality <- ncol(dataset)
  nObservations <- nrow(dataset)

  if(class(gmin) != "numeric" || class(gmax) != "numeric") {
    stop("Class of gmin and gmin should be numeric.")
  }

  if (gmax < gmin) {
    stop("gmax cannot be less than gmin.")
  }

  if (class(nChains) != "numeric") {
    stop("nChains should be a positive integer of class numeric specifying
      the number of Markov chains.")
  }

  if(nChains < 3) {
    cat("nChains is less than 3. Note: recommended minimum number of
      chains is 3.")
  }

  if(class(nIterations) != "numeric") {
   stop("nIterations should be a positive integer of class numeric,
      specifying the number of iterations for each Markov chain
      (including warmup).")
  }

  if(class(nIterations) == "numeric" && nIterations < 40) {
    stop("nIterations argument should be greater than 40.")
  }

  if(all(membership != "none") && class(membership) != "numeric") {
    stop("membership should be a numeric vector containing the
      cluster membership. Otherwise, leave as 'none'.")
  }

  if(all(membership != "none") && length(membership) != nObservations) {
    stop("membership should be a numeric vector, where length(membership)
      should equal the number of observations. Otherwise, leave as 'none'.")
  }

  if(all(membership != "none") &&
      all((diff(sort(unique(membership))) == 1) != TRUE) ) {
    stop("Cluster memberships in the membership vector
      are missing a cluster, e.g. 1, 3, 4, 5, 6 is missing cluster 2.")
  }

  if (gmax > nObservations) {
    stop("gmax cannot be larger than nrow(dataset).")
  }

  if (class(initMethod) == "character"){
    if(initMethod != "kmeans" & initMethod != "random" & initMethod != "medoids" & initMethod != "medoids" & initMethod != "clara" & initMethod != "fanny") {
    stop("initMethod should of class character, specifying
      either: kmeans, random, medoids, clara, or fanny.")
    }
  } else if (class(initMethod) != "character") {
    stop("initMethod should of class character, specifying
      either: kmeans, random, medoids, clara, or fanny.")
  }

  if (class(nInitIterations) != "numeric") {
    stop("nInitIterations should be positive integer or zero, specifying
      the number of initialization runs to be considered.")
  }

  if (class(normalize) != "character") {
    stop("normalize should be a string of class character specifying
      if normalization should be performed.")
  }

  if (normalize != "Yes" && normalize != "No") {
    stop("normalize should be a string indicating Yes or No, specifying
       if normalization should be performed.")
  }

  # Check numNodes and calculate the number of cores and initiate cluster
  if(is.na(numNodes) == TRUE) {
    cl <- parallel::makeCluster(parallel::detectCores() - 1)
  } else if (class(numNodes) == "numeric") {
    cl <- parallel::makeCluster(numNodes)
  } else {
    stop("numNodes should be a positive integer indicating the number of
        nodes to be used from the local machine. Else leave as NA, so the
        numNodes will be determined by
        parallel::makeCluster(parallel::detectCores() - 1).")
  }

  # Remove rows with all zeros, if present
  if(length(which(apply(dataset, 1, function(x) all(x == 0)) == TRUE)) != 0) {
    cat("\nDataset row(s)",
      c(which(apply(dataset, 1, function(x) all(x == 0)) == TRUE)),
      "will be removed as this/these contain(s) all zeros")

    if(class(membership) != "character") {
      membership <- membership[- c(which(apply(dataset, 1, function(x)
        all(x == 0)) == TRUE))]
    }

    dataset <- dataset[- c(which(apply(dataset, 1, function(x)
      all(x == 0)) == TRUE)), ]
    nObservations <- nrow(dataset)
  }

  if(class(membership) == "character") {
    membership <- "Not provided"
  }

  # Calculating normalization factors
  if(normalize == "Yes") {
    normFactors <- log(as.vector(edgeR::calcNormFactors(as.matrix(dataset),
      method = "TMM")))
  } else if(normalize == "No") {
    normFactors <- rep(0, dimensionality)
  } else{
    stop("normalize should be 'Yes' or 'No' ")
  }

  # Construct a Stan model
  # Stan model was developed with help of Sanjeena Dang
  # <sdang@math.binghamton.edu>
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
  mplnParallelCode <- function(g) {
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
      "mplnParallelCode",
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
                                     fun = mplnParallelCode,
                                     g = gmin:gmax)
  cat("\nDone parallel code.")
  parallel::stopCluster(cl)

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

  class(RESULTS) <- "mplnParallel"
  return(RESULTS)
  # [END]
}


#' Model-Based Clustering Using MPLN via Non-Parallel Performance
#'
#' Performs clustering using mixtures of multivariate Poisson-log
#' normal (MPLN) distribution and model selection using AIC, AIC3,
#' BIC and ICL. No internal parallelization is performed.
#'
#' @param dataset A dataset of class matrix and type integer such that
#'    rows correspond to observations and columns correspond to variables.
#' @param membership A numeric vector of length nrow(dataset) containing the
#'    cluster membership of each observation. If not available,
#'    leave as "none".
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
#'    "kmeans", "random", "medoids", "clara", or "fanny". Default is "kmeans"
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
#' # Clustering
#' # mplnResults <- mplnNonParallel(dataset = sampleData$dataset,
#' #                      membership = sampleData$trueMembership,
#' #                      gmin = 1,
#' #                      gmax = 2,
#' #                      nChains = 3,
#' #                      nIterations = 1000,
#' #                      initMethod = "kmeans",
#' #                      nInitIterations = 2,
#' #                      normalize = "Yes")
#'
#' @author {Anjali Silva, \email{anjali.silva@uhnresearch.ca}, Sanjeena Dang,
#'          \email{sdang@math.binghamton.edu}. }
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
#' @importFrom edgeR calcNormFactors
#' @import MASS
#' @importFrom mclust unmap
#' @importFrom mclust map
#' @importFrom mvtnorm rmvnorm
#' @import Rcpp
#' @importFrom rstan sampling
#' @importFrom rstan stan_model
#' @import parallel
#' @import stats
#' @importFrom utils tail
#' @importFrom utils write.csv
#'
mplnNonParallel <- function(dataset, membership = "none",
  gmin = 1, gmax = 2,
  nChains = 3, nIterations = 1000,
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

  dimensionality <- ncol(dataset)
  nObservations <- nrow(dataset)

  if(class(gmin) != "numeric" || class(gmax) != "numeric") {
    stop("Class of gmin and gmin should be numeric.")
  }

  if (gmax < gmin) {
    stop("gmax cannot be less than gmin.")
  }

  if (class(nChains) != "numeric") {
    stop("nChains should be a positive integer of class numeric specifying
      the number of Markov chains.")
  }

  if(nChains < 3) {
    cat("nChains is less than 3. Note: recommended minimum number of
      chains is 3.")
  }

  if(class(nIterations) != "numeric") {
    stop("nIterations should be a positive integer of class numeric,
      specifying the number of iterations for each Markov chain
      (including warmup).")
  }

  if(class(nIterations) == "numeric" && nIterations < 40) {
    stop("nIterations argument should be greater than 40.")
  }

  if(all(membership != "none") && class(membership) != "numeric") {
    stop("membership should be a numeric vector containing the
      cluster membership. Otherwise, leave as 'none'.")
  }

  if(all(membership != "none") && length(membership) != nObservations) {
    stop("membership should be a numeric vector, where length(membership)
      should equal the number of observations. Otherwise, leave as 'none'.")
  }

  if(all(membership != "none") &&
      all((diff(sort(unique(membership))) == 1) != TRUE) ) {
    stop("Cluster memberships in the membership vector
      are missing a cluster, e.g. 1, 3, 4, 5, 6 is missing cluster 2.")
  }

  if (gmax > nObservations) {
    stop("gmax cannot be larger than nrow(dataset).")
  }

  if (class(initMethod) == "character"){
    if(initMethod != "kmeans" & initMethod != "random" & initMethod != "medoids" & initMethod != "medoids" & initMethod != "clara" & initMethod != "fanny") {
      stop("initMethod should of class character, specifying
        either: kmeans, random, medoids, clara, or fanny.")
    }
    } else if (class(initMethod) != "character") {
      stop("initMethod should of class character, specifying
        either: kmeans, random, medoids, clara, or fanny.")
  }

  if (class(nInitIterations) != "numeric") {
    stop("nInitIterations should be positive integer or zero, specifying
      the number of initialization runs to be considered.")
  }

  if (class(normalize) != "character") {
    stop("normalize should be a string of class character specifying
      if normalization should be performed.")
  }

  if (normalize != "Yes" && normalize != "No") {
    stop("normalize should be a string indicating Yes or No, specifying
      if normalization should be performed.")
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
  # Stan model was developed with help of Sanjeena Dang
  # <sdang@math.binghamton.edu>
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
  nonParallelRun <- list()

  # Constructing non parallel code
  for (gmodel in seq_along(1:(gmax - gmin + 1))) {

    if(length(1:(gmax - gmin + 1)) == gmax) {
      clustersize <- gmodel
    } else if(length(1:(gmax - gmin + 1)) < gmax) {
      clustersize <- seq(gmin, gmax, 1)[gmodel]
    }

    if(nInitIterations != 0) {
      # cat("\nRunning initialization for G =", clustersize)
      initializeruns <- initializationRun(gmodel = clustersize,
        dataset = dataset,
        init_method = initMethod,
        init_iterations = nInitIterations,
        n_chain = nChains,
        numb_iterations = nIterations,
        initialization = NA,
        normalizefactors = normFactors,
        mod = mod)
      # cat("\nInitialization done for G =", clustersize)
      # cat("\nRunning clustering for G =", clustersize)
      nonParallelRun[[gmodel]] <- mplnCluster(dataset = dataset,
        z = NA,
        G = clustersize,
        nChains = nChains,
        nIterations = nIterations,
        initialization = initializeruns,
        normalizefac = normFactors,
        mod = mod)
      # cat("\nClustering done for G =", clustersize)
    } else if(nInitIterations == 0) {
      # cat("\nNo initialization done for G =", clustersize)
      # cat("\nRunning clustering for G =", clustersize)
      nonParallelRun[[gmodel]] <- mplnCluster(dataset = dataset,
        z = mclust::unmap(stats::kmeans(log(dataset + 1/3),
          clustersize)$cluster),
        G = clustersize,
        nChains = nChains,
        nIterations = nIterations,
        initialization = NA,
        normalizefac = normFactors,
        mod = mod)
      # cat("\nClustering done for G =", clustersize)
    }
  }


  BIC <- ICL <- AIC <- AIC3 <- Djump <- DDSE <- k <- ll <- vector()

  for(g in seq_along(1:(gmax - gmin + 1))) {
    # save the final log-likelihood
    ll[g] <- unlist(utils::tail(nonParallelRun[[g]]$loglikelihood, n = 1))

    k[g] <- calcParameters(numberG = g, dimensionality = dimensionality)

    if (g == max(1:(gmax - gmin + 1))) { # starting model selection
      bic <- BICFunction(ll = ll,
        k = k,
        n = nObservations,
        run = nonParallelRun,
        gmin = gmin,
        gmax = gmax,
        parallel = FALSE)

      icl <- ICLFunction(bIc = bic,
        gmin = gmin,
        gmax = gmax,
        run = nonParallelRun,
        parallel = FALSE)

      aic <- AICFunction(ll = ll,
        k = k,
        run = nonParallelRun,
        gmin = gmin,
        gmax = gmax,
        parallel = FALSE)

      aic3 <- AIC3Function(ll = ll,
        k = k,
        run = nonParallelRun,
        gmin = gmin,
        gmax = gmax,
        parallel = FALSE)
    }
  }

  # for Djump and DDSE
  if((gmax - gmin + 1) > 10 ) {
    # adapted based on HTSCluster package 2.0.8 (25 Oct 2016)
    PMM <- nonParallelRun
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
      all_results = nonParallelRun,
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
      normalization_factors = normFactors,
      gmin = gmin,
      gmax = gmax,
      initalization_method = initMethod,
      all_results = nonParallelRun,
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

  class(RESULTS) <- "mplnNonParallel"
  return(RESULTS)
  # [END]
}


# Log likelihood calculation
calcLikelihood <- function(z, PI, dataset,
  mu_g, G, Sig_g,
  thetaStan, normFactors) {
  nObservations <- nrow(dataset)
  like <- matrix(NA, nrow = nObservations, ncol = G)
  for (g in seq_along(1:G)) {
    for (i in seq_along(1:nObservations)) {
      x <- thetaStan[[g]][i, ]
      dimensionality <- ncol(dataset)
      like[i, g] <- (z[i, g] * (log(PI[g]) +
          t(dataset[i, ]) %*% (x + normFactors) -
          sum(exp(x + normFactors)) -
          sum(lfactorial(dataset[i, ])) -
          dimensionality / 2 * log(2 * pi) -
          1 / 2 * log(det(Sig_g[((g - 1) *
              dimensionality + 1):(g * dimensionality), ])) -
          0.5 * t(x-mu_g[g, ]) %*%
          solve(Sig_g[((g - 1) *
              dimensionality + 1):(g * dimensionality), ]) %*%
          (x - mu_g[g, ])))
    }
  }
  loglike <- sum(rowSums(like))
  return(loglike)
  # Developed by Anjali Silva
}


# Parameter calculation
calcParameters <- function(numberG, dimensionality) {

  muPara <- dimensionality * numberG

  sigmaPara <- (dimensionality * ((dimensionality + 1) / 2)) * numberG
  # because if you have numberG-1 parameters,
  # you can do 1-these to get the last one

  piPara <- numberG - 1

  # total parameters
  paraTotal <- muPara + sigmaPara + piPara

  return(paraTotal)
  # Developed by Anjali Silva
}

# Zvalue calculation
calcZvalue <- function(theta_Stan, dataset, G,
  mu_g, Sig_g, PI,
  normalizefactors) {
  dimensionality <- ncol(dataset)
  nObservations <- nrow(dataset)
  forz <- matrix(NA,ncol = G, nrow = nObservations)

  for (g in seq_along(1:G)) {
    for (i in seq_along(1:nObservations)) {
      x <- theta_Stan[[g]][i, ]
      # for zig calculation (the numerator part)
      forz[i, g] <- PI[g] * exp(t(dataset[i, ]) %*%
          (x + normalizefactors) -
          sum(exp(x + normalizefactors)) -
          sum(lfactorial(dataset[i, ])) -
          dimensionality / 2 * log(2 * pi) - 1 / 2 *
          log(det(Sig_g[((g - 1) *
              dimensionality + 1):(g*dimensionality), ])) -
          0.5 * t(x-mu_g[g, ]) %*%
          solve(Sig_g[((g - 1) *
              dimensionality + 1):(g * dimensionality), ]) %*%
          (x - mu_g[g, ]))
    }

    # check which forz == 0 and rowSums(forz)==0 and which of these
    # have both equalling to 0 (because 0/0 =NaN)
    if (G == 1) {
      errorpossible <- Reduce(intersect, list(which(forz == 0),
        which(rowSums(forz) == 0)))
      zvalue <- forz / rowSums(forz)
      zvalue[errorpossible, ] <- 1
    }else {
      zvalue <- forz / rowSums(forz)
    }
  }

  # check which forz == 0 and rowSums(forz)==0 and which of these
  # have both equalling to 0 (because 0/0 =NaN)
  if (G == 1) {
    errorpossible <- Reduce(intersect,
      list(which(forz == 0),
        which(rowSums(forz) == 0)))
    zvalue <- forz / rowSums(forz)
    zvalue[errorpossible, ] <- 1
  }else {
    zvalue <- forz / rowSums(forz)
  }
  return(zvalue)
  # Developed by Anjali Silva
}


# Calling clustering
callingClustering <- function(data, gmin, gmax, nChains,
  nIterations = NA, initMethod=NA,
  nInitIterations = NA,
  normFactors, model) {
  ptm_inner = proc.time()

  for (gmodel in seq_along(1:(gmax - gmin + 1))) {

    if(length(1:(gmax - gmin + 1)) == gmax) {
      clustersize <- gmodel
    } else if(length(1:(gmax - gmin + 1)) < gmax) {
      clustersize <- seq(gmin, gmax, 1)[gmodel]
    }

    if(nInitIterations != 0) {
      # cat("\nRunning initialization for G =", clustersize)
      initializeruns <- initializationRun(gmodel = clustersize,
        dataset = data,
        init_method = initMethod,
        init_iterations = nInitIterations,
        n_chain = nChains,
        numb_iterations = nIterations,
        initialization = NA,
        normalizefactors = normFactors,
        mod = model)
      # cat("\nInitialization done for G =", clustersize)
      # cat("\nRunning clustering for G =", clustersize)
      allruns <- mplnCluster(dataset = data,
        z = NA,
        G = clustersize,
        nChains = nChains,
        nIterations = nIterations,
        initialization = initializeruns,
        normalizefac = normFactors,
        mod = model)
      # cat("\nClustering done for G =", clustersize)
    } else if(nInitIterations == 0) {
      # cat("\nNo initialization done for G =", clustersize)
      # cat("\nRunning clustering for G =", clustersize)
      allruns <- mplnCluster(dataset = data,
        z = mclust::unmap(stats::kmeans(log(data + 1/3),
          clustersize)$cluster),
        G = clustersize,
        nChains = nChains,
        nIterations = nIterations,
        initialization = NA,
        normalizefac = normFactors,
        mod = model)
      # cat("\nClustering done for G =", clustersize)
    }
  }

  final_inner <- proc.time() - ptm_inner

  RESULTS <- list(gmin = gmin,
    gmax = gmax,
    initalization_method = initMethod,
    all_results = allruns,
    total_time = final_inner)

  class(RESULTS) <- "mplnCallingClustering"
  return(RESULTS)
  # Developed by Anjali Silva
}

# Initialization
initializationRun <- function(gmodel, dataset, init_method,
  init_iterations, n_chain,
  numb_iterations,
  initialization = NA,
  normalizefactors, mod) {

  z <- init_runs <- list()
  logL_init <- vector()
  nObservations <- nrow(dataset)
  dimensionality <- ncol(dataset)

  for(iterations in seq_along(1:init_iterations)) {
    # setting seed, to ensure if multiple iterations are selected by
    # user, then each run will give a different result.
    set.seed(iterations)
    if (init_method == "kmeans" | is.na(init_method)) {
      z[[iterations]] <- mclust::unmap(stats::kmeans(log(dataset + 1/3),
        gmodel)$cluster)
    } else if (init_method == "random") {
      if(gmodel == 1) { # generating z if g=1
        z[[iterations]] <- as.matrix(rep.int(1, times = nObservations),
          ncol = gmodel,
          nrow = nObservations)
      } else { # generating z if g>1
        z_conv = 0
        while(! z_conv) {
          # ensure that dimension of z is same as G (i.e.,
          # if one column contains all 0s, then generate z again)
          z[[iterations]] <- t(stats::rmultinom(nObservations, size = 1,
            prob = rep(1 / gmodel, gmodel)))
          if(length(which(colSums(z[[iterations]]) > 0)) == gmodel) {
            z_conv = 1
          }
        }
      }
    }else if (init_method == "medoids") {
      z[[iterations]] <- mclust::unmap(cluster::pam(log(dataset + 1/3),
        k = gmodel)$cluster)
    }else if (init_method == "clara") {
      z[[iterations]] <- mclust::unmap(cluster::clara(log(dataset + 1/3),
        k = gmodel)$cluster)
    }else if (init_method == "fanny") {
      z[[iterations]] <- mclust::unmap(cluster::fanny(log(dataset + 1/3),
        k = gmodel)$cluster)
    }

    init_runs[[iterations]] <- mplnCluster(dataset = dataset,
      z = z[[iterations]],
      G = gmodel,
      nChains = n_chain,
      nIterations = numb_iterations,
      initialization = "init",
      normalizefac = normalizefactors,
      mod = mod)
    logL_init[iterations] <-
      unlist(utils::tail((init_runs[[iterations]]$loglikelihood), n = 1))
  }

  initialization <- init_runs[[which(logL_init ==
      max(logL_init, na.rm = TRUE))[1]]]
  return(initialization)
  # Developed by Anjali Silva
}


# AIC calculation
AICFunction <- function(ll, k, run, gmin, gmax, parallel = TRUE) {
  AIC <- - 2 * ll + 2 * k
  AICmodel <- seq(gmin, gmax, 1)[grep(min(AIC,na.rm = TRUE), AIC)]
  if(isTRUE(parallel) == "FALSE"){
    # if non parallel run
    AICmodel_labels <- run[[grep(min(AIC,na.rm = TRUE), AIC)]]$clusterlabels
  }else{
    # if parallel run
    AICmodel_labels <- run[[grep(min(AIC,na.rm = TRUE), AIC)]]$all_results$clusterlabels
  }
  AICMessage <- NA

  if (max(AICmodel_labels) != AICmodel) {
    AICmodel <- max(AICmodel_labels)
    AICMessage <- "Spurious or empty cluster resulted."
  }

  AICresults<-list(allAICvalues = AIC,
    AICmodelselected = AICmodel,
    AICmodelselected_labels = AICmodel_labels,
    AICMessage = AICMessage)
  class(AICresults) <- "AIC"
  return(AICresults)
}


# AIC3 calculation
AIC3Function <- function(ll, k, run, gmin, gmax, parallel = TRUE) {
  AIC3 <- - 2 * ll + 3 * k
  AIC3model <- seq(gmin, gmax, 1)[grep(min(AIC3,na.rm = TRUE), AIC3)]
  if(isTRUE(parallel) == "FALSE"){
    # if non parallel run
    AIC3model_labels <- run[[grep(min(AIC3,na.rm = TRUE), AIC3)]]$clusterlabels
  } else{
    # if parallel run
    AIC3model_labels <- run[[grep(min(AIC3,na.rm = TRUE), AIC3)]]$all_results$clusterlabels
  }
  AIC3Message <- NA

  if (max(AIC3model_labels) != AIC3model) {
    AIC3model <- max(AIC3model_labels)
    AIC3Message <- "Spurious or empty cluster resulted."
  }
  AIC3results <- list(allAIC3values = AIC3,
    AIC3modelselected = AIC3model,
    AIC3modelselected_labels = AIC3model_labels,
    AIC3Message = AIC3Message)
  class(AIC3results) <- "AIC3"
  return(AIC3results)
}


# BIC calculation
BICFunction <- function(ll, k, n, run, gmin, gmax, parallel = TRUE) {
  BIC <- - 2 * ll + (k * log(n))
  BICmodel <- seq(gmin, gmax, 1)[grep(min(BIC, na.rm = TRUE), BIC)]
  if(isTRUE(parallel) == "FALSE") {
    # if non parallel run
    BICmodel_labels <- run[[grep(min(BIC, na.rm = TRUE),
      BIC)]]$clusterlabels
  } else {
    # if parallel run
    BICmodel_labels <- run[[grep(min(BIC, na.rm = TRUE),
      BIC)]]$all_results$clusterlabels
  }
  BICMessage <- NA

  if (max(BICmodel_labels) != BICmodel) {
    BICmodel <- max(BICmodel_labels)
    BICMessage <- "Spurious or empty cluster resulted."
  }

  BICresults <- list(allBICvalues = BIC,
    BICmodelselected = BICmodel,
    BICmodelselected_labels = BICmodel_labels,
    BICMessage = BICMessage)
  class(BICresults) <- "BIC"
  return(BICresults)
}

# ICL calculation
ICLFunction <- function(bIc, gmax, gmin, run, parallel = TRUE) {
  ICL <- vector()
  for (g in 1:(gmax - gmin + 1)) {
    if(isTRUE(parallel) == "FALSE") {
      # if non parallel run
      z <- run[[g]]$probaPost
      mapz <- mclust::unmap(run[[g]]$clusterlabels)
    } else {
      # if parallel run
      z <- run[[g]]$all_results$probaPost
      mapz <- mclust::unmap(run[[g]]$all_results$clusterlabels)
    }
    forICL <- function(g){sum(log(z[which(mapz[, g] == 1), g]))}
    ICL[g] <- bIc$allBICvalues[g] + sum(sapply(1:ncol(mapz), forICL))
  }
  ICLmodel <- seq(gmin, gmax, 1)[grep(min(ICL, na.rm = TRUE), ICL)]

  if(isTRUE(parallel) == "FALSE") {
    # if non parallel run
    ICLmodel_labels <- run[[grep(min(ICL, na.rm = TRUE),
      ICL)]]$clusterlabels
  } else {
    # if parallel run
    ICLmodel_labels <- run[[grep(min(ICL, na.rm = TRUE),
      ICL)]]$all_results$clusterlabels
  }

  ICLMessage <- NA

  if (max(ICLmodel_labels) != ICLmodel) {
    ICLmodel <- max(ICLmodel_labels)
    ICLMessage <- "Spurious or empty cluster resulted."
  }

  ICLresults <- list(allICLvalues = ICL,
    ICLmodelselected = ICLmodel,
    ICLmodelselected_labels = ICLmodel_labels,
    ICLMessage = ICLMessage)
  class(ICLresults) <- "ICL"
  return(ICLresults)
  # Developed by Anjali Silva
}


# Clustering function
mplnCluster <- function(dataset, z, G, nChains, nIterations,
  initialization, normalizefac, mod) {

  dimensionality <- ncol(dataset)
  nObservations <- nrow(dataset)

  # for convergence calculation
  norm_mu_outer <- norm_sigma_outer <- vector()
  median_mu_outer <- median_sigma_outer <- list()
  # for saving mu and sigma values
  mu_all_outer <- sigma_all_outer <- list()
  obs <- PI <- logL <- vector()
  it_outer <- 2 # the starting value of interation for outer loop
  conv_outer <- 0


  if (all(is.na(initialization)) == TRUE || all(initialization == "init")) {
    # mean for both t and normal distribution
    mu_all_outer[[1]] <- mu_g <- matrix(log(mean(dataset)),
      ncol = dimensionality, nrow = G)
    # sig for sigma of t distribtuion
    sigma_all_outer[[1]] <- Sig_g <- do.call("rbind",
      rep(list(cov(log(dataset + 1)) * dimensionality), G))
  } else {
    mu_all_outer[[1]] <- mu_g <- initialization$finalmu
    sigma_all_outer[[1]] <- Sig_g <- initialization$finalsigma
    z <- initialization$probaPost
  }

  while(! conv_outer) {
    for(g in seq_along(1:G)) {
      obs[g] <- sum(z[, g]) # number of observations in each group
      PI[g] <- obs[g] / nObservations  # obtain probability of each group
    }

    theta_Stan <- E_theta2 <- list()
    rstan_results <- stanRun(model = mod,
      gmin = 1,
      gmax = G,
      dataset = dataset,
      mu_all_outer = mu_all_outer,
      it_outer = it_outer,
      sigma_all_outer = sigma_all_outer,
      numb_iterations = nIterations,
      n_chain = nChains,
      normalizefacs = normalizefac)

    fit <- rstan_results$fitrstan
    nIterations <- rstan_results$numb_iterations

    for (g in seq_along(1:G)) {
      tt <- as.matrix(fit[[g]])
      theta_Stan[[g]] <- matrix(NA, nrow = nObservations,
        ncol = dimensionality)
      E_theta2[[g]] <- list()

      for (i in seq_along(1:nObservations)) {
        zz <- c(1:(dimensionality - 1)) * nObservations + i
        theta_mat <- tt[, c(i, zz)]
        theta_Stan[[g]][i, ] <- colMeans(theta_mat)
        E_theta2[[g]][[i]] <- z[i, g] * t(tt[, c(i,zz)]) %*% tt[, c(i,zz)]/
          ((0.5 * nIterations) * nChains)
      }

      mu_g[g, ] <- colSums(z[, g] * theta_Stan[[g]]) / sum(z[, g])
      Sig_g[((g - 1) * dimensionality + 1):(g * dimensionality), ] <-
        Reduce("+", E_theta2[[g]]) / sum(z[, g]) - mu_g[g, ] %*% t(mu_g[g, ])
    }

    mu_all_outer[[it_outer]] <- mu_g
    sigma_all_outer[[it_outer]] <- Sig_g

    logL[it_outer] <- calcLikelihood(z = z,
      PI = PI,
      dataset = dataset,
      mu_g = mu_all_outer[[it_outer]],
      G = G,
      Sig_g = sigma_all_outer[[it_outer]],
      thetaStan = theta_Stan,
      normFactors = normalizefac)

    # convergence of outer loop
    norm_mu_outer[it_outer] <- norm((mu_all_outer[[it_outer]] -
        mu_all_outer[[it_outer - 1]]),
      type = "F")
    norm_sigma_outer[it_outer] <- norm(sigma_all_outer[[it_outer]] -
        sigma_all_outer[[it_outer - 1]],
      type="F")
    median_mu_outer[[it_outer]] <- median(norm_mu_outer, na.rm = TRUE)
    median_sigma_outer[[it_outer]] <- median(norm_sigma_outer, na.rm = TRUE)
    # par(mfrow=c(1,2))
    # plot(norm_mu_outer, main=paste0("Norm outer mean, G=", G),
    # type="l", ylab="median(norm_mu_outer)", xlab="iterations")
    # plot(norm_sigma_outer, main=paste0("Norm outer sigma, G=", G),
    # type="l", ylab="median(norm_sigma_outer)", xlab="iterations")

    threshold_outer <- 2
    if(it_outer > (threshold_outer + 1)) {

      # cat("\nMedian difference of mean and sigma in outer loop respectively ",
      # c(abs(median_mu_outer[[it_outer-threshold_outer]]-
      # median_mu_outer[[it_outer]])))
      if( ( (abs(median_mu_outer[[it_outer - threshold_outer]] -
          median_mu_outer[[it_outer]]) < 5) &&
          (abs(median_sigma_outer[[it_outer - threshold_outer]] -
              median_sigma_outer[[it_outer]]) < 5) ) || it_outer > 100) {
        # cat("\nConvergence of mu and sigma at outer loop
        # iteration ", it_outer)
        # take out absolute value
        programclust <- vector()
        programclust <- mclust::map(z)

        # checking for spurious clusters and getting rid of them
        # keep<-as.numeric(names(which(table(programclust)>5)))
        # if ( (length(keep) !=length(unique(programclust))) &&
        # (length(keep) !=0) ){
        #  z<-as.matrix(z[,keep])
        #  z<-z/rowSums(z)
        #  programclust<-map(z)
        #}

        # checking for empty clusters
        J <- 1:ncol(z)
        K <- as.logical(match(J, sort(unique(programclust)), nomatch = 0))
        if(length(J[! K]) > 0) { # J[!K] tells which are empty clusters
          z <- z[, - J[! K]]
          programclust <- mclust::map(z)
        }
        conv_outer <- 1
      }
    }

    # if running for initialization, need to stop after 10 iterations
    if(it_outer == 10 && all(is.na(initialization) != TRUE)) {
      if(all(initialization == "init")) {
        programclust <- vector()
        programclust <- mclust::map(z)
        conv_outer <- 1
      }
    }

    if(conv_outer != 1) { # only update until convergence, not after
      z <- calcZvalue(theta_Stan = theta_Stan,
        dataset = dataset,
        G = G,
        mu_g = mu_g,
        Sig_g = Sig_g,
        PI = PI,
        normalizefactors = normalizefac)
      it_outer <- it_outer + 1 # updating outer loop iteration
      nIterations <- nIterations + 10
    }
  } # end of outer loop

  results <- list(finalmu = mu_all_outer[[it_outer]] +
      matrix(rep(normalizefac,
        nrow(mu_all_outer[[it_outer]])),
        byrow=TRUE, ncol =
          ncol(mu_all_outer[[it_outer]])),
    finalsigma = sigma_all_outer[[it_outer]],
    allmu = lapply(mu_all_outer, function(x)
      (x + matrix(rep(normalizefac,
        nrow(mu_all_outer[[it_outer]])),
        byrow = TRUE,
        ncol = ncol(mu_all_outer[[it_outer]])))),
    allsigma = sigma_all_outer,
    clusterlabels = programclust,
    iterations = it_outer,
    proportion = PI,
    loglikelihood = logL,
    probaPost = z,
    stanresults = fit)

  class(results) <- "MPLNcluster"
  return(results)
  # Developed by Anjali Silva
}

# Calculate and remove rows with zeros
removeZeroCounts <- function(dataset, output = FALSE) {

  zeroSUMrows <- which(rowSums(dataset) == 0)
  nbrZeroSUMrows <- length(zeroSUMrows)

  if (output) {
    write.csv(rownames(dataset[zeroSUMrows, ]),
      file=paste0("zeroSUMrows_",
        format(Sys.time(), "%d%b%Y"),".csv") )
  }

  if (nbrZeroSUMrows > 0) {
    dataset <- dataset[- zeroSUMrows, ]
    cat(paste(nbrZeroSUMrows,
      "row(s) removed from the dataset because sum = 0",'\n'))
  }


  RESULTS <- list(dataset = dataset)
  class(RESULTS) <- "MPLN_ZerosRemoved"
  return(RESULTS)
  # Developed by Anjali Silva
}

# Stan sampling
stanRun <- function(model, gmin, gmax, dataset,
  mu_all_outer, it_outer,
  sigma_all_outer, numb_iterations,
  n_chain = n_chain, normalizefacs) {

  fitrstan <- list()
  dimensionality <- ncol(dataset)

  for (g in seq_along(gmin:gmax)) {
    data1 = list(d = ncol(dataset),
      N = nrow(dataset),
      y = dataset,
      mu = mu_all_outer[[it_outer-1]][g, ],
      Sigma = sigma_all_outer[[it_outer-1]][((g - 1) *
          dimensionality + 1):(g*dimensionality), ],
      normfactors = as.vector(normalizefacs))
    stanproceed <- 0
    try = 1

    while (! stanproceed) {

      # cat("\nRstan generating sample at outer iteration",
      # it_outer, "for g: ",g , "try: ", try)
      # cat("\nNumber of iterations is", numb_iterations, "\n")
      fitrstan[[g]] <- rstan::sampling(object = model,
        data = data1,
        iter = numb_iterations,
        chains = n_chain,
        verbose = FALSE,
        refresh = -1)

      if (all(rstan::summary(fitrstan[[g]])$summary[, "Rhat"] < 1.1) == TRUE &&
          all(rstan::summary(fitrstan[[g]])$summary[, "n_eff"] > 100) == TRUE) {
        stanproceed <- 1
      } else if(all(rstan::summary(fitrstan[[g]])$summary[,"Rhat"] < 1.1) != TRUE ||           all(rstan::summary(fitrstan[[g]])$summary[, "n_eff"] > 100) != TRUE) {
        if(try == 10) { # stop after 10 attempts
          stanproceed = 1
        }
        numb_iterations = numb_iterations + 100
        try = try + 1
      }
    }
  }

  results <- list(fitrstan = fitrstan,
    numb_iterations = numb_iterations)
  class(results) <- "RStan"
  return(results)
  # Developed by Anjali Silva
}
