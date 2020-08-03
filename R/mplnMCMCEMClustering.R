#' Clustering Using MPLN With MCMC-EM Via Parallel Performance
#'
#' Performs clustering using mixtures of multivariate Poisson-log
#' normal (MPLN) distribution with Markov chain Monte Carlo
#' expectation-maximization algorithm (MCMC-EM) for parameter
#' estimation. Coarse grain parallelization is employed, such that
#' when a range of components/clusters (g = 1,...,G) are considered, each
#' G is run on a different processor. This can be performed because
#' each component/cluster size is independent from another. All
#' components/clusters in the range to be tested have been parallelized
#' to run on a seperate core using the *parallel* R package. The number of
#' nodes to be used for clustering can be specified or calculated using
#' *parallel::detectCores() - 1*. Model selection is performed using
#' AIC, AIC3, BIC and ICL.
#'
#' @param dataset A dataset of class matrix and type integer such that
#'    rows correspond to observations and columns correspond to variables.
#'    The dataset have dimensions n x d, where n is the total number of
#'    observations and d is the dimensionality. If rowSums are zero,
#'    these rows will be removed prior to cluster analysis.
#' @param membership A numeric vector of length nrow(dataset) containing the
#'    cluster membership of each observation. If not available,
#'    leave as "none".
#' @param gmin A positive integer specifying the minimum number of components
#'    to be considered in the clustering run.
#' @param gmax A positive integer, >= gmin, specifying the maximum number of
#'    components to be considered in the clustering run.
#' @param nChains A positive integer specifying the number of Markov chains.
#'    Default is 3, the recommended minimum number.
#' @param nIterations A positive integer specifying the number of iterations
#'    for each MCMC chain (including warmup). The value should be greater than
#'    40. The upper limit will depend on size of dataset.
#' @param initMethod An algorithm for initialization. Current options are
#'    "kmeans", "random", "medoids", "clara", or "fanny". Default is "kmeans".
#' @param nInitIterations A positive integer or zero, specifying the number
#'    of initialization runs to be performed. This many runs, each with 10
#'    iterations, will be performed via MPLNClust and values from the run with
#'    highest log-likelihood will be used as initialization values. Default is 0.
#' @param normalize A string with options "Yes" or "No" specifying
#'     if normalization should be performed. Currently, normalization factors
#'     are calculated using TMM method of edgeR package. Default is "Yes".
#' @param numNodes A positive integer indicating the number of nodes to be
#'     used from the local machine to run the clustering algorithm. Else
#'     leave as NA, so default will be detected as
#'     parallel::makeCluster(parallel::detectCores() - 1).
#'
#' @return Returns an S3 object of class mplnMCMCParallel with results.
#' \itemize{
#'   \item dataset - The input dataset on which clustering is performed.
#'   \item dimensionality - Dimensionality of the input dataset.
#'   \item normalizationFactors - A vector of normalization factors used
#'      for input dataset.
#'   \item gmin - Minimum number of components considered in the clustering
#'      run.
#'   \item gmax - Maximum number of components considered in the clustering
#'      run.
#'   \item initalizationMethod - Method used for initialization.
#'   \item allResults - A list with all results.
#'   \item logLikelihood - A vector with value of final log-likelihoods for
#'      each cluster size.
#'   \item numbParameters - A vector with number of parameters for each
#'      cluster size.
#'   \item trueLabels - The vector of true labels, if provided by user.
#'   \item ICLresults - A list with all ICL model selection results.
#'   \item BICresults - A list with all BIC model selection results.
#'   \item AICresults - A list with all AIC model selection results.
#'   \item AIC3results - A list with all AIC3 model selection results.
#'   \item slopeHeuristics - If more than 10 models are considered, slope heuristic
#'      results.
#'   \item DjumpModelSelected - If more than 10 models are considered, slope heuristic
#'      results.
#'   \item DDSEModelSelected - If more than 10 models are considered, slope heuristic
#'      results.
#'   \item totalTime - Total time used for clustering and model selection.
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
#' # sampleData <- MPLNClust::mplnDataGenerator(nObservations = 100,
#' #                                 dimensionality = 6,
#' #                                 mixingProportions = c(0.79, 0.21),
#' #                                 mu = rbind(trueMu1, trueMu2),
#' #                                 sigma = rbind(trueSigma1, trueSigma2),
#' #                                 produceImage = "No")
#'
#' # Clustering
#' # mplnResults <- MPLNClust::mplnMCMCParallel(dataset = sampleData$dataset,
#' #                                 membership = sampleData$trueMembership,
#' #                                 gmin = 1,
#' #                                 gmax = 2,
#' #                                 nChains = 3,
#' #                                 nIterations = 1000,
#' #                                 initMethod = "kmeans",
#' #                                 nInitIterations = 2,
#' #                                 normalize = "Yes")
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
mplnMCMCParallel <- function(dataset,
                             membership = "none",
                             gmin = 1,
                             gmax = 2,
                             nChains = 3,
                             nIterations = 1000,
                             initMethod = "kmeans",
                             nInitIterations = 0,
                             normalize = "Yes",
                             numNodes = NA) {

  # Keeping track of time
  ptm <- proc.time()

  # Performing checks of user input
  if (typeof(dataset) != "double" & typeof(dataset) != "integer") {
    stop("Dataset type needs to be integer.")
  }

  if (any((dataset %% 1 == 0) == FALSE)) {
    stop("Dataset should be a matrix of counts.")
  }

  if (class(dataset) != "matrix") {
    stop("Dataset needs to be a matrix.")
  }

  if (any(colSums(dataset) <= 0)) {
    stop("Column sums cannot be less than or equal to 0. Double check dataset.")
  }

  dimensionality <- ncol(dataset)
  nObservations <- nrow(dataset)

  if(class(gmin) != "numeric" || class(gmax) != "numeric") {
    stop("Class of gmin and gmin should be numeric.")
  }

  if (gmax < gmin) {
    stop("gmax cannot be less than gmin.")
  }

  if (gmax > nObservations) {
    stop("gmax cannot be larger than nrow(dataset).")
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

  if(all(membership != "none") &&
      all((diff(sort(unique(membership))) == 1) != TRUE) ) {
    stop("Cluster memberships in the membership vector
      are missing a cluster, e.g. 1, 3, 4, 5, 6 is missing cluster 2.")
  }

  if(all(membership != "none") && length(membership) != nObservations) {
    stop("membership should be a numeric vector, where length(membership)
      should equal the number of observations. Otherwise, leave as 'none'.")
  }

  # Remove rows with only zeros, if present
  removezeros <- removeZeroCounts(dataset = dataset, membership = membership)
  dataset <- removezeros$dataset
  membership <- removezeros$membership
  dimensionality <- ncol(dataset)
  nObservations <- nrow(dataset)

  if (class(initMethod) == "character"){
    initMethodsUsed <- c("kmeans", "random", "medoids", "clara", "fanny")
    if(all((initMethod == initMethodsUsed) == FALSE)) {
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
    varlist = c(
      "nInitIterations",
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

  parallelRun <- list()
  cat("\nRunning parallel code now.")
  parallelRun <- parallel::clusterMap(cl = cl,
    fun = mplnParallelCode,
    g = gmin:gmax)
  cat("\nDone parallel code.")
  parallel::stopCluster(cl)

  BIC <- ICL <- AIC <- AIC3 <- Djump <- DDSE <- nParameters <- logLikelihood <- vector()

  for(g in seq_along(1:(gmax - gmin + 1))) {
    # save the final log-likelihood
    logLikelihood[g] <- unlist(tail(parallelRun[[g]]$allResults$logLikelihood,
      n = 1))

    if(length(1:(gmax - gmin + 1)) == gmax) {
      clustersize <- g
    } else if(length(1:(gmax - gmin + 1)) < gmax) {
      clustersize <- seq(gmin, gmax, 1)[g]
    }

    nParameters[g] <- calcParameters(numberG = clustersize,
      dimensionality = dimensionality)

    if (g == max(1:(gmax - gmin + 1))) { # starting model selection
      bic <- BICFunction(logLikelihood = logLikelihood,
                         nParameters = nParameters,
                         nObservations = nObservations,
                         clusterRunOutput = parallelRun,
                         gmin = gmin,
                         gmax = gmax,
                         parallel = TRUE)

      icl <- ICLFunction(logLikelihood = logLikelihood,
                         nParameters = nParameters,
                         nObservations = nObservations,
                         gmin = gmin,
                         gmax = gmax,
                         clusterRunOutput = parallelRun,
                         parallel = TRUE)

      aic <- AICFunction(logLikelihood = logLikelihood,
                          nParameters = nParameters,
                          clusterRunOutput = parallelRun,
                          gmin = gmin,
                          gmax = gmax,
                          parallel = TRUE)

      aic3 <- AIC3Function(logLikelihood = logLikelihood,
                          nParameters = nParameters,
                          clusterRunOutput = parallelRun,
                          gmin = gmin,
                          gmax = gmax,
                          parallel = TRUE)
    }
  }

  # for Djump and DDSE
  if((gmax - gmin + 1) > 10 ) {
    # adapted based on HTSCluster package 2.0.8 (25 Oct 2016)
    Kchoice <- gmin:gmax # all the cluster/component sizes considered
    message("Note: diagnostic plots for results corresponding
      to model selection via slope heuristics (Djump and DDSE)
      should be examined to ensure that sufficiently complex
      models have been considered.")
    mat <- cbind(Kchoice, nParameters / nObservations, nParameters / nObservations, - logLikelihood)
    ResCapushe <- capushe ::capushe(mat, nObservations)
    DDSEmodel<- ResCapushe@DDSE@model
    Djumpmodel<- ResCapushe@Djump@model
    final <- proc.time() - ptm

    RESULTS <- list(dataset = dataset,
                    dimensionality = dimensionality,
                    normalizationFactors = normFactors,
                    gmin = gmin,
                    gmax = gmax,
                    initalizationMethod = initMethod,
                    allResults = parallelRun,
                    logLikelihood = logLikelihood,
                    numbParameters = nParameters,
                    trueLabels = membership,
                    ICLresults = icl,
                    BICresults = bic,
                    AICresults = aic,
                    AIC3results = aic3,
                    slopeHeuristics = ResCapushe,
                    DjumpModelSelected = ResCapushe@Djump@model,
                    DDSEModelSelected = ResCapushe@DDSE@model,
                    totalTime = final)

  } else {
    final <- proc.time() - ptm
    RESULTS <- list(dataset = dataset,
                    dimensionality = dimensionality,
                    normalizationFactors = normFactors,
                    gmin = gmin,
                    gmax = gmax,
                    initalizationMethod = initMethod,
                    allResults = parallelRun,
                    logLikelihood = logLikelihood,
                    numbParameters = nParameters,
                    trueLabels = membership,
                    ICLresults = icl,
                    BICresults = bic,
                    AICresults = aic,
                    AIC3results = aic3,
                    slopeHeuristics = "Not used",
                    totalTime = final)
  }

  class(RESULTS) <- "mplnMCMCParallel"
  return(RESULTS)
}



#' Clustering Using MPLN  With MCMC-EM Via Non-Parallel Performance
#'
#' Performs clustering using mixtures of multivariate Poisson-log
#' normal (MPLN) distribution with Markov chain Monte Carlo
#' expectation-maximization algorithm (MCMC-EM) for parameter
#' estimation. No internal parallelization, thus code is run in
#' serial. Model selection is performed using AIC, AIC3, BIC
#' and ICL.
#'
#' @param dataset A dataset of class matrix and type integer such that
#'    rows correspond to observations and columns correspond to variables.
#'    The dataset have dimensions n x d, where n is the total number of
#'    observations and d is the dimensionality. If rowSums are zero,
#'    these rows will be removed prior to cluster analysis.
#' @param membership A numeric vector of length nrow(dataset) containing the
#'    cluster membership of each observation. If not available,
#'    leave as "none".
#' @param gmin A positive integer specifying the minimum number of components
#'    to be considered in the clustering run.
#' @param gmax A positive integer, >= gmin, specifying the maximum number of
#'    components to be considered in the clustering run.
#' @param nChains A positive integer specifying the number of Markov chains.
#'    Default is 3, the recommended minimum number.
#' @param nIterations A positive integer specifying the number of iterations
#'    for each chain (including warmup). The value should be greater than 40.
#'    The upper limit will depend on size of dataset.
#' @param initMethod An algorithm for initialization. Current options are
#'    "kmeans", "random", "medoids", "clara", or "fanny". Default is "kmeans"
#' @param nInitIterations A positive integer or zero, specifying the number
#'    of initialization runs to be performed. This many runs, each with 10
#'    iterations, will be performed via MPLNClust and values from the run with
#'    highest log-likelihood will be used as initialization values. Default is 0.
#' @param normalize A string with options "Yes" or "No" specifying
#'     if normalization should be performed. Currently, normalization factors
#'     are calculated using TMM method of edgeR package. Default is "Yes".
#'
#' @return Returns an S3 object of class mplnMCMCNonParallel with results.
#' \itemize{
#'   \item dataset - The input dataset on which clustering is performed.
#'   \item dimensionality - Dimensionality of the input dataset.
#'   \item normalizationFactors - A vector of normalization factors used
#'      for input dataset.
#'   \item gmin - Minimum number of components considered in the clustering
#'      run
#'   \item gmax - Maximum number of components considered in the clustering
#'      run
#'   \item initalizationMethod - Method used for initialization.
#'   \item allResults - A list with all results.
#'   \item logLikelihood - A vector with value of final log-likelihoods for
#'      each cluster size.
#'   \item numbParameters - A vector with number of parameters for each
#'      cluster size.
#'   \item trueLabels - The vector of true labels, if provided by user.
#'   \item ICLresults - A list with all ICL model selection results.
#'   \item BICresults - A list with all BIC model selection results.
#'   \item AICresults - A list with all AIC model selection results.
#'   \item AIC3results - A list with all AIC3 model selection results.
#'   \item slopeHeuristics - If more than 10 models are considered, slope heuristic
#'      results.
#'   \item DjumpModelSelected - If more than 10 models are considered, slope heuristic
#'      results.
#'   \item DDSEModelSelected - If more than 10 models are considered, slope heuristic
#'      results.
#'   \item totalTime - Total time used for clustering and model selection.
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
#' # mplnResults <- mplnMCMCNonParallel(dataset = sampleData$dataset,
#' #                                    membership = sampleData$trueMembership,
#' #                                    gmin = 1,
#' #                                    gmax = 2,
#' #                                    nChains = 3,
#' #                                    nIterations = 1000,
#' #                                    initMethod = "kmeans",
#' #                                    nInitIterations = 2,
#' #                                    normalize = "Yes")
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
mplnMCMCNonParallel <- function(dataset,
                                membership = "none",
                                gmin = 1,
                                gmax = 2,
                                nChains = 3,
                                nIterations = 1000,
                                initMethod = "kmeans",
                                nInitIterations = 0,
                                normalize = "Yes") {
  ptm <- proc.time()

  # Performing checks
  if (typeof(dataset) != "double" & typeof(dataset) != "integer") {
    stop("Dataset type needs to be integer.")
  }

  if (any((dataset %% 1 == 0) == FALSE)) {
    stop("Dataset should be a matrix of counts.")
  }

  if (class(dataset) != "matrix") {
    stop("Dataset needs to be a matrix.")
  }

  if (any(colSums(dataset) <= 0)) {
    stop("Column sums cannot be less than or equal to 0. Double check dataset.")
  }

  dimensionality <- ncol(dataset)
  nObservations <- nrow(dataset)

  if(class(gmin) != "numeric" || class(gmax) != "numeric") {
    stop("Class of gmin and gmin should be numeric.")
  }

  if (gmax < gmin) {
    stop("gmax cannot be less than gmin.")
  }

  if (gmax > nObservations) {
    stop("gmax cannot be larger than nrow(dataset).")
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

  if(all(membership != "none") &&
      all((diff(sort(unique(membership))) == 1) != TRUE) ) {
    stop("Cluster memberships in the membership vector
      are missing a cluster, e.g. 1, 3, 4, 5, 6 is missing cluster 2.")
  }

  if(all(membership != "none") && length(membership) != nObservations) {
    stop("membership should be a numeric vector, where length(membership)
      should equal the number of observations. Otherwise, leave as 'none'.")
  }

  # Remove rows with zeros, if present
  removezeros <- removeZeroCounts(dataset = dataset,
                                  membership = membership)
  dataset <- removezeros$dataset
  membership <- removezeros$membership
  dimensionality <- ncol(dataset)
  nObservations <- nrow(dataset)

  if (class(initMethod) == "character"){
    initMethodsUsed <- c("kmeans", "random", "medoids", "clara", "fanny")
    if(all((initMethod == initMethodsUsed) == FALSE)) {
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
  } else {
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

  mod <- rstan::stan_model(model_code = stancode,
    verbose = FALSE)
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
                                          initMethod = initMethod,
                                          initIterations = nInitIterations,
                                          nChain = nChains,
                                          numbIterations = nIterations,
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
                                              z = mclust::unmap(stats::kmeans(log(dataset + 1 / 3),
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


  BIC <- ICL <- AIC <- AIC3 <- Djump <- DDSE <- nParameters <- logLikelihood <- vector()

  for(g in seq_along(1:(gmax - gmin + 1))) {
    # save the final log-likelihood
    logLikelihood[g] <- unlist(utils::tail(nonParallelRun[[g]]$logLikelihood, n = 1))

    if(length(1:(gmax - gmin + 1)) == gmax) {
      clustersize <- g
    } else if(length(1:(gmax - gmin + 1)) < gmax) {
      clustersize <- seq(gmin, gmax, 1)[g]
    }

    nParameters[g] <- calcParameters(numberG = clustersize, dimensionality = dimensionality)

    if (g == max(1:(gmax - gmin + 1))) { # starting model selection
      bic <- BICFunction(logLikelihood = logLikelihood,
                        nParameters = nParameters,
                        nObservations = nObservations,
                        clusterRunOutput = nonParallelRun,
                        gmin = gmin,
                        gmax = gmax,
                        parallel = FALSE)

      icl <- ICLFunction(logLikelihood = logLikelihood,
                        nParameters = nParameters,
                        nObservations = nObservations,
                        gmin = gmin,
                        gmax = gmax,
                        clusterRunOutput = nonParallelRun,
                        parallel = FALSE)

      aic <- AICFunction(logLikelihood = logLikelihood,
                        nParameters = nParameters,
                        clusterRunOutput = nonParallelRun,
                        gmin = gmin,
                        gmax = gmax,
                        parallel = FALSE)

      aic3 <- AIC3Function(logLikelihood = logLikelihood,
                          nParameters = nParameters,
                          clusterRunOutput = nonParallelRun,
                          gmin = gmin,
                          gmax = gmax,
                          parallel = FALSE)
    }
  }

  # for Djump and DDSE
  if((gmax - gmin + 1) > 10 ) {
    # adapted based on HTSCluster package 2.0.8 (25 Oct 2016)
    Kchoice <- gmin:gmax # all the cluster/component sizes considered
    message("Note: diagnostic plots for results corresponding
      to model selection via slope heuristics (Djump and DDSE)
      should be examined to ensure that sufficiently complex
      models have been considered.")
    mat <- cbind(Kchoice, nParameters / nObservations, nParameters / nObservations, - logLikelihood)
    ResCapushe <- capushe ::capushe(mat, nObservations)
    DDSEmodel<- ResCapushe@DDSE@model
    Djumpmodel<- ResCapushe@Djump@model
    final <- proc.time() - ptm

    RESULTS <- list(dataset = dataset,
                    dimensionality = dimensionality,
                    normalizationFactors = normFactors,
                    gmin = gmin,
                    gmax = gmax,
                    initalizationMethod = initMethod,
                    allResults = nonParallelRun,
                    logLikelihood = logLikelihood,
                    numbParameters = nParameters,
                    trueLabels = membership,
                    ICLresults = icl,
                    BICresults = bic,
                    AICresults = aic,
                    AIC3results = aic3,
                    slopeHeuristics = ResCapushe,
                    DjumpModelSelected = ResCapushe@Djump@model,
                    DDSEModelSelected = ResCapushe@DDSE@model,
                    totalTime = final)

  } else {# end of Djump and DDSE
    final <- proc.time() - ptm
    RESULTS <- list(dataset = dataset,
                    dimensionality = dimensionality,
                    normalizationFactors = normFactors,
                    gmin = gmin,
                    gmax = gmax,
                    initalizationMethod = initMethod,
                    allResults = nonParallelRun,
                    logLikelihood = logLikelihood,
                    numbParameters = nParameters,
                    trueLabels = membership,
                    ICLresults = icl,
                    BICresults = bic,
                    AICresults = aic,
                    AIC3results = aic3,
                    slopeHeuristics = "Not used",
                    totalTime = final)
  }

  class(RESULTS) <- "mplnMCMCNonParallel"
  return(RESULTS)
}



# Log likelihood calculation
calcLikelihood <- function(z,
                           PI,
                           dataset,
                           muG,
                           G,
                           sigmaG,
                           thetaStan,
                           normFactors) {

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
          0.5 * log(det(sigmaG[((g - 1) *
              dimensionality + 1):(g * dimensionality), ])) -
          0.5 * t(x - muG[g, ]) %*%
          solve(sigmaG[((g - 1) *
              dimensionality + 1):(g * dimensionality), ]) %*%
          (x - muG[g, ])))
    }
  }
  loglike <- sum(rowSums(like))
  return(loglike)
  # Developed by Anjali Silva
}



# Parameter calculation
calcParameters <- function(numberG,
                           dimensionality) {

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
calcZvalue <- function(thetaStan,
                       dataset,
                       G,
                       muG,
                       sigmaG,
                       PI,
                       normalizefactors) {

  dimensionality <- ncol(dataset)
  nObservations <- nrow(dataset)
  forz <- matrix(NA,ncol = G, nrow = nObservations)

  for (g in seq_along(1:G)) {
    for (i in seq_along(1:nObservations)) {
      x <- thetaStan[[g]][i, ]
      # for zig calculation (the numerator part)
      forz[i, g] <- PI[g] * exp(t(dataset[i, ]) %*%
          (x + normalizefactors) -
          sum(exp(x + normalizefactors)) -
          sum(lfactorial(dataset[i, ])) -
          dimensionality / 2 * log(2 * pi) - 0.5 *
          log(det(sigmaG[((g - 1) *
              dimensionality + 1):(g*dimensionality), ])) -
          0.5 * t(x - muG[g, ]) %*%
          solve(sigmaG[((g - 1) *
              dimensionality + 1):(g * dimensionality), ]) %*%
          (x - muG[g, ]))
    }

    # check which forz == 0 and rowSums(forz)==0 and which of these
    # have both equalling to 0 (because 0/0 =NaN)
    if (G == 1) {
      errorpossible <- Reduce(intersect,
        list(which(forz == 0),
          which(rowSums(forz) == 0)))
      zvalue <- forz / rowSums(forz)
      zvalue[errorpossible, ] <- 1
    } else {
      zvalue <- forz / rowSums(forz)
    }
  }

  return(zvalue)
  # Developed by Anjali Silva
}



# Calling clustering
callingClustering <- function(data,
                              gmin,
                              gmax,
                              nChains,
                              nIterations = NA,
                              initMethod = NA,
                              nInitIterations = NA,
                              normFactors,
                              model) {

  ptmInner = proc.time()

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
                                          initMethod = initMethod,
                                          initIterations = nInitIterations,
                                          nChain = nChains,
                                          numbIterations = nIterations,
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
                             z = mclust::unmap(stats::kmeans(log(data + 1 / 3),
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

  finalInner <- proc.time() - ptmInner

  RESULTS <- list(gmin = gmin,
                  gmax = gmax,
                  initalizationMethod = initMethod,
                  allResults = allruns,
                  totalTime = finalInner)

  class(RESULTS) <- "mplnCallingClustering"
  return(RESULTS)
  # Developed by Anjali Silva
}



# Initialization
initializationRun <- function(gmodel,
                              dataset,
                              initMethod,
                              initIterations,
                              nChain,
                              numbIterations,
                              initialization = NA,
                              normalizefactors,
                              mod) {

  z <- initRuns <- list()
  logLinit <- vector()
  nObservations <- nrow(dataset)
  dimensionality <- ncol(dataset)

  # Internal function for random initialization
  randomInitfunction <- function(gmodel, nObservations) {
    if(gmodel == 1) { # generating z if g=1
      z <- as.matrix(rep.int(1, times = nObservations),
                          ncol = gmodel,
                          nrow = nObservations)
    } else { # generating z if g>1
      zValueConv <- 0
      while(! zValueConv) {
        # ensure that dimension of z is same as G (i.e.,
        # if one column contains all 0s, then generate z again)
        z <- t(stats::rmultinom(nObservations, size = 1,
                                     prob = rep(1 / gmodel, gmodel)))
        if(length(which(colSums(z) > 0)) == gmodel) {
          zValueConv <- 1
        }
      }
    }
    return(z)
  }


  for(iterations in seq_along(1:initIterations)) {
    # setting seed, to ensure if multiple iterations are selected by
    # user, then each run will give a different result.
    set.seed(iterations)
    if (initMethod == "kmeans" | is.na(initMethod)) {
      z[[iterations]] <- mclust::unmap(stats::kmeans(log(dataset + 1 / 3),
        gmodel)$cluster)
    } else if (initMethod == "random") {
      z[[iterations]] <- randomInitfunction(gmodel = gmodel, nObservations = nObservations)

    } else if (initMethod == "medoids") {
      z[[iterations]] <- mclust::unmap(cluster::pam(log(dataset + 1 / 3),
        k = gmodel)$cluster)
      # if z generated has less columns than gmodel, then use random initialization
      if(ncol(z[[iterations]]) < gmodel) {
        z[[iterations]] <- randomInitfunction(gmodel = gmodel, nObservations = nObservations)
      }

    } else if (initMethod == "clara") {
      z[[iterations]] <- mclust::unmap(cluster::clara(log(dataset + 1 / 3),
        k = gmodel)$cluster)
      # if z generated has less columns than gmodel, then use random initialization
      if(ncol(z[[iterations]]) < gmodel) {
        z[[iterations]] <- randomInitfunction(gmodel = gmodel, nObservations = nObservations)
      }

    } else if (initMethod == "fanny") {
      z[[iterations]] <- mclust::unmap(cluster::fanny(log(dataset + 1 / 3),
        k = gmodel)$cluster)
      # if z generated has less columns than gmodel, then use random initialization
      if(ncol(z[[iterations]]) < gmodel) {
        z[[iterations]] <- randomInitfunction(gmodel = gmodel, nObservations = nObservations)
      }
    }

    initRuns[[iterations]] <- mplnCluster(dataset = dataset,
                                           z = z[[iterations]],
                                           G = gmodel,
                                           nChains = nChain,
                                           nIterations = numbIterations,
                                           initialization = "init",
                                           normalizefac = normalizefactors,
                                           mod = mod)
    logLinit[iterations] <-
      unlist(utils::tail((initRuns[[iterations]]$logLikelihood), n = 1))
  }

  initialization <- initRuns[[which(logLinit == max(logLinit, na.rm = TRUE))[1]]]
  return(initialization)
  # Developed by Anjali Silva
}



# Clustering function
mplnCluster <- function(dataset,
                        z,
                        G,
                        nChains,
                        nIterations,
                        initialization,
                        normalizefac,
                        mod) {

  dimensionality <- ncol(dataset)
  nObservations <- nrow(dataset)

  # for convergence calculation
  normMuOuter <- normSigmaOuter <- vector()
  medianMuOuter <- medianSigmaOuter <- list()
  # for saving mu and sigma values
  muAllOuter <- sigmaAllOuter <- list()
  obs <- PI <- logL <- vector()
  nIterationsOuter <- 2 # the starting value of interation for outer loop
  convOuter <- 0


  if (all(is.na(initialization)) == TRUE || all(initialization == "init")) {
    # mean for both t and normal distribution
    muAllOuter[[1]] <- muG <- matrix(log(mean(dataset)),
      ncol = dimensionality,
      nrow = G)
    # sig for sigma of t distribtuion
    sigmaAllOuter[[1]] <- sigmaG <- do.call("rbind",
      rep(list(cov(log(dataset + 1)) *
          dimensionality), G))
  } else {
    muAllOuter[[1]] <- muG <- initialization$finalmu
    sigmaAllOuter[[1]] <- sigmaG <- initialization$finalsigma
    z <- initialization$probaPost
  }

  while(! convOuter) {
    for(g in seq_along(1:G)) {
      obs[g] <- sum(z[, g]) # number of observations in each group
      PI[g] <- obs[g] / nObservations  # obtain probability of each group
    }

    thetaStan <- ExpTheta2 <- list()
    rstanResults <- stanRun(model = mod,
                             gmin = 1,
                             gmax = G,
                             dataset = dataset,
                             muAllOuter = muAllOuter,
                             nIterationsOuter = nIterationsOuter,
                             sigmaAllOuter = sigmaAllOuter,
                             numbIterations = nIterations,
                             nChain = nChains,
                             normalizefacs = normalizefac)

    fit <- rstanResults$fitrstan
    nIterations <- rstanResults$numbIterations

    for (g in seq_along(1:G)) {
      tt <- as.matrix(fit[[g]])
      thetaStan[[g]] <- matrix(NA, nrow = nObservations, ncol = dimensionality)
      ExpTheta2[[g]] <- list()

      for (i in seq_along(1:nObservations)) {
        zz <- c(1:(dimensionality - 1)) * nObservations + i
        thetaMatrix <- tt[, c(i, zz)]
        thetaStan[[g]][i, ] <- colMeans(thetaMatrix)
        ExpTheta2[[g]][[i]] <- z[i, g] * t(tt[, c(i, zz)]) %*% tt[, c(i, zz)]/
          ((0.5 * nIterations) * nChains)
      }

      muG[g, ] <- colSums(z[, g] * thetaStan[[g]]) / sum(z[, g])
      sigmaG[((g - 1) * dimensionality + 1):(g * dimensionality), ] <-
        Reduce("+", ExpTheta2[[g]]) / sum(z[, g]) - muG[g, ] %*% t(muG[g, ])
    }

    muAllOuter[[nIterationsOuter]] <- muG
    sigmaAllOuter[[nIterationsOuter]] <- sigmaG

    logL[nIterationsOuter] <- calcLikelihood(z = z,
                                    PI = PI,
                                    dataset = dataset,
                                    muG = muAllOuter[[nIterationsOuter]],
                                    G = G,
                                    sigmaG = sigmaAllOuter[[nIterationsOuter]],
                                    thetaStan = thetaStan,
                                    normFactors = normalizefac)

    # convergence of outer loop
    normMuOuter[nIterationsOuter] <- norm((muAllOuter[[nIterationsOuter]] -
                                           muAllOuter[[nIterationsOuter - 1]]),
                                           type = "F")
    normSigmaOuter[nIterationsOuter] <- norm(sigmaAllOuter[[nIterationsOuter]] -
                                             sigmaAllOuter[[nIterationsOuter - 1]],
                                             type = "F")
    medianMuOuter[[nIterationsOuter]] <- median(normMuOuter, na.rm = TRUE)
    medianSigmaOuter[[nIterationsOuter]] <- median(normSigmaOuter, na.rm = TRUE)
    # par(mfrow = c(1, 2))
    # plot(normMuOuter, main = paste0("Norm outer mean, G=", G),
    # type = "l", ylab = "median(normMuOuter)", xlab = "iterations")
    # plot(normSigmaOuter, main = paste0("Norm outer sigma, G = ", G),
    # type = "l", ylab = "median(normSigmaOuter)", xlab = "iterations")

    thresholdOuter <- 2
    if(nIterationsOuter > (thresholdOuter + 1)) {

      # cat("\nMedian difference of mean and sigma in outer loop respectively ",
      # c(abs(medianMuOuter[[nIterationsOuter - thresholdOuter]]-
      # medianMuOuter[[nIterationsOuter]])))
      if( ( (abs(medianMuOuter[[nIterationsOuter - thresholdOuter]] -
          medianMuOuter[[nIterationsOuter]]) < 5) &&
          (abs(medianSigmaOuter[[nIterationsOuter - thresholdOuter]] -
              medianSigmaOuter[[nIterationsOuter]]) < 5) ) || nIterationsOuter > 100) {
        # cat("\nConvergence of mu and sigma at outer loop
        # iteration ", nIterationsOuter)
        # take out absolute value
        programclust <- vector()
        programclust <- mclust::map(z)

        # checking for spurious clusters and getting rid of them
        # keep <- as.numeric(names(which(table(programclust) > 5)))
        # if ((length(keep) != length(unique(programclust))) &&
        # (length(keep) != 0)) {
        #  z <- as.matrix(z[ , keep])
        #  z <- z / rowSums(z)
        #  programclust <- map(z)
        #}

        # checking for empty clusters
        J <- 1:ncol(z)
        K <- as.logical(match(J, sort(unique(programclust)), nomatch = 0))
        if(length(J[! K]) > 0) { # J[!K] tells which are empty clusters
          z <- z[, - J[! K]]
          programclust <- mclust::map(z)
        }
        convOuter <- 1
      }
    }

    # if running for initialization, need to stop after 10 iterations
    if(nIterationsOuter == 10 && all(is.na(initialization) != TRUE)) {
      if(all(initialization == "init")) {
        programclust <- vector()
        programclust <- mclust::map(z)
        convOuter <- 1
      }
    }

    if(convOuter != 1) { # only update until convergence, not after
      z <- calcZvalue(thetaStan = thetaStan,
                      dataset = dataset,
                      G = G,
                      muG = muG,
                      sigmaG = sigmaG,
                      PI = PI,
                      normalizefactors = normalizefac)
      nIterationsOuter <- nIterationsOuter + 1 # updating outer loop iteration
      nIterations <- nIterations + 10
    }
  } # end of outer loop

  results <- list(finalmu = muAllOuter[[nIterationsOuter]] +
                    matrix(rep(normalizefac,
                      nrow(muAllOuter[[nIterationsOuter]])),
                      byrow = TRUE, ncol =
                        ncol(muAllOuter[[nIterationsOuter]])),
                  finalsigma = sigmaAllOuter[[nIterationsOuter]],
                  allmu = lapply(muAllOuter, function(x)
                    (x + matrix(rep(normalizefac,
                      nrow(muAllOuter[[nIterationsOuter]])),
                      byrow = TRUE,
                      ncol = ncol(muAllOuter[[nIterationsOuter]])))),
                  allsigma = sigmaAllOuter,
                  clusterlabels = programclust,
                  iterations = nIterationsOuter,
                  proportion = PI,
                  logLikelihood = logL,
                  probaPost = z,
                  stanresults = fit)

  class(results) <- "MPLNcluster"
  return(results)
  # Developed by Anjali Silva
}



# Calculate and remove rows with zeros
#' @author {Anjali Silva, \email{anjali.silva@uhnresearch.ca}}
removeZeroCounts <- function(dataset,
                             membership = "none") {

  zeroSUMrows <- which(rowSums(dataset) == 0)

  if (length(zeroSUMrows) > 0 && class(membership) == "numeric") {
    dataset <- dataset[- zeroSUMrows, ]
    membership <- membership[- zeroSUMrows]
  } else if(length(zeroSUMrows) > 0 && all(membership == "none")) {
    dataset <- dataset[- zeroSUMrows, ]
    membership <- "none"
  }

  RESULTS <- list(dataset = dataset,
                  membership = membership)
  class(RESULTS) <- "MPLNZerosRemoved"
  return(RESULTS)
}



# Stan sampling
stanRun <- function(model,
                    gmin,
                    gmax,
                    dataset,
                    muAllOuter,
                    nIterationsOuter,
                    sigmaAllOuter,
                    numbIterations,
                    nChain = nChain,
                    normalizefacs) {

  fitrstan <- list()
  dimensionality <- ncol(dataset)

  for (g in seq_along(gmin:gmax)) {
    data1 <- list(d = ncol(dataset),
                  N = nrow(dataset),
                  y = dataset,
                  mu = muAllOuter[[nIterationsOuter - 1]][g, ],
                  Sigma = sigmaAllOuter[[nIterationsOuter - 1]][((g - 1) *
                      dimensionality + 1):(g*dimensionality), ],
                  normfactors = as.vector(normalizefacs))
    stanproceed <- 0
    try <- 1

    while (! stanproceed) {

      # cat("\nRstan generating sample at outer iteration",
      # nIterationsOuter, "for g: ",g , "try: ", try)
      # cat("\nNumber of iterations is", numbIterations, "\n")
      fitrstan[[g]] <- rstan::sampling(object = model,
                                       data = data1,
                                       iter = numbIterations,
                                       chains = nChain,
                                       verbose = FALSE,
                                       refresh = -1)

      # Checking convergence
      if (all(rstan::summary(fitrstan[[g]])$summary[, "Rhat"] < 1.1) == TRUE &&
          all(rstan::summary(fitrstan[[g]])$summary[, "n_eff"] > 100) == TRUE) {
        stanproceed <- 1
      } else if(all(rstan::summary(fitrstan[[g]])$summary[, "Rhat"] < 1.1) != TRUE ||
          all(rstan::summary(fitrstan[[g]])$summary[, "n_eff"] > 100) != TRUE) {
        if(try == 10) { # Stop after 10 attempts
          stanproceed <- 1
        }
        numbIterations <- numbIterations + 100
        try <- try + 1
      }
    }
  }

  results <- list(fitrstan = fitrstan,
                  numbIterations = numbIterations)
  class(results) <- "MPLNRStan"
  return(results)
  # Developed by Anjali Silva
}

# [END]
