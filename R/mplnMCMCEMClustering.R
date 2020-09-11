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
#'      results as obtained via capushe::capushe().
#'   \item DjumpModelSelected - If more than 10 models are considered, slope heuristic
#'      results as obtained via capushe::capushe().
#'   \item DDSEModelSelected - If more than 10 models are considered, slope heuristic
#'      results as obtained via capushe::capushe().
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
#' # sampleData <- MPLNClust::mplnDataGenerator(nObservations = 40,
#' #                                             dimensionality = 6,
#' #                                             mixingProportions = c(0.79, 0.21),
#' #                                             mu = rbind(trueMu1, trueMu2),
#' #                                             sigma = rbind(trueSigma1, trueSigma2),
#' #                                             produceImage = "No")
#'
#' # Clustering
#' # mplnResults <- MPLNClust::mplnMCMCParallel(dataset = sampleData$dataset,
#' #                                             membership = sampleData$trueMembership,
#' #                                             gmin = 1,
#' #                                             gmax = 1,
#' #                                             nChains = 3,
#' #                                             nIterations = 400,
#' #                                             initMethod = "kmeans",
#' #                                             nInitIterations = 0,
#' #                                             normalize = "Yes",
#' #                                             numNodes = 2)
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
#' Arlot, S., Brault, V., Baudry, J., Maugis, C., and Michel, B. (2016).
#' capushe: CAlibrating Penalities Using Slope HEuristics. R package version 1.1.1.
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
#' Robinson, M.D., and Oshlack, A. (2010). A scaling normalization method for differential
#' expression analysis of RNA-seq data. \emph{Genome Biology} 11, R25.
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

  if (is.matrix(dataset) != TRUE) {
    stop("Dataset needs to be a matrix.")
  }

  if (any(colSums(dataset) <= 0)) {
    stop("Column sums cannot be less than or equal to 0. Double check dataset.")
  }

  dimensionality <- ncol(dataset)
  nObservations <- nrow(dataset)

  if(is.numeric(gmin) != TRUE || is.numeric(gmax) != TRUE) {
    stop("Class of gmin and gmin should be numeric.")
  }

  if (gmax < gmin) {
    stop("gmax cannot be less than gmin.")
  }

  if (gmax > nObservations) {
    stop("gmax cannot be larger than nrow(dataset).")
  }

  if (is.numeric(nChains) != TRUE) {
    stop("nChains should be a positive integer of class numeric specifying
      the number of Markov chains.")
  }

  if(nChains < 3) {
    cat("nChains is less than 3. Note: recommended minimum number of
      chains is 3.")
  }

  if(is.numeric(nIterations) != TRUE) {
    stop("nIterations should be a positive integer of class numeric,
      specifying the number of iterations for each Markov chain
      (including warmup).")
  }

  if(is.numeric(nIterations) == TRUE && nIterations < 40) {
    stop("nIterations argument should be greater than 40.")
  }

  if(all(membership != "none") && is.numeric(membership) != TRUE) {
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

  if (is.character(initMethod) == TRUE) {
    initMethodsUsed <- c("kmeans", "random", "medoids", "clara", "fanny")
    if(all((initMethod == initMethodsUsed) == FALSE)) {
      stop("initMethod should of class character, specifying
        either: kmeans, random, medoids, clara, or fanny.")
    }
    } else if (is.character(initMethod) != TRUE) {
      stop("initMethod should of class character, specifying
        either: kmeans, random, medoids, clara, or fanny.")
  }

  if (is.numeric(nInitIterations) != TRUE) {
    stop("nInitIterations should be positive integer or zero, specifying
      the number of initialization runs to be considered.")
  }

  if (is.character(normalize) != TRUE) {
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
  } else if (is.numeric(numNodes) == TRUE) {
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
#'      results as obtained via capushe::capushe().
#'   \item DjumpModelSelected - If more than 10 models are considered, slope heuristic
#'      results as obtained via capushe::capushe().
#'   \item DDSEModelSelected - If more than 10 models are considered, slope heuristic
#'      results as obtained via capushe::capushe().
#'   \item totalTime - Total time used for clustering and model selection.
#' }
#'
#' @examples
#' # Not run
#' # Generating simulated data
#'
#' # trueMu1 <- c(6.5, 6, 6, 6, 6, 6)
#' # trueMu2 <- c(2, 2.5, 2, 2, 2, 2)
#'
#' # trueSigma1 <- diag(6) * 2
#' # trueSigma2 <- diag(6)
#'
#' # sampleData <- MPLNClust::mplnDataGenerator(nObservations = 40,
#' #                                 dimensionality = 6,
#' #                                 mixingProportions = c(0.79, 0.21),
#' #                                 mu = rbind(trueMu1, trueMu2),
#' #                                 sigma = rbind(trueSigma1, trueSigma2),
#' #                                 produceImage = "No")
#'
#' # Clustering
#' # mplnResults <- MPLNClust::mplnMCMCNonParallel(dataset = sampleData$dataset,
#' #                                               membership = sampleData$trueMembership,
#' #                                               gmin = 1,
#' #                                               gmax = 1,
#' #                                               nChains = 3,
#' #                                               nIterations = 700,
#' #                                               initMethod = "kmeans",
#' #                                               nInitIterations = 0,
#' #                                               normalize = "Yes")
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
#' Arlot, S., Brault, V., Baudry, J., Maugis, C., and Michel, B. (2016).
#' capushe: CAlibrating Penalities Using Slope HEuristics. R package version 1.1.1.
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
#' Robinson, M.D., and Oshlack, A. (2010). A scaling normalization method for differential
#' expression analysis of RNA-seq data. \emph{Genome Biology} 11, R25.
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
#' @import cluster
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

  if (is.matrix(dataset) != TRUE) {
    stop("Dataset needs to be a matrix.")
  }

  if (any(colSums(dataset) <= 0)) {
    stop("Column sums cannot be less than or equal to 0. Double check dataset.")
  }

  dimensionality <- ncol(dataset)
  nObservations <- nrow(dataset)

  if(is.numeric(gmin) != TRUE || is.numeric(gmax) != TRUE) {
    stop("Class of gmin and gmin should be numeric.")
  }

  if (gmax < gmin) {
    stop("gmax cannot be less than gmin.")
  }

  if (gmax > nObservations) {
    stop("gmax cannot be larger than nrow(dataset).")
  }

  if (is.numeric(nChains) != TRUE) {
    stop("nChains should be a positive integer of class numeric specifying
      the number of Markov chains.")
  }

  if(nChains < 3) {
    cat("nChains is less than 3. Note: recommended minimum number of
      chains is 3.")
  }

  if(is.numeric(nIterations) != TRUE) {
    stop("nIterations should be a positive integer of class numeric,
      specifying the number of iterations for each Markov chain
      (including warmup).")
  }

  if(is.numeric(nIterations) == TRUE && nIterations < 40) {
    stop("nIterations argument should be greater than 40.")
  }

  if(all(membership != "none") && is.numeric(membership) != TRUE) {
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

  if (is.character(initMethod) == TRUE) {
    initMethodsUsed <- c("kmeans", "random", "medoids", "clara", "fanny")
    if(all((initMethod == initMethodsUsed) == FALSE)) {
      stop("initMethod should of class character, specifying
        either: kmeans, random, medoids, clara, or fanny.")
    }
    } else if (is.character(initMethod) != TRUE) {
      stop("initMethod should of class character, specifying
        either: kmeans, random, medoids, clara, or fanny.")
  }


  if (is.numeric(nInitIterations) != TRUE) {
    stop("nInitIterations should be positive integer or zero, specifying
      the number of initialization runs to be considered.")
  }

  if (is.character(normalize) != TRUE) {
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

  if (length(zeroSUMrows) > 0 && is.numeric(membership) == TRUE) {
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








#' Clustering Using MPLN Via Variational-EM
#'
#' Performs clustering using mixtures of multivariate Poisson-log
#' normal (MPLN) distribution with variational expectation-maximization
#' (EM) for parameter estimation. Model selection is performed using
#' AIC, AIC3, BIC and ICL. If more than 10 models are considered, Djump
#' and DDSE is also applied for model selection. No internal
#' parallelization, thus code is run in serial.
#'
#' @param dataset A dataset of class matrix and type integer such that
#'    rows correspond to observations and columns correspond to variables.
#'    The dataset have dimensions n x d, where n is the total number of
#'    observations and d is the dimensionality. If rowSums are zero, these
#'    rows will be removed prior to cluster analysis.
#' @param membership A numeric vector of length nrow(dataset) containing the
#'    cluster membership of each observation. If not available,
#'    leave as "none".
#' @param gmin A positive integer specifying the minimum number of components
#'    to be considered in the clustering run.
#' @param gmax A positive integer, >= gmin, specifying the maximum number of
#'    components to be considered in the clustering run.
#' @param initMethod An algorithm for initialization. Current options are
#'    "kmeans", "random", "medoids", "clara", or "fanny". Default is "kmeans".
#' @param nInitIterations A positive integer or zero, specifying the number
#'    of initialization runs to be performed. This many runs, each with 10
#'    iterations, will be performed via MPLNClust and values from the run with
#'    highest log-likelihood will be used as initialization values. Default is 2.
#' @param normalize A string with options "Yes" or "No" specifying
#'     if normalization should be performed. Currently, normalization factors
#'     are calculated using TMM method of edgeR package. Default is "Yes".
#'
#' @return Returns an S3 object of class mplnVariational with results.
#' \itemize{
#'   \item dataset - The input dataset on which clustering is performed.
#'   \item dimensionality - Dimensionality of the input dataset.
#'   \item normalizationFactors - A vector of normalization factors used
#'      for input dataset.
#'   \item gmin - Minimum number of components/clusters considered in the clustering
#'      run.
#'   \item gmax - Maximum number of components/clusters considered in the clustering
#'      run.
#'   \item initalizationMethod - Method used for initialization.
#'   \item allResults - A list with all results.
#'   \item logLikelihood - A vector with value of final log-likelihoods for
#'      each component/cluster size.
#'   \item numbParameters - A vector with number of parameters for each
#'      component/cluster size.
#'   \item trueLabels - The vector of true labels, if provided by user.
#'   \item ICLresults - A list with all ICL model selection results.
#'   \item BICresults - A list with all BIC model selection results.
#'   \item AICresults - A list with all AIC model selection results.
#'   \item AIC3results - A list with all AIC3 model selection results.
#'   \item slopeHeuristics - If more than 10 models are considered, slope heuristic
#'      results as obtained via capushe::capushe().
#'   \item DjumpModelSelected - If more than 10 models are considered, slope heuristic
#'      results as obtained via capushe::capushe().
#'   \item DDSEModelSelected - If more than 10 models are considered, slope heuristic
#'      results as obtained via capushe::capushe().
#'   \item totalTime - Total time used for clustering and model selection.
#' }
#'
#' @examples
#' # Generating simulated data
#'
#'  trueMu1 <- c(6.5, 6, 6, 6, 6, 6)
#'  trueMu2 <- c(2, 2.5, 2, 2, 2, 2)
#'
#'  trueSigma1 <- diag(6) * 2
#'  trueSigma2 <- diag(6)
#'
#'  sampleData <- MPLNClust::mplnDataGenerator(nObservations = 1000,
#'                                             dimensionality = 6,
#'                                             mixingProportions = c(0.79, 0.21),
#'                                             mu = rbind(trueMu1, trueMu2),
#'                                             sigma = rbind(trueSigma1, trueSigma2),
#'                                             produceImage = "No")
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
#' Arlot, S., Brault, V., Baudry, J., Maugis, C., and Michel, B. (2016).
#' capushe: CAlibrating Penalities Using Slope HEuristics. R package version 1.1.1.
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
#' Robinson, M.D., and Oshlack, A. (2010). A scaling normalization method for differential
#' expression analysis of RNA-seq data. \emph{Genome Biology} 11, R25.
#'
#' Schwarz, G. (1978). Estimating the dimension of a model. \emph{The Annals of Statistics}
#' 6.
#'
#' Silva, A. et al. (2019). A multivariate Poisson-log normal mixture model
#' for clustering transcriptome sequencing data. \emph{BMC Bioinformatics} 20.
#' \href{https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-2916-0}{Link}
#'
#' Subedi, S., and R. Browne (2020). A parsimonious family of multivariate Poisson-lognormal
#' distributions for clustering multivariate count data. arXiv preprint arXiv:2004.06857.
#' \href{https://arxiv.org/pdf/2004.06857.pdf}{Link}
#'
#' @export
#' @importFrom capushe capushe
#' @importFrom edgeR calcNormFactors
#' @importFrom mclust unmap
#' @importFrom mclust map
#' @import stats
#' @import cluster
#'
mplnVariational <- function(dataset,
                            membership = "none",
                            gmin,
                            gmax,
                            initMethod = "kmeans",
                            nInitIterations = 2,
                            normalize = "Yes") {

  initialTime <- proc.time()

  # Performing checks
  if (typeof(dataset) != "double" & typeof(dataset) != "integer") {
    stop("Dataset type needs to be integer.")
  }

  if (any((dataset %% 1 == 0) == FALSE)) {
    stop("Dataset should be a matrix of counts.")
  }

  if (is.matrix(dataset) != TRUE) {
    stop("Dataset needs to be a matrix.")
  }

  if (any(colSums(dataset) <= 0)) {
    stop("Column sums cannot be less than or equal to 0. Double check dataset.")
  }

  dimensionality <- ncol(dataset)
  nObservations <- nrow(dataset)

  if(is.numeric(gmin) != TRUE || is.numeric(gmax) != TRUE) {
    stop("Class of gmin and gmin should be numeric.")
  }

  if (gmax < gmin) {
    stop("gmax cannot be less than gmin.")
  }

  if (gmax > nObservations) {
    stop("gmax cannot be larger than nrow(dataset).")
  }

  if(all(membership != "none") && is.numeric(membership) != TRUE) {
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

  if (is.character(initMethod) == TRUE) {
    initMethodsUsed <- c("kmeans", "random", "medoids", "clara", "fanny")
    if(all((initMethod == initMethodsUsed) == FALSE)) {
      stop("initMethod should of class character, specifying
        either: kmeans, random, medoids, clara, or fanny.")
    }
  } else if (is.character(initMethod) != TRUE) {
    stop("initMethod should of class character, specifying
        either: kmeans, random, medoids, clara, or fanny.")
  }

  if (is.numeric(nInitIterations) != TRUE) {
    stop("nInitIterations should be positive integer or zero, specifying
      the number of initialization runs to be considered.")
  }

  if (is.character(normalize) != TRUE) {
    stop("normalize should be a string of class character specifying
      if normalization should be performed.")
  }

  if ((normalize != "Yes") && (normalize != "No")) {
    stop("normalize should be a string indicating Yes or No, specifying
      if normalization should be performed.")
  }


  if(normalize == "Yes") {
    normFactors <- as.vector(edgeR::calcNormFactors(as.matrix(dataset),
                                                    method = "TMM"))
  } else if(normalize == "No") {
    normFactors <- rep(1, dimensionality)
  } else{
    stop("normalize should be 'Yes' or 'No'.")
  }

  clusterResults <- list() # to save cluster output
  for (gmodel in seq_along(1:(gmax - gmin + 1))) {

    if(length(1:(gmax - gmin + 1)) == gmax) {
      clustersize <- gmodel
    } else if(length(1:(gmax - gmin + 1)) < gmax) {
      clustersize <- seq(gmin, gmax, 1)[gmodel]
    }
    cat("\n Running for g =", clustersize)
    clusterResults[[gmodel]] <- varMPLNClustering(dataset = dataset,
                                                  initMethod = initMethod,
                                                  nInitIterations = nInitIterations,
                                                  G = clustersize,
                                                  maxIterations = 1000,
                                                  normFactors = normFactors)
  }

  names(clusterResults) <- paste0(rep("G=", length(seq(gmin, gmax, 1))), seq(gmin, gmax, 1))

  BIC <- ICL <- AIC <- AIC3 <- Djump <- DDSE <- nParameters <- logLikelihood <- vector()

  for(g in seq_along(1:(gmax - gmin + 1))) {
    # save the final log-likelihood

    if(length(1:(gmax - gmin + 1)) == gmax) {
      clustersize <- g
    } else if(length(1:(gmax - gmin + 1)) < gmax) {
      clustersize <- seq(gmin, gmax, 1)[g]
    }

    logLikelihood[g] <- unlist(utils::tail(clusterResults[[g]]$logLikelihood, n = 1))

    nParameters[g] <- calcParameters(numberG = clustersize, dimensionality = dimensionality)

    if (g == max(1:(gmax - gmin + 1))) { # starting model selection
      bic <- BICFunction(logLikelihood = logLikelihood,
                         nParameters = nParameters,
                         nObservations = nObservations,
                         clusterRunOutput = clusterResults,
                         gmin = gmin,
                         gmax = gmax,
                         parallel = FALSE)

      icl <- ICLFunction(logLikelihood = logLikelihood,
                         nParameters = nParameters,
                         nObservations = nObservations,
                         gmin = gmin,
                         gmax = gmax,
                         clusterRunOutput = clusterResults,
                         parallel = FALSE)

      aic <- AICFunction(logLikelihood = logLikelihood,
                         nParameters = nParameters,
                         clusterRunOutput = clusterResults,
                         gmin = gmin,
                         gmax = gmax,
                         parallel = FALSE)

      aic3 <- AIC3Function(logLikelihood = logLikelihood,
                           nParameters = nParameters,
                           clusterRunOutput = clusterResults,
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
    ResCapushe <- capushe::capushe(mat, nObservations)
    DDSEmodel<- ResCapushe@DDSE@model
    Djumpmodel<- ResCapushe@Djump@model

    finalTime <- proc.time() - initialTime

    RESULTS <- list(dataset = dataset,
                    dimensionality = dimensionality,
                    normalizationFactors = normFactors,
                    gmin = gmin,
                    gmax = gmax,
                    initalizationMethod = initMethod,
                    allResults = clusterResults,
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
                    totalTime = finalTime)

  } else {
    finalTime <- proc.time() - initialTime

    RESULTS <- list(dataset = dataset,
                    dimensionality = dimensionality,
                    normalizationFactors = normFactors,
                    gmin = gmin,
                    gmax = gmax,
                    initalizationMethod = initMethod,
                    allResults = clusterResults,
                    logLikelihood = logLikelihood,
                    numbParameters = nParameters,
                    trueLabels = membership,
                    ICLresults = icl,
                    BICresults = bic,
                    AICresults = aic,
                    AIC3results = aic3,
                    slopeHeuristics = "Not used",
                    totalTime = finalTime)
  }

  class(RESULTS) <- "mplnVariational"
  return(RESULTS)
}






varMPLNClustering <- function(dataset,
                              initMethod,
                              nInitIterations,
                              G,
                              maxIterations,
                              normFactors) {

  nObservations <- nrow(dataset)
  dimensionality <- ncol(dataset)

  # If no intialization is requested by user
  if (nInitIterations == 0) {

    # basic initialization performed for parameters
    mu <- sigma <- isigma <- list()
    m <- S <- list() # variational parameters

    ###Other intermediate items initialized
    Sk <- array(0, c(dimensionality, dimensionality, G))
    mPreviousValue <- GX <- dGX <- zSValue <- list()


    zValue <- mclust::unmap(stats::kmeans(log(dataset + 1 / 6),
                                          centers = G, nstart = 100)$cluster )

    piG <- colSums(zValue) / nObservations

    for (g in 1:G) {
      obs <- which(zValue[ , g] == 1)
      # starting value for mu
      mu[[g]] <- colMeans(log(dataset[obs, ] + 1 / 6))
      # starting value for sample covariance matrix
      sigma[[g]] <- var(log(dataset[obs, ] + 1 / 6))
      # starting value for inverse of sample covariance matrix
      # If the inverse is not present for covariance matrix, handle that
      isigma[[g]] <- tryCatch(solve(sigma[[g]]), error = function(err) NA)
      if(all(is.na(isigma[[g]]))) {
        isigma[[g]] <- solve(sigma[[g]] + 0.001) # if error with inverse
      }
    }

    for (g in 1:G) {
      mPreviousValue[[g]] <- m[[g]] <- log(dataset + 1 / 6)
      # mPreviousValue = starting value for m
      S[[g]] <- list() # starting value for S
      for (i in 1:nObservations) {
        S[[g]][[i]] <- diag(dimensionality) * 0.000000001
      }
    }
  } else if (nInitIterations != 0) {
    # if initialization is requested by user

    initializationResults <- varMPLNInitialization(dataset = dataset,
                                                   numbG = G,
                                                   initMethod = initMethod,
                                                   nInitIterations = nInitIterations,
                                                   normFactors = normFactors)

    mu <- initializationResults$mu
    sigma <- initializationResults$sigma
    isigma <- initializationResults$isigma
    mPreviousValue <- m <- initializationResults$m # variational parameters
    S <- initializationResults$S # variational parameters
    zValue <- initializationResults$zValue
    piG <- colSums(zValue) / nObservations

    ###Other intermediate items initialized
    Sk <- array(0, c(dimensionality, dimensionality, G))
    GX <- dGX <- zSValue <- list()

  }






  # Start clustering
  itOuter <- 1
  aloglik <- logLikelihood <- NULL
  checks <- aloglik[c(1, 2, 3)] <- 0


  while (checks == 0) {


    for (g in 1:G) {
      GX[[g]] <- list()
      dGX[[g]] <- list()
      zSValue[[g]] <- list()
      for (i in 1:nObservations) {
        dGX[[g]][[i]] <- diag(exp((log(normFactors) + mPreviousValue[[g]][i, ]) +
                                    0.5 * diag(S[[g]][[i]])), dimensionality) + isigma[[g]]
        S[[g]][[i]] <- solve(dGX[[g]][[i]]) # update S
        zSValue[[g]][[i]] <- zValue[i, g] * S[[g]][[i]] # will be used for updating sample covariance matrix
        GX[[g]][[i]] <- dataset[i, ] - exp(mPreviousValue[[g]][i, ] +
                                             log(normFactors) + 0.5 * diag(S[[g]][[i]])) -
          (isigma[[g]]) %*% (mPreviousValue[[g]][i, ] - mu[[g]])
        m[[g]][i , ] <- mPreviousValue[[g]][i , ] + S[[g]][[i]] %*% GX[[g]][[i]] # update m
      }

      mPreviousValue[[g]] <- m[[g]]

      # Updating mu
      mu[[g]] <- colSums(zValue[ , g] * m[[g]]) / sum(zValue[ , g])

      # Updating sample covariance
      muMatrix <- matrix(rep(mu[[g]], nObservations), nrow = nObservations, byrow = TRUE)
      res <- m[[g]] - muMatrix
      # Calculate weighted covariance
      temp <- stats::cov.wt(res, wt = zValue[ , g], center = FALSE, method = "ML")
      Sk[ , , g] <- temp$cov
      sigma[[g]] <- temp$cov + Reduce("+", zSValue[[g]]) / sum(zValue[ , g])
      isigma[[g]] <- solve(sigma[[g]])
    }

    piG <- colSums(zValue) / nObservations

    # Matrix containing normaization factors
    normFactorsAsMatrix <- matrix(rep(normFactors, nObservations), nrow = nObservations, byrow = TRUE)

    # A useful function
    zvalueFuncTerm <- function(x, y = isigma[[g]]) {
      sum(diag(x %*% y))
    }

    # # Zvalue calculation
    forz <- matrix(NA, ncol = G, nrow = nObservations)

    for (g in 1:G) {
      # exp(m_igj + 0.5 S_ig,jj)
      two <- rowSums(exp(m[[g]] + log(normFactorsAsMatrix) +
                           0.5 * matrix(unlist(lapply(S[[g]], diag)), ncol = dimensionality, byrow = TRUE)))
      five <- 0.5 * unlist(lapply(S[[g]], zvalueFuncTerm)) # trace(isigma x S_ig)
      six <- 0.5 * log(unlist(lapply(S[[g]], det))) # 0.5 * log|S_ig|

      # For zig calculation (the numerator part)
      forz[ , g] <- piG[g] *
        exp(rowSums(m[[g]] * dataset) - two - rowSums(lfactorial(dataset)) +
              rowSums(log(normFactorsAsMatrix) * dataset) -
              0.5 * mahalanobis(m[[g]], center = mu[[g]], cov = sigma[[g]]) -
              five + six - 0.5 * log(det(sigma[[g]])) - dimensionality / 2)
    }


    # Calculate zValue value
    # check which forz == 0 and rowSums(forz)==0 and which of these
    # have both equalling to 0 (because 0/0 =NaN)
    if (G == 1) {
      errorpossible <- Reduce(intersect,
                              list(which(forz == 0),
                                   which(rowSums(forz) == 0)))
      forz[errorpossible] <- 1e-100
      zvalue <- forz / rowSums(forz)
    } else {

      # check for error, if rowsums are zero
      rowSumsZero <- which(rowSums(forz) == 0)
      if(length(rowSumsZero) > 1) {
        forz[rowSumsZero, ] <- mclust::unmap(stats::kmeans(log(dataset + 1 / 6),
                                                           centers = G,
                                                           nstart = 100)$cluster)[rowSumsZero, ]
        zvalue <- forz / rowSums(forz)
      } else {
        zvalue <- forz / rowSums(forz)
      }
    }


    # Calculate log-likelihood
    logLikelihood[itOuter] <- sum(log(rowSums(forz)))

    # Stopping criterion
    if (itOuter > 2) {

      if ((logLikelihood[itOuter - 1] - logLikelihood[itOuter - 2]) == 0) {
        checks <-1
      } else {
        # Aitken stopping criterion
        termAitkens <- (logLikelihood[itOuter]- logLikelihood[itOuter - 1]) /
          (logLikelihood[itOuter - 1] - logLikelihood[itOuter - 2])
        term2Aitkens <- (1 / (1 - termAitkens) * (logLikelihood[itOuter] -
                                                    logLikelihood[itOuter - 1]))
        aloglik[itOuter] <- logLikelihood[itOuter - 1] + term2Aitkens
        if (abs(aloglik[itOuter] - logLikelihood[itOuter - 1]) < 0.001) {
          # If this critera, as per Böhning et al., 1994 is achieved
          # convergence is achieved
          checks <- 1
        } else {
          checks <- checks}
      }
    }
    # print(itOuter)
    itOuter <- itOuter + 1
    if (itOuter == maxIterations) {
      checks <- 1
    }
  }

  # Naming parameters
  names(mu) <- names(sigma) <- paste0(rep("G=", G), 1:G)

  # Saving results for output
  Results <- list(piG = piG,
                  mu = mu,
                  sigma = sigma,
                  probaPost = zValue,
                  clusterlabels = mclust::map(zValue),
                  logLikelihood = logLikelihood)

  class(Results) <- "varMPLNClustering"
  return(Results)
}



varMPLNInitialization <- function(dataset,
                                  numbG,
                                  initMethod,
                                  nInitIterations,
                                  normFactors) {

  dimensionality <- ncol(dataset)
  nObservations <- nrow(dataset)

  zValue <- initRuns <- list()
  logLinit <- vector()

  # Internal function for random initialization
  randomInitfunction <- function(numbG, nObservations) {
    if(numbG == 1) { # generating zValue if g=1
      zValue <- as.matrix(rep.int(1, times = nObservations),
                          ncol = numbG,
                          nrow = nObservations)
    } else { # generating zValue if g>1
      zValueConv <- 0
      while(! zValueConv) {
        # ensure that dimension of zValue is same as G (i.e.,
        # if one column contains all 0s, then generate zValue again)
        zValue <- t(stats::rmultinom(nObservations, size = 1,
                                     prob = rep(1 / numbG, numbG)))
        if(length(which(colSums(zValue) > 0)) == numbG) {
          zValueConv <- 1
        }
      }
    }
    return(zValue)
  }

  for(iterations in seq_along(1:nInitIterations)) {
    # setting seed, to ensure if multiple iterations are selected by
    # user, then each run will give a different result.
    set.seed(iterations)
    if (initMethod == "kmeans" | is.na(initMethod)) {
      zValue[[iterations]] <- mclust::unmap(stats::kmeans(log(dataset + 1 / 6),
                                                          centers = numbG, nstart = 100)$cluster )
      # if z generated has less columns than numbG, then use random initialization
      if(ncol(zValue[[iterations]]) < numbG) {
        zValue[[iterations]] <- randomInitfunction(numbG = numbG, nObservations = nObservations)
      }

    } else if (initMethod == "random") {
      zValue[[iterations]] <- randomInitfunction(numbG = numbG, nObservations = nObservations)

    } else if (initMethod == "medoids") {
      zValue[[iterations]] <- mclust::unmap(cluster::pam(log(dataset + 1 / 3),
                                                         k = numbG,  cluster.only = TRUE))
      # if z generated has less columns than numbG, then use random initialization
      if(ncol(zValue[[iterations]]) < numbG) {
        zValue[[iterations]] <- randomInitfunction(numbG = numbG, nObservations = nObservations)
      }

    } else if (initMethod == "clara") {
      zValue[[iterations]] <- mclust::unmap(cluster::clara(log(dataset + 1 / 3),
                                                           k = numbG)$cluster)
      # if z generated has less columns than numbG, then use random initialization
      if(ncol(zValue[[iterations]]) < numbG) {
        zValue[[iterations]] <- randomInitfunction(numbG = numbG, nObservations = nObservations)
      }

    } else if (initMethod == "fanny") {
      zValue[[iterations]] <- mclust::unmap(cluster::fanny(log(dataset + 1 / 3),
                                                           k = numbG, memb.exp = numbG, cluster.only = TRUE)$clustering)
      # if z generated has less columns than numbG, then use random initialization
      if(ncol(zValue[[iterations]]) < numbG) {
        zValue[[iterations]] <- randomInitfunction(numbG = numbG, nObservations = nObservations)
      }
    }


    initRuns[[iterations]] <- varMPLNInitClustering(dataset = dataset,
                                                    G = numbG,
                                                    zValue = zValue[[iterations]],
                                                    normFactors = normFactors,
                                                    maxIterations = 10)
    # maxIterations set to 10 for initialization

    logLinit[iterations] <-
      unlist(utils::tail((initRuns[[iterations]]$logLikelihood), n = 1))
  }

  # select the initialization run with highest loglikelihood
  initializationResults <- initRuns[[which(logLinit == max(logLinit, na.rm = TRUE))[1]]]

  class(initializationResults) <- "varMPLNInitialization"
  return(initializationResults)
}


varMPLNInitClustering <- function(dataset,
                                  G,
                                  zValue,
                                  normFactors,
                                  maxIterations = 10) {
  # if running for initialization, need to stop after 10 iterations

  dimensionality <- ncol(dataset)
  nObservations <- nrow(dataset)


  # basic initialization performed for parameters
  mu <- sigma <- isigma <- list()
  m <- S <- list() # variational parameters

  ###Other intermediate items initialized
  Sk <- array(0, c(dimensionality, dimensionality, G))
  mPreviousValue <- GX <- dGX <- zSValue <- list()


  piG <- colSums(zValue) / nObservations

  for (g in 1:G) {
    obs <- which(zValue[ , g] == 1)
    mu[[g]] <- colMeans(log(dataset[obs, ] + 1 / 6)) # starting value for mu
    # starting value for sample covariance matrix
    sigma[[g]] <- cov(log(dataset[obs, ] + 1 / 6))
    # starting value for inverse of sample covariance matrix
    # If the inverse is not present for covariance matrix, handle that
    isigma[[g]] <- tryCatch(solve(sigma[[g]]), error = function(err) NA)
    if(all(is.na(isigma[[g]]))) {
      isigma[[g]] <- solve(sigma[[g]] + 0.001) # if error with inverse
    }

  }

  for (g in 1:G) {
    # mPreviousValue = starting value for m
    mPreviousValue[[g]] <- m[[g]] <- log(dataset + 1 / 6)
    S[[g]] <- list() # starting value for S
    for (i in 1:nObservations) {
      S[[g]][[i]] <- diag(dimensionality) * 0.000000001
    }
  }

  # Start clustering
  itOuter <- 1
  aloglik <- logLikelihood <- NULL
  checks <- aloglik[c(1, 2, 3)] <- 0


  while (checks == 0) {


    for (g in 1:G) {
      GX[[g]] <- list()
      dGX[[g]] <- list()
      zSValue[[g]] <- list()
      for (i in 1:nObservations) {
        dGX[[g]][[i]] <- diag(exp((log(normFactors) + mPreviousValue[[g]][i, ]) +
                                    0.5 * diag(S[[g]][[i]])), dimensionality) + isigma[[g]]
        S[[g]][[i]] <- solve(dGX[[g]][[i]]) # update S
        # will be used for updating sample covariance matrix
        zSValue[[g]][[i]] <- zValue[i, g] * S[[g]][[i]]
        GX[[g]][[i]] <- dataset[i, ] - exp(mPreviousValue[[g]][i, ] +
                                             log(normFactors) + 0.5 * diag(S[[g]][[i]])) -
          (isigma[[g]]) %*% (mPreviousValue[[g]][i, ] - mu[[g]])
        m[[g]][i , ] <- mPreviousValue[[g]][i , ] + S[[g]][[i]] %*% GX[[g]][[i]] # update m
      }

      mPreviousValue[[g]] <- m[[g]]

      # Updating mu
      mu[[g]] <- colSums(zValue[ , g] * m[[g]]) / sum(zValue[ , g])

      # Updating sample covariance
      muMatrix <- matrix(rep(mu[[g]], nObservations), nrow = nObservations, byrow = TRUE)
      res <- m[[g]] - muMatrix
      # Calculate weighted covariance
      temp <- stats::cov.wt(res, wt = zValue[ , g], center = FALSE, method = "ML")
      Sk[ , , g] <- temp$cov
      sigma[[g]] <- temp$cov + Reduce("+", zSValue[[g]]) / sum(zValue[ , g])
      isigma[[g]] <- solve(sigma[[g]])
    }


    piG <- colSums(zValue) / nObservations

    # Matrix containing normaization factors
    normFactorsAsMatrix <- matrix(rep(normFactors, nObservations), nrow = nObservations, byrow = TRUE)


    # A useful function
    zvalueFuncTerm <- function(x, y = isigma[[g]]) {
      sum(diag(x %*% y))
    }

    # # Zvalue calculation
    forz <- matrix(NA, ncol = G, nrow = nObservations)

    for (g in 1:G) {
      # exp(m_igj + 0.5 S_ig,jj)
      two <- rowSums(exp(m[[g]] + log(normFactorsAsMatrix) +
                           0.5 * matrix(unlist(lapply(S[[g]], diag)), ncol = dimensionality, byrow = TRUE)))
      five <- 0.5 * unlist(lapply(S[[g]], zvalueFuncTerm)) # trace(isigma x S_ig)
      six <- 0.5 * log(unlist(lapply(S[[g]], det))) # 0.5 * log|S_ig|

      # For zig calculation (the numerator part)
      forz[ , g] <- piG[g] *
        exp(rowSums(m[[g]] * dataset) - two - rowSums(lfactorial(dataset)) +
              rowSums(log(normFactorsAsMatrix) * dataset) -
              0.5 * mahalanobis(m[[g]], center = mu[[g]], cov = sigma[[g]]) -
              five + six - 0.5 * log(det(sigma[[g]])) - dimensionality / 2)
    }

    # Calculate zValue value
    # check which forz == 0 and rowSums(forz)==0 and which of these
    # have both equalling to 0 (because 0/0 =NaN)
    if (G == 1) {
      errorpossible <- Reduce(intersect,
                              list(which(forz == 0),
                                   which(rowSums(forz) == 0)))
      forz[errorpossible] <- 1e-100
      zvalue <- forz / rowSums(forz)
    } else {

      # check for error, if rowsums are zero
      rowSumsZero <- which(rowSums(forz) == 0)
      if(length(rowSumsZero) > 1) {
        forz[rowSumsZero, ] <- mclust::unmap(stats::kmeans(log(dataset + 1 / 6),
                                                           centers = G,
                                                           nstart = 100)$cluster)[rowSumsZero, ]
        zvalue <- forz / rowSums(forz)
      } else {
        zvalue <- forz / rowSums(forz)
      }
    }



    # Calculate log-likelihood
    logLikelihood[itOuter] <- sum(log(rowSums(forz)))


    # Stopping criterion
    if (itOuter > 2) {

      if ((logLikelihood[itOuter - 1] - logLikelihood[itOuter - 2]) == 0) {
        checks <-1
      } else {
        # Aitken stopping criterion
        termAitkens <- (logLikelihood[itOuter]- logLikelihood[itOuter - 1]) /
          (logLikelihood[itOuter - 1] - logLikelihood[itOuter - 2])
        term2Aitkens <- (1 / (1 - termAitkens) * (logLikelihood[itOuter] - logLikelihood[itOuter - 1]))
        aloglik[itOuter] <- logLikelihood[itOuter - 1] + term2Aitkens
        if (abs(aloglik[itOuter] - logLikelihood[itOuter - 1]) < 0.001) {
          # If this critera, as per Böhning et al., 1994 is achieved
          # convergence is achieved
          checks <- 1
        } else {
          checks <- checks}
      }
    }
    # print(itOuter)
    itOuter <- itOuter + 1
    if (itOuter == maxIterations) {
      checks <- 1
    }
  }

  # Naming parameters
  names(mu) <- names(sigma) <- paste0(rep("G=", G), 1:G)

  # Saving results for output
  Results <- list(piG = piG,
                  mu = mu,
                  sigma = sigma,
                  isigma = isigma,
                  zValue = zValue,
                  m = m,
                  S = S,
                  clusterlabels = mclust::map(zValue),
                  logLikelihood = logLikelihood)

  class(Results) <- "varMPLNInitClustering"
  return(Results)
}






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







#' Visualize Clustered Results Via MPLN
#'
#' A function to visualize data and clustering results obtained
#' from a mixtures of multivariate Poisson-log normal (MPLN) model.
#' The function produces heatmaps and line plots of input data.
#' If cluster membership of observations are provided, this
#' information will be indicated in the figures. If a matrix of
#' probabilities for the observations belonging to each cluster
#' is provided, the option to produce a barplot of probabilities
#' is also available.
#'
#'
#' @param dataset A dataset of class matrix and type integer such that
#'    rows correspond to observations and columns correspond to variables.
#' @param plots A character string indicating which plots to be produced.
#'    Options are 'heatmaps', 'lines', 'bar', and 'all'. Default is 'all'.
#' @param probabilities A matrix of size N x C, such that rows correspond
#'    to N observations and columns correspond to C clusters. Each row
#'    should sum to 1. Default is NA.
#' @param clusterMembershipVector A numeric vector of length nrow(dataset)
#'    containing the cluster membership of each observation as generated by
#'    mpln(). Default is NA.
#' @param LinePlotColours Character string indicating if the line plots
#'    should be multicoloured or monotone, in black. Options are
#'    'multicolour' or 'black'. Default is 'black'.
#' @param printPlot Logical indicating if plot(s) should be saved in local
#'    directory. Default TRUE. Options TRUE or FALSE.
#' @param fileName Unique character string indicating the name for the plot
#'    being generated. Default is Plot_date, where date is obtained from
#'    date().
#' @param format Character string indicating the format of the image to
#'    be produced. Default 'pdf'. Options 'pdf' or 'png'.
#'
#' @return Plotting function provides the possibility for line and heatmap
#'    plots.
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
#'  simulatedCounts <- MPLNClust::mplnDataGenerator(nObservations = 70,
#'                                       dimensionality = 6,
#'                                       mixingProportions = c(0.79, 0.21),
#'                                       mu = rbind(trueMu1, trueMu2),
#'                                       sigma = rbind(trueSigma1, trueSigma2),
#'                                       produceImage = "No")
#'
#'  # Clustering data
#'  MPLNClustResults <- MPLNClust::mplnVariational(
#'                               dataset = as.matrix(simulatedCounts$dataset),
#'                               membership = "none",
#'                               gmin = 1,
#'                               gmax = 2,
#'                               initMethod = "kmeans",
#'                               nInitIterations = 1,
#'                               normalize = "Yes")
#'
#'  # Visualize data
#'  MPLNVisuals <- MPLNClust::mplnVisualize(dataset = simulatedCounts$dataset,
#'                                          plots = 'all',
#'                                          probabilities =
#'                                          MPLNClustResults$allResults$`G=2`$probaPost,
#'                                          clusterMembershipVector =
#'                                          MPLNClustResults$allResults$`G=2`$clusterlabels,
#'                                          fileName = 'TwoClusterModel',
#'                                          printPlot = FALSE,
#'                                          format = 'png')
#'
#' @author Anjali Silva, \email{anjali.silva@uhnresearch.ca}
#'
#' @export
#' @import graphics
#' @import ggplot2
#' @importFrom grDevices png
#' @importFrom grDevices pdf
#' @importFrom grDevices dev.off
#' @importFrom RColorBrewer brewer.pal.info
#' @importFrom RColorBrewer brewer.pal
#' @importFrom randomcoloR distinctColorPalette
#' @importFrom pheatmap pheatmap
#' @importFrom gplots heatmap.2
#' @importFrom gplots redgreen
mplnVisualize <- function(dataset,
                          plots = 'all',
                          probabilities = NA,
                          clusterMembershipVector = NA,
                          fileName = paste0('Plot_',date()),
                          LinePlotColours = "black",
                          printPlot = TRUE,
                          format = 'pdf') {

  # Checking user input
  if (typeof(dataset) != "double" & typeof(dataset) != "integer") {
    stop("\n Dataset type needs to be integer")
  }

  if (is.matrix(dataset) != TRUE) {
    stop("\n Dataset needs to be a matrix")
  }

  if (is.character(plots) == TRUE) {
    plotsMethodsPossible<- c("all", "bar", "lines", "heatmaps")
    if(all((plots == plotsMethodsPossible) == FALSE)) {
      stop("initMethod should of class character, specifying
             either: all, bar, lines, heatmaps.")
    }
  } else if (is.character(plots) != TRUE) {
    stop("initMethod should of class character, specifying
             either: all, bar, lines, heatmaps.")
  }

  if (is.logical(clusterMembershipVector) == TRUE) {
    cat("\n clusterMembershipVector is not provided.")
    clusterMembershipVector <- rep(1, nrow(dataset))

  } else if (is.numeric(clusterMembershipVector) == TRUE) {
    if (nrow(dataset) != length(clusterMembershipVector)) {
      stop("\n length(clusterMembershipVector) should match
          nrow(dataset)")
    }
  }



  if (is.logical(probabilities) == TRUE) {
    cat("\n Probabilities are not provided. Barplot of probabilities will not be produced.")
  } else if (is.matrix(probabilities) == TRUE) {
    if (nrow(probabilities) != length(clusterMembershipVector)) {
      stop("\n length(probabilities) should match nrow(dataset)")
    }
    if (any(rowSums(probabilities) >= 1.01)) {
      stop("\n rowSums(probabilities) reveals at least
          one observation has probability != 1.")
    }
    if (any(rowSums(probabilities) <= 0.99)) {
      stop("\n rowSums(probabilities) reveals at least
          one observation has probability != 1.")
    }
  }


  # Obtaining path to save images
  pathNow <- getwd()

  # Saving cluster membership for each observation
  DataPlusLabs <- cbind(dataset, clusterMembershipVector)
  ordervector <- anothervector <- list()

  # Divide observations into each cluster based on membership
  for (i in 1:max(clusterMembershipVector)) {
    ordervector[[i]] <- which(DataPlusLabs[,
                                           ncol(dataset) + 1] == i)
    # divide observations as an integer based on cluster membership
    anothervector[[i]] <- rep(i,
                              length(which(DataPlusLabs[,
                                                        ncol(dataset) + 1] == i)))
  }

  vec <- unlist(ordervector) # put observations in order of cluster membership
  colorsvector <- unlist(anothervector) # put all details together as integers

  # Setting the colours
  if(max(clusterMembershipVector) > 17) {
    qualColPals <- RColorBrewer::brewer.pal.info[brewer.pal.info$category == 'qual', ]
    coloursBarPlot <- unlist(mapply(RColorBrewer::brewer.pal,
                                    qualColPals$maxcolors,
                                    rownames(qualColPals)))
  } else {
    coloursBarPlot <- c('#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6',
                        '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324',
                        '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1',
                        '#000075', '#808080')
  }


  # empty plots
  heatmapOne <- heatmapTwo <- linePlots <- barPlot <- NULL

  if (plots == 'all' || plots == 'heatmaps') {

    # Heatmap 1
    if (printPlot == TRUE) {

      if (format == 'png') {
        grDevices::png(paste0(pathNow, "/heatmap1_", fileName, ".png"))
      } else {
        grDevices::pdf(paste0(pathNow, "/heatmap1_", fileName, ".pdf"))
      }

      gplots::heatmap.2(as.matrix(dataset[vec, ]),
                        dendrogram = "column",
                        trace = "none",
                        scale = "row",
                        Rowv = FALSE, col = rev(gplots::redgreen(75)),
                        RowSideColor = coloursBarPlot[colorsvector],
                        labRow = FALSE,
                        main = paste("Clustering results, G =",
                                     max(clusterMembershipVector)))
      graphics::par(xpd = TRUE)
      graphics::legend(xpd = TRUE, x = -0.1, y = 0,
                       legend = paste0("Cluster ", unique(colorsvector)),
                       col = unique(coloursBarPlot[colorsvector]),
                       lty = 1,
                       lwd = 5,
                       cex =.8, horiz = FALSE)

      grDevices::dev.off()
    }


    # Heatmap 2
    # Only produced if less than 17 clusters
    if(max(clusterMembershipVector) < 18) {
      # Defining annotation row
      annotation_row = data.frame(Cluster = factor(clusterMembershipVector[vec]))
      if(is.null(rownames(dataset)) == TRUE) {
        rownames(dataset)  = paste("Gene", c(1:nrow(dataset[vec, ])))
        rownames(annotation_row) = rownames(dataset[vec, ])
      } else {
        rownames(annotation_row) = rownames(dataset[vec, ])
      }

      # Define row annotation colours
      heatMap2RowAnnotation <- c("1" = coloursBarPlot[1], "2" = coloursBarPlot[2],
                                 "3" = coloursBarPlot[3], "4" = coloursBarPlot[4],
                                 "5" = coloursBarPlot[5], "6" = coloursBarPlot[6],
                                 "7" = coloursBarPlot[7], "8" = coloursBarPlot[8],
                                 "9" = coloursBarPlot[9], "10" = coloursBarPlot[10],
                                 "11" = coloursBarPlot[11], "12" = coloursBarPlot[12],
                                 "13" = coloursBarPlot[13], "14" = coloursBarPlot[14],
                                 "15" = coloursBarPlot[15], "16" = coloursBarPlot[16],
                                 "17" = coloursBarPlot[17])

      # Show row names or not based on dataset size
      if(nrow(dataset) < 50){
        showLabels = TRUE
      } else {
        showLabels = FALSE
      }


      if (printPlot == TRUE) {
        if (format == 'png') {
          grDevices::png(paste0(pathNow, "/heatmap2_", fileName, ".png"))
        } else {
          grDevices::pdf(paste0(pathNow, "/heatmap2_", fileName, ".pdf"))
        }

        heatmapFunctionTwo(dataset = dataset,
                           vec = vec,
                           showLabels = showLabels,
                           heatMap2RowAnnotation = heatMap2RowAnnotation,
                           annotation_row = annotation_row,
                           clusterMembershipVector = clusterMembershipVector)
        grDevices::dev.off()
      }


      heatmapTwo <- heatmapFunctionTwo(dataset = dataset,
                                       vec = vec,
                                       showLabels = showLabels,
                                       heatMap2RowAnnotation = heatMap2RowAnnotation,
                                       annotation_row = annotation_row,
                                       clusterMembershipVector = clusterMembershipVector)
    }
  }

  if (plots == 'all' || plots == 'lines') {
    # Line Plots

    if (LinePlotColours == "multicolour") {
      linePlots <- list()
      for(cluster in unique(clusterMembershipVector)) {

        # Save how many observations below to each cluster size,
        # given by 'cluster'
        if (length(which(DataPlusLabs[, ncol(dataset) + 1] == cluster)) == 1) {
          toPlot2 <- as.matrix(DataPlusLabs[which(DataPlusLabs[,
                                                                ncol(dataset) + 1] == cluster), c(1:ncol(dataset))],
                                ncol = ncol(dataset))
          rownames(toPlot2) <- names(which(DataPlusLabs[, ncol(dataset) + 1] == cluster))
        } else if (length(which(DataPlusLabs[, ncol(dataset) + 1] == cluster)) > 1) {
          toPlot2 <- as.matrix(DataPlusLabs[which(DataPlusLabs[,
                                                                ncol(dataset) + 1] == cluster), c(1:ncol(dataset))],
                                ncol = ncol(dataset))
        }

        # Save column mean in last row
        toplot1 <- rbind(log(toPlot2 + 1), colMeans(log(toPlot2 + 1)))
        # If discontinunity is needed between samples (e.g. for 6 samples)
        # toplot1_space=cbind(toplot1[,c(1:3)],rep(NA,nrow(toPlot2)+1),
        # toplot1[,c(4:6)])


        if (printPlot == TRUE) {
          if (format == 'png') {
            grDevices::png(paste0(pathNow, "/LinePlot_Cluster", cluster,
                                  "_", fileName, ".png"))
          } else {
            grDevices::pdf(paste0(pathNow, "/LinePlot_Cluster", cluster,
                                  "_", fileName, ".pdf"))
          }

          linePlotMultiCol(dataset = dataset,
                           toplot1 = toplot1,
                           toPlot2 = toPlot2,
                           coloursBarPlot = coloursBarPlot,
                           cluster = cluster)
          grDevices::dev.off()
        }

        linePlots[[cluster]] <- linePlotMultiCol(dataset = dataset,
                                                 toplot1 = toplot1,
                                                 toPlot2 = toPlot2,
                                                 coloursBarPlot = coloursBarPlot,
                                                 cluster = cluster)
      }
    } else if (LinePlotColours == "black") {
      linePlots <- list()
      for(cluster in unique(clusterMembershipVector)) {


        # Save how many observations below to each cluster size,
        # given by 'cluster'
        if (length(which(DataPlusLabs[, ncol(dataset) + 1] == cluster)) == 1) {
          toPlot2 <- t(as.matrix(DataPlusLabs[which(DataPlusLabs[,
                                                                  ncol(dataset) + 1] == cluster), c(1:ncol(dataset))],
                                  ncol = ncol(dataset)))
          rownames(toPlot2) <- names(which(DataPlusLabs[, ncol(dataset) + 1] == cluster))
        } else if (length(which(DataPlusLabs[, ncol(dataset) + 1] == cluster)) > 1) {
          toPlot2 <- as.matrix(DataPlusLabs[which(DataPlusLabs[,
                                                                ncol(dataset) + 1] == cluster), c(1:ncol(dataset))],
                                ncol = ncol(dataset))
        }

        # Save column mean in last row
        toplot1 <- rbind(log(toPlot2 + 1), colMeans(log(toPlot2 + 1)))
        # If discontinunity is needed between samples (e.g. for 6 samples)
        # toplot1_space=cbind(toplot1[,c(1:3)],rep(NA,nrow(toPlot2)+1),
        # toplot1[,c(4:6)])

        if (printPlot == TRUE) {
          if (format == 'png') {
            grDevices::png(paste0(pathNow, "/LinePlot_Cluster", cluster,
                                  "_", fileName, ".png"))
          } else {
            grDevices::pdf(paste0(pathNow, "/LinePlot_Cluster", cluster,
                                  "_", fileName, ".pdf"))
          }
          linePlotMonoCol(dataset = dataset,
                          toplot1 = toplot1,
                          toPlot2 = toPlot2,
                          cluster = cluster)
          grDevices::dev.off()
        }

        linePlots[[cluster]] <- linePlotMonoCol(dataset = dataset,
                                                toplot1 = toplot1,
                                                toPlot2 = toPlot2,
                                                cluster = cluster)
      }
    }
  }

  if (plots == 'all' || plots == 'bar') {

    if(is.logical(probabilities) == TRUE){
      stop("\n probabilities should be provided to make bar plot.")
    }

    # Bar plot
    tableProbabilities <- as.data.frame(cbind(Sample = c(1:nrow(probabilities)),
                                              Cluster = mclust::map(probabilities),
                                              probabilities))

    names(tableProbabilities) <- c("Sample", "Cluster",
                                   paste0("P", rep(1:(ncol(tableProbabilities)-2))))

    tableProbabilitiesMelt <- reshape::melt(tableProbabilities,
                                            id.vars = c("Sample","Cluster"))

    if (printPlot == TRUE) {
      barPlot <- barPlotFunction(tableProbabilitiesMelt = tableProbabilitiesMelt,
                                 coloursBarPlot = coloursBarPlot,
                                 probabilities = probabilities)
      ggplot2::ggsave(paste0(pathNow,"/barplot_", fileName,".",format))
    }

    barPlot <- barPlotFunction(tableProbabilitiesMelt = tableProbabilitiesMelt,
                               coloursBarPlot = coloursBarPlot,
                               probabilities = probabilities)
  }

  return(list(heatmapOne,
              heatmapTwo,
              linePlots,
              barPlot))
}




heatmapFunctionTwo <- function(dataset,
                               vec,
                               showLabels,
                               heatMap2RowAnnotation,
                               annotation_row,
                               clusterMembershipVector) {
  pheatmapPlot <- pheatmap::pheatmap(as.matrix(dataset[vec, ]), show_colnames = TRUE,
                                     show_rownames = showLabels,
                                     labels_col = colnames(dataset),
                                     annotation_row = annotation_row,
                                     annotation_colors = list(Cluster = heatMap2RowAnnotation[
                                       sort(unique(clusterMembershipVector))]),
                                     fontface = "italic", legend = TRUE, scale ="row",
                                     border_color = "black", cluster_row = FALSE,
                                     cluster_col = FALSE,
                                     color =  rev(gplots::redgreen(1000)) )
  return(pheatmapPlot)
}


linePlotMultiCol <- function(dataset,
                             toplot1,
                             toPlot2,
                             coloursBarPlot,
                             cluster) {
  linePlotMultiCol <- graphics::matplot(t(toplot1), type = "l", pch = 1,
                                        col = c(rep(coloursBarPlot[cluster], nrow(toPlot2)), 7),
                                        xlab = "Samples", ylab = "Expression (log counts)", cex = 1,
                                        lty = c(rep(2, nrow(toPlot2)), 1),
                                        lwd = c(rep(3, nrow(toPlot2)), 4),
                                        xaxt = "n", xlim = c(1, ncol(toplot1)),
                                        main = paste("Cluster ", cluster))
  linePlotMultiCol <- linePlotMultiCol + axis(1, at = c(1:ncol(dataset)), labels = colnames(dataset))
  return(linePlotMultiCol)
}


linePlotMonoCol <- function(dataset,
                            toplot1,
                            toPlot2,
                            cluster) {
  linePlotMonoCol <- graphics::matplot(t(toplot1), type = "l", pch = 1,
                                       col = c(rep(1, nrow(toPlot2)), 7),
                                       xlab = "Samples", ylab = "Expression (log counts)", cex = 1,
                                       lty = c(rep(2, nrow(toPlot2)), 1),
                                       lwd = c(rep(3, nrow(toPlot2)), 4),
                                       xaxt = "n", xlim = c(1, ncol(toplot1)),
                                       main = paste("Cluster ", cluster))
  linePlotMonoCol <- linePlotMonoCol + axis(1, at = c(1:ncol(dataset)), labels = colnames(dataset))
  return(linePlotMonoCol)
}


barPlotFunction <- function(tableProbabilitiesMelt,
                            coloursBarPlot,
                            probabilities) {

  variable <- value <- Sample <- NULL

  if(is.data.frame(tableProbabilitiesMelt) != TRUE) {
    stop("tableProbabilitiesMelt should be a data frame")
  }

  if(is.character(coloursBarPlot) != TRUE) {
    stop("coloursBarPlot should be character")
  }

  if(is.matrix(probabilities) != TRUE) {
    stop("probabilities should be a matrix")
  }

  barPlot <- ggplot2::ggplot(data = tableProbabilitiesMelt,
                             ggplot2::aes(fill = variable, y = value, x = Sample))

  barPlot <- barPlot + ggplot2::geom_bar(position = "fill", stat = "identity") +
    scale_fill_manual(values = coloursBarPlot,
                      name = "Cluster") + theme_bw() +
    theme(text = element_text(size = 10),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(face = "bold"),
          axis.text.y = element_text(face = "bold")) +
    coord_cartesian(ylim = c(0, 1), xlim = c(1, nrow(probabilities))) +
    labs(x = "Observation") +
    scale_y_continuous(name = "Posterior probability", limits = c(0: 1))
  return(barPlot)
}







#' Generating Data Using Mixtures of MPLN
#'
#' This function simulates data from a mixture of MPLN model.
#'
#' @param nObservations A positive integer indicating the number of observations for
#'    the dataset.
#' @param dimensionality A positive integer indicating the dimensionality for the
#'    dataset.
#' @param mixingProportions A numeric vector that length equal to the number of total
#'    components, indicating the proportion of each component. Vector content should
#'    sum to 1.
#' @param mu A matrix of size (dimensionality x number of components), indicating
#'    the mean for each component. See example.
#' @param sigma A matrix of size ((dimensionality * number of components) x
#'    dimensionality), indicating the covariance matrix for each component.
#'    See example.
#' @param produceImage A character string indicating whether or not to
#'    produce an image. Options "Yes" or "No". Image will be produced as
#'    'Pairs plot of log-transformed data.png" in the current working
#'    directory.
#' @param ImageName A character string indicating name for image, if
#'    produceImage is set to "Yes". Default is "TwoComponents".
#'
#' @return Returns an S3 object of class mplnDataGenerator with results.
#' \itemize{
#'   \item dataset - Simulated dataset.
#'   \item trueMembership -A numeric vector indicating the membership of
#'      each observation.
#'   \item probaPost - A matrix indicating the posterior probability that
#'      each observation belong to the component/cluster.
#'   \item truenormfactors - A numeric vector indicating the true
#'      normalization factors used for adjusting the library sizes.
#'   \item observations - Number of observations in the simulated dataset.
#'   \item dimensionality - Dimensionality of the simulated dataset.
#'   \item mixingProportions - A numeric vector indicating the mixing
#'      proportion of each component.
#'   \item mu - True mean used for the simulated dataset.
#'   \item sigma - True covariances used for the simulated dataset.
#'}
#'
#' @examples
#' # Generating simulated data
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
#'                                            produceImage = "No",
#'                                            ImageName = "TwoComponents")
#'
#' @author Anjali Silva, \email{anjali.silva@uhnresearch.ca}
#'
#' @references
#' Aitchison, J. and C. H. Ho (1989). The multivariate Poisson-log normal distribution.
#' \emph{Biometrika} 76.
#'
#' Silva, A. et al. (2019). A multivariate Poisson-log normal mixture model
#' for clustering transcriptome sequencing data. \emph{BMC Bioinformatics} 20.
#' \href{https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-2916-0}{Link}
#'
#' @export
#' @import stats
#' @importFrom mvtnorm rmvnorm
#' @importFrom edgeR calcNormFactors
#' @importFrom grDevices dev.off
#' @importFrom grDevices png
#' @importFrom graphics pairs
mplnDataGenerator <- function(nObservations,
                              dimensionality,
                              mixingProportions,
                              mu,
                              sigma,
                              produceImage = "No",
                              ImageName = "sampleName") {

  # Checking user input
  if(is.numeric(nObservations) != TRUE) {
    stop("nObservations argument should be of class numeric.")
  }

  if(is.numeric(dimensionality) != TRUE) {
    stop("dimensionality argument should be of class numeric.")
  }

  if(is.numeric(mixingProportions) != TRUE) {
    stop("mixingProportions argument should be a vector of class numeric.")
  }

  if (sum(mixingProportions) != 1) {
    stop("mixingProportions argument should be a vector that sum to 1.")
  }

  if(is.matrix(mu) != TRUE) {
    stop("mu argument should be of class matrix.")
  }

  if(ncol(mu) != dimensionality) {
    stop("mu should be a matrix, which has number of columns equalling dimensionality.")
  }

  if(nrow(mu) != length(mixingProportions)) {
    stop("mu should be a matrix, which has number of rows equalling number of components.")
  }

  if(is.matrix(sigma) != TRUE) {
    stop("sigma argument should be of class matrix.")
  }

  if(ncol(sigma) != dimensionality) {
    stop("sigma should be a matrix, which has number of columns equalling
      dimensionality.")
  }

  if(nrow(sigma) != dimensionality * length(mixingProportions)) {
    stop("sigma should be a matrix, which has number of rows equalling
      (dimensionality * number of components).")
  }

  if (produceImage == "Yes" && is.character(ImageName) != TRUE) {
    stop("ImageName should be a character string of class character.")
  }

  # Begin calculations - generate z
  z <- t(stats::rmultinom(nObservations, size = 1, mixingProportions))

  y <- theta <- nG <- vector("list", length = length(mixingProportions))

  # For visualization only
  theta2 <- matrix(NA, ncol = dimensionality, nrow = nObservations)

  for (i in seq_along(1:length(mixingProportions))) {
    nG[[i]] <- which(z[, i] == 1)
    theta[[i]] <- mvtnorm::rmvnorm(n = length(nG[[i]]),
                                   mean = mu[i, ],
                                   sigma = sigma[((i - 1) *
                                                    dimensionality + 1):(i * dimensionality), ])
    theta2[nG[[i]], ] <- mvtnorm::rmvnorm(n = length(nG[[i]]),
                                          mean = mu[i, ],
                                          sigma = sigma[((i  -1) *
                                                           dimensionality + 1):(i * dimensionality), ])
  }

  y <- matrix(NA, ncol = dimensionality, nrow = nObservations)
  for (i in seq_along(1:nObservations)) {
    for (j in seq_along(1:dimensionality)) {
      y[i, j] <- stats::rpois(1, exp(theta2[i, j]))
    }
  }

  norms <- log(edgeR::calcNormFactors(y))

  # Generating counts with norm factors
  y2 <- matrix(NA, ncol = dimensionality, nrow = nObservations)
  for (i in seq_along(1:nObservations)) {
    for (j in seq_along(1:dimensionality)) {
      y2[i, j] <- stats::rpois(n = 1, lambda = exp(theta2[i, j] + norms[j]))
    }
  }

  if (produceImage == "Yes") {
    # Obtaining path to save images
    pathNow <- getwd()
    grDevices::png(paste0(pathNow, "/PairsPlot_", ImageName,".png"))
    graphics::pairs(log(y2 + 1 / 100), col = mclust::map(z) + 1,
                    main = "Pairs plot of log-transformed data")
    grDevices::dev.off()
  }

  results <- list(dataset = y2,
                  trueMembership = mclust::map(z),
                  probaPost = z,
                  trueNormFactors = norms,
                  observations = nObservations,
                  dimensionality = dimensionality,
                  mixingProportions = mixingProportions,
                  trueMu = mu,
                  trueSigma = sigma)

  class(results) <- "mplnDataGenerator"
  return(results)
}

# [END]
