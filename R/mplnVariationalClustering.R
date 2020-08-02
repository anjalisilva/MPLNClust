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
#' # sampleData <- mplnDataGenerator(nObservations = 1000,
#' #                                 dimensionality = 6,
#' #                                 mixingProportions = c(0.79, 0.21),
#' #                                 mu = rbind(trueMu1, trueMu2),
#' #                                 sigma = rbind(trueSigma1, trueSigma2),
#' #                                 produceImage = "No")
#'
#' # Clustering
#' # mplnResults <- mplnVariational(dataset = sampleData$dataset,
#' #                                 membership = sampleData$trueMembership,
#' #                                 gmin = 1,
#' #                                 gmax = 2,
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

  if (class(initMethod) == "character") {
    initMethodsUsed <- c("kmeans", "random", "medoids", "clara", "fanny")
    if(all((initMethod == initMethodsUsed) == FALSE)) {
      stop("initMethod should of class character, specifying
        either: kmeans, random, medoids, clara, or fanny.")
    }
  } else if ((class(initMethod)) != "character") {
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
    logLikelihood[g] <- unlist(utils::tail(clusterResults[[g]]$logLikelihood, n = 1))

    nParameters[g] <- calcParameters(numberG = g, dimensionality = dimensionality)

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
      mu[[g]] <- colMeans(log(dataset[obs, ] + 1 / 6)) # starting value for mu
      sigma[[g]] <- var(log(dataset[obs, ] + 1 / 6)) # starting value for sample covariance matrix
      isigma[[g]] <- solve(sigma[[g]])
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
      zvalue <- forz / rowSums(forz)
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
        if (abs(aloglik[itOuter] - logLikelihood[itOuter - 1]) < 0.01) {
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
    isigma[[g]] <- solve(sigma[[g]])
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
      zvalue <- forz / rowSums(forz)
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
        if (abs(aloglik[itOuter] - logLikelihood[itOuter - 1]) < 0.01) {
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


# [END]
