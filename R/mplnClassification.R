#' Classification Using MPLN Via Variational-EM
#'
#' Performs classification using mixtures of multivariate Poisson-log
#' normal (MPLN) distribution with variational expectation-maximization
#' (EM) for parameter estimation. Model selection is performed using
#' AIC, AIC3, BIC and ICL. No internal parallelization, thus code is
#' run in serial.
#'
#' @param dataset A dataset of class matrix and type integer such that
#'    rows correspond to observations and columns correspond to variables.
#'    The dataset have dimensions n x d, where n is the total number of
#'    observations and d is the dimensionality. If rowSums are zero, these
#'    rows will be removed prior to cluster analysis.
#' @param membership A numeric vector of length nrow(dataset) containing the
#'    cluster membership of each observation. For observations with unknown
#'    membership, assign value zero. E.g., for a dataset with 10 observations,
#'    and 2 known groups, c(1, 1, 1, 2, 2, 0, 0, 0, 1, 2).
#' @param gmin A positive integer specifying the minimum number of components
#'    to be considered in the clustering run. The value should be equal to
#'    or greater than max(membership).
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
#' # Example 1
#' trueMu1 <- c(6.5, 6, 6, 6, 6, 6)
#' trueMu2 <- c(2, 2.5, 2, 2, 2, 2)
#'
#' trueSigma1 <- diag(6) * 2
#' trueSigma2 <- diag(6)
#'
#' # Generating simulated data
#' sampleData <- MPLNClust::mplnDataGenerator(nObservations = 1000,
#'                      dimensionality = 6,
#'                      mixingProportions = c(0.79, 0.21),
#'                      mu = rbind(trueMu1, trueMu2),
#'                      sigma = rbind(trueSigma1, trueSigma2),
#'                      produceImage = "No")
#'
#' # Classification
#' membershipInfo <- sampleData$trueMembership
#' length(membershipInfo) # 1000 observations
#'
#' # Assume membership of 200 of 1000 observations were unknown
#' set.seed(1234)
#' randomNumb <- sample(1:length(membershipInfo), 200, replace = FALSE)
#' membershipInfo[randomNumb] <- 0
#' table(membershipInfo)
#' #   0   1   2
#' # 200 593 207
#'
#' # Run for g = 2:3 groups
#' # mplnClassificationResults <- MPLNClust::mplnVarClassification(
#' #                                dataset = sampleData$dataset,
#' #                                membership = membershipInfo,
#' #                                gmin = 2,
#' #                                gmax = 3,
#' #                                initMethod = "kmeans",
#' #                                nInitIterations = 2,
#' #                                normalize = "Yes")
#' #names(mplnClassificationResults)
#'
#'
#'
#'
#'
#' # Example 2
#' # Use an external dataset
#' if (requireNamespace("MBCluster.Seq", quietly = TRUE)) {
#' library(MBCluster.Seq)
#' data("Count")
#' dim(Count) # 1000    8
#'
#' # Clustering subset of data
#' subsetCountData <- as.matrix(Count[c(1:500), ])
#' mplnResultsMBClust <- MPLNClust::mplnVariational(
#'                             dataset = subsetCountData,
#'                             membership = "none",
#'                             gmin = 1,
#'                             gmax = 3,
#'                             initMethod = "kmeans",
#'                             nInitIterations = 2,
#'                             normalize = "Yes")
#' names(mplnResultsMBClust)
#'
#'
#' # Classification
#' # Using 500 labels from clustering above, classify rest of 500 observations
#' membershipInfo <- c(mplnResultsMBClust$BICresults$BICmodelSelectedLabels,
#'                     rep(0, 500))
#' max(mplnResultsMBClust$BICresults$BICmodelSelectedLabels) # 2
#' # Assume there maybe 2 to 4 underlying groups
#' # mplnClassificationResults <- MPLNClust::mplnVarClassification(
#' #                               dataset = Count,
#' #                               membership = membershipInfo,
#' #                               gmin = 1,
#' #                               gmax = 4,
#' #                               initMethod = "kmeans",
#' #                               nInitIterations = 2,
#' #                               normalize = "Yes")
#' # names(mplnClassificationResults)
#'
#'
#'}
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
#' @importFrom edgeR calcNormFactors
#' @importFrom mclust unmap
#' @importFrom mclust map
#' @import stats
#' @import cluster
#'
mplnVarClassification <- function(dataset,
                                  membership,
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

        if(is.numeric(membership) != TRUE) {
          stop("membership should be a numeric vector containing the
            known cluster memberships and unknowns as 0s. E.g., for a
            dataset with 10 observations and 2 known groups,
            c(1, 1, 1, 2, 2, 0, 0, 0, 1, 2).")
        }

        if(max(membership) > gmin) {
          stop("The value of gmin cannot be less than max(membership).")
        }

        if(all((diff(sort(unique(membership))) == 1) != TRUE) ) {
          stop("Cluster memberships in the membership vector
            are missing a cluster. E.g., 0, 1, 3, 4, is missing cluster 2.")
        }

        if(length(membership) != nObservations) {
          stop("membership should be a numeric vector, where length(membership)
            should equal the number of observations. E.g., for a
            dataset with 10 observations and 2 known groups,
            c(1, 1, 1, 2, 2, 0, 0, 0, 1, 2).")
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

        classificationResults <- list() # to save classification output
        for (gmodel in seq_along(1:(gmax - gmin + 1))) {

          if(length(1:(gmax - gmin + 1)) == gmax) {
            clustersize <- gmodel
          } else if(length(1:(gmax - gmin + 1)) < gmax) {
            clustersize <- seq(gmin, gmax, 1)[gmodel]
          }
          cat("\n Running for g =", clustersize)
          classificationResults[[gmodel]] <- varMPLNClassification(
                                                        dataset = dataset,
                                                        initMethod = initMethod,
                                                        nInitIterations = nInitIterations,
                                                        membership = membership,
                                                        G = clustersize,
                                                        maxIterations = 1000,
                                                        normFactors = normFactors)
          # cat("\n Legnth of classificationResults is ", length(classificationResults))
        }

        names(classificationResults) <- paste0(rep("G=", length(seq(gmin, gmax, 1))),
                                               seq(gmin, gmax, 1))

        BIC <- ICL <- AIC <- AIC3 <- Djump <- DDSE <- nParameters <- logLikelihood <- vector()

        for(g in seq_along(1:(gmax - gmin + 1))) {
          # save the final log-likelihood

          if(length(1:(gmax - gmin + 1)) == gmax) {
            clustersize <- g
          } else if(length(1:(gmax - gmin + 1)) < gmax) {
            clustersize <- seq(gmin, gmax, 1)[g]
          }

          if(all(is.na(classificationResults[[g]])) != TRUE) {
            logLikelihood[g] <- unlist(utils::tail(classificationResults[[g]]$logLikelihood, n = 1))
          } else if(all(is.na(classificationResults[[g]])) == TRUE){
            logLikelihood[g] <- NA
          }

          if(all(is.na(classificationResults[[g]])) != TRUE) {
            nParameters[g] <- calcParameters(numberG = clustersize, dimensionality = dimensionality)
          } else if(all(is.na(classificationResults[[g]])) == TRUE){
            nParameters[g] <- NA
          }

          if (g == max(1:(gmax - gmin + 1))) { # starting model selection
            bic <- BICFunction(logLikelihood = logLikelihood,
                               nParameters = nParameters,
                               nObservations = nObservations,
                               clusterRunOutput = classificationResults,
                               gmin = gmin,
                               gmax = gmax,
                               parallel = FALSE)
            # cat("\n Done BIC")

            icl <- ICLFunction(logLikelihood = logLikelihood,
                               nParameters = nParameters,
                               nObservations = nObservations,
                               gmin = gmin,
                               gmax = gmax,
                               clusterRunOutput = classificationResults,
                               parallel = FALSE)
            # cat("\n Done ICL")

            aic <- AICFunction(logLikelihood = logLikelihood,
                               nParameters = nParameters,
                               clusterRunOutput = classificationResults,
                               gmin = gmin,
                               gmax = gmax,
                               parallel = FALSE)
            # cat("\n Done AIC")

            aic3 <- AIC3Function(logLikelihood = logLikelihood,
                                 nParameters = nParameters,
                                 clusterRunOutput = classificationResults,
                                 gmin = gmin,
                                 gmax = gmax,
                                 parallel = FALSE)
            # cat("\n Done AIC3")
          }
        }



        finalTime <- proc.time() - initialTime

        RESULTS <- list(dataset = dataset,
                          dimensionality = dimensionality,
                          normalizationFactors = normFactors,
                          gmin = gmin,
                          gmax = gmax,
                          initalizationMethod = initMethod,
                          allResults = classificationResults,
                          logLikelihood = logLikelihood,
                          numbParameters = nParameters,
                          trueLabels = membership,
                          ICLresults = icl,
                          BICresults = bic,
                          AICresults = aic,
                          AIC3results = aic3,
                          slopeHeuristics = "Not used",
                          totalTime = finalTime)


        class(RESULTS) <- "mplnVarClassification"
        return(RESULTS)
 }




varMPLNClassification <- function(dataset,
                                  initMethod,
                                  nInitIterations,
                                  membership,
                                  G,
                                  maxIterations,
                                  normFactors) {

  nObservations <- nrow(dataset)
  dimensionality <- ncol(dataset)

  # If no intialization is requested by user
  if (nInitIterations == 0) {

    # basic initialization performed for parameters
    mu <- sigma <- isigma <- list()
    mi <- mj <- Si <- Sj <- list() # variational parameters

    ###Other intermediate items initialized
    Ski <- Skj <- array(0, c(dimensionality, dimensionality, G))
    miPreviousValue <- mjPreviousValue <- list()
    GXi <- GXj <- dGXi <- dGXj <- zSiValue <- zSjValue <- list()


    # Initialize z
    classIndicator  <- (membership == 0)
    zValue <- matrix(0, nObservations, G)
    for (i in 1:nObservations) {
      # zjg (for non-labeled) generate values randomly
      if(classIndicator[i]) {zValue[i, sample(1:G, 1, replace=FALSE)] <- 1 } # zjg (for non-labeled)
      else{zValue[i, membership[i]] <- 1} # zig (for labeled)
    }

    nObservationsi <- length(which(membership != 0)) # labeled observations
    nObservationsj <- length(which(membership == 0)) # unlabeled observations
    zValuei <- zValue[which(membership != 0), ] # labeled observations
    zValuej <- zValue[which(membership == 0), ] # unlabeled observations

    piG <- (colSums(zValuei) + colSums(zValuej)) / nObservations

    for (g in 1:G) {
      obs <- which(zValue[ , g] == 1)
      if(length(obs) > 1) {
        # Initialize mu
        mu[[g]] <- colMeans(log(dataset[obs, ] + 1 / 6))
        # Initialize sample covariance matrix
        sigma[[g]] <- cov(log(dataset[obs, ] + 1 / 6))
        # Initialize inverse of sample covariance matrix
        # If the inverse is not present for covariance matrix, handle that
        isigma[[g]] <- tryCatch(solve(sigma[[g]]), error = function(err) NA)
        if(all(is.na(isigma[[g]]))) {
          isigma[[g]] <- diag(ncol(dataset[obs, ])) # if error with inverse
        }
      } else if(length(obs) == 1) {
        # Initialize mu
        mu[[g]] <- log(dataset[obs, ] + 1 / 6)
        # Initialize sample covariance matrix
        sigma[[g]] <- diag(ncol(dataset))
        # Initialize inverse of sample covariance matrix
        # If the inverse is not present for covariance matrix, handle that
        isigma[[g]] <- tryCatch(solve(sigma[[g]]), error = function(err) NA)
        if(all(is.na(isigma[[g]]))) {
          isigma[[g]] <- diag(ncol(dataset[obs, ])) # if error with inverse
        }
      }
    }


    for (g in 1:G) {
      # mPreviousValue = starting value for m
      miPreviousValue[[g]] <- mi[[g]] <- log(dataset[which(membership != 0), ] + 1 / 6)
      mjPreviousValue[[g]] <- mj[[g]] <- log(dataset[which(membership == 0), ] + 1 / 6)

      # starting value for S
      Si[[g]] <- Sj[[g]] <- list()
      for (i in 1:nObservationsi) {
        Si[[g]][[i]] <- diag(dimensionality) * 0.000000001
      }
      for (i in 1:nObservationsj) {
        Sj[[g]][[i]] <-  diag(dimensionality) * 0.000000001
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

    miPreviousValue <- mi <- initializationResults$mi # variational parameters
    mjPreviousValue <- mj <- initializationResults$mj # variational parameters
    Si <- initializationResults$Si # variational parameters
    Sj <- initializationResults$Sj # variational parameters
    zValuei <- initializationResults$zValuei
    zValuej <- initializationResults$zValuej

    nObservationsi <- length(which(membership != 0)) # labeled observations
    nObservationsj <- length(which(membership == 0)) # unlabeled observations

    piG <- (colSums(zValuei) + colSums(zValuej)) / nObservations

    ###Other intermediate items initialized
    Ski <- Skj <- array(0, c(dimensionality, dimensionality, G))
    GXi <- GXj <- dGXi <- dGXj <- zSiValue <- zSjValue <- list()
  }






  # Start clustering
  itOuter <- 1
  aloglik <- logLikelihood <- NULL
  checks <- aloglik[c(1, 2, 3)] <- 0


  while (checks == 0) {


    for (g in 1:G) {
      GXi[[g]] <- GXj[[g]] <- list()
      dGXi[[g]] <- dGXj[[g]] <- list()
      zSiValue[[g]] <- zSjValue[[g]] <- list()

      # for i
      for (i in 1:nObservationsi) {
        dGXi[[g]][[i]] <- diag(exp((log(normFactors) + miPreviousValue[[g]][i, ]) +
                                    0.5 * diag(Si[[g]][[i]])), dimensionality) + isigma[[g]]
        # update S; will be used for updating sample covariance matrix
        Si[[g]][[i]] <- tryCatch(solve(dGXi[[g]][[i]]), error = function(err) NA)
         if(all(is.na(Si[[g]][[i]]))) {
           cat("\n", i)
           Si[[g]][[i]] <- diag(ncol(dGXi[[g]][[i]])) # if error with inverse
         } else {
           Si[[g]][[i]] <- solve(dGXi[[g]][[i]]) # update Si
         }
        zSiValue[[g]][[i]] <- zValuei[i, g] * Si[[g]][[i]] # for updating sample covariance matrix
        GXi[[g]][[i]] <- dataset[i, ] - exp(miPreviousValue[[g]][i, ] +
                                             log(normFactors) + 0.5 * diag(Si[[g]][[i]])) -
                                             (isigma[[g]]) %*% (miPreviousValue[[g]][i, ] - mu[[g]])
        mi[[g]][i , ] <- miPreviousValue[[g]][i , ] + Si[[g]][[i]] %*% GXi[[g]][[i]] # update mi
      }

      # for j
      for (i in 1:nObservationsj) {
        dGXj[[g]][[i]] <- diag(exp((log(normFactors) + mjPreviousValue[[g]][i, ]) +
                                     0.5 * diag(Sj[[g]][[i]])), dimensionality) + isigma[[g]]
        # update S; will be used for updating sample covariance matrix
        Sj[[g]][[i]] <- tryCatch(solve(dGXj[[g]][[i]]), error = function(err) NA)
         if(all(is.na(Sj[[g]][[i]]))) {
           cat("\n", i)
           Sj[[g]][[i]] <- diag(ncol(dGXj[[g]][[i]])) # if error with inverse
         } else {
           Sj[[g]][[i]] <- solve(dGXj[[g]][[i]]) # update Sj
         }
        zSjValue[[g]][[i]] <- zValuej[i, g] * Sj[[g]][[i]] # for updating sample covariance matrix
        GXj[[g]][[i]] <- dataset[i, ] - exp(mjPreviousValue[[g]][i, ] +
                                              log(normFactors) + 0.5 * diag(Sj[[g]][[i]])) -
          (isigma[[g]]) %*% (mjPreviousValue[[g]][i, ] - mu[[g]])
        mj[[g]][i , ] <- mjPreviousValue[[g]][i , ] + Sj[[g]][[i]] %*% GXj[[g]][[i]] # update mj
      }

      miPreviousValue[[g]] <- mi[[g]]
      mjPreviousValue[[g]] <- mj[[g]]

      # Updating mu
      mu[[g]] <- (colSums(zValuei[, g] * mi[[g]]) + colSums(zValuej[, g] * mj[[g]])) /
                      (sum(zValuei[ , g]) + sum(zValuej[ , g]))

      # Updating sample covariance
      muMatrixi <- matrix(rep(mu[[g]], nObservationsi), nrow = nObservationsi, byrow = TRUE)
      muMatrixj <- matrix(rep(mu[[g]], nObservationsj), nrow = nObservationsj, byrow = TRUE)
      resi <- mi[[g]] - muMatrixi
      resj <- mj[[g]] - muMatrixj
      # Calculate weighted covariance
      tempi <- stats::cov.wt(resi, wt = zValuei[ , g], center = FALSE, method = "ML")
      tempj <- stats::cov.wt(resj, wt = zValuej[ , g], center = FALSE, method = "ML")
      Ski[ , , g] <- tempi$cov
      Skj[ , , g] <- tempj$cov
      sigma[[g]] <- ((tempi$cov + Reduce("+", zSiValue[[g]])) +
                        (tempj$cov + Reduce("+", zSjValue[[g]]))) /
                            (sum(zValuei[ , g]) + sum(zValuej[ , g]))
      isigma[[g]] <- solve(sigma[[g]])
    }

    piG <- (colSums(zValuei) + colSums(zValuej)) / nObservations


    # Matrix containing normaization factors
    normFactorsAsMatrixi <- matrix(rep(normFactors, nObservationsi),
                                  nrow = nObservationsi, byrow = TRUE)
    normFactorsAsMatrixj <- matrix(rep(normFactors, nObservationsj),
                                   nrow = nObservationsj, byrow = TRUE)

    # A useful function
    zvalueFuncTerm <- function(x, y = isigma[[g]]) {
      sum(diag(x %*% y))
    }

    # # Zvalue calculation
    forzi <- matrix(NA, ncol = G, nrow = nObservationsi)
    forzj <- matrix(NA, ncol = G, nrow = nObservationsj)

    for (g in 1:G) {
      # exp(m_igj + 0.5 S_ig,jj)
      twoi <- rowSums(exp(mi[[g]] + log(normFactorsAsMatrixi) +
                           0.5 * matrix(unlist(lapply(Si[[g]], diag)),
                                        ncol = dimensionality, byrow = TRUE)))
      twoj <- rowSums(exp(mj[[g]] + log(normFactorsAsMatrixj) +
                           0.5 * matrix(unlist(lapply(Sj[[g]], diag)),
                                        ncol = dimensionality, byrow = TRUE)))
      fivei <- 0.5 * unlist(lapply(Si[[g]], zvalueFuncTerm)) # trace(isigma x S_ig)
      fivej <- 0.5 * unlist(lapply(Sj[[g]], zvalueFuncTerm)) # trace(isigma x S_ig)
      sixi <- 0.5 * log(unlist(lapply(Si[[g]], det))) # 0.5 * log|S_ig|
      sixj <- 0.5 * log(unlist(lapply(Sj[[g]], det))) # 0.5 * log|S_ig|

      # For zig calculation (the numerator part)
      forzi[ , g] <- piG[g] *
        exp(rowSums(mi[[g]] * dataset[which(membership != 0), ]) -
              twoi - rowSums(lfactorial(dataset[which(membership != 0), ])) +
              rowSums(log(normFactorsAsMatrixi) * dataset[which(membership != 0), ]) -
              0.5 * mahalanobis(mi[[g]], center = mu[[g]], cov = sigma[[g]]) -
              fivei + sixi - 0.5 * log(det(sigma[[g]])) - (dimensionality / 2))


      # For zjg calculation (the numerator part)
      forzj[ , g] <- piG[g] *
        exp(rowSums(mj[[g]] * dataset[which(membership == 0), ]) -
              twoj - rowSums(lfactorial(dataset[which(membership == 0), ])) +
              rowSums(log(normFactorsAsMatrixj) * dataset[which(membership == 0), ]) -
              0.5 * mahalanobis(mj[[g]], center = mu[[g]], cov = sigma[[g]]) -
              fivej + sixj - 0.5 * log(det(sigma[[g]])) - (dimensionality / 2))
    }


    # Calculate zValue value
    # check which forz == 0 and rowSums(forz) == 0 and which of these
    # have both equalling to 0 (because 0/0 = NaN)
    if (G == 1) {
      errorpossiblei <- Reduce(intersect,
                              list(which(forzi == 0),
                                   which(rowSums(forzi) == 0)))
      forzi[errorpossiblei] <- 1e-100
      zvaluei <- forzi / rowSums(forzi)

      errorpossiblej <- Reduce(intersect,
                               list(which(forzj == 0),
                                    which(rowSums(forzj) == 0)))
      forzj[errorpossiblej] <- 1e-100
      zvaluej <- forzj / rowSums(forzj)
    } else {

      # check for error, if rowsums are zero
      rowSumsZeroi <- which(rowSums(forzi) == 0)
      if(length(rowSumsZeroi) > 0) {
        forzi[rowSumsZeroi, ] <- mclust::unmap(stats::kmeans(log(dataset + 1 / 6),
                                                            centers = G,
                                                            nstart = 100)$cluster)[rowSumsZeroi, ]
        zvaluei <- forzi / rowSums(forzi)
      } else {
        zvaluei <- forzi / rowSums(forzi)
      }

      rowSumsZeroj <- which(rowSums(forzj) == 0)
      if(length(rowSumsZeroj) > 0) {
        forzj[rowSumsZeroj, ] <- mclust::unmap(stats::kmeans(log(dataset + 1 / 6),
                                                             centers = G,
                                                             nstart = 100)$cluster)[rowSumsZeroj, ]
        zvaluej <- forzj / rowSums(forzj)
      } else {
        zvaluej <- forzj / rowSums(forzj)
      }


    }

    # if z generated has less clusters than numbG, then use random initialization
    # checkClusters <- 0
    # if(length(unique(mclust::map(zvalue))) < G) {
    #   while(! checkClusters) {
    #     zvalue <- randomInitfunction(gmodel = G, nObservations = nObservations)
    #     if(length(unique(mclust::map(zvalue))) == G) {
    #       checkClusters <- 1
    #       # cat("\n checkClusters", checkClusters)
    #     }
    #   }
    # }

    zValue <- matrix(0, nObservations, G)
    zValue[which(membership != 0),] <- zValuei
    zValue[which(membership == 0),] <- zValuej


    # if z generated has less clusters than numbG, then NA
    if(length(unique(mclust::map(zValue))) < G) {
      cat("\n length(unique(mclust::map(zvalue))) < G for G =", G)
      checks <- 2

    } else { # if z generated clusters == numbG

      # Calculate log-likelihood
      logLikelihood[itOuter] <- sum(log(rowSums(forzi) + rowSums(forzj)))

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
  }

  if(checks == 1) {
    # Naming parameters
    names(mu) <- names(sigma) <- paste0(rep("G=", G), 1:G)

    # Saving results for output
    Results <- list(piG = piG,
                    mu = mu,
                    sigma = sigma,
                    probaPost = zvalue,
                    clusterlabels = mclust::map(zvalue),
                    logLikelihood = logLikelihood)

    class(Results) <- "varMPLNClassification"
    return(Results)
  } else if(checks == 2) {
    # if z generated has less clusters than numbG, then NA
    return(NA)
  }

}





varMPLNClassificationInit <- function(dataset,
                                      numbG,
                                      initMethod,
                                      nInitIterations,
                                      membership,
                                      normFactors) {

  dimensionality <- ncol(dataset)
  nObservations <- nrow(dataset)

  zValue <- zValueKnown <- initRuns <- list()
  logLinit <- vector()

  # Initialize z with known memberships
  classIndicator  <- (membership == 0)
  zValueKnown[[1]] <- matrix(0, nObservations, G)
  for (i in 1:nObservations) {
    if(classIndicator[i]) {zValueKnown[[1]][i, ] <- 0} # zjg (for non-labeled)
    else{zValueKnown[[1]][i, membership[i]] <- 1} # zig (for labeled)
  }

  for(iterations in seq_along(1:nInitIterations)) {

    # the known memberships will always remain same
    zValueKnown[[i]] <- zValueKnown[[1]]

    # setting seed, to ensure if multiple iterations are selected by
    # user, then each run will give a different result.
    set.seed(iterations)
    if (initMethod == "kmeans" | is.na(initMethod)) {
      zValue[[iterations]] <- mclust::unmap(stats::kmeans(log(dataset + 1 / 6),
                                                          centers = numbG, nstart = 100)$cluster )
      # if z generated has less clusters than numbG, then use random initialization
      checkClusters <- 0
      if(length(unique(mclust::map(zValue[[iterations]]))) < numbG) {
        while(! checkClusters) {
          zValue[[iterations]] <- randomInitfunction(gmodel = numbG, nObservations = nObservations)
          if(length(unique(mclust::map(zValue[[iterations]]))) == numbG) {
            checkClusters <- 1
            # cat("\n checkClusters", checkClusters)
          }
        }
      }

      # assing estimated values only to the observations requiring classification
      zValueKnown[[i]][which(rowSums(zValueKnown) == 0), ] <-
        zValue[[iterations]][which(rowSums(zValueKnown) == 0), ]

    } else if (initMethod == "random") {
      zValue[[iterations]] <- randomInitfunction(gmodel = numbG, nObservations = nObservations)

      # assing estimated values only to the observations requiring classification
      zValueKnown[[i]][which(rowSums(zValueKnown) == 0), ] <-
        zValue[[iterations]][which(rowSums(zValueKnown) == 0), ]

    } else if (initMethod == "medoids") {
      zValue[[iterations]] <- mclust::unmap(cluster::pam(log(dataset + 1 / 3),
                                                         k = numbG,  cluster.only = TRUE))
      # if z generated has less clusters than numbG, then use random initialization
      checkClusters <- 0
      if(length(unique(mclust::map(zValue[[iterations]]))) < numbG) {
        while(! checkClusters) {
          zValue[[iterations]] <- randomInitfunction(gmodel = numbG, nObservations = nObservations)
          if(length(unique(mclust::map(zValue[[iterations]]))) == numbG) {
            checkClusters <- 1
            # cat("\n checkClusters", checkClusters)
          }
        }
      }

      # assing estimated values only to the observations requiring classification
      zValueKnown[[i]][which(rowSums(zValueKnown) == 0), ] <-
        zValue[[iterations]][which(rowSums(zValueKnown) == 0), ]

    } else if (initMethod == "clara") {
      zValue[[iterations]] <- mclust::unmap(cluster::clara(log(dataset + 1 / 3),
                                                           k = numbG)$cluster)
      # if z generated has less clusters than numbG, then use random initialization
      checkClusters <- 0
      if(length(unique(mclust::map(zValue[[iterations]]))) < numbG) {
        while(! checkClusters) {
          zValue[[iterations]] <- randomInitfunction(gmodel = numbG, nObservations = nObservations)
          if(length(unique(mclust::map(zValue[[iterations]]))) == numbG) {
            checkClusters <- 1
            # cat("\n checkClusters", checkClusters)
          }
        }
      }

      # assing estimated values only to the observations requiring classification
      zValueKnown[[i]][which(rowSums(zValueKnown) == 0), ] <-
        zValue[[iterations]][which(rowSums(zValueKnown) == 0), ]

    } else if (initMethod == "fanny") {
      zValue[[iterations]] <- mclust::unmap(cluster::fanny(log(dataset + 1 / 3),
                                                           k = numbG, memb.exp = numbG, cluster.only = TRUE)$clustering)
      # if z generated has less clusters than numbG, then use random initialization
      checkClusters <- 0
      if(length(unique(mclust::map(zValue[[iterations]]))) < numbG) {
        while(! checkClusters) {
          zValue[[iterations]] <- randomInitfunction(gmodel = numbG, nObservations = nObservations)
          if(length(unique(mclust::map(zValue[[iterations]]))) == numbG) {
            checkClusters <- 1
            # cat("\n checkClusters", checkClusters)
          }
        }
      }

      # assing estimated values only to the observations requiring classification
      zValueKnown[[i]][which(rowSums(zValueKnown) == 0), ] <-
        zValue[[iterations]][which(rowSums(zValueKnown) == 0), ]
    }


    initRuns[[iterations]] <- varMPLNInitClustering(dataset = dataset,
                                                    G = numbG,
                                                    zValue = zValueKnown[[iterations]],
                                                    normFactors = normFactors,
                                                    maxIterations = 10)

    # maxIterations set to 10 for initialization
    logLinit[iterations] <-
      unlist(utils::tail((initRuns[[iterations]]$logLikelihood), n = 1))
  }

  # select the initialization run with highest loglikelihood
  initializationResults <- initRuns[[which(logLinit == max(logLinit, na.rm = TRUE))[1]]]

  class(initializationResults) <- "varMPLNClassificationInit"
  return(initializationResults)
}




varMPLNInitClustering <- function(dataset,
                                  G,
                                  zValue,
                                  normFactors,
                                  membership,
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
    if(length(obs) > 1) {
      mu[[g]] <- colMeans(log(dataset[obs, ] + 1 / 6)) # starting value for mu
      # starting value for sample covariance matrix
      sigma[[g]] <- cov(log(dataset[obs, ] + 1 / 6))
      # starting value for inverse of sample covariance matrix
      # If the inverse is not present for covariance matrix, handle that
      isigma[[g]] <- tryCatch(solve(sigma[[g]]), error = function(err) NA)
      if(all(is.na(isigma[[g]]))) {
        isigma[[g]] <- diag(ncol(dataset[obs, ])) # if error with inverse
      }
    } else if(length(obs) == 1) {
      mu[[g]] <- log(dataset[obs, ] + 1 / 6) # starting value for mu
      # starting value for sample covariance matrix
      sigma[[g]] <- diag(ncol(dataset))
      # starting value for inverse of sample covariance matrix
      # If the inverse is not present for covariance matrix, handle that
      isigma[[g]] <- tryCatch(solve(sigma[[g]]), error = function(err) NA)
      if(all(is.na(isigma[[g]]))) {
        isigma[[g]] <- diag(ncol(dataset[obs, ])) # if error with inverse
      }
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
        # update S; will be used for updating sample covariance matrix
        S[[g]][[i]] <- tryCatch(solve(dGX[[g]][[i]]), error = function(err) NA)
        if(all(is.na(S[[g]][[i]]))) {
          S[[g]][[i]] <- diag(ncol(dGX[[g]][[i]])) # if error with inverse
        }


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
      if(length(rowSumsZero) > 0) {
        forz[rowSumsZero, ] <- mclust::unmap(stats::kmeans(log(dataset + 1 / 6),
                                                           centers = G,
                                                           nstart = 100)$cluster)[rowSumsZero, ]
        zvalue <- forz / rowSums(forz)
      } else {
        zvalue <- forz / rowSums(forz)
      }

      # if z generated has less clusters than numbG, then use random initialization
      checkClusters <- 0
      if(length(unique(mclust::map(zvalue))) < G) {
        while(! checkClusters) {
          zvalue <- randomInitfunction(gmodel = G, nObservations = nObservations)
          if(length(unique(mclust::map(zvalue))) == G) {
            checkClusters <- 1
            # cat("\n checkClusters", checkClusters)
          }
        }
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
                  zValue = zvalue,
                  m = m,
                  S = S,
                  clusterlabels = mclust::map(zvalue),
                  logLikelihood = logLikelihood)

  class(Results) <- "varMPLNInitClustering"
  return(Results)
}


