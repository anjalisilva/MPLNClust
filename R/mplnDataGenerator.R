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
#' sampleData <- mplnDataGenerator(nObservations = 100,
#'                                 dimensionality = 6,
#'                                 mixingProportions = c(0.79, 0.21),
#'                                 mu = rbind(trueMu1, trueMu2),
#'                                 sigma = rbind(trueSigma1, trueSigma2),
#'                                 produceImage = "Yes",
#'                                 ImageName = "TwoComponents")
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
  if(class(nObservations) != "numeric") {
    stop("nObservations argument should be of class numeric.")
  }

  if(class(dimensionality) != "numeric") {
    stop("dimensionality argument should be of class numeric.")
  }

  if(class(mixingProportions) != "numeric") {
    stop("mixingProportions argument should be a vector of class numeric.")
  }

  if (sum(mixingProportions) != 1) {
    stop("mixingProportions argument should be a vector that sum to 1.")
  }

  if(class(mu) != "matrix") {
    stop("mu argument should be of class matrix.")
  }

  if(ncol(mu) != dimensionality) {
    stop("mu should be a matrix, which has number of columns equalling dimensionality.")
  }

  if(nrow(mu) != length(mixingProportions)) {
    stop("mu should be a matrix, which has number of rows equalling number of components.")
  }

  if(class(sigma) != "matrix") {
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

  if (produceImage == "Yes" && class(ImageName) != "character") {
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



