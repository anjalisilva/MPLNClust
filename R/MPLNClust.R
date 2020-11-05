#' MPLNClust: Clustering via mixtures of multivariate Poisson-log normal distribution
#'
#' `MPLNClust` is an R package for performing clustering using mixtures of
#' multivariate Poisson-log normal (MPLN) distribution proposed by
#' Silva et al., 2019. It was developed for count data, with clustering of
#' RNA sequencing data as a motivation. However, the vector of normalization
#' factors can be relaxed and clustering method may be applied to other
#' types of count data.
#'
#' @section MPLNClust functions:
#' The MPLNClust package provides 10 functions:
#' \itemize{
#'   \item mplnVariational
#'   \item runMPLNClust
#'   \item mplnMCMCParallel
#'   \item mplnMCMCNonParallel
#'   \item mplnVisualize
#'   \item mplnDataGenerator
#'   \item AICFunction
#'   \item BICFunction
#'   \item AIC3Function
#'   \item ICLFunction
#' }
#' For a quick introduction to MPLNClust see the vignettes.
#'
#' @author {Anjali Silva, \email{anjali.silva@uhnresearch.ca}, Sanjeena Dang,
#'          \email{sdang@math.binghamton.edu}. }
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
#' sampleData <- MPLNClust::mplnDataGenerator(nObservations = 1000,
#'                                             dimensionality = 6,
#'                                             mixingProportions = c(0.79, 0.21),
#'                                             mu = rbind(trueMu1, trueMu2),
#'                                             sigma = rbind(trueSigma1, trueSigma2),
#'                                             produceImage = "No")
#'
#' # Clustering via mplnVariational
#' mplnResults <- MPLNClust::mplnVariational(dataset = sampleData$dataset,
#'                                           membership = sampleData$trueMembership,
#'                                           gmin = 1,
#'                                           gmax = 2,
#'                                           initMethod = "kmeans",
#'                                           nInitIterations = 2,
#'                                           normalize = "Yes")
#' \dontrun{
#' # Clustering via mplnMCMCParallel
#' mplnResults <- MPLNClust::mplnMCMCParallel(dataset = sampleData$dataset,
#'                                              membership = sampleData$trueMembership,
#'                                              gmin = 1,
#'                                              gmax = 1,
#'                                              nChains = 3,
#'                                              nIterations = 400,
#'                                              initMethod = "kmeans",
#'                                              nInitIterations = 0,
#'                                              normalize = "Yes",
#'                                              numNodes = 2)
#'
#' # Clustering via mplnMCMCNonParallel
#' mplnResults <- MPLNClust::mplnMCMCNonParallel(dataset = sampleData$dataset,
#'                                                membership = sampleData$trueMembership,
#'                                                gmin = 1,
#'                                                gmax = 1,
#'                                                nChains = 3,
#'                                                nIterations = 700,
#'                                                initMethod = "kmeans",
#'                                                nInitIterations = 0,
#'                                                normalize = "Yes")
#' }
#' @references
#' Silva, A. et al. (2019). A multivariate Poisson-log normal mixture model
#' for clustering transcriptome sequencing data. \emph{BMC Bioinformatics} 20.
#' \href{https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-2916-0}{Link}
#'
#' Subedi, S., and R. Browne (2020). A parsimonious family of multivariate Poisson-lognormal
#' distributions for clustering multivariate count data. arXiv preprint arXiv:2004.06857.
#' \href{https://arxiv.org/pdf/2004.06857.pdf}{Link}
#'
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
#' @docType package
#' @name MPLNClust
NULL
#> NULL
