#' Launch Shiny App For Package MPLNClust
#'
#' A function that launches the shiny app for this package.
#' The shiny app permit to perform clustering using mixtures
#' of MPLN via variational-EM. Model selection is performed
#' using AIC, AIC3, BIC and ICL.
#'
#' @return No return value but open up a shiny page.
#'
#' @examples
#' \dontrun{
#' runMPLNClust()
#' }
#'
#' @author Anjali Silva, \email{anjali.silva@uhnresearch.ca}
#'
#' @references
#' Silva, A. et al. (2019). A multivariate Poisson-log normal mixture model
#' for clustering transcriptome sequencing data. \emph{BMC Bioinformatics} 20.
#' \href{https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-2916-0}{Link}
#'
#' Subedi, S., and R. Browne (2020). A parsimonious family of multivariate Poisson-lognormal
#' distributions for clustering multivariate count data. arXiv preprint arXiv:2004.06857.
#' \href{https://arxiv.org/pdf/2004.06857.pdf}{Link}
#'
#' @export
#' @importFrom shiny runApp
runMPLNClust <- function() {
  appDir <- system.file("shiny-scripts",
                        package = "MPLNClust")
  shiny::runApp(appDir, display.mode = "normal")
  return()
}

