#' Launch the shiny app for package MPLNClust
#'
#' A function that launches the shiny app for this package.
#' The shiny app permit to perform clustering using mixtures
#' of MPLN via variational-EM. The code has been placed in
#' \code{./inst/shiny-scripts}.
#'
#' @return No return value but open up a shiny page.
#'
#' @examples
#' \dontrun{
#' runMPLNClust()
#' }
#'
#' @export
#' @importFrom shiny runApp

runMPLNClust <- function() {
  appDir <- system.file("shiny-scripts",
                        package = "MPLNClust")
  shiny::runApp(appDir, display.mode = "normal")
  return()
}

