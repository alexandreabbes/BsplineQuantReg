#' Package load message
#'
#' @keywords internal
.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "ConstrainedQuantileSplines loaded.\n",
    "Polynomial and spline functions are reimplemented for consistency.\n",
    "Use test_karlin_simple() for a demo."
  )
}
