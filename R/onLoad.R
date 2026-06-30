#' Package load message
#'
#' @keywords internal
.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "BsplineQuantReg loaded.\n",
    "Polynomial and spline functions are reimplemented for consistency.\n",
    "Use test_karlin_simple() for a demo."
  )
}
#' Package load message
#'
#' @keywords internal
.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "BsplineQuantReg (version ", utils::packageVersion("BsplineQuantReg"), ")\n",
    "Constrained Quantile Regression with Cubic B-Splines\n",
    "Available on CRAN: https://cran.r-project.org/package=BsplineQuantReg\n",
    "Use demo() to see examples: demo(package = 'BsplineQuantReg')"
  )
}
