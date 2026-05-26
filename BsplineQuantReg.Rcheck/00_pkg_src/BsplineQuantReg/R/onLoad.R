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
    "BsplineQuantReg(version 0.1.0-beta)\n",
    "Quantile regression with splines under shape constraints\n",
    "With independent B-splines tools"
#    "Use demo() to see examples: demo(package = 'BsplineQuantReg')"
  )
}
