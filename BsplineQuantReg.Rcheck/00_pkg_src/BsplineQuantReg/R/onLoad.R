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
    "This is a BETA version. The API may change in future releases.\n",
    "Use demo() to see examples: demo(package = 'BsplineQuantReg')"
  )
}
