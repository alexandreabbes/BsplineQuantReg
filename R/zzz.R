#' Get package version
#'
#' @return Current version string
#' @export
package_version <- function() {
  as.character(utils::packageVersion("BsplineQuantReg"))
}

#' Check if package is beta version
#'
#' @return TRUE if beta version
#' @export
is_beta <- function() {
  grepl("beta", as.character(utils::packageVersion("BsplineQuantReg")))
}
