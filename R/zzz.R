#' @importFrom utils install.packages packageVersion

.onLoad <- function(libname, pkgname) {
  # Définir le solveur par défaut
  if (is.null(getOption("CVXR.solver"))) {
    options(CVXR.solver = "ECOS")
  }
}

#' Vérifier et configurer le solveur CVXR
#'
#' @param solver Nom du solveur ('OSQP', 'ECOS', 'SCS', 'CLARABEL')
#' @param install Si TRUE, tente d'installer le solveur manquant
#' @return Logical indiquant si le solveur est disponible
#' @export
setup_solver <- function(solver = "OSQP", install = FALSE) {
  solver_pkg <- switch(toupper(solver),
                       ECOS = "ECOSolveR",
                       SCS = "scs",
                       OSQP = "osqp",
                       CLARABEL = "CLARABEL")

  if (!requireNamespace(solver_pkg, quietly = TRUE)) {
    if (install) {
      message("Installing ", solver_pkg, "...")
      if (solver_pkg == "CLARABEL") {
        message("CLARABEL requires Rust compiler.")
        message("Try: install.packages('CLARABEL', type = 'binary')")
        message("Or: install.packages('osqp') for an alternative.")
      }
      install.packages(solver_pkg)
    }
    return(FALSE)
  }

  options(CVXR.solver = toupper(solver))
  message("CVXR solver set to: ", toupper(solver))
  return(TRUE)
}
