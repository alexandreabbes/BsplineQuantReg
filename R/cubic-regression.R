#' Constrained quantile regression with cubic splines
#'
#' Performs quantile regression using cubic B-splines, with optional
#' monotonicity constraints (via Karlin-Studden) and convexity constraints.
#'
#' @param xtab Predictor vector (x)
#' @param ytab Response vector (y)
#' @param knot Knot vector or number of knot (quantiles are then used)
#' @param tau Quantile (between 0 and 1)
#' @param monot Monotonicity constraint vector per interval:
#'        1 = increasing, -1 = decreasing, 0 = unconstrained. If scalar, repeated.
#' @param convcons Convexity constraint vector per knot:
#'        1 = convex, -1 = concave, 0 = unconstrained. If scalar, repeated.
#' @param der3cons Constraint on the 3rd derivative (on esach intervall:
#'        -1: négative, 0: no constraint, 1: positive constraint
#' @param solver CVXR solver to use (default = "CLARABEL")
#' @param weight Observation weights (default = 1 for all)
#' @param verbose boolean FALSE (default) or TRUE.
#' @return A list containing:
#'   \item{coefficients}{B-spline coefficients (including y mean)}
#'   \item{degree}{Spline degree (always 3)}
#'   \item{ext_knot}{ext_knot vector used}
#'   \item{int_knot}{knot vector}
#' @examples
#' #optional set.seed(42)
#' x <- seq(0, 1, length=100)
#' y <- 2*x + sin(6*pi*x)/2 + rnorm(100, 0, 0.05)
#' knot <- quantile(x, probs=seq(0,1,length.out=10))
#'
#' # Median quantile regression without constraints
#' fit <- SplineConstQuantRegBs3(x, y, knot, tau=0.5)
#'
#' # With increasing monotonicity constraint
#' fit_monot <- SplineConstQuantRegBs3(x, y, knot, tau=0.5, monot=1)
#'
#' # With convexity constraint
#' fit_convex <- SplineConstQuantRegBs3(x, y, knot, tau=0.5, convcons=1)
#'
#' @seealso
#' Related R packages:
#' \itemize{
#'   \item \code{quantreg} - Quantile regression with linear programming
#'   \item \code{cobs} - Constrained B-sines (linear and quadratic only)
#' }
#'
#' Other implementations:
#' \itemize{
#'   \item MATLAB/Python versions: \url{https://github.com/alexandreabbes/Constrained-Quantile-Regression-with-cubic-splines}
#' }
#' @references
#' \itemize{
#'   \item Abbes, A. (2025). \emph{Quantile regression with cubic polynomial splines under shape constraints with applications}
#'         . Zenodo.
#'         \doi{10.5281/zenodo.16999784}
#'   \item de Boor, C. (1978). \emph{A Practical Guide to Splines}. Springer-Verlag.
#'         \doi{10.1007/978-1-4612-6333-3}
#'   \item Karlin, S., & Studden, W. J. (1966). \emph{Tchebycheff Systems: With
#'         Applications in Analysis and Statistics}. Interscience.
#'   \item Koenker, R., & Bassett, G. (1978). Regression Quantiles.
#'         \emph{Econometrica}, 46(1), 33-50. \doi{10.2307/1913643}
#'   \item Koenker, R. (2025). quantreg: Quantile Regression. R package version 5.99.
#'         \url{https://CRAN.R-project.org/package=quantreg}
#'   \item Ng, P., & Maechler, M. (2024). cobs: Constrained B-Splines.
#'         R package version 1.3-8. \url{https://CRAN.R-project.org/package=cobs}
#' }
#' @export


SplineCubicQuant<- function(xtab, ytab, knot, tau,
                                   monot = 0,
                                   convcons=0,
                                   der3cons=0,
                                   solver = "CLARABEL", weight = NULL,
                                   verbose=FALSE)
{
  if (is.null(weight)) {
    weight <- rep(1, length(xtab))
  }

  ordre <- order(xtab)
  xtab <- xtab[ordre]
  ytab <- ytab[ordre]
  weight <- weight[ordre]

  n <- length(xtab)

  if (length(knot) == 1 && is.numeric(knot))
  {
    kn <- knot - 1
    knot <- quantile(xtab, probs = (0:kn)/(kn))
  }

  kn <- length(knot) - 1
  if (verbose) {
    message("knot:", knot, "\n")
  }

  if (verbose) {
  message("Monotonicity constraints (Karlin):", monot, "\n")
  }
  boundary_knot <- range(knot)
  degree=3
  N=length(knot)+3-1
  # Calcul des coefficients normalises des derivees
  deriv_spline <- bspline_to_deriv_coeffs_cubic(knot, degree = 3,x_values=xtab,verbose=verbose)
  deriv_coeffs <-deriv_spline$d1 # kn x N x 3
  deriv_coeffs2<-deriv_spline$d2 # (kn+1) x N
  deriv_coeffs3<-deriv_spline$d3 # kn x N
  B<-deriv_spline$d0
  B=t(B)
  y_mean <- mean(ytab)
  ytab_centered <- ytab - y_mean

  alpha <- Variable(N)
  #Contraintes monotones
  # Variables auxiliaires z
  residuals <- ytab_centered - (B %*% alpha)

  u_plus <- pos(residuals)
  u_minus <- pos(-residuals)
  weighted_loss <- sum(weight * (tau * u_plus + (1 - tau) * u_minus))

  objective <- Minimize(weighted_loss)

  constraints <- list()
  z_vars <- list()

  if (any(monot != 0)) {
    if (length(monot) == 1) {
      monot <- rep(monot, kn)
    }
    if (length(monot)<(kn)){
      message("Not enough monotonicity constraints, completing with 0")
      monot<-c(monot,rep(0,kn-length(convcons))) }

    for (i in 1:(kn)) {
      if (monot[i] != 0) {
        z_vars[[i]] <- Variable(1, name = paste0("z", i))
        a_coef=sum(deriv_coeffs[i,,1]*alpha) *monot[i]
        b_coef=sum(deriv_coeffs[i,,2]*alpha) *monot[i]
        c_coef=sum(deriv_coeffs[i,,3]*alpha) *monot[i]
        #a*x^2+b*x+c
        CK<-apply_karlin_quadratic(a_coef,b_coef,c_coef,z_vars[[i]])
        constraints<-c(constraints,CK)
      }
    }
  }

  #"contraintes convexes

  # eliminate the null (unconstrained) case
  if (any(convcons !=0)){
    if (length(convcons) == 1) {
      convcons <- rep(convcons, (kn+1) )}
    if (length(convcons)<(kn+1)){
      message("Not enough convexity constraints, completing with 0")
      convcons<-c(convcons,rep(0,kn+1-length(convcons))) }
    if (verbose){
      message("Convexity constraints (Linear):", convcons, "\n")
      }
    CV<-list((convcons*(deriv_coeffs2 %*% alpha))>=0)
    # Very simple, only use the sign of
    # the second derivatives at the knot.
    constraints<-c(constraints,CV)
  }

#3rd derivative constraints
  if (any(der3cons!=0)){
    #if (dim(deriv_coeffs3)[2]!=N){
      #prepare the values for matrix mult.      #calculation
     # deriv_coeff3=t(deriv_coeffs3)}
    if (length(der3cons) == 1){der3cons<-rep(der3cons,kn)}
    if (length(der3cons) < kn ){
      message("not enough 3rd order constraints, completing with 0")
      der3cons=c(der3cons,rep(0,kn-length(der3cons))) }
    if (verbose) {
      message("3rd order derivative constraints (constant):", der3cons, "\n")}
    for (j in 1:kn) {
      sig <- der3cons[j]
      if (sig != 0) {
        # Utiliser sum() au lieu de vdot
        d3 <- sum(deriv_coeffs3[j, ] * alpha)
        DER3 <- apply_linear_constraint(d3, sig)
        constraints <- c(constraints, DER3)
        if (verbose) {
          message("  Constraint on interval ", j, ": ", sig, " * d3 >= 0")
        }
      }
    }
    }

  # Solve the problem
  problem <- Problem(objective, constraints)

  result <- NULL
  fallback_result <- NULL
  solvers_to_try <- c(solver, "CLARABEL", "OSQP", "ECOS", "SCS")
  solvers_to_try <- unique(solvers_to_try)  # Supprimer les doublons

  for (s in solvers_to_try) {
    if (verbose) cat("Trying solver:", s, "\n")

    # Vérifier la disponibilité de ECOS
    if (s == "ECOS" && !requireNamespace("ECOSolveR", quietly = TRUE)) {
      if (verbose) cat("  ECOS not available (ECOSolveR missing)\n")
      next
    }

    result <- tryCatch({
      opt_val <- psolve(problem, solver = toupper(s), verbose = verbose)
      list(
        value = opt_val,
        status = status(problem),
        alpha_value = value(alpha)
      )
    }, error = function(e) {
      if (verbose) cat("  Failed:", e$message, "\n")
      NULL
    })

    # Vérifier si le solveur a réussi
    if (!is.null(result) && !is.null(result$alpha_value)) {
      if (result$status == "optimal") {
        if (verbose) cat("  Solver succeeded with optimal status:", s, "\n")
        break  # OK, on sort de la boucle
      } else if (result$status == "optimal_inaccurate") {
        if (verbose) cat("  Solver returned optimal_inaccurate:", s, "\n")
        fallback_result <- result
        # Continuer à essayer d'autres solveurs pour un meilleur résultat
      } else {
        if (verbose) cat("  Solver returned non-optimal status:", result$status, "\n")
        fallback_result <- result
      }
    }
  }

  # Après la boucle, vérifier le résultat
  if (is.null(result) || is.null(result$alpha_value)) {
    # Utiliser le fallback si disponible
    if (!is.null(fallback_result)) {
      result <- fallback_result
      if (verbose) cat("Using fallback result with status:", result$status, "\n")
    } else {
      warning("Optimisation did not converge with any available solver")
      return(NULL)
    }
  }

  # Si le résultat est optimal_inaccurate, on peut quand même l'utiliser avec un avertissement
  if (result$status == "optimal_inaccurate") {
    warning("Solution may be inaccurate. Try another solver or adjust settings.")
  }

  alpha_val <- result$alpha_value + y_mean
  result$y_mean <- y_mean

  if (verbose) {
    message(" Statut:", result$status, "\n",
            "Valeur objectif:", result$value, "\n",
            "Coefficients alpha (range):", range(alpha_val), "\n")
  }


  return(list(
    coefficients = alpha_val,
    degree=3,
    knot=knot,
    result=result
  ))
}

#' Wrapper function for Cubic spline quantile regression
#'
#' This is a convenience wrapper for SplineCubicQuant
#' with convention of previous versions  SplineConstQuantRegBs3.
#' @inheritParams SplineCubicQuant
#' @return Same as SplineCubicQuant
#' @export

SplineConstQuantRegBs3<-function(xtab, ytab, knot, tau,
                                         monot = 0,
                                         convcons=0,
                                         der3cons=0,
                                         solver = "CLARABEL", weight = NULL,
                                         verbose=FALSE)
          {SplineCubicQuant(xtab, ytab, knot, tau,
          monot = monot,
          convcons=convcons,
          der3cons=der3cons,
          solver = solver, weight = weight,
          verbose=verbose)
          }

