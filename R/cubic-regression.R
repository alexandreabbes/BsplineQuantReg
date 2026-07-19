# SplineConstQuantRegBs3, apply_karlin_constraints

#' Local Wrapper function to apply_karlin_quadratic
apply_karlin_constraints<-function(p2, p1, p0, z0,verbose=FALSE)
  {
  apply_karlin_quadratic(p2=p2, p1=p1, p0=p0, z0=z0,sign=1,verbose=verbose)
  }

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
    knot <- quantile(xtab, probs = seq(0, 1, length.out = kn + 1))
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

  int_knot=knot[2:kn]

  # Calcul des coefficients normalises des derivees
  deriv_spline <- bspline_to_deriv_coeffs_pp(knot, degree = 3,x_values=xtab,verbose=verbose)
  deriv_coeffs <-deriv_spline$d1
  deriv_coeffs2<-deriv_spline$d2
  deriv_coeffs3<-deriv_spline$d3
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
    for (i in 1:(kn)) {
      if (monot[i] != 0) {
        z_vars[[i]] <- Variable(1, name = paste0("z", i))
        a_coef=sum(deriv_coeffs[i,,1]*alpha) *monot[i]
        b_coef=sum(deriv_coeffs[i,,2]*alpha) *monot[i]
        c_coef=sum(deriv_coeffs[i,,3]*alpha) *monot[i]
        #a*x^2+b*x+c
        CK<-apply_karlin_constraints(a_coef,b_coef,c_coef,z_vars[[i]],verbose=verbose)
        constraints<-c(constraints,CK)
      }
    }
  }

  #"contraintes convexes
  # eliminate the null (unconstrained) case
  if (any(convcons !=0)){
    if (length(convcons) == 1) {
      convcons <- rep(convcons, (kn+1))}
    CV<-list((convcons*(deriv_coeffs2 %*% alpha))>=0)
    # Very simple, only use the sign of
    # the second derivatives at the knot.
    constraints<-c(constraints,CV)
  }
#3rd derivative constraints
  if (any(der3cons!=0)){
    if (length(der3cons)==1){der3cons<-rep(der3cons,kn)}

    DER3<-list((der3cons*(deriv_coeffs3%*%alpha))>=0)
    constraints<-c(constraints,DER3)
  }
  problem <- Problem(objective, constraints)

  result <- NULL
  solvers_to_try <- c(solver, "CLARABEL", "OSQP", "ECOS", "SCS")

  for (s in unique(solvers_to_try)) {
    if (verbose) {
      message("attempt with  solver:", s, "\n")
      }
    result <- tryCatch(
      psolve(problem, solver = toupper(s), verbose=verbose),
    error = function(e) {warning("Missed:", e$message, "\n")
      NULL} )
    if (!is.null(result) && !is.null(value(alpha))) {
      if (verbose) {
        message("Solveur succeeded:", s, "\n")
      }
      break}
  }

  if (is.null(result) || is.null(value(alpha))) {
    warning("Optimisation did not converge with any available solver")
    return(NULL)
  }

  alpha_val <- value(alpha)+y_mean
  if (verbose) {
    message(" Statut:", result$status, "\n",
            "Valeur objectif:", result$value, "\n",
             "Coefficients alpha (range):", range(alpha_val), "\n")
}

  return(list(
    #spline = spline_result,
    coefficients = alpha_val,
    degree=3,
    #basis_matrix = B,
    knot = knot,
    int_knot = knot
  ))
}

#' Wrapper function for Cubic spline quantile regression
#'
#' This is a convenience wrapper for SplineCubicQuant
#' with convention of previous versions  SplineConstQuantRegBs3.
#'
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
          verbose=verbose)}

