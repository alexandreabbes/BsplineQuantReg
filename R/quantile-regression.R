#' Unified Quantile Regression with B-Splines of Any Degree (1 to 4)
#'
#' This function provides a unified interface for quantile regression with
#' B-splines of degrees 1 to 4, with shape constraints (monotonicity,
#' convexity, third derivative).
#'
#' @param xtab Predictor vector (x)
#' @param ytab Response vector (y)
#' @param knot Knot vector or number of knots (quantiles are then used)
#' @param tau Quantile (between 0 and 1)
#' @param degree Spline degree: 1 (linear), 2 (quadratic), 3 (cubic), or 4 (quartic)
#' @param monot Monotonicity constraint vector per interval:
#'        1 = increasing, -1 = decreasing, 0 = unconstrained. If scalar, repeated.
#' @param convcons Convexity constraint vector:
#'        For degree 1: not available (ignored)
#'        For degree 2: per interval (1 = convex, -1 = concave)
#'        For degree 3: per knot (1 = convex, -1 = concave)
#'        For degree 4: per interval (1 = convex, -1 = concave)
#' @param der3cons Third derivative constraint:
#'        For degree 1: not available (third derivative = 0)
#'        For degree 2: not available (third derivative = 0)
#'        For degree 3: per interval (1 = positive, -1 = negative)
#'        For degree 4: per knot (1 = positive, -1 = negative)
#' @param solver 'CVXR' solver to use (default = 'CLARABEL')
#' @param callable render the final container object callable:
#'  'y=Bspline(x)' or 'y=Bspline(x,Bvalues)' for evaluation at x
#' @param weight Observation weights (default = 1 for all)
#' @param verbose logical; if TRUE, print progress messages
#' @param type_reg 'quantile' or 'mean_square' type of regression
#' @return A list containing coefficients, degree, and knots
#' @examples
#' # Generate data
#' x <- seq(0, 1, length.out = 100)
#' y <- 2*x + sin(6*pi*x)/2 + rnorm(100, 0, 0.05)
#' knot <- quantile(x, probs = seq(0, 1, length.out = 10))
#'
#' # Linear spline (degree 1) with increasing constraint
#' fit_lin <- quantile_spline(x, y, knot, tau = 0.5, degree = 1, monot = 1)
#'
#' # Quadratic spline (degree 2) with convexity constraint
#' fit_quad <- quantile_spline(x, y, knot, tau = 0.5, degree = 2, convcons = 1)
#'
#' # Cubic spline (degree 3) with both constraints
#' fit_cubic <- quantile_spline(x, y, knot, tau = 0.5, degree = 3,
#'                              monot = 1, convcons = 1)
#'
#' # Quartic spline (degree 4) with all constraints
#' fit_quart <- quantile_spline(x, y, knot, tau = 0.5, degree = 4,
#'                              monot = 1, convcons = 1, der3cons = 1)
#' @references
#' \itemize{
#'   \item Abbes, A. (2025). \emph{Quantile regression with cubic polynomial splines
#'         under shape constraints with applications}. Zenodo.
#'         \doi{10.5281/zenodo.16999784}
#'   \item Karlin, S., & Studden, W. J. (1966). \emph{Tchebycheff Systems: With
#'         Applications in Analysis and Statistics}. Interscience.
#' }
#' @seealso
#' \code{\link{SplineLinearQuant}}, \code{\link{SplineQuadraticQuant}},
#' \code{\link{SplineCubicQuant}}, \code{\link{SplineQuarticQuant}}
#' @export
quantile_spline <- function(xtab,
                            ytab,
                            knot=NULL,
                            tau = 0.5,
                            degree = 3,
                            monot = 0,
                            convcons = 0,
                            der3cons = 0,
                            solver = "CLARABEL",
                            weight = NULL,
                            verbose = FALSE,
                            callable = TRUE,
                            type_reg='quantile') {
  if (is.null(knot)){knot=c(min(xtab),max(xtab))}

    # Validate degree
  if (degree < 0 || degree > 4) {
    stop("degree must be between 1 and 4. Received: ", degree)
  }
  if (degree == 0) {
    fit <- SplineConstantQuant(
      xtab,
      ytab,
      knot,
      tau = tau,
      monot = monot,
      solver = solver,
      weight = weight,
      verbose = verbose,
      type_reg=type_reg
    )}

  # Dispatch to the appropriate function based on degree
  if (degree == 1) {
    fit <- SplineLinearQuant(
      xtab,
      ytab,
      knot,
      tau = tau,
      monot = monot,
      solver = solver,
      weight = weight,
      verbose = verbose,
      type_reg=type_reg
    )
  } else if (degree == 2) {
    fit <- SplineQuadraticQuant(
      xtab,
      ytab,
      knot,
      tau = tau,
      monot = monot,
      convcons = convcons,
      solver = solver,
      weight = weight,
      verbose = verbose,
      type_reg=type_reg
    )
  } else if (degree == 3) {
    fit <- SplineCubicQuant(
      xtab,
      ytab,
      knot,
      tau = tau,
      monot = monot,
      convcons = convcons,
      der3cons = der3cons,
      solver = solver,
      weight = weight,
      verbose = verbose,
      type_reg=type_reg
    )
  } else if (degree == 4) {
    fit <- SplineQuarticQuant(
      xtab,
      ytab,
      knot,
      tau = tau,
      monot = monot,
      convcons = convcons,
      der3cons = der3cons,
      solver = solver,
      weight = weight,
      verbose = verbose,
      type_reg=type_reg
    )
  }

  # Return callable object if requested

  result = make_spline(fit, callable = callable)

  return(result)
}
