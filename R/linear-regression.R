# Linear Spline Quantile Regression with Shape Constraints
# Author: Alexandre Abbes
#' Quantile regression with linear splines and monotonicity constraints
#'
#' Performs quantile regression using linear B-splines with monotonicity constraints.
#' For linear splines, the derivative is constant on each interval, so monotonicity
#' is a simple linear constraint on the derivative.
#'
#' @param xtab Predictor vector (x)
#' @param ytab Response vector (y)
#' @param knot Knot vector or number of knots
#' @param tau Quantile (between 0 and 1)
#' @param monot Monotonicity constraint vector per interval:
#'        1 = increasing, -1 = decreasing, 0 = unconstrained
#' @param solver CVXR solver to use (default = "OSQP")
#' @param weight Observation weights (default = 1 for all)
#' @param verbose logical; if TRUE, print progress messages
#' @param type_reg 'quantile' or 'mean_square' type of regression
#' i.e form of the objective.
#' @return A list containing coefficients, degree, and knots
#' @export
SplineLinearQuant <- function(xtab,
                              ytab,
                              knot=NULL,
                              tau=0.5,
                              monot = 0,
                              solver = "CLARABEL",
                              weight = NULL,
                              verbose = FALSE,
                              type_reg='quantile') {
  if (is.null(knot)){knot=c(min(xtab),max(xtab))}
  if (is.null(weight)) {
    weight <- rep(1, length(xtab))
  }

  # Sort data
  ordre <- order(xtab)
  xtab <- xtab[ordre]
  ytab <- ytab[ordre]
  weight <- weight[ordre]

  n <- length(xtab)

  # Handle knots
  if (length(knot) == 1 && is.numeric(knot)) {
    kn <- knot - 1
    knot <- as.numeric(quantile(xtab, probs = seq(0, 1, length.out = kn + 1)))
  }

  kn <- length(knot) - 1
  degree <- 1
  N <- kn + degree  # Number of basis functions for linear = kn + 1

  if (verbose) {
    message("=== Linear Spline Quantile Regression (degree = 1) ===")
    message(sprintf("Knots: %d, Basis functions: %d", kn, N))
  }

  # Handle monotonicity constraints
  if (length(monot) == 1) {
    monot <- rep(monot, kn)
  }

  if (verbose) {
    message("Monotonicity constraints:", paste(monot, collapse = " "))
  }

  # Build B-spline basis and derivative coefficients
  deriv_data <- bspline_to_deriv_coeffs_lin(knot,
                                            degree = 1,
                                            x_values = xtab,
                                            verbose = verbose)

  B <- deriv_data$d0
  B <- t(B)  # Design matrix: n x N

  # Derivative coefficients (constant per interval)
  deriv1_coeffs <- deriv_data$d1  # [kn, N]

  # Center data
  y_mean <- mean(ytab)
  ytab_centered <- ytab - y_mean

  # Optimization variables
  alpha <- Variable(N)

  # Objective function
  residuals <- ytab_centered - B %*% alpha
  if (type_reg=='mean_square'){
    weighted_loss <- norm2(residuals)
  }
   else{
   u_plus <- pos(residuals)
   u_minus <- pos(-residuals)
   weighted_loss <- sum(weight * (u_plus*tau+u_minus*(1-tau)))
   }
  objective <- Minimize(weighted_loss)

  constraints <- list()

  # Monotonicity constraints (linear on derivative)
  if (any(monot != 0)) {
    for (i in 1:kn) {
      if (monot[i] != 0) {
        # Derivative on interval i is constant
        deriv_val <- sum(alpha * deriv1_coeffs[i, ])

        if (monot[i] > 0) {
          constraints <- c(constraints, list(deriv_val >= 0))
        } else {
          constraints <- c(constraints, list(deriv_val <= 0))
        }
      }
    }
  }

  # Solve the problem
  problem <- Problem(objective, constraints)

  result <- NULL
  solvers_to_try <- c(solver, "HIGHS", "CLARABEL", "GUROBI", "OSQP", "ECOS", "SCS")

  for (s in unique(solvers_to_try)) {
    if (verbose)
      cat("Trying solver:", s, "\n")

    # Use new CVXR syntax: psolve() for optimal value
    result <- tryCatch({
      # Solve the problem with new syntax
      opt_val <- psolve(problem, solver = toupper(s), verbose = verbose)

      # Create a result list compatible with old expectations
      list(
        value = opt_val,
        status = status(problem),
        alpha_value = value(alpha)
      )
    }, error = function(e) {
      if (verbose)
        cat("Failed:", e$message, "\n")
      NULL
    })
  }

  if (is.null(result) || is.null(value(alpha))) {
    warning("Optimization did not converge with any solver")
    return(NULL)
  }

  alpha_val <- result$alpha_value + y_mean
  result$y_mean <- y_mean

  if (verbose) {
    message("Status:", result$status)
    message("Objective value:", result$value)
  }

  # Return results
  return(list(
    coeff = alpha_val,
    degree = degree,
    knot = knot,
    result = result
  ))
}

#' Wrapper function for linear spline quantile regression
#'
#' This is a convenience wrapper for SplineLinearQuant.
#'
#' @inheritParams SplineLinearQuant
#' @return Same as SplineLinearQuant
#' @export
SplineConstQuantRegBs1 <- function(xtab,
                                   ytab,
                                   knot,
                                   tau,
                                   monot = 0,
                                   solver = "CLARABEL",
                                   weight = NULL,
                                   verbose = FALSE) {
  SplineLinearQuant(
    xtab,
    ytab,
    knot,
    tau,
    monot = monot,
    solver = solver,
    weight = weight,
    verbose = verbose
  )
}
