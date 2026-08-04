# Quadratic Spline Quantile Regression with Shape Constraints
# Based on the Python implementation in quantile_reg.py
# Author: Alexandre Abbes

#' Quantile regression with quadratic splines and shape constraints
#'
#' Performs quantile regression using quadratic B-splines with:
#' - Monotonicity: Karlin-Studden constraints on the derivative (linear)
#' - Convexity: Karlin-Studden constraints on the second derivative (constant)
#'
#' @param xtab Predictor vector (x)
#' @param ytab Response vector (y)
#' @param knot Knot vector or number of knots
#' @param tau Quantile (between 0 and 1)
#' @param monot Monotonicity constraint vector per interval:
#'        1 = increasing, -1 = decreasing, 0 = unconstrained.
#'        Although constraints are set at knots; we apply them
#'        on each interval, at both extremities.
#'        kn+1 knots, kn constraints taken into account.
#' @param convcons Convexity constraint vector per interval:
#'        1 = convex, -1 = concave, 0 = unconstrained.
#'        kn constraints are considered for kn+1 knots.
#' @param solver 'CVXR' solver to use (default = 'CLARABEL')
#' @param weight Observation weights (default = 1 for all)
#' @param verbose logical; if TRUE, print progress messages
#' @param type_reg 'quantile' or 'mean_square' type of regression,
#' @return A list containing coefficients, degree, and knots
#' @export
SplineQuadraticQuant <- function(xtab,
                                 ytab,
                                 knot=NULL,
                                 tau=0.5,
                                 monot = 0,
                                 convcons = 0,
                                 solver = "CLARABEL",
                                 weight = NULL,
                                 verbose = FALSE,
                                 type_reg='quantile'){
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
    knot <- quantile(xtab, probs = seq(0, 1, length.out = kn + 1))
  }

  kn <- length(knot) - 1
  degree <- 2
  N <- kn + degree  # Number of basis functions for quadratic = kn + 2

  if (verbose) {
    message("=== Quadratic Spline Quantile Regression (degree = 2) ===")
    message(sprintf("Knots: %d, Basis functions: %d", kn, N))
  }

  # Handle constraints
  if (any(monot != 0)) {
    if (length(monot) == 1) {
      monot <- rep(monot, kn)
    }
    if (length(monot) < (kn+1)) {
      message("Not enough monotonicity constraints, completing with 0")
      monot <- c(monot, rep(0, kn+1 - length(convcons)))
    }
  }
  # handle convexity constraints
  if (any(convcons != 0)) {
    if (length(convcons) == 1) {
      convcons <- rep(convcons, (kn))
    }
    if (length(convcons) < (kn + 1)) {
      message("Not enough convexity constraints, completing with 0")
      convcons <- c(convcons, rep(0, kn - length(convcons)))
    }
    if (verbose) {
      message("Convexity constraints (Linear):", convcons, "\n")
    }
  }


  if (verbose) {
    message("Monotonicity constraints:", paste(monot, collapse = ","))
    message("Convexity constraints:", paste(convcons, collapse = ","))
  }

  # Build B-spline basis and derivative coefficients
  deriv_data <- bspline_to_deriv_coeffs_quad(knot,
                                             degree = 2,
                                             x_values = xtab,
                                             verbose = verbose)

  B <- deriv_data$d0
  B <- t(B)  # Design matrix: n x N

  # Derivative coefficients
  deriv1_coeffs <- deriv_data$d1  # [kn, N, 2] for linear derivative
  deriv2_coeffs <- deriv_data$d2  # [kn, N] for constant second derivative

  # Center data
  y_mean <- mean(ytab)
  ytab_centered <- ytab - y_mean

  # Optimization variables
  alpha <- Variable(N)
  z_vars <- list()

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

  # 1. Monotonicity constraints
  # Derivative is linear: P'(u) = a*u + b
  if (any(monot != 0)) {
    for (i in 1:kn) {
      if (monot[i] != 0) {
        # Coefficients of the derivative on interval i
        # deriv1_coeffs[i, j, ] = [a, b] for a*u + b
        a_coef <- sum(alpha * deriv1_coeffs[i, , 1])
        b_coef <- sum(alpha * deriv1_coeffs[i, , 2])

        # Apply sign of monotonicity
        s <- monot[i]
        if (verbose) {
          message("applying monotonicity ",
                  s,
                  " on intervall No",
                  i,
                  "[",
                  knot[i],
                  ",",
                  knot[i + 1],
                  "]")
        }
        if (s > 0) {
          # P'(u) >= 0 on [0,1] for linear: min(P'(0), P'(1)) >= 0
          constraints <- c(constraints, list(b_coef >= 0))
          constraints <- c(constraints, list(a_coef + b_coef >= 0))
        } else {
          # P'(u) <= 0 on [0,1]: max(P'(0), P'(1)) <= 0
          constraints <- c(constraints, list(b_coef <= 0))
          constraints <- c(constraints, list(a_coef + b_coef <= 0))
        }
      }
    }
  }

  # 2. Convexity constraints (second derivative constant)
  #"contraintes convexes



  # Second derivative is constant: P''(u) = c
  if (any(convcons != 0)) {
    for (i in 1:kn) {
      if (convcons[i] != 0) {
        # Second derivative value on interval i
        if (verbose) {
          message("applying convexity ",
                  convcons[i],
                  " at knot No",
                  i,
                  ":",
                  knot[i])
        }
        s2_val <- sum(alpha * deriv2_coeffs[i, ])

        if (convcons[i] > 0) {
          constraints <- c(constraints, list(s2_val >= 0))
        } else {
          constraints <- c(constraints, list(s2_val <= 0))
        }
      }
    }
  }

  # Solve the problem
  problem <- Problem(objective, constraints)

  result <- NULL
  solvers_to_try <- c(solver, "CLARABEL", "HIGHS", "OSQP", "ECOS", "SCS")

  for (s in unique(solvers_to_try)) {
    if (verbose)
      cat("Trying solver:", s, "\n")

    # Use new 'CVXR' syntax: psolve() for optimal value
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

#' Wrapper function for quadratic spline quantile regression
#' Previous version notation style
#' This is a convenience wrapper for SplineQuadraticQuant.
#'
#' @inheritParams SplineQuadraticQuant
#' @return Same as SplineQuadraticQuant
#' @export
SplineConstQuantRegBs2 <- function(xtab,
                                   ytab,
                                   knot,
                                   tau,
                                   monot = 0,
                                   convcons = 0,
                                   solver = "CLARABEL",
                                   weight = NULL,
                                   verbose = FALSE) {
  SplineQuadraticQuant(
    xtab,
    ytab,
    knot,
    tau,
    monot = monot,
    convcons = convcons,
    solver = solver,
    weight = weight,
    verbose = verbose
  )
}
