# Quartic Spline Quantile Regression with Shape Constraints
# Based on the Python implementation in quantile_reg.py
# Author: Alexandre Abbes
# Updated for CVXR new syntax (psolve, value, status)

#' Quantile regression with quartic splines and shape constraints
#'
#' Performs quantile regression using quartic (degree 4) B-splines with
#' monotonicity (Karlin-Studden constraints on the cubic derivative),
#' convexity (Karlin-Studden constraints on the quadratic second derivative),
#' and third derivative constraints (linear at knot).
#'
#' @param xtab Predictor vector (x)
#' @param ytab Response vector (y)
#' @param knot Knot vector or number of knot
#' @param tau Quantile (between 0 and 1)
#' @param monot Monotonicity constraint vector per interval:
#'        1 = increasing, -1 = decreasing, 0 = unconstrained
#' @param convcons Convexity constraint vector per interval:
#'        1 = convex, -1 = concave, 0 = unconstrained
#' @param der3cons Third derivative constraint vector at knot:
#'        1 = positive third derivative, -1 = negative, 0 = unconstrained
#' @param solver CVXR solver to use (default = "OSQP")
#' @param weight Observation weights (default = 1 for all)
#' @param verbose Logical; if TRUE, print progress messages
#' @return A list containing coefficients, degree, and knot
#' @export
SplineQuarticQuant <- function(xtab, ytab, knot, tau,
                               monot = 0,
                               convcons = 0,
                               der3cons = 0,
                               solver = "CLARABEL",
                               weight = NULL,
                               verbose = FALSE) {

  if (is.null(weight)) {
    weight <- rep(1, length(xtab))
  }

  # Sort data
  ordre <- order(xtab)
  xtab <- xtab[ordre]
  ytab <- ytab[ordre]
  weight <- weight[ordre]

  n <- length(xtab)

  # Handle knot
  if (length(knot) == 1 && is.numeric(knot)) {
    kn <- knot - 1
    knot <- quantile(xtab, probs = seq(0, 1, length.out = kn + 1))
  }

  kn <- length(knot) - 1
  degree <- 4
  N <- kn + degree  # Number of basis functions for quartic = kn + 4

  if (verbose) {
    cat("=== Quartic Spline Quantile Regression (degree = 4) ===\n")
    cat(sprintf("knot: %d, Basis functions: %d\n", kn, N))
  }

  # Handle constraints
  if (length(monot) == 1) {
    monot <- rep(monot, kn)
  }

  if (length(convcons) == 1) {
    convcons <- rep(convcons, kn)
  }

  if (length(der3cons) == 1) {
    der3cons <- rep(der3cons, kn + 1)
  }

  if (verbose) {
    cat("Monotonicity constraints:", monot, "\n")
    cat("Convexity constraints:", convcons, "\n")
    cat("3rd derivative constraints:", der3cons, "\n")
  }

  # Build B-spline basis and derivative coefficients
  deriv_data <- bspline_to_deriv_coeffs_quart(knot, degree = 4, x_values = xtab)

  B <- deriv_data$d0
  B <- t(B)  # Design matrix: n x N

  # Derivative coefficients
  deriv1_coeffs <- deriv_data$d1  # [kn, N, 4] for cubic derivative
  deriv2_coeffs <- deriv_data$d2  # [kn, N, 3] for quadratic second derivative
  deriv3_knot <- deriv_data$d3   # [kn+1, N] for linear third derivative at knot

  # Center data
  y_mean <- mean(ytab)
  ytab_centered <- ytab - y_mean

  # Optimization variables
  alpha <- Variable(N)
  z_vars <- list()

  # Objective function
  residuals <- ytab_centered - B %*% alpha
  weighted_loss <- sum(weight * (tau * pos(residuals) + (1 - tau) * pos(-residuals)))
  objective <- Minimize(weighted_loss)

  constraints <- list()

  # 1. Monotonicity constraints (Karlin-Studden on cubic derivative)
  # Derivative is cubic: P'(u) = a*u^3 + b*u^2 + c*u + d
  if (any(monot != 0)) {
    # Create z variables for cubic constraints (one per interval)
    z0_vars <- list()
    z1_vars <- list()

    for (i in 1:kn) {
      if (monot[i] != 0) {
        z0_vars[[i]] <- Variable(1, name = paste0("z0_monot_", i))
        z1_vars[[i]] <- Variable(1, name = paste0("z1_monot_", i))

        # Coefficients of the derivative on interval i
        # deriv1_coeffs[i, j, ] = [a3, a2, a1, a0] for a3*u^3 + a2*u^2 + a1*u + a0
        a3 <- sum(alpha * deriv1_coeffs[i, , 1])
        a2 <- sum(alpha * deriv1_coeffs[i, , 2])
        a1 <- sum(alpha * deriv1_coeffs[i, , 3])
        a0 <- sum(alpha * deriv1_coeffs[i, , 4])

        # Apply sign of monotonicity
        s <- sign(monot[i])

        # Karlin-Studden constraints for cubic
        karlin_constr <- apply_karlin_cubic(a3, a2, a1, a0,
                                            z0_vars[[i]], z1_vars[[i]],
                                            sign = s)
        constraints <- c(constraints, karlin_constr)
      }
    }
  }

  # 2. Convexity constraints (Karlin-Studden on quadratic second derivative)
  # Second derivative is quadratic: P''(u) = a*u^2 + b*u + c
  if (any(convcons != 0)) {
    z_conv_vars <- list()

    for (i in 1:kn) {
      if (convcons[i] != 0) {
        z_conv_vars[[i]] <- Variable(1, name = paste0("z_conv_", i))

        # Coefficients of the second derivative on interval i
        # deriv2_coeffs[i, j, ] = [a2, a1, a0] for a2*u^2 + a1*u + a0
        a2 <- sum(alpha * deriv2_coeffs[i, , 1])
        a1 <- sum(alpha * deriv2_coeffs[i, , 2])
        a0 <- sum(alpha * deriv2_coeffs[i, , 3])

        # Apply sign of convexity
        s <- sign(convcons[i])

        # Karlin-Studden constraints for quadratic
        karlin_constr <- apply_karlin_quadratic(a2, a1, a0,
                                                z_conv_vars[[i]],
                                                sign = s)
        constraints <- c(constraints, karlin_constr)
      }
    }
  }

  # 3. Third derivative constraints (linear at knot)
  # For quartic splines, third derivative is affine on each interval
  # We impose sign constraints at knot
  if (any(der3cons != 0)) {
    if (verbose) { cat(c("On applique les contraintes d'ordre 3:", der3cons,"\n"))}
    for (i in 1:(kn + 1)) {
      if (der3cons[i] != 0) {
        # Third derivative value at knot i
        s3_val <- sum(alpha * deriv3_knot[i, ])
        lin_constr <- apply_linear_constraint(s3_val, sign = der3cons[i])
        constraints <- c(constraints, lin_constr)
      }
    }
  }

  # Solve the problem using new CVXR syntax
  problem <- Problem(objective, constraints)

  result <- NULL
  solvers_to_try <- c(solver, "CLARABEL", "OSQP", "ECOS", "SCS")

  for (s in unique(solvers_to_try)) {
    if (verbose) cat("Trying solver:", s, "\n")

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
      if (verbose) cat("Failed:", e$message, "\n")
      NULL
    })

    if (!is.null(result) && !is.null(result$alpha_value)) {
      if (verbose) cat("Solver succeeded:", s, "\n")
      break
    }
  }

  if (is.null(result) || is.null(result$alpha_value)) {
    warning("Optimization did not converge with any solver")
    return(NULL)
  }

  alpha_val <- result$alpha_value + y_mean
  result$y_mean<-y_mean
  if (verbose) {
    cat("Status:", result$status, "\n")
    cat("Objective value:", result$value, "\n")
  }

  # Return results
  return(list(
    coefficients = alpha_val,
    degree = degree,
    knot = t(knot),
    result<-result
  ))
}

#' Wrapper function for quartic spline quantile regression
#'
#' This is a convenience wrapper for SplineQuarticQuant
#' that follows the same naming convention as SplineConstQuantRegBs3.
#'
#' @inheritParams SplineQuarticQuant
#' @return Same as SplineQuarticQuant
#' @export
SplineConstQuantRegBs4 <- function(xtab, ytab, knot, tau,
                                   monot = 0,
                                   convcons = 0,
                                   der3cons = 0,
                                   solver = "CLARABEL",
                                   weight = NULL,
                                   verbose = FALSE) {
  SplineQuarticQuant(xtab, ytab, knot, tau,
                     monot = monot,
                     convcons = convcons,
                     der3cons = der3cons,
                     solver = solver,
                     weight = weight,
                     verbose = verbose)
}
