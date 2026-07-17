# Quartic Spline Quantile Regression with Shape Constraints
# Based on the Python implementation in quantile_reg.py
# Author: Alexandre Abbes

#' Apply Karlin-Studden constraints for a cubic polynomial
#'
#' For a cubic polynomial p(u) = a*u^3 + b*u^2 + c*u + d on [0,1],
#' this function applies the Karlin-Studden SOCP constraints
#' to ensure p(u) >= 0 or p(u) <= 0 on [0,1].
#'
#' @param p3 Coefficient of u^3
#' @param p2 Coefficient of u^2
#' @param p1 Coefficient of u
#' @param p0 Constant term
#' @param z0 Auxiliary SOCP variable for cubic
#' @param z1 Auxiliary SOCP variable for cubic
#' @param sign Sign of the constraint (+1 for >= 0, -1 for <= 0)
#' @return List of CVXR constraints
#' @keywords internal
apply_karlin_cubic <- function(p3, p2, p1, p0, z0, z1, sign = 1) {
  constraints <- list()

  # z0 >= 0, z1 >= 0
  constraints <- c(constraints, list(z0 >= 0))
  constraints <- c(constraints, list(z1 >= 0))

  # Apply sign transformation if negative
  if (sign < 0) {
    p3 <- -p3
    p2 <- -p2
    p1 <- -p1
    p0 <- -p0
  }

  # Karlin-Studden characterization for cubic polynomials
  # p(u) = a*u^3 + b*u^2 + c*u + d
  # There exist x0, x1, x2, y0, y1, y2 such that:
  # d = y0
  # c = 2*y1 + x0 - y0
  # b = y2 + 2*x1 - 2*y1
  # a = x2 - y2
  # with SOC constraints:
  # (x0+x2, x0-x2, 2*x1) in Q3
  # (y0+y2, y0-y2, 2*y1) in Q3

  # Variables for the representation
  x0 <- Variable(1)
  x1 <- Variable(1)
  x2 <- Variable(1)
  y0 <- Variable(1)
  y1 <- Variable(1)
  y2 <- Variable(1)

  # Equations (6a)-(6d) from Karlin-Studden
  constraints <- c(constraints, list(p0 == y0))
  constraints <- c(constraints, list(p1 == 2*y1 + x0 - y0))
  constraints <- c(constraints, list(p2 == y2 + 2*x1 - 2*y1))
  constraints <- c(constraints, list(p3 == x2 - y2))

  # SOC constraints (6e)-(6f)
  vec_x <- vstack(x0 - x2, 2*x1)
  constraints <- c(constraints, list(x0 + x2 >= p_norm(vec_x, 2)))

  vec_y <- vstack(y0 - y2, 2*y1)
  constraints <- c(constraints, list(y0 + y2 >= p_norm(vec_y, 2)))

  return(constraints)
}

#' Apply Karlin-Studden constraints for a quadratic polynomial
#'
#' For a quadratic polynomial p(u) = a*u^2 + b*u + c on [0,1],
#' this function applies the Karlin-Studden SOCP constraints
#' to ensure p(u) >= 0 or p(u) <= 0 on [0,1].
#'
#' @param p2 Coefficient of u^2
#' @param p1 Coefficient of u
#' @param p0 Constant term
#' @param z0 Auxiliary SOCP variable
#' @param sign Sign of the constraint (+1 for >= 0, -1 for <= 0)
#' @return List of CVXR constraints
#' @keywords internal
apply_karlin_quadratic <- function(p2, p1, p0, z0, sign = 1) {
  constraints <- list()

  # z0 >= 0
  constraints <- c(constraints, list(z0 >= 0))

  if (sign < 0) {
    p2 <- -p2
    p1 <- -p1
    p0 <- -p0
  }

  # Karlin-Studden characterization for quadratic polynomials
  # p(u) = a*u^2 + b*u + c
  # Condition: there exists z0 >= 0 such that
  # (p0 + p2 + z0, p0 - p2 - z0, p1 - z0) in Q3

  K1_x <- p0 + p2 + z0
  K1_y <- p0 - p2 - z0
  K2 <- p1 - z0

  vec <- vstack(K1_y, K2)
  constraints <- c(constraints, list(K1_x >= p_norm(vec, 2)))

  return(constraints)
}

#' Apply linear constraints at knots for the third derivative
#'
#' For a quartic spline, the third derivative is affine (linear) on each interval.
#' This function applies sign constraints at the knots.
#'
#' @param const_value Value of the third derivative at a knot
#' @param sign Sign of the constraint (+1 for >= 0, -1 for <= 0)
#' @return List of CVXR constraints
#' @keywords internal
apply_linear_constraint <- function(const_value, sign = 1) {
  constraints <- list()

  if (sign > 0) {
    constraints <- c(constraints, list(const_value >= 0))
  } else if (sign < 0) {
    constraints <- c(constraints, list(const_value <= 0))
  }

  return(constraints)
}

#' Convert quartic B-spline to derivative coefficients
#'
#' Computes normalized first, second, and third derivative coefficients
#' for quartic B-splines on each interval.
#'
#' @param knots Knot vector (effective partition)
#' @param degree Spline degree (should be 4)
#' @param xvalues Evaluation points for design matrix
#' @return A list containing:
#'   \item{d0}{Design matrix (if xvalues provided)}
#'   \item{d1}{First derivative coefficients [a3, a2, a1, a0] for each interval}
#'   \item{d2}{Second derivative coefficients [a2, a1, a0] for each interval}
#'   \item{d3}{Third derivative values at knots (linear constraints)}
#' @export
bspline_to_deriv_coeffs_quart <- function(knots, degree = 4, xvalues = 0) {

  kn <- length(knots) - 1
  N <- kn + degree

  # Extended knot vector
  sn <- c(rep(knots[1], degree), knots, rep(knots[kn + 1], degree))

  # Build B-spline basis
  BB <- Bspline_base(sn, degree = degree)
  basis <- BB$base

  # Number of intervals in extended notation
  n_intervals <- length(sn) - 1

  # Derivative coefficients
  # First derivative: cubic on each interval -> [a3, a2, a1, a0] for a3*u^3 + a2*u^2 + a1*u + a0
  deriv1_coeffs <- array(0, dim = c(kn, N, 4))

  # Second derivative: quadratic on each interval -> [a2, a1, a0] for a2*u^2 + a1*u + a0
  deriv2_coeffs <- array(0, dim = c(kn, N, 3))

  # Third derivative at knots: linear constraints
  # For quartic splines, third derivative is affine on each interval
  # We evaluate at knots for linear constraints
  deriv3_knots <- array(0, dim = c(kn + 1, N))

  for (j in 1:N) {
    for (nu in (degree + 1):(kn + degree)) {
      h <- sn[nu + 1] - sn[nu]

      # basis[j, nu, ] = [a4, a3, a2, a1, a0] for polynomial on interval
      # P(u) = a0 + a1*u + a2*u^2 + a3*u^3 + a4*u^4
      a0 <- basis[j, nu, 5]
      a1 <- basis[j, nu, 4]
      a2 <- basis[j, nu, 3]
      a3 <- basis[j, nu, 2]
      a4 <- basis[j, nu, 1]

      # First derivative: P'(u) = a1 + 2*a2*u + 3*a3*u^2 + 4*a4*u^3
      # Normalized: a1 + (2*a2*h)*u + (3*a3*h^2)*u^2 + (4*a4*h^3)*u^3
      deriv1_coeffs[nu - degree, j, ] <- c(
        4 * a4 * h^3,
        3 * a3 * h^2,
        2 * a2 * h,
        a1
      )

      # Second derivative: P''(u) = 2*a2 + 6*a3*u + 12*a4*u^2
      # Normalized: (2*a2) + (6*a3*h)*u + (12*a4*h^2)*u^2
      deriv2_coeffs[nu - degree, j, ] <- c(
        12 * a4 * h^2,
        6 * a3 * h,
        2 * a2
      )
    }

    # Third derivative at knots (for linear constraints)
    # For quartic splines, third derivative is affine on each interval
    # We evaluate at knots for linear constraints
    for (i in 1:(kn + 1)) {
      # Evaluate at knot i
      if (i <= kn) {
        # Use polynomial on interval i
        nu <- degree + i
        h <- sn[nu + 1] - sn[nu]
        a3 <- basis[j, nu, 2]
        a4 <- basis[j, nu, 1]
        # P'''(x) = 6*a3 + 24*a4*(x - t_k)
        # At knot i (x = t_k): P'''(t_k) = 6*a3
        deriv3_knots[i, j] <- 6 * a3
      } else {
        # Last knot: use polynomial on last interval
        nu <- degree + kn
        a3 <- basis[j, nu, 2]
        deriv3_knots[i, j] <- 6 * a3
      }
    }
  }

  # Design matrix if xvalues provided
  if (length(xvalues) != 1) {
    yvalues <- bs_direct(BB, xvalues)
  } else {
    yvalues <- 0
  }

  return(list(
    d0 = yvalues,
    d1 = deriv1_coeffs,
    d2 = deriv2_coeffs,
    d3 = deriv3_knots
  ))
}

#' Quantile regression with quartic splines and shape constraints
#'
#' Performs quantile regression using quartic (degree 4) B-splines with
#' monotonicity (Karlin-Studden constraints on the cubic derivative),
#' convexity (Karlin-Studden constraints on the quadratic second derivative),
#' and third derivative constraints (linear at knots).
#'
#' @param xtab Predictor vector (x)
#' @param ytab Response vector (y)
#' @param knots Knot vector or number of knots
#' @param tau Quantile (between 0 and 1)
#' @param monot Monotonicity constraint vector per interval:
#'        1 = increasing, -1 = decreasing, 0 = unconstrained
#' @param convcons Convexity constraint vector per interval:
#'        1 = convex, -1 = concave, 0 = unconstrained
#' @param der3cons Third derivative constraint vector at knots:
#'        1 = positive third derivative, -1 = negative, 0 = unconstrained
#' @param solver CVXR solver to use (default = "OSQP")
#' @param weight Observation weights (default = 1 for all)
#' @param verbose Logical; if TRUE, print progress messages
#' @return A list containing coefficients, degree, and knots
#' @export
SplineQuarticQuant <- function(xtab, ytab, knots, tau,
                               monot = 0,
                               convcons = 0,
                               der3cons = 0,
                               solver = "OSQP",
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

  # Handle knots
  if (length(knots) == 1 && is.numeric(knots)) {
    kn <- knots - 1
    knots <- quantile(xtab, probs = seq(0, 1, length.out = kn + 1))
  }

  kn <- length(knots) - 1
  degree <- 4
  N <- kn + degree  # Number of basis functions for quartic = kn + 4

  if (verbose) {
    cat("=== Quartic Spline Quantile Regression (degree = 4) ===\n")
    cat(sprintf("Knots: %d, Basis functions: %d\n", kn, N))
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
  deriv_data <- bspline_to_deriv_coeffs_quart(knots, degree = 4, xvalues = xtab)

  B <- deriv_data$d0
  B <- t(B)  # Design matrix: n x N

  # Derivative coefficients
  deriv1_coeffs <- deriv_data$d1  # [kn, N, 4] for cubic derivative
  deriv2_coeffs <- deriv_data$d2  # [kn, N, 3] for quadratic second derivative
  deriv3_knots <- deriv_data$d3   # [kn+1, N] for linear third derivative at knots

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

  # 3. Third derivative constraints (linear at knots)
  # For quartic splines, third derivative is affine on each interval
  # We impose sign constraints at knots
  if (any(der3cons != 0)) {
    for (i in 1:(kn + 1)) {
      if (der3cons[i] != 0) {
        # Third derivative value at knot i
        s3_val <- sum(alpha * deriv3_knots[i, ])
        lin_constr <- apply_linear_constraint(s3_val, sign = der3cons[i])
        constraints <- c(constraints, lin_constr)
      }
    }
  }

  # Solve the problem
  problem <- Problem(objective, constraints)

  result <- NULL
  solvers_to_try <- c(solver, "OSQP", "ECOS", "SCS")

  for (s in unique(solvers_to_try)) {
    if (verbose) cat("Trying solver:", s, "\n")
    result <- tryCatch({
      solve(problem, solver = toupper(s), verbose = FALSE)
    }, error = function(e) {
      if (verbose) cat("Failed:", e$message, "\n")
      NULL
    })

    if (!is.null(result) && !is.null(result$getValue(alpha))) {
      if (verbose) cat("Solver succeeded:", s, "\n")
      break
    }
  }

  if (is.null(result) || is.null(result$getValue(alpha))) {
    warning("Optimization did not converge with any solver")
    return(NULL)
  }

  alpha_val <- result$getValue(alpha) + y_mean

  if (verbose) {
    cat("Status:", result$status, "\n")
    cat("Objective value:", result$value, "\n")
  }

  # Return results
  return(list(
    coefficients = alpha_val,
    degree = degree,
    knots = knots,
    int_knots= knots[2:kn],
    y_mean = y_mean,
    status = result$status,
    value = result$value
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
SplineConstQuantRegBs4 <- function(xtab, ytab, knots, tau,
                                   monot = 0,
                                   convcons = 0,
                                   der3cons = 0,
                                   solver = "OSQP",
                                   weight = NULL,
                                   verbose = FALSE) {
  SplineQuarticQuant(xtab, ytab, knots, tau,
                     monot = monot,
                     convcons = convcons,
                     der3cons = der3cons,
                     solver = solver,
                     weight = weight,
                     verbose = verbose)
}
