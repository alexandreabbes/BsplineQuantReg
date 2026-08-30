# Polynomial Spline Quantile Regression with Shape Constraints
# Based on the Python implementation in quantile_reg.py
# Author: Alexandre Abbes
# Updated for 'CVXR' new syntax (psolve, value, status)

#' Quantile regression with general polynomial splines and shape constraints
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
#' @param degree the degree of the polynomial spline degree>=5
#' @param convcons not yet implemented Convexity constraint vector per interval:
#'        1 = convex, -1 = concave, 0 = unconstrained
#' @param der3cons Not yet implemented Third derivative constraint vector at knot:
#'        1 = positive third derivative, -1 = negative, 0 = unconstrained
#' @param solver 'CVXR' solver to use (default = "OSQP")
#' @param weight Observation weights (default = 1 for all)
#' @param verbose Logical; if TRUE, print progress messages
#' @param type_reg 'quantile' or 'mean_square' type of regression
#' @return A list containing coefficients, degree, and knot
#' @export
SplinePolynQuant <- function(xtab,
                               ytab,
                               knot=NULL,
                               convcons=0,
                                der3cons=0,
                               tau=0.5,
                               degree=5,
                               solver = "CLARABEL",
                               weight = NULL,
                               verbose = FALSE,
                               type_reg='quantile')
{
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

  # Handle knot
  if (length(knot) == 1 && is.numeric(knot)) {
    kn <- knot - 1
    knot <- as.numeric(quantile(xtab, probs = seq(0, 1, length.out = kn + 1)))
  }

  kn <- length(knot) - 1

  N <- kn + degree  # Number of basis functions for quartic = kn + 4

  if (verbose) {
    cat("=== General  Spline Quantile Regression (degree = ",degree,") ===\n")
    cat(sprintf("knot: %d, Basis functions: %d\n", kn, N))
  }


  # Build B-spline basis and derivative coefficients
  sn <- c(rep(knot[1], degree), knot, rep(knot[kn + 1], degree))

  BB <- Bspline_base(sn, degree = degree)
  B<-bs_direct(BB,xtab)
  B <- t(B)  # Design matrix: n x N
  # Center data
  y_mean <- mean(ytab)
  ytab_centered <- ytab - y_mean

  # Optimization variables
  alpha <- Variable(N)
  z_vars <- list()

  # Objective function
  residuals <- ytab_centered - B %*% alpha

  if (type_reg=='mean_square'){
    weighted_loss <- p_norm(residuals,2)
  }
  else{#if (type_reg=='quantile'){
    #quantile
    u_plus <- pos(residuals)
    u_minus <- pos(-residuals)
    weighted_loss <- sum(weight * (u_plus*tau+u_minus*(1-tau)))
  }
  objective <- Minimize(weighted_loss)



# Handle constraints
  constraints <- list()
if (degree==5){
  #"contraintes convexes

  # eliminate the null (unconstrained) case
  if (any(convcons != 0)) {
    if (length(convcons) == 1) {
      convcons <- rep(convcons, (kn + 1))
    }
    if (length(convcons) < (kn + 1)) {
      message("Not enough convexity constraints, completing with 0")
      convcons <- c(convcons, rep(0, kn + 1 - length(convcons)))
    }
    if (verbose) {
      message("Convexity constraints (Linear):", convcons, "\n")
    }
    }

  #3rd derivative constraints
  if (any(der3cons != 0)) {
    #if (dim(deriv_coeffs3)[2]!=N){
    #prepare the values for matrix mult.      #calculation
    # deriv_coeff3=t(deriv_coeffs3)}
    if (length(der3cons) == 1) {
      der3cons <- rep(der3cons, kn)
    }
    if (length(der3cons) < kn+1) {
      message("not enough 3rd order constraints, completing with 0")
      der3cons = c(der3cons, rep(0, kn+1 - length(der3cons)))
    }
    if (verbose) {
      message("3rd order derivative constraints (constant):",
              der3cons,
              "\n")
    }
  }



  deriv_data <- bspline_to_deriv_coeffs_general(BB)

  # Derivative coefficients
  #deriv1_coeffs <- deriv_data$d1  # [kn, N, 4] for cubic derivative
  deriv2_coeffs <- deriv_data$d2  # [kn, N, 3] for quadratic second derivative
  deriv3_coeffs <- deriv_data$d3   # [kn+1, N] for linear third derivative at knot




  # 1. convexity constraints (Karlin-Studden on cubic derivative)
  # Second Derivative is cubic: P'(u) = a*u^3 + b*u^2 + c*u + d
  if (any(convcons!= 0)) {
    cat("constraints are applyed")
    # Create z variables for cubic constraints (one per interval)
    z0_vars <- list()
    z1_vars <- list()

    for (i in 1:kn) {
      if (verbose){cat("Applying convexity on itervall ",i,"\n") }
      if (convcons[i] != 0) {
        z0_vars[[i]] <- Variable(1, name = paste0("z0_convcons_", i))
        z1_vars[[i]] <- Variable(1, name = paste0("z1_convcons_", i))

        #Coefficients of the derivative on interval i
        #deriv2_coeffs[i, j, ] <- [a3, a2, a1, a0]  #for a3*u^3 + a2*u^2 + a1*u + a0
        a3 <- sum(alpha * deriv2_coeffs[i, , 1])
        a2 <- sum(alpha * deriv2_coeffs[i, , 2])
        a1 <- sum(alpha * deriv2_coeffs[i, , 3])
        a0 <- sum(alpha * deriv2_coeffs[i, , 4])

        # Apply sign of monotonicity
        s <- sign(convcons[i])

        # Karlin-Studden constraints for cubic
        karlin_constr <- apply_karlin_cubic(a3, a2, a1, a0, z0_vars[[i]], z1_vars[[i]], sign = s)
        constraints <- c(constraints, karlin_constr)
      }
    }
  }

  #

  # 2. Third derivative constraints (Karlin-Studden on quadratic second derivative)
  # Second derivative is quadratic: P''(u) = a*u^2 + b*u + c
  if (any(der3cons != 0)) {
    z_der3_vars <- list()

    for (i in 1:kn) {
       if (der3cons[i] != 0) {
         z_der3_vars[[i]] <- Variable(1, name = paste0("z_conv_", i))

         # Coefficients of the second derivative on interval i
         # deriv2_coeffs[i, j, ] = [a2, a1, a0] for a2*u^2 + a1*u + a0
         a2 <- sum(alpha * deriv3_coeffs[i, , 1])
         a1 <- sum(alpha * deriv3_coeffs[i, , 2])
         a0 <- sum(alpha * deriv3_coeffs[i, , 3])

  #       # Apply sign of convexity
         s <- sign(der3cons[i])
  #
  #       # Karlin-Studden constraints for quadratic
         karlin_constr <- apply_karlin_quadratic(a2, a1, a0, z_der3_vars[[i]], sign = s)
         constraints <- c(constraints, karlin_constr)
       }
     }
   }
}
  # Solve the problem using new 'CVXR' syntax

  problem <- Problem(objective, constraints)

  result <- NULL
  solvers_to_try <- c(solver, "CLARABEL", "HIGHS", "ECOS", "OSQP", "SCS", "MOSEK")

  for (s in unique(solvers_to_try)) {
    if (verbose)
      cat("Trying solver:", s, "\n")

    # Use new 'CVXR' syntax: psolve()
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

    if (!is.null(result) && !is.null(result$alpha_value)) {
      if (verbose)
        cat("Solver succeeded:", s, "\n")
      break
    }
  }

  if (is.null(result) || is.null(result$alpha_value)) {
    warning("Optimization did not converge with any solver")
    return(NULL)
  }

  alpha_val <- result$alpha_value + y_mean
  result$y_mean <- y_mean
  if (verbose) {
    cat("Status:", result$status, "\n")
    cat("Objective value:", result$value, "\n")
  }

  # Return results
  return(list(
    coeff = alpha_val,
    degree = degree,
    knot = t(knot),
    result <- result
  ))
}

