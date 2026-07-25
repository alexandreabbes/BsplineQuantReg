#' Apply Karlin-Studden constraints for a cubic polynomial
#'
#' For a cubic polynomial p(u) = a*u^3 + b*u^2 + c*u + d on [0,1],
#' this function applies the 'Karlin-Studden SOCP constraints'
#' to ensure p(u) >= 0 or p(u) <= 0 on [0,1].
#'
#' @param p3 Coefficient of u^3
#' @param p2 Coefficient of u^2
#' @param p1 Coefficient of u
#' @param p0 Constant term
#' @param z0 Auxiliary SOCP variable for cubic
#' @param z1 Auxiliary SOCP variable for cubic
#' @param sign Sign of the constraint (+1 for >= 0, -1 for <= 0)
#' @return List of 'CVXR' constraints
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
#' @return List of 'CVXR' constraints
#' @keywords internal*

apply_karlin_quadratic <- function(p2, p1, p0, z0, sign = 1)
  {
  constraints <- list()

  # z0 >= 0
  constraints <- c(constraints, list(z0 >= 0))

  if (sign < 0) {
    p2 <- -p2
    p1 <- -p1
    p0 <- -p0
  }

   #Karlin-Studden characterization for quadratic polynomials)
   #p(u) = a*u^2 + b*u + c
   #Condition: there exists z0 >= 0 such that
  # (p0 + p2 + z0, p0 - p2 - z0, p1 - z0) in Q3)

  K1_x <- p0 + p2 + z0
  K1_y <- p0 - p2 - z0
  K2 <- p1 - z0

  vec <- vstack(K1_y, K2)
  constraints <- c(constraints, list(K1_x >= p_norm(vec, 2)))

  return(constraints)
}

#' Apply linear constraints at knot for the third derivative
#'
#' For a quartic spline, the third derivative is affine (linear) on each interval.
#' This function applies sign constraints at the knot.
#'
#' @param const_value Value of the third derivative at a knot
#' @param sign Sign of the constraint (+1 for >= 0, -1 for <= 0)
#' @return List of 'CVXR' constraints
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
