#' normalized first and second derivative coefficients on each interval.
#' for a cubic B-spline basis
#'
#' @param tn Knot vector (effective partition, not extended. OPtionnal if Bsbasis
#' is given)
#' @param degree Spline degree (default = 3)
#' @param x_values Evaluation points for design matrix (0 = no evaluation)
#' @param Bsbasis Bspline basis structure if already computed. Optional if knots are given.
#' @param verbose boolean FALSE (default) or TRUE.
#' @return A list containing (kn: Nb intervals, N=kn+3: Nb basis functions):
#'   \item{d0}{Design matrix (if x_values provided)}
#'   \item{d1}{First derivative coefficients [a3, a2, a1] for each interval: (kn x N x 3) array}
#'   \item{d2}{Second derivative values at knot: ((kn+1) x N) matrix}
#'   \item{d3}{Third derivative values at knots: (kn x N) matrix}
#' @export
bspline_to_deriv_coeffs_cubic <- function(tn=NULL,
                                          degree = 3,
                                          x_values = 0,
                                          Bsbasis=NULL,
                                          verbose = FALSE) {
if (!is.null(Bsbasis)){
  tn<-Bsbasis$knot
  degree=Bsbasis$degree}
else {if (is.null(tn)){
  cat('Provide at least the knots')
  return()
}
  # create  basis with create.bspline.basis
  kn <- length(tn) - 1
  sn = c(tn[1] * rep(1, degree), tn, tn[kn + 1] * rep(1, degree)) # this is the extended knot sequence
  Bsbasis <- Bspline_base(sn, degree)
}
  basis <- Bsbasis$base

  N <- Bsbasis$n_splines
  if (verbose) {
    message("Number of  basis functions", N, "\n")
  }
  # Matrix of normalised coefficient derivativs
  deriv_coeffs <- array(0, dim = c(kn, N, 3))
  deriv2_val <- array(0, dim = c(kn + 1, N))
  deriv3_val <- array(0, dim = c(kn, N))
  for (j in 1:N) {
    for (nu in (degree + 1):(kn + degree))
    {
      h = sn[nu + 1] - sn[nu]
      c3 <- basis[j, nu, 1]
      c2 <- basis[j, nu, 2]
      a3 <- 3 * c3 * h^2
      a2 <- 2 * basis[j, nu, 2] * h
      a1 <- basis[j, nu, 3]
      # coeffs_poly est [a3, a2, a1, a0] a0+a1*x+a_2*x^2+a3*x^3
      deriv_coeffs[nu - degree, j, ] <- c(a3, a2, a1)
      deriv2_val[nu - degree, j] <- c2
      deriv3_val[nu - degree, j] <- c3 #up to a factor 6, but the sign is the same.
    }

    # for the last knot the second deriv is an affine function
    # p+m(t-t_{kn-1}) h is the last intervall space
    p <- 2 * basis[j, nu, 2]
    m <- 6 * basis[j, nu, 1]
    deriv2_val[nu - degree + 1, j] = p + m * h
  }
  if (length(x_values) != 1) {
    yvalues = bs_direct(Bsbasis, x_values)
  }
  else {
    yvalues = 0
  }
  return(list(
    d0 = yvalues,
    d1 = deriv_coeffs,
    d2 = deriv2_val,
    d3 = deriv3_val
  ))
}



#' Computes normalized first, second, and third derivative coefficients
#' for quartic B-splines on each interval.
#'
#' @param knot Knot vector (effective partition)
#' @param degree Spline degree (should be 4)
#' @param x_values Evaluation points for design matrix
#' @return A list containing:
#'   \item{d0}{Design matrix (if x_values provided)}
#'   \item{d1}{First derivative coefficients [a3, a2, a1, a0] for each interval}
#'   \item{d2}{Second derivative coefficients [a2, a1, a0] for each interval}
#'   \item{d3}{Third derivative values at knot (linear constraints)}
#' @export
bspline_to_deriv_coeffs_quart <- function(knot,
                                          degree = 4,
                                          x_values = 0) {
  kn <- length(knot) - 1
  N <- kn + degree

  # Extended knot vector
  sn <- c(rep(knot[1], degree), knot, rep(knot[kn + 1], degree))

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

  # Third derivative at knot: linear constraints
  # For quartic splines, third derivative is affine on each interval
  # We evaluate at knot for linear constraints
  deriv3_knot <- array(0, dim = c(kn + 1, N))

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
      deriv1_coeffs[nu - degree, j, ] <- c(4 * a4 * h^3, 3 * a3 * h^2, 2 * a2 * h, a1)

      # Second derivative: P''(u) = 2*a2 + 6*a3*u + 12*a4*u^2
      # Normalized: (2*a2) + (6*a3*h)*u + (12*a4*h^2)*u^2
      deriv2_coeffs[nu - degree, j, ] <- c(12 * a4 * h^2, 6 * a3 * h, 2 * a2)
    }

    # Third derivative at knot (for linear constraints)
    # For quartic splines, third derivative is affine on each interval
    # We evaluate at knot for linear constraints
    for (i in 1:(kn + 1)) {
      # Evaluate at knot i
      if (i <= kn) {
        # Use polynomial on interval i
        nu <- degree + i
        h <- sn[nu + 1] - sn[nu]
        a3 <- basis[j, nu, 2]
        a4 <- basis[j, nu, 1]
        # P'''(x) = 6*a3 + 24*a4*(x - t_k)
        # At left knot i (x = t_k): P'''(t_k) = 6*a3
        deriv3_knot[i, j] <- 6 * a3
      } else {
        # Last knot: use polynomial on last interval
        #P'''(t_{kn+1}) = 6*a3 + 24*a4*(t_{kn+1} - t_k)
        nu <- degree + kn
        h <- sn[nu + 1] - sn[nu]
        a3 <- basis[j, nu, 2]
        a4 <- basis[j, nu, 1]
        deriv3_knot[i, j] <- 6 * a3 + 24 * a4 * h
      }
    }
  }

  # Design matrix if x_values provided
  if (length(x_values) != 1) {
    yvalues <- bs_direct(BB, x_values)
  } else {
    yvalues <- 0
  }

  return(list(
    d0 = yvalues,
    d1 = deriv1_coeffs,
    d2 = deriv2_coeffs,
    d3 = deriv3_knot
  ))
}



#' quadratic B-spline  derivative coefficients
#'
#' Computes normalized first and second derivative coefficients
#' for quadratic B-splines on each interval.
#'
#' @param tn Knot vector (effective partition)
#' @param degree Spline degree (should be 2)
#' @param x_values Evaluation points for design matrix (0 = no evaluation)
#' @param verbose logical; if TRUE, print progress messages
#' @return A list containing:
#'   \item{d0}{Design matrix (if x_values provided)}
#'   \item{d1}{First derivative coefficients [a1, a0] for each interval}
#'   \item{d2}{Second derivative values (constant per interval)}
#' @export
bspline_to_deriv_coeffs_quad <- function(tn,
                                         degree = 2,
                                         x_values = 0,
                                         verbose = FALSE) {
  kn <- length(tn) - 1
  N <- kn + degree  # Number of basis functions for quadratic = kn + 2

  # Extended knot vector
  sn <- c(rep(tn[1], degree), tn, rep(tn[kn + 1], degree))

  # Build B-spline basis
  BB <- Bspline_base(sn, degree = degree)
  basis <- BB$base

  if (verbose) {
    message("Quadratic B-spline basis: ",
            N,
            " basis functions, ",
            kn,
            " intervals")
  }

  # Derivative coefficients
  # First derivative: linear on each interval -> [a1, a0] for a1*u + a0
  deriv1_coeffs <- array(0, dim = c(kn, N, 2))

  # Second derivative: constant on each interval
  deriv2_coeffs <- array(0, dim = c(kn, N))

  for (j in 1:N) {
    for (nu in (degree + 1):(kn + degree)) {
      h <- sn[nu + 1] - sn[nu]

      # basis[j, nu, ] = [a2, a1, a0] for polynomial on interval
      # P(u) = a0 + a1*u + a2*u^2
      a1 <- basis[j, nu, 2]
      a2 <- basis[j, nu, 1]

      # First derivative: P'(u) = a1 + 2*a2*u with u in [sn_nu,sn_{nu+1}]
      # Normalized: a1 + (2*a2*h)*u on u in [0,1]
      deriv1_coeffs[nu - degree, j, ] <- c(2 * a2 * h, a1)

      # Second derivative: P''(u) = 2*a2 (constant)
      deriv2_coeffs[nu - degree, j] <- 2 * a2
    }
  }

  # Design matrix if x_values provided
  if (length(x_values) != 1) {
    yvalues <- bs_direct(BB, x_values)
  } else {
    yvalues <- 0
  }

  return(list(d0 = yvalues, d1 = deriv1_coeffs, d2 = deriv2_coeffs))
}


#' Derivative coefficients for linear B-spline
#' Computes normalized first derivative coefficients for linear B-splines.
#' For linear splines, the derivative is constant on each interval.
#'
#' @param tn Knot vector (effective partition)
#' @param degree Spline degree (should be 1)
#' @param x_values Evaluation points for design matrix (0 = no evaluation)
#' @param verbose logical; if TRUE, print progress messages
#' @return A list containing:
#'   \item{d0}{Design matrix (if x_values provided)}
#'   \item{d1}{First derivative values (constant per interval)}
#' @export
bspline_to_deriv_coeffs_lin <- function(tn,
                                        degree = 1,
                                        x_values = 0,
                                        verbose = FALSE) {
  kn <- length(tn) - 1
  N <- kn + degree  # Number of basis functions for linear = kn + 1

  # Extended knot vector
  sn <- c(rep(tn[1], degree), tn, rep(tn[kn + 1], degree))

  # Build B-spline basis
  BB <- Bspline_base(sn, degree = degree)
  basis <- BB$base

  if (verbose) {
    message("Linear B-spline basis: ",
            N,
            " basis functions, ",
            kn,
            " intervals")
  }

  # First derivative: constant on each interval
  # deriv1_coeffs[interval, basis] = constant value
  deriv1_coeffs <- array(0, dim = c(kn, N))

  for (j in 1:N) {
    for (nu in (degree + 1):(kn + degree)) {
      h <- sn[nu + 1] - sn[nu]

      # basis[j, nu, ] = [a1, a0] for polynomial on interval
      # P(u) = a0 + a1*u
      a1 <- basis[j, nu, 1]

      # First derivative: P'(u) = a1 (constant)
      # Normalized: a1/h on [0,1]
      deriv1_coeffs[nu - degree, j] <- a1 / h
    }
  }

  # Design matrix if x_values provided
  if (length(x_values) != 1) {
    yvalues <- bs_direct(BB, x_values)
  } else {
    yvalues <- 0
  }

  return(list(d0 = yvalues, d1 = deriv1_coeffs))
}
