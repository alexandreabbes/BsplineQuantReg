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
                                          verbose = FALSE)
  {
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
      a3 <- 3 * c3 * h^3
      a2 <- 2 * basis[j, nu, 2] * h^2
      a1 <- basis[j, nu, 3]*h
      # coeffs_poly est [a3, a2, a1, a0] a0+a1*x+a_2*x^2+a3*x^3
      deriv_coeffs[nu - degree, j, ] <- c(a3, a2, a1)
      deriv2_val[nu - degree, j] <- c2*h^2
      deriv3_val[nu - degree, j] <- c3*h #up to a factor 6, but the sign is the same.
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
  #deriv1_coeff <-Spline_der_knot(BB,der=1)
  deriv1_coeffs <- array(0, dim = c(kn, N))

  for (j in 1:N) {
    for (nu in (degree + 1):(kn + degree)) {
      h <- sn[nu + 1] - sn[nu]

      # basis[j, nu, ] = [a1, a0] for polynomial on interval
      # P(u) = a0 + a1*u
      a1 <- basis[j, nu, 1]

      # First derivative: P'(u) = a1 (constant)
      # Normalized: a1*h on [0,1]
      deriv1_coeffs[nu - degree, j] <- a1*h
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



#' Derivative Coefficients of a General B-Spline Basis
#'
#' For a B-spline basis of arbitrary degree \eqn{d \ge 1} built on a knot
#' vector \code{knot}, this helper returns the polynomial coefficients of the
#' second and third derivatives of every basis function on every knot
#' interval. These coefficients are the building blocks used by the
#' Karlin–Studden shape constraints in \code{\link{SplinePolynQuant}}.
#'
#' On each interval \eqn{[t_i, t_{i+1}]}, every B-spline basis function
#' \eqn{B_j} is a polynomial of degree \eqn{d}. Writing it in the local
#' normalized variable \eqn{u \in [0,1]} as
#' \SplinePolynQuantdeqn{B_j(u) = a_0 + a_1 u + a_2 u^2 + a_3 u^3 + a_4 u^4,}
#' the second and third derivatives are
#' \deqn{B_j''(u) = 2a_2 + 6a_3 h\,u + 12a_4 h^2 u^2,}
#' \deqn{B_j'''(u) = 6a_3 h + 24a_4 h^2 u,}
#' where \eqn{h = t_{i+1} - t_i} is the width of the interval. The function
#' returns the coefficients of these two polynomials, stored so that
#' \code{d2[i, j, ]} and \code{d3[i, j, ]} can be contracted with the spline
#' coefficient vector \code{alpha} to obtain the second and third derivative
#' of the fitted spline on interval \code{i}.
#'
#' @param BB A list describing the B-spline basis, as produced by
#'   \code{\link{Bspline_base}}. It must contain at least the elements
#'   \code{knot} (the interior + boundary knot vector) and \code{degree}
#'   (the spline degree).
#' @param x_values Numeric vector only used for evaluation points of the non-differentiated
#' Basis.
#' @return A named list with two arrays:
#' \describe{
#'   \item{\code{d2}}{Numeric array of dimension \code{c(kn, N, 4)}. Entry
#'     \code{d2[i, j, ]} gives the coefficients
#'     \code{c(a3, a2, a1, a0)} of the second derivative of basis function
#'     \eqn{j} on interval \eqn{i}, ordered from highest to lowest power of
#'     \eqn{u}: \eqn{a_3 u^3 + a_2 u^2 + a_1 u + a_0}.}
#'   \item{\code{d3}}{Numeric array of dimension \code{c(kn, N, 3)}. Entry
#'     \code{d3[i, j, ]} gives the coefficients \code{c(a2, a1, a0)} of the
#'     third derivative of basis function \eqn{j} on interval \eqn{i}, i.e.
#'     \eqn{a_2 u^2 + a_1 u + a_0}.}
#' }
#' Here \code{kn = length(knot) - 1} is the number of knot intervals and
#' \code{N = kn + degree} is the number of B-spline basis functions.
#'
#' @details
#' The derivatives are expressed in the **normalized** local variable
#' \eqn{u \in [0,1]} on each interval, with the interval width \eqn{h}
#' folded into the coefficients. This is the form expected by
#' \code{\link{apply_karlin_cubic}} and \code{\link{apply_karlin_quadratic}},
#' which impose sign or convexity constraints by testing sign changes of the
#' derivative polynomial on \eqn{[0,1]}.
#'
#' The first derivative \code{d1} is intentionally **not** returned by this
#' helper: it is not needed by the Karlin–Studden constraints implemented in
#' the current version of the package. Monotonicity constraints are handled
#' by the lower-degree constructors (\code{SplineLinearQuant},
#' \code{SplineQuadraticQuant}, \code{SplineCubicQuant},
#' \code{SplineQuarticQuant}) through their own derivative machinery.
#'
#' @seealso
#' \code{\link{Bspline_base}}, \code{\link{Bspline_base_deriv}},
#' \code{\link{apply_karlin_cubic}}, \code{\link{apply_karlin_quadratic}},
#' \code{\link{SplinePolynQuant}}
#'
#' @examples
#' ## ---- Example 1: cubic spline on [0, 1] with 3 interior knots ----
#' knot <- seq(0, 1, length.out = 5)          # 0, 0.25, 0.5, 0.75, 1
#' BB   <- Bspline_base(knot, degree = 3)
#' dc   <- bspline_to_deriv_coeffs_general(BB)
#'
#' # d2 has dimension c(kn, N, 4) = c(4, 7, 4)
#' dim(dc$d2)
#' # d3 has dimension c(kn, N, 3) = c(4, 7, 3)
#' dim(dc$d3)
#'
#' ## ---- Example 2: check the second derivative at the left knot ----
#' # B_j''(u = 0) is the constant term a0 of d2[i, j, ]
#' i <- 1; j <- 3
#' B_j_second_deriv_at_left_knot <- dc$d2[i, j, 4]
#'
#' ## ---- Example 3: quartic basis on an irregular knot grid ----
#' knot <- c(0, 0.1, 0.3, 0.6, 1)
#' BB   <- Bspline_base(knot, degree = 4)
#' dc   <- bspline_to_deriv_coeffs_general(BB)
#' stopifnot(
#'   identical(dim(dc$d2), c(4L, 8L, 4L)),
#'   identical(dim(dc$d3), c(4L, 8L, 3L))
#' )
#'
#' ## ---- Example 4: sign of the third derivative on interval 2 ----
#' # For a linear third derivative, the sign on the whole interval is
#' # determined by the two endpoint values of d3[i, j, ]:
#' i <- 2; j <- 1
#' a2 <- dc$d3[i, j, 1]; a1 <- dc$d3[i, j, 2]; a0 <- dc$d3[i, j, 3]
#' value_at_u0 <- a0
#' value_at_u1 <- a2 + a1 + a0
#' c(left = value_at_u0, right = value_at_u1)
#'
#' @export
bspline_to_deriv_coeffs_general <- function(BB,
                                            x_values = 0)
{
  knot   <- BB$knot
  degree <- BB$degree
  kn     <- length(knot) - 1
  N      <- kn + degree

  # Extended knot vector
  sn <- c(rep(knot[1], degree), knot, rep(knot[kn + 1], degree))

  # Build B-spline basis (first derivative basis)
  D1BB  <- Bspline_base_deriv(BB, der = 1)
  basis <- D1BB$base

  # Second derivative: cubic on each interval -> [a3, a2, a1, a0]
  deriv2_coeffs <- array(0, dim = c(kn, N, 4))

  # Third derivative: quadratic on each interval -> [a2, a1, a0]
  deriv3_coeffs <- array(0, dim = c(kn, N, 3))

  for (j in 1:N) {
    for (nu in (degree + 1):(kn + degree)) {
      h <- sn[nu + 1] - sn[nu]

      a0 <- basis[j, nu, 5]
      a1 <- basis[j, nu, 4]
      a2 <- basis[j, nu, 3]
      a3 <- basis[j, nu, 2]
      a4 <- basis[j, nu, 1]

      deriv2_coeffs[nu - degree, j, ] <- c(4 * a4 * h^3, 3 * a3 * h^2, 2 * a2 * h, a1)
      deriv3_coeffs[nu - degree, j, ] <- c(12 * a4 * h^2, 6 * a3 * h, 2 * a2)
    }
  }

  # Design matrix if x_values provided
  if (length(x_values) != 1) {
    yvalues <- bs_direct(BB, x_values)
  } else {
    yvalues <- 0
  }

  return(list(d0 = yvalues, d1 = deriv1_coeffs))
  return(list(yvalues,
    d2 = deriv2_coeffs,
    d3 = deriv3_coeffs
  ))
}
