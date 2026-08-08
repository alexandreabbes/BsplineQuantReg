#Bspline_deriv, Spline_der_knot

#' Compute derivative coefficients of a B-spline
#'
#' Given a B-spline of given degree, coefficients, knots, compute the coefficients
#' of the derivative of order 'der' in the B-spline basis of degree d-der.
#' Uses the standard algorithm for differentiating bsplines
#'
#'
#' @param bspline callable_spline or non_callable_spline class list
#'  or function, or list with variable names 'degree', 'knot' ,'coeff'.
#' @param der Derivative order (default = 1)
#' @param callable Boolean, if TRUE returns a callable Bspline. If not provided
#' @param verbose Boolean gives output messages when TRUE
#' return the same class as input : 'callable' or 'non_callable'
#' @return Coefficients of the derivative B-spline (degree = degree - der)
#' @export

Bspline_deriv <- function(bspline,
                          der = 1,
                          callable = NULL,
                          verbose = FALSE) {
  if (der == 0) {
    diff_spline <- make_spline(bspline)
    if (verbose)
      cat("No change, nul derivative")
  } else{
    if (inherits(bspline, "callable_spline")) {
      params <- get_parameters(bspline)
      coeff <- params$coeff
      knot <- params$knot
      degree <- params$degree
      result <- params$result
      if (is.null(callable)) {
        callable <- TRUE
      }
    } else {
      coeff <- bspline$coeff
      knot <- bspline$knot
      degree <- bspline$degree
      result <- bspline$result
      if (is.null(callable)) {
        callable <- FALSE
      }
    }
    if (der > degree) {
      if (verbose) {
        message("Derivative order exceeds degree")
      }
      coeff <- rep(0, length(knot) - 1)
      d <- 0
    }

    # Start with original coefficients

    co <- coeff

    for (j in 1:der) {
      # Number of coefficients decreases by 1 at each derivative

      d <- degree - j + 1
      N <- length(co)

      co_new <- numeric(N - 1)
      s <- c(rep(knot[1], d), knot, rep(rev(knot)[1], d))
      # Extended knots for the current degree
      # t_i for i = 1..n+d+1 (standard knot sequence)
      # For a B-spline of degree d with n basis functions,
      # the knot sequence has length n + d + 1
      # Here we need t_{i} and t_{i+d} for i = 1..n-1

      # For the derivative formula, we need t_i and t_{i+d}
      # where i goes from 1 to n-1 (new basis functions)
      for (i in (1):(N - 1)) {
        # s_i is the i-th knot (1-indexed)
        # t_{i+d} is the (i+d)-th knot
        # We need the effective knots
        denom <- (if (is.na(s[i + d + 1])) {
          0
        } else{
          s[i + d + 1]
        }) - (if (is.na(s[i + 1])) {
          0
        } else{
          s[i + 1]
        })
        #usual formula      #s[i+k-j]-s[i] k= initial deg

        if (denom == 0) {
          co_new[i] <- 0
        } else {
          co_new[i] <- d * (co[i + 1] - co[i]) / denom
        }
      }

      # Update for next derivative
      co <- co_new

    }

    diff_list = list(
      coeff = co,
      degree = degree - der,
      knot = knot,
      result = result
    )
    diff_spline <- make_spline(diff_list, callable = callable)
  }
  return(diff_spline)
}

#' Derivatives at knot of a B-spline
#'
#' Computes derivative values of a B-spline at knot (efficient because it
#' directly uses polynomial coefficients).
#'
#' @param Bsbase Object returned by \code{Bspline_base}
#' @param der Derivative order (default = 1)
#' @return Matrix of derivative values (n_splines x n_knot)
#' @export

Spline_der_knot <- function(Bsbase, der = 1)
  #compute the values of a derivatives only at the knot
  #(simple, it only uses the coefficients)
{
  coeff = Bsbase$base
  nsplines = Bsbase$n_splines
  tn = Bsbase$knot
  kn = length(tn) - 1
  m = Bsbase$degree

  if (der > m) {
    Der2_knot = array(data = 0, c(nsplines, kn))
  }
  else{
    #Der2_knot=coeff[,,(der+1)]*factorial(der) #in increasing pol notation
    Der2_knot = coeff[, , (m - der + 1)] * factorial(der) #in decreasing notation
    #computation of the last value
    h = tn[kn] - tn[kn - 1]
    for (j in 1:nsplines)
    {
      p_kn_der = polyderiv(coeff[j, kn + m - 1, ], der)

      v_kn = poly_eval(p_kn_der, h)
      Der2_knot[j, kn + m] = v_kn
    }
  }
  return(t(Der2_knot))
}

