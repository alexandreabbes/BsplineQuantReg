# polymul, polyadd, polyderiv, poly_eval, reduce_pol
# change_polynomial_base_taylor,
# makpp, evalpp, print.callable_pp, print.non_callable_pp


#' Polynomial multiplication
#'
#' Multiplies two polynomials represented by coefficients in decreasing power order.
#'
#' @param p1 First polynomial (coefficient vector, decreasing powers)
#' @param p2 Second polynomial (coefficient vector, decreasing powers)
#' @param ord Unused (compatibility parameter)
#' @param verbose boolean FALSE (default) or TRUE.
#' @return Coefficient vector of the product polynomial (decreasing powers)
#' @examples
#' # (1 + x) * (1 + x) = 1 + 2x + x^2
#' polymul(c(1, 1), c(1, 1)) # returns c(1, 2, 1)
#' @export

polymul <- function(p1,
                    p2,
                    ord = 0,
                    verbose = FALSE)
{
  if ((sum(abs(p1)) == 0) || (sum(abs(p2)) == 0)) {
    return(0)
  }
  else
  {
    p1 = rev(p1)#invert for decreasing convention
    p2 = rev(p2)
    deg1 <- length(p1) - 1
    deg2 <- length(p2) - 1
    res <- numeric(deg1 + deg2 + 1)
    for (i in 0:deg1) {
      for (j in 0:deg2)
      {
        res[i + j + 1] <- res[i + j + 1] + p1[i + 1] * p2[j + 1]
      }
    }
    #reduce pol
    res = rev(res)#in standard notation
    res = reduce_pol(res, verbose)
    return(res)
  }
}


#' Polynomial addition
#'
#' Adds two polynomials represented by coefficients in decreasing power order.
#'
#' @param p1 First polynomial (coefficient vector)
#' @param p2 Second polynomial (coefficient vector)
#' @param verbose boolean FALSE (default) or TRUE.
#' @return Coefficient vector of the sum
#' @examples
#' polyadd(c(1, 1), c(1, -1)) # returns c(2, 0)
#' @export

polyadd <- function(p1, p2, verbose = FALSE) {
  p1 = rev(p1)
  p2 = rev(p2)
  l1 <- length(p1)
  l2 <- length(p2)
  if (l1 > l2) {
    p = p2
    p2 = p1
    p1 = p
    l = l2
    l2 = l1
    l1 = l
  }
  p1 = c(p1, rep(0, (l2 - l1)))
  sp1p2 = p1 + p2
  res = rev(p1 + p2)
  res = reduce_pol(res, verbose)
  return(res)
}


#' Change polynomial basis (Taylor expansion)
#'
#' Converts a polynomial expressed in the basis (t-a)^k to its representation
#' in the basis (t-b)^k using Taylor's formula.
#' P(t) = sum c_k (t-a)^k
#' P(t) = sum c"_k (t-b)^k
#" c'_j = P^{(j)}(b)/j!
#' @param coeffs_a Coefficients in basis centered at a (decreasing powers)
#' @param a Original expansion point
#' @param b New expansion point
#' @return Coefficients in basis centered at b
#' @export
change_polynomial_base_taylor <- function(coeffs_a, a, b)
{
  #pol is in std decreasing notation
  coeffs_a = rev(coeffs_a)
  n <- length(coeffs_a) - 1


  coeffs_b <- numeric(n + 1)

  for (j in 0:n) {
    # Calculer P^{(j)}(b) = sum_{k=j}^{n} c_k * k!/(k-j)! * (b-a)^{k-j}
    deriv_val <- 0
    for (k in j:n) {
      if (abs(coeffs_a[k + 1]) > 1e-10) {
        deriv_val <- deriv_val + coeffs_a[k + 1] *
          factorial(k) / factorial(k - j) *
          (b - a)^(k - j)
      }
    }
    coeffs_b[j + 1] <- deriv_val / factorial(j)
  }

  return(rev(coeffs_b)) #return in std Not
}

#' Reduce polynomial
#'
#' Removes leading zeros from a polynomial coefficient vector.
#'
#' @param p Polynomial coefficient vector(coefficients are in decreasing order)
#' @param verbose Boolean FALSE (default) or TRUE.
#' @return Reduced vector (without leading zeros)
#' @examples
#' reduce_pol(c(0,0, 1, 1))
#' @export
reduce_pol <- function(p, verbose = FALSE) {
  l = length(p)
  k = 1
  while (p[k] == 0 & k < l) {
    k = k + 1
  }
  if (verbose)
  {
    message("removed ", k, " useless zeroes to ", p)
  }
  return(p[k:l])
}
#' Evaluate polynomial
#'
#' Evaluates a polynomial at one or more points.
#'
#' @param p Coefficient vector (decreasing powers)
#' @param xvalues vector at which to evaluate the polynomial
#' @return Vector of same length as xvalues : polynomial values at the points xvalues
#' @examples
#' # P(x) = 1 + x + x^2
#' poly_eval(c(1, 1, 1), c(0, 1, 2)) # returns c(1, 3, 7)
#' @export

poly_eval <- function(p, xvalues) {
  p = rev(p)
  #we evaluate the values in the convention p=c(p0,p1,p2)
  #represent the polynomial p0+p1x+p2*x^2
  y = c()
  d = length(p) - 1
  for (x in xvalues) {
    val = 0
    for (i in 0:d) {
      val = val + p[i + 1] * x^i
    }
    y = c(y, val)
  }
  return(y)
}

#' Polynomial derivative
#' Computes the derivative of order \code{der} of a polynomial.
#' @param p Coefficient vector (decreasing powers)
#' @param der Derivative order (default = 1)
#' @return Coefficients of the derivative polynomial
#' @examples
#' # P(x) = x^2 -> P'(x) = 2x
#' polyderiv(c(1, 0, 0), 1) # returns c(2, 0)
#' @export

polyderiv <- function(p, der = 1) {
  if (der == 0)   {
    q = p
  }
  else {
    if (length(p) == 1) {
      q = 0
    } else{
      l = length(p)
      #p=rev(p)#reverse
      A = array(data = 0, c(l, l))
      for (i in 1:(l - 1)) {
        A[i + 1, (i)] = (l - i)
      }
      D = A
      if (der > 1) {
        for (i in 2:der) {
          D = A %*% D
        }
      } #compute the d-th power of A
      q = D %*% p
      q = q[(1 + der):(l)]
    }

    return(q)
  }
}


#' Evaluates a piecewise polynomial (PP form)  function at given points.
#'
#' @param p List with components \code{ext_knot} (ext_knot) and \code{coeff}
#' @param x_values Vector of evaluation points
#' @return Function values at the requested points
#' @export
#' @keywords internal

evalpp <- function(p, x_values) {
  #this evaluates a polynomial p under the pp form,
  #p if given with its knot and the local coefficients
  #This funciton is independent from the order convention
  #for polynomials
  #The x_values out of the knot give 0 in the corresponding yvalues
  if (inherits(p,'callable_pp')){p<-get_parameters(p)}
  tn = p$knot
  coeff = p$coeff
  kn = length(tn) - 1 #number of intervals
  n_values = length(x_values)
  degree<-p$degree
  pval <- c()
  xvalues<-as.vector(x_values)
  if (kn == 1) {
    #Only one piece
    pval <- poly_eval(coeff, x_values - tn[1])
  } else if (kn>1)
  {
    if (x_values[1] < tn[1]) {
      message("Some x values smaler than first knot, extrapolating")
      #values before the first knot
      pre_k = x_values[(x_values < tn[1])]
      h = pre_k - tn[1]
      if (degree>0){
        poly_loc <- coeff[1, ] # extrapolate using the first piece
        pval <- poly_eval(poly_loc, h)
      } else if (degree==0){
        pval=rep(coeff[1],length(h))
      }}

    for (i in 1:(kn))
    {
      xval = x_values[(x_values >= tn[i]) & (x_values < tn[i + 1])]
      if (degree>0)
      {poly_loc <- coeff[i, ]} else{
        poly_loc <- coeff[i]
      }

      # reverse our convention to match polyval convention
      #pval<-c(pval,polyval(p=rev(poly_loc),xval) )
      h = xval - tn[i]
      #shit to fit the local basis
      yval <- poly_eval(poly_loc, h)
      pval <- c(pval, yval)
    }


    xval = x_values[x_values == tn[kn + 1]]
    if (length(xval) > 0) {
      h = tn[kn + 1] - tn[kn] # if the last knot is in x_values
      yval<-poly_eval(poly_loc, h)
      pval = c(pval,yval )

    }
    if (x_values[n_values] > tn[kn + 1]) {
      message("Some x values greater than last knot, extrapolating")
      # values after the last knot
      post_k = x_values[(x_values > tn[kn + 1])]
      h = post_k - tn[kn]
      if(degree>0){
        poly_loc <- coeff[kn, ] # extrapolate using the last piece
        yval<- poly_eval(poly_loc, h)
        pval <- c(pval, yval)
      }
      else if (degree==0){

        pval<-c(pval,rep(coeff[kn],length(h) ))}
    }
  }

  return(pval)
}


#' Build a piecewise polynomial (PP) form
#'
#' Creates a PP structure from polynomial coefficients and knots.
#' If callable = TRUE, returns a function that evaluates the PP.
#'
#' @param coeff coefficients matrix or pp polynomial. Coefficient matrix (kn x (degree+1))
#' @param tn a Knot vector of length kn+1. Needed if coefficient is not pp
#' @param callable Boolean; if TRUE, returns a callable function
#' @param verbose Boolean
#' @return A PP object or a callable function
#' @export
makpp <- function(coeff,
                  tn=NULL,
                  callable = FALSE,
                  verbose = FALSE) {
  if ( inherits(coeff,'callable_pp')){
     if (callable) {
      if(verbose){
        cat("Already a 'callable_pp'")}

        return(coeff) }
    else{ pp<-get_parameters(coeff)
    coeff<-pp$coeff
    tn<-pp$knot
    if (verbose)
    cat("transform 'non_callable_pp' to 'callable_pp'\n")}}

    if ( inherits(coeff,'non_callable_pp')){
      { if (!callable) {
        if(verbose)
          cat("Already a 'non_callable_pp'\n ")

        return(coeff) }
        else{ pp<-coeff
        coeff<-pp$coeff
        tn<-pp$knot
        if (verbose)
          cat("transform 'non_callable_pp' to 'callable_pp'\n")
        }}
    }
#  coeff<-as.array(coeff)


  if (length(tn) == 2) {
    kn <-1 #only one intervals
    degree<-length(coeff)-1 # only one polynopial
  } else if (!is.null(dim(coeff))) {
    kn <- dim(coeff)[1]
    degree <- dim(coeff)[2] - 1
  } else if (any(dim(coeff)==1) || is.null(dim(coeff))){
    degree = 0
    kn <- length(coeff)
  }

  if (length(tn) != (kn + 1)) {
    stop("length of coefficients and number of knots do not match")
  }

  # Creer l'objet PP
  pp_obj <- list(coeff = coeff,
                 knot = tn,
                 degree = degree)
  class(pp_obj) <- "non_callable_pp"

  if (callable) {
    # Retourner une fonction callable
    eval_func <- function(x_values) {
      evalpp(pp_obj, x_values)
    }

    # Ajouter les parametres comme attributs
    attr(eval_func, "coeff") <- coeff
    attr(eval_func, "knot") <- tn
    attr(eval_func, "degree") <- degree
    attr(eval_func, "pp_obj") <- pp_obj

    class(eval_func) <- c("callable_pp", "function")
    return(eval_func)
  }

  return(pp_obj)
}


#' Print method for callable_pp objects
#'
#' @param x A callable_pp object
#' @param ... Additional arguments
#' @export
print.callable_pp <- function(x, ...) {
  # Get the attributes
  degree <- attr(x, "degree")
  knot <- attr(x, "knot")
  coeff <- attr(x, "coeff")
  pp_obj <- attr(x, "pp_obj")

  # Si pas d'attributs, essayer de les extraire de l'environnement
  if (is.null(degree) || is.null(knot) || is.null(coeff)) {
    env <- environment(x)
    degree <- env$degree %||% attr(x, "degree")
    knot <- env$knot %||% attr(x, "knot")
    coeff <- env$coeff %||% attr(x, "coeff")
  }

  cat("Callable Piecewise Polynomial (PP) Object\n")
  cat("==========================================\n")
  cat("  Degree:", degree %||% "unknown", "\n")
  cat("  Intervals:", if (!is.null(knot))
    length(knot) - 1
    else
      "unknown", "\n")
  cat("  Knots:", if (!is.null(knot))
    length(knot)
    else
      "unknown", "\n")
  if (!is.null(knot)) {
    cat("  Knot range: [", round(min(knot), 4), ", ", round(max(knot), 4), "]\n")
  }
  if (!is.null(coeff)) {
    cat("  coefficients:",
        if (is.matrix(coeff))
          paste(dim(coeff), collapse = " x ")
        else
          length(coeff),
        "\n")
  }
  cat("\n  Usage: pp(x_values) or evalpp(pp, x_values)\n")
  invisible(x)
}

#' Print method for non_callable_pp objects
#'
#' @param x A non_callable_pp object
#' @param ... Additional arguments
#' @export
print.non_callable_pp <- function(x, ...) {
  coeff <- x$coeff
  knot <- x$knot
  degree <- x$degree
  n_intervals <- length(knot) - 1
  cat("Piecewise Polynomial (PP) (non-callable)\n")
  cat("================================\n")
  cat("  $degree:", degree, "\n")
  cat("  $knots:", length(knot), knot, "\n")
  if (!is.null(dim(coeff))) {
    cat(" coefficients dimension:",
        dim(coeff)[1],
        "x",
        dim(coeff)[2],
        "\n")
    # Afficher les coefficients (troncated if too numerous)
    if (n_intervals <= 5 && dim(coeff)[2] <= 4) {
      cat("\n  $coeff:\n")
      for (i in 1:n_intervals) {
        cat("    Interval", i, ":", paste(round(coeff[i, ], 4), collapse = ", "), "\n")
      }
    } else {
      cat("\n  First interval coefficients:",
          paste(round(coeff[1, ], 4), collapse = ", "),
          "\n")
      if (n_intervals > 1) {
        cat("  Last interval coefficients: ",
            paste(round(coeff[n_intervals, ], 4), collapse = ", "),
            "\n")
      }
    }
  } else {
    cat("  Coefficients length:", length(coeff), "\n")
  }
  cat("\n  Usage: pp_eval(pp, x_values)\n")
  invisible(x)
}
