# bs_direct, spline_eval, view_spline


#' Evaluate a B-spline
#'
#' Evaluates a spline (linear combination of B-splines) at given points.
#'
#' @param Bspline Spline object (list with coefficients on the 'Bspline' basis, degree, extended knot)
#' A 'Bspline' can also be rendered callable witk 'make_spline(BSpline)'
#' @param x_values Vector of evaluation points. By default, evaluation is calculated
#' at the knots
#' @param der integer 0 (defalut)...degree is the order of derivative
#' @param Bvalues : the values of the 'Bspline' basis evaluated at the values x
#' this increases the speed by avoiding re-doing the same calculations
#' of the basis for each spline with the same knots.
#' @param verbose boollean, if TRUE output more messages
#' @return vector of same length as x_values, with Spline values at the requested points
#' @examples
#'{ # Create and evaluate a spline
#' # This function is also used when a 'Bspline' object is rendered callable as a function
#' # with make_spline()
#' tn<-c(0,1,2,3,4,5)
#' x_values<-seq(0,5,length=100)
#' Bspline=list(degree=3, knot=tn,coeffs=runif(length(tn)+3-1))
#' y <- spline_eval(Bspline, x_values)
#' #alternatively compute the base before to optimize if needed
#' sn <- c(0,0,0,0,1,2,3,4,5,5,5,5)
#' Bsbase <- Bspline_base(sn, degree=3)
#' Bvalues <-bs_direct(Bsbase,x_values)
#' y <- spline_eval(Bspline=Bspline,x_values=x_values, Bvalues=Bvalues )}
#' @export

spline_eval <- function(Bspline,
                        x_values = NULL,
                        der = 0,
                        Bvalues = NULL,
                        verbose = FALSE)
{
  if (is.function(Bspline) && inherits(Bspline, "callable_spline")) {
    # Already callable, retrieve the parameters
    param <- get_parameters(Bspline)
    if(verbose){print(param)}
    Bspline<-param
  }

  knot = Bspline$knot #vector of effective knots
  degree = Bspline$degree
  coeff = Bspline$coeff
  if (verbose) {
    cat("Evaluation de la spline:\n", print(Bspline))
  }
  if (is.null(x_values)) {
    message("no x_value given, compute at the knots")
    x_values <- knot
  }
  #Bvalues=bs(x_values,knot=knot,degree)  "can be used instead of the following lines
  #to accelerate the calculations
  # values at knots by default
  if (is.null(Bvalues))
    # if the values of the basis
    #are not provided, or do not match the spline coefficients
    #then compute them
  {
    t1 = knot[1]
    tkn = rev(knot)[1] #last knot
    sn = c(rep(t1, degree), knot, rep(tkn, degree)) #extended knot partition
    D_BB = Bspline_base(sn,
                        degree = degree,
                        der = der,
                        verbose = verbose) # first compute the Basis
    #BB=Bspline_base(sn,degree=degree,verbose=verbose)
    #then evaluate the basis as the values
    #D_BB<-Bspline_base_deriv(BB,der=der,verbose=verbose)
    D_Bvalues <- bs_direct(D_BB, x_values = x_values, verbose = verbose)
  } else{
    if (inherits(Bvalues, "bspline_basis")) {
      #case basis
      BB <- Bvalues
      #We trust the dimensions, knots, are compatible
      if (der == 0) {
        D_Bvalues = bs_direct(BB, x_values)
      } else{
        D_BB <- Bspline_base_deriv(BB, der)
        D_Bvalues <- D_BB$base
        D_Bvalues <- bs_direct(D_BB, x_values)
      }
    }
    else{
      #case values
      D_Bvalues <- Bvalues
    }

    n1 = dim(D_Bvalues)[1]
    if (n1 != length(coeff)) {
      D_Bvalues <- t(D_Bvalues)
    }
  }
  yvalues <- (t(D_Bvalues)) %*% coeff

  return(yvalues)
}

#' Direct evaluation of a B-spline basis
#'
#' Computes the values of all B-spline basis functions at given points.
#'
#' @param Basis Object returned by \code{Bspline_base}
#' @param x_values Vector of evaluation points
#' @param verbose set to TRUE increases verbosity
#' @return Matrix of basis function values (n_splines x length(x_values))
#' @export

bs_direct <- function(Basis,
                      x_values = NULL,
                      verbose = FALSE)
{
  #Calcule les valeurs d'une base 'Basis'
  #comme la fonction 'bs' de 'R', mais en utilisant la
  # base calculee sous PP-forme : coeff des polynomes sur la base locale.
  #indep. du choix de la notation croissant/decroissant
  if (is.null(x_values))
  {
    message("no xvalue given, computing at knots")
    return(Spline_der_knot(Basis))
  }

  n_values = length(x_values)
  knot = Basis$knot  # knot
  kn = length(knot) - 1
  d = Basis$degree
  diff <- Basis$deriv_order
  nsplines = Basis$n_splines
  base = Basis$base
  if (verbose) {
    message(
      c(
        "Knots:  ",
        paste(knot, collapse = " , ") ,
        "\n degree of the basis: ",
        d,
        "\n Nb splines in the  basis : ",
        nsplines
      )
    )
  }
  bb=array(0,dim=c(nsplines,nsplines-d-diff,d+1))
  if (d >= 0) {
    yvalues = array(data = 0, c(nsplines, n_values))
    bb[,,] = base[, (d + 1 + diff):(nsplines), ] #only keep the effective pieces
    for (j in 1:nsplines)
      #go through splines of the base
    {
      if (kn == 1) {
        # In case of a single piece
        if (verbose) {
          print("only one piece")
        }
      #  p = makpp(bb[j, ], tn = knot)
        if (verbose) {
          print("function number", j, "is", p)
        }

       # yvalues[j, ] <- evalpp(p, x_values)
      }
#      else {

        p = makpp(bb[j, , ], tn = knot)

        yvalues[j, ] <- evalpp(p, x_values)
 #     }
    }
  }


  return(yvalues)
}


#' Convert a B-spline to Piecewise Polynomial (PP) form
#'
#' Transforms a B-spline object (callable or non-callable) into a piecewise
#' polynomial representation. The resulting PP object contains the polynomial
#' coefficients for each interval between knots.
#'
#' @param Bspline A B-spline object (list, non_callable_spline, or callable_spline)
#'        containing at least 'coeff', 'degree', and 'knot'.
#' @param Bsbasis Optional pre-computed B-spline basis object from Bspline_base().
#'        If NULL, the basis is computed automatically.
#' @param callable Logical; if TRUE, returns a callable function for evaluation.
#'        If FALSE (default), returns a non_callable_pp object. If a callable spline is given
#'        as input, then by default is return a callable pp.
#' @param verbose Logical; if TRUE, print progress messages.
#'
#' @return A PP object (piecewise polynomial) with class:
#'         - "callable_pp" if callable = TRUE
#'         - "non_callable_pp" if callable = FALSE
#'         The object can be evaluated with evalpp() or directly if callable.
#'
#' @examples
#' \dontrun{
#'
#' # Create a B-spline
#' sn <- c(0,0,0,0,1,2,3,4,5,5,5,5)
#' basis <- Bspline_base(sn, degree = 3)
#' basis$coeff <- runif(basis$n_splines)
#'
#' # Convert to PP (non-callable)
#' pp <- Bsplinetopp(basis, callable = FALSE)
#' y <- evalpp(pp, seq(0, 5, length.out = 100))
#'
#' # Convert to PP (callable)
#' pp_call <- Bsplinetopp(basis, callable = TRUE)
#' y <- pp_call(seq(0, 5, length.out = 100))
#' }
#'
#' @seealso \code{\link{makpp}}, \code{\link{evalpp}}, \code{\link{Bspline_base}}
#' @export

Bsplinetopp <- function(Bspline,
                        Bsbasis = NULL,
                        callable = FALSE,
                        verbose = FALSE) {
  #Convert a B-spline to a PP-polynomial
  if (inherits(Bspline, "callable_spline")) {
    Bspline = get_parameters(Bspline)
  }
  coeff = Bspline$coeff
  degree = Bspline$degree
  tn <- Bspline$knot
  sn <- c(rep(tn[1], degree), tn, rep(rev(tn)[1], degree)) #extended knots
  if (!is.null(Bspline$base)) {
    BB <- Bspline$base
  }
  if (is.null(Bsbasis) && is.null(Bspline$base)) {
    Bsbasis <- Bspline_base(sn, degree = degree, verbose = verbose)
    BB <- Bsbasis$base
  }
  PP <- array(0 , dim = c(length(tn) - 1, degree + 1))
  for (nu in 1:(length(tn) - 1)) {
    PP[nu, ] <- t(BB[, nu + degree , ]) %*% coeff
  }
  PP <- makpp(PP, tn, callable = callable)
}





#' Visualize a B-spline functions basis
#' Visualize a B-spline basis
#' Plots all basis functions of a B-spline basis with improved styling.
#' @param BB Object returned by \code{Bspline_base}
#' @param x_values Vector of evaluation points for plotting.
#'        If length is 1 (default = 0), generates 200 points in the knot range.
#' @param view_knot Logical; if TRUE, adds vertical lines at knot positions.
#' @param add_knot logical, same value as 'view_knot' for retro compatibility
#' @param main Title to print over the plot.
#' @return No return value, called for side effects (generates a plot)
#' @export
view_basis <- function(BB,
                       x_values = 0,
                       main = NULL,
                       view_knot = TRUE,
                       add_knot = NULL
                       ) {
  add_knot <- view_knot # for compatibility with 2.3 version
  if (length(x_values) == 1) {
    k <- range(BB$knot)
    margin <- 0.01 * diff(k)
    x_values <- seq(k[1] - margin, k[2] + margin, length.out = 200)
  }

  yvalues <- bs_direct(BB, x_values)

  n_splines <- nrow(yvalues)

  # Palette de couleurs automatique
  #colors <- rainbow(n_splines)
  if (is.null(main))
    main <- paste("B-spline Basis (degree", BB$degree, ")")
  # Graphique
  matplot(
    x_values,
    t(yvalues),
    type = "l",
    lwd = 1.5,
    xlab = "x",
    ylab = "Basis values",
    main = main,
    lty = 1
  )
  #col = colors, lty = 1)

  # Ajouter les nœuds
  # Ajouter les nœuds
  if (view_knot) {
    abline(
      v = BB$knot,
      col = "red",
      lty = 2,
      lwd = 0.8
    )
    grid(col = "gray90", lty = 1)
  }

  # Légende
  #legend("topright", legend = c("Basis", "Knots"),
  #       col = c("black", "red"), lty = c(1, 2), cex = 0.8)
}


