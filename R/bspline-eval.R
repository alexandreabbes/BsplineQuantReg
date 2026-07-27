# bs_direct, makpp, evalpp, spline_eval, view_spline
# Spline_der_knot, bspline_to_deriv_coeffs_pp
# make_spline, print.callable_spline

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
#' Bspline=list(degree=3, knot=tn,coefficients=runif(length(tn)+3-1))
#' y <- spline_eval(Bspline, x_values)
#' #alternatively compute the base before to optimize if needed
#' sn <- c(0,0,0,0,1,2,3,4,5,5,5,5)
#' Bsbase <- Bspline_base(sn, degree=3)
#' Bvalues <-bs_direct(Bsbase,x_values)
#' y <- spline_eval(Bspline=Bspline,x_values=x_values, Bvalues=Bvalues )}
#' @export

spline_eval<-function(Bspline, x_values=NULL,der=0, Bvalues=NULL,verbose=FALSE)
{
  if (is.function(Bspline) && inherits(Bspline, "callable_spline") ) {
    # Already callable, retrieve the parameters
    Bspline<-get_parameters(Bspline)
  }


  knot=Bspline$knot #vector of effective knots
  degree=Bspline$degree
  coeff=Bspline$coeff
  if (is.null(x_values)){message("no x_value given, compute at the knots")
    x_values<-knot}
  #Bvalues=bs(x_values,knot=knot,degree)  "can be used instead of the following lines
  #to accelerate the calculations
  # values at knots by default
  if (is.null(Bvalues))
      # if the values of the basis
      #are not provided, or do not match the spline coefficients
      #then compute them
      {t1=knot[1]
      tkn=rev(knot)[1] #last knot
      sn=c(rep(t1,degree),knot,rep(tkn,degree) ) #extended knot partition
      BB=Bspline_base(sn,degree=degree) # first compute the Basis
      # then evaluate the basis as the values
      D_BB<-Bspline_base_deriv(BB,der=der,verbose=verbose)
      D_Bvalues<-bs_direct(D_BB,x_values=x_values,verbose=verbose)
      }else{
    if (inherits(Bvalues,"bspline_basis")){ #case basis
    BB<-Bvalues
  #We trust the dimensions, knots, are compatible
    if (der==0){
      D_Bvalues=bs_direct(BB,x_values)
    }else{
      D_BB<-Bspline_base_deriv(BB,der)
    D_Bvalues<-D_BB$base
    D_Bvalues<-bs_direct(D_BB,x_values)}
    }
    else{ #case values
      D_Bvalues<-Bvalues}

    n1=dim(D_Bvalues)[1]
    if (n1!=length(coeff)){D_Bvalues<-t(D_Bvalues)}
  }
   yvalues<-(t(D_Bvalues))%*%coeff

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

bs_direct<-function(Basis,x_values=NULL,verbose=FALSE)
{
  #Calcule les valeurs d'une base 'Basis'
  #comme la fonction 'bs' de 'R', mais en utilisant la
  # base calculee sous PP-forme : coeff des polynomes sur la base locale.
  #indep. du choix de la notation croissant/decroissant
  if (is.null(x_values))
    {message("no xvalue given, computing at knots")
    return(Spline_der_knot(Basis))
           }

  n_values=length(x_values)
  knot=Basis$knot  # knot
  kn=length(knot)-1
  d=Basis$degree
  diff<-Basis$deriv_order
  nsplines=Basis$n_splines
  if(verbose){message(c("Noeuds:  ", paste(knot, collapse=" , ") ,
                        "\n degré de la base: ",d,
                        "\n Nb splines dans la base : ",nsplines
  ))}

  yvalues=array(data=0,c(nsplines,n_values))
    if (d>0){
      bb=Basis$base[,(d+1+diff):(nsplines),] #only keep the effective pieces
    for (j in 1:nsplines) #go through splines of the base
      {
      if (kn==1){ # In case of a single piece
        if (verbose){print("only one piece")}
        p=makpp(bb[j,],tn=knot)
        if (verbose){print("function number",j,"is",p)}
        yvalues[j,]<-evalpp(p,x_values)}
      else {
        p=makpp(bb[j,,],tn=knot)
        yvalues[j,]<-evalpp(p,x_values)
        }
      }
      }
      if (d==0){
        bb=Basis$base
        for (j in 1:nsplines)
          {
          p=makpp(bb[j,],tn=knot)
          yvalues[j,]<-evalpp(p,x_values)
          }
      }

   return(yvalues)
  }


#' Evaluates a piecewise polynomial (PP form)  function at given points.
#'
#' @param p List with components \code{ext_knot} (ext_knot) and \code{coeff}
#' @param x_values Vector of evaluation points
#' @return Function values at the requested points
#' @export
#' @keywords internal

evalpp<-function(p,x_values){
  #this evaluates a polynomial p under the pp form,
  #p if given with its knot and the local coefficients
  #This funciton is independent from the order convention
  #for polynomials
  #The x_values out of the knot give 0 in the corresponding yvalues

  tn=p$knot
  coeff=p$coefficients
  kn=length(tn)-1 #number of intervals
  n_values=length(x_values)

  pval<-c()
  if(kn==1){#Only one piece
    pval<-poly_eval(coeff,x_values-tn[1])
  }else{
  if (!is.null(dim(coeff))){#if the degree is not 0
    if (x_values[1]<tn[1]){ message("Some x values smaler than first knot, extrapolating")
      #values before the first knot
      pre_k=x_values[(x_values<tn[1])]
      poly_loc<-coeff[1,] # extrapolate using the first piece
      h=pre_k-tn[1]
      pval<-poly_eval(poly_loc,h)
      }

    for (i in 1:(kn))
    {
    xval=x_values[(x_values>=tn[i]) & (x_values<tn[i+1])]
    poly_loc<-coeff[i,]
    # reverse our convention to match polyval convention
    #pval<-c(pval,polyval(p=rev(poly_loc),xval) )
    h=xval-tn[i]
    #shit to fit the local basis
    yval<-poly_eval(poly_loc,h)
    pval<-c(pval,yval)
    }}
  if (is.null(dim(coeff)))# if degree=0
    {
    for (i in 1:(kn)) {
    xval=x_values[(x_values>=tn[i]) & (x_values<tn[i+1])]
    poly_loc<-coeff[i]
    # reverse our convention to match polyval convention
    #pval<-c(pval,polyval(p=rev(poly_loc),xval) )
    h=xval-tn[i]
    #shit to fit the local basis
    yval<-poly_eval(poly_loc,h)

    pval<-c(pval,yval)

}}
  xval=x_values[x_values==tn[kn+1]]
  if (length(xval)>0){
  h=tn[kn+1]-tn[kn] # if the last knot is in x_values
  pval=c(pval,poly_eval(poly_loc,h))

  }
  if (x_values[n_values]>tn[kn+1]){message("Some x values greater than last knot, extrapolating")
    # values after the last knot
    post_k=x_values[(x_values>tn[kn+1])]
    poly_loc<-coeff[kn,] # extrapolate using the last piece
    h=post_k-tn[kn]
    pval<-c(pval,poly_eval(poly_loc,h))
  }
  }

  return(pval)
}


#' Build a piecewise polynomial (PP) form
#'
#' Creates a PP structure from polynomial coefficients and knots.
#' If callable = TRUE, returns a function that evaluates the PP.
#'
#' @param coefficients Coefficient matrix (kn x (degree+1))
#' @param tn a Knot vector of length kn+1.
#' @param callable Boolean; if TRUE, returns a callable function
#' @param verbose Boolean
#' @return A PP object or a callable function
#' @export
makpp <- function(coefficients, tn, callable = FALSE, verbose=FALSE) {
  # Déterminer kn
  if (length(tn) == 2) {
    kn <- length(tn) - 1
  } else if (!is.null(dim(coefficients))) {
    kn <- dim(coefficients)[1]
    degree <-dim(coefficients)[2]-1
  } else {
    degree=0
    kn <- length(coefficients)
  }

  if (length(tn) != (kn + 1)) {
    stop("length of coefficients and number of knots do not match")
  }

  # Créer l'objet PP
  pp_obj <- list(
    coefficients = coefficients,
    knot = tn,
    degree = degree
  )
  class(pp_obj) <- "non_callable_pp"

  if (callable) {
    # Retourner une fonction callable
    eval_func <- function(x_values) {
      evalpp(pp_obj, x_values)
    }

    # Ajouter les paramètres comme attributs
    attr(eval_func, "coefficients") <- coefficients
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
  # Récupérer les attributs
  degree <- attr(x, "degree")
  knot <- attr(x, "knot")
  coeff <- attr(x, "coefficients")
  pp_obj <- attr(x, "pp_obj")

  # Si pas d'attributs, essayer de les extraire de l'environnement
  if (is.null(degree) || is.null(knot) || is.null(coeff)) {
    env <- environment(x)
    degree <- env$degree %||% attr(x, "degree")
    knot <- env$knot %||% attr(x, "knot")
    coeff <- env$coeff %||% attr(x, "coefficients")
  }

  cat("Callable Piecewise Polynomial (PP) Object\n")
  cat("==========================================\n")
  cat("  Degree:", degree %||% "unknown", "\n")
  cat("  Intervals:", if (!is.null(knot)) length(knot) - 1 else "unknown", "\n")
  cat("  Knots:", if (!is.null(knot)) length(knot) else "unknown", "\n")
  if (!is.null(knot)) {
    cat("  Knot range: [", round(min(knot), 4), ", ", round(max(knot), 4), "]\n")
  }
  if (!is.null(coeff)) {
    cat("  Coefficients:", if (is.matrix(coeff)) paste(dim(coeff), collapse=" x ") else length(coeff), "\n")
  }
  cat("\n  Usage: pp(x_values) or pp_eval(pp, x_values)\n")
  invisible(x)
}

#' Print method for non_callable_pp objects
#'
#' @param x A non_callable_pp object
#' @param ... Additional arguments
#' @export
print.non_callable_pp <- function(x, ...) {
  coeff <- x$coefficients
  knot <- x$knot
  degree <- x$degree
  n_intervals<-length(knot)-1
  cat("Piecewise Polynomial (PP) (non-callable)\n")
  cat("================================\n")
  cat("  $degree:", degree, "\n")
  cat("  $knots:", length(knot),knot, "\n" )
  if (!is.null(dim(coeff))) {
    cat(" Coefficients dimension:", dim(coeff)[1], "x", dim(coeff)[2], "\n")
    # Afficher les coefficients (tronqués si trop nombreux)
    if (n_intervals <= 5 && dim(coeff)[2] <= 4) {
      cat("\n  $coefficients:\n")
      for (i in 1:n_intervals) {
        cat("    Interval", i, ":", paste(round(coeff[i, ], 4), collapse = ", "), "\n")
      }
    } else {
      cat("\n  First interval coefficients:", paste(round(coeff[1, ], 4), collapse = ", "), "\n")
      if (n_intervals > 1) {
        cat("  Last interval coefficients: ", paste(round(coeff[n_intervals, ], 4), collapse = ", "), "\n")
      }
    }
  } else {
    cat("  Coefficients length:", length(coeff), "\n")
  }
  cat("\n  Usage: pp_eval(pp, x_values)\n")
  invisible(x)
}



#' Visualize a B-spline functions basis
#'
#' Plots all  functions of a B-spline basis.
#'
#' @param BB Object returned by \code{Bspline_base}
#' @param x_values Vector of evaluation points for plotting (by default 100 points are computed in the knot range)
#' @return No return value, called for side effects (generates a plot)
#' @export
#'
view_basis<-function(BB,x_values=0)
{
  if (length(x_values)==1){
    k=range(BB$ext_knot)
    x_values=(k[1]:(k[2]*100))/100}

  yvalues=bs_direct(BB,x_values)

  matplot(x_values, t(yvalues))
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

Spline_der_knot<-function(Bsbase,der=1)
  #compute the values of a derivatives only at the knot
  #(simple, it only uses the coefficients)
{
  coeff=Bsbase$base
  nsplines=Bsbase$n_splines
  tn=Bsbase$ext_knot
  kn=length(tn)-1
  m=Bsbase$degree

  if (der>m){
    Der2_knot=array(data=0,c(nsplines,kn))
  }
  else{
    #Der2_knot=coeff[,,(der+1)]*factorial(der) #in increasing pol notation
    Der2_knot=coeff[,,(m-der+1)]*factorial(der) #in decreasing notation
    #computation of the last value
    h=tn[kn]-tn[kn-1]
    for (j in 1:nsplines)
    {
      p_kn_der=polyderiv(coeff[j,kn+m-1,],der)

      v_kn=poly_eval(p_kn_der,h)
      Der2_knot[j,kn+m]=v_kn
    }
  }
  return(t(Der2_knot))
}






#' Create a callable spline object
#'
#' Transforms a list container 'Bspline' that may result from
#'  any of the regression functions associated to 'quantile_spline' into a
#' callable function that is the spline function with given knots and
#' coefficients on the corresponding B-spline basis. It can be evaluated
#' at any point of the knot range transfering the call to the 'spline_eval()' function,
#' while preserving access to parameters.
#' Outside the knots range, it extrapolates the side part of the spline.
#' @param Bspline A list containing at least (knot, coefficients, degree)
#'        like the one returned by quantile_spline or one of the degree-specific functions.
#'        If it already contains a 'spline' element, returns the object unchanged.
#' @param verbose Boolean. If TRUE : output the attributes of the 'Bspline'
#' @param callable is a Boolean flag that allows to return a simple
#'        spline list with class 'non_callable_spline'.
#' @return Either a callable function can be called as `Bspline(x)`, with
#'         class ('callable_spline'), and all parameters as attributes
#'         or the list of parameters with an additional 'non_callable_spline' classe.
#'
#' @export
make_spline <- function(Bspline,verbose=FALSE,callable=TRUE){
  if (!callable) {#erase the callable status
    if (inherits(Bspline, "non_callable_spline") )
      { #non callable from non callable
      message ("Already a non_callable_spline object")
      if (verbose) {print(Bspline)}
      return(Bspline)
      }
    else if(inherits(Bspline, "callable_spline")){
      NCspline_obj<-get_parameters(Bspline)
      class(NCspline_obj)<-c("non_callable_spline","list")
      if(verbose){print(NCspline_obj)}
        return(NCspline_obj)
      }
      else {#nothing to non callable
        for (i in (1:3)){if (!is.element(c("knot","degree","coeff")[i],names(Bspline) ))
          {message(c("missing ", c("knot","degree","coeff")[i]," in input"))
          return(NULL)
          }}
        NCspline_obj<-Bspline
        class(NCspline_obj)<-c("non_callable_spline","list")
        return(NCspline_obj)
        if(verbose){print(NCspline_obj)}
      }
    }

  if (callable)
    {
    #Callable from callable
        if (is.function(Bspline) && inherits(Bspline, "callable_spline"))
          {
    # Déjà une callable, l'utiliser directement
    message("Already a callable spline object")
    return(Bspline)
    if (verbose){print(Bspline)}}
    else
  {  # Callable from non_callable
    #Extract components
    for (i in (1:3)){if (!is.element(list("knot","degree","coeff")[i],names(Bspline) ))
    {message(c("missing ", c("knot","degree","coeff")[i]," in input"))
      return(NULL)
    }}
  coeff <- Bspline$coeff
  degree <- Bspline$degree
  knot <- Bspline$knot
  result <- Bspline$result # keep other parameters
  # Create the callable function
  spline_func <- function(x_values){
    # Build the spline object structure expected by spline_eval
    spline_obj <- list(
      coeff = coeff,
      degree = degree,
      knot = knot,
      result=result
    )
    # Evaluate the spline
    spline_eval(Bspline, x_values)
  }
  # Attach parameters as attributes (accessible via attr())
  attr(spline_func, "degree") <- degree
  attr(spline_func, "knot") <- knot
  attr(spline_func, "coeff") <- coeff

  # Store the full 'Bspline' as an attribute
  attr(spline_func, "result") <- result

  # Set class for print method
  class(spline_func) <- c("callable_spline", "function")
  if (verbose){print(spline_func)}
  return(spline_func)
  }
  }
  }


#' Get parameters from a callable spline or a callable pp object
#'
#' @param x A callable spline object
#' @return A list with degree, knots, coefficients, and 'Bspline' object
#' @export
get_parameters <- function(x) {
  if (!inherits(x, "callable_spline") &&
      !inherits(x,"callable_pp")) {
    message("Object is neither a callable spline
            or a callable pp")
    return(NULL)
  }
  else
  list(
    degree = attr(x, "degree"),
    knot = attr(x, "knot"),
    coeff = attr(x, "coeff"),
    result = attr(x, "result")
  )
  }



#' Print method for callable spline
#'
#' @param x A callable spline object
#' @param ... Additional arguments
#' @export
print.callable_spline <- function(x, ...) {
  cat("callable_spline Object\n")
  cat("======================\n")
  cat("  Degree:", attr(x, "degree"), "\n")
  knot <- attr(x, "knot")
  coeff <- attr(x, "coeff")
  result <- attr(x, "result")
  if (!is.null(knot)) {
    cat("  Knots (", length(knot), "): ",knot,"\n")
  } else {
    cat("  Knots: NULL\n")
  }
  cat("  Coefficients (", length(coeff),"): ", coeff, "\n")
  cat("for result, type get_parameter(spline)\n")
    invisible(x)
}

#' Print method for non_callable_spline results
#'
#' @param x Result object from quantile_spline when callable=FALSE
#' @param ... Additional arguments
#' @export
print.non_callable_spline <- function(x, ...) {
  cat("Non callable Spline List \n")
  cat("======================== \n")
  cat(" $degree: ", x$degree, "\n")
  cat(" $knot : [",paste(x$knot,collapse = ", " ) ,"]", "(", length(x$knot), "knots ) \n")
  cat(" $coeff (rounded 10^(-7)) :[",paste(round(x$coeff,7),collapse = ", " ) , "]  ( dim Basis is ",length(x$coeff),")\n")
  if (!is.null(x$result)){cat(" $result is non NULL")}
  invisible(x)
}
