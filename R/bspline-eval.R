# bs_direct, makpp, evalpp, spline_eval, view_spline
# Spline_der_knot, bspline_to_deriv_coeffs_pp

#' Evaluate a B-spline
#'
#' Evaluates a spline (linear combination of B-splines) at given points.
#'
#' @param Bspline Spline object (list with coefficients on the Bspline basis, degree, extended knot)
#' @param x_values Vector of evaluation points. By default, evaluation is calculated
#' at the knots
#' @param Bvalues : the values of the Bspline basis evaluated at the values x
#' this increases the speed by avoiding re-doing the same calculations
#' of the basis for each spline with the same knots.
#' @return vector of same length as x_values, with Spline values at the requested points
#' @examples
#'{ # Create and evaluate a spline
#' # This function is also used when a Bspline object is rendered callable as a function
#' # with meake_spline()
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

spline_eval<-function(Bspline, x_values=NULL, Bvalues=NULL)
  #Bspline has a new type R container,
  #designed by its coefficients, the degree
  #and the knot,
  #It is independent from the polynomial notation order
{
  knot=Bspline$knot #vector of effective knots
  degree=Bspline$degree
  coeff=Bspline$coefficients
  #Bvalues=bs(x_values,knot=knot,degree)  "can be used instead of the following lines
  #to accelerate the calculations
  if (is.null(x_values)){x_values<-knots}# values at knots by default
  if (is.null(Bvalues))
      # if the values of the basis
      #are not provided, or do not match the spline coefficients
      #then compute them
      {t1=knot[1]
      tkn=rev(knot)[1] #last knot
      sn=c(rep(t1,degree),knot,rep(tkn,degree) ) #extended knot partition
      BB=Bspline_base(sn,degree=degree) # first compute the Basis
      Bvalues=t(bs_direct(BB,x_values))} # then evaluate the basis as the values}
  db=dim(Bvalues)[2]
  if (!db==length(coeff)){Bvalues=t(Bvalues)}
  yvalues<-(Bvalues)%*%coeff
  return(yvalues)
}



#' Direct evaluation of a B-spline basis
#'
#' Computes the values of all B-spline basis functions at given points.
#'
#' @param Basis Object returned by \code{Bspline_base}
#' @param x_values Vector of evaluation points
#' @return Matrix of basis function values (n_splines x length(x_values))
#' @export

bs_direct<-function(Basis,x_values,verbose=FALSE)
{
  #Calcule les valeurs d'une base Bspline.
  #comme la fonction bs de R, mais en utilisant la
  # base calculee sous PP-forme : coeff des polynomes sur la base locale.
  #indep. du choix de la notation croissant/decroissant
  n_values=length(x_values)
  knot=Basis$knot  # knot
  kn=length(knot)-1
  d=Basis$degree
  nsplines=Basis$n_splines
  if(verbose){print("La base est:\n",Basis)}
  yvalues=array(data=0,c(nsplines,n_values))
    if (d>0){
      bb=Basis$base[,(d+1):(nsplines),] #only keep the effective pieces
    for (j in 1:nsplines) #go through splines of the base
    { print(Basis)
      print(c("coef",bb ,"knot",knot))
    # In case of a single piece
      if (kn==1){
        p=makpp(bb[j,],tn=knot)
        yvalues[j,]<-evalpp(p,x_values)}
        else {p=makpp(bb[j,,],tn=knot)
        yvalues[j,]<-evalpp(p,x_values)
        }}
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
#' Evaluate a piecewise polynomial (PP) form
#'
#' Evaluates a piecewise polynomial function at given points.
#'
#' @param p List with components \code{ext_knot} (ext_knot) and \code{coeff}
#' @param x_values Vector of evaluation points
#' @return Function values at the requested points
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
  if (!is.null(dim(coeff))){
    for (i in 1:(kn))
  {
    xval=x_values[(x_values>=tn[i]) & (x_values<tn[i+1])]
    poly_loc<-coeff[i,]
    # reverse our convention to match polyval convention
    #pval<-c(pval,polyval(p=rev(poly_loc),xval) )
    h=xval-tn[i]
    #shit to fit the local basis
    pval<-c(pval,poly_eval(poly_loc,h))
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
    pval<-c(pval,poly_eval(poly_loc,h))
}}
  xval=x_values[x_values==tn[kn+1]]
  if (length(xval)>0){
  h=tn[kn+1]-tn[kn] # if the last knot is in x_values
  pval=c(pval,poly_eval(poly_loc,h))

  }}
  return(pval)
}

#' Build a piecewise polynomial (PP) form
#'
#' Creates a PP structure from polynomial coefficients and knot.
#'
#' @param coeff Coefficient matrix (kn x (degree+1))
#' @param tn Knot vector of length kn+1
#' @return List with components \code{coefficients} and \code{knot}
#' @keywords internal

makpp<-function(coeff,tn){
  #coeff is an array of dim: kn,(d+1)
  #kn=length(tn)
  #this is independent from the notation convention order
  #particular case one single piece : length(tn)=2
  if (length(tn)==2){kn=(length(tn)-1)}
  if ((length(tn)!=2) & (!is.null(dim(coeff)) ) )#if the degree is >0 or more than 1 piece
             {kn=dim(coeff)[1]
               o=dim(coeff)[2]}
  if ((length(tn)!=2) & (is.null(dim(coeff)))){kn<-length(coeff)}# case degree=0

  if (length(tn) != (kn+1)){
    stop("length of coeff and number of knot do not match")

  }
  else{
    return(list(coefficients=(coeff),knot=tn))
  }
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
  kn=length(tn)
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
#' Transforms a spline regression result into a callable function that
#' can be evaluated at any point, while preserving access to parameters.
#'
#' @param result A list returned by quantile_spline or one of the degree-specific functions
#' @return A function that can be called as `result(x)` and has attributes
#'         for `degree`, `knot`, `coefficients`, and `status`
#' @export
make_spline <- function(result) {
  # Extract components
  coeff <- result$coefficients
  deg <- result$degree
  knots <- result$knot
  # Create the callable function
  spline_func <- function(x_values, Bvalues=NULL) {
    # Build the spline object structure expected by spline_eval
    spline_obj <- list(
      coefficients = coeff,
      degree = deg,
      knot = knots
    )
    class(spline_obj) <- "Bspline"

    # Evaluate the spline
    spline_eval(spline_obj, x_values, Bvalues=Bvalues)
  }

  # Attach parameters as attributes (accessible via attr())
  attr(spline_func, "degree") <- deg
  attr(spline_func, "knot") <- knots
  attr(spline_func, "coefficients") <- coeff
  attr(spline_func, "status") <- result$status
  attr(spline_func, "value") <- result$value
  attr(spline_func, "y_mean") <- result$y_mean

  # Store the full result as an attribute
  attr(spline_func, "result") <- result

  # Set class for print method
  class(spline_func) <- c("callable_spline", "function")

  return(spline_func)
}

#' Print method for callable spline
#'
#' @param x A callable spline object
#' @param ... Additional arguments
#' @export
print.callable_spline <- function(x, ...) {
  cat("Callable Spline Object\n")
  cat("  Degree:", attr(x, "degree"), "\n")
  cat("  Knots:", length(attr(x, "knot")), "knots\n")
  cat("  Coefficients:", length(attr(x, "coefficients")), "\n")
  cat("  Status:", attr(x, "status"), "\n")
  cat("  Usage: my_spline(x) to evaluate\n")
  invisible(x)
}
