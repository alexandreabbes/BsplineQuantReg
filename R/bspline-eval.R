# bs_direct, makpp, evalpp, spline_eval, view_spline
# Spline_der_knots, bspline_to_deriv_coeffs_pp

#' Evaluate a B-spline
#'
#' Evaluates a spline (linear combination of B-splines) at given points.
#'
#' @param Bspline Spline object (list with coefficients, degree, knots)
#' @param xvalues Vector of evaluation points
#' @return Spline values at the requested points
#' @examples
#' # Create and evaluate a spline
#' sn <- c(0,0,0,0,1,2,3,4,5,5,5,5)
#' basis <- Bspline_base(sn, degree=3)
#' basis$coefficients <- runif(basis$n_splines)
#' y <- spline_eval(basis, seq(0,5,length=100))
#' @export

spline_eval<-function(Bspline,xvalues)
  #Bspline has a new type R container,
  #désign by its coefficients, the degree
  #and the knots
  #It is independent from the polynomial notation order
{
  knots=Bspline$int_knots #interior knots
  degree=Bspline$degree
  coefficients=Bspline$coefficients
  #Bvalues=bs(xvalues,knots=knots,degree)  "can be used instead of the following lines
  #to accelerate the calculations
  t1=knots[1]
  tkn=rev(knots)[1] #last knot

  sn=c(rep(t1,degree),knots,rep(tkn,degree) ) #extended knot partition
  BB=Bspline_base(sn,degree=degree)
  Bvalues=t(bs_direct(BB,xvalues))
  N=length(knots)+degree-1
  yvalues=Bvalues[,1:N]%*%coefficients
  return(yvalues)
}



#' Direct evaluation of a B-spline basis
#'
#' Computes the values of all B-spline basis functions at given points.
#'
#' @param Basis Object returned by \code{Bspline_base}
#' @param xvalues Vector of evaluation points
#' @return Matrix of basis function values (n_splines x length(xvalues))
#' @export

bs_direct<-function(Basis,xvalues)
{
  #Calcule les valeurs d'une base Bspline.
  #comme la fonction bs de R, mais en utilisant la
  # base calculée sous PP-forme : coeff des polynômes sur la base locale.
  #indep. du choix de la notation croissant/decroissant
  n=length(xvalues)
  int_knots=Basis$int_knots#interior knots
  kn=length(knots)
  degree=Basis$degree
  nsplines=Basis$n_splines
  bb=Basis$base[,(degree+1):(nsplines),]
  #internal bspline knots
  yvalues=array(data=0,c(nsplines,n))

  for (j in 1:nsplines)
  {
    p=makpp(bb[j,,],tn=c(int_knots))
    yvalues[j,]<-evalpp(p,xvalues)
  }
  return(yvalues)
}


#' Evaluate a piecewise polynomial (PP) form
#'
#' Evaluates a piecewise polynomial function at given points.
#'
#' @param p List with components \code{knots} (knots) and \code{coefficients}
#' @param xvalues Vector of evaluation points
#' @return Function values at the requested points
#' @keywords internal

evalpp<-function(p,xvalues){
  #this evaluates a polynomial p under the pp form,
  #p if given with its knots and the local coefficients
  #This funciton is independent from the order convention   for polynomials
  tn=p$knots
  coeff=p$coefficients
  kn=length(tn)
  n=length(xvalues)

  pval<-c()
  for (i in 1:(kn-1))
  {
    xval=xvalues[(xvalues>=tn[i]) & (xvalues<tn[i+1])]
    poly_loc<-coeff[i,]
    # reverse our convention to match polyval convention
    #pval<-c(pval,polyval(p=rev(poly_loc),xval) )
    h=xval-tn[i]
    pval<-c(pval,poly_eval(poly_loc,h))
  }
  xval=xvalues[xvalues=tn[kn]]
  #pval=c(pval,polyval(p=rev(poly_loc),xval)) #if use of R convention
  h=tn[kn]-tn[kn-1]
  pval=c(pval,poly_eval(poly_loc,h)) #use our convention for polynomial
  return(pval)
}

#' Build a piecewise polynomial (PP) form
#'
#' Creates a PP structure from polynomial coefficients and knots.
#'
#' @param coef Coefficient matrix (kn x (degree+1))
#' @param tn Knot vector of length kn+1
#' @return List with components \code{coefficients} and \code{knots}
#' @keywords internal

makpp<-function(coef,tn){
  #coef is an array of dim: kn,(d+1)
  #kn=length(tn)
  #this is independent from the notation convention order
  kn=dim(coef)[1]
  o=dim(coef)[2]
  if (length(tn) != (kn+1)){
    return("length of coef and number of knots do not match")
    break
  }
  else{
    return(list(coefficients=(coef),knots=tn))
  }
}

#' Visualize a B-spline functions basis
#'
#' Plots all basis functions of a B-spline.
#'
#' @param Bspline Object returned by \code{Bspline_base}
#' @param xvalues Vector of evaluation points for plotting (by default 100 points are computed in the knot range)
#' @export
#'
view_basis<-function(Bspline,xvalues=0)
{
  if (xvalues==0){
    k=range(Bspline$knots)
    xvalues=(k[1]:(k[2]*100))/100}

  yvalues=bs_direct(Bspline,xvalues)

  matplot(xvalues, t(yvalues))
}


#' Derivatives at knots of a B-spline
#'
#' Computes derivative values of a B-spline at knots (efficient because it
#' directly uses polynomial coefficients).
#'
#' @param Bspline Object returned by \code{Bspline_base}
#' @param der Derivative order (default = 1)
#' @return Matrix of derivative values (n_splines x n_knots)
#' @export

Spline_der_knots<-function(Bspline,der=1)
  #compute the values of a derivatives only at the knots
  #(simple, it only uses the coefficients)
{
  coeff=Bspline$base
  nsplines=Bspline$n_splines
  tn=Bspline$int_knots
  kn=length(tn)
  m=Bspline$degree
  if (der>m){
    Der2_knots=array(data=0,c(nsplines,kn))
  }
  else{
    #Der2_knots=coeff[,,(der+1)]*factorial(der) #in increasing pol notation
    Der2_knots=coeff[,,(m-der+1)]*factorial(der) #in decreasing notation
    #computation of the last value
    h=tn[kn]-tn[kn-1]
    for (j in 1:nsplines)
    {
      p_kn_der=polyderiv(coeff[j,kn+m-1,],der)

      v_kn=poly_eval(p_kn_der,h)
      Der2_knots[j,kn+m]=v_kn
    }
  }
  return(t(Der2_knots))
}

