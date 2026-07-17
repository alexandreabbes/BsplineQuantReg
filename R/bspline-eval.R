# bs_direct, makpp, evalpp, spline_eval, view_spline
# Spline_der_knot, bspline_to_deriv_coeffs_pp

#' Evaluate a B-spline
#'
#' Evaluates a spline (linear combination of B-splines) at given points.
#'
#' @param Bspline Spline object (list with coefficients on the Bspline basis, degree, extended knot)
#' @param xvalues Vector of evaluation points
#' @return vector of same length as xvalues, with Spline values at the requested points
#' @examples
#'{ # Create and evaluate a spline
#' sn <- c(0,0,0,0,1,2,3,4,5,5,5,5)
#' basis <- Bspline_base(sn, degree=3)
#' basis$coefficients <- runif(basis$n_splines)
#' y <- spline_eval(basis, seq(0,5,length=100))}
#' @export

spline_eval<-function(Bspline,xvalues)
  #Bspline has a new type R container,
  #designed by its coefficients, the degree
  #and the knot
  #It is independent from the polynomial notation order
{
  knot=Bspline$knot #interior knot
  degree=Bspline$degree
  coeff=Bspline$coefficients
  #Bvalues=bs(xvalues,knot=knot,degree)  "can be used instead of the following lines
  #to accelerate the calculations
  t1=knot[1]
  tkn=rev(knot)[1] #last knot
  sn=c(rep(t1,degree),knot,rep(tkn,degree) ) #extended knot partition
  BB=Bspline_base(sn,degree=degree)
  Bvalues=t(bs_direct(BB,xvalues))
  N=length(knot)+degree-1
  yvalues=Bvalues[,1:N]%*%coeff
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
  # base calculee sous PP-forme : coeff des polynomes sur la base locale.
  #indep. du choix de la notation croissant/decroissant
  n_values=length(xvalues)
  knot=Basis$knot  # knot
  kn=length(knot)-1
  d=Basis$degree
  nsplines=Basis$n_splines
  if (d>0){bb=Basis$base[,(d+1):(nsplines),]}
  if (d==0){bb=Basis$base}
  yvalues=array(data=0,c(nsplines,n_values))
  print(c('les noeuds de la spline',knot))
  print(Basis$n_splines)
  for (j in 1:nsplines)
  {
    p=makpp(bb[j,,],tn=knot)
    yvalues[j,]<-evalpp(p,xvalues)
  }
  return(yvalues)
}


#' Evaluate a piecewise polynomial (PP) form
#'
#' Evaluates a piecewise polynomial function at given points.
#'
#' @param p List with components \code{ext_knot} (ext_knot) and \code{coeff}
#' @param xvalues Vector of evaluation points
#' @return Function values at the requested points
#' @keywords internal

evalpp<-function(p,xvalues){
  #this evaluates a polynomial p under the pp form,
  #p if given with its knot and the local coefficients
  #This funciton is independent from the order convention
  #for polynomials
  #The xvalues out of the knot give 0 in the corresponding yvalues
  tn=p$knot
  coeff=p$coefficients
  kn=length(tn)
  n_values=length(xvalues)

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
  xval=xvalues[xvalues==tn[kn]]
  #pval=c(pval,polyval(p=rev(poly_loc),xval)) #if use of R convention
  h=tn[kn]-tn[kn-1]
  pval=c(pval,poly_eval(poly_loc,h)) #use our convention for polynomial
  #zero1=rep(0,length(xvalues<tn[1]))
  #zero2=rep(0,length(xvalues>tn[kn]))
  #return(c(zero1,pval,zero2))
  return(pval)
}

#' Build a piecewise polynomial (PP) form
#'
#' Creates a PP structure from polynomial coefficients and knot.
#'
#' @param coef Coefficient matrix (kn x (degree+1))
#' @param tn Knot vector of length kn+1
#' @return List with components \code{coefficients} and \code{knot}
#' @keywords internal

makpp<-function(coeff,tn){
  #coeff is an array of dim: kn,(d+1)
  #kn=length(tn)
  #this is independent from the notation convention order
  kn=dim(coeff)[1]
  o=dim(coeff)[2]
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
#' @param xvalues Vector of evaluation points for plotting (by default 100 points are computed in the knot range)
#' @return No return value, called for side effects (generates a plot)
#' @export
#'
view_basis<-function(BB,xvalues=0)
{
  if (length(xvalues)==1){
    k=range(BB$ext_knot)
    xvalues=(k[1]:(k[2]*100))/100}

  yvalues=bs_direct(BB,xvalues)

  matplot(xvalues, t(yvalues))
}


#' Derivatives at knot of a B-spline
#'
#' Computes derivative values of a B-spline at knot (efficient because it
#' directly uses polynomial coefficients).
#'
#' @param Bspline Object returned by \code{Bspline_base}
#' @param der Derivative order (default = 1)
#' @return Matrix of derivative values (n_splines x n_knot)
#' @export

Spline_der_knot<-function(Bspline,der=1)
  #compute the values of a derivatives only at the knot
  #(simple, it only uses the coefficients)
{
  coeff=Bspline$base
  nsplines=Bspline$n_splines
  tn=Bspline$ext_knot
  kn=length(tn)
  m=Bspline$degree
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

