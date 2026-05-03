# Omega, Bspline_base, Bspline_deriv



#' Omega function for De Boor recursion
#'
#' Computes the affine element used in the recursive De Boor algorithm
#' for B-spline construction.
#'
#' @param s Extended knot vector
#' @param j Knot index
#' @param o Spline order (degree + 1)
#' @return Vector of coefficients [alpha, beta] representing (alpha*t + beta)
#' @keywords internal

Omega<-function(s,j,o)#t: knots in the base t-t[j]

{
  if(s[j]==s[j+o-1]){w<-c(0,0)}
  else
  {alpha=1/(s[j+o-1]-s[j])
  beta=(-s[j]*alpha)
  #w<-c(beta,alpha) } # This is the affine element wjl of deBoor
  w<-c(alpha,beta)#MatNot
  #wik=(t-s_i)/(s_{i+k-1}-s_{i})=alpha*t+beta
  print(w)
  return(w)
  }
}


#' Build B-spline basis in piecewise polynomial form
#' Computes local polynomial coefficients for each B-spline basis function
#' on each interval. Polynomials are expressed in the canonical basis
#' Uses De Boor's recursion formula.
#' @param sn Extended knot vector (including endpoint repetitions ; This means if t0..t_{kn} it the set of knots
#' then sn should be given as a vector with  "degree" times t_0 and t_{kn} at the begining
#'  and the ends.  its length is number of intervals+1+2*degree)
#' @param degree B-spline degree (default = 3 for cubic)
#' @param der Derivative order (0 = original basis)
#' @return A list containing:
#'   \item{base}{Coefficients in the local bases, in the form of an 3-d array [j,nu,coeff], j : the number of the spline in the basis,
#'   nu: the number of the interval in the extended notation, coeff : the coefficients in decreasing order
#'   convention on the local bases (t-s_nu)^l, l=3..1
#'   base[j,,] is a matrix of piecewise polynomial function compatible with the pp-form.}
#'   \item{base0}{Coefficients in canonical basis (centered at 0) (1, t, t^2, t^3) centered at the interval origin.}
#'   \item{knots}{Extended knot vector}
#'   \item{int_knots}{Internal knots (effective partition including ends)}
#'   \item{degree}{Spline degree}
#'   \item{n_splines}{Number of basis functions}
#'   \item{deriv_order}{Applied derivative order}
#' @examples
#' sn <- c(0,0,0,0,1,2,3,4,5,5,5,5)
#' basis <- Bspline_base(sn, degree=3)
#' x=(0:(5*100))/100
#' y=bs_direct(basis,x)
#' matplot(x,y)
#' #or simple:
#' #view_basis(basis)
#'
#' @export

Bspline_base<-function(sn,degree=3,der=0)
{
  tn=sn[(degree+1):(length(sn)-degree)] #effective knots partition
  kn=length(tn)-1 # tn is the interior knots without the extended partition.

  n_intervals<-kn+2*degree #Nb extended intervals
  n_splines<-kn+degree
  B<-array(0,dim=c((degree+1),n_splines,n_intervals,(degree+1))) # B is the initial B-spline basis : piecewise constant
  #for (i in (degree+1):(kn+degree)){B[1,i,i,1]<-1} # initialisation of the splines in increasing convention
  #with degree 0 wich has dimension kn
  # 1: degree 0, kn = Nb interval; kn= Nb elements in the basis of degree 0

  for (i in (degree+1):(kn+degree)){B[1,i,i,degree+1]<-1}#in decreasing convention
  if (degree>0){
    for (o in (2:(degree+1))){
      #k : dimension of local basis=deg+1
      for (j in (1:(n_splines))) {
        #go through the elements of the basis
        for (nu in ((degree+1):(n_intervals))) #go through the  pieces of the spline of order l
        {
          #Bjnu<-B[(o-1),j,nu,1:(o-1)] #inc. convention
          Bjnu<-reduce_pol(B[(o-1),j,nu,(degree+1-(o-1)):(degree+1)])#dec. convention
          if ((j+1)>n_splines){
            Bjpnu=0
            wjp1nu=0}
          else
          {
            #Bjp1nu=B[(o-1),(j+1),nu,1:(o-1)] #inc. convention
            Bjp1nu=reduce_pol( B[(o-1),(j+1),nu,(degree-o+2):(degree+1)]) #decr. convention
            #wjp1o=polyadd(c(1,0),-Omega(sn,(j+1),o))}#inc convent.
            wjp1o=polyadd(c(0,1),-Omega(sn,(j+1),o))}

          wjo=Omega(sn,j,o)

          term1=polymul(Bjnu,wjo)
          term2=polymul(Bjp1nu,wjp1o)

          sumterm=polyadd(term1,term2)

          #B[o,j,nu,(1:o)]<-sumterm
          B[o,j,nu,(degree-o+2):(degree+1)]<-sumterm
        }
      }
    }}

  #changement de base
  # because the coefficients of each polynomial piece
  #are computed on the canonical basis 1,x,x^2...x^degree
  Bn=array(0,dim=dim(B))
  for (o in 1:(degree+1)){
    for (j in 1:n_splines)
    {for (nu in 1:n_intervals)
    {
      Bn[o,j,nu,]<-change_polynomial_base_taylor(B[o,j,nu,],0,sn[nu])
    }}}


  Base0=B[(degree+1),,,]
  BaseL=Bn[(degree+1),,,]
  Base0=round(Base0,10)
  BaseL=round(BaseL,10)

  if (der!=0){
    Base0=Bspline_deriv(Base0,der = der)
    BaseL=Bspline_deriv(BaseL,der=der)
  }
  return(list(base =BaseL ,base0=Base0, knots = sn, int_knots=tn, degree = degree, n_splines = (n_splines),deriv_order=der ) )
  #  return (B)
}



#' Differentiate a B-spline basis
#'
#' Computes the basis of order \code{der} derivatives of a B-spline basis.
#'
#' @param bspline Object returned by \code{Bspline_base}
#' @param der Derivative order
#' @return A list similar to \code{Bspline_base} for the derivative basis
#' @export


Bspline_deriv<-function(bspline,der=2){
  #computes the derivative fo a Bspline basis
  Bn=bspline$base
  B0=bspline$base0
  n_splines=bspline$n_splines
  knots=bspline$knots

  degree=bspline$degree
  degree_der=max(degree-der,0)
  NS=length(knots)-1 #Nb extended intervals
  Bn_der=array(dim=c(n_splines,NS,max((degree_der+1),1) ),0)

  B0_der=array(dim=dim(Bn_der),0)

  for (j in 1:n_splines){
    for (nu in 4:NS){
      p=polyderiv(Bn[j,nu,],der)
      q=polyderiv(B0[j,nu,],der)
      if (!is.null(p)){
        print(c(j,nu))
        Bn_der[j,nu,]=p
        B0_der[j,nu,]=q
      }
    }
  }

  return(list(base =Bn_der, base0=B0_der, knots = knots, int_knots=bspline$int_knots, degree = degree_der, n_splines = (n_splines) ) )
}
