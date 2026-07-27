# Omega, Bspline_base, Bspline_deriv

#' Omega function for De Boor recursion
#'
#' Computes the affine element used in the recursive De Boor algorithm
#' for B-spline construction. for a given extended knot partition s(1),..
#' Omega(j,l)(x)=(s(j-x))/(s(j+l-1)-s(j)\), with l: order of the spline
#'
#' @param s Extended knot vector
#' @param j Knot index
#' @param o Spline order (degree + 1)
#' @return (Side function) A Vector of coefficients [alpha, beta] representing (alpha*t + beta)
#' @keywords internal
Omega<-function(s,j,o)#t: knot in the base t-t[j]

{
  if(s[j]==s[j+o-1]){w<-c(0,0)}
  else
  {alpha=1/(s[j+o-1]-s[j])
  beta=(-s[j]*alpha)
  w<-c(alpha,beta)#MatNot
  return(w)
  }
}

#' Build B-spline basis in piecewise polynomial form
#' Computes local polynomial coefficients for each B-spline basis function
#' on each interval. Polynomials are expressed in the canonical basis
#' Uses De Boor's recursion formula.
#' @param sn Extended knot vector (including endpoint repetitions)
#' This means if t0..tkn it the set of knot
#' then sn should be given as a vector with  "degree" times t_0 and t_kn
#' at the beginning
#'  and the ends.  its length is number of intervals+1+2*degree.
#' @param degree 'B-spline' degree (default = 3 for cubic)
#' @param der Derivative order (0 = original basis)
#' @param verbose boolean FALSE (default) or TRUE.
#' @return A list containing:
#'   \item{base}{Coefficients in the local bases, in the form of an 3-d array [j,nu,coeff],
#'    j : the number of the spline in the basis,
#'   nu: the number of the interval in the extended notation,
#'   coeff : the coefficients in decreasing order
#'   Decresing (matlab) convention on the local bases (t-s_nu)^l, l=3.2.1
#'   base[j,,] is a matrix of piecewise polynomial function compatible with the pp-form.}
#'   \item{base0}{Coefficients in canonical basis (centered at 0) (1, t, t^2, t^3) centered at the interval origin.}
#'   \item{ext_knot}{Extended knot vector}
#'   \item{knot}{Effective knot partition including ends)}
#'   \item{degree}{Spline degree}
#'   \item{n_splines}{Number of basis functions}
#'   \item{deriv_order}{Applied derivative order}
#' @examples
#' sn <- c(0,0,0,0,1,2,3,4,5,5,5,5)
#' basis <- Bspline_base(sn, degree=3)
#' x=(0:(5*100))/100
#' y=bs_direct(basis,x)
#' matplot(x,t(y))
#' #or simple:
#' #view_basis(basis)
#' @export
#'
#'
Bspline_base<-function(sn,degree=3,der=0,verbose=FALSE)
{ if (verbose) {message("Constructing the B-spline basis, \n Degree = ", degree, " \n Extended knots partition: ", sn)}

  tn <- tryCatch({
    sn[(degree + 1):(length(sn) - degree)]
  }, error = function(e) {
    stop("Error extracting knots: ", e$message)
  }) #effective knot partition
  kn=length(tn)-1 # tn is the list of knot without the extended partition.
  n_ext_intervals<-kn+2*degree #Nb extended intervals
  n_splines<-kn+degree
  B<-array(0,dim=c((degree+1),n_splines,n_ext_intervals,(degree+1))) # B is the initial B-spline basis : piecewise constant

  for (i in (degree+1):(kn+degree)){B[1,i,i,degree+1]<-1}#in decreasing convention
  if (degree>0){
    for (o in (2:(degree+1)))#first Bspline base
      #to be computed: order2=degree 1.
      {
      #k : dimension of local basis=deg+1
      for (j in (1:(n_splines))) {
        #go through the elements of the basis
        for (nu in ((degree+1):(n_ext_intervals))) #go through the  pieces of the spline of order l
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
    {for (nu in 1:n_ext_intervals)
    {
      Bn[o,j,nu,]<-change_polynomial_base_taylor(B[o,j,nu,],0,sn[nu])
    }}}


  Base0=B[(degree+1),,,]  #basis x^i at each knot
  BaseL=Bn[(degree+1),,,] #Local basis (x-t_nu)^i at knot nu
  Base0=round(Base0,10)
  BaseL=round(BaseL,10)
  basis<-list(base =BaseL ,
              base0=Base0,
              ext_knot = sn,
              knot=tn,
              degree = degree,
              n_splines = (n_splines),
              deriv_order=0 )
  class(basis)<-c('bspline_basis',class(basis))
  if (der==0){
    return(basis)
  }else
  {return(Bspline_base_deriv(basis,der))}

}


#' Differentiate a B-spline basis
#'
#' Computes the basis of order \code{der} derivatives of a B-spline basis.
#'
#' @param Bsbasis Object returned by \code{Bspline_base} or an equivalent
#'  coherent list with enough parameters.
#' @param der Derivative order
#' @param verbose boolean FALSE (default) or TRUE.
#' @return A list similar to \code{Bspline_base} for the derivative basis
#' @export
Bspline_base_deriv<-function(Bsbasis,der=1,verbose=FALSE){
  #computes the derivative for a Bspline basis
  if (is.null(der)){der<-0}
  Bn=Bsbasis$base
  B0=Bsbasis$base0

  n_splines=Bsbasis$n_splines
  ext_knot=Bsbasis$ext_knot
  der0<-Bsbasis$deriv_order

  degree=Bsbasis$degree
  degree_der=max(degree-der,0)
  NS=length(ext_knot)-1 #Nb extended intervals
  Bn_der=array(dim=c(n_splines,NS,max((degree_der+1),1) ),0)

  B0_der=array(dim=dim(Bn_der),0)

  if (der<=degree){
    #otherwise no work is needed, only zeros

  for (j in 1:n_splines){
    for (nu in (degree+1):NS){
      p=polyderiv(Bn[j,nu,],der)
      q=polyderiv(B0[j,nu,],der)
      if (!is.null(p)){

        Bn_der[j,nu,]=p
        B0_der[j,nu,]=q
      }
    }
  }}

  D_basis<-list(base =Bn_der,
              base0=B0_der,
              ext_knot = ext_knot,
              knot=Bsbasis$knot,
              degree = degree_der,
              n_splines = (n_splines),
              deriv_order=der+der0)
  class(D_basis)<-c('bspline_basis',class(D_basis))
  return(D_basis)
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

Bsplinetopp <- function(Bspline,Bsbasis=NULL,callable=FALSE,verbose=FALSE) {
  #Convert a B-spline to a PP-polynomial
  if (inherits(Bspline,"callable_spline")){
    Bspline=get_parameters(Bspline)
    callable=TRUE}
    coeff=Bspline$coeff
    degree=Bspline$degree
    tn<-Bspline$knot
    sn<-c(rep(tn[1],degree),tn, rep( rev(tn)[1],degree )) #extended knots
  if (!is.null(Bspline$base)){BB<-Bspline$base}
  if (is.null(Bsbasis) && is.null(Bspline$base)) {
    Bsbasis<-Bspline_base(sn,degree=degree,verbose=verbose)
    BB<-Bsbasis$base
  }
    PP<-array( 0 , dim=c(length(tn)-1, degree+1) )
    for (nu in 1:(length(tn)-1)){
      PP[nu, ]<-t(BB[ , nu+degree ,  ])%*%coeff
    }
    PP<-makpp(PP,tn,callable=callable)
  }


