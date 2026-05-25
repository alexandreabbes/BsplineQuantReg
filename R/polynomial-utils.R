# polymul, polyadd, polyderiv, poly_eval, reduce_pol
# change_polynomial_base_taylor

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

polymul <- function(p1, p2,ord=0,verbose=FALSE)
{
  if ((sum(abs(p1))==0) || (sum(abs(p2))==0)){ return(0)}
  else
  {
    p1=rev(p1)#invert for decreasing convention
    p2=rev(p2)
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
    res=rev(res)#in standard notation
    res=reduce_pol(res,verbose)
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

polyadd<- function(p1, p2,verbose) {#in stnd notation
  p1=rev(p1)
  p2=rev(p2)
  l1 <- length(p1)
  l2 <- length(p2)
  if (l1>l2) {p=p2
  p2=p1
  p1=p
  l=l2
  l2=l1
  l1=l}
  p1=c(p1,rep(0,(l2-l1)))
  sp1p2=p1+p2
  res=rev(p1+p2)
  res=reduce_pol(res,verbose)
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
  coeffs_a=rev(coeffs_a)
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
#' @param p Polynomial coefficient vector(coef in decreasing order)
#' @param verbose boolean FALSE (default) or TRUE.
#' @return Reduced vector (without leading zeros)
#' @examples
#' reduce_pol((c(0, 1, 1)) # returns c(1, 1) #since 0x^2+x+1=x+1
#' @export
reduce_pol<-function(p,verbose=FALSE){
  l=length(p)
  k=1
  while (p[k]==0 & k<l){k=k+1}
  if (vebose){message("removed ",k," useless zeroes")}
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

poly_eval<-function(p,xvalues){
  p=rev(p)
  #we evaluate the values in the convention p=c(p0,p1,p2)
  #represent the polynomial p0+p1x+p2*x^2
  y=c()
  d=length(p)-1
  for (x in xvalues){
    val=0
    for (i in 0:d){
      val=val+p[i+1]*x^i
    }
    y=c(y,val)
  }
  return(y)
}

#' Polynomial derivative
#'
#' Computes the derivative of order \code{der} of a polynomial.
#'
#' @param p Coefficient vector (decreasing powers)
#' @param der Derivative order (default = 1)
#' @return Coefficients of the derivative polynomial
#' @examples
#' # P(x) = x^2 -> P'(x) = 2x
#' polyderiv(c(1, 0, 0), 1) # returns c(2, 0)
#' @export

polyderiv<-function(p,der=1) # this is a working equivalent to polyder function,
  #because the order higher than 1 do not work  in the R function polyder
{#It is here in Std Not
  if (der==0)   {return(p)}
  else {
    l=length(p)
    #p=rev(p)#reverse
    A=array(data=0,c(l,l))
    for (i in 1:(l-1)){A[i+1,(i)]=(l-i)}
    D=A
    if (der>1){for (i in 2:der){D=A%*%D}} #compute the d-th power of A
    q=D%*%p
    q=q[(1+der):(l)]
    #return(rev(q)) # reverse back to match  stnd convention
    return(q)
  }
}
