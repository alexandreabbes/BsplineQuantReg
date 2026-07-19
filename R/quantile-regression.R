# SplineConstQuantRegBs3, apply_karlin_constraints

#' Convert B-spline to derivative coefficients
#'
#' Converts a B-spline basis to normalized first and second derivative
#' coefficients on each interval.
#'
#' @param tn Knot vector (effective partition, not extended)
#' @param degree Spline degree (default = 3)
#' @param x_values Evaluation points for design matrix (0 = no evaluation)
#' @param verbose boolean FALSE (default) or TRUE.
#' @return A list containing:
#'   \item{d0}{Design matrix (if x_values provided)}
#'   \item{d1}{First derivative coefficients [a3, a2, a1] for each interval}
#'   \item{d2}{Second derivative values at knot}
#' @export
bspline_to_deriv_coeffs_pp <- function(tn,degree = 3,x_values=0, verbose=FALSE) {

  # create  basis with create.bspline.basis
  kn <- length(tn) - 1
  # Nombre correct de fonctions de base: kn + degree+1
  nbasis <- kn + degree

  norder <- degree + 1  # 4 pour cubique

  sn=c(tn[1]*rep(1,degree),tn,tn[kn+1]*rep(1,degree)) # this is the extended knot sequence
  BB <- Bspline_base(sn,degree)
  basis<-BB$base

  N <-BB$n_splines
  if (verbose) {
    message("Number of  basis functions", N, "\n")
  }
  # Matrix of normalised coefficient derivativs
  deriv_coeffs <- array(0, dim = c(kn, N, 3))
  deriv2_val<-array(0, dim = c(kn+1, N))
  deriv3_val<-array(0, dim = c(kn, N))
  for (j in 1:N){
    for (nu in (degree+1):(kn+degree))
    {
      h=sn[nu+1]-sn[nu]
      c3<-basis[j,nu,1]
      c2<-basis[j,nu,2]
      a3<-3*c3*h^2
      a2<-2*basis[j,nu,2]*h
      a1<-basis[j,nu,3]
      # coeffs_poly est [a3, a2, a1, a0] a0+a1*x+a_2*x^2+a3*x^3
      deriv_coeffs[nu-degree,j,]<-c(a3,a2,a1)
      deriv2_val[nu-degree,j]<-c2
      deriv3_val[nu-degree,j]<-c3 #up to a factor 6, but the sign is the same.
    }

    # for the last knot the second deriv is an affine function
    # p+m(t-t_{kn-1}) h is the last intervall space
    p<-2*basis[j,nu,2]
    m<-6*basis[j,nu,1]
    deriv2_val[nu-degree+1,j]=p+m*h
  }
  if (length(x_values)!=1){yvalues=bs_direct(BB,x_values)}
  else {yvalues=0}
  return(list(d0=yvalues,d1=deriv_coeffs, d2=deriv2_val,d3=deriv3_val))
}

#' Local Wrapper function to apply_karlin_quadratic
apply_karlin_constraints<-function(p2, p1, p0, z0,verbose=FALSE)
  {
  apply_karlin_quadratic(p2=p2, p1=p1, p0=p0, z0=z0,sign=1,verbose=verbose)
  }

#' Constrained quantile regression with cubic splines
#'
#' Performs quantile regression using cubic B-splines, with optional
#' monotonicity constraints (via Karlin-Studden) and convexity constraints.
#'
#' @param xtab Predictor vector (x)
#' @param ytab Response vector (y)
#' @param knot Knot vector or number of knot (quantiles are then used)
#' @param tau Quantile (between 0 and 1)
#' @param monot Monotonicity constraint vector per interval:
#'        1 = increasing, -1 = decreasing, 0 = unconstrained. If scalar, repeated.
#' @param convcons Convexity constraint vector per knot:
#'        1 = convex, -1 = concave, 0 = unconstrained. If scalar, repeated.
#' @param solver CVXR solver to use (default = "CLARABEL")
#' @param weight Observation weights (default = 1 for all)
#' @param verbose boolean FALSE (default) or TRUE.
#' @return A list containing:
#'   \item{coefficients}{B-spline coefficients (including y mean)}
#'   \item{degree}{Spline degree (always 3)}
#'   \item{ext_knot}{ext_knot vector used}
#'   \item{int_knot}{knot vector}
#' @examples
#' #optional set.seed(42)
#' x <- seq(0, 1, length=100)
#' y <- 2*x + sin(6*pi*x)/2 + rnorm(100, 0, 0.05)
#' knot <- quantile(x, probs=seq(0,1,length.out=10))
#'
#' # Median quantile regression without constraints
#' fit <- SplineConstQuantRegBs3(x, y, knot, tau=0.5)
#'
#' # With increasing monotonicity constraint
#' fit_monot <- SplineConstQuantRegBs3(x, y, knot, tau=0.5, monot=1)
#'
#' # With convexity constraint
#' fit_convex <- SplineConstQuantRegBs3(x, y, knot, tau=0.5, convcons=1)
#'
#' @seealso
#' Related R packages:
#' \itemize{
#'   \item \code{quantreg} - Quantile regression with linear programming
#'   \item \code{cobs} - Constrained B-sines (linear and quadratic only)
#' }
#'
#' Other implementations:
#' \itemize{
#'   \item MATLAB/Python versions: \url{https://github.com/alexandreabbes/Constrained-Quantile-Regression-with-cubic-splines}
#' }
#' @references
#' \itemize{
#'   \item Abbes, A. (2025). \emph{Quantile regression with cubic polynomial splines under shape constraints with applications}
#'         . Zenodo.
#'         \doi{10.5281/zenodo.16999784}
#'   \item de Boor, C. (1978). \emph{A Practical Guide to Splines}. Springer-Verlag.
#'         \doi{10.1007/978-1-4612-6333-3}
#'   \item Karlin, S., & Studden, W. J. (1966). \emph{Tchebycheff Systems: With
#'         Applications in Analysis and Statistics}. Interscience.
#'   \item Koenker, R., & Bassett, G. (1978). Regression Quantiles.
#'         \emph{Econometrica}, 46(1), 33-50. \doi{10.2307/1913643}
#'   \item Koenker, R. (2025). quantreg: Quantile Regression. R package version 5.99.
#'         \url{https://CRAN.R-project.org/package=quantreg}
#'   \item Ng, P., & Maechler, M. (2024). cobs: Constrained B-Splines.
#'         R package version 1.3-8. \url{https://CRAN.R-project.org/package=cobs}
#' }
#' @export


SplineCubicQuant<- function(xtab, ytab, knot, tau,
                                   monot = 0,
                                   convcons=0,
                                   der3cons=0,
                                   solver = "CLARABEL", weight = NULL,
                                   verbose=FALSE)
{

  if (is.null(weight)) {
    weight <- rep(1, length(xtab))
  }

  ordre <- order(xtab)
  xtab <- xtab[ordre]
  ytab <- ytab[ordre]
  weight <- weight[ordre]

  n <- length(xtab)

  if (length(knot) == 1 && is.numeric(knot))
  {
    kn <- knot - 1
    knot <- quantile(xtab, probs = seq(0, 1, length.out = kn + 1))
  }

  kn <- length(knot) - 1
  if (verbose) {
    message("knot:", knot, "\n")
  }


  if (verbose) {
  message("Monotonicity constraints (Karlin):", monot, "\n")
  }
  boundary_knot <- range(knot)
  degree=3
  N=length(knot)+3-1

  int_knot=knot[2:kn]

  # Calcul des coefficients normalises des derivees
  deriv_spline <- bspline_to_deriv_coeffs_pp(knot, degree = 3,x_values=xtab,verbose=verbose)
  deriv_coeffs <-deriv_spline$d1
  deriv_coeffs2<-deriv_spline$d2
  deriv_coeffs3<-deriv_spline$d3
  B<-deriv_spline$d0
  B=t(B)
  y_mean <- mean(ytab)
  ytab_centered <- ytab - y_mean

  alpha <- Variable(N)
  #Contraintes monotones
  # Variables auxiliaires z
  residuals <- ytab_centered - (B %*% alpha)

  u_plus <- pos(residuals)
  u_minus <- pos(-residuals)
  weighted_loss <- sum(weight * (tau * u_plus + (1 - tau) * u_minus))

  objective <- Minimize(weighted_loss)

  constraints <- list()
  z_vars <- list()

  if (any(monot != 0)) {
    if (length(monot) == 1) {
      monot <- rep(monot, kn)
    }
    for (i in 1:(kn)) {
      if (monot[i] != 0) {
        z_vars[[i]] <- Variable(1, name = paste0("z", i))
        a_coef=sum(deriv_coeffs[i,,1]*alpha) *monot[i]
        b_coef=sum(deriv_coeffs[i,,2]*alpha) *monot[i]
        c_coef=sum(deriv_coeffs[i,,3]*alpha) *monot[i]
        #a*x^2+b*x+c
        CK<-apply_karlin_constraints(a_coef,b_coef,c_coef,z_vars[[i]],verbose=verbose)
        constraints<-c(constraints,CK)
      }
    }
  }

  #"contraintes convexes
  # eliminate the null (unconstrained) case
  if (any(convcons !=0)){
    if (length(convcons) == 1) {
      convcons <- rep(convcons, (kn+1))}
    print(convcons)
    CV<-list((convcons*(deriv_coeffs2 %*% alpha))>=0) # Very simple, only use the second derivatives at the knot.
    constraints<-c(constraints,CV)
  }
#3rd derivative constraints
  if (any(der3cons!=0)){
    der3cons<-rep(der3cons,kn)
    DER3<-list(der3cons*alpha*deriv_coeffs3>0)
  }
  problem <- Problem(objective, constraints)

  result <- NULL
  solvers_to_try <- c(solver, "CLARABEL", "OSQP", "ECOS", "SCS")

  for (s in unique(solvers_to_try)) {
    if (verbose) {
      message("attempt with  solver:", s, "\n")
      }
    result <- tryCatch(
      psolve(problem, solver = toupper(s), verbose=verbose),
    error = function(e) {warning("Missed:", e$message, "\n")
      NULL} )
    if (!is.null(result) && !is.null(value(alpha))) {
      if (verbose) {
        message("Solveur succeeded:", s, "\n")
      }
      break}
  }

  if (is.null(result) || is.null(value(alpha))) {
    warning("Optimisation did not converge with any available solver")
    return(NULL)
  }

  alpha_val <- value(alpha)+y_mean
  if (verbose) {
    message(" Statut:", result$status, "\n",
            "Valeur objectif:", result$value, "\n",
             "Coefficients alpha (range):", range(alpha_val), "\n")
}

  return(list(
    #spline = spline_result,
    coefficients = alpha_val,
    degree=3,
    #basis_matrix = B,
    knot = knot,
    int_knot = knot
  ))
}

#' Wrapper function for Cubic spline quantile regression
#'
#' This is a convenience wrapper for SplineCubicQuant
#' with convention of previous versions  SplineConstQuantRegBs3.
#'
#' @inheritParams SplineCubicQuant
#' @return Same as SplineCubicQuant
#' @export

SplineConstQuantRegBs3<-function(xtab, ytab, knot, tau,
                                         monot = 0,
                                         convcons=0,
                                         solver = "CLARABEL", weight = NULL,
                                         verbose=FALSE)
          {SplineCubicQuant(xtab, ytab, knot, tau,
          monot = monot,
          convcons=convcons,
          solver = solver, weight = weight,
          verbose=verbose)}
