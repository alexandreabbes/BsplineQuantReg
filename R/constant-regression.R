#''
#' Performs quantile/mean regression using constant B-splines with monotonicity constraints.
#' No solver is needed, since for the degree 0, the B-splines in the basis are
#' independent piece-wise constants.
#'
#' @param xtab Predictor vector (x)
#' @param ytab Response vector (y)
#' @param knot Knot vector or number of knots
#' @param tau Quantile (between 0 and 1)
#' @param monot Monotonicity constraint vector: Monotonicity is expressed at knots
#'        1 || k>0 = increasing, -1 ||k<0 = decreasing, 0 = unconstrained.
#'        kn-1 constraints are taken into account, kn=length(knot)-1
#'        If monot is a single number, it is repaeated as uniform coonstraint.
#'        if length(monot)<kn-1 then it is completed with zeros.
#' @param weight Observation weights (default = 1 for all)
#' @param verbose logical; if TRUE, print progress messages
#' @param type_reg 'quantile' or 'mean_square' type of regression
#' i.e form of the objective.
#' @return A list containing coefficients, degree, and knots
#'
#' @export
SplineConstantQuant <- function(xtab,
                              ytab,
                              knot=NULL,
                              tau=0.5,
                              monot = 0,
                              weight = NULL,
                              verbose = FALSE,
                              type_reg='quantile') {
  if (is.null(knot)){knot=c(min(xtab),max(xtab))}
  if (is.null(weight)) {
    weight <- rep(1, length(xtab))
  }
(0)

  # Sort data
  ordre <- order(xtab)
  xtab <- xtab[ordre]
  ytab <- ytab[ordre]
  weight <- weight[ordre]

  n <- length(xtab)

  # Handle knots
  if (length(knot) == 1 && is.numeric(knot)) {
    kn <- knot - 1
    knot <- as.numeric(quantile(xtab, probs = seq(0, 1, length.out = kn + 1)))
  }

  kn <- length(knot) - 1
  degree <- 1
  N <- kn + degree  # Number of basis functions for linear = kn + 1
#manage constraints

  if (length(monot) == 1) {
    monot <- rep(monot, kn)
  }

  if (length(monot)<kn) monot<-c(monot,rep(0,kn-length(monot)))

  if (verbose) {
    message("=== Constant Quantile Regression (degree = 1) ===")
    message(sprintf("Knots: %d, Basis functions: %d", kn, N))
    if (type_reg=='mean_square') message('Mean square regression')
    else message('Quantile regression')
  }



  if (verbose) {
    message("Monotonicity constraints:", paste(monot, collapse = " "))
  }
  fonct<-function(x,tau, type){
#    print(type_reg)
     if (type=='mean_square'){
       return(mean(x))}
    else {
    return(quantile(x,tau))}

}

  coeff=rep(0,kn)
  sn<-knot
  B0=Bspline_base(sn,degree=0)
  for (i in 1:kn){
      local_ytab<-ytab[ knot[i]<xtab & xtab<knot[i+1]]
      local_weight<-weight[ knot[i]<xtab & xtab<knot[i+1]]
      local_weighted_ytab<-rep(local_ytab,local_weight)
      s=sign(monot[i])
      if (s!=0 & i>1)
      coeff[i]<-s*max(s*fonct(local_weighted_ytab,tau),s*coeff[i-1])
      else
      coeff[i]<-fonct(local_weighted_ytab,tau,type=type_reg)
      }

    Bspline<-list(
    coeff = coeff,
    degree = 0,
    knot = knot,
    result=NULL)
}



