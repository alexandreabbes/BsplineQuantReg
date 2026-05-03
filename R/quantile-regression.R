# rhotau, SplineCubicQuantBspkn3, apply_karlin_constraints




#' Convert B-spline to derivative coefficients
#'
#' Converts a B-spline basis to normalized first and second derivative
#' coefficients on each interval.
#'
#' @param tn Knot vector (effective partition)
#' @param degree Spline degree (default = 3)
#' @param xvalues Evaluation points for design matrix (0 = no evaluation)
#' @return A list containing:
#'   \item{d0}{Design matrix (if xvalues provided)}
#'   \item{d1}{First derivative coefficients [a1, a2, a3] for each interval}
#'   \item{d2}{Second derivative values at knots}
#' @export
bspline_to_deriv_coeffs_pp <- function(tn,degree = 3,xvalues=0) {

  # Créer la base avec create.bspline.basis
  kn <- length(tn) - 1

  # Nombre correct de fonctions de base: kn + degree+1
  nbasis <- kn + degree

  norder <- degree + 1  # 4 pour cubique

  # Créer la base avec create.bspline.basis
  sn=c(tn[1]*rep(1,degree),tn,tn[kn+1]*rep(1,degree)) # this is the extended knots sequence
  BB <- Bspline_base(sn,degree)
  basis<-BB$base

  N <-BB$n_splines

  cat("Nombre de fonctions de base", N, "\n")

  # Matrice des coefficients pour la dérivée première normalisée

  deriv_coeffs <- array(0, dim = c(kn, N, 3))
  deriv2_val<-array(0, dim = c(kn+1, N))

  for (j in 1:N){
    for (nu in (degree+1):(kn+degree))
    {

      h=sn[nu+1]-sn[nu]
      a3<-3*basis[j,nu,1]*h^2
      a2<-2*basis[j,nu,2]*h
      a1<-basis[j,nu,3]
      c1<-2*basis[j,nu,2]

      # coeffs_poly est [a0, a1, a2, a3] a0+a1*x+a_2*x^2...
      deriv_coeffs[nu-degree,j,]=c(a1,a2,a3)
      deriv2_val[nu-degree,j]=c1
    }

    # for the last knot the second deriv is an affine function
    # c1+c2(t-t_{kn-1}) h is the last intervall space
    c2=6*basis[j,nu,4]
    deriv2_val[nu-degree+1,j]=c1+c2*h
  }
  if (length(xvalues)!=1){yvalues=bs_direct(BB,xvalues)}
  else {yvalues=0}
  return(list(d0=yvalues,d1=deriv_coeffs, d2=deriv2_val))
}

#' Karlin-Studden constraints for monotonicity
#'
#' Applies Karlin-Studden SOCP constraints to ensure monotonicity of a
#' quadratic polynomial on the interval [0,1].
#'
#' @param p2 Coefficient of u^2
#' @param p1 Coefficient of u
#' @param p0 Constant term
#' @param z0 Auxiliary SOCP variable
#' @return List of CVXR constraints
#' @export
apply_karlin_constraints <- function(p2, p1, p0, z0) {
  # P2, p1, p0 sont les coefficients du polynôme quadratique: p2*u^2 + p1*u + p0
  # Dans la notation de l'article

  constraints <- list()
  constraints <- c(constraints, list(z0 >= 0))

  K1_vec <- vstack(p0 + p2 - z0,p1-z0)
  K2_vec <- (p0 +p2+ z0)

  constraints <- c(constraints, list(K2_vec >= p_norm(K1_vec, 2)))

  return(constraints)
}

#' Régression quantile avec splines cubiques - Version avec contraintes de Karlin
#' pour la monotonie
#' et contraintes de convexite aux noeuds.
#'


#' Constrained quantile regression with cubic splines
#'
#' Performs quantile regression using cubic B-splines, with optional
#' monotonicity constraints (via Karlin-Studden) and convexity constraints.
#'
#' @param xtab Predictor vector (x)
#' @param ytab Response vector (y)
#' @param knots Knot vector or number of knots (quantiles are then used)
#' @param tau Quantile (between 0 and 1)
#' @param monot Monotonicity constraint vector per interval:
#'        1 = increasing, -1 = decreasing, 0 = unconstrained. If scalar, repeated.
#' @param convcons Convexity constraint vector per knot:
#'        1 = convex, -1 = concave, 0 = unconstrained. If scalar, repeated.
#' @param solver CVXR solver to use (default = "CLARABEL")
#' @param weight Observation weights (default = 1 for all)
#' @return A list containing:
#'   \item{coefficients}{B-spline coefficients (including y mean)}
#'   \item{degree}{Spline degree (always 3)}
#'   \item{knots}{Knot vector used}
#'   \item{int_knots}{Same as knots (compatibility)}
#' @examples
#' set.seed(42)
#' x <- seq(0, 1, length=100)
#' y <- 2*x + sin(6*pi*x)/2 + rnorm(100, 0, 0.05)
#' knots <- quantile(x, probs=seq(0,1,length.out=10))
#'
#' # Median quantile regression without constraints
#' fit <- SplineCubicQuantBspkn3(x, y, knots, tau=0.5)
#'
#' # With increasing monotonicity constraint
#' fit_monot <- SplineCubicQuantBspkn3(x, y, knots, tau=0.5, monot=1)
#'
#' # With convexity constraint
#' fit_convex <- SplineCubicQuantBspkn3(x, y, knots, tau=0.5, convcons=1)
#' @export


SplineConstQuantRegBs3 <- function(xtab, ytab, knots, tau,
                                   monot = 0,
                                   convcons=0,
                                   solver = "CLARABEL", weight = NULL)
{

  if (is.null(weight)) {
    weight <- rep(1, length(xtab))
  }

  ordre <- order(xtab)
  xtab <- xtab[ordre]
  ytab <- ytab[ordre]
  weight <- weight[ordre]

  n <- length(xtab)

  if (length(knots) == 1 && is.numeric(knots))
  {
    kn <- knots - 1
    knots <- quantile(xtab, probs = seq(0, 1, length.out = kn + 1))
  }

  kn <- length(knots) - 1

  cat("knots:", knots, "\n")

  if (length(monot) == 1) {
    monot <- rep(monot, kn)
  }
  cat("Monotonicity constraints (Karlin):", monot, "\n")
  boundary_knots <- range(knots)
  degree=3
  N=length(knots)+3-1
  #B <- bs(xtab, knots = knots, degree = degree)
  #B=B[,1:N]

  #        Boundary.knots = boundary_knots, intercept = TRUE)
  #kn=length(knots)-1
  #N <- kn+degree
  int_knots=knots[2:kn]
  #cat("Nombre de fonctions de base:", N, "\n")

  # Calcul des coefficients normalisés des dérivées
  deriv_spline <- bspline_to_deriv_coeffs_pp(knots, degree = 3,xvalues=xtab)
  deriv_coeffs <-deriv_spline$d1
  deriv_coeffs2<-deriv_spline$d2
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
    for (i in 1:(kn)) {
      if (monot[i] != 0) {
        z_vars[[i]] <- Variable(1, name = paste0("z", i))
        a_coef=sum(deriv_coeffs[i,,1]*alpha) *monot[i]
        b_coef=sum(deriv_coeffs[i,,2]*alpha) *monot[i]
        c_coef=sum(deriv_coeffs[i,,3]*alpha) *monot[i]

        CK<-apply_karlin_constraints(a_coef,b_coef,c_coef,z_vars[[i]])
        constraints<-c(constraints,CK)
      }
    }
  }

  #"contraintes convexes
  if (length(convcons) == 1) {
    convcons <- rep(convcons, (kn+1))
  }
  # eliminate the null (unconstrained) case

  if (any(convcons !=0)){
    CV<-list(convcons*(deriv_coeffs2 %*% alpha)>=0) # Very simple, only use the second derivatives at the knots.
    constraints<-c(constraints,CV)
  }

  problem <- Problem(objective, constraints)

  result <- NULL
  solvers_to_try <- c(solver, "CLARABEL", "OSQP", "ECOS", "SCS")

  for (s in unique(solvers_to_try)) {
    cat("Tentative avec solveur:", s, "\n")
    result <- tryCatch({
      solve(problem, solver = toupper(s), verbose = FALSE)
    }, error = function(e) {
      cat("Échec:", e$message, "\n")
      NULL
    })

    if (!is.null(result) && !is.null(result$getValue(alpha))) {
      cat("Solveur réussi:", s, "\n")
      break
    }
  }

  if (is.null(result) || is.null(result$getValue(alpha))) {
    warning("L'optimisation n'a pas convergé avec aucun solveur")
    return(NULL)
  }

  alpha_val <- result$getValue(alpha)+y_mean

  #  cat("Statut:", result$status, "\n")
  #  cat("Valeur objectif:", result$value, "\n")
  #  cat("Coefficients alpha (range):", range(alpha_val), "\n")


  return(list(
    #spline = spline_result,
    coefficients = alpha_val,
    degree=3,
    #basis_matrix = B,
    knots = knots,
    int_knots = knots
  ))
}




