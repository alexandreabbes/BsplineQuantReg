#' Comprehensive test function
#'
#' Runs quantile regression tests with and without constraints,
#' and displays results. Demo function.
#' @param verbose boolean FALSE (default) or TRUE.
#' @param seed (default=NULL) value for the random generator.
#' @return No return value, produces plots.
#' @examples
#' test_karlin_simple()
#' @export
test_karlin_simple <- function(verbose=FALSE,seed=NULL) {
  # Store par
  old <- par(mfrow = c(2,2))
  # Restore par
  on.exit(par(old))

  if(!is.null(seed))
    {    set.seed(seed)}
  n_points <- 50
  xtab <- (0:n_points)/n_points

  # simple pscillating data
  #ytab <- -3 * xtab +sin(3*2*xtab*3.14)+ 0.2 * rnorm(n_points)
  ytab <- 2* xtab + 0.5 * sin(6 * pi * xtab) + 0.05 * rnorm(n_points+1)
  #ytab<-xtab*(1-xtab)
  kn <- 12

  n=7
  monot=c(rep(1,n),rep(0,(12-n)))
  #  monot=0
  knots <- quantile(xtab, probs = seq(0, 1, length.out = kn + 1))
  res<-SplineConstQuantRegBs3(xtab, ytab, knots, tau = 0.9,
                              monot = 0, solver = "OSQP")
  if (verbose) {
    message("\n===  PARTIAL INCREASING TEST ===\n")
  }

  res_croissant <- SplineConstQuantRegBs3(xtab, ytab, knots, tau = 0.5,
                                          monot = monot, convcons=0,solver = "OSQP")

  if (verbose) {
    message("\n=== DECREASING TEST ===\n")}

  res_decroissant <- SplineConstQuantRegBs3(xtab, ytab, knots, tau = 0.5,
                                            monot = -1, solver = "OSQP")
  if (verbose) {
    message("\n=== CONVEX  TEST ===\n")}
  res_convexe <- SplineConstQuantRegBs3(xtab, ytab, knots,monot=0,convcons=1, tau = 0.5)
  res_convexe <- SplineConstQuantRegBs3(xtab, ytab, knots,monot=0,convcons=1, tau = 0.5)
  #                                                   monot = -1, solver = "OSQP")
  # Visualisation
  par(mfrow = c(2, 2))
  x_eval <- seq(0, 1, length.out = 200)
  y_sans=spline_eval(res,x_eval)
  y_croiss=spline_eval(res_croissant,x_eval)
  y_decroiss=spline_eval(res_decroissant,x_eval)
  y_convexe=spline_eval(res_convexe,x_eval)

  # Croissant
  plot(xtab, ytab, pch = 16, cex = 0.5, col = "black",
       main = "increasing constraint (On the first inter-knots)")
  lines(x_eval,y_croiss,lwd=2)

  #  view_spline
  #lines(x_eval, y_croiss, col = "red", lwd = 2)

  abline(v = knots, col = "blue", lty = 2)

  # Decroissant
  plot(xtab, ytab, pch = 16, cex = 0.5, col = "black",
       main = "decreasing constraint")

  lines(x_eval, y_decroiss, col = "red", lwd = 2)

  abline(v = knots, col = "blue", lty = 2)


  # Convexe
  plot(xtab, ytab, pch = 16, cex = 0.5, col = "black",
       main = "convex constraint")
  if (!is.null(res_convexe)) {
    lines(x_eval, y_convexe, col = "red", lwd = 2)
  }
  abline(v = knots, col = "blue", lty = 2)

  # sans contraintes
  plot(xtab, ytab, pch = 16, cex = 0.5, col = "black",
       main = "without constraint ")
  lines(x_eval,y_sans,lwd=2)
  abline(v = knots, col = "blue", lty = 2)

}
#
