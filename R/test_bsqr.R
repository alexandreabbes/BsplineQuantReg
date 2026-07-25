#' Comprehensive test function
#'
#' Runs quantile regression tests with and without constraints,
#' and displays results. Demo function.
#' This example also shows that solution may be evaluated (extrapolation)
#' out of the knot range
#' @param verbose boolean FALSE (default) or TRUE.
#' @param seed (default=NULL) value for the random generator.
#' @param degree 1, 2, 3 or 4 (default 3). The degree of the regression spline.
#' @export
test_karlin_simple<-function(degree=3,seed=NULL,verbose=FALSE)
{
  old <- par(mfrow = c(2,2))
  # Restore par
  on.exit(par(old))

if (!exists('verbose')){verbose=FALSE}
if (!exists('degree')){degree=3}
if (!exists('seed')){seed=NULL}

if (verbose) {message("Simple Test. degree=:" , degree)}
  # Store par



  if(!is.null(seed))
  {    set.seed(seed)}

  n_points <- 50
  xtab=((0:(n_points))/(n_points))

  # simple pscillating data
  #ytab<- -3 * xtab +sin(3*2*xtab*3.14)+ 0.2 * rnorm(n_points+1)
  #ytab= <- 2* xtab + 0.5 * sin(6 * pi * xtab) + 0.1 * rnorm(n_points+1)
  ytab<-xtab*(1-xtab)+0.1 * rnorm(n_points+1)
  kn <-8

  n=7

  monot=c(rep(1,n),rep(0,(12-n)))
  knot <-((3:(kn+3))/(kn+5) )
  res<-quantile_spline(xtab=xtab, ytab=ytab, knot, tau = 0.5,
                              monot = 0, degree=degree,verbose=verbose)
  if (verbose) {
    message("\n===  PARTIAL INCREASING TEST ===\n
            Monotone constraints list =", monot)
  }

  res_croissant <- quantile_spline(xtab=xtab, ytab=ytab, knot=knot, tau = 0.5,degree=degree,
                                          monot = monot, convcons=0, verbose=verbose)

  if (verbose) {
    message("\n=== DECREASING TEST ===\n")}

  res_decroissant <- quantile_spline(degree=degree,xtab=xtab, ytab=ytab, knot, tau = 0.5,
                                            monot = -1, verbose=verbose)
  if (verbose) {
    message("\n=== CONVEX  TEST ===\n")}
  res_convexe <- quantile_spline(xtab=xtab, ytab=ytab, knot,degree=degree,monot=0,convcons=-1, tau = 0.5, verbose=verbose)

der3cons=rep(1,length(knot))
der3cons[1:7]<-0
if (verbose){message("3rd order constraints list:",der3cons)}
    res_3rd <- quantile_spline(xtab, ytab, knot,degree=degree,monot=0,convcons=0,der3cons=der3cons, tau = 0.5, verbose=verbose)

  # Visualisation
  par(mfrow = c(3, 2))
  x_eval <- seq(0, 1, length.out = 200)

  y_sans=spline_eval(res,x_eval)
  y_croiss=spline_eval(res_croissant,x_eval)
  y_decroiss=spline_eval(res_decroissant,x_eval)
  y_convexe=spline_eval(res_convexe,x_eval)
  y_3rd=spline_eval(res_3rd,x_eval)

  plot(xtab, ytab, pch = 16, cex = 0.5, col = "black",
       main = c("Mixed mootone :\n increasing constraints
                \n on the first inter-knots"))
  lines(x_eval,y_croiss,lwd=2)

  #  view_spline
  #lines(x_eval, y_croiss, col = "red", lwd = 2)

  abline(v = knot, col = "blue", lty = 2)

  # Decroissant
  plot(xtab, ytab, pch = 16, cex = 0.5, col = "black",
       main = "decreasing constraint \n (inside the knots range)")

  lines(x_eval, y_decroiss, col = "red", lwd = 2)

  abline(v = knot, col = "blue", lty = 2)


  # Convexe
  plot(xtab, ytab, pch = 16, cex = 0.5, col = "black",
       main = "convex constraint\n inside the knots range")
  if (!is.null(res_convexe)) {
    lines(x_eval, y_convexe, col = "red", lwd = 2)
  }
  abline(v = knot, col = "blue", lty = 2)

  # sans contraintes
  plot(xtab, ytab, pch = 16, cex = 0.5, col = "black",
       main = "without constraint ")
  lines(x_eval,y_sans,lwd=2)
  abline(v = knot, col = "blue", lty = 2)

#der3>0
  plot(xtab, ytab, pch = 16, cex = 0.5, col = "black",
      main = "3rd derivative >0 \n on the second half \n of the knots range")
  lines(x_eval,y_3rd,lwd=2)
  abline(v = knot, col = "blue", lty = 2)

  (par(old))
}
