# Test: Quartic Spline Quantile Regression
# Compare with cubic version

library(BsplineQuantReg)

#test_quartic <- function() {
  set.seed(42)
  n <- 100
  x <- seq(0, 1, length.out = n)

  # Quartic function with noise (ideal for quartic splines)
  y <- 4 * x^4 - 3 * x^2 + 0.5 * sin(4 * pi * x) + 0.2 * rnorm(n)

  kn <- 6
  knot <- (0:(kn))/(kn)
  knot<-knot
  cat("=== Comparing Quartic vs Cubic Splines ===\n\n")

  # Quartic spline (unconstrained)
  cat("Fitting quartic spline (unconstrained)...\n")
  fit_quart <- SplineConstQuantRegBs4(x, y, knot, tau = 0.5, verbose = TRUE)

  # Quartic spline with monotonicity
  cat("\nFitting quartic spline (increasing)...\n")
  fit_quart_monot <- SplineConstQuantRegBs4(x, y, knot, tau = 0.5,
                                            monot = 1, verbose = TRUE)

  # Quartic spline with convexity
  cat("\nFitting quartic spline (convex)...\n")
  fit_quart_conv <- SplineConstQuantRegBs4(x, y, knot, tau = 0.5,
                                           convcons = 1, verbose = TRUE)

  # Quartic spline with third derivative constraint
  cat("\nFitting quartic spline (positive 3rd derivative)...\n")
  fit_quart_d3 <- SplineConstQuantRegBs4(x, y, knot, tau = 0.5,
                                         der3cons = 1, verbose = TRUE)

  # Cubic spline for comparison
  cat("\nFitting cubic spline (unconstrained)...\n")
  fit_cubic <- SplineConstQuantRegBs3(x, y, knot, tau = 0.5, verbose = FALSE)

  # Evaluation
  x_eval <- seq(0, 1, length.out = 200)
  y_quart <- spline_eval(fit_quart, x_eval)
  y_quart_monot <- spline_eval(fit_quart_monot, x_eval)
  y_quart_conv <- spline_eval(fit_quart_conv, x_eval)
  y_quart_d3 <- spline_eval(fit_quart_d3, x_eval)
  y_cubic <- spline_eval(fit_cubic, x_eval)

  # Plot
  par(mfrow = c(2, 3), mar = c(4, 4, 4, 2))

  plot(x, y, pch = 16, cex = 0.4, col = "black", main = "Quartic (unconstrained)")
  lines(x_eval, y_quart, col = "blue", lwd = 2)
  abline(v = knot, col = "blue", lty = 2, lwd = 0.5)

  plot(x, y, pch = 16, cex = 0.4, col = "black", main = "Quartic (increasing)")
  lines(x_eval, y_quart_monot, col = "darkgreen", lwd = 2)
  abline(v = knot, col = "blue", lty = 2, lwd = 0.5)

  plot(x, y, pch = 16, cex = 0.4, col = "black", main = "Quartic (convex)")
  lines(x_eval, y_quart_conv, col = "purple", lwd = 2)
  abline(v = knot, col = "blue", lty = 2, lwd = 0.5)

  plot(x, y, pch = 16, cex = 0.4, col = "black", main = "Quartic (3rd deriv >= 0)")
  lines(x_eval, y_quart_d3, col = "orange", lwd = 2)
  abline(v = knot, col = "blue", lty = 2, lwd = 0.5)

  plot(x, y, pch = 16, cex = 0.4, col = "black", main = "Cubic (unconstrained)")
  lines(x_eval, y_cubic, col = "red", lwd = 2)

#  return(list(quart = fit_quart, quart_monot = fit_quart_monot,
#              quart_conv = fit_quart_conv, quart_d3 = fit_quart_d3,
#              cubic = fit_cubic))
#}

#test_quartic()
