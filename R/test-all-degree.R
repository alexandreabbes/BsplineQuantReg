# Test: Compare all spline degrees (1, 2, 3, 4)

#library(BsplineQuantReg)

test_all_degrees <- function() {
  set.seed(42)
  n <- 150
  x <- seq(0, 1, length.out = n)

  # Test function with different characteristics
  y <- 2 * x^4 + 0.5 * sin(4 * pi * x) + 0.05 * rnorm(n)

  kn <- 8
  knots <- quantile(x, probs = seq(0, 1, length.out = kn + 1))

  cat("=== Comparing Splines of Different Degrees ===\n\n")

  # Fit all degrees with increasing constraint
  cat("Fitting with increasing constraint...\n")

  fit1 <- SplineConstQuantRegBs1(x, y, knots, tau = 0.5,
                                 monot = 1, verbose = FALSE)
  cat("  Degree 1 (linear)... done\n")

  fit2 <- SplineConstQuantRegBs2(x, y, knots, tau = 0.5,
                                 monot = 1, verbose = FALSE)
  cat("  Degree 2 (quadratic)... done\n")

  fit3 <- SplineConstQuantRegBs3(x, y, knots, tau = 0.5,
                                 monot = 1, verbose = FALSE)
  cat("  Degree 3 (cubic)... done\n")

  fit4 <- SplineConstQuantRegBs4(x, y, knots, tau = 0.5,
                                 monot = 1, verbose = FALSE)
  cat("  Degree 4 (quartic)... done\n")

  # Evaluation
  x_eval <- seq(0, 1, length.out = 200)
  y1 <- spline_eval(fit1, x_eval)
  y2 <- spline_eval(fit2, x_eval)
  y3 <- spline_eval(fit3, x_eval)
  y4 <- spline_eval(fit4, x_eval)

  # Plot
  par(mfrow = c(2, 3), mar = c(4, 4, 4, 2))

  plot(x, y, pch = 16, cex = 0.4, col = "gray", main = "Data")

  plot(x, y, pch = 16, cex = 0.4, col = "gray", main = "Linear (deg 1)")
  lines(x_eval, y1, col = "blue", lwd = 2)
  abline(v = knots, col = "blue", lty = 2, lwd = 0.5)

  plot(x, y, pch = 16, cex = 0.4, col = "gray", main = "Quadratic (deg 2)")
  lines(x_eval, y2, col = "green", lwd = 2)
  abline(v = knots, col = "blue", lty = 2, lwd = 0.5)

  plot(x, y, pch = 16, cex = 0.4, col = "gray", main = "Cubic (deg 3)")
  lines(x_eval, y3, col = "red", lwd = 2)
  abline(v = knots, col = "blue", lty = 2, lwd = 0.5)

  plot(x, y, pch = 16, cex = 0.4, col = "gray", main = "Quartic (deg 4)")
  lines(x_eval, y4, col = "purple", lwd = 2)
  abline(v = knots, col = "blue", lty = 2, lwd = 0.5)

  cat("\nDemo completed.\n")
}

test_all_degrees()
