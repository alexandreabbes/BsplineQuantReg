# tests/test-quantile-spline.R

test_quantile_spline <- function() {
  set.seed(42)
  n <- 30
  x <- seq(0, 1, length.out = n)
  y <- 2*x^3 + 0.5*sin(4*pi*x) + 0.45*rnorm(n)
  knot <- quantile(x, probs = seq(0, 1, length.out = 8))

  cat("=== Testing unified quantile_spline function ===\n\n")

  # Test all degrees
  degrees <- 1:4
  fits <- list()

  for (d in degrees) {
    cat("Degree", d, "...")
    fits[[d]] <- quantile_spline(x, y, knot, tau = 0.5, degree = d,
                                 monot = 1, verbose = FALSE, callable=TRUE)
    cat(" done\n")
  }

  # Evaluation and plot
  x_eval <- seq(0, 1, length.out = 200)
  colors <- c("blue", "green", "red", "purple")
  labels <- c("Linear (deg 1)", "Quadratic (deg 2)",
              "Cubic (deg 3)", "Quartic (deg 4)")

  par(mfrow = c(2, 2), mar = c(4, 4, 4, 2))

  for (d in degrees) {
    #y_fit <- spline_eval(fits[[d]], x_eval)
    y_fit <- fits[[d]](x_eval)
    plot(x, y, pch = 16, cex = 0.4, col = "black",
         main = labels[d], xlab = "x", ylab = "y")
    lines(x_eval, y_fit, col = colors[d], lwd = 2)
    abline(v = knot, col = "blue", lty = 2, lwd = 0.5)
  }

  cat("\nDemo completed.\n")
}

test_quantile_spline()
