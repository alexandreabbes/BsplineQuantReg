# Demo: Logistic Curve Regression (He & Ng, 1999)
# This demo shows quantile regression on a logistic/sigmoid function
# Author: Alexandre Abbes

library(BsplineQuantReg)

cat("========================================\n")
cat("Demo: Logistic Curve Quantile Regression\n")
cat("========================================\n\n")

# Generate data from a logistic/sigmoid function
set.seed(42)
n_points <- 200
x <- seq(-5, 5, length.out = n_points)

# Logistic (sigmoid) function: P(x) = 1 / (1 + exp(-x))
true_logistic <- function(x) 1 / (1 + exp(-x))

# Add heteroscedastic noise (more noise in the middle)
noise <- rnorm(n_points, 0, 0.05) * (1 + 0.5 * exp(-(x^2)/4))
y <- true_logistic(x) + noise

# Knots (equally spaced quantiles)
kn <- 12
knots <- quantile(x, probs = seq(0, 1, length.out = kn + 1))

cat(sprintf("Number of points: %d\n", n_points))
cat(sprintf("Number of intervals: %d\n", kn))
cat(sprintf("Number of basis functions: %d\n\n", kn + 3))

# Fit multiple quantiles
tau_values <- c(0.1, 0.25, 0.5, 0.75, 0.9)
colors <- c("red", "orange", "blue", "green", "purple")
fits <- list()

for (i in seq_along(tau_values)) {
  tau <- tau_values[i]
  cat(sprintf("Fitting quantile tau = %.2f...\n", tau))
  fits[[i]] <- SplineConstQuantRegBs3(x, y, knots, tau = tau,
                                      monot = 0, convcons = 0)
}

# Evaluation
x_eval <- seq(-5, 5, length.out = 300)
y_true <- true_logistic(x_eval)

# Plot
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))

# Plot 1: Data and fits
plot(x, y, pch = 16, cex = 0.4, col = "gray",
     xlab = "x", ylab = "y",
     main = "Quantile Regression on Logistic Curve")
lines(x_eval, y_true, col = "black", lwd = 2, lty = 2)

for (i in seq_along(tau_values)) {
  y_fit <- spline_eval(fits[[i]], x_eval)
  lines(x_eval, y_fit, col = colors[i], lwd = 1.5)
}

legend("topleft", legend = c("True", paste("tau =", tau_values)),
       col = c("black", colors), lty = c(2, rep(1, 5)),
       lwd = 2, cex = 0.6)

# Plot 2: With monotonicity constraint (increasing)
cat("\n=== Fitting with monotonicity constraint ===\n")
fit_monot <- SplineConstQuantRegBs3(x, y, knots, tau = 0.5,
                                    monot = 1, convcons = 0)
y_monot <- spline_eval(fit_monot, x_eval)

plot(x, y, pch = 16, cex = 0.4, col = "gray",
     xlab = "x", ylab = "y",
     main = "Median Regression with Monotonicity")
lines(x_eval, y_true, col = "black", lwd = 2, lty = 2)
lines(x_eval, y_monot, col = "blue", lwd = 2)

# Check monotonicity
dy <- diff(y_monot) / diff(x_eval)
cat(sprintf("Minimum derivative (should be >= 0): %.6f\n", min(dy)))
if(min(dy) >= -1e-6) {
  cat(" Monotonicity constraint satisfied!\n")
} else {
  cat(" Monotonicity constraint violated\n")
}

legend("topleft", legend = c("True", "Monotonic fit (tau=0.5)"),
       col = c("black", "blue"), lty = c(2, 1), lwd = 2,cex=0.6)

cat("\nDemo completed.\n")
