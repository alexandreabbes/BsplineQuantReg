# Demo: Logistic Curve Regression (He & Ng, 1999)
# This demo shows quantile regression on a logistic/sigmoid function
# Author: Alexandre Abbes
#
# Constraints demonstrated:
# 1. Multiple quantiles (0.1, 0.25, 0.5, 0.75, 0.9)
# 2. Monotonicity constraint (increasing)
# 3. Convexity constraint (convex on left, concave on right)
# 4. Combined monotonicity + convexity

library(BsplineQuantReg)

cat("========================================\n")
cat("Demo: Logistic Curve Quantile Regression\n")
cat("========================================\n\n")

# Save original graphical parameters
oldpar <- par(mfrow = c(2, 2), mar = c(4, 4, 4, 2))

# Generate data from a logistic/sigmoid function
#set.seed(42)
n_points <- 200
x <- seq(-5, 5, length.out = n_points)

# Logistic (sigmoid) function: P(x) = 1 / (1 + exp(-x))
true_logistic <- function(x) 1 / (1 + exp(-x))

# Add heteroscedastic noise (more noise in the middle)
noise <- rnorm(n_points, 0, 0.05) * (1 + 0.5 * exp(-(x^2)/4))
y <- true_logistic(x) + noise

# Knots (equally spaced quantiles)
kn <- 6
knots <- quantile(x, probs = seq(0, 1, length.out = kn + 1))

# Define convexity constraints:
# Convex (second derivative >= 0) on left half (x < 0)
# Concave (second derivative <= 0) on right half (x > 0)
convex <- rep(0, kn + 1)
for (i in 1:(kn + 1)) {
  if (knots[i] < 0) {
    convex[i] <- 1      # convex on left (x < 0)
  } else if (knots[i] > 0) {
    convex[i] <- -1     # concave on right (x > 0)
  } else {
    convex[i] <- 0      # at x = 0, no constraint
  }
}

cat(sprintf("Number of points: %d\n", n_points))
cat(sprintf("Number of intervals: %d\n", kn))
cat(sprintf("Number of basis functions: %d\n\n", kn + 3))

cat("Knots positions:\n")
for (i in 1:(kn + 1)) {
  cat(sprintf("  Knot %d: x = %.2f", i, knots[i]))
  if (knots[i] < 0) {
    cat(" (convex region)\n")
  } else if (knots[i] > 0) {
    cat(" (concave region)\n")
  } else {
    cat(" (inflection point)\n")
  }
}
cat("\n")

# Fit multiple quantiles
tau_values <- c(0.1, 0.25, 0.5, 0.75, 0.9)
colors <- c("red", "orange", "blue", "green", "purple")
fits <- list()

cat("Fitting multiple quantiles...\n")
for (i in seq_along(tau_values)) {
  tau <- tau_values[i]
  cat(sprintf("  tau = %.2f...\n", tau))
  fits[[i]] <- SplineConstQuantRegBs3(x, y, knots, tau = tau,
                                      monot = 0, convcons = 0)
}

# Evaluation grid
x_eval <- seq(-5, 5, length.out = 300)
y_true <- true_logistic(x_eval)

# ============================================================
# Plot 1: Multiple quantiles
# ============================================================
plot(x, y, pch = 16, cex = 0.4, col = "black",
     xlab = "x", ylab = "y",
     main = "1. Quantile Regression\nMultiple Quantiles")
lines(x_eval, y_true, col = "black", lwd = 2, lty = 2)

for (i in seq_along(tau_values)) {
  y_fit <- spline_eval(fits[[i]], x_eval)
  lines(x_eval, y_fit, col = colors[i], lwd = 1.5)
}

# Add knots
abline(v = knots, col = "blue", lty = 2, lwd = 0.5)
legend("topleft", legend = c("True", paste("tau =", tau_values)),
       col = c("black", colors), lty = c(2, rep(1, 5)),
       lwd = 1.5, cex = 0.6)

# ============================================================
# Plot 2: Monotonicity constraint (increasing)
# ============================================================
cat("\n=== Fitting with monotonicity constraint ===\n")
fit_monot <- SplineConstQuantRegBs3(x, y, knots, tau = 0.5,
                                    monot = 1, convcons = 0)
y_monot <- spline_eval(fit_monot, x_eval)

plot(x, y, pch = 16, cex = 0.4, col = "black",
     xlab = "x", ylab = "y",
     main = "2. Median Regression\nwith Monotonicity Constraint")
lines(x_eval, y_true, col = "black", lwd = 2, lty = 2)
lines(x_eval, y_monot, col = "blue", lwd = 2)
abline(v = knots, col = "blue", lty = 2, lwd = 0.5)

# Check monotonicity
dy <- diff(y_monot) / diff(x_eval)
cat(sprintf("  Minimum derivative: %.6f (should be >= 0)\n", min(dy)))
if (min(dy) >= -1e-6) {
  cat("   Monotonicity constraint satisfied!\n")
} else {
  cat("   Monotonicity constraint violated\n")
}

legend("topleft", legend = c("True", "Monotonic (tau=0.5)"),
       col = c("black", "blue"), lty = c(2, 1), lwd = 2, cex = 0.6)

# ============================================================
# Plot 3: Convexity constraint (convex left, concave right)
# ============================================================
cat("\n=== Fitting with convexity/concavity constraint ===\n")
fit_conv <- SplineConstQuantRegBs3(x, y, knots, tau = 0.5,
                                   monot = 0, convcons = convex)
y_conv <- spline_eval(fit_conv, x_eval)

plot(x, y, pch = 16, cex = 0.4, col = "black",
     xlab = "x", ylab = "y",
     main = "3. Median Regression\nwith Convex/Concave Constraint")
lines(x_eval, y_true, col = "black", lwd = 2, lty = 2)
lines(x_eval, y_conv, col = "darkgreen", lwd = 2)
abline(v = knots, col = "blue", lty = 2, lwd = 0.5)
abline(v = 0, col = "red", lty = 2, lwd = 2)

# Add shading for convex/concave regions
rect(-5, -0.2, 0, 1.2, col = rgb(0, 1, 0, 0.1), border = NA)
rect(0, -0.2, 5, 1.2, col = rgb(1, 0, 0, 0.1), border = NA)
text(-2.5, 1.0, "Convex\n(2nd deriv >= 0)", col = "darkgreen", cex = 0.7)
text(2.5, 1.0, "Concave\n(2nd deriv <= 0)", col = "darkred", cex = 0.7)

legend("topleft", legend = c("True", "Convex/Concave (tau=0.5)"),
       col = c("black", "darkgreen"), lty = c(2, 1), lwd = 2, cex = 0.6)

# ============================================================
# Plot 4: Combined monotonicity + convexity
# ============================================================
cat("\n=== Fitting with monotonicity + convexity ===\n")
fit_both <- SplineConstQuantRegBs3(x, y, knots, tau = 0.5,
                                   monot = 1, convcons = convex)
y_both <- spline_eval(fit_both, x_eval)

plot(x, y, pch = 16, cex = 0.4, col = "black",
     xlab = "x", ylab = "y",
     main = "4. Median Regression\nMonotonicity + Convex/Concave")
lines(x_eval, y_true, col = "black", lwd = 2, lty = 2)
lines(x_eval, y_both, col = "purple", lwd = 2)
abline(v = knots, col = "blue", lty = 2, lwd = 0.5)
abline(v = 0, col = "red", lty = 2, lwd = 2)

# Add shading
rect(-5, -0.2, 0, 1.2, col = rgb(0, 1, 0, 0.1), border = NA)
rect(0, -0.2, 5, 1.2, col = rgb(1, 0, 0, 0.1), border = NA)

legend("topleft", legend = c("True", "Both constraints (tau=0.5)"),
       col = c("black", "purple"), lty = c(2, 1), lwd = 2, cex = 0.6)

# ============================================================
# Summary
# ============================================================
cat("\n========================================\n")
cat("Summary\n")
cat("========================================\n")
cat("The logistic curve is naturally increasing and S-shaped:\n")
cat("  - Left side (x < 0): convex (accelerating growth)\n")
cat("  - Right side (x > 0): concave (decelerating growth)\n")
cat("  - Inflection point at x = 0\n")
cat("\nConstraint configurations:\n")
cat("  Plot 1: Multiple quantiles (0.1 to 0.9)\n")
cat("  Plot 2: Monotonicity only\n")
cat("  Plot 3: Convexity only (convex left, concave right)\n")
cat("  Plot 4: Both monotonicity + convexity\n")
cat("\nDemo completed.\n")

# Restore graphical parameters
par(oldpar)
