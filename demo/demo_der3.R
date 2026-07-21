# Demo: Third Derivative Constraints
# Author: Alexandre Abbes
#
# This test checks the third derivative constraints for cubic (degree 3)
# and quartic (degree 4) splines on a highly perturbated curve.
#
# Usage:
#   - Set degree = 3 for cubic, degree = 4 for quartic
#   - Or run with: degree <- 3; source("tests/test_der3.R")

oldpar <- par(mfrow = c(2,2))


library(BsplineQuantReg)

cat("==============================================\n")
cat("Demo: Third Derivative Constraints\n")
cat("==============================================\n\n")


# ============================================================
# 1. Set degree (default = 3)
# ============================================================

if (!exists("degree")) {
  degree <- 3
}

cat("========================================\n")
cat("Test: Third Derivative Constraints\n")
cat(sprintf("Degree: %d (%s)\n", degree, ifelse(degree == 3, "Cubic", "Quartic")))
cat("========================================\n\n")

# ============================================================
# 2. Generate highly perturbated data
# ============================================================

set.seed(123)
n <- 150
x <- seq(0, 1, length.out = n)

# Function with multiple changes in curvature
true_function <- function(x) {
  3 * x^3 - 2 * x^2 + 0.5 * sin(8 * pi * x) + 0.2 * cos(12 * pi * x)
}

# Add significant noise
noise_sd <- 0.3
y <- true_function(x) + noise_sd * rnorm(n)

# Knots (more knots for flexibility)
kn <- 12
knots <- quantile(x, probs = seq(0, 1, length.out = kn + 1))

cat(sprintf("Number of points: %d\n", n))
cat(sprintf("Number of knots: %d\n", length(knots)))
cat(sprintf("Noise standard deviation: %.2f\n\n", noise_sd))

# ============================================================
# 3. Define third derivative constraints
# ============================================================

# Third derivative = 0 means cubic polynomial (constant third derivative for cubic)
# Third derivative > 0 means third derivative positive
# Third derivative < 0 means third derivative negative

cat("Third derivative constraint configurations:\n")
cat("  1. Positive (increasing curvature)\n")
cat("  2. Negative (decreasing curvature)\n")
cat("  3. Mixed (positive then negative)\n\n")

# ============================================================
# 4. Fit models
# ============================================================

cat(sprintf("Fitting models (degree = %d)...\n\n", degree))

# --- Unconstrained ---
cat("  Unconstrained...")
fit_uncon <- quantile_spline(x, y, knots, tau = 0.5, degree = degree,
                             monot = 0, convcons = 0, der3cons = 0,
                             callable = TRUE, verbose = FALSE)
cat(" done\n")

# --- Positive third derivative ---
if (degree == 3) {
  der3_positive <- rep(1, kn)
  der3_negative <- rep(-1, kn)
  der3_mixed <- c(rep(1, floor(kn/2)), rep(-1, kn - floor(kn/2)))
} else {
  # degree == 4
  der3_positive <- rep(1, kn + 1)
  der3_negative <- rep(-1, kn + 1)
  der3_mixed <- c(rep(1, floor((kn+1)/2)), rep(-1, (kn+1) - floor((kn+1)/2)))
}

cat("  Positive third derivative...")
fit_pos <- quantile_spline(x, y, knots, tau = 0.5, degree = degree,
                           monot = 0, convcons = 0, der3cons = der3_positive,
                           callable = TRUE, verbose = FALSE)
cat(" done\n")

# --- Negative third derivative ---
cat("  Negative third derivative...")
fit_neg <- quantile_spline(x, y, knots, tau = 0.5, degree = degree,
                           monot = 0, convcons = 0, der3cons = der3_negative,
                           callable = TRUE, verbose = FALSE)
cat(" done\n")

# --- Mixed third derivative ---
cat("  Mixed third derivative...")
fit_mixed <- quantile_spline(x, y, knots, tau = 0.5, degree = degree,
                             monot = 0, convcons = 0, der3cons = der3_mixed,
                             callable = TRUE, verbose = FALSE)
cat(" done\n\n")

# ============================================================
# 5. Evaluate fits
# ============================================================

x_eval <- seq(0, 1, length.out = 300)
y_true <- true_function(x_eval)

y_uncon <- fit_uncon(x_eval)
y_pos <- fit_pos(x_eval)
y_neg <- fit_neg(x_eval)
y_mixed <- fit_mixed(x_eval)

# ============================================================
# 6. Compute and verify third derivatives
# ============================================================

# Function to compute third derivative numerically
compute_third_deriv <- function(y, x) {
  h <- diff(x)[1]
  n <- length(y)
  # Central differences for interior points
  d3 <- numeric(n)
  for (i in 3:(n-2)) {
    d3[i] <- (-y[i-2] + 2*y[i-1] - 2*y[i+1] + y[i+2]) / (2 * h^3)
  }
  # Boundary approximations
  d3[1] <- d3[3]
  d3[2] <- d3[3]
  d3[n-1] <- d3[n-2]
  d3[n] <- d3[n-2]
  return(d3)
}

d3_pos <- compute_third_deriv(y_pos, x_eval)
d3_neg <- compute_third_deriv(y_neg, x_eval)
d3_mixed <- compute_third_deriv(y_mixed, x_eval)

# Verification
cat("=== Verification of Third Derivative Constraints ===\n")
cat(sprintf("Degree %d (%s):\n", degree, ifelse(degree == 3, "Cubic", "Quartic")))
cat(sprintf("  Positive:   min d3 = %.4f (should be >= 0)\n", min(d3_pos, na.rm = TRUE)))
cat(sprintf("  Negative:   max d3 = %.4f (should be <= 0)\n", max(d3_neg, na.rm = TRUE)))
cat(sprintf("  Mixed:      min d3 on first half = %.4f, max on second half = %.4f\n",
            min(d3_mixed[1:150], na.rm = TRUE), max(d3_mixed[151:300], na.rm = TRUE)))

# Check if constraints are satisfied
cat("\nConstraint satisfaction:\n")
cat(sprintf("  Positive: %s\n", ifelse(min(d3_pos, na.rm = TRUE) >= -1e-4, "OK", "FAILED")))
cat(sprintf("  Negative: %s\n", ifelse(max(d3_neg, na.rm = TRUE) <= 1e-4, "OK", "FAILED")))
cat(sprintf("  Mixed:    %s (first half >= 0, second half <= 0)\n",
            ifelse(min(d3_mixed[1:150], na.rm = TRUE) >= -1e-4 &&
                     max(d3_mixed[151:300], na.rm = TRUE) <= 1e-4, "OK", "FAILED")))

# ============================================================
# 7. Visualization
# ============================================================

# Save original par
oldpar <- par(mfrow = c(2, 2), mar = c(4, 4, 4, 2))

# Plot 1: All fits
plot(x, y, pch = 16, cex = 0.4, col = "gray",
     xlab = "x", ylab = "y",
     main = sprintf("Degree %d - Third Derivative Constraints", degree))
lines(x_eval, y_true, col = "black", lwd = 2, lty = 2)
lines(x_eval, y_uncon, col = "red", lwd = 1.5)
lines(x_eval, y_pos, col = "blue", lwd = 1.5)
lines(x_eval, y_neg, col = "darkgreen", lwd = 1.5)
lines(x_eval, y_mixed, col = "purple", lwd = 1.5)
abline(v = knots, col = "blue", lty = 3, lwd = 0.5)
legend("topleft", legend = c("True", "Unconstrained", "Positive", "Negative", "Mixed"),
       col = c("black", "red", "blue", "darkgreen", "purple"),
       lty = c(2, 1, 1, 1, 1), lwd = 2, cex = 0.7)

# Plot 2: Third derivatives - Positive constraint
plot(x_eval, d3_pos, type = "l", col = "blue", lwd = 2,
     xlab = "x", ylab = "Third Derivative",
     main = sprintf("Degree %d - Third Derivatives (Positive)", degree))
abline(h = 0, col = "red", lty = 2)
grid()

# Plot 3: Third derivatives - Negative constraint
plot(x_eval, d3_neg, type = "l", col = "darkgreen", lwd = 2,
     xlab = "x", ylab = "Third Derivative",
     main = sprintf("Degree %d - Third Derivatives (Negative)", degree))
abline(h = 0, col = "red", lty = 2)
grid()

# Plot 4: Third derivatives - Mixed constraint
plot(x_eval, d3_mixed, type = "l", col = "purple", lwd = 2,
     xlab = "x", ylab = "Third Derivative",
     main = sprintf("Degree %d - Third Derivatives (Mixed)", degree))
abline(h = 0, col = "red", lty = 2)
abline(v = 0.5, col = "orange", lty = 3, lwd = 2)
text(0.25, max(d3_mixed, na.rm = TRUE) * 0.8, "Positive", col = "blue")
text(0.75, min(d3_mixed, na.rm = TRUE) * 0.8, "Negative", col = "darkgreen")
grid()

# Restore par
par(oldpar)

# ============================================================
# 8. Summary
# ============================================================

cat("\n========================================\n")
cat("Summary\n")
cat("========================================\n")
cat(sprintf("Degree %d (%s) third derivative constraints:\n",
            degree, ifelse(degree == 3, "Cubic", "Quartic")))
if (degree == 3) {
  cat("  - der3cons is applied per interval (length = kn)\n")
} else {
  cat("  - der3cons is applied per knot (length = kn + 1)\n")
}
cat("\nThe constraints force the third derivative to be positive, negative,\n")
cat("or mixed, controlling the curvature changes in the fitted function.\n")

par(oldpar)

cat("\nDemo completed.\n")
