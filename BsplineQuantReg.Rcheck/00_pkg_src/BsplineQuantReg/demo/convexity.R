# Demo: Convexity and Concavity Constraints
# Author: Alexandre Abbes
#
# This demo shows:
# 1. Convexity constraint (second derivative >= 0)
# 2. Concavity constraint (second derivative <= 0)
# 3. Partial convexity (only on the right half)

library(BsplineQuantReg)

cat("========================================\n")
cat("Demo: Convexity & Concavity Constraints\n")
cat("========================================\n\n")

set.seed(42)
n_points <- 150
xtab <- seq(-2, 2, length.out = n_points)

# Function: x^3 - 3x (convex for x > 0, concave for x < 0)
true_function <- function(x) x^3 - 3*x

# Add noise
ytab <- true_function(xtab) + 0.2 * rnorm(n_points)

# Knots
kn <- 12
knots <- quantile(xtab, probs = seq(0, 1, length.out = kn + 1))

cat("Data: cubic function x^3 - 3x with noise\n")
cat(sprintf("Number of points: %d\n", n_points))
cat(sprintf("Number of intervals: %d\n\n", kn))

# Fit models
cat("Fitting models...\n")
cat("  - Unconstrained...")
fit_uncon <- SplineConstQuantRegBs3(xtab, ytab, knots, tau = 0.5,
                                    monot = 0, convcons = 0)
cat(" done\n")

cat("  - Convexity constraint (everywhere)...")
fit_convex <- SplineConstQuantRegBs3(xtab, ytab, knots, tau = 0.5,
                                     monot = 0, convcons = 1)
cat(" done\n")

cat("  - Concavity constraint (everywhere)...")
fit_concave <- SplineConstQuantRegBs3(xtab, ytab, knots, tau = 0.5,
                                      monot = 0, convcons = -1)
cat(" done\n")

# Partial convexity (only on right half, x >= 0)
convcons_partial <- rep(0, kn + 1)
for (i in 1:(kn + 1)) {
  if (knots[i] >= 0) {
    convcons_partial[i] <- 1
  }
}
cat("  - Partial convexity (x >= 0 only)...")
fit_partial <- SplineConstQuantRegBs3(xtab, ytab, knots, tau = 0.5,
                                      monot = 0, convcons = convcons_partial)
cat(" done\n\n")

# Evaluation
x_eval <- seq(-2, 2, length.out = 300)
y_true <- true_function(x_eval)
y_uncon <- spline_eval(fit_uncon, x_eval)
y_convex <- spline_eval(fit_convex, x_eval)
y_concave <- spline_eval(fit_concave, x_eval)
y_partial <- spline_eval(fit_partial, x_eval)

# Create plots - 2x2 layout with clean labels
par(mfrow = c(2, 2), mar = c(4, 4, 3, 2))

# Plot 1: Unconstrained
plot(xtab, ytab, pch = 16, cex = 0.4, col = "lightgray",
     xlab = "x", ylab = "y", main = "Unconstrained Fit")
lines(x_eval, y_true, col = "black", lwd = 2, lty = 2)
lines(x_eval, y_uncon, col = "red", lwd = 2)
legend("bottomleft", legend = c("True", "Unconstrained"),
       col = c("black", "red"), lty = c(2, 1), lwd = 2, cex = 0.6)

# Plot 2: Convexity everywhere
plot(xtab, ytab, pch = 16, cex = 0.4, col = "lightgray",
     xlab = "x", ylab = "y", main = "Convexity Constraint (everywhere)")
lines(x_eval, y_true, col = "black", lwd = 2, lty = 2)
lines(x_eval, y_convex, col = "blue", lwd = 2)
legend("bottomleft", legend = c("True", "Convex"),
       col = c("black", "blue"), lty = c(2, 1), lwd = 2, cex = 0.7)

# Plot 3: Concavity everywhere
plot(xtab, ytab, pch = 16, cex = 0.4, col = "lightgray",
     xlab = "x", ylab = "y", main = "Concavity Constraint (everywhere)")
lines(x_eval, y_true, col = "black", lwd = 2, lty = 2)
lines(x_eval, y_concave, col = "darkgreen", lwd = 2)
legend("bottomleft", legend = c("True", "Concave"),
       col = c("black", "darkgreen"), lty = c(2, 1), lwd = 2, cex = 0.7)

# Plot 4: Partial convexity (x >= 0 only)
plot(xtab, ytab, pch = 16, cex = 0.4, col = "lightgray",
     xlab = "x", ylab = "y", main = "Partial Convexity (x >= 0 only)")
lines(x_eval, y_true, col = "black", lwd = 2, lty = 2)
lines(x_eval, y_partial, col = "purple", lwd = 2)
abline(v = 0, col = "red", lty = 2, lwd = 1.5)
text(0.15, -1.5, "Convex region", col = "red", cex = 0.7)
text(-0.8, 1.5, "Unconstrained", col = "red", cex = 0.7)
legend("bottomleft", legend = c("True", "Partial convex"),
       col = c("black", "purple"), lty = c(2, 1), lwd = 2, cex = 0.7)

cat("========================================\n")
cat("Demo completed.\n")
cat("========================================\n")
