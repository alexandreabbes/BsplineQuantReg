# Demo: Comparison of Spline Degrees (1 to 4)
# Author: Alexandre Abbes
#
# This demo shows the effect of spline degree on quantile regression
# with shape constraints. It compares degrees 1 (linear) to 4 (quartic)
# on the same dataset with monotonicity and convexity constraints.



library(BsplineQuantReg)

cat("========================================\n")
cat("Demo: Comparison of Spline Degrees\n")
cat("========================================\n\n")

# Save original graphical parameters
oldpar <- par(mfrow = c(2, 2), mar = c(4, 4, 4, 2))

# Generate test data
set.seed(42)
n <- 200
x <- seq(0, 1, length.out = n)

# Function with increasing trend and curvature
true_function <- function(x) {
  2 * x + 0.5 * sin(3 * pi * x) + 0.3 * x^2
}

# Add noise
y <- true_function(x) + 0.7 * rnorm(n)

# Knots
knot <- quantile(x, probs = seq(0, 1, length.out = 10))

cat(sprintf("Number of points: %d\n", n))
cat(sprintf("Number of knot: %d\n", length(knot)))
cat(sprintf("Tau: 0.5 (median)\n\n"))

# Fit all degrees with increasing constraint
cat("Fitting models with increasing constraint...\n")

fit1 <- quantile_spline(x, y, knot, tau = 0.5, degree = 1,
                        monot = 1, verbose = FALSE, callable = TRUE)
cat("  Degree 1 (linear)... done\n")

fit2 <- quantile_spline(x, y, knot, tau = 0.5, degree = 2, monot = 1, verbose = FALSE, callable = TRUE)
cat("  Degree 2 (quadratic)... done\n")

fit3 <- quantile_spline(x, y, knot, tau = 0.5, degree = 3,
                        monot = 1, verbose = FALSE, callable = TRUE)
cat("  Degree 3 (cubic)... done\n")

fit4 <- quantile_spline(x, y, knot, tau = 0.5, degree = 4,
                        monot = 1, verbose = FALSE, callable = TRUE)
cat("  Degree 4 (quartic)... done\n")

# Fit all degrees with convexity constraint
cat("\nFitting models with convexity constraint...\n")

fit1_conv <- quantile_spline(x, y, knot, tau = 0.5, degree = 1,
                             convcons = 1, verbose = FALSE, callable = TRUE)
cat("  Degree 1 (linear)... done\n")

fit2_conv <- quantile_spline(x, y, knot, tau = 0.5, degree = 2,
                             convcons = 1, verbose = FALSE, callable = TRUE)
cat("  Degree 2 (quadratic)... done\n")

fit3_conv <- quantile_spline(x, y, knot, tau = 0.5, degree = 3,
                             convcons = 1, verbose = FALSE, callable = TRUE)
cat("  Degree 3 (cubic)... done\n")

fit4_conv <- quantile_spline(x, y, knot, tau = 0.5, degree = 4,
                             convcons = 1, verbose = FALSE, callable = TRUE)
cat("  Degree 4 (quartic)... done\n")

# Evaluation
x_eval <- seq(0, 1, length.out = 300)
y_true <- true_function(x_eval)

y1 <- fit1(x_eval)
y2 <- fit2(x_eval)
y3 <- fit3(x_eval)
y4 <- fit4(x_eval)

y1_conv <- fit1_conv(x_eval)
y2_conv <- fit2_conv(x_eval)
y3_conv <- fit3_conv(x_eval)
y4_conv <- fit4_conv(x_eval)

# Colors
colors <- c("blue", "darkgreen", "red", "purple")
degree_labels <- c("Linear (deg 1)", "Quadratic (deg 2)", "Cubic (deg 3)", "Quartic (deg 4)")

# ============ PLOT 1: Increasing constraint ============
plot(x, y, pch = 16, cex = 0.4, col = "gray",
     xlab = "x", ylab = "y",
     main = "Monotonicity Constraint (increasing)")
lines(x_eval, y_true, col = "black", lwd = 2, lty = 2)

for (i in 1:4) {
  y_fit <- list(y1, y2, y3, y4)[[i]]
  lines(x_eval, y_fit, col = colors[i], lwd = 2)
}

legend("topleft", legend = c("True", degree_labels),
       col = c("black", colors), lty = c(2, rep(1, 4)),
       lwd = 2, cex = 0.7)

abline(v = knot, col = "blue", lty = 3, lwd = 0.5)

# ============ PLOT 2: Convexity constraint ============
plot(x, y, pch = 16, cex = 0.4, col = "gray",
     xlab = "x", ylab = "y",
     main = "Convexity Constraint")
lines(x_eval, y_true, col = "black", lwd = 2, lty = 2)

for (i in 1:4) {
  y_fit <- list(y1_conv, y2_conv, y3_conv, y4_conv)[[i]]
  lines(x_eval, y_fit, col = colors[i], lwd = 2)
}

legend("topleft", legend = c("True", degree_labels),
       col = c("black", colors), lty = c(2, rep(1, 4)),
       lwd = 2, cex = 0.7)

abline(v = knot, col = "blue", lty = 3, lwd = 0.5)

# ============ PLOT 3: Comparison of degrees (unconstrained) ============
fit1_uncon <- quantile_spline(x, y, knot, tau = 0.5, degree = 1,
                              verbose = FALSE, callable = TRUE)
fit2_uncon <- quantile_spline(x, y, knot, tau = 0.5, degree = 2,
                              verbose = FALSE, callable = TRUE)
fit3_uncon <- quantile_spline(x, y, knot, tau = 0.5, degree = 3,
                              verbose = FALSE, callable = TRUE)
fit4_uncon <- quantile_spline(x, y, knot, tau = 0.5, degree = 4,
                              verbose = FALSE, callable = TRUE)

y1_uncon <- fit1_uncon(x_eval)
y2_uncon <- fit2_uncon(x_eval)
y3_uncon <- fit3_uncon(x_eval)
y4_uncon <- fit4_uncon(x_eval)

plot(x, y, pch = 16, cex = 0.4, col = "gray",
     xlab = "x", ylab = "y",
     main = "Unconstrained")
lines(x_eval, y_true, col = "black", lwd = 2, lty = 2)

for (i in 1:4) {
  y_fit <- list(y1_uncon, y2_uncon, y3_uncon, y4_uncon)[[i]]
  lines(x_eval, y_fit, col = colors[i], lwd = 2)
}

legend("topleft", legend = c("True", degree_labels),
       col = c("black", colors), lty = c(2, rep(1, 4)),
       lwd = 2, cex = 0.7)

abline(v = knot, col = "blue", lty = 3, lwd = 0.5)

# ============ PLOT 4: Convergence with degree (increasing constraint) ============
# Show all four fits together
plot(x, y, pch = 16, cex = 0.4, col = "gray",
     xlab = "x", ylab = "y",
     main = "All Degrees with Monotonicity Constraint")
lines(x_eval, y_true, col = "black", lwd = 2, lty = 2)

lines(x_eval, y1, col = colors[1], lwd = 2, lty = 2)
lines(x_eval, y2, col = colors[2], lwd = 2, lty = 3)
lines(x_eval, y3, col = colors[3], lwd = 2, lty = 1)
lines(x_eval, y4, col = colors[4], lwd = 2, lty = 4)

legend("topleft", legend = c("True", degree_labels),
       col = c("black", colors), lty = c(2, 2, 3, 1, 4),
       lwd = 2, cex = 0.6)

abline(v = knot, col = "blue", lty = 3, lwd = 0.5)

# ============ SUMMARY ============
cat("\n========================================\n")
cat("Summary\n")
cat("========================================\n")
cat("Observations:\n")
cat("  - Linear splines (deg 1) are piecewise linear\n")
cat("  - Quadratic splines (deg 2) have continuous derivative\n")
cat("  - Cubic splines (deg 3) have continuous second derivative\n")
cat("  - Quartic splines (deg 4) have continuous third derivative\n")
cat("\nHigher degrees provide more smoothness but can overfit.\n")
cat("The Karlin-Studden constraints work for all degrees.\n")

# Restore graphical parameters
par(oldpar)

cat("\nDemo completed.\n")
