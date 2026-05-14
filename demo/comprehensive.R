# Demo: Comprehensive Test of Quantile Regression Constraints
# Author: Alexandre Abbes
#
# This demo runs a complete test of all constraint types:
# - Partial monotonicity (increasing on first intervals only)
# - Full monotonicity (decreasing everywhere)
# - Convexity constraint
# - Unconstrained fit
#
# The results are displayed in a 2x2 panel for comparison.

library(BsplineQuantReg)

cat("==============================================\n")
cat("Demo: Comprehensive Quantile Regression Tests\n")
cat("==============================================\n\n")

# Generate synthetic data
set.seed(42)
n_points <- 100
xtab <- seq(0, 1, length.out = n_points + 1)  # 51 points from 0 to 1

# Data: increasing trend with oscillation + noise
# y = 2x + 0.5*sin(6 pi x) + noise
ytab <- 2 * xtab + 0.5 * sin(6 * pi * xtab) + 0.2 * rnorm(n_points + 1)

# Number of intervals
kn <- 12
cat(sprintf("Number of observations: %d\n", length(xtab)))
cat(sprintf("Number of intervals (kn): %d\n", kn))
cat(sprintf("Number of basis functions: %d\n\n", kn + 3))

# Define knots (quantile-based)
knots <- quantile(xtab, probs = seq(0, 1, length.out = kn + 1))

# Define partial monotonicity: increasing only on first 7 intervals
n_constrained <- 7
monot_partial <- c(rep(1, n_constrained), rep(0, kn - n_constrained))

cat("Constraint configurations:\n")
cat("-------------------------\n")
cat("1. Unconstrained:        no constraints\n")
cat("2. Partial increasing:   increasing on intervals 1-7 only\n")
cat("3. Full decreasing:      decreasing everywhere\n")
cat("4. Convexity:            second derivative >= 0 everywhere\n\n")

cat("Fitting models...\n")

# Fit 1: Unconstrained (tau = 0.9, 0.5, 0.1)
cat("  - Unconstrained (tau = 0.9, 0.5, 0.1)")
res_uncon1 <- SplineConstQuantRegBs3(xtab, ytab, knots, tau = 0.9,
                                    monot = 0, convcons = 0, solver = "OSQP")
res_uncon2 <- SplineConstQuantRegBs3(xtab, ytab, knots, tau = 0.1,
                                    monot = 0, convcons = 0, solver = "OSQP")
res_uncon3 <- SplineConstQuantRegBs3(xtab, ytab, knots, tau = 0.5,
                                    monot = 0, convcons = 0, solver = "OSQP")

cat(" done\n")

# Fit 2: Partial increasing (tau = 0.5)
cat("  - Partial increasing (tau = 0.5, 0.1, intervals 1-7)...")
res_croissant1 <- SplineConstQuantRegBs3(xtab, ytab, knots, tau = 0.5,
                                        monot = monot_partial, convcons = 0,
                                        solver = "OSQP")
res_croissant2 <- SplineConstQuantRegBs3(xtab, ytab, knots, tau = 0.1,
                                         monot = monot_partial, convcons = 0,
                                         solver = "OSQP")
cat(" done\n")

# Fit 3: Full decreasing (tau = 0.5)
cat("  - Full decreasing (tau = 0.5)...")
res_decroissant <- SplineConstQuantRegBs3(xtab, ytab, knots, tau = 0.5,
                                          monot = -1, convcons = 0,
                                          solver = "OSQP")
cat(" done\n")

# Fit 4: Convexity constraint (tau = 0.5)
cat("  - Convexity constraint (tau = 0.5)...")
res_convexe <- SplineConstQuantRegBs3(xtab, ytab, knots, tau = 0.5,
                                      monot = 0, convcons = 1,
                                      solver = "OSQP")
cat(" done\n\n")

# Evaluate all fits on a fine grid
x_eval <- seq(0, 1, length.out = 200)
y_uncon1 <- spline_eval(res_uncon1, x_eval)
y_uncon2 <- spline_eval(res_uncon2, x_eval)
y_uncon3 <- spline_eval(res_uncon3, x_eval)
y_croiss1 <- spline_eval(res_croissant1, x_eval)
y_croiss2 <- spline_eval(res_croissant2, x_eval)
y_decroiss <- spline_eval(res_decroissant, x_eval)
y_convexe <- spline_eval(res_convexe, x_eval)

# Compute derivatives for verification
compute_derivative <- function(y, x) {
  diff(y) / diff(x)
}

deriv_croiss <- compute_derivative(y_croiss, x_eval)
deriv_decroiss <- compute_derivative(y_decroiss, x_eval)

cat("Verification of constraints:\n")
cat("----------------------------\n")
cat(sprintf("  Partial increasing: min derivative = %.4f (should be >= 0 on constrained region)\n",
            min(deriv_croiss[1:140])))  # first 140 points ~ first 7 intervals
cat(sprintf("  Full decreasing:    max derivative = %.4f (should be <= 0)\n",
            max(deriv_decroiss)))
cat("\n")

# Create plots
par(mfrow = c(2, 2), mar = c(4, 4, 4, 2))

# Plot 1: Partial increasing constraint
plot(xtab, ytab, pch = 16, cex = 0.5, col = "black",
     xlab = "x", ylab = "y",
     main = "Partial Increasing Constraint\n(first 7 intervals only)\n tau=0.5, 0.1")
lines(x_eval, y_croiss1, col = "blue", lwd = 2)
lines(x_eval, y_croiss2, col = "blue", lwd = 2)
abline(v = knots, col = "blue", lty = 2, lwd = 0.5)
# Highlight the constrained region
abline(v = knots[n_constrained + 1], col = "red", lty = 2, lwd = 2)
text(knots[n_constrained + 1] + 0.02, max(ytab) - 0.2,
     "Constrained region", col = "red", cex = 0.7, srt = 90)
grid()

# Plot 2: Full decreasing constraint
plot(xtab, ytab, pch = 16, cex = 0.5, col = "black",
     xlab = "x", ylab = "y",
     main = "Full Decreasing Constraint\n(everywhere)")
lines(x_eval, y_decroiss, col = "darkgreen", lwd = 2)
abline(v = knots, col = "blue", lty = 2, lwd = 0.5)
grid()

# Plot 3: Convexity constraint
plot(xtab, ytab, pch = 16, cex = 0.5, col = "black",
     xlab = "x", ylab = "y",
     main = "Convexity Constraint\n(second derivative >= 0)")
lines(x_eval, y_convexe, col = "purple", lwd = 2)
abline(v = knots, col = "blue", lty = 2, lwd = 0.5)
grid()

# Plot 4: Unconstrained (tau = 0.9)
plot(xtab, ytab, pch = 16, cex = 0.5, col = "black",
     xlab = "x", ylab = "y",
     main = "Unconstrained\n(tau = 0.9, , 0.5, 0.1)")
lines(x_eval, y_uncon1, col = "red", lwd = 2)
lines(x_eval, y_uncon2, col = "red", lwd = 2)
lines(x_eval, y_uncon3, col = "red", lwd = 2)
abline(v = knots, col = "blue", lty = 2, lwd = 0.5)
grid()

cat("==============================================\n")
cat("Demo completed. Check the plots for comparison.\n")
cat("==============================================\n")
