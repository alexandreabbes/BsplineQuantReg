# Demo: Monotonicity constraints with Karlin-Studden
# Author: Alexandre Abbes
# This demo shows quantile regression with increasing/decreasing constraints

# Load the package
library(BsplineQuantReg)

oldpar <- par(mfrow = c(2,2))

# Generate data
#set.seed(42)
n_points <- 100
xtab <- seq(0, 1, length.out = n_points)

# Increasing function with noise
ytab <- 2 * xtab + 0.5 * sin(4 * pi * xtab) + 0.05 * rnorm(n_points)

# Knots
kn <- 10
knots <- quantile(xtab, probs = seq(0, 1, length.out = kn + 1))

# Unconstrained fit
fit_uncon <- SplineConstQuantRegBs3(xtab, ytab, knots, tau = 0.5,
                                    monot = 0, convcons = 0)

# Increasing constraint
fit_inc <- SplineConstQuantRegBs3(xtab, ytab, knots, tau = 0.5,
                                  monot = 1, convcons = 0)

# Decreasing constraint
fit_dec <- SplineConstQuantRegBs3(xtab, ytab, knots, tau = 0.5,
                                  monot = -1, convcons = 0)

# Plot results
x_eval <- seq(0, 1, length.out = 200)
y_uncon <- spline_eval(fit_uncon, x_eval)
y_inc <- spline_eval(fit_inc, x_eval)
y_dec <- spline_eval(fit_dec, x_eval)

par(mfrow = c(2, 2))
plot(xtab, ytab, pch = 16, cex = 0.5, col = "gray", main = "Data")
plot(xtab, ytab, pch = 16, cex = 0.5, col = "gray", main = "Unconstrained")
lines(x_eval, y_uncon, col = "red", lwd = 2)
plot(xtab, ytab, pch = 16, cex = 0.5, col = "gray", main = "Increasing")
lines(x_eval, y_inc, col = "blue", lwd = 2)
plot(xtab, ytab, pch = 16, cex = 0.5, col = "gray", main = "Decreasing")
lines(x_eval, y_dec, col = "green", lwd = 2)


par(oldpar)
cat("\nDemo completed. The increasing fit should have non-negative slope.\n")
