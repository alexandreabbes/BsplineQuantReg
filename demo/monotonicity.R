# Demo: Monotonicity constraints with Karlin-Studden
# Author: Alexandre Abbes
# This demo shows quantile regression with increasing/decreasing constraints

# Load the package
library(BsplineQuantReg)

oldpar <- par(mfrow = c(2,2))

# Generate data
#set.seed(42)
n_points <- 50
xtab <- seq(0, 1, length.out = n_points)

# Increasing function with noise
ytab <- 2 * xtab + 0.5 * sin(4 * pi * xtab) + 0.2 * rnorm(n_points)

# Knots
kn <- 10
knots <- quantile(xtab, probs = seq(0, 1, length.out = kn + 1))

# Unconstrained fit
fit_uncon <- SplineConstQuantRegBs3(xtab, ytab, knots, tau = 0.5,
                                    monot = 0, convcons = 0)

# Increasing constraint
fit_inc2 <- SplineConstQuantRegBs3(xtab, ytab, knots, tau = 0.5,
                                  monot = 1, convcons = 0)
fit_inc1 <- SplineConstQuantRegBs3(xtab, ytab, knots, tau = 0.1,
                                   monot = 1, convcons = 0)
fit_inc3 <- SplineConstQuantRegBs3(xtab, ytab, knots, tau = 0.9,
                                   monot = 1, convcons = 0)
# Decreasing constraint
fit_dec <- SplineConstQuantRegBs3(xtab, ytab, knots, tau = 0.5,
                                  monot = -1, convcons = 0)

# Plot results
x_eval <- seq(0, 1, length.out = 200)
y_uncon <- spline_eval(fit_uncon, x_eval)
y_inc2 <- spline_eval(fit_inc2, x_eval)
y_inc3 <- spline_eval(fit_inc3, x_eval)
y_inc1 <- spline_eval(fit_inc1, x_eval)
y_dec <- spline_eval(fit_dec, x_eval)

par(mfrow = c(2, 2))
plot(xtab, ytab, pch = 16,  col = "black", main = "Data")
plot(xtab, ytab, pch = 16,  col = "black", main = "Unconstrained")
lines(x_eval, y_uncon, col = "red", lwd = 2)
plot(xtab, ytab, pch = 16,  col = "black", main = "Increasing")
lines(x_eval, y_inc1,col='yellow', lwd = 2)
lines(x_eval, y_inc2,col='red', lwd = 2)
lines(x_eval, y_inc3,col='brown', lwd = 2)
legend("topleft", legend = c(paste("tau =", c(0.1,0.5,0.9))),
       col = colors, lty = c(rep(1, 3)),
       lwd = 1.5, cex = 0.6)

plot(xtab, ytab, pch = 16,  col = "black", main = "Decreasing")
lines(x_eval, y_dec, col = "green", lwd = 2)


par(oldpar)
cat("\nDemo completed. The increasing fit should have non-negative slope.\n")
