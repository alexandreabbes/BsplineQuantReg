# Demo: Logistic Curve Regression (He & Ng, 1999)
# This demo shows quantile regression on a logistic/sigmoid function
# Author: Alexandre Abbes

library(BsplineQuantReg)

cat("========================================\n")
cat("Demo: Logistic Curve Quantile Regression\n")
cat("========================================\n\n")

oldpar <- par(mfrow = c(2,2))

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

cat(sprintf("Number of points: %d\n", n_points))
cat(sprintf("Number of intervals: %d\n", kn))
cat(sprintf("Number of basis functions: %d\n\n", kn + 3))

# Fit multiple quantiles
tau_values <- c(0.1, 0.25, 0.5, 0.75, 0.9)
colors <- c("red", "orange", "blue", "green", "purple")
fits <- list()


convex=rep(0,  (kn+1))
for (i in  1:(kn+1) )
  {
  if (i<6) {convex[i]<-1}
  if  (i>6) {convex[i<--1]}
}

fit_monot_conv <- SplineConstQuantRegBs3(x, y, knots, tau = 0.5,
                                    monot = 1, convcons = convex)
y_monot_conv <- spline_eval(fit_monot_conv, x_eval)

plot(x, y, pch = 16, cex = 0.4, col = "black",
     xlab = "x", ylab = "y",
     main = "Median Regression with Monotonicity \n and convexity")
lines(x_eval, y_true, col = "black", lwd = 2, lty = 2)
lines(x_eval, y_monot_conv, col = "blue", lwd = 2)




par(oldpar)


cat("\nDemo completed.\n")
