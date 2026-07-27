# Demo: B-spline Derivatives
# Author: Alexandre Abbes
#
# This demo shows how to compute and visualize derivatives of B-splines
# using the Bspline_deriv and Bspline_base_deriv functions.

library(BsplineQuantReg)

oldpar <- par(mfrow = c(2,2))
cat("========================================\n")
cat("Demo: B-spline Derivatives\n")
cat("========================================\n\n")

# ============================================================
# 1. Create a B-spline basis
# ============================================================

cat("1. Creating B-spline basis...\n")

# Knots
knots <- c(0, 1, 2, 3, 4, 5)
degree <- 3
sn <- c(rep(knots[1], degree), knots, rep(rev(knots)[1], degree))

# Build basis
basis <- Bspline_base(sn, degree = degree)
n_splines <- basis$n_splines
cat(sprintf("  Degree: %d\n", degree))
cat(sprintf("  Number of basis functions: %d\n", n_splines))
cat(sprintf("  Knots: %s\n", paste(knots, collapse = ", ")))

# Create a random spline
set.seed(42)
coeff <- rnorm(n_splines, mean = 0, sd = 2)
basis$coeff <- coeff

cat("\n")

# ============================================================
# 2. Visualize the original basis
# ============================================================

cat("2. Computing the original B-spline basis...\n")

# Ouvrir une nouvelle fenetre graphique pour la base
#dev.new()
#view_basis(basis)
#title("Original B-spline Basis (degree 3)")

#cat("  Original basis plotted.\n\n")

# ============================================================
# 3. Compute derivatives using Bspline_deriv
# ============================================================

cat("3. Computing derivatives with Bspline_deriv...\n")

# Original spline
spline_orig <- basis
class(spline_orig) <- "non_callable_spline"

# Derivatives
der1 <- Bspline_deriv(spline_orig, der = 1)
der2 <- Bspline_deriv(spline_orig, der = 2)
der3 <- Bspline_deriv(spline_orig, der = 3)

cat(sprintf("  Original: degree %d, coefficients %d\n",
            spline_orig$degree, length(spline_orig$coeff)))
cat(sprintf("  First derivative: degree %d, coefficients %d\n",
            der1$degree, length(der1$coeff)))
cat(sprintf("  Second derivative: degree %d, coefficients %d\n",
            der2$degree, length(der2$coeff)))
cat(sprintf("  Third derivative: degree %d, coefficients %d\n",
            der3$degree, length(der3$coeff)))

cat("\n")

# ============================================================
# 4. Compute derivatives using Bspline_base_deriv
# ============================================================

cat("4. Computing derivatives with Bspline_base_deriv...\n")

# Derivative bases
basis_der1 <- Bspline_base_deriv(basis, der = 1)
basis_der2 <- Bspline_base_deriv(basis, der = 2)
basis_der3 <- Bspline_base_deriv(basis, der = 3)

cat(sprintf("  Derivative basis 1: degree %d\n", basis_der1$degree))
cat(sprintf("  Derivative basis 2: degree %d\n", basis_der2$degree))
cat(sprintf("  Derivative basis 3: degree %d\n", basis_der3$degree))

cat("\n")

# ============================================================
# 5. Visualize derivative bases
# ============================================================

cat("5. Visualizing derivative bases...\n")

# Ouvrir une nouvelle fenetre pour les bases derivees
dev.new()
par(mfrow = c(2, 2), mar = c(4, 4, 4, 2))

# Base originale
view_basis(basis)
title("Original Basis (deg 3)")

# Base derivee 1
view_basis(basis_der1)
title("First Derivative Basis (deg 2)")

# Base derivee 2
view_basis(basis_der2)
title("Second Derivative Basis (deg 1)")

# Base derivee 3
view_basis(basis_der3)
title("Third Derivative Basis (deg 0)")

cat("  Derivative bases plotted.\n\n")

# Restaurer les parametres graphiques
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2))

# ============================================================
# 6. Evaluate and plot the spline and its derivatives
# ============================================================

cat("6. Evaluating and plotting the spline and its derivatives...\n")

x_eval <- seq(0, 5, length.out = 300)

# Evaluer la spline et ses derivees
y_orig <- spline_eval(spline_orig, x_eval)
y_der1 <- spline_eval(spline_orig, x_eval, der = 1)
y_der2 <- spline_eval(spline_orig, x_eval, der = 2)
y_der3 <- spline_eval(spline_orig, x_eval, der = 3)

# Evaluer egalement avec les coefficients derives
y_der1_coeff <- spline_eval(der1, x_eval)
y_der2_coeff <- spline_eval(der2, x_eval)
y_der3_coeff <- spline_eval(der3, x_eval)

# Verifier la coherence des deux methodes
max_diff1 <- max(abs(y_der1 - y_der1_coeff), na.rm = TRUE)
max_diff2 <- max(abs(y_der2 - y_der2_coeff), na.rm = TRUE)
max_diff3 <- max(abs(y_der3 - y_der3_coeff), na.rm = TRUE)

cat(sprintf("  Max difference between methods (1st derivative): %.2e\n", max_diff1))
cat(sprintf("  Max difference between methods (2nd derivative): %.2e\n", max_diff2))
cat(sprintf("  Max difference between methods (3rd derivative): %.2e\n", max_diff3))
cat("max diff are :", paste(c(max_diff1,max_diff2,max_diff3),collapse="  "),"\n")
cat("  Both methods are consistent up to", max(max_diff1,max_diff2,max_diff3),"\n")
cat("\n")

# Ouvrir une nouvelle fenetre pour les graphiques
dev.new()
oldpar <- par(mfrow = c(2, 3), mar = c(4, 4, 4, 2))

# Plot 1: Original spline
plot(x_eval, y_orig, type = "l", col = "blue", lwd = 2,
     xlab = "x", ylab = "y",
     main = "Original Spline (deg 3)")
abline(v = knots, col = "red", lty = 2, lwd = 0.5)
grid()

# Plot 2: First derivative
plot(x_eval, y_der1, type = "l", col = "darkgreen", lwd = 2,
     xlab = "x", ylab = "y'",
     main = "First Derivative (deg 2)")
abline(v = knots, col = "red", lty = 2, lwd = 0.5)
abline(h = 0, col = "gray", lty = 3)
grid()

# Plot 3: Second derivative
plot(x_eval, y_der2, type = "l", col = "purple", lwd = 2,
     xlab = "x", ylab = "y''",
     main = "Second Derivative (deg 1)")
abline(v = knots, col = "red", lty = 2, lwd = 0.5)
abline(h = 0, col = "gray", lty = 3)
grid()

# Plot 4: Third derivative
plot(x_eval, y_der3, type = "l", col = "orange", lwd = 2,
     xlab = "x", ylab = "y'''",
     main = "Third Derivative (deg 0)")
abline(v = knots, col = "red", lty = 2, lwd = 0.5)
abline(h = 0, col = "gray", lty = 3)
grid()

# Plot 5: Comparison of derivatives (1st derivative)
plot(x_eval, y_der1, type = "l", col = "darkgreen", lwd = 2,
     xlab = "x", ylab = "y'",
     main = "1st Derivative: Both Methods")
lines(x_eval, y_der1_coeff, col = "red", lwd = 2, lty = 2)
legend("topright", legend = c("spline_eval(der=1)", "Bspline_deriv"),
       col = c("darkgreen", "red"), lty = c(1, 2), lwd = 2, cex = 0.7)
grid()

# Plot 6: Comparison of derivatives (2nd derivative)
plot(x_eval, y_der2, type = "l", col = "purple", lwd = 2,
     xlab = "x", ylab = "y''",
     main = "2nd Derivative: Both Methods")
lines(x_eval, y_der2_coeff, col = "red", lwd = 2, lty = 2)
legend("topright", legend = c("spline_eval(der=2)", "Bspline_deriv"),
       col = c("purple", "red"), lty = c(1, 2), lwd = 2, cex = 0.7)
grid()

# Restaurer les parametres graphiques
par(oldpar)

cat("  Plots created.\n\n")

# ============================================================
# 7. Callable derivatives
# ============================================================

cat("7. Creating callable derivatives...\n")

# Make the original spline callable
spline_call <- make_spline(spline_orig, callable = TRUE)

# Compute derivatives (they inherit callable status)
der1_call <- Bspline_deriv(spline_call, der = 1, callable = TRUE)
der2_call <- Bspline_deriv(spline_call, der = 2, callable = TRUE)

# Evaluate directly
y_der1_call <- der1_call(x_eval)
y_der2_call <- der2_call(x_eval)

cat("  Callable derivatives created.\n")
cat("  Usage: der1_call(x) or der2_call(x)\n\n")

# ============================================================
# 8. Summary
# ============================================================

cat("========================================\n")
cat("Summary\n")
cat("========================================\n")
cat("B-spline derivatives can be computed in two conceptually (but equivalent) different ways:\n")
cat("  1. Bspline_deriv(): computes coefficients of derivative B-spline
    on a natural b-spline basis. It can handle both class : 'callable_spline' and 'non_callable_spline'\n")
cat("  2. Bspline_base_deriv(): differentiate the basis functions and keep
    the same coefficient on this derivative basis\n")
cat("The functions are all compatible. the funcitons 'spline_eval()', and 'Bspline_basis()' can also receive a derivative parameter.\n")
cat("Basis visualization: 'view_basis()' displays the B-spline basis functions\n")

cat("\nAll Methods are consistent and produce the same results.\n")
cat("Demo completed.\n")

par(oldpar)
