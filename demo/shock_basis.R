# inst/demo/singularity.R
# ============================================================================
# Demo: Discontinuity and Derivability Singularity
# ============================================================================
# This demo shows how B-spline knot multiplicity controls regularity:
#   - Multiplicity d+1 -> C^-1 : discontinuity (jump)
#   - Multiplicity d   -> C^0  : continuous but non-differentiable (kink, like |x|)
#   - Multiplicity 1   -> C^(d-1) : fully smooth
#
# The knot vectors are built dynamically to work with any degree.
# The corresponding B-spline basis is displayed to show how the
# singularities propagate from the basis to the fit.
# ============================================================================

oldpar <- par(mfrow = c(2,2))


# ---------------------------------------------------------------------------
# Data with both a jump and a kink
# ---------------------------------------------------------------------------
n <- 300
x <- seq(-1, 2, length.out = n)

y <- numeric(n)
for (i in seq_along(x)) {
  if (x[i] < 0)      y[i] <- 0.5 + abs(x[i]) + 0.1 * rnorm(1)
  else if (x[i] < 1) y[i] <- 0.5 + x[i] + 0.1 * rnorm(1)
  else               y[i] <- 2.0 + 0.3 * (x[i] - 1) + 0.1 * rnorm(1)
}


if (!exists("degree")) {
  degree <- 3
}else {degree<-max(degree,1)}


if (degree < 1) stop("Degree must be >= 1")


m_kink <- degree       # Multiplicity for C^0 (kink)
m_jump <- degree + 1   # Multiplicity for C^-1 (jump)

cat("=== Demo with degree =", degree, "===\n")
cat("Kink (x = 0) multiplicity =", m_kink, "(C^0)\n")
cat("Jump (x = 1) multiplicity =", m_jump, "(C^-1)\n\n")

# ---------------------------------------------------------------------------
# Knot vectors
# ---------------------------------------------------------------------------

# (a) No multiplicity: smooth everywhere
knot_a <- c(-1, 0, 0.5, 1, 2)

# (b) Kink at 0, jump at 1
knot_b <- c(-1, rep(0, m_kink), 0.5, rep(1, m_jump), 2)

# (c) Only jump at 1, smooth kink at 0
knot_c <- c(-1, 0, 0.5, rep(1, m_jump), 2)

# Helper to build extended knot sequence (repeated ends)
build_extended_knot <- function(knot, degree) {
  t0 <- knot[1]
  tn <- knot[length(knot)]
  c(rep(t0, degree), knot, rep(tn, degree))
}

# ---------------------------------------------------------------------------
# Fits
# ---------------------------------------------------------------------------

fit_a <- tryCatch(
  quantile_spline(x, y, knot = knot_a, tau = 0.5, degree = degree,
                  monot = 0, convcons = 0, der3cons = 0,
                  solver = "ECOS", callable = TRUE),
  error = function(e) NULL)

fit_b <- tryCatch(
  quantile_spline(x, y, knot = knot_b, tau = 0.5, degree = degree,
                  monot = 0, convcons = 0, der3cons = 0,
                  solver = "ECOS", callable = TRUE),
  error = function(e) NULL)

fit_c <- tryCatch(
  quantile_spline(x, y, knot = knot_c, tau = 0.5, degree = degree,
                  monot = 0, convcons = 0, der3cons = 0,
                  solver = "ECOS", callable = TRUE),
  error = function(e) NULL)

safe_eval <- function(fit, xv) {
  if (is.null(fit)) return(rep(NA, length(xv)))
  tryCatch(fit(xv), error = function(e) rep(NA, length(xv)))
}

# ---------------------------------------------------------------------------
# Build B-spline bases
# ---------------------------------------------------------------------------

x_eval <- seq(-1, 2, length.out = 500)

build_basis <- function(knot, degree) {
  sn <- build_extended_knot(knot, degree)
  BB <- tryCatch(Bspline_base(sn, degree = degree), error = function(e) NULL)
  if (is.null(BB)) return(NULL)
  Bv <- tryCatch(bs_direct(BB, x_eval), error = function(e) NULL)
  list(BB = BB, Bv = Bv, sn = sn)
}

basis_a <- build_basis(knot_a, degree)
basis_b <- build_basis(knot_b, degree)
basis_c <- build_basis(knot_c, degree)

# ---------------------------------------------------------------------------
# Helper: plot basis with singularity markers
# ---------------------------------------------------------------------------
plot_basis <- function(basis, title, kink_x = 0, jump_x = 1, colors = NULL) {
  if (is.null(basis) || is.null(basis$Bv)) {
    plot.new(); title(main = paste(title, "(basis unavailable)"))
    return(invisible())
  }

  Bv <- basis$Bv
  nsplines <- nrow(Bv)

  if (is.null(colors)) {
    colors <- rainbow(nsplines, s = 0.7, v = 0.8)
  }

  ylim <- range(Bv, na.rm = TRUE)
  ylim <- ylim + c(-0.05, 0.15) * diff(ylim)

  matplot(x_eval, t(Bv),
          type = "l", lty = 1, lwd = 1.5, col = colors,
          xlab = "x", ylab = "Basis value",
          main = title, ylim = ylim)

  # Mark singularities
  abline(v = kink_x, col = "blue", lty = 2, lwd = 1.5)
  abline(v = jump_x, col = "red", lty = 2, lwd = 1.5)

  # Count knots and multiplicities
  info <- paste0("degree=", degree,
                 " | n_splines=", nsplines)
  mtext(info, side = 3, line = 0.2, cex = 0.7, col = "gray40")

  legend("topright",
         legend = c("kink (C^0)", "jump (C^-1)"),
         col = c("blue", "red"), lty = 2, lwd = 1.5,
         cex = 0.7, bg = "white")
}

# ---------------------------------------------------------------------------
# Visualization: 3 rows (fit) x 2 columns (fit + basis)
# ---------------------------------------------------------------------------
par(mfrow = c(3, 2), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))

# --- Row 1 : smooth fit + smooth basis ---
plot(x, y, pch = 16, cex = 0.25, col = "gray60",
     main = paste0("(a) Smooth fit (d=", degree, ")"),
     xlab = "x", ylab = "y")
lines(x_eval, safe_eval(fit_a, x_eval), col = "darkgreen", lwd = 2)
abline(v = c(0, 1), col = c("blue", "red"), lty = 2)
legend("topleft", paste0("C^", degree - 1, " everywhere (smooths features)"),
       col = "darkgreen", lty = 1, lwd = 2, cex = 0.7)

plot_basis(basis_a, paste0("(a') Basis: no multiplicity"))

# --- Row 2 : kink + jump fit + singular basis ---
plot(x, y, pch = 16, cex = 0.25, col = "gray60",
     main = paste0("(b) Kink (m=", m_kink, ") + Jump (m=", m_jump, ")"),
     xlab = "x", ylab = "y")
lines(x_eval, safe_eval(fit_b, x_eval), col = "purple", lwd = 2)
abline(v = c(0, 1), col = c("blue", "red"), lty = 2)
legend("topleft", "C^0 at 0, C^-1 at 1 (correct)",
       col = "purple", lty = 1, lwd = 2, cex = 0.7)

plot_basis(basis_b, paste0("(b') Basis: kink (m=", m_kink, ") + jump (m=", m_jump, ")"))

# --- Row 3 : jump only fit + basis with jump ---
plot(x, y, pch = 16, cex = 0.25, col = "gray60",
     main = paste0("(c) Smooth kink + Jump (m=", m_jump, ")"),
     xlab = "x", ylab = "y")
lines(x_eval, safe_eval(fit_c, x_eval), col = "darkorange", lwd = 2)
abline(v = c(0, 1), col = c("blue", "red"), lty = 2)
legend("topleft", paste0("C^", degree - 1, " at 0, C^-1 at 1"),
       col = "darkorange", lty = 1, lwd = 2, cex = 0.7)

plot_basis(basis_c, paste0("(c') Basis: smooth kink + jump (m=", m_jump, ")"))

mtext(paste0("Regularity control with knot multiplicity (degree = ", degree, ")"),
      outer = TRUE, cex = 0.9, font = 2)

# ---------------------------------------------------------------------------
# Console summary
# ---------------------------------------------------------------------------
cat("=== Regularity control with knot multiplicity (degree =", degree, ") ===\n")
cat("Multiplicity m -> regularity C^(d-m)\n")
cat("  m = 1        : C^", degree - 1, " (smooth)\n", sep = "")
if (degree >= 2)
  cat("  m = 2        : C^", degree - 2, " (curvature jumps)\n", sep = "")
if (degree >= 3)
  cat("  m = 3        : C^", degree - 3, " (third derivative jumps)\n", sep = "")
if (degree >= 1)
  cat("  m = ", degree, "       : C^0  (kink, like |x|)\n", sep = "")
cat("  m = ", degree + 1, "       : C^-1 (discontinuity/jump)\n\n", sep = "")

rss <- function(fit) {
  if (is.null(fit)) return(NA)
  sum((y - fit(x))^2, na.rm = TRUE)
}
cat("Residual sum of squares:\n")
cat("  (a) Smooth everywhere :", round(rss(fit_a), 3), "\n")
cat("  (b) Kink + jump       :", round(rss(fit_b), 3), "\n")
cat("  (c) Jump only         :", round(rss(fit_c), 3), "\n\n")

# Summary of basis dimensions
cat("B-spline basis dimensions (n_splines):\n")
if (!is.null(basis_a$Bv)) cat("  (a) Smooth            :", nrow(basis_a$Bv), "\n")
if (!is.null(basis_b$Bv)) cat("  (b) Kink + jump       :", nrow(basis_b$Bv), "\n")
if (!is.null(basis_c$Bv)) cat("  (c) Jump only         :", nrow(basis_c$Bv), "\n")


par(oldpar)
cat("\nDemo completed.\n")

