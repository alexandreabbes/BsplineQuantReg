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
# ============================================================================

oldpar <- par(mfrow = c(2,2))




# ---------------------------------------------------------------------------
# Data with both a jump and a kink
# ---------------------------------------------------------------------------
n <- 300
x <- seq(-1, 2, length.out = n)

# Kink at x = 0 (like |x|), jump at x = 1
y <- numeric(n)
for (i in seq_along(x)) {
  if (x[i] < 0)      y[i] <- 0.5 + abs(x[i]) + 0.1 * rnorm(1)
  else if (x[i] < 1) y[i] <- 0.5 + x[i] + 0.1 * rnorm(1)
  else               y[i] <- 2.0 + 0.3 * (x[i] - 1) + 0.1 * rnorm(1)
}

# Degree can be set to any value >= 1
# For the demo to work, the number of intervals must be sufficient
# (at least degree + 1 intervals)

if (!exists("degree")) {
  degree <- 3
}else {degree<-max(degree,1)}



# Helper to build a knot vector with given multiplicities at specified points
# intervals: matrix with columns (left, right) for each interval
# special: named list of multiplicities at specific points
build_demo_knot <- function(intervals, special = list()) {
  # Build regular knots from intervals
  knots <- sort(unique(c(intervals)))
  # Apply multiplicities
  result <- c()
  for (k in knots) {
    if (as.character(k) %in% names(special)) {
      result <- c(result, rep(k, special[[as.character(k)]]))
    } else {
      result <- c(result, k)
    }
  }
  return(result)
}

# ---------------------------------------------------------------------------
# Check degree feasibility
# ---------------------------------------------------------------------------
if (degree < 1) stop("Degree must be >= 1")

# Multiplicity needed for a kink (C^0) = degree
# Multiplicity needed for a jump (C^-1) = degree + 1
m_kink <- degree
m_jump <- degree + 1

cat("=== Demo with degree =", degree, "===\n")
cat("Kink (x = 0) multiplicity =", m_kink, "(C^0)\n")
cat("Jump (x = 1) multiplicity =", m_jump, "(C^-1)\n\n")

# ---------------------------------------------------------------------------
# Three fits with increasing regularity control
# ---------------------------------------------------------------------------

# (a) No multiplicity: smooth everywhere (poor fit)
knot_a <- c(-1, 0, 0.5, 1, 2)
fit_a <- tryCatch(
  quantile_spline(x, y, knot = knot_a, tau = 0.5, degree = degree,
                  monot = 0, convcons = 0, der3cons = 0,
                  solver = "ECOS", callable = TRUE),
  error = function(e) NULL
)

# (b) Kink at 0 (mult = degree), jump at 1 (mult = degree + 1)
knot_b <- c(-1,
            rep(0, m_kink),
            0.5,
            rep(1, m_jump),
            2)
fit_b <- tryCatch(
  quantile_spline(x, y, knot = knot_b, tau = 0.5, degree = degree,
                  monot = 0, convcons = 0, der3cons = 0,
                  solver = "ECOS", callable = TRUE),
  error = function(e) NULL
)

# (c) Only jump at 1, smooth kink at 0
knot_c <- c(-1,
            0,
            0.5,
            rep(1, m_jump),
            2)
fit_c <- tryCatch(
  quantile_spline(x, y, knot = knot_c, tau = 0.5, degree = degree,
                  monot = 0, convcons = 0, der3cons = 0,
                  solver = "ECOS", callable = TRUE),
  error = function(e) NULL
)

# Fallback: if any fit failed, replace with NA function
safe_eval <- function(fit, xv) {
  if (is.null(fit)) return(rep(NA, length(xv)))
  tryCatch(fit(xv), error = function(e) rep(NA, length(xv)))
}

# ---------------------------------------------------------------------------
# Visualization (single page, 2x2)
# ---------------------------------------------------------------------------
x_eval <- seq(-1, 2, length.out = 500)

par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

# Data
plot(x, y, pch = 16, cex = 0.3, col = "gray60",
     main = "Data: kink (x=0) + jump (x=1)",
     xlab = "x", ylab = "y")
abline(v = c(0, 1), col = c("blue", "red"), lty = 2, lwd = 2)
legend("topleft", c("kink", "jump"),
       col = c("blue", "red"), lty = 2, lwd = 2, cex = 0.8)

# (a) Smooth fit
plot(x, y, pch = 16, cex = 0.3, col = "gray60",
     main = paste0("Smooth (no multiplicity, d=", degree, ")"),
     xlab = "x", ylab = "y")
lines(x_eval, safe_eval(fit_a, x_eval), col = "darkgreen", lwd = 2)
abline(v = c(0, 1), col = c("blue", "red"), lty = 2)
legend("topleft", paste0("C^", degree - 1, " everywhere (smooths features)"),
       col = "darkgreen", lty = 1, lwd = 2, cex = 0.75)

# (b) Kink + jump
plot(x, y, pch = 16, cex = 0.3, col = "gray60",
     main = paste0("Kink (m=", m_kink, ") + Jump (m=", m_jump, ")"),
     xlab = "x", ylab = "y")
lines(x_eval, safe_eval(fit_b, x_eval), col = "purple", lwd = 2)
abline(v = c(0, 1), col = c("blue", "red"), lty = 2)
legend("topleft", "C^0 at 0, C^-1 at 1 (correct)",
       col = "purple", lty = 1, lwd = 2, cex = 0.75)

# (c) Only jump
plot(x, y, pch = 16, cex = 0.3, col = "gray60",
     main = paste0("Smooth kink + Jump (m=", m_jump, ")"),
     xlab = "x", ylab = "y")
lines(x_eval, safe_eval(fit_c, x_eval), col = "darkorange", lwd = 2)
abline(v = c(0, 1), col = c("blue", "red"), lty = 2)
legend("topleft", paste0("C^", degree - 1, " at 0, C^-1 at 1"),
       col = "darkorange", lty = 1, lwd = 2, cex = 0.75)

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
cat("  (c) Jump only         :", round(rss(fit_c), 3), "\n")


par(oldpar)

cat("\nDemo completed.\n")
