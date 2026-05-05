# Demo: Temperature Anomaly Trend Analysis
# This demo shows temperature data with various monotonicity constraints
# Data source: Global temperature anomalies (1880-1992)
#
# Constraint scenarios:
# 1. Unconstrained fit
# 2. Full increasing monotonicity (too restrictive)
# 3. Partial increasing only after 1970
# 4. Mixed constraints: decreasing (1945-1970) + increasing elsewhere

library(ConstrainedQuantileSplines)

cat("========================================\n")
cat("Demo: Temperature Anomaly Trend Analysis\n")
cat("========================================\n\n")

# Create temperature data (1880-1992)
years <- 1880:1992
temperature <- c(
  -0.32, -0.32, -0.40, -0.39, -0.65, -0.43, -0.40, -0.52, -0.30, -0.12,
  -0.40, -0.42, -0.39, -0.45, -0.35, -0.36, -0.19, -0.14, -0.37, -0.22,
  0.00, -0.08, -0.24, -0.36, -0.49, -0.27, -0.19, -0.43, -0.29, -0.30,
  -0.29, -0.29, -0.28, -0.23, -0.04, -0.02, -0.24, -0.42, -0.35, -0.16,
  -0.17, -0.09, -0.13, -0.16, -0.14, -0.14,  0.10, -0.03,  0.03, -0.18,
  -0.06,  0.04,  0.02, -0.13,  0.03, -0.06,  0.02,  0.13,  0.13, -0.03,
  0.15,  0.12,  0.10,  0.04,  0.11, -0.04,  0.01,  0.13, -0.01, -0.06,
  -0.14, -0.02,  0.04,  0.14, -0.07, -0.06, -0.17,  0.10,  0.10,  0.05,
  -0.01,  0.08,  0.02,  0.02, -0.26, -0.16, -0.09, -0.02, -0.12,  0.03,
  0.04, -0.11, -0.07,  0.19, -0.07, -0.05, -0.22,  0.16,  0.09,  0.14,
  0.28,  0.39,  0.07,  0.29,  0.11,  0.11,  0.16,  0.32,  0.35,  0.25,
  0.47,  0.41,  0.13
)

cat(sprintf("Years: %d to %d (%d observations)\n",
            min(years), max(years), length(years)))
cat(sprintf("Temperature range: [%.2f, %.2f] °C\n\n",
            min(temperature), max(temperature)))

# Normalize years to [0,1] for numerical stability
xtab <- (years - min(years)) / (max(years) - min(years))
ytab <- temperature

# Define knots at specific years
target_knots <- c(1917, 1936, 1945, 1970)
knots_years <- sort(c(min(years), target_knots, max(years)))
knots <- (knots_years - min(years)) / (max(years) - min(years))

cat("Knots at years:", knots_years, "\n\n")

# Create evaluation grid
x_eval <- seq(0, 1, length.out = 300)
years_eval <- min(years) + x_eval * (max(years) - min(years))

# Helper function to get interval index for a given year
get_interval_idx <- function(year, knots_years) {
  which(knots_years <= year)[length(which(knots_years <= year))]
}

# ============================================================
# Model 1: Unconstrained median regression
# ============================================================
cat("=== 1. Unconstrained median regression ===\n")
fit_uncon <- SplineConstQuantRegBs3(xtab, ytab, knots, tau = 0.5,
                                    monot = 0, convcons = 0)
y_uncon <- spline_eval(fit_uncon, x_eval)

# ============================================================
# Model 2: Full monotonicity (increasing everywhere)
# ============================================================
cat("\n=== 2. Full increasing monotonicity (everywhere) ===\n")
fit_full_inc <- SplineConstQuantRegBs3(xtab, ytab, knots, tau = 0.5,
                                       monot = 1, convcons = 0)
y_full_inc <- spline_eval(fit_full_inc, x_eval)

# ============================================================
# Model 3: Partial monotonicity (only after 1970)
# ============================================================
break_year_1970 <- 1970
monot_partial <- rep(0, length(knots) - 1)
for (i in 1:length(monot_partial)) {
  interval_end <- knots_years[i + 1]
  if (interval_end > break_year_1970) {
    monot_partial[i] <- 1
  }
}
cat("\n=== 3. Partial monotonicity (increasing only after 1970) ===\n")
cat("   Constrained intervals:", which(monot_partial == 1), "\n")
fit_partial <- SplineConstQuantRegBs3(xtab, ytab, knots, tau = 0.5,
                                      monot = monot_partial, convcons = 0)
y_partial <- spline_eval(fit_partial, x_eval)

# ============================================================
# Model 4: Mixed constraints
# Decreasing between 1945 and 1970, increasing elsewhere
# ============================================================
break_year_start <- 1945
break_year_end <- 1970

monot_mixed <- rep(0, length(knots) - 1)
for (i in 1:length(monot_mixed)) {
  interval_start <- knots_years[i]
  interval_end <- knots_years[i + 1]

  if (interval_start >= break_year_end) {
    monot_mixed[i] <- 1      # increasing after 1970
  } else if (interval_end <= break_year_start) {
    monot_mixed[i] <- 1      # increasing before 1945
  } else if (interval_start >= break_year_start && interval_end <= break_year_end) {
    monot_mixed[i] <- -1     # decreasing between 1945 and 1970
  } else if (interval_start < break_year_start && interval_end > break_year_start) {
    monot_mixed[i] <- 1      # increasing (partial interval before 1945)
  } else if (interval_start < break_year_end && interval_end > break_year_end) {
    monot_mixed[i] <- 1      # increasing (partial interval after 1970)
  }
}

cat("\n=== 4. Mixed constraints ===\n")
cat("   Decreasing between 1945 and 1970\n")
cat("   Increasing elsewhere\n")
cat("   Constraint vector:", monot_mixed, "\n")

fit_mixed <- SplineConstQuantRegBs3(xtab, ytab, knots, tau = 0.5,
                                    monot = monot_mixed, convcons = 0)
y_mixed <- spline_eval(fit_mixed, x_eval)

# ============================================================
# Model 5: Multiple quantiles with mixed constraints
# ============================================================
cat("\n=== 5. Multiple quantiles (0.1, 0.5, 0.9) with mixed constraints ===\n")
tau_multi <- c(0.1, 0.5, 0.9)
fits_multi <- list()
for (i in seq_along(tau_multi)) {
  fits_multi[[i]] <- SplineConstQuantRegBs3(xtab, ytab, knots,
                                            tau = tau_multi[i],
                                            monot = monot_mixed,
                                            convcons = 0)
  cat(sprintf("   tau = %.1f done\n", tau_multi[i]))
}

# Compute derivatives for analysis
compute_derivative <- function(y, x) {
  diff(y) / diff(x)
}

deriv_mixed <- compute_derivative(y_mixed, years_eval)
deriv_full <- compute_derivative(y_full_inc, years_eval)

# Find indices for analysis
idx_1945 <- which.min(abs(years_eval - 1945))
idx_1970 <- which.min(abs(years_eval - 1970))

cat("\n=== Derivative Analysis ===\n")
cat(sprintf("Full increasing:       min derivative = %.4f (should be >= 0)\n",
            min(deriv_full)))
cat(sprintf("Mixed constraints:     min derivative = %.4f\n", min(deriv_mixed)))
cat(sprintf("                       derivative between 1945-1970: range [%.4f, %.4f]\n",
            min(deriv_mixed[idx_1945:idx_1970]),
            max(deriv_mixed[idx_1945:idx_1970])))

# ============================================================
# Visualization
# ============================================================
par(mfrow = c(2, 2), mar = c(4, 4, 4, 3))

# Plot 1: Unconstrained fit
plot(years, temperature, pch = 16, cex = 0.6, col = "gray",
     xlab = "Year", ylab = "Temperature Anomaly (°C)",
     main = "1. Unconstrained Median (tau = 0.5)")
lines(years_eval, y_uncon, col = "red", lwd = 2)
abline(v = knots_years, col = "blue", lty = 2, lwd = 0.5)
abline(v = c(1945, 1970), col = "orange", lty = 3, lwd = 1.5)
abline(h = 0, col = "black", lty = 3, lwd = 0.5)
grid()

# Plot 2: Full monotonicity (too restrictive)
plot(years, temperature, pch = 16, cex = 0.6, col = "gray",
     xlab = "Year", ylab = "Temperature Anomaly (°C)",
     main = "2. Full Increasing Constraint (too restrictive)")
lines(years_eval, y_full_inc, col = "blue", lwd = 2)
abline(v = knots_years, col = "blue", lty = 2, lwd = 0.5)
abline(v = c(1945, 1970), col = "orange", lty = 3, lwd = 1.5)
abline(h = 0, col = "black", lty = 3, lwd = 0.5)
text(1948, -0.1, "Forced increasing\nbut data shows cooling", col = "blue", cex = 0.7)
grid()

# Plot 3: Mixed constraints (decreasing 1945-1970, increasing elsewhere)
plot(years, temperature, pch = 16, cex = 0.6, col = "gray",
     xlab = "Year", ylab = "Temperature Anomaly (°C)",
     main = "3. Mixed Constraints")
lines(years_eval, y_mixed, col = "darkgreen", lwd = 2)
abline(v = knots_years, col = "blue", lty = 2, lwd = 0.5)
abline(v = c(1945, 1970), col = "orange", lty = 2, lwd = 2)
abline(h = 0, col = "black", lty = 3, lwd = 0.5)

# Add shaded region for decreasing constraint
rect(1945, -0.8, 1970, 0.6, col = rgb(1, 0.5, 0, 0.1), border = NA)
text(1957, 0.5, "Decreasing constraint", col = "orange", cex = 0.8)
text(1900, 0.4, "Increasing", col = "darkgreen", cex = 0.7)
text(1985, 0.4, "Increasing", col = "darkgreen", cex = 0.7)
grid()

# Plot 4: Multiple quantiles with mixed constraints
colors <- c("orange", "red", "darkred")
plot(years, temperature, pch = 16, cex = 0.6, col = "gray",
     xlab = "Year", ylab = "Temperature Anomaly (°C)",
     main = "4. Quantile Regression with Mixed Constraints")
for (i in seq_along(tau_multi)) {
  y_fit <- spline_eval(fits_multi[[i]], x_eval)
  lines(years_eval, y_fit, col = colors[i], lwd = 2,
        lty = ifelse(tau_multi[i] == 0.5, 1, 2))
}
abline(v = knots_years, col = "blue", lty = 2, lwd = 0.5)
abline(v = c(1945, 1970), col = "orange", lty = 2, lwd = 2)
abline(h = 0, col = "black", lty = 3, lwd = 0.5)
rect(1945, -0.8, 1970, 0.6, col = rgb(1, 0.5, 0, 0.1), border = NA)
legend("topleft", legend = c("tau = 0.1", "tau = 0.5", "tau = 0.9"),
       col = colors, lty = c(2, 1, 2), lwd = 2, cex = 0.8)
grid()

# ============================================================
# Conclusion
# ============================================================
cat("\n========================================\n")
cat("Conclusion\n")
cat("========================================\n")
cat("The temperature data shows:\n")
cat("  - A general warming trend, especially after 1970\n")
cat("  - A cooling period between 1945 and 1970\n")
cat("  - Full increasing constraint ignores this cooling period\n")
cat("  - Mixed constraints (decreasing 1945-1970, increasing elsewhere)\n")
cat("    better reflect the actual pattern\n")
cat("  - Quantile regression reveals changing distribution over time\n")
cat("\nThis demonstrates that flexible constraint patterns\n")
cat("are more realistic than global monotonicity.\n")
cat("\nDemo completed.\n")

