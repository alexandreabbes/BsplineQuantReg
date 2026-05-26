# Demo: Temperature Anomaly Trend Analysis
# Data source: Global temperature anomalies (1880-1992)
#
# Constraints tested:
# 1. Unconstrained fit
# 2. Uniform increasing constraint (everywhere)
# 3. Mixed: Decreasing only between 1956-1973, unconstrained elsewhere

library(BsplineQuantReg)
oldpar <- par(mfrow = c(2,2))

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
cat(sprintf("Temperature range: [%.2f, %.2f]  C\n\n",
            min(temperature), max(temperature)))

# Use raw years (no normalization)
xtab <- years
ytab <- temperature

# More knots for better flexibility
knots_years <- c(1880, 1900, 1917, 1936, 1956, 1965, 1973, 1980, 1992)
knots <- knots_years

cat("Knots at years:", knots_years, "\n")
cat(sprintf("Number of intervals: %d\n\n", length(knots) - 1))

# Create evaluation grid
years_eval <- seq(min(years), max(years), length.out = 300)

# ============================================================
# Model 1: Unconstrained median regression
# ============================================================
cat("=== 1. Unconstrained median regression ===\n")
fit_uncon <- SplineConstQuantRegBs3(xtab, ytab, knots, tau = 0.5,
                                    monot = 0, convcons = 0)
y_uncon <- spline_eval(fit_uncon, years_eval)

# ============================================================
# Model 2: Uniform increasing constraint (everywhere)
# ============================================================
cat("\n=== 2. Uniform increasing constraint (everywhere) ===\n")
fit_uniform_inc <- SplineConstQuantRegBs3(xtab, ytab, knots, tau = 0.5,
                                          monot = 1, convcons = 0)
y_uniform_inc <- spline_eval(fit_uniform_inc, years_eval)

# ============================================================
# Model 3: Mixed constraints
# Decreasing only between 1956-1973, unconstrained elsewhere
# ============================================================
decr_start <- 1956
decr_end <- 1973

# Build constraint vector for each interval
monot_mixed <- rep(0, length(knots) - 1)  # 0 = unconstrained

for (i in 1:length(monot_mixed)) {
  interval_start <- knots[i]
  interval_end <- knots[i + 1]

  # Interval is entirely within 1956-1973
  if (interval_start >= decr_start && interval_end <= decr_end) {
    monot_mixed[i] <- -1  # decreasing
  }
  # Interval partially overlaps - apply decreasing if mostly inside
  else if (interval_start < decr_start && interval_end > decr_start &&
           interval_end <= decr_end) {
    # Interval starts before, ends inside cooling period
    monot_mixed[i] <- -1  # decreasing
  }
  else if (interval_start >= decr_start && interval_start < decr_end &&
           interval_end > decr_end) {
    # Interval starts inside, ends after cooling period
    monot_mixed[i] <- -1  # decreasing
  }
  else if (interval_start < decr_start && interval_end > decr_end) {
    # Interval spans the entire cooling period
    monot_mixed[i] <- -1  # decreasing
  }
  # All other intervals: unconstrained (0)
}

cat("\n=== 3. Mixed constraints ===\n")
cat(sprintf("   Decreasing only between %d and %d\n", decr_start, decr_end))
cat("   Unconstrained elsewhere\n\n")
cat("   Constraint vector by interval:\n")
for (i in 1:length(monot_mixed)) {
  constraint_type <- ifelse(monot_mixed[i] == 1, "increasing",
                            ifelse(monot_mixed[i] == -1, "decreasing", "unconstrained"))
  cat(sprintf("     Interval %d (%d-%d): %s\n",
              i, knots[i], knots[i+1], constraint_type))
}

fit_mixed <- SplineConstQuantRegBs3(xtab, ytab, knots, tau = 0.5,
                                    monot = monot_mixed, convcons = 0)
y_mixed <- spline_eval(fit_mixed, years_eval)

# ============================================================
# Model 4: Multiple quantiles with mixed constraints
# ============================================================
cat("\n=== 4. Multiple quantiles (0.1, 0.5, 0.9) with mixed constraints ===\n")
tau_multi <- c(0.1, 0.5, 0.9)
fits_multi <- list()
for (i in seq_along(tau_multi)) {
  fits_multi[[i]] <- SplineConstQuantRegBs3(xtab, ytab, knots,
                                            tau = tau_multi[i],
                                            monot = monot_mixed,
                                            convcons = 0)
  cat(sprintf("   tau = %.1f done\n", tau_multi[i]))
}

# ============================================================
# Visualization
# ============================================================
par(mfrow = c(2, 2), mar = c(4, 4, 4, 3))

# Plot 1: Unconstrained fit
plot(years, temperature, pch = 16, cex = 0.6, col = "gray",
     xlab = "Year", ylab = "Temperature Anomaly (C)",
     main = "1. Unconstrained Median (tau = 0.5)")
lines(years_eval, y_uncon, col = "red", lwd = 2)
abline(v = knots, col = "blue", lty = 2, lwd = 0.5)
abline(v = c(decr_start, decr_end), col = "orange", lty = 2, lwd = 1.5)
abline(h = 0, col = "black", lty = 3, lwd = 0.5)
grid()

# Plot 2: Uniform increasing constraint
plot(years, temperature, pch = 16, cex = 0.6, col = "gray",
     xlab = "Year", ylab = "Temperature Anomaly (C)",
     main = "2. Uniform Increasing Constraint (everywhere)")
lines(years_eval, y_uniform_inc, col = "blue", lwd = 2)
abline(v = knots, col = "blue", lty = 2, lwd = 0.5)
abline(v = c(decr_start, decr_end), col = "orange", lty = 2, lwd = 1.5)
abline(h = 0, col = "black", lty = 3, lwd = 0.5)
grid()

# Plot 3: Mixed constraints (decreasing only 1956-1973)
plot(years, temperature, pch = 16, cex = 0.6, col = "gray",
     xlab = "Year", ylab = "Temperature Anomaly (C)",
     main = "3. Mixed Constraints (decreasing 1956-1973 only)")
lines(years_eval, y_mixed, col = "darkgreen", lwd = 2)
abline(v = knots, col = "blue", lty = 2, lwd = 0.5)
abline(v = c(decr_start, decr_end), col = "orange", lty = 2, lwd = 2)
abline(h = 0, col = "black", lty = 3, lwd = 0.5)

# Shade the decreasing constraint region
rect(decr_start, -0.8, decr_end, 0.5, col = rgb(1, 0.5, 0, 0.1), border = NA)
text(1964, 0.4, "Decreasing\nconstraint", col = "orange", cex = 0.8)
grid()

# Plot 4: Multiple quantiles with mixed constraints
colors <- c("orange", "red", "darkred")
plot(years, temperature, pch = 16, cex = 0.6, col = "gray",
     xlab = "Year", ylab = "Temperature Anomaly (C)",
     main = "4. Quantile Regression with Mixed Constraints")
for (i in seq_along(tau_multi)) {
  y_fit <- spline_eval(fits_multi[[i]], years_eval)
  lines(years_eval, y_fit, col = colors[i], lwd = 2,
        lty = ifelse(tau_multi[i] == 0.5, 1, 2))
}
abline(v = knots, col = "blue", lty = 2, lwd = 0.5)
abline(v = c(decr_start, decr_end), col = "orange", lty = 2, lwd = 2)
abline(h = 0, col = "black", lty = 3, lwd = 0.5)
rect(decr_start, -0.8, decr_end, 0.5, col = rgb(1, 0.5, 0, 0.1), border = NA)
legend("topleft", legend = c("tau = 0.1", "tau = 0.5", "tau = 0.9"),
       col = colors, lty = c(2, 1, 2), lwd = 2, cex = 0.8)
grid()

par(oldpar)

# ============================================================
# Conclusion
# ============================================================
cat("\n========================================\n")
cat("Conclusion\n")
cat("========================================\n")
cat("The temperature data shows:\n")
cat("  - A general warming trend, especially after 1970\n")
cat("  - A cooling period between 1956 and 1973\n")
cat("  - Uniform increasing constraint forces monotonicity everywhere\n")
cat("    and may not capture the cooling period well\n")
cat("  - Mixed constraints (decreasing only 1956-1973, unconstrained elsewhere)\n")
cat("    allow the fit to reflect both the cooling period and the overall warming\n")
cat("\nDemo completed.\n")
