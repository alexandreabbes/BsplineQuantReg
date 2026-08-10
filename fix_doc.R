# verify_docs.R
verify_documentation <- function() {
  cat("=== Verifying BsplineQuantReg Documentation ===\n\n")

  # Check each function
  functions <- list(
    bspline_to_deriv_coeffs_cubic = list(
      expected = c("tn", "degree", "x_values", "Bsbasis", "verbose"),
      file = "man/bspline_to_deriv_coeffs_cubic.Rd"
    ),
    makpp = list(
      expected = c("coeff", "tn", "callable", "verbose"),
      file = "man/makpp.Rd"
    ),
    view_basis = list(
      expected = c("BB", "x_values", "main", "view_knot", "add_knot"),
      file = "man/view_basis.Rd"
    )
  )

  issues <- 0

  for (func_name in names(functions)) {
    cat("\nChecking", func_name, "...\n")

    # Get actual arguments
    if (exists(func_name)) {
      actual <- names(formals(func_name))
      cat("  Actual arguments:", paste(actual, collapse = ", "), "\n")

      expected <- functions[[func_name]]$expected
      cat("  Expected arguments:", paste(expected, collapse = ", "), "\n")

      # Check for mismatches
      if (!identical(actual, expected)) {
        cat("  ❌ MISMATCH!\n")
        missing <- setdiff(expected, actual)
        extra <- setdiff(actual, expected)
        if (length(missing) > 0) {
          cat("    Missing:", paste(missing, collapse = ", "), "\n")
        }
        if (length(extra) > 0) {
          cat("    Extra:", paste(extra, collapse = ", "), "\n")
        }
        issues <- issues + 1
      } else {
        cat("  ✅ OK\n")
      }
    } else {
      cat("  ⚠️ Function not found\n")
    }
  }

  cat("\n=== Summary ===\n")
  if (issues == 0) {
    cat("✅ All documentation matches function signatures!\n")
  } else {
    cat("❌ Found", issues, "documentation issue(s)\n")
  }
}

# Run verification
verify_documentation()
