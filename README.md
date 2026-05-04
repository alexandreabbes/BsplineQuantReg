# BsplineQuantReg
## Citation

If you use this package in your research, please cite:

```bibtex
@Article{Abbes2026,
  author  = {Alexandre Abbes},
  title   = {Quantile regression with cubic polynomial splines under shape constraints with applications},
  year    = {2025},
  doi     = {10.5281/zenodo.16999784}
}

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.16999784.svg)](https://doi.org/10.5281/zenodo.16999784)

## Beta Version Notice

This package is currently in **beta** (version 0.1.0-beta). 

### What to expect:
- ✅ Core functionality (quantile regression, monotonicity, convexity) is stable
- ⚠️ API may change based on user feedback
- 📝 Documentation is being refined
- 🐛 Please report issues on GitHub

## Related Packages

### R Packages

| Package | Description | Constraint Type | Spline Degree |
|---------|-------------|-----------------|---------------|
| **ConstrainedQuantileSplines** (this package) | Quantile regression with Karlin-Studden constraints | Monotonicity, Convexity | Cubic |
| [quantreg](https://CRAN.R-project.org/package=quantreg) | Classical quantile regression | None (linear programming) | Linear |
| [cobs](https://CRAN.R-project.org/package=cobs) | Constrained B-splines | Monotonicity, Convexity | Linear, Quadratic |

### Comparison with cobs

The `cobs` package (Constrained B-Splines) is the closest to this package, but with key differences:

| Feature | ConstrainedQuantileSplines | cobs |
|---------|---------------------------|------|
| Spline degree | Cubic (degree 3) | Linear, Quadratic |
| Constraint method | Karlin-Studden SOCP | Traditional constraints |
| Convexity | ✅ Yes | ✅ Yes |
| Monotonicity | ✅ Yes | ✅ Yes |
| Quantile regression | ✅ Yes | ✅ Yes |
| Partial constraints | ✅ Yes (per interval) | Limited |
| Polynomial coefficient export | ✅ Yes | ❌ No |

### When to use this package vs cobs

- **Use `cobs`** : For linear or quadratic splines, simpler constraints
- **Use this package** : For cubic splines testing , Karlin-Studden exact constraints, polynomial coefficient extraction

## ⚠️ Performance Notice

**This R package is currently intended for demonstration, prototyping,  and educational purposes only.**

Due to the current implementation (pure R with CVXR), the package is **significantly slower** than its Python counterpart. Cubic B-spline quantile regression with constraints involves solving SOCP problems, and the R implementation does not yet leverage optimized linear algebra libraries.
[Python version](https://github.com/alexandreabbes/Constrained-Quantile-Regression-with-cubic-splines) 

### Future Improvements

We plan to improve performance in future releases by:
- Linking with faster optimization libraries (OSQP, Gurobi)
- Implementing more efficient SOCP solvers
- Optimizing the B-spline basis computation

**The Python version remains the recommended choice for production use.**

### Installation for beta testing:

```r
# Install from GitHub
devtools::install_github("alexandreabbes/Constrained-Quantile-Regression-with-cubic-splines", ref = "R")
