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
- **Use this package** : For cubic splines, Karlin-Studden exact constraints, polynomial coefficient extraction

### Installation for beta testing:

```r
# Install from GitHub
devtools::install_github("alexandreabbes/Constrained-Quantile-Regression-with-cubic-splines", ref = "R")
