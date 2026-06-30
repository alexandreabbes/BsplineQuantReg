# BsplineQuantReg
## Citation

If you use this package in your research, please cite:

```bibtex
@Article{Abbes2026,
  author  = {Alexandre Abbes},
  title   = {Quantile regression with cubic polynomial splines under shape constraints with applications},
  year    = {2026},
  doi     = {10.5281/zenodo.17427913}
}

# BsplineQuantReg

[![CRAN status](https://www.r-pkg.org/badges/version/BsplineQuantReg)](https://cran.r-project.org/package=BsplineQuantReg)
[![DOI](https://doi.org/10.5281/zenodo.17427913)](https://doi.org/10.5281/zenodo.17427913)
[![CRAN downloads](https://cranlogs.r-pkg.org/badges/BsplineQuantReg)](https://cran.r-project.org/package=BsplineQuantReg)

**Constrained Quantile Regression with Cubic B-Splines**

This package is now available on CRAN! 🎉

### What to expect:



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
- Optimizing the B-spline basis computation
- Implement the quartic version already available in Python.
- Improve the API based on user feedback


🐛 Please report issues on GitHub

**The Python version remains the recommended choice for production use.**

## Installation

### From CRAN (stable, recommended)
```r
install.packages("BsplineQuantReg")

### Install from GitHub (development version)
pak::pak("alexandreabbes/BsplineQuantReg")
# or
devtools::install_github("alexandreabbes/BsplineQuantReg")

###System Requirements
##Linux Users

##On Linux systems, the packages CVXR and CLARABEL require the Rust compiler and Cargo package manager to be installed.
#Ubuntu/Debian:
# Install Rust and Cargo
```bash (sudo or root)
apt-get install rust
apt-get install cargo
