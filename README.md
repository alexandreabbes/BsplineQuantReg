# BsplineQuantReg

[![CRAN status](https://www.r-pkg.org/badges/version/BsplineQuantReg)](https://cran.r-project.org/package=BsplineQuantReg)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17427913.svg)](https://doi.org/10.5281/zenodo.17427913)
[![CRAN downloads](https://cranlogs.r-pkg.org/badges/BsplineQuantReg)](https://cran.r-project.org/package=BsplineQuantReg)

Constrained Quantile Regression with B-Splines (Degrees 1 to 4)

This package is available on CRAN.

```
## Features

- Quantile regression for any tau in (0,1)
- B-splines of degree 1 to 4 (linear to quartic)
- Monotonicity constraints (increasing or decreasing)
- Convexity constraints (convex or concave)
- Third derivative constraints for cubic and quartic splines
- Karlin-Studden SOCP formulation for rigorous shape constraints
- Partial constraints (per interval or per knot)
- Polynomial coefficient export
- Conversion between PP-form and B-spline representation
- Callable and non-callable spline objects



## Installation

### From CRAN (stable, recommended)

```r
install.packages("BsplineQuantReg")
```

### From GitHub (development version)

```r
devtools::install_github("alexandreabbes/BsplineQuantReg")
```

## System Requirements

### Linux Users

On Linux systems, the packages `CVXR` and `CLARABEL` require the Rust compiler and Cargo package manager to be installed.

#### Ubuntu/Debian:
```bash
sudo apt-get install cargo rustc
```

#### Fedora and other linux dist.

```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -
source ~/.cargo/env
```

#### Verify installation:

```bash
rustc --version
cargo --version
```

After installing Rust and Cargo, restart R and install the package:

```r
install.packages("BsplineQuantReg")
```

### Windows Users

Windows users do not need to install Rust separately. The package uses pre-compiled binaries available on CRAN.

## Graphical Interface (Shiny)

This package includes an interactive Shiny interface that allows you to manage most of the functions without writing code.

###Launch the GUI
```r
library(BsplineQuantRegGui)
run_gui()
```

The interface will open in your default browser.
###Features of the GUI

    Load data (CSV, Excel, built-in datasets)

    Configure spline parameters (degree, knots)

    Apply shape constraints (monotonicity, convexity, third derivative)

    Define constraints per region interactively

    Run quantile regression with various solvers

    Visualize results with interactive plots

    Export reproducible R code

    Run built-in demos

###Docker Deployment
A Docker image is available for easy deployment:
bash
bash```
docker pull ghcr.io/alexandreabbes/bsplinequantreggui:latest
docker run -p 3838:3838 ghcr.io/alexandreabbes/bsplinequantreggui:latest
```
Then open http://localhost:3838 in your browser.


### Installation

```r
# Install from GitHub
remotes::install_github("alexandreabbes/bsplinequantreggui")

# Or with pak
pak::pkg_install("alexandreabbes/bsplinequantreggui")
```



## Related Packages

### R Packages

| Package | Description | Constraint Type | Spline Degree |
|---------|-------------|-----------------|---------------|
| BsplineQuantReg (this package) | Quantile regression with Karlin-Studden constraints | Monotonicity, Convexity | 1 to 4 |
| quantreg | Classical quantile regression | None (linear programming) | Linear |
| cobs | Constrained B-splines | Monotonicity, Convexity | Linear, Quadratic |

### Comparison with cobs

The `cobs` package (Constrained B-Splines with linear or quadratic splines)
is the closest to this package.

### Python package
Python version: https://pypi.org/project/BsplineQuantRegpy/

### Matlab core
You may take a look here for a glimpse at the initial matlab code.
https://github.com/alexandreabbes/Constrained-Quantile-Regression-with-cubic-splines/tree/matlab
`
## Performance Notice
This R package is intended for demonstration, prototyping, and educational purposes.
Due to the current implementation (pure R with CVXR) the package is  almost 5 times slower than its Python counterpart (benchmark test), but much (10 times 
maybe) faster than the matlab one.
B-spline quantile regression with constraints involves solving SOCP problems, and the R implementation does not yet leverage optimized linear algebra libraries.


### Future Improvements
- Optimize the B-spline basis computation
- Improve the API based on user feedback

## Getting Started

```r
library(BsplineQuantReg)

# Generate sample data
set.seed(42)
x <- seq(0, 1, length.out = 100)
y <- 2*x + 0.5*sin(6*pi*x) + 0.05*rnorm(100)
knots <- quantile(x, probs = seq(0, 1, length.out = 10))

# Quantile regression with cubic spline and increasing constraint
fit <- SplineCubicQuant(x, y, knots, tau = 0.5, monot = 1)

# Evaluate the spline
x_eval <- seq(0, 1, length.out = 200)
#y_eval <- spline_eval(fit, x_eval) # deprecated now
y_eval <- fit(x_eval) # the fit is now callable
```

## Demos

```r
# List available demos
demo(package = "BsplineQuantReg")

# Run a specific demo
demo("comprehensive", package = "BsplineQuantReg")
demo("temperature", package = "BsplineQuantReg")
```


## Citation

If you use this package in your research, please cite:

```bibtex
@Article{Abbes2025,
  author  = {Alexandre Abbes},
  title   = {Constrained Quantile Regression with Cubic B-Splines under Shape Constraints},
  year    = {2025},
  doi     = {10.5281/zenodo.17427913}
}

## Bug Reports

Please report issues on GitHub: https://github.com/alexandreabbes/BsplineQuantReg/issues
```
