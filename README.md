# FastKRR

<!-- badges: start -->
<!-- badges: end -->

The goal of FastKRR is to fit a Kernel Ridge regression estimator based on a set of data $(x_i, y_i)^n_{i = 1}$.

**Dependencies:** Rcpp, RcppArmadillo, CVST, parsnip  
This package uses **CVST** (GPL ≥ 2). Overall license: **GPL (≥ 3)**.

## Installation

You can install the development version of FastKRR from [GitHub](https://github.com/jang-miyoung-041/FastKRR) with:

``` r
# install.packages("pak")
pak::pak("jang-miyoung-041/FastKRR")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(FastKRR)
set.seed(1)

# example data set
n = 1000; d = 1
rho = 1
X = matrix(runif(n*d, 0, 1), nrow = n, ncol = d)
y = as.vector(sin(2*pi*rowMeans(X)^3) + rnorm(n, 0, 0.1))

# model fitting - exact
model_exact = fastkrr(X, y, kernel = "gaussian", rho = rho, opt = "exact")

# model fitting - pivoted
model_pivoted = fastkrr(X, y, kernel = "gaussian", rho = rho, opt = "pivoted")

# model fitting - nystrom
model_nystrom = fastkrr(X, y, kernel = "gaussian", rho = rho, opt = "nystrom")

# model fitting - rff
model_rff = fastkrr(X, y, kernel = "gaussian", rho = rho, opt = "rff")

# predict
new_n = 1500
new_x = matrix(runif(new_n*d, 0, 1), nrow = new_n, ncol = d)
new_y = as.vector(sin(2*pi*rowMeans(new_x)^3) + rnorm(new_n, 0, 0.1))
pred = pred_krr(model_exact, new_x)
crossprod(pred, new_y) / new_n
```

