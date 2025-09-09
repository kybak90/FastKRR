
<!-- README.md is generated from README.Rmd. Please edit that file -->

# FastKRR

<!-- badges: start -->
<!-- badges: end -->

The ‘FastKRR’ implements its core computational operations in C++ via
‘RcppArmadillo’, enabling faster performance than pure R, improved
numerical stability, and parallel execution with OpenMP where available.
On systems without OpenMP support, the package automatically falls back
to single-threaded execution with no user configuration required. For
efficient model selection, it integrates with ‘CVST’ to provide
sequential-testing cross-validation that identifies competitive
hyperparameters without exhaustive grid search. The package offers a
unified interface for exact kernel ridge regression and three scalable
approximations—Nyström, Pivoted Cholesky, and Random Fourier
Features—allowing analyses with substantially larger sample sizes than
are feasible with exact KRR. It also integrates with the ‘tidymodels’
ecosystem via the ‘parsnip’ model specification ‘krr_reg’, the S3 method
‘tunable.krr_reg()’, and the direct fitting helper ‘fit_krr()’.

**Dependencies:** Rcpp, RcppArmadillo, CVST, parsnip  
This package uses **CVST** (GPL ≥ 2). Overall license: **GPL (≥ 2)**.

- Authorss:
  - Gyeongmin Kim, Sungshin Women’s University,
    <rlarudals0824@gmail.com>
  - Seyoung Lee, Sungshin Women’s University, <sudang0404@gmail.com>
  - Miyoung Jang, Sungshin Women’s University, <miyoung9072@gmail.com>
  - Kwan-Young Bak, professor at Sungshin Women’s University,
    <kybak@sungshin.ac.kr>,
    [ORCID:0000-0002-4541-160X](https://orcid.org/0000-0002-4541-160X%7D%7BORCID:0000-0002-4541-160X)

## Installation

You can install the development version of FastKRR from
[GitHub](https://github.com/kybak90/FastKRR) with:

``` r
# install.packages("pak")
pak::pak("kybak90/FastKRR")
#> ! Using bundled GitHub PAT. Please add your own PAT using `gitcreds::gitcreds_set()`.
#> ℹ Loading metadata database
#> ✔ Loading metadata database ... done
#> 
#> 
#> → Will install 36 packages.
#> → Will update 1 package.
#> → Will download 33 CRAN packages (45.53 MB), cached: 3 (2.53 MB).
#> → Will download 1 package with unknown size.
#> + cli                   3.6.5    🔧 ⬇ (1.46 MB)
#> + CVST                  0.2-3     ⬇ (85.44 kB)
#> + dplyr                 1.1.4    🔧 ⬇ (1.60 MB)
#> + farver                2.1.2    🔧 ⬇ (1.97 MB)
#> + FastKRR       0.1.0 → 0.1.0    👷🏻‍♂️🔧 ⬇ (GitHub: 7b27019)
#> + generics              0.1.4     ⬇ (81.27 kB)
#> + ggplot2               3.5.2     ⬇ (4.96 MB)
#> + globals               0.17.0    ⬇ (126.67 kB)
#> + glue                  1.8.0    🔧 ⬇ (172.79 kB)
#> + gtable                0.3.6     ⬇ (224.21 kB)
#> + hardhat               1.4.1    
#> + isoband               0.2.7    🔧 ⬇ (1.87 MB)
#> + kernlab               0.9-33   🔧 ⬇ (2.31 MB)
#> + labeling              0.4.3     ⬇ (60.95 kB)
#> + lifecycle             1.0.4     ⬇ (124.48 kB)
#> + magrittr              2.0.3    🔧 ⬇ (232.43 kB)
#> + parsnip               1.3.2    
#> + pillar                1.11.0    ⬇ (658.36 kB)
#> + pkgconfig             2.0.3     ⬇ (18.21 kB)
#> + prettyunits           1.2.0     ⬇ (156.51 kB)
#> + purrr                 1.0.4    🔧 ⬇ (560.62 kB)
#> + R6                    2.6.1     ⬇ (86.61 kB)
#> + RColorBrewer          1.1-3     ⬇ (53.10 kB)
#> + Rcpp                  1.1.0    🔧 ⬇ (3.36 MB)
#> + RcppArmadillo         14.6.0-1 🔧 ⬇ (1.69 MB)
#> + rlang                 1.1.6    🔧 ⬇ (1.88 MB)
#> + scales                1.4.0     ⬇ (867.77 kB)
#> + sparsevctrs           0.3.4    🔧
#> + stringi               1.8.7    🔧 ⬇ (14.77 MB)
#> + stringr               1.5.1     ⬇ (312.86 kB)
#> + tibble                3.3.0    🔧 ⬇ (688.24 kB)
#> + tidyr                 1.3.1    🔧 ⬇ (1.32 MB)
#> + tidyselect            1.2.1     ⬇ (224.10 kB)
#> + utf8                  1.2.6    🔧 ⬇ (209.04 kB)
#> + vctrs                 0.6.5    🔧 ⬇ (1.88 MB)
#> + viridisLite           0.4.2     ⬇ (1.30 MB)
#> + withr                 3.0.2     ⬇ (222.20 kB)
#> ℹ Getting 33 pkgs (45.53 MB) and 1 pkg with unknown size, 3 (2.53 MB) cached
#> ✔ Cached copy of FastKRR 0.1.0 (source) is the latest build
#> ✔ Got RColorBrewer 1.1-3 (aarch64-apple-darwin20) (53.26 kB)
#> ✔ Got CVST 0.2-3 (aarch64-apple-darwin20) (85.69 kB)
#> ✔ Got R6 2.6.1 (aarch64-apple-darwin20) (86.61 kB)
#> ✔ Got cli 3.6.5 (aarch64-apple-darwin20) (1.46 MB)
#> ✔ Got RcppArmadillo 14.6.0-1 (aarch64-apple-darwin20) (1.69 MB)
#> ✔ Got gtable 0.3.6 (aarch64-apple-darwin20) (224.01 kB)
#> ✔ Got isoband 0.2.7 (aarch64-apple-darwin20) (1.87 MB)
#> ✔ Got pkgconfig 2.0.3 (aarch64-apple-darwin20) (18.36 kB)
#> ✔ Got utf8 1.2.6 (aarch64-apple-darwin20) (209.04 kB)
#> ✔ Got Rcpp 1.1.0 (aarch64-apple-darwin20) (3.36 MB)
#> ✔ Got kernlab 0.9-33 (aarch64-apple-darwin20) (2.31 MB)
#> ✔ Got glue 1.8.0 (aarch64-apple-darwin20) (172.91 kB)
#> ✔ Got lifecycle 1.0.4 (aarch64-apple-darwin20) (124.40 kB)
#> ✔ Got generics 0.1.4 (aarch64-apple-darwin20) (81.27 kB)
#> ✔ Got scales 1.4.0 (aarch64-apple-darwin20) (867.77 kB)
#> ✔ Got purrr 1.0.4 (aarch64-apple-darwin20) (560.62 kB)
#> ✔ Got tibble 3.3.0 (aarch64-apple-darwin20) (688.24 kB)
#> ✔ Got withr 3.0.2 (aarch64-apple-darwin20) (222.42 kB)
#> ✔ Got viridisLite 0.4.2 (aarch64-apple-darwin20) (1.30 MB)
#> ✔ Got magrittr 2.0.3 (aarch64-apple-darwin20) (232.44 kB)
#> ✔ Got tidyr 1.3.1 (aarch64-apple-darwin20) (1.32 MB)
#> ✔ Got vctrs 0.6.5 (aarch64-apple-darwin20) (1.89 MB)
#> ✔ Got globals 0.17.0 (aarch64-apple-darwin20) (126.67 kB)
#> ✔ Got tidyselect 1.2.1 (aarch64-apple-darwin20) (224.10 kB)
#> ✔ Got stringr 1.5.1 (aarch64-apple-darwin20) (312.86 kB)
#> ✔ Got labeling 0.4.3 (aarch64-apple-darwin20) (61.27 kB)
#> ✔ Got prettyunits 1.2.0 (aarch64-apple-darwin20) (156.34 kB)
#> ✔ Got rlang 1.1.6 (aarch64-apple-darwin20) (1.88 MB)
#> ✔ Got pillar 1.11.0 (aarch64-apple-darwin20) (658.36 kB)
#> ✔ Got ggplot2 3.5.2 (aarch64-apple-darwin20) (4.96 MB)
#> ✔ Got dplyr 1.1.4 (aarch64-apple-darwin20) (1.60 MB)
#> ✔ Got farver 2.1.2 (aarch64-apple-darwin20) (1.97 MB)
#> ✔ Got stringi 1.8.7 (aarch64-apple-darwin20) (14.77 MB)
#> ✔ Installed FastKRR 0.1.0 (github::kybak90/FastKRR@7b27019) (69ms)
#> ✔ Installed CVST 0.2-3  (74ms)
#> ✔ Installed R6 2.6.1  (82ms)
#> ✔ Installed RColorBrewer 1.1-3  (86ms)
#> ✔ Installed cli 3.6.5  (95ms)
#> ✔ Installed dplyr 1.1.4  (106ms)
#> ✔ Installed RcppArmadillo 14.6.0-1  (154ms)
#> ✔ Installed Rcpp 1.1.0  (161ms)
#> ✔ Installed farver 2.1.2  (81ms)
#> ✔ Installed generics 0.1.4  (57ms)
#> ✔ Installed ggplot2 3.5.2  (36ms)
#> ✔ Installed globals 0.17.0  (59ms)
#> ✔ Installed glue 1.8.0  (59ms)
#> ✔ Installed gtable 0.3.6  (32ms)
#> ✔ Installed hardhat 1.4.1  (33ms)
#> ✔ Installed isoband 0.2.7  (33ms)
#> ✔ Installed kernlab 0.9-33  (31ms)
#> ✔ Installed labeling 0.4.3  (29ms)
#> ✔ Installed lifecycle 1.0.4  (30ms)
#> ✔ Installed magrittr 2.0.3  (31ms)
#> ✔ Installed parsnip 1.3.2  (33ms)
#> ✔ Installed pillar 1.11.0  (33ms)
#> ✔ Installed pkgconfig 2.0.3  (30ms)
#> ✔ Installed prettyunits 1.2.0  (55ms)
#> ✔ Installed purrr 1.0.4  (32ms)
#> ✔ Installed rlang 1.1.6  (62ms)
#> ✔ Installed scales 1.4.0  (63ms)
#> ✔ Installed sparsevctrs 0.3.4  (32ms)
#> ✔ Installed stringr 1.5.1  (18ms)
#> ✔ Installed stringi 1.8.7  (75ms)
#> ✔ Installed tibble 3.3.0  (40ms)
#> ✔ Installed tidyr 1.3.1  (34ms)
#> ✔ Installed tidyselect 1.2.1  (32ms)
#> ✔ Installed utf8 1.2.6  (32ms)
#> ✔ Installed vctrs 0.6.5  (57ms)
#> ✔ Installed viridisLite 0.4.2  (45ms)
#> ✔ Installed withr 3.0.2  (18ms)
#> ✔ 1 pkg + 42 deps: upd 1, added 36, dld 33 (45.53 MB) [12.4s]
```

<!-- ## Parallelization -->
<!-- Some functions in **FastKRR** support parallel computation via OpenMP. -->
<!-- - On Windows and most Linux systems, OpenMP is available by default, and computations will use multiple threads. -->
<!-- - On macOS, if OpenMP is not installed, the package will run in single-threaded mode automatically. -->
<!-- No special installation steps are required. -->

## Example

This is a basic example of fitting a Gaussian kernel ridge regression
estimator to a dataset ${(x_i, y_i)}_{i=1}^n$ using (1) exact
computation, (2) the Nyström approximation, (3) pivoted Cholesky
decomposition, and (4) random Fourier features.

``` r
library(FastKRR)

# example data set
set.seed(1)
n = 1000; d = 1
rho = 1
X = matrix(runif(n*d, 0, 1), nrow = n, ncol = d)
y = as.vector(sin(2*pi*rowMeans(X)^3) + rnorm(n, 0, 0.1))

# model fitting - exact
model_exact = fastkrr(X, y, kernel = "gaussian", rho = rho, opt = "exact", verbose = FALSE)

# model fitting - pivoted
model_pivoted = fastkrr(X, y, kernel = "gaussian", rho = rho, opt = "pivoted", verbose = FALSE)

# model fitting - nystrom
model_nystrom = fastkrr(X, y, kernel = "gaussian", rho = rho, opt = "nystrom", verbose = FALSE)

# model fitting - rff
model_rff = fastkrr(X, y, kernel = "gaussian", rho = rho, opt = "rff", verbose = FALSE)


# predict
new_n = 500
new_x = matrix(runif(new_n*d, 0, 1), nrow = new_n, ncol = d)
new_y = as.vector(sin(2*pi*rowMeans(new_x)^3) + rnorm(new_n, 0, 0.1))

pred_exact = pred_krr(model_exact, new_x)
pred_pivoted = pred_krr(model_pivoted, new_x)
pred_nystrom = pred_krr(model_nystrom, new_x)
pred_rff = pred_krr(model_rff, new_x)
```

The visualization of the fitted results is shown below.

    #> Warning: package 'ggplot2' was built under R version 4.3.3

<img src="man/figures/README-unnamed-chunk-3-1.png" width="100%" />
