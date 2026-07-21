## Resubmission

This is a new version of an existing CRAN package (0.1.2 -> 0.2.0). This release includes major functionality enhancements, computational updates (e.g., REML), and interface improvements implemented during the revision of our paper for The R Journal, addressing feedback from the editors and reviewers.

## Changes in version 0.2.0

- Added a restricted maximum likelihood (REML) option (`selection_method = "REML"`) as an efficient, C++-implemented alternative to cross-validation-based tuning of the regularization parameter, available for both the exact method and all three low-rank approximations (Nystrom, pivoted Cholesky, RFF).
- `fastkrr()` now takes a `data.frame` as its primary input (via `data` and `response` arguments), instead of requiring separate `X`/`y` matrices.
- Added an `na.rm` argument to `fastkrr()` (and the tidymodels engine) to explicitly handle missing values, instead of silently producing `NA` outputs or an uninformative error.
- The `lambda` argument now accepts four input forms (`NULL`, a single value, a length-2 `(min, max)` vector for REML, or a grid of length >= 3 for cross-validation), with validation that depends on the chosen `selection_method`.
- Replaced attribute-based (`attr()`/`attributes()`) access to model components throughout the package with standard list-field access (e.g., `model$coefficients`, `model$fitted.values`), so that fitted model objects behave like ordinary R lists rather than requiring undocumented `attr()` calls. This affected `coef.krr()`, `error.krr()`, `param.krr()`, `plot.krr()`, `print.krr()`, and `summary.krr()`.
- Fixed `coef.krr()`, which previously only printed coefficients to the console and did not return a value; it now returns the coefficient vector as expected of a standard `coef()` method.
- Exposed the Cholesky/low-rank factor matrices directly as components of the returned model object (`chol_factor` for the exact method; `approx_factor` for the low-rank approximations), and added `K_approx` as a returned field exclusively in `approx_kernel()`.
- Added a `data_new` argument to `error()` to compute prediction MSE on held-out data, in addition to training MSE.
- Added a comprehensive `testthat` test suite covering the main fitting, prediction, and utility functions.
- Various documentation clarifications and manual updates.

## R CMD check results

0 errors | 0 warnings | 0 notes
