#' Coef method for fitted Kernel Ridge Regression models
#'
#' @description
#' Displays key information from a fitted Kernel Ridge Regression (KRR) model,
#' including the original call, first few coefficients, a 6×6 block of the
#' kernel (or approximated kernel) matrix, and the main kernel options.
#'
#' @param x An S3 object of class \code{krr}, typically returned by
#'   \code{\link{fastkrr}}.
#' @param ... Additional arguments (currently ignored).
#'
#' @return A human-readable summary of the fitted KRR model to the console.
#'
#' @seealso \code{\link{fastkrr}}, \code{\link{print.approx_kernel}},
#'   \code{\link{print.kernel_matrix}}
#'
#' @examples
#' # Data setting
#' set.seed(1)
#' lambda = 1e-4
#' d = 1
#' n = 50
#' rho = 1
#' X = matrix(runif(n*d, 0, 1), nrow = n, ncol = d)
#' y = as.vector(sin(2*pi*rowMeans(X)^3) + rnorm(n, 0, 0.1))
#'
#' # Example: exact
#' model = fastkrr(X, y,
#'                 kernel = "gaussian", opt = "exact",
#'                 rho = rho, lambda = 1e-4)
#' model = fastkrr(X, y,
#'                 kernel = "gaussian", opt = "nystrom",
#'                 rho = rho, lambda = 1e-4)
#' model = fastkrr(X, y,
#'                 kernel = "gaussian", opt = "rff",
#'                 rho = rho, lambda = 1e-4)
#' class(model)
#'
#' coef(model)
#' @export
coef.krr = function(x, ...){
  model = x
  cat("Call:\n")
  print(attributes(model)$call)
  cat("\n")

  cat("Coefficients:\n")
  coefs = as.vector(attr(model, "coefficients"))
  n_show = min(6, length(coefs))
  cat("  ")
  cat(coefs[1:n_show])
  if (length(coefs) > n_show) cat(" ...")
}
