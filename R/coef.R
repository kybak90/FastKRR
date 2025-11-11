#' Coef method for fitted Kernel Ridge Regression models
#'
#' @description
#' Displays the main coefficient information from a fitted Kernel Ridge Regression (KRR) model,
#' including the original function call and the first few estimated coefficients.
#' The type of coefficient reported depends on the kernel approximation method:
#' for \code{opt = "exact"}, \code{"nystrom"}, or \code{"pivoted"}, the coefficients represent
#' \eqn{\alpha}; for \code{opt = "rff"}, they represent the coefficient (\eqn{\beta}).
#'
#' @param object An S3 object of class \code{krr}, typically returned by
#'   \code{\link{fastkrr}}.
#' @param ... Additional arguments (currently ignored).
#'
#' @return A human-readable summary of the fitted KRR model to the console.
#'
#' @seealso \code{\link{fastkrr}}
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
#'
#' @rdname coef.krr
#' @method coef krr
#' @export
coef.krr = function(object, ...){
  model = object
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
