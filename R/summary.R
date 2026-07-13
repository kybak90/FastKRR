#' Summary method for fitted Kernel Ridge Regression models
#'
#' @description
#' Displays key information from a fitted Kernel Ridge Regression (KRR) model,
#' including the original call, first few coefficients, a 6×6 block of the
#' kernel (or approximated kernel) matrix, and the main kernel options.
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
#' X = matrix(runif(n*d, 0, 1), nrow = n, ncol = d); colnames(X) = paste0("X", seq_len(d))
#' y = sin(2 * pi * rowMeans(X)^3) + rnorm(n, mean = 0, sd = 0.1)
#'
#' data = data.frame(X, y = y)
#'
#' # Example: exact
#' model = fastkrr(data = data, response = "y",
#'                 kernel = "gaussian", opt = "exact",
#'                 rho = rho, selection_method = "fastCV")
#'
#' summary(model)
#'
#' @export
summary.krr = function(object, ...){
  model = object

  invisible(param(model))

  cat("Training error:\n")
  print(error(model))
}
