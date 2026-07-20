#' Summary method for fitted Kernel Ridge Regression models
#'
#' @description
#' Computes and displays a comprehensive summary of a fitted Kernel Ridge Regression (KRR) model,
#' including the original function call, fitted hyperparameters (via \code{\link{param}}),
#' and the final training mean squared error (via \code{\link{error}}).
#'
#' @param object An S3 object of class \code{krr}, typically returned by
#'   \code{\link{fastkrr}}.
#' @param ... Additional arguments (currently ignored).
#'
#' @return Invisibly returns the hyperparameter list or model statistics after
#'   printing the summary panel to the console.
#'
#' @seealso \code{\link{fastkrr}}, \code{\link{param.krr}}, \code{\link{error.krr}}
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
#'                  kernel = "gaussian", opt = "exact",
#'                  rho = rho, lambda = lambda)
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
