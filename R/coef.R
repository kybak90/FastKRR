#' Coef method for fitted Kernel Ridge Regression models
#'
#' @description
#' Extracts the estimated coefficients from a fitted Kernel Ridge Regression (KRR) model.
#' The type of coefficient reported depends on the kernel approximation method:
#' for \code{opt = "exact"}, \code{"nystrom"}, or \code{"pivoted"}, the coefficients represent
#' the dual weights \eqn{\alpha}; for \code{opt = "rff"}, they represent the primal weights \eqn{\beta}.
#'
#' @param object An S3 object of class \code{krr}, typically returned by
#'   \code{\link{fastkrr}}.
#' @param ... Additional arguments (currently ignored).
#'
#' @return A numeric vector of the estimated model coefficients (\eqn{\alpha} or \eqn{\beta}).
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
#'                  kernel = "gaussian", opt = "exact",
#'                  rho = rho, lambda = 1e-4)
#'
#' coef(model)
#'
#' @rdname coef.krr
#' @method coef krr
#' @export
coef.krr = function(object, ...){
  return(as.vector(object$coefficients))
}
