#' Compute Model Error (Generic)
#'
#' @description
#' Generic function for computing model error.
#'
#' @param x An object to compute model error for.
#' @param ... Additional arguments passed to methods.
#'
#' @return A numeric value or class-specific result.
#' @export
error = function(x, ...) {
  UseMethod("error")
}

#' Compute Model Error for Kernel Ridge Regression Models
#'
#' @description
#' Computes the model error for kernel ridge regression (\code{"krr"}) objects.
#' Returns the mean squared error (MSE) between the observed responses
#' and the fitted values stored in the object.
#'
#' @param x An object of class \code{"krr"}.
#' @param ... Additional arguments (ignored).
#'
#' @return A numeric value giving the mean squared error (MSE).
#'
#' @details
#' This method computes:
#' \deqn{\text{MSE} = \frac{1}{n} \sum_i (y_i - \hat{y}_i)^2}
#' where \code{y} and \code{fitted.values} are stored in the \code{"krr"} object attributes.
#'
#' @seealso \code{\link{summary.krr}}, \code{\link{plot.krr}}, \code{\link{predict.krr}}
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
#' model = fastkrr(X, y, kernel = "gaussian", lambda = 0.001)
#' error(model)
#'
#'
#' @export
error.krr = function(x, ...) {
  model = x
  as.numeric(
    crossprod(attributes(model)$y - attributes(model)$fitted.values) /
      length(attr(model, "y"))
  )
}
