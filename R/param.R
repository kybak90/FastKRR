#' Extract/print hyperparameters of fitted models
#'
#' @description
#' \code{param} is a generic S3 function that displays (and invisibly returns)
#' model hyperparameters. Methods are provided for \code{"krr"} objects.
#'
#' @param x An object.
#' @param ... Additional arguments passed to methods.
#' @return A named list of hyperparameters (invisibly); may print side effects.
#' @keywords internal
#' @export
param = function(x, ...) UseMethod("param")

#' @rdname param
#' @export
param.default = function(x, ...) {
  stop(sprintf("No 'param' method for objects of class: %s", paste(class(x), collapse=", ")))
}

#' Displays hyperparameters of fitted Kernel Ridge Regression models
#'
#' @description
#' Displays (and invisibly returns) the hyperparameters actually used by a
#' fitted object. For \code{krr} objects returned by \code{\link{fastkrr}},
#' this prints a concise hyperparameter panel (e.g., \code{rho}, \code{lambda},
#' \code{m}, \code{eps}, \code{n_threads}, \code{d}).
#'
#' @details
#' \strong{Pivoted approximation note}: When \code{opt = "pivoted"}, the
#' effective number of pivots \code{m} used during the approximation may be
#' smaller than the user-specified \code{m} because the algorithm can stop
#' early based on \code{eps}. If you want to confirm the initial \code{m}
#' that you set, please see the printed \emph{Call} (the original function
#' call shows your input arguments).
#'
#' @param x An object of class \code{"krr"}, typically returned by
#'   \code{\link{fastkrr}}.
#' @param ... Additional arguments.
#'
#' @return
#' Prints a human-readable panel to the console and invisibly returns a
#' named list of class \code{"krr_params"} containing the extracted hyperparameters:
#' \itemize{
#'   \item \code{kernel}: Kernel type used ("gaussian" or "laplace").
#'   \item \code{opt}: Kernel approximation method.
#'   \item \code{selection_method}: Lambda tuning method.
#'   \item \code{rho}: Kernel scaling parameter.
#'   \item \code{lambda}: Regularization parameter.
#'   \item \code{m}: Rank or number of random features used for approximation.
#'   \item \code{eps}: Tolerance parameter (for pivoted Cholesky).
#'   \item \code{n_threads}: Number of parallel threads used.
#' }
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
#' model = fastkrr(data = data, response = "y",
#'                  kernel="gaussian", opt="pivoted",
#'                  rho=1, lambda=lambda, n_threads = 1)
#'
#' # Inspect hyperparameters
#' param(model)
#'
#' @export
param.krr = function(x, ...) {
  model = x

  keys = c("kernel", "opt", "selection_method",
           "rho", "lambda", "m", "eps", "n_threads")
  present = intersect(keys, names(model))

  out = model[present]
  class(out) = c("krr_params", "list")

  # Call
  cat("Call:\n")
  print(model$call)
  cat("\n")

  # Parameters
  if (length(present) > 0) {
    df = as.data.frame(model[present], optional = TRUE)
    row.names(df) = ""
    cat("Parameters:\n")
    print(df)
    cat("\n")
  }

  invisible(out)
}
