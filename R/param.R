#' Extract/print hyperparameters of fitted models
#'
#' @description
#' \code{"param()"} is a generic S3 function that displays (and invisibly returns)
#' model hyperparameters. Methods are provided for \code{"krr"} objects.
#'
#' @param x An object.
#' @param ... Additional arguments passed to methods.
#' @return A named list of hyperparameters (invisibly); may print side effects.
#' @export
param = function(x, ...) UseMethod("param")

#' @rdname param
#' @export
param.default = function(x, ...) {
  stop(sprintf("No 'param' method for objects of class: %s", paste(class(x), collapse=", ")))
}

#' Param method for fitted Kernel Ridge Regression models
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
#' Prints a human-readable panel and returns (invisibly) a named list
#' of hyperparameters.
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
#' model = fastkrr(X, y, kernel="gaussian", opt="nystrom",
#'               rho=1, lambda=1e-4, m=200, n_threads=4, fastcv=FALSE)
#'
#'
#' class(model)
#' param(model)
#'
#' @export
param.krr = function(x, ...) {
  model = x
  at = attributes(model)

  keys = c("rho","lambda","m","eps", "n_threads")
  present = intersect(keys, names(at))

  out = at[present]
  class(out) = c("krr_params", "list")

  # Call
  cat("Call:\n")
  print(at$call)
  cat("\n")

  # core options -> table
  core_keys = c("rho","lambda", "eps", "m", "n_threads")
  core_keys = intersect(core_keys, present)
  if (length(core_keys) > 0) {
    df = as.data.frame(at[core_keys], optional = TRUE)
    row.names(df) <- ""
    cat("Parameters :\n")
    print(df)
    cat("\n")
  }

  invisible(out)
}
