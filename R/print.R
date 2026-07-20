#' Print method for fitted Kernel Ridge Regression models
#'
#' @description
#' Displays a concise summary of key information from a fitted Kernel Ridge Regression (KRR) model,
#' including the original function call and the main hyperparameter options used during fitting.
#'
#' @param x An S3 object of class \code{krr}, typically returned by
#'   \code{\link{fastkrr}}.
#' @param ... Additional arguments (currently ignored).
#'
#' @return Invisibly returns the input \code{krr} object after printing the summary to the console.
#'
#' @seealso \code{\link{fastkrr}}, \code{\link{print.approx_kernel}}
#'
#' @examples
#' # Data setting
#' set.seed(1)
#' lambda = 1e-4
#' d = 1
#' n = 50
#' rho = 1
#' X = matrix(runif(n*d, 0, 1), nrow = n, ncol = d)
#' y = sin(2 * pi * rowMeans(X)^3) + rnorm(n, mean = 0, sd = 0.1)
#'
#' data = data.frame(X, y = y)
#'
#' # Example: exact
#' model = fastkrr(data = data, response = "y",
#'                  kernel = "gaussian", opt = "exact",
#'                  rho = rho, lambda = lambda)
#'
#' print(model)
#' @export
print.krr = function(x, ...) {
  model = x

  cat("Call:\n")
  print(model$call)
  cat("\n")

  keys = c("kernel", "opt", "selection_method",
           "rho", "lambda", "m", "eps", "n_threads")
  present = intersect(keys, names(model))

  if (length(present) > 0) {
    df = as.data.frame(model[present], optional = TRUE)

    cat("Parameters:\n")
    print(df, row.names = FALSE)
    cat("\n")
  }

  invisible(model)
}


#' Print method for approximated kernel matrices
#'
#' @description
#' Displays the key algorithmic choices and hyperparameter options used to construct
#' an approximated kernel matrix object.
#'
#' @details
#' The function summarizes critical metadata attributes stored within the
#' \code{approx_kernel} object, including the approximation strategy (\code{opt}),
#' kernel type, approximation degree (\code{m}), numerical tolerance (\code{eps}),
#' and the number of parallel computational threads used.
#'
#' @param x An S3 object created by \code{\link{approx_kernel}}.
#' @param ... Additional arguments (currently ignored).
#'
#' @return Invisibly returns the input \code{approx_kernel} object after printing the options.
#'
#' @seealso \code{\link{approx_kernel}}, \code{\link{print.krr}}
#'
#' @examples
#' # Data setting
#' set.seed(1)
#' d = 1
#' n = 100
#' rho = 1
#' X = matrix(runif(n*d, 0, 1), nrow = n, ncol = d)
#'
#' # Example: nystrom
#' K_nystrom = approx_kernel(X = X, opt = "nystrom", rho = rho, n_threads = 1)
#'
#' print(K_nystrom)
#' @export
print.approx_kernel = function(x, ...){
  approx_object = x
  cat("Call:\n")
  print(approx_object$call)

  cat("\n")
  idx = names(approx_object) %in% c("opt", "kernel", "d", "rho", "m", "eps", "n_threads")
  df = as.data.frame(approx_object[idx])
  row.names(df) = ""
  cat("Options:\n")
  print(df)
}
