#' Print method for fitted Kernel Ridge Regression models
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
#' @seealso \code{\link{fastkrr}}, \code{\link{print.approx_kernel}}
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
#'                 rho = rho, lambda = 1e-4)
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
#' Displays the approximated kernel matrix and key options used to construct it.
#'
#' @details
#' The function prints the stored approximated kernel matrix (top-left 6x6)
#' and summarizes options such as the approximation method (\code{opt}), approximaion degree (\code{m}),
#' numerical tolerance (\code{eps}), and number of threads used (\code{n_threads}).
#'
#' @param x An S3 object created by \code{\link{approx_kernel}}.
#' @param ... Additional arguments (currently ignored).
#'
#' @return An approximated kernel matrix and its associated options.
#'
#' @seealso \code{\link{approx_kernel}, \link{print.krr}}
#'
#' @examples
#' # Data setting
#' set.seed(1)
#' d = 2
#' n = 100
#' rho = 1
#' X = matrix(runif(n*d, 0, 1), nrow = n, ncol = d)
#'
#' # Example: nystrom
#' K_nystrom = approx_kernel(X = X, opt = "nystrom", d = d, rho = rho, n_threads = 1)
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
