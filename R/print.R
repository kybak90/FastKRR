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
#' class(model)
#'
#' print(model)
#' @export
print.krr = function(x, ...){
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

  cat("\n\n")
  if(attributes(model)$opt == "exact"){
    cat("Kernel matrix\n")
    print(attributes(model)$K) ## print.kernel_matrix()
    cat("\n")

    idx = names(attributes(model)) %in% c("kernel", "opt", "rho", "fastcv", "n_threads", "m")
    df = as.data.frame(attributes(model)[idx])
    row.names(df) = ""

    cat("Options:\n")
    print(df)

  }else{
    cat("Approximated kernel matrix\n")
    print(attributes(model)$K_approx) ## print.kernel_matrix()
    cat("\n")

    idx = names(attributes(model)) %in% c("kernel", "opt", "rho", "fastcv", "n_threads", "m")
    df = as.data.frame(attributes(model)[idx])
    row.names(df) = ""

    cat("Options:\n")
    print(df)
  }
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
#' @seealso \code{\link{approx_kernel}, \link{print.krr}}, \code{\link{print.kernel_matrix}}
#'
#' @examples
#' # Data setting
#' set.seed(1)
#' d = 1
#' n = 1000
#' m = 50
#' rho = 1
#' X = matrix(runif(n*d, 0, 1), nrow = n, ncol = d)
#' y = as.vector(sin(2*pi*rowMeans(X)^3) + rnorm(n, 0, 0.1))
#' K = make_kernel(X, kernel = "gaussian", rho = rho)
#'
#' # Example: nystrom
#' K_nystrom = approx_kernel(K = K, opt = "nystrom", m = m, d = d, rho = rho, n_threads = 1)
#' class(K_nystrom)
#'
#' print(K_nystrom)
#' @export
print.approx_kernel = function(x, ...){
  approx_object = x
  cat("Call:\n")
  print(attributes(approx_object)$call)
  cat("\n")
  print(attributes(approx_object)$K_approx) ## print.kernel_matrix()
  cat("\n")


  idx = names(attributes(approx_object)) %in% c("opt", "m", "eps", "n_threads")
  df = as.data.frame(attributes(approx_object)[idx])
  row.names(df) = ""

  cat("Options:\n")
  print(df)
}


#' Print method for kernel matrices
#'
#' @description
#' Displays the top-left 6×6 portion of a kernel or approximated kernel matrix
#' for quick inspection.
#'
#' @param x An object of class \code{kernel_matrix}, which may represent
#'   either an exact kernel matrix (from \code{\link{make_kernel}} or
#'   \code{\link{fastkrr}}) or an approximated kernel matrix (from
#'   \code{\link{approx_kernel}}).
#' @param ... Additional arguments (currently ignored).
#'
#' @return A top-left 6x6 block of the kernel matrix to the console.
#'
#' @seealso \code{\link{approx_kernel}}, \code{\link{fastkrr}}, \code{\link{print.approx_kernel}}, \code{\link{print.krr}}
#'
#' @examples
#' # data setting
#' set.seed(1)
#' n = 1000 ; d = 1
#' m = 100
#' rho = 1
#' X = matrix(runif(n*d, 0, 1), nrow = n, ncol = d)
#' y = as.vector(sin(2*pi*X^3) + rnorm(n, 0, 0.1))
#'
#'
#' # Example for fastkrr
#' fit_pivoted = fastkrr(X, y,
#'                       kernel = "gaussian", opt = "pivoted",
#'                       m = 100, fastcv = TRUE, verbose = FALSE)
#'
#' class(attr(fit_pivoted, "K"))
#' print(class(attr(fit_pivoted, "K")))
#'
#' class(attr(fit_pivoted, "K_approx"))
#' print(class(attr(fit_pivoted, "K_approx")))
#'
#'
#' # Example for make_kernel
#' K = make_kernel(X, kernel = "gaussian", rho = rho)
#'
#' class(K)
#' print(K)
#'
#'
#' # Example for make_kernel
#' K_rff = approx_kernel(X = X, opt = "rff", kernel = "gaussian",
#'                       d = d, rho = rho, n_threads = 1, m = 100)
#'
#' class(attr(K_rff, "K_approx"))
#' print(attr(K_rff, "K_approx"))
#' @export
print.kernel_matrix = function(x, ...){
  ker_mat = x
  cat("Showing top-left 6x6:\n")
  print(ker_mat[1:6, 1:6])
}
