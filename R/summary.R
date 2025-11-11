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
#' X = matrix(runif(n*d, 0, 1), nrow = n, ncol = d)
#' y = as.vector(sin(2*pi*rowMeans(X)^3) + rnorm(n, 0, 0.1))
#'
#' # Example: exact
#' model = fastkrr(X, y,
#'                 kernel = "gaussian", opt = "exact",
#'                 rho = rho, lambda = 1e-4)
#' class(model)
#'
#' summary(model)
#'
#' @export
summary.krr = function(object, ...){
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
    print(attributes(model)$K_approx) ## summary.kernel_matrix()
    cat("\n")

    idx = names(attributes(model)) %in% c("kernel", "opt", "rho", "fastcv", "n_threads", "m")
    df = as.data.frame(attributes(model)[idx])
    row.names(df) = ""

    cat("Options:\n")
    print(df)
  }
}
