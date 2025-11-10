#' Summary method for fitted Kernel Ridge Regression models
#'
#' @description
#' Displays key information from a fitted Kernel Ridge Regression (KRR) model,
#' including the original call, first few coefficients, a 6×6 block of the
#' kernel (or approximated kernel) matrix, and the main kernel options.
#'
#'
#'
#'
#'
#' @export
summary.krr = function(x, ...){
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
