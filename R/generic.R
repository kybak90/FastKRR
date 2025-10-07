#' @export
print.krr = function(model, ...){
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
    cat("Kernel matrix (showing top-left 6x6):\n")
    print(attributes(model)$K) ## print.kernel_matrix
    cat("\n")

    idx = names(attributes(model)) %in% c("kernel", "opt", "rho", "fastcv", "n_threads", "m")
    df = as.data.frame(attributes(model)[idx])
    row.names(df) = ""

    cat("Options:\n")
    print(df)

  }else{
    cat("Approximated kernel matrix (showing top-left 6x6):\n")
    print(attributes(model)$K_approx) ## print.kernel_matrix
    cat("\n")

    idx = names(attributes(model)) %in% c("kernel", "opt", "rho", "fastcv", "n_threads", "m")
    df = as.data.frame(attributes(model)[idx])
    row.names(df) = ""

    cat("Options:\n")
    print(df)
  }
}

#' @export
print.approx_kernel = function(approx_mat, ...){
  cat("Call:\n")
  print(attributes(approx_mat)$call)
  cat("\n")
  print(attributes(approx_mat)$K_approx) ## print.kernel_matrix
  cat("\n")


  idx = names(attributes(approx_mat)) %in% c("opt", "m", "eps", "n_threads")
  df = as.data.frame(attributes(approx_mat)[idx])
  row.names(df) = ""

  cat("Options:\n")
  print(df)
}

#' @export
print.kernel_matrix = function(ker_mat, ...){
  cat("Showing top-left 6x6:\n")
  print(ker_mat[1:6, 1:6])
}
