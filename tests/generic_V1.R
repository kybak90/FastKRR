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
    print(attributes(model)$K[1:6, 1:6])
    cat("\n")

    cat("Options:\n")
    cat(" ", "kernel=", attributes(model)$kernel, ", ",
        "opt=", attributes(model)$opt, ", ",
        "rho=", attributes(model)$rho, ", ",
        "fastcv=", attributes(model)$fastcv, ", ",
        "n_threads=", attributes(model)$n_threads,
        sep="")
  }else{
    cat("Approximated kernel matrix (showing top-left 6x6):\n")
    print(attributes(model)$K_approx[1:6, 1:6])
    cat("\n")
    cat("Options:\n")
    cat(" ", "kernel=", attributes(model)$kernel, ", ",
        "opt=", attributes(model)$opt, ", ",
        "rho=", attributes(model)$rho, ", ",
        "fastcv=", attributes(model)$fastcv, ", ",
        "n_threads=", attributes(model)$n_threads, ", ",
        "m=", attributes(model)$m,
        sep="")
  }
}

print.kernel_matrix = function(ker_mat, ...){
  if(!is.list(ker_mat)){
    cat("Showing top-left 6x6:\n")
    print(ker_mat[1:6, 1:6])
  }else{
    cat("Showing top-left 6x6:\n")
    print(attributes(ker_mat)$K_approx[1:6, 1:6])
  }
}
