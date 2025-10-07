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
print(model_exact)
print(model_pivoted)
print(model_rff)
print(model_nystrom)

set.seed(1)
n = 1000
rho = 1
X = runif(n, 0, 1)
y = sin(2*pi*X^3) + rnorm(n, 0, 0.1)


model_exact = fastkrr(X, y, kernel = "gaussian", rho = rho, opt = "exact", fastcv = TRUE, verbose = FALSE)
model_pivoted = fastkrr(X, y, kernel = "gaussian", rho = rho, opt = "pivoted", fastcv = TRUE, verbose = FALSE)
model_rff = fastkrr(X, y, kernel = "gaussian", rho = rho, opt = "rff", fastcv = TRUE, verbose = FALSE)
model_nystrom = fastkrr(X, y, kernel = "gaussian", rho = rho, opt = "nystrom", fastcv = TRUE, verbose = FALSE)

print.kernel_matrix = function(ker_mat, ...){
  if(!is.list(ker_mat)){
    cat("Showing top-left 6x6:\n")
    print(ker_mat[1:6, 1:6])
  }else{
    cat("Showing top-left 6x6:\n")
    print(attributes(ker_mat)$K_approx[1:6, 1:6])
  }
}
print(attr(model_exact, "K"))
print(attr(model_pivoted, "K"))
print(attr(model_rff, "K"))
print(attr(model_nystrom, "K"))



ker_mat = K_nystrom
idx = names(attributes(ker_mat)) %in% c("opt", "m", "eps", "n_threads")
df = as.data.frame(attributes(ker_mat)[idx])
row.names(df) = ""
print(df)

names(attributes(K_pivoted))
names(attributes(K_nystrom))
names(attributes(K_rff))
c("opt", "m", "eps", "n_threads")





model = model_exact
idx = names(attributes(model)) %in% c("kernel", "opt", "rho", "fastcv", "n_threads", "m")
df = as.data.frame(attributes(model)[idx])
row.names(df) = ""

cat("Options:\n")
print(df)
