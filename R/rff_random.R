# Generate random Fourier features parameters
rff_random = function(m, d = d, rho = 1, kernel = "gaussian"){
  W = matrix(0, nrow = m, ncol = d)
  b = rep(0, m)

  if(kernel == 'gaussian'){
    W = matrix(rnorm(m * d, 0, sqrt(2 * rho)), nrow = m)
  }else if(kernel == "laplace"){
    W = matrix(rcauchy(m * d, location = 0, scale = 1/sqrt(2 * rho)), nrow = m)
  }

  b = runif(m, 0, 2 * pi)

  return(list("W" = W,
              "b" = b))
}


rff_kernel = function(X,
                      W = NULL,
                      b = NULL,
                      m = NULL,
                      rho = 1,
                      kernel = "gaussian")
{
  n = nrow(X)
  d = ncol(X)

  # m default: n/10 (min: 1)
  if (is.null(m)) m = n / 10 * log(d + 5)
  m = as.integer(max(1, floor(m)))

  if (is.null(W) || is.null(b)) {
    # Temporarily set global variable `d` because rff.random references it
    old_d_exists = exists("d", envir = .GlobalEnv, inherits = FALSE)
    if (old_d_exists) old_d = get("d", envir = .GlobalEnv)
    on.exit({
    }, add = TRUE)

    rb = rff_random(m = m, d = d, rho = rho, kernel = kernel)
    W = rb$W
    b = rb$b
  } else {

    if (ncol(W) != d) stop("ncol(W) must equal ncol(X).")
    if (length(b) != nrow(W)) stop("length(b) must equal nrow(W).")
    m = nrow(W)
  }

  Z = make_Z(X, W, b)

  K_approx = tcrossprod(Z)

  return(list(
    K_approx = K_approx,
    m = m,
    W = W,
    b = b
  ))
}
