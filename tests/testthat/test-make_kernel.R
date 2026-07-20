test_that("make_kernel constructs symmetric kernel matrix when X_new is NULL", {
  set.seed(123)
  n <- 30
  d <- 2
  X <- matrix(runif(n * d), nrow = n, ncol = d)
  rho <- 1.0

  # 1. Gaussian kernel symmetric matrix
  K_gauss <- make_kernel(X = X, kernel = "gaussian", rho = rho)
  expect_true(is.matrix(K_gauss))
  expect_equal(dim(K_gauss), c(n, n))
  expect_equal(K_gauss, t(K_gauss)) # Symmetry check
  expect_equal(diag(K_gauss), rep(1.0, n)) # K(x_i, x_i) = exp(0) = 1

  # 2. Laplace kernel symmetric matrix
  K_laplace <- make_kernel(X = X, kernel = "laplace", rho = rho, n_threads = 1)
  expect_true(is.matrix(K_laplace))
  expect_equal(dim(K_laplace), c(n, n))
  expect_equal(K_laplace, t(K_laplace))
  expect_equal(diag(K_laplace), rep(1.0, n))
})

test_that("make_kernel constructs rectangular cross-kernel matrix when X_new is provided", {
  set.seed(123)
  n <- 30
  n_new <- 15
  d <- 2
  X <- matrix(runif(n * d), nrow = n, ncol = d)
  X_new <- matrix(runif(n_new * d), nrow = n_new, ncol = d)
  rho <- 1.0

  # 1. Gaussian cross-kernel matrix (n_new x n)
  K_cross_gauss <- make_kernel(X = X, X_new = X_new, kernel = "gaussian", rho = rho)
  expect_true(is.matrix(K_cross_gauss))
  expect_equal(dim(K_cross_gauss), c(n_new, n))

  # 2. Laplace cross-kernel matrix (n_new x n)
  K_cross_laplace <- make_kernel(X = X, X_new = X_new, kernel = "laplace", rho = rho, n_threads = 1)
  expect_true(is.matrix(K_cross_laplace))
  expect_equal(dim(K_cross_laplace), c(n_new, n))
})
