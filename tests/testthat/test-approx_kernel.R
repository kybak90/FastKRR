test_that("approx_kernel computes Nystrom approximation correctly", {
  # Set up reproducible synthetic data
  set.seed(123)
  n <- 40
  d <- 2
  X <- matrix(runif(n * d), nrow = n, ncol = d)
  m <- 10
  rho <- 1.0

  # Execute Nystrom kernel approximation
  res <- approx_kernel(X = X, opt = "nystrom", kernel = "gaussian", m = m, rho = rho, n_threads = 1)

  # Check object structure, dimensions, and output components
  expect_s3_class(res, "approx_kernel")
  expect_equal(res$opt, "nystrom")
  expect_equal(dim(res$K_approx), c(n, n))
  expect_equal(dim(res$approx_factor), c(n, m))
  expect_equal(res$m, m)
  expect_equal(res$rho, rho)
  expect_false(is.null(res$n_threads))
})

test_that("approx_kernel computes Pivoted Cholesky approximation and handles dynamic eps", {
  set.seed(123)
  n <- 40
  d <- 2
  X <- matrix(runif(n * d), nrow = n, ncol = d)
  m <- 15
  rho <- 1.0

  # Test 1: Gaussian kernel defaults to eps = 1e-6
  res_gauss <- approx_kernel(X = X, opt = "pivoted", kernel = "gaussian", m = m, rho = rho)
  expect_equal(res_gauss$opt, "pivoted")
  expect_equal(dim(res_gauss$K_approx), c(n, n))
  expect_equal(nrow(res_gauss$approx_factor), n)
  expect_equal(res_gauss$eps, 1e-6)

  # Test 2: Laplace kernel defaults to eps = 1e-4
  res_laplace <- approx_kernel(X = X, opt = "pivoted", kernel = "laplace", m = m, rho = rho)
  expect_equal(res_laplace$eps, 1e-4)
})

test_that("approx_kernel handles Random Fourier Features (RFF) with auto & user-supplied W, b", {
  set.seed(123)
  n <- 30
  d <- 3
  X <- matrix(runif(n * d), nrow = n, ncol = d)
  m <- 12
  rho <- 1.0

  # Case 1: Automatically generated W and b
  res_auto <- approx_kernel(X = X, opt = "rff", kernel = "gaussian", m = m, rho = rho, n_threads = 1)
  expect_equal(res_auto$opt, "rff")
  expect_equal(dim(res_auto$K_approx), c(n, n))
  expect_equal(dim(res_auto$approx_factor), c(n, m))
  expect_equal(dim(res_auto$W), c(m, d))
  expect_equal(length(res_auto$b), m)
  expect_false(res_auto$used_supplied_Wb)

  # Case 2: User-supplied W and b
  custom_W <- matrix(rnorm(m * d), nrow = m, ncol = d)
  custom_b <- runif(m, 0, 2 * pi)

  res_supplied <- approx_kernel(X = X, opt = "rff", kernel = "gaussian",
                                rho = rho, W = custom_W, b = custom_b, n_threads = 1)
  expect_true(res_supplied$used_supplied_Wb)
  expect_equal(res_supplied$m, m)
  expect_equal(res_supplied$W, custom_W)
  expect_equal(res_supplied$b, custom_b)
})

test_that("approx_kernel throws expected errors on invalid inputs", {
  set.seed(123)
  X <- matrix(runif(20), nrow = 10, ncol = 2)

  # Error when X is missing
  expect_error(approx_kernel(X = NULL, rho = 1), "Argument 'X' must be provided")

  # Error when rho is missing
  expect_error(approx_kernel(X = X, opt = "nystrom"), "'rho' must be provided")

  # Error when only one of W or b is provided in RFF mode
  custom_W <- matrix(rnorm(10), nrow = 5, ncol = 2)
  expect_error(
    approx_kernel(X = X, opt = "rff", rho = 1, W = custom_W, b = NULL),
    "provide both 'W' and 'b' or neither"
  )
})
