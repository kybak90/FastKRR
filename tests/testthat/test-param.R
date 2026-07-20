test_that("param.default throws an informative error for unsupported classes", {
  # Expect error when passing unsupported object type (e.g., numeric vector or lm)
  expect_error(
    param("invalid_object"),
    "No 'param' method for objects of class: character"
  )

  dummy_list <- list(a = 1)
  class(dummy_list) <- "unsupported_class"
  expect_error(
    param(dummy_list),
    "No 'param' method for objects of class: unsupported_class"
  )
})

test_that("param.krr extracts hyperparameters and returns a krr_params object invisibly", {
  set.seed(123)
  n <- 30
  d <- 2
  X <- matrix(runif(n * d), nrow = n, ncol = d)
  colnames(X) <- c("x1", "x2")
  y <- sin(X[, 1]) + rnorm(n, sd = 0.1)
  df <- data.frame(X, y = y)

  # Fit model using fastkrr
  fit <- fastkrr(
    data = df,
    response = "y",
    kernel = "gaussian",
    opt = "pivoted",
    rho = 1.0,
    lambda = 1e-4,
    n_threads = 1
  )

  # Check invisible return value structure and class
  res <- param(fit)
  expect_s3_class(res, "krr_params")
  expect_type(res, "list")
  expect_equal(res$kernel, "gaussian")
  expect_equal(res$opt, "pivoted")
  expect_equal(res$rho, 1.0)
  expect_equal(res$lambda, 1e-4)
})

test_that("param.krr prints human-readable output to the console", {
  set.seed(123)
  n <- 20
  d <- 1
  X <- matrix(runif(n * d), nrow = n, ncol = d)
  y <- X[, 1] * 2 + rnorm(n, sd = 0.05)
  df <- data.frame(X, y = y)

  fit <- fastkrr(
    data = df,
    response = "y",
    kernel = "gaussian",
    opt = "exact",
    lambda = 1e-3
  )

  # Verify console print output contains Call and Parameters sections
  expect_output(param(fit), "Call:")
  expect_output(param(fit), "Parameters:")
})
