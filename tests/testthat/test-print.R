test_that("print.krr displays model details and invisibly returns the original object", {
  set.seed(123)
  n <- 30
  d <- 2
  X <- matrix(runif(n * d), nrow = n, ncol = d)
  colnames(X) <- c("x1", "x2")
  y <- sin(X[, 1]) + rnorm(n, sd = 0.1)
  df <- data.frame(X, y = y)

  # Fit a FastKRR model
  fit <- fastkrr(
    data = df,
    response = "y",
    kernel = "gaussian",
    opt = "exact",
    rho = 1.0,
    lambda = 1e-3
  )

  # 1. Check console output text patterns
  expect_output(print(fit), "Call:")
  expect_output(print(fit), "Parameters:")

  # 2. Verify that print.krr invisibly returns the input krr object
  res <- print(fit)
  expect_s3_class(res, "krr")
  expect_equal(res, fit)
})

test_that("print.approx_kernel displays approximation details to console", {
  set.seed(123)
  n <- 30
  d <- 2
  X <- matrix(runif(n * d), nrow = n, ncol = d)

  # Construct an approximated kernel object
  K_approx <- approx_kernel(X = X, opt = "nystrom", m = 10, rho = 1.0, n_threads = 1)

  # Verify console output text patterns
  expect_output(print(K_approx), "Call:")
  expect_output(print(K_approx), "Options:")

  # Capture output dataframe returned by print.approx_kernel
  res_df <- print(K_approx)
  expect_s3_class(res_df, "data.frame")
  expect_true("opt" %in% names(res_df))
})
