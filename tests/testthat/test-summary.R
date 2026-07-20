test_that("summary.krr prints model summary and returns training error", {
  # Set up reproducible synthetic data
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

  # Verify console output sections from param() and error() calls
  expect_output(summary(fit), "Call:")
  expect_output(summary(fit), "Parameters:")
  expect_output(summary(fit), "Training error:")

  # Capture the return value of summary() which returns training error
  res_err <- summary(fit)

  # Validate that the returned value is a positive numeric training error
  expect_type(res_err, "double")
  expect_true(is.numeric(res_err))
  expect_gt(res_err, 0)
  expect_equal(res_err, as.numeric(error(fit)))
})
