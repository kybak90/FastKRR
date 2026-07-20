test_that("coef.krr extracts model coefficients correctly as a numeric vector", {
  # Set seed for reproducibility
  set.seed(123)
  n <- 30
  d <- 2
  X <- matrix(runif(n * d), nrow = n, ncol = d)
  colnames(X) <- c("x1", "x2")
  y <- sin(X[, 1]) + rnorm(n, sd = 0.1)
  df <- data.frame(X, y = y)

  # Fit a FastKRR model
  model <- fastkrr(data = df, response = "y", kernel = "gaussian",
                   opt = "exact", rho = 1.0, lambda = 1e-3)

  # Extract coefficients using coef()
  extracted_coefs <- coef(model)

  # Validate output type, dimension, and structure
  expect_type(extracted_coefs, "double")
  expect_true(is.vector(extracted_coefs))
  expect_equal(length(extracted_coefs), length(model$coefficients))
  expect_equal(extracted_coefs, as.vector(model$coefficients))
})

test_that("coef.krr works across different approximation options", {
  set.seed(123)
  n <- 30
  d <- 2
  X <- matrix(runif(n * d), nrow = n, ncol = d)
  colnames(X) <- c("x1", "x2")
  y <- sin(X[, 1]) + rnorm(n, sd = 0.1)
  df <- data.frame(X, y = y)

  # Test for RFF option where coefficients represent beta weights
  model_rff <- fastkrr(data = df, response = "y", kernel = "gaussian",
                       opt = "rff", m = 10, rho = 1.0, lambda = 1e-3)

  extracted_coefs_rff <- coef(model_rff)

  expect_type(extracted_coefs_rff, "double")
  expect_true(is.vector(extracted_coefs_rff))
  expect_equal(length(extracted_coefs_rff), 10)
})
