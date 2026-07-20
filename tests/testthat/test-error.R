test_that("error.krr computes training MSE correctly when data_new is NULL", {
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

  # Compute training MSE
  train_mse <- error(model)

  # Manual expectation check
  expected_mse <- mean((df$y - model$fitted.values)^2)

  # Validate output type, structure, and numeric correctness
  expect_type(train_mse, "double")
  expect_length(train_mse, 1)
  expect_equal(train_mse, expected_mse)
})

test_that("error.krr computes prediction MSE (PMSE) correctly when data_new is provided", {
  set.seed(123)
  n <- 30
  d <- 2
  X <- matrix(runif(n * d), nrow = n, ncol = d)
  colnames(X) <- c("x1", "x2")
  y <- sin(X[, 1]) + rnorm(n, sd = 0.1)
  df <- data.frame(X, y = y)

  model <- fastkrr(data = df, response = "y", kernel = "gaussian",
                   opt = "exact", rho = 1.0, lambda = 1e-3)

  # Construct new test data frame
  new_n <- 20
  new_X <- matrix(runif(new_n * d), nrow = new_n, ncol = d)
  colnames(new_X) <- c("x1", "x2")
  new_y <- sin(new_X[, 1]) + rnorm(new_n, sd = 0.1)
  new_df <- data.frame(new_X, y = new_y)

  # Compute prediction MSE (PMSE)
  pmse <- error(model, data_new = new_df)

  # Manual prediction expectation check
  preds <- predict(model, newdata = as.matrix(new_X))
  expected_pmse <- mean((new_y - preds)^2)

  expect_type(pmse, "double")
  expect_length(pmse, 1)
  expect_equal(pmse, expected_pmse)
})

test_that("error.krr throws expected errors on invalid data_new inputs", {
  set.seed(123)
  n <- 20
  d <- 2
  X <- matrix(runif(n * d), nrow = n, ncol = d)
  colnames(X) <- c("x1", "x2")
  y <- sin(X[, 1]) + rnorm(n, sd = 0.1)
  df <- data.frame(X, y = y)

  model <- fastkrr(data = df, response = "y", kernel = "gaussian",
                   opt = "exact", rho = 1.0, lambda = 1e-3)

  # 1. Error when data_new is not a data.frame
  expect_error(
    error(model, data_new = matrix(1:10, nrow = 5)),
    "data_new must be a data.frame"
  )

  # 2. Error when data_new contains NA values
  na_df <- data.frame(x1 = c(1, NA), x2 = c(2, 3), y = c(1, 2))
  expect_error(
    error(model, data_new = na_df),
    "data_new contains missing values"
  )

  # 3. Error when the response column is missing in data_new
  wrong_col_df <- data.frame(x1 = runif(5), x2 = runif(5), target = runif(5))
  expect_error(
    error(model, data_new = wrong_col_df),
    "Response variable not found in data_new"
  )
})
