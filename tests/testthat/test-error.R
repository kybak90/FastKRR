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

test_that("error.krr supports x_new/y_new and stored response names", {
  set.seed(123)
  n <- 30
  X <- matrix(runif(n * 2), ncol = 2)
  y <- sin(X[, 1]) + rnorm(n, sd = 0.1)
  df <- data.frame(X, y = y)
  X_new <- matrix(runif(20 * 2), ncol = 2)
  y_new <- sin(X_new[, 1]) + rnorm(20, sd = 0.1)
  df_new <- data.frame(X_new, y = y_new)

  m_df <- fastkrr(df, "y", opt = "exact", lambda = 1e-3, verbose = FALSE)
  m_xy <- fastkrr(X, y, opt = "exact", lambda = 1e-3, verbose = FALSE)
  expected <- mean((y_new - predict(m_df, X_new))^2)

  # data.frame model: data_new and x_new/y_new both work
  expect_equal(error(m_df, data_new = df_new), expected)
  expect_equal(error(m_df, x_new = X_new, y_new = y_new), expected)

  # matrix model: x_new/y_new work, data_new gives a clear message
  expect_equal(error(m_xy, x_new = X_new, y_new = y_new), expected)
  expect_equal(error(m_xy), mean((y - m_xy$fitted.values)^2))
  expect_error(error(m_xy, data_new = df_new), "use x_new and y_new")

  # response passed through a variable is resolved to the column name
  resp <- "y"
  m_var <- fastkrr(df, response = resp, opt = "exact", lambda = 1e-3, verbose = FALSE)
  expect_equal(error(m_var, data_new = df_new), expected)

  # invalid combinations
  expect_error(error(m_df, data_new = df_new, x_new = X_new, y_new = y_new), "not both")
  expect_error(error(m_df, x_new = X_new), "Both x_new and y_new")
  expect_error(error(m_df, x_new = X_new, y_new = y_new[-1]), "must match")
})
