test_that("fastkrr fits correctly across all computation options (exact, nystrom, pivoted, rff)", {
  set.seed(123)
  n <- 30
  d <- 2
  X <- matrix(runif(n * d), nrow = n, ncol = d)
  colnames(X) <- c("x1", "x2")
  y <- sin(X[, 1]) + rnorm(n, sd = 0.1)
  df <- data.frame(X, y = y)

  # 1. Exact option & S3 object structure verification
  fit_exact <- fastkrr(data = df, response = "y", kernel = "gaussian", opt = "exact", lambda = 1e-4)
  expect_s3_class(fit_exact, "krr")
  expect_type(fit_exact$coefficients, "double")
  expect_equal(length(fit_exact$coefficients), n)
  expect_equal(length(fit_exact$fitted.values), n)
  expect_false(is.null(fit_exact$chol_factor))

  # 2. Nystrom option
  fit_nys <- fastkrr(data = df, response = "y", kernel = "gaussian", opt = "nystrom", lambda = 1e-4, m = 10)
  expect_equal(length(fit_nys$fitted.values), n)
  expect_false(is.null(fit_nys$approx_factor))

  # 3. Pivoted Cholesky option
  fit_piv <- fastkrr(data = df, response = "y", kernel = "gaussian", opt = "pivoted", lambda = 1e-4, m = 10)
  expect_equal(length(fit_piv$fitted.values), n)
  expect_false(is.null(fit_piv$approx_factor))

  # 4. RFF option
  fit_rff <- fastkrr(data = df, response = "y", kernel = "gaussian", opt = "rff", lambda = 1e-4, m = 10)
  expect_equal(length(fit_rff$fitted.values), n)
  expect_equal(length(fit_rff$coefficients), 10) # m-dimensional weight vector
  expect_false(is.null(fit_rff$approx_factor))
})

test_that("predict.krr handles both missing newdata and out-of-sample predictions", {
  set.seed(123)
  n <- 30
  d <- 2
  X <- matrix(runif(n * d), nrow = n, ncol = d)
  colnames(X) <- c("x1", "x2")
  y <- sin(X[, 1]) + rnorm(n, sd = 0.1)
  df <- data.frame(X, y = y)

  fit <- fastkrr(data = df, response = "y", kernel = "gaussian", opt = "exact", lambda = 1e-4)

  # Returns fitted values when newdata is missing
  expect_equal(predict(fit), fit$fitted.values)

  # Evaluates predictions on new input design matrix
  new_X <- matrix(runif(10 * d), nrow = 10, ncol = d)
  colnames(new_X) <- c("x1", "x2")

  preds <- predict(fit, new_X)
  expect_equal(length(preds), 10)
  expect_false(anyNA(preds))

  # Evaluates predictions for RFF option
  fit_rff <- fastkrr(data = df, response = "y", kernel = "gaussian", opt = "rff", lambda = 1e-4, m = 8)
  preds_rff <- predict(fit_rff, new_X)
  expect_equal(length(preds_rff), 10)
  expect_false(anyNA(preds_rff))
})

test_that("lambda specification and selection_method behave correctly", {
  set.seed(123)
  n <- 25
  d <- 2
  X <- matrix(runif(n * d), nrow = n, ncol = d)
  colnames(X) <- c("x1", "x2")
  y <- sin(X[, 1]) + rnorm(n, sd = 0.1)
  df <- data.frame(X, y = y)

  # 1. REML tuning with length-2 min/max vector
  fit_reml <- fastkrr(data = df, response = "y", kernel = "gaussian", opt = "exact",
                      selection_method = "REML", lambda = c(1e-8, 1e-2))
  expect_type(fit_reml$lambda, "double")
  expect_length(fit_reml$lambda, 1)

  # 2. Error when passing length-2 vector to CV methods
  expect_error(
    fastkrr(data = df, response = "y", selection_method = "exactCV", lambda = c(1e-4, 1e-2)),
    "only supported when selection_method = 'REML'"
  )

  # 3. Error when passing length >= 3 grid to REML
  expect_error(
    fastkrr(data = df, response = "y", selection_method = "REML", lambda = c(1e-4, 1e-3, 1e-2)),
    "REML does not support a lambda grid of length >= 3"
  )
})

test_that("na.rm argument handles missing values appropriately", {
  set.seed(123)
  n <- 30
  d <- 2
  X <- matrix(runif(n * d), nrow = n, ncol = d)
  colnames(X) <- c("x1", "x2")
  y <- sin(X[, 1]) + rnorm(n, sd = 0.1)

  # Inject missing value
  X[5, 1] <- NA
  df <- data.frame(X, y = y)

  # Throws error when na.rm = FALSE
  expect_error(
    fastkrr(data = df, response = "y", kernel = "gaussian", opt = "exact", lambda = 1e-4, na.rm = FALSE),
    "Missing values"
  )

  # Removes missing rows gracefully when na.rm = TRUE
  expect_no_error({
    fit_na <- fastkrr(data = df, response = "y", kernel = "gaussian", opt = "exact", lambda = 1e-4, na.rm = TRUE)
  })
  expect_equal(fit_na$removed_row_idx, 5)
  expect_equal(length(fit_na$fitted.values), n - 1)
})
