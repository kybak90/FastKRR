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

test_that("fastkrr accepts data.frame, (X, y), and (x = X, y = y) with identical results", {
  set.seed(1)
  n <- 60
  d <- 2
  X <- matrix(runif(n * d), nrow = n, ncol = d)
  y <- sin(2 * pi * rowMeans(X)) + rnorm(n, sd = 0.1)
  df <- data.frame(X, y = y)
  X_new <- matrix(runif(10 * d), ncol = d)

  for (opt in c("exact", "pivoted", "nystrom", "rff")) {
    for (kernel in c("gaussian", "laplace")) {
      set.seed(2); f_df <- fastkrr(df, "y", kernel = kernel, opt = opt, lambda = 1e-3, verbose = FALSE)
      set.seed(2); f_xy <- fastkrr(X, y, kernel = kernel, opt = opt, lambda = 1e-3, verbose = FALSE)
      set.seed(2); f_nm <- fastkrr(x = X, y = y, kernel = kernel, opt = opt, lambda = 1e-3, verbose = FALSE)

      expect_equal(f_xy$coefficients, f_df$coefficients)
      expect_equal(f_nm$coefficients, f_df$coefficients)
      expect_equal(predict(f_xy, X_new), predict(f_df, X_new))
      expect_equal(unname(f_xy$x), unname(f_df$x))
    }
  }

  # lambda selection gives the same result
  f_df <- fastkrr(df, "y", lambda = 10^seq(-6, -2, length.out = 5), verbose = FALSE)
  f_xy <- fastkrr(X, y, lambda = 10^seq(-6, -2, length.out = 5), verbose = FALSE)
  expect_equal(f_xy$lambda, f_df$lambda)

  # interface and response are stored
  expect_equal(f_df$interface, "data.frame")
  expect_equal(f_df$response, "y")
  expect_equal(f_xy$interface, "matrix")
  expect_null(f_xy$response)
  expect_equal(colnames(f_xy$x), c("X1", "X2"))
})

test_that("fastkrr matrix interface accepts vectors and data frames and validates input", {
  set.seed(1)
  n <- 40
  x1 <- runif(n)
  y <- sin(2 * pi * x1) + rnorm(n, sd = 0.1)

  # numeric vector X (d = 1) is treated as a one-column matrix
  f_vec <- fastkrr(x1, y, lambda = 1e-3, verbose = FALSE)
  f_mat <- fastkrr(matrix(x1, ncol = 1), y, lambda = 1e-3, verbose = FALSE)
  expect_equal(f_vec$coefficients, f_mat$coefficients)
  expect_equal(ncol(f_vec$x), 1)

  # numeric data frame X
  f_dfx <- fastkrr(data.frame(a = x1), y, lambda = 1e-3, verbose = FALSE)
  expect_equal(f_dfx$coefficients, f_mat$coefficients)
  expect_equal(colnames(f_dfx$x), "a")

  # a predictor named "y" does not clash with the response
  f_y <- fastkrr(data.frame(y = x1), y, lambda = 1e-3, verbose = FALSE)
  expect_equal(f_y$coefficients, f_mat$coefficients)

  # input errors
  expect_error(fastkrr(data.frame(x1, y = y)), "response is missing")
  expect_error(fastkrr(x1, y[-1], lambda = 1e-3), "nrow(x) must match length(y)", fixed = TRUE)
  expect_error(fastkrr(x1, factor(y > 0), lambda = 1e-3), "single column name of data or a numeric response vector")
  expect_error(fastkrr(data.frame(a = x1, b = letters[seq_len(n)]), y, lambda = 1e-3), "only numeric columns")
  expect_error(fastkrr(x = x1, lambda = 1e-3), "Both x and y must be supplied")
  expect_error(fastkrr(data.frame(x1, y = y), "y", x = x1, y = y), "either data/response or x/y")
  expect_error(fastkrr(x1, "y", lambda = 1e-3), "data must be a data.frame")
})

test_that("na.rm works the same way for the matrix interface", {
  set.seed(1)
  n <- 30
  X <- matrix(runif(n * 2), ncol = 2)
  y <- sin(2 * pi * rowMeans(X)) + rnorm(n, sd = 0.1)
  X[3, 1] <- NA
  y[7] <- NA
  df <- data.frame(X, y = y)

  expect_error(fastkrr(X, y, lambda = 1e-3, verbose = FALSE), "Missing values")
  f_df <- suppressMessages(fastkrr(df, "y", lambda = 1e-3, verbose = FALSE, na.rm = TRUE))
  f_xy <- suppressMessages(fastkrr(X, y, lambda = 1e-3, verbose = FALSE, na.rm = TRUE))
  expect_equal(f_xy$removed_row_idx, c(3L, 7L))
  expect_equal(f_xy$coefficients, f_df$coefficients)
})
