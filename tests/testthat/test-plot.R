test_that("plot.krr returns a valid ggplot object for 1D input", {
  # Generate 1D synthetic dataset
  set.seed(123)
  n <- 30
  X <- matrix(runif(n), nrow = n, ncol = 1)
  colnames(X) <- "x1"
  y <- sin(2 * pi * X[, 1]) + rnorm(n, sd = 0.1)
  df <- data.frame(X, y = y)

  # Fit model
  fit <- fastkrr(data = df, response = "y", kernel = "gaussian",
                 opt = "exact", lambda = 1e-3)

  # Test plot with points (default)
  p1 <- plot(fit, show_points = TRUE)
  expect_s3_class(p1, "ggplot")
  expect_equal(length(p1$layers), 2) # geom_point + geom_line

  # Test plot without points
  p2 <- plot(fit, show_points = FALSE)
  expect_s3_class(p2, "ggplot")
  expect_equal(length(p2$layers), 1) # geom_line only
})

test_that("plot.krr throws informative errors on invalid inputs or d >= 2", {
  set.seed(123)

  # Case 1: Input is not a 'krr' object
  expect_error(
    plot.krr("not_a_krr_object"),
    "Input must be a fitted model of class 'krr'."
  )

  # Case 2: Multi-dimensional input matrix (d >= 2)
  n <- 30
  d <- 2
  X_2d <- matrix(runif(n * d), nrow = n, ncol = d)
  colnames(X_2d) <- c("x1", "x2")
  y <- X_2d[, 1] + X_2d[, 2] + rnorm(n, sd = 0.1)
  df_2d <- data.frame(X_2d, y = y)

  fit_2d <- fastkrr(data = df_2d, response = "y", kernel = "gaussian",
                    opt = "exact", lambda = 1e-3)

  # Expect error for multivariate inputs
  expect_error(
    plot(fit_2d),
    "plot.krr\\(\\) currently supports only 1D input."
  )
})
