test_that("krr_reg initializes correctly and validates mode", {
  skip_if_not_installed("parsnip")

  # 1. Valid specification creation
  spec <- krr_reg(
    kernel = "gaussian",
    opt = "exact",
    penalty = 0.01,
    rho = 1.0
  )
  expect_s3_class(spec, "krr_reg")
  expect_s3_class(spec, "model_spec")
  expect_equal(spec$mode, "regression")

  # 2. Invalid mode error handling check
  expect_error(
    krr_reg(mode = "classification"),
    "`mode` should be 'regression'"
  )
})

test_that("tunable.krr_reg returns appropriate parameter information", {
  skip_if_not_installed("parsnip")
  skip_if_not_installed("generics")

  spec <- krr_reg()
  tb <- generics::tunable(spec)

  expect_s3_class(tb, "tbl_df")
  expect_equal(tb$name, "penalty")
  expect_equal(tb$component, "krr_reg")
})

test_that("update.krr_reg modifies parameters correctly with fresh = FALSE/TRUE", {
  skip_if_not_installed("parsnip")

  spec <- krr_reg(kernel = "gaussian", opt = "exact", rho = 1.0)

  # 1. Update individual parameters (fresh = FALSE)
  spec_updated <- update(
    spec,
    kernel = "laplace",
    opt = "nystrom",
    m = 20,
    eps = 1e-4,
    n_threads = 2,
    rho = 2.0,
    penalty = 0.05,
    na.rm = TRUE,
    fresh = FALSE
  )

  expect_s3_class(spec_updated, "krr_reg")
  expect_equal(rlang::quo_get_expr(spec_updated$args$kernel), "laplace")
  expect_equal(rlang::quo_get_expr(spec_updated$args$opt), "nystrom")

  # 2. Update with fresh = TRUE
  spec_fresh <- update(spec_updated, kernel = "gaussian", fresh = TRUE)
  expect_true("kernel" %in% names(spec_fresh$args))
  expect_false("opt" %in% names(spec_fresh$args)) # previous args cleared
})

test_that("make_krr_reg and helper functions handle registration idempotently", {
  skip_if_not_installed("parsnip")

  # Ensure model is registered
  make_krr_reg()

  # Calling again should do nothing and return NULL safely
  res <- make_krr_reg()
  expect_null(res)

  # Check helper function without breaking parsnip environment
  expect_invisible(.reset_parsnip_model("non_existent_model_xyz"))
})

test_that("fastkrr_fit_wrapper handles dataframe conversions and outcome name collisions", {
  set.seed(123)
  n <- 25
  X <- matrix(runif(n * 2), nrow = n, ncol = 2)
  y <- rnorm(n)

  # 1. Normal matrix fitting via wrapper
  fit_res <- fastkrr_fit_wrapper(
    x = X,
    y = y,
    kernel = "gaussian",
    opt = "exact",
    rho = 1.0,
    lambda = 0.01
  )
  expect_s3_class(fit_res, "krr")

  # 2. Name collision test (dataframe already containing `.outcome` column)
  df_x <- data.frame(x1 = runif(n), .outcome = runif(n))
  fit_col_test <- fastkrr_fit_wrapper(
    x = df_x,
    y = y,
    kernel = "gaussian",
    opt = "exact",
    rho = 1.0,
    lambda = 0.01
  )
  expect_s3_class(fit_col_test, "krr")
})

test_that("parsnip workflow integration (fit and predict) functions seamlessly", {
  skip_if_not_installed("parsnip")

  # Ensure model specification is cleanly registered before fitting
  make_krr_reg()

  set.seed(123)
  df_train <- data.frame(
    x1 = runif(30),
    x2 = runif(30),
    y  = rnorm(30)
  )

  df_test <- data.frame(
    x1 = runif(5),
    x2 = runif(5)
  )

  # Set parsnip model specification
  krr_spec <- krr_reg(
    kernel = "gaussian",
    opt = "exact",
    penalty = 0.01,
    rho = 1.0
  ) %>%
    parsnip::set_engine("fastkrr") %>%
    parsnip::set_mode("regression")

  # Fit parsnip model
  fit_parsnip <- parsnip::fit(krr_spec, y ~ x1 + x2, data = df_train)
  expect_s3_class(fit_parsnip, "model_fit")

  # Predict via parsnip
  preds <- predict(fit_parsnip, new_data = df_test)
  expect_s3_class(preds, "tbl_df")
  expect_equal(nrow(preds), 5)
  expect_true(".pred" %in% names(preds))
})
