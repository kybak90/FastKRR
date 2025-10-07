#' Kernel Ridge Regression
#'
#' @description
#' `krr_reg()` defines a Kernel Ridge Regression (KRR) model specification
#' for use with the tidymodels ecosystem via \pkg{parsnip}. This spec can be
#' paired with the \code{"fastkrr"} engine implemented in this package to fit
#' exact or kernel approximation (Nyström, Pivoted Cholesky, Random Fourier Features) within
#' \pkg{recipes}/\pkg{workflows} pipelines.
#'
#'
#' @param mode A single string; only `"regression"` is supported.
#' @param kernel Kernel matrix \eqn{K} has two kinds of Kernel ("gaussian", "laplace").
#' @param opt Method for constructing or approximating :
#'  \describe{
#'   \item{\code{"exact"}}{Construct the full kernel matrix
#'   \eqn{K \in \mathbb{R}^{n\times n}} using design matrix \eqn{X}.}
#'   \item{\code{"nystrom"}}{Construct a low-rank approximation of
#'       the kernel matrix \eqn{K \in \mathbb{R}^{n \times n}}
#'       using the Nyström approximation.}
#'   \item{\code{"pivoted"}}{Construct a low-rank approximation of
#'       the kernel matrix \eqn{K \in \mathbb{R}^{n \times n}}
#'       using Pivoted Cholesky decomposition.}
#' \item{\code{"rff"}}{Use Random Fourier Features to construct a feature map
#'   \eqn{Z \in \mathbb{R}^{n \times m}} (with \eqn{m} random features) so that
#'   \eqn{K \approx Z Z^\top}. Here, \eqn{m} is the number of features.}

#'  }
#' @param m  Approximation rank(number of random features) used for the low-rank kernel approximation.
#' @param eps Tolerance parameter used only in \code{"pivoted"}
#'   for stopping criterion of the Pivoted Cholesky decomposition.
#' @param n_threads Number of parallel threads. It is applied only for
#'  \code{opt = "nystrom"} or \code{opt = "rff"}, and for the
#'   Laplace kernel (\code{kernel = "laplace"}).
#' @param rho  Scaling parameter of the kernel(\eqn{\rho}).
#' @param penalty Regularization parameter.
#' @param fastcv If \code{TRUE}, accelerated cross-validation is
#'   performed via sequential testing (early stopping) as implemented in the \pkg{CVST} package.
#'
#' @return A parsnip model specification of class \code{"krr_reg"}.
#'
#' @examples
#' \donttest{
#' if (all(vapply(
#'   c("parsnip","stats","modeldata"),
#'   requireNamespace, quietly = TRUE, FUN.VALUE = logical(1)
#' ))) {
#' library(tidymodels)
#' library(parsnip)
#' library(stats)
#' library(modeldata)
#'
#' # Data analysis
#' data(ames)
#' ames = ames %>% mutate(Sale_Price = log10(Sale_Price))
#'
#' set.seed(502)
#' ames_split = initial_split(ames, prop = 0.80, strata = Sale_Price)
#' ames_train = training(ames_split) # dim (2342, 74)
#' ames_test  = testing(ames_split) # dim (588, 74)
#'
#' # Model spec
#' krr_spec = krr_reg(kernel = "gaussian", opt = "exact",
#'                    m = 50, eps = 1e-6, n_threads = 4,
#'                    rho = 1, penalty = tune()) %>%
#'  set_engine("fastkrr") %>%
#'  set_mode("regression")
#'
#' # Define rec
#' rec = recipe(Sale_Price ~ Longitude + Latitude, data = ames_train)
#'
#' # workflow
#' wf = workflow() %>%
#'   add_recipe(rec) %>%
#'   add_model(krr_spec)
#'
#' # Define hyper-parameter grid
#' param_grid = grid_regular(
#'   dials::penalty(range = c(-10, -3)),
#'   levels = 5
#' )
#'
#' # CV setting
#' set.seed(123)
#' cv_folds = vfold_cv(ames_train, v = 5, strata = Sale_Price)
#'
#' # Tuning
#' tune_results = tune_grid(
#'   wf,
#'   resamples = cv_folds,
#'   grid = param_grid,
#'   metrics = metric_set(rmse),
#'   control = control_grid(verbose = TRUE, save_pred = TRUE)
#' )
#'
#' # Result check
#' collect_metrics(tune_results)
#'
#' # Select best parameter
#' best_params = select_best(tune_results, metric = "rmse")
#'
#' # Finalized model spec using best parameter
#' final_spec = finalize_model(krr_spec, best_params)
#' final_wf = workflow() %>%
#'   add_recipe(rec) %>%
#'   add_model(final_spec)
#'
#' # Finalized fitting using best parameter
#' final_fit = final_wf %>% fit(data = ames_train)
#'
#' # Prediction
#' predict(final_fit, new_data = ames_test)
#' print(best_params)
#'
#' }}
#'
#' @export
krr_reg = function(mode = "regression", kernel = NULL, opt = NULL, eps = NULL,
                    n_threads = NULL, m = NULL, rho = NULL, penalty = NULL, fastcv = NULL) {
  if (mode != "regression") {
    rlang::abort("`mode` should be 'regression'.")
  }
  args = list(
    kernel  = rlang::enquo(kernel),
    opt     = rlang::enquo(opt),
    m       = rlang::enquo(m),
    eps      = rlang::enquo(eps),
    n_threads = rlang::enquo(n_threads),
    rho     = rlang::enquo(rho),
    penalty = rlang::enquo(penalty),
    fastcv  = rlang::enquo(fastcv)
  )
  parsnip::new_model_spec(
    "krr_reg",
    args     = args,
    eng_args = NULL,
    mode     = mode,
    method   = NULL,
    engine   = NULL
  )
}

#' Expose tunable parameters for `krr_reg`
#'
#' @description
#' Supplies a tibble of tunable arguments for `krr_reg()`.
#'
#' @param x A `krr_reg` model specification.
#' @param ... Not used; included for S3 method compatibility.
#'
#' @return A tibble (one row per tunable parameter)
#'    with columns 'name', 'call_info', 'source',
#'   'component', and 'component_id'.
#'
#' @importFrom generics tunable
#' @importFrom tibble tibble
#' @export
tunable.krr_reg = function(x, ...) {
  tibble::tibble(
    name = "penalty",
    call_info = list(list(pkg = "dials", fun = "penalty", range = c(-10, -3))),
    source = "model_spec",
    component = "krr_reg",
    component_id = "main"
  )
}


# ===================================================
# helper: remove an existing parsnip model registration
# ===================================================
.reset_parsnip_model = function(model) {
  if (!requireNamespace("parsnip", quietly = TRUE)) return(invisible())
  get_model_env = get("get_model_env", envir = asNamespace("parsnip"))
  env = get_model_env()
  if (exists(model, envir = env, inherits = FALSE)) {
    rm(list = model, envir = env)
  }
}

# =========================
# parsnip registration
# =========================
make_krr_reg = function(){
  # Idempotent guard: do nothing if already registered
  get_model_env = get("get_model_env", envir = asNamespace("parsnip"))
  env = get_model_env()
  if (exists("krr_reg", envir = env, inherits = FALSE)) {
    return(invisible(NULL))
  }

  # 1) model/mode
  parsnip::set_new_model("krr_reg")
  parsnip::set_model_mode(model = "krr_reg", mode = "regression")

  # 2) engine
  parsnip::set_model_engine(model = "krr_reg", mode = "regression", eng = "fastkrr")
  parsnip::set_dependency(model = "krr_reg", eng = "fastkrr", pkg = "FastKRR")

  # 3) parameter mapping
  parsnip::set_model_arg(
    model    = "krr_reg",
    eng      = "fastkrr",
    parsnip  = "penalty",
    original = "lambda",
    func     = list(pkg = "dials", fun = "penalty"),# func     = list(pkg = "base", fun = "identity"),
    has_submodel = FALSE
  )

  for (nm in c("kernel", "opt", "m", "eps", "n_threads", "rho", "fastcv")){
    parsnip::set_model_arg(
      model    = "krr_reg",
      eng      = "fastkrr",
      parsnip  = nm,
      original = nm,
      func     = list(pkg = "base", fun = "identity"),
      has_submodel = FALSE
    )
  }

  # 4) fit function mapping
  #    fastkrr(x, y, kernel, opt, m, rho, lambda, fastcv)
  parsnip::set_fit(
    model = "krr_reg",
    eng   = "fastkrr",
    mode  = "regression",
    value = list(
      interface = "matrix",
      protect   = c("x", "y"),
      func      = c(pkg = "FastKRR", fun = "fastkrr"),
      args      = list(
        x      = rlang::expr(x),
        y      = rlang::expr(y),
        kernel = rlang::expr(kernel),
        opt    = rlang::expr(opt),
        m   = rlang::expr(m),
        eps  = rlang::expr(eps),
        n_threads  = rlang::expr(n_threads),
        rho    = rlang::expr(rho),
        lambda = rlang::expr(penalty),
        fastcv = rlang::expr(fastcv)
      ),
      defaults = list()
    )
  )

  # 5) prediction function mapping
  #    predict(model, newdata) -> numeric vector
  parsnip::set_pred(
    model = "krr_reg",
    eng   = "fastkrr",
    mode  = "regression",
    type  = "numeric",
    value = list(
      pre  = NULL,
      post = NULL,
      func = c(fun = "predict"),
      args = list(
        object = rlang::expr(object$fit),  ### model   = rlang::expr(object$fit)
        newdata = rlang::expr(new_data)
      )
    )
  )

  # 6) encoding
  parsnip::set_encoding(
    model = "krr_reg",
    eng   = "fastkrr",
    mode  = "regression",
    options = list(
      predictor_indicators = "none",
      compute_intercept    = FALSE,
      remove_intercept     = FALSE,
      allow_sparse_x       = FALSE
    )
  )
}

# =========================
# package load hook — force reset and re-register on every load
# =========================
.onLoad = function(libname, pkgname) {
  if (!requireNamespace("parsnip", quietly = TRUE)) return(invisible())


  # .reset_parsnip_model("krr_reg")

  make_krr_reg()

}



# =========================
# optional: spec updater
# =========================
#' @importFrom stats update
#' @exportS3Method update krr_reg
update.krr_reg = function(object, parameters = NULL,
                           kernel = NULL, opt = NULL, m = NULL, eps = NULL,
                           n_threads = NULL, rho = NULL, fastcv = NULL, penalty = NULL,
                           fresh = FALSE, ...) {
  if (requireNamespace("parsnip", quietly = TRUE) &&
      "update.model_spec" %in% getNamespaceExports("parsnip")) {
    return(stats::update(
      object,
      kernel = kernel, opt = opt, m = m, eps = eps, n_threads = n_threads,
      rho = rho, fastcv = fastcv, penalty = penalty,
      parameters = parameters, fresh = fresh, ...
    ))
  }
  args_new = rlang::enquos(
    kernel=kernel, opt=opt, m=m, eps=eps, n_threads=n_threads,
    rho=rho, fastcv=fastcv, penalty=penalty,
    .ignore_empty="all"
  )
  object$args = if (fresh) args_new else utils::modifyList(object$args, args_new)
  object
}




utils::globalVariables(c(
  "x","y","opt","m", "eps","n_threads","rho","penalty","fastcv","object","new_data"
))
