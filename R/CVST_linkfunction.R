#' Predict responses for new data using fitted KRR model
#'
#' Generates predictions from a fitted Kernel Ridge Regression (KRR) model
#' for new data.
#'
#' @param object A S3 object of class \code{krr} created by \code{\link{fastkrr}}.
#' @param newdata New design matrix or data frame containing new observations
#'                for which predictions are to be made. If \code{newdata} is missing, the function returns fitted values.
#' @param ... Additional arguments (currently ignored).
#'
#'
#' @return A numeric vector of predicted values corresponding to \code{newdata} or fitted values.
#'
#' @seealso \code{\link{fastkrr}}, \code{\link{make_kernel}}
#'
#' @examples
#' # Data setting
#' n = 30
#' d = 2
#'
#' X = matrix(runif(n*d, 0, 1), nrow = n, ncol = d); colnames(X) = paste0("X", seq_len(d))
#' y = sin(2 * pi * rowMeans(X)^3) + rnorm(n, mean = 0, sd = 0.1)
#'
#' data = data.frame(X, y = y)
#'
#' lambda = 1e-4
#' rho = 1
#'
#' # Fitting model: pivoted
#' model = fastkrr(data = data, response = "y",
#'                 kernel = "gaussian", rho = rho, lambda = lambda, opt = "pivoted")
#'
#' # Predict
#' new_n = 50
#' new_x = matrix(runif(new_n*d, 0, 1), nrow = new_n, ncol = d)
#' new_y = as.vector(sin(2*pi*rowMeans(new_x)^3) + rnorm(new_n, 0, 0.1))
#'
#' pred = predict(model, new_x)
#' crossprod(pred - new_y) / new_n
#'
#' @importFrom stats predict
#' @export
predict.krr = function(object, newdata, ...){
  if(missing(newdata)) return(object$fitted.values)

  if(object$opt == "rff"){
    x_new = newdata
    W = object$W
    b = object$b
    coef = object$coefficients
    Z_new = make_Z(x_new, W = W, b = b)

    return(as.vector(Z_new %*% coef))
  }else{
    x_new = newdata
    x = object$x
    kernel = object$kernel
    rho = object$rho
    coef = object$coefficients

    K_new = make_kernel(x, x_new, kernel = kernel, rho = rho)

    return(as.vector(K_new %*% coef))
  }
}


#' Fit kernel ridge regression using exact or approximate methods
#'
#' This function performs kernel ridge regression (KRR) in high-dimensional
#' settings. The regularization parameter \eqn{\lambda} can be supplied by
#' the user or selected automatically using cross-validation or
#' restricted maximum likelihood (REML). For scalability,
#' three different kernel approximation strategies are supported (Nyström approximation,
#' Pivoted Cholesky decomposition, Random Fourier Features(RFF)), and kernel matrix
#' can be computed using two methods(Gaussian kernel, Laplace kerenl).
#'
#' @param data A data frame containing the data point variables and response
#' variable.
#' @param response A character string specifying the name of the response
#' variable in \code{data}.
#' @param kernel Kernel type either "gaussian"or "laplace".
#' @param rho Scaling parameter of the kernel(\eqn{\rho}),  specified by the user.
#'   Defaults to \code{1}.
#' \deqn{\text{Gaussian kernel : } \mathcal{K}(x, x') = \exp(-\rho \| x - x'\|^2_2)}
#' \deqn{\text{Laplace kernel : } \mathcal{K}(x, x') = \exp(-\rho \| x - x'\|_1)}
#' @param m Approximation rank(number of random features) used for the low-rank kernel approximation.
#'   If not provided by the user, it defaults to
#'   \deqn{\lceil n^{1/3} \cdot \log(d + 5) \rceil,}
#'   where \eqn{n = nrow(X)} and \eqn{d = ncol(X)}.
#' @param eps Tolerance parameter used only in \code{"pivoted"}
#'   for stopping criterion of the Pivoted Cholesky decomposition.
#' @param lambda Regularization parameter. If \code{NULL}, the penalty parameter
#'   is chosen automatically via \pkg{CVST} package or REML. If not provided, the argument is set to a
#'   kernel-specific grid of 100 values: \eqn{[10^{-10}, 10^{-3}]} for Gaussian,
#'   \eqn{[10^{-5}, 10^{-2}]} for Laplace.
#' @param opt Method for constructing or approximating :
#'  \describe{
#'   \item{\code{"exact"}}{Construct the full kernel matrix
#'   \eqn{K \in \mathbb{R}^{n\times n}} using design martix \eqn{X}.}
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
#' @param n_threads Number of parallel threads.
#'   The default is 4. If the system does not support 4 threads,
#'   it automatically falls back to 1 thread.
#'   Parallelization (implemented in C++) is one of the main advantages
#'   of this package and is applied only for \code{opt = "nystrom"} or \code{opt = "rff"}, and for the
#'   Laplace kernel (\code{kernel = "laplace"}).
#' @param selection_method Method used to select \eqn{\lambda} when a grid or \code{NULL} is
#'   passed. One of:
#'   \describe{
#'     \item{\code{"exactCV"}}{Full cross-validation via \pkg{CVST} (default).}
#'     \item{\code{"fastCV"}}{Accelerated sequential-testing CV via \pkg{CVST}.}
#'     \item{\code{"REML"}}{Restricted Maximum Likelihood.}
#'   }
#' @param verbose If TRUE, detailed progress and cross-validation
#' results are printed to the console. If FALSE, suppresses
#' intermediate output and only returns the final result.
#'
#' @details
#' The function performs several input checks and automatic adjustments:
#'
#' \itemize{
#'   \item \code{kernel} must be either \code{gaussian} or \code{laplace}.
#'   \item \code{opt} must be one of \code{"exact"}, \code{"pivoted"},
#'     \code{"nystrom"}, or \code{"rff"}.
#'   \item If \code{m} is \code{NULL}, it defaults to
#'     \deqn{\lceil n^{1/3} \cdot \log(d + 5) \rceil}
#'     where \eqn{n = nrow(X)} and \eqn{d = ncol(X)}.
#'     Otherwise, \code{m} must be a positive integer.
#'   \item \code{rho} must be a positive real number (default is 1).
#'   \item \code{lambda} can be specified in three ways:
#'     \enumerate{
#'       \item A positive numeric scalar, in which case the model is fitted with
#'          this single value.
#'       \item A numeric vector (length >= 3) of positive values used as a tuning grid;
#'          selection is performed by \pkg{CVST} cross-validation (sequential testing if
#'          \code{fastcv = TRUE}).
#'       \item \code{NULL}: use a default grid (internal setting) and tune \code{lambda}
#'         via \pkg{CVST} or REML}
#'
#'   \item \code{n_threads}: Number of threads for parallel computation.
#'     Default is \code{4}. If the system has <= 3 available processors,
#'     it uses \code{1}.
#'}
#'
#'
#' @return
#'
#' \itemize{
#'   \item{\code{coefficients}: Estimated coefficient vector \eqn{\mathbb{R}^{n}}.
#'   Accessible via \code{model$coefficients}.}
#'   \item{\code{fitted.values}: Fitted values \eqn{\mathbb{R}^{n}}.
#'   Accessible via \code{model$fitted.values}.}
#'   \item{\code{opt}: Kernel approximation option.
#'   One of \code{"exact"}, \code{"pivoted"}, \code{"nystrom"}, \code{"rff"}.}
#'   \item{\code{kernel}: Kernel used (\code{"gaussian"} or \code{"laplace"}).}
#'   \item{\code{x}: Input design matrix.}
#'   \item{\code{y}: Response vector.}
#'   \item{\code{lambda}: Regularization parameter. If \code{NULL}, tuned
#'   by cross-validation via \pkg{CVST} or REML.}
#'   \item{\code{rho}: Additional user-specified hyperparameter.}
#'   \item{\code{selection_method}: Tunning method for select hyperparmeter lambda}
#'   \item{\code{call}: Ahe matched function call used to create the object.}
#'   \item{\code{n_threads}: Number of threads used for parallelization.}
#' }
#'
#' Additional components depend on the value of \code{opt}:
#'
#' \subsection{opt = \dQuote{exact}}{
#' \itemize{
#'   \item{\code{K}: The full kernel matrix \eqn{K \in \mathbb{R}^{n \times n}}.}
#'   \item{\code{chol_factor}: Lower triangular Cholesky factor \eqn{L \in \mathbb{R}^{n \times n}}
#'     of \eqn{K + n\lambda I}; satisfies \eqn{K + n\lambda I = L L^\top}.}
#' }}
#'
#'
#' \subsection{opt = \dQuote{nystrom}}{
#' \itemize{
#'   \item{\code{m}: Kernel approximation degree.}
#'   \item{\code{approx_factor}: The method provides a low-rank approximation to the kernel matrix
#'     \eqn{R \in \mathbb{R}^{n \times m}} obtained via
#'     Nyström approximation; satisfies \eqn{K \approx R R^\top}.}
#' }}
#'
#' \subsection{opt = \dQuote{pivoted}}{
#' \itemize{
#'   \item{\code{m}: Kernel pproximation degree.}
#'   \item{\code{approx_factor}: The method provides a low-rank approximation to the kernel matrix
#'     \eqn{PR \in \mathbb{R}^{ n \times m}} obtained via
#'     Pivoted Cholesky decomposition; satisfies \eqn{K \approx PR\,(PR)^\top}.}
#'   \item{\code{eps}: Numerical tolerance used for early stopping in the pivoted Cholesky decomposition.}
#' }}
#'
#' \subsection{opt = \dQuote{rff}}{
#' \itemize{
#'   \item{\code{m}: Number of random features.}
#'   \item{\code{approx_factor}: Random Fourier Feature matrix \eqn{Z \in \mathbb{R}^{n \times m}} with
#'     \eqn{Z_{ij} = z_j(x_i) = \sqrt{2/m}\cos(\omega_j^\top x_i + b_j), \quad j = 1, \cdots, m,}
#'     so that \eqn{K \approx Z Z^\top}.}
#'   \item{\code{W}: Random frequency matrix \eqn{\omega \in \mathbb{R}^{m \times d}}
#'       (row \eqn{j} is \eqn{\omega_j^\top \in \mathbb{R}^d}), drawn i.i.d. from the spectral density of the chosen kernel:
#'       \itemize{
#'         \item {Gaussian: \eqn{\omega_{jk} \sim \mathcal{N}(0, 2\gamma)} (e.g., \eqn{\gamma=1/\ell^2}).}
#'         \item {Laplace: \eqn{\omega_{jk} \sim \mathrm{Cauchy}(0, 1/\sigma)} i.i.d.
#'       }}}
#'   \item{\code{b} Random phase vector \eqn{b \in \mathbb{R}^m}, i.i.d. \eqn{\mathrm{Unif}[0,\,2\pi]}.}
#' }}
#'
#'
#' @examples
#' # Data setting
#' set.seed(1)
#' lambda = 1e-4
#' d = 3
#' rho = 1
#' n = 50
#'
#' X = matrix(runif(n*d, 0, 1), nrow = n, ncol = d); colnames(X) = paste0("X", seq_len(d))
#' y = sin(2 * pi * rowMeans(X)^3) + rnorm(n, mean = 0, sd = 0.1)
#'
#' data = data.frame(X, y = y)
#'
#' # Exapmle: pivoted cholesky
#' model = fastkrr(data = data, response = "y",  kernel = "gaussian",
#'                 opt = "pivoted", rho = rho, lambda = 1e-4)
#'
#' # Example: nystrom
#' model = fastkrr(data = data, response = "y", kernel = "gaussian",
#'                 opt = "nystrom", rho = rho, lambda = 1e-4)
#'
#' # Example: random fourier features
#' model = fastkrr(data = data, response = "y", kernel = "gaussian",
#'                 opt = "rff", rho = rho, lambda = 1e-4)
#'
#' # Example: Laplace kernel
#' model = fastkrr(data = data, response = "y",  kernel = "laplace",
#'                 opt = "nystrom", n_threads = 1, rho = rho)
#'
#' @export
fastkrr = function(data, response,
                   kernel = "gaussian", # c(gaussian, laplace)
                   opt = "exact",  # c(exact, pivoted, nystrom, rff)
                   m = NULL,
                   eps = 1e-6,
                   rho = 1,
                   lambda = NULL,
                   selection_method = "exactCV",  # c(exactCV, fastCV, REML)
                   n_threads = 4,
                   verbose =  TRUE,
                   na.rm = FALSE)
{
  call = match.call()

  # Data frame
  if(!is.data.frame(data)) stop("data must be a data.frame")
  removed_row_idx = integer(0)
  if (anyNA(data)){
    na_cols = names(data)[colSums(is.na(data)) > 0]
    if(!na.rm){
      stop("Missing values (NA) found in: ",
           paste(na_cols, collapse = ", "),
           "\nSet na.rm = TRUE to remove rows containing missing values.",
           call. = FALSE)
    }else{
      removed_row_idx = which(!stats::complete.cases(data))
      data = stats::na.omit(data)
      message("Removed ", length(removed_row_idx),
              " row", if (length(removed_row_idx) != 1) "s" else "",
              " containing missing values.\n",
              "Removed row indices: ",
              paste(removed_row_idx, collapse = ", "), ".")
    }
  }
  if(!(length(response) == 1 && is.character(response) &&
       response %in% colnames(data)))
    stop("response must be a single column name in data.")

  # Kernel and approximation option
  if (!kernel %in% c("gaussian", "laplace"))
    stop("kernel must be one of 'gaussian', 'laplace'")
  if (!opt %in% c("exact", "pivoted", "nystrom", "rff"))
    stop("opt must be one of 'exact', 'pivoted', 'nystrom', 'rff'")
  if (!selection_method %in% c("exactCV", "fastCV", "REML"))
    stop("selection_method must be one of 'exactCV', 'fastCV', 'REML'")

  y = data[[response]]
  x = as.matrix(data[, setdiff(names(data), response), drop = FALSE])

  # Parameter setting
  if (eps <= 0)
    stop("eps must be a positive real number")
  if (is.null(m))
    m = as.integer(max(1, ceiling(nrow(x)^(1/3) * log(ncol(x) + 5))))
  else if(m <= 0)
    stop("m must be a positive integer")
  if (rho <= 0)
    stop("rho must be a positive real number")


  # Complexity parameter setting
  if(is.null(lambda)){
    lambda = if (kernel == "gaussian") seq(1e-10, 1e-3, len = 100) else seq(1e-5, 1e-2, len = 100) # vector
  }else if(is.vector(lambda) && all(lambda > 0) && length(lambda) >= 3){
    lambda = lambda # vector
  }else if(is.numeric(lambda) && length(lambda) == 1 && lambda > 0){
    lambda = lambda # scalar
  }else{
    stop("lambda must be a positive number or a numeric vector of positive real numbers with length greater than 3")
  }


  # Adjust number of threads
  if ((opt %in% c("nystrom", "rff")) || (kernel == "laplace") || (selection_method == "REML")){
    max_threads = get_num_procs()
    if (max_threads <= 3)
      n_threads = 1
    else
      n_threads = min(n_threads, max_threads - 1)
  }else{
    n_threads = 1
  }

  n = nrow(x)
  d = ncol(x)
  m = as.integer(m)
  #rate = round(m / n, 5)
  idx_ny = NULL  # will be set during REML nystrom selection; reused in fitting

  # Lambda selection
  if(length(lambda) > 1){
    if(selection_method == "REML"){
      lambda = if(opt == "exact"){
        K_tmp = make_kernel(x, kernel = kernel, rho = rho, n_threads = n_threads)
        reml_exact(K_tmp, y, lambda, n_threads = n_threads)
      }else if(opt == "pivoted"){
        PR_tmp = pchol_kernel(x, rho = rho, kernel = kernel, m = m, eps = eps, verbose = FALSE)$PR
        reml_lowrank(PR_tmp, y, lambda, n_threads = n_threads)
      }else if(opt == "nystrom"){
        idx_ny = sample(seq_len(n), m)
        K_nm_tmp = make_kernel(x[idx_ny, , drop = FALSE], x,
                               kernel = kernel, rho = rho, n_threads = n_threads)
        K_mm_tmp = make_kernel(x[idx_ny, , drop = FALSE],
                               kernel = kernel, rho = rho, n_threads = n_threads)
        R_tmp = nystrom_kernel(K_mm_tmp, K_nm_tmp, m_in = m, n_threads = n_threads)$R
        reml_lowrank(R_tmp, y, lambda, n_threads = n_threads)
      }else if(opt == "rff"){
        rand_tmp = rff_random(m = m, rho = rho, d = d, kernel = kernel)
        Z_tmp = make_Z(x, rand_tmp$W, rand_tmp$b, n_threads = n_threads)
        reml_lowrank(Z_tmp, y, lambda, n_threads = n_threads)
      }

    }else{
      if(selection_method == "fastCV")
        CVSTmethod = function(...) CVST::fastCV(...,
                                                setup = CVST::constructCVSTModel(),
                                                verbose = verbose)
      else
        CVSTmethod = function(...) CVST::CV(..., verbose = verbose)

      data = CVST::constructData(x, y)
      if(opt == "exact"){
        ojct = CVST::constructLearner(krr_fit_exact, krr_pred)
        param_sets = CVST::constructParams(kernel = kernel, rho = rho, lambda = lambda,
                                           n_threads = n_threads, verbose = verbose)
      }else if(opt == "pivoted"){
        ojct = CVST::constructLearner(krr_fit_pivoted, krr_pred)
        param_sets = CVST::constructParams(kernel = kernel, m = m, eps = eps,
                                           rho = rho, lambda = lambda,
                                           n_threads = n_threads, verbose = verbose)
      }else if(opt == "nystrom"){
        ojct = CVST::constructLearner(krr_fit_nystrom, krr_pred)
        param_sets = CVST::constructParams(kernel = kernel, m = m,
                                           rho = rho, lambda = lambda,
                                           n_threads = n_threads, verbose = verbose)
      }else if(opt == "rff"){
        ojct = CVST::constructLearner(krr_fit_rff, krr_pred_rff)
        param_sets = CVST::constructParams(kernel = kernel, m = m,
                                           rho = rho, lambda = lambda,
                                           n_threads = n_threads, verbose = verbose)
      }
      best_param = CVSTmethod(data, ojct, param_sets)
      lambda = best_param[[1]]$lambda
    }
  }


  # Fitting
  result_values = list()
  class(result_values) = "krr"
  result_values$removed_row_idx = removed_row_idx

  if(opt == "rff"){
    rand_set = rff_random(m = m, rho = rho, d = d, kernel = kernel)
    rslt = rff(x, y, rand_set$W, rand_set$b, lambda = lambda, n_threads = n_threads)

    result_values$coefficients = rslt$coefficients
    result_values$fitted.values = as.vector(rslt$fitted.values)
    result_values$opt = opt
    result_values$kernel = kernel
    result_values$x = x
    result_values$y = y
    result_values$lambda = lambda
    result_values$rho = rho
    result_values$n_threads = n_threads
    result_values$selection_method = selection_method
    result_values$call = call
    result_values$m = m
    result_values$approx_factor = rslt$Z
    result_values$W = rslt$W
    result_values$b = rslt$b
    return(result_values)

  }else if(opt == "exact"){
    K = make_kernel(x, kernel = kernel, rho = rho, n_threads = n_threads)

    rslt = solve_chol(K + diag(n * lambda, n), y)

    result_values$coefficients = rslt$coefficients
    result_values$fitted.values = as.vector(K %*% rslt$coefficients)
    result_values$opt = opt
    result_values$kernel = kernel
    result_values$x = x
    result_values$y = y
    result_values$lambda = lambda
    result_values$rho = rho
    result_values$n_threads = n_threads
    result_values$selection_method = selection_method
    result_values$call = call
    result_values$K = K
    result_values$chol_factor = rslt$chol_factor
    return(result_values)

  }else if(opt == "pivoted"){
    rslt = pchol(x, y, lambda = lambda, rho = rho,
                 kernel = kernel, m = m, verbose = verbose)

    result_values$coefficients = rslt$coefficients
    result_values$fitted.values = as.vector(rslt$fitted.values)
    result_values$opt = opt
    result_values$kernel = kernel
    result_values$x = x
    result_values$y = y
    result_values$lambda = lambda
    result_values$rho = rho
    result_values$n_threads = n_threads
    result_values$selection_method = selection_method
    result_values$call = call
    result_values$m = rslt$m
    result_values$approx_factor = rslt$PR
    result_values$eps = eps
    return(result_values)

  }else if(opt == "nystrom"){
    if(is.null(idx_ny)) idx_ny = sample(seq_len(n), m)

    K_nm = make_kernel(x[idx_ny, , drop = FALSE], x,
                       kernel = kernel, rho = rho, n_threads = n_threads)

    K_mm = make_kernel(x[idx_ny, , drop = FALSE],
                       kernel = kernel, rho = rho, n_threads = n_threads)

    rslt = nystrom(K_nm = K_nm, K_mm = K_mm, m = m, y = y,
                   lambda = lambda, n_threads = n_threads)

    result_values$coefficients = as.vector(rslt$coefficients)
    result_values$fitted.values = as.vector(rslt$fitted.values)
    result_values$opt = opt
    result_values$kernel = kernel
    result_values$x = x
    result_values$y = y
    result_values$lambda = lambda
    result_values$rho = rho
    result_values$n_threads = n_threads
    result_values$selection_method = selection_method
    result_values$call = call

    result_values$m = rslt$m
    result_values$approx_factor = rslt$R
    return(result_values)
  }
}

krr_fit_exact = function(data, param) {
  x = data$x
  y = data$y
  n = nrow(x)

  lambda = as.numeric(param$lambda)

  K = make_kernel(x, kernel = param$kernel, rho = param$rho, n_threads = param$n_threads)
  rslt <- solve_chol(K + diag(n * lambda, n), y)

  return(list(data = data, kernel = param$kernel,
              coefficients = rslt$coefficients,
              rho = param$rho, lambda = param$lambda,
              n_threads = param$n_threads))

}

krr_fit_pivoted = function(data, param) {
  x = data$x
  y = data$y
  n = nrow(x)
  m = as.integer(min(n, param$m))  # as.integer(n * param$rate)
  lambda = as.numeric(param$lambda)

  rslt = pchol(x, y, lambda = lambda, rho = param$rho,
               kernel = param$kernel, m = m,
               eps = param$eps, verbose = FALSE)
  coefficients = rslt$coefficients

  return(list(data = data, kernel = param$kernel,
              coefficients = coefficients,
              rho = param$rho, lambda = param$lambda,
              verbose = param$verbose,
              n_threads = param$n_threads))

}


krr_fit_nystrom = function(data, param) {
  x = data$x
  y = data$y
  n = nrow(x)
  d = ncol(x)

  m = as.integer(min(n, param$m))  # as.integer(n * param$rate)

  lambda = as.numeric(param$lambda)

  idx_ny = sample(n, m)

  K_nm = make_kernel(x[idx_ny, , drop = FALSE], x,
                     kernel = param$kernel, rho = param$rho, n_threads = param$n_threads)

  K_mm = make_kernel(x[idx_ny, , drop = FALSE], x[idx_ny, , drop = FALSE],
                     kernel = param$kernel, rho = param$rho, n_threads = param$n_threads)

  rslt = nystrom(K_nm = K_nm, K_mm = K_mm, m = m, y = y,
                 lambda = lambda, n_threads = param$n_threads)

  coefficients = rslt$coefficients

  return(list(data = data, kernel = param$kernel,
              coefficients = coefficients,
              rho = param$rho, lambda = param$lambda,
              n_threads = rslt$n_threads))

}

krr_fit_rff = function(data, param) {
  x = data$x
  y = data$y
  n = nrow(x)
  d = ncol(x)

  lambda = as.numeric(param$lambda)
  rand_set = rff_random(m = as.integer(min(n, param$m)), d = d, param$rho, kernel = param$kernel)

  rslt = rff(x, y, rand_set$W, rand_set$b, lambda = lambda, n_threads = param$n_threads)
  coefficients = rslt$coefficients

  return(list(data = data,
              kernel = param$kernel,
              coefficients = rslt$coefficients,
              Z = rslt$Z,
              m = rslt$m,
              W = rslt$W,
              b = rslt$b,
              rho = param$rho, lambda = param$lambda,
              n_threads = rslt$n_threads))

}


krr_pred = function(model, newdata) {
  x = model$data$x
  x_new = newdata$x
  rho = model$rho
  kernel = model$kernel
  opt = model$opt

  K_new = make_kernel(x, x_new, kernel = kernel, rho = rho, n_threads = model$n_threads)


  return(as.vector(K_new %*% model$coefficients))
}


krr_pred_rff = function(model, newdata) {
  x_new = newdata$x
  W = model$W
  b = model$b
  coef = model$coefficients
  Z_new = make_Z(x_new, W = W, b = b)

  return(as.vector(Z_new %*% coef))
}
