#' Compute low-rank approximations(Nyström, Pivoted Cholesky, RFF)
#'
#' This function compute low-rank kernel approximation \eqn{\tilde{K} \in \mathbb{R}^{n \times n}}using three methods:
#' Nyström approximation, Pivoted Cholesky decomposition, and
#' Random Fourier Features (RFF).
#'
#' @param K Exact Kernel matrix \eqn{K \in \mathbb{R}^{n \times n}}. Used in \code{"nystrom"} and \code{"pivoted"}.
#' @param X Design matrix \eqn{X \in \mathbb{R}^{n \times d}}. Only required for \code{"rff"}.
#' @param opt Method for constructing or approximating :
#'  \describe{
#'   \item{\code{"nystrom"}}{Construct a low-rank approximation of
#'       the kernel matrix \eqn{K \in \mathbb{R}^{n \times n}}
#'       using the Nyström approximation.}
#'   \item{\code{"pivoted"}}{Construct a low-rank approximation of
#'       the kernel matrix \eqn{K \in \mathbb{R}^{n \times n}}
#'       using Pivoted Cholesky decomposition.}
#'   \item{\code{"rff"}}{Construct a low-rank approximation of
#'       the kernel matrix \eqn{K \in \mathbb{R}^{n \times n}}
#'       using Random Fourier Features (RFF).}
#'  }
#' @param kernel Kernel type either "gaussian"or "laplace".
#' @param m Approximation rank (number of random features) for the
#'   low-rank kernel approximation. If not specified, the recommended
#'   choice is
#'   \deqn{\lceil n \cdot \log(d + 5) / 10 \rceil}
#'   where \eqn{X} is design matrix, \eqn{n = nrow(X)} and \eqn{d = ncol(X)}.
#' @param d Design matrix's dimension (\eqn{d = ncol(X)}).
#' @param rho Scaling parameter of the kernel (\eqn{\rho}), specified by the user.
#' @param eps Tolerance parameter used only in \code{"pivoted"}
#'   for stopping criterion of the Pivoted Cholesky decomposition.
#' @param W Random frequency matrix \eqn{\omega \in \mathbb{R}^{m \times d}}
#' @param b Random phase vector \eqn{b \in \mathbb{R}^m}, i.i.d. \eqn{\mathrm{Unif} [ 0,\,2\pi ]}.
#' @param n_threads Number of parallel threads.
#'   The default is 4. If the system does not support 4 threads,
#'   it automatically falls back to 1 thread. It is applied only for \code{opt = "nystrom"} or \code{opt = "rff"}
#'   , and for the Laplace kernel (\code{kernel = "laplace"}).
#'
#' @details
#' Requirements and what to supply:
#'
#' \describe{
#'   \item{Common}{
#'
#'     \itemize{
#'       \item \code{d} and \code{rho} must be provided (non-\code{NULL}).
#'     }
#'   }
#'
#'   \item{\code{nystrom} / \code{pivoted}}{
#'
#'     \itemize{
#'       \item Require a precomputed kernel matrix \code{K}; error if \eqn{K} is \code{NULL}.
#'       \item If \code{m} is \code{NULL}, use \eqn{\lceil n \cdot \log(d + 5) / 10 \rceil}.
#'       \item For \code{"pivoted"}, a tolerance \code{eps} is used; the decomposition stops early
#'       when the next pivot (residual diagonal) drops below \code{eps}.
#'     }
#'   }
#'
#'   \item{\code{rff}}{
#'
#'     \itemize{
#'       \item \eqn{K} must be \code{NULL} (not used) and \code{X} must be provided with \code{d = ncol(X)}.
#'       \item To pre-supply random features, provide both
#'             \code{W} (random frequency matrix \eqn{\omega \in \mathbb{R}^{m \times d}}) and
#'             \code{b} (random phase vector \eqn{b \in \mathbb{R}^{m}}); error if only one is supplied or shapes mismatch.
#'       \item When \code{W} and \code{b} are supplied by \code{\link{rff_random}}.
#'     }
#'   }
#' }
#'
#' @return
#' A list containing the results of the low-rank kernel approximation.
#' The exact structure depends on the chosen method:
#'
#' \describe{
#'   \item{\code{"nystrom"}}{
#'     \itemize{
#'       \item \code{K_approx}: Approximated kernel matrix (\eqn{K \in \mathbb{R}^{n \times n}})
#'       from the Nyström approximation.
#'       \item \code{m}: Approximation rank used for the low-rank kernel approximation.
#'       \item \code{n_threads}: Number of threads used in the computation.
#'       \item \code{method}: The kernel approximation method actually used. The string \code{"nystrom"}.
#'     }
#'   }
#'
#'   \item{\code{"pivoted"}}{
#'     \itemize{
#'       \item \code{K_approx}: Approximated kernel matrix (\eqn{K \in \mathbb{R}^{n \times n}})
#'             from the Pivoted Cholesky decomposition.
#'       \item \code{m}: Effective rank actually used.
#'             This value is at most the requested \code{m} and
#'             may be smaller if early stopping is triggered by \code{eps}
#'       \item \code{method}: The kernel approximation method actually used. The string \code{"pivoted"}.
#'     }
#'   }
#'
#'   \item{\code{"rff"}}{
#'     \itemize{
#'       \item \code{K_approx}: Approximated kernel matrix (\eqn{K \in \mathbb{R}^{n \times n}}).
#'       \item \code{method}: The kernel approximation method actually used. The string \code{"rff"}.
#'       \item \code{m}: Number of random features used.
#'       \item \code{d}: Input design matrix's dimension.
#'       \item \code{rho}: Scaling parameter of the kernel.
#'       \item \code{W}: Random frequency matrix (\eqn{m \times d}).
#'       \item \code{b}: Random phase vector (\eqn{m}).
#'       \item \code{used_supplied_Wb}: Logical; \code{TRUE} if user-supplied
#'             \code{W}, \code{b} were used, \code{FALSE} otherwise.
#'       \item \code{n_threads}: Number of threads used in the computation.
#'     }
#'   }
#' }
#'
#'
#'
#' @examples
#' # data setting
#' set.seed(1)
#' d = 1
#' n = 1000
#' m = 50
#' X = matrix(runif(n*d, 0, 1), nrow = n, ncol = d)
#' y = as.vector(sin(2*pi*rowMeans(X)^3) + rnorm(n, 0, 0.1))
#' K = make_kernel(X, kernel = "gaussian", rho = 1)
#'
#' # Example: RFF approximation
#' rff_pars = rff_random(m = m, d = d, rho = 1, kernel = "gaussian")
#' K_rff = approx_kernel(X = X, opt = "rff", kernel = "gaussian",
#'                       m = m, d = d, rho = 1,
#'                       W = rff_pars$W, b = rff_pars$b,
#'                       n_threads = 1)
#'
#' # Exapmle: Nystrom approximation
#' K_nystrom = approx_kernel(K = K, opt = "nystrom",
#'                           m = m, d = d, rho = 1,
#'                           n_threads = 1)
#'
#' # Example: Pivoted Cholesky approximation
#' K_pivoted = approx_kernel(K = K, opt = "pivoted",
#'                           m = m, d = d, rho = 1)
#' @export
approx_kernel = function(K = NULL, X = NULL,
                         opt = c("nystrom", "pivoted", "rff"),
                         kernel = c("gaussian", "laplace"),
                         m = NULL, d, rho, eps = 1e-6,
                         W = NULL, b = NULL, n_threads = 4) {
  call   = match.call()
  opt    = match.arg(opt)
  kernel = match.arg(kernel)


  # Adjust number of threads only for heavy computation cases
  if ((opt %in% c("nystrom", "rff")) || (kernel == "laplace")){
    max_threads = get_num_procs()
    if (max_threads <= 3)
      n_threads = 1
    else
      n_threads = min(n_threads, max_threads - 1)
  }else{
    n_threads = 1
  }

  if (missing(d) || is.null(d))   stop("'d' must be provided and non-NULL.")
  if (missing(rho) || is.null(rho)) stop("'rho' must be provided and non-NULL.")


  result_values = list()
  class(result_values) = "approx_kernel"
  attr(result_values, "call") = call
  attr(result_values, "opt") = opt

  if (opt %in% c("nystrom", "pivoted")) {
    if (is.null(K)) stop("For opt='", opt, "', argument 'K' must be provided (not NULL).")
    n = nrow(K)
    if (is.null(m)) m = n / 10 * log(d + 5)
    m = as.integer(max(1, min(n, floor(m))))

    if (opt == "nystrom") {
      rslt = nystrom_kernel(K, m, n_threads = n_threads)
      attr(result_values, "K_approx") = rslt$K_approx
      class(attr(result_values, "K_approx")) = "kernel_matrix"
      attr(result_values, "m") = rslt$m
      attr(result_values, "n_threads") = rslt$n_threads

      return(result_values)
    } else {
      rslt = pchol_kernel(K, m, eps = eps)
      attr(result_values, "K_approx") = rslt$K_approx
      class(attr(result_values, "K_approx")) = "kernel_matrix"
      attr(result_values, "m") = rslt$rank
      attr(result_values, "eps") = rslt$eps
      return(result_values)
    }
  }

  if (!is.null(K)) stop("For opt='rff', argument 'K' must be NULL (not used).")
  if (is.null(X))  stop("For opt='rff', argument 'X' must be provided (not NULL).")

  n   = nrow(X)
  d_x = ncol(X)
  if (d != d_x) stop("For opt='rff', provided d must equal ncol(X).")

  if (!is.null(W) || !is.null(b)) {
    if (is.null(W) || is.null(b)) stop("For opt='rff', provide both 'W' and 'b' or neither.")
    if (!is.matrix(W)) stop("'W' must be a numeric matrix (m x d).")
    if (ncol(W) != d) stop("ncol(W) must equal d (and ncol(X)).")
    if (length(b) != nrow(W)) stop("length(b) must equal nrow(W).")

    m_used = nrow(W)
    Z = make_Z(X, W, b, n_threads = n_threads )
    K_approx = tcrossprod(Z)

    attr(result_values, "K_approx") = tcrossprod(Z)
    class(attr(result_values, "K_approx")) = "kernel_matrix"
    attr(result_values, "m") = m_used
    attr(result_values, "d") = d
    attr(result_values, "rho") = rho
    attr(result_values, "W") = W
    attr(result_values, "b") = b
    attr(result_values, "used_supplied_Wb") = TRUE
    attr(result_values, "n_threads") = n_threads

    return(result_values)
  }

  if (is.null(m)) m = n / 10 * log(d + 5)
  m = as.integer(max(1, floor(m)))

  old_d_exists = exists("d", envir = .GlobalEnv, inherits = FALSE)
  if (old_d_exists) old_d = get("d", envir = .GlobalEnv)
  on.exit({
  }, add = TRUE)

  rb = rff_random(m = m, rho = rho, d = d, kernel = kernel)
  W = rb$W; b = rb$b
  Z = make_Z(X, W, b, n_threads = n_threads)

  attr(result_values, "K_approx") = tcrossprod(Z)
  class(attr(result_values, "K_approx")) = "kernel_matrix"
  attr(result_values, "m") = m_used
  attr(result_values, "d") = d
  attr(result_values, "rho") = rho
  attr(result_values, "W") = W
  attr(result_values, "b") = b
  attr(result_values, "used_supplied_Wb") = FALSE
  attr(result_values, "n_threads") = n_threads

  return(result_values)
}
