#' Compute low-rank approximations (Nyström, Pivoted Cholesky, RFF)
#'
#' Computes low-rank kernel approximation \eqn{\tilde{K} \in \mathbb{R}^{n \times n}} using three methods:
#' Nyström approximation, Pivoted Cholesky decomposition, and
#' Random Fourier Features (RFF).
#'
#' @param X Design matrix \eqn{X \in \mathbb{R}^{n \times d}}.
#' @param opt Method for constructing or approximating:
#'   \describe{
#'     \item{\code{"nystrom"}}{Construct a low-rank approximation of
#'         the kernel matrix \eqn{K \in \mathbb{R}^{n \times n}}
#'         using the Nyström approximation.}
#'     \item{\code{"pivoted"}}{Construct a low-rank approximation of
#'         the kernel matrix \eqn{K \in \mathbb{R}^{n \times n}}
#'         using Pivoted Cholesky decomposition.}
#'     \item{\code{"rff"}}{Construct a low-rank approximation of
#'         the kernel matrix \eqn{K \in \mathbb{R}^{n \times n}}
#'         using Random Fourier Features (RFF).}
#'   }
#' @param kernel Kernel type either "gaussian" or "laplace".
#' @param m Approximation rank (number of random features) for the
#'   low-rank kernel approximation. If not specified, the recommended
#'   choice is
#'   \eqn{\lceil n^{1/2} \cdot \log(d + 5) \rceil}
#'   where \eqn{X} is design matrix, \eqn{n = nrow(X)} and \eqn{d = ncol(X)}.
#' @param rho Scaling parameter of the kernel (\eqn{\rho}), specified by the user.
#' @param eps Tolerance parameter used only in \code{"pivoted"}
#'   for stopping criterion of the Pivoted Cholesky decomposition.
#'   If \code{NULL}, it is dynamically adjusted based on the smoothness of the chosen kernel:
#'   defaults to \code{1e-6} for the infinitely smooth Gaussian kernel to ensure precision,
#'   and \code{1e-4} for the less smooth Laplace kernel.
#' @param W Random frequency matrix \eqn{\omega \in \mathbb{R}^{m \times d}}
#' @param b Random phase vector \eqn{b \in \mathbb{R}^m}, i.i.d. \eqn{\mathrm{Unif} [ 0, 2\pi ]}.
#' @param n_threads Number of parallel threads.
#'   If \code{NULL}, defaults to half of the available system processors.
#'   It automatically falls back to 1 thread if the system has 3 or fewer processors.
#'   Applied only for \code{opt = "nystrom"}, \code{opt = "rff"}, or \code{kernel = "laplace"}.
#'
#' @details
#' Requirements and what to supply:
#'
#' \strong{Common}
#'
#' \itemize{
#'   \item \code{rho} must be provided (non-\code{NULL}).
#' }
#'
#' \strong{nystrom / pivoted}
#'
#' \itemize{
#'   \item If \code{m} is \code{NULL}, use \eqn{\lceil n^{1/2} \cdot \log(d + 5) \rceil}.
#'   \item For \code{"pivoted"}, a tolerance \code{eps} is used; the decomposition stops early
#'   when the next pivot (residual diagonal) drops below \code{eps}.
#' }
#'
#' \strong{rff}
#'
#' \itemize{
#'   \item The function automatically generates
#'         \code{W} (random frequency matrix \eqn{\omega \in \mathbb{R}^{m \times d}}) and
#'         \code{b} (random phase vector \eqn{b \in \mathbb{R}^{m}}).
#'   \item If the user provides them manually, both \code{W} and \code{b} must be specified and their dimensions must be compatible.
#' }
#'
#' @return
#' \itemize{
#'   \item \code{call}: The matched function call used to create the object.
#'   \item \code{opt}: The kernel approximation method actually used (\code{"nystrom", "pivoted", "rff"}).
#'   \item \code{K_approx}: The fully reconstructed \eqn{n \times n} approximated kernel matrix.
#'   \item \code{approx_factor}: Low-rank component matrix (\eqn{R} for Nyström, \eqn{PR} for Pivoted Cholesky, \eqn{Z} for RFF).
#'   \item \code{m}: Kernel approximation degree.
#'   \item \code{rho}: Scaling parameter of the kernel.
#' }
#'
#' Additional components depend on the value of opt:
#'
#' \strong{nystrom}
#' \itemize{
#'   \item \code{n_threads}: Number of threads used in the computation.
#' }
#'
#' \strong{pivoted}
#' \itemize{
#'   \item \code{eps}: Numerical tolerance used for early stopping in the
#'                      pivoted Cholesky decomposition.
#' }
#'
#' \strong{rff}
#' \itemize{
#'   \item \code{d}: Input design matrix's dimension.
#'   \item \code{W}: \eqn{m \times d} Random frequency matrix.
#'   \item \code{b}: Random phase \eqn{m}-vector.
#'   \item \code{used_supplied_Wb}: Logical; \code{TRUE} if user-supplied
#'         \code{W}, \code{b} were used, \code{FALSE} otherwise.
#'   \item \code{n_threads}: Number of threads used in the computation.
#' }
#'
#' @examples
#' # Data setting
#' set.seed(1)
#' d = 1
#' n = 100
#' m = 50
#' X = matrix(runif(n*d, 0, 1), nrow = n, ncol = d)
#'
#' # Example: Nystrom approximation
#' K_nystrom = approx_kernel(X = X, opt = "nystrom",
#'                            m = m, rho = 1, n_threads = 1)
#'
#' # Example: Pivoted Cholesky approximation
#' K_pivoted = approx_kernel(X = X, opt = "pivoted",
#'                            m = m, rho = 1)
#' # Example: RFF approximation
#' K_rff = approx_kernel(X = X, opt = "rff", kernel = "gaussian",
#'                        m = m, rho = 1, n_threads = 1)
#'
#' @export
approx_kernel = function(X = NULL,
                         opt = c("nystrom", "pivoted", "rff"),
                         kernel = c("gaussian", "laplace"),
                         m = NULL, rho, eps = NULL,
                         W = NULL, b = NULL, n_threads = NULL) {
  call   = match.call()
  opt    = match.arg(opt)
  kernel = match.arg(kernel)

  if (is.null(X)) stop("Argument 'X' must be provided (not NULL).", call. = FALSE)
  if (!is.matrix(X)) stop("Argument 'X' must be a numeric matrix.", call. = FALSE)
  d = ncol(X)

  # Adjust eps for opt=="pivoted"
  if (opt == "pivoted") {
    if (is.null(eps)) {
      eps = if (kernel == "gaussian") 1e-6 else 1e-4
    } else {
      if (eps <= 0) stop("eps must be a positive real number", call. = FALSE)
    }
  }

  # Adjust number of threads only for heavy computation cases
  max_threads = get_num_procs()
  if (is.null(n_threads)) n_threads = as.integer(max(1, max_threads %/% 2))
  if ((opt %in% c("nystrom", "rff")) || (kernel == "laplace")){
    n_threads = if (max_threads <= 3) 1 else min(n_threads, max_threads - 1)
  }else{
    n_threads = 1
  }

  if (missing(rho) || is.null(rho)) stop("'rho' must be provided and non-NULL.", call. = FALSE)


  result_values = list()
  class(result_values) = "approx_kernel"
  result_values$call = call
  result_values$opt = opt
  result_values$kernel = kernel

  if (opt %in% c("nystrom", "pivoted")) {

    n = nrow(X)
    if (is.null(m)) m = n^(1/2) * log(d + 5)
    m = as.integer(max(1, min(n, floor(m))))

    if (opt == "nystrom") {
      idx_ny = sample(seq_len(nrow(X)), m)

      K_nm = make_kernel(X[idx_ny, , drop = FALSE], X,
                         kernel = kernel, rho = rho, n_threads = n_threads)

      K_mm = make_kernel(X[idx_ny, , drop = FALSE], X[idx_ny, , drop = FALSE],
                         kernel = kernel, rho = rho, n_threads = n_threads)

      rslt = nystrom_kernel(K_nm = K_nm, K_mm = K_mm, m, n_threads = n_threads)

      result_values$K_approx = tcrossprod(rslt$R)
      result_values$approx_factor = rslt$R
      result_values$m = rslt$m
      result_values$rho = rho
      result_values$n_threads = rslt$n_threads

      return(result_values)
    } else {
      rslt = rslt = pchol_kernel(X, rho = rho, kernel = kernel, m = m, eps = eps)

      result_values$K_approx = tcrossprod(rslt$PR)
      result_values$approx_factor = rslt$PR
      result_values$m = rslt$rank
      result_values$rho = rho
      result_values$eps = rslt$eps
      return(result_values)
    }
  }



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

    result_values$approx_factor = Z
    result_values$K_approx = K_approx
    result_values$m = m_used
    result_values$d = d
    result_values$rho = rho
    result_values$W = W
    result_values$b = b
    result_values$used_supplied_Wb = TRUE
    result_values$n_threads = n_threads

    return(result_values)
  }

  if (is.null(m)) m = n^(1/2) * log(d + 5)
  m = as.integer(max(1, floor(m)))

  old_d_exists = exists("d", envir = .GlobalEnv, inherits = FALSE)
  if (old_d_exists) old_d = get("d", envir = .GlobalEnv)
  on.exit({
  }, add = TRUE)

  rb = rff_random(m = m, rho = rho, d = d, kernel = kernel)
  W = rb$W; b = rb$b
  Z = make_Z(X, W, b, n_threads = n_threads)

  result_values$approx_factor = Z
  result_values$K_approx = tcrossprod(Z)
  result_values$m = m
  result_values$d = d
  result_values$rho = rho
  result_values$W = W
  result_values$b = b
  result_values$used_supplied_Wb = FALSE
  result_values$n_threads = n_threads

  return(result_values)
}
