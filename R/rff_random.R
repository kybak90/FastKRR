#' Generate random Fourier features parameters
#'
#' This function generates random frequency matrix \eqn{\omega} and offsets
#' random phase vector \eqn{b}
#' used in the Random Fourier Features (RFF) method.
#
#' @details
#' \strong{Notation.} Let \eqn{X \in \mathbb{R}^{n \times d}} denote the design matrix
#' of inputs on which RFF will be applied downstream (with \eqn{n = nrow(X)},
#' \eqn{d = ncol(X)}). A commonly recommended choice for the number of random
#' features is
#' \deqn{\left\lceil n \cdot \frac{\log(d + 5)}{10} \right\rceil.}
#'
#'
#' @param m Number of random features \eqn{m} used in the RFF
#'   approximation. If not specified, the recommended choice is
#'   \deqn{\left\lceil n \cdot \frac{\log(d + 5)}{10} \right\rceil}
#'   where \eqn{X} is design matrix, \eqn{n = nrow(X)} and \eqn{d = ncol(X)}.
#' @param d Design matrix's dimension (\eqn{d = ncol(X)}).
#' @param rho Kernel width parameter (\eqn{\rho}), specified by the user.
#'   Controls the scale of the kernel function. Defaults to \code{1}.
#' @param kernel Kernel matrix \eqn{K \in \mathbb{R}^{n \times n}} has two kinds of Kernel ("gaussian", "laplace").
#'
#' @return A list with components:
#' \describe{
#'   \item{w}{Random frequency matrix \eqn{\omega \in \mathbb{R}^{d \times m}},
#'     sampled i.i.d. from a Cauchy distribution.}
#'   \item{b}{Random phase vector \eqn{b \in \mathbb{R}^m}, i.i.d. \eqn{\mathrm{Unif} [ 0,\,2\pi ]}.}
#' }
#'
#' @examples
#' set.seed(1)
#' lambda = 1e-4
#' d = 1
#' rho = 1
#' n = 1000
#' X = matrix(runif(n*d, 0, 1), nrow = n, ncol = d)
#' m = ceiling(n* log(d + 5)/ 10)
#' y = as.vector(sin(2*pi*rowMeans(X)^3) + rnorm(n, 0, 0.1))
#' rv = rff_random(m = m, d = d, rho = 1, kernel = "gaussian")
#' str(rv)
#'
#' @export
rff_random = function(m, d = d, rho = 1, kernel = "gaussian"){
  W = matrix(0, nrow = m, ncol = d)
  b = rep(0, m)

  if(kernel == 'gaussian'){
    W = matrix(rnorm(m * d, 0, sqrt(2 * rho)), nrow = m)
  }else if(kernel == "laplace"){
    W = matrix(rcauchy(m * d, location = 0, scale = 1/sqrt(2 * rho)), nrow = m)
  }

  b = runif(m, 0, 2 * pi)

  return(list("W" = W,
              "b" = b))
}


rff_kernel = function(X,
                      W = NULL,
                      b = NULL,
                      m = NULL,
                      rho = 1,
                      kernel = "gaussian")
{
  n = nrow(X)
  d = ncol(X)

  # m default: n/10 (min: 1)
  if (is.null(m)) m = max(1, floor(n / 10))

  if (is.null(W) || is.null(b)) {
    # Temporarily set global variable `d` because rff.random references it
    old_d_exists = exists("d", envir = .GlobalEnv, inherits = FALSE)
    if (old_d_exists) old_d = get("d", envir = .GlobalEnv)
    on.exit({
    }, add = TRUE)

    rb = rff_random(m = m, d = d, rho = rho, kernel = kernel)
    W = rb$W
    b = rb$b
  } else {

    if (ncol(W) != d) stop("ncol(W) must equal ncol(X).")
    if (length(b) != nrow(W)) stop("length(b) must equal nrow(W).")
    m = nrow(W)
  }

  Z = make_Z(X, W, b)

  K_approx = Z %*% t(Z)

  return(list(
    K_approx = K_approx,
    m = m,
    W = W,
    b = b
  ))
}
