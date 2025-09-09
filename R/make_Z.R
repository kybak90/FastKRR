#'  Random Fourier Feature matrix \eqn{Z} construction for given datasets
#'
#' Constructs a Random Fourier Feature matrix \eqn{Z \in \mathbb{R}^{n \times m}} with
#' \eqn{Z_{ij} = z_j(x_i) = \sqrt{2/m}\cos(\omega_j^\top x_i + b_j), \quad j = 1, \cdots, m,
#'  \quad i = 1, \cdots, n}.
#' Implemented in C++ via RcppArmadillo.
#'
#' @param X Design matrix \eqn{X \in \mathbb{R}^{n \times d}} (rows  \eqn{x_i \in \mathbb{R}^d}).
#' @param W  Random frequency matrix \eqn{\omega \in \mathbb{R}^{m \times d}}
#'       (row \eqn{j} is \eqn{\omega_j^\top \in \mathbb{R}^d}), drawn i.i.d. from the spectral density of the chosen kernel:
#'       \itemize{
#'         \item Gaussian: \eqn{\omega_{jk} \sim \mathcal{N}(0, 2\gamma)} (e.g., \eqn{\gamma=1/\ell^2}).
#'         \item Laplace: \eqn{\omega_{jk} \sim \mathrm{Cauchy}(0, 1/\sigma)} i.i.d.
#'       }
#' @param b Random phase vector \eqn{b \in \mathbb{R}^m}, i.i.d. \eqn{\mathrm{Unif}[0,\,2\pi]}.
#'
#' @return Random Fourier Feature matrix \eqn{Z}
#'
#'
#' @examples
#' set.seed(1)
#' lambda = 1e-4
#' d = 1
#' rho = 1
#' n = 50
#' X = matrix(runif(n*d, 0, 1), nrow = n, ncol = d)
#' m = ceiling(n* log(d + 5)/ 10)
#' y = as.vector(sin(2*pi*rowMeans(X)^3) + rnorm(n, 0, 0.1))
#' rv = rff_random(m = m, d = d, rho = 1, kernel = "gaussian")
#' str(rv)
#'
#' Z = make_Z(X, rv$W, rv$b)
#' str(Z)
#'
#' @name make_Z
#' @export make_Z
#' @useDynLib FastKRR, .registration = TRUE
#' @importFrom Rcpp evalCpp
NULL
