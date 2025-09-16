#'  Kernel matrix \eqn{K} construction for given datasets
#'
#' Constructs a kernel matrix \eqn{K \in \mathbb{R}^{n \times n'}} given two
#' datasets \eqn{X \in \mathbb{R}^{n \times d}} and \eqn{X' \in \mathbb{R}^{n' \times d}},
#' where \eqn{x_i \in \mathbb{R}^d} and \eqn{x'_j \in \mathbb{R}^d} denote the
#' i-th and j-th rows of \eqn{X} and \eqn{X'}, respectively, and
#' \eqn{K_{ij}=\mathcal{K}(x_i, x'_j)} for a user-specified kernel.
#' Implemented in C++ via RcppArmadillo.
#'
#' @param X Design matrix \eqn{X \in \mathbb{R}^{n \times d}} (rows  \eqn{x_i \in \mathbb{R}^d}).
#' @param X_new  Second matrix \eqn{X' \in \mathbb{R}^{n' \times d}} (rows  \eqn{x'_j \in \mathbb{R}^d}).
#' If omitted, \eqn{X' = X} and \eqn{n' = n}.
#' @param kernel Kernel type; one of \code{"gaussian"} or \code{"laplace"}.
#' @param rho Kernel width parameter (\eqn{\rho > 0}).
#' @param n_threads Number of parallel threads.
#'   The default is 4. If the system does not support 4 threads,
#'   it automatically falls back to 1 thread.
#'   Parallelization (implemented in C++) is one of the main advantages
#'   of this package and is applied only for \code{"laplace"} kernels.
#'
#' @details
#' Gaussian: \deqn{\mathcal{K}(x_i,x_j)=\exp\!\big(-\rho\|x_i-x_j\|_2^2\big)}
#' Laplace:  \deqn{\mathcal{K}(x_i,x_j)=\exp\!\big(-\rho\|x_i-x_j\|_1\big)}
#'
#' @return
#' If \code{X_new} is \code{NULL}, a symmetric matrix \eqn{K_{ij} = \mathcal{K}(x_i,x_j), K \in \mathbb{R}^{n \times n}}.
#' Otherwise, a matrix \eqn{K'_{ij} = \mathcal{K}(x_i,x'_j), K' \in \mathbb{R}^{n' \times n}}.
#'
#' @examples
#' set.seed(1)
#' lambda = 1e-4
#' d = 1
#' rho = 1
#' n = 1000
#' X = matrix(runif(n*d, 0, 1), nrow = n, ncol = d)
#'
#' # second matrix
#' new_n = 1500
#' new_X = matrix(runif(new_n*d, 0, 1), nrow = new_n, ncol = d)
#'
#' # make kernel : Gaussian kernel
#' K = make_kernel(X, kernel = "gaussian", rho = rho)
#' new_K = make_kernel(X, new_X, kernel = "gaussian", rho = rho)
#'
#' # make kernel : Laplace kernel
#' K = make_kernel(X, kernel = "laplace", rho = rho, n_threads = 1)
#' new_K = make_kernel(X, new_X, kernel = "laplace", rho = rho, n_threads = 1)
#'
#' @name make_kernel
#' @export make_kernel
#' @useDynLib FastKRR, .registration = TRUE
#' @importFrom Rcpp evalCpp
NULL
