// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// K(i,j) = exp(-rho * ||x_i - x_j||^2)
inline double gaussian_entry(const arma::mat& X, int i, int j, double rho) {
  arma::rowvec diff = X.row(i) - X.row(j);
  return std::exp(-rho * arma::accu(diff % diff));
}

// K(i,j) = exp(-rho * ||x_i - x_j||_1)
inline double laplace_entry(const arma::mat& X, int i, int j, double rho) {
  arma::rowvec diff = X.row(i) - X.row(j);
  return std::exp(-rho * arma::accu(arma::abs(diff)));
}

inline double kernel_entry(const arma::mat& X, int i, int j,
                           double rho, const std::string& kernel) {
  if (kernel == "gaussian") return gaussian_entry(X, i, j, rho);
  if (kernel == "laplace")  return laplace_entry(X, i, j, rho);
  Rcpp::stop("Unknown kernel: use 'gaussian' or 'laplace'");
  return 0.0;
}

// [[Rcpp::export]]
Rcpp::List pchol_kernel(const arma::mat& X,
                        double rho = 1.0,
                        std::string kernel = "gaussian",
                        Rcpp::Nullable<int> m = R_NilValue,
                        double eps = 1e-6,
                        bool verbose = true) {
  int n = X.n_rows;
  int r = m.isNull() ? n : Rcpp::as<int>(m);
  if (!m.isNull() && (r <= 0 || r > n)) {
    Rcpp::Rcout << "Reset m=n\n";
    r = n;
  }

  arma::mat R(n, r, arma::fill::zeros);
  arma::ivec piv(n);
  for (int i = 0; i < n; ++i) piv(i) = i;

  // K(i,i) = exp(0) = 1
  arma::vec schur_diag(n, arma::fill::ones);
  int rank = 0;

  for (int i = 0; i < r; ++i) {
    double sqsum = arma::accu(schur_diag.subvec(i, n - 1));
    if (sqsum <= eps) break;

    arma::uword max_idx = i + schur_diag.subvec(i, n - 1).index_max();
    if (max_idx != static_cast<arma::uword>(i)) {
      std::swap(piv(i), piv(max_idx));
      std::swap(schur_diag(i), schur_diag(max_idx));
      for (int k = 0; k < i; ++k)
        std::swap(R(i, k), R(max_idx, k));
    }

    double pivot_val = schur_diag(i);
    if (pivot_val <= 0.0) break;
    R(i, i) = std::sqrt(pivot_val);
    rank++;

    for (int j = i + 1; j < n; ++j) {
      double kij = kernel_entry(X, piv(j), piv(i), rho, kernel);

      double dot = 0.0;
      for (int k = 0; k < i; ++k)
        dot += R(j, k) * R(i, k);

      R(j, i) = (kij - dot) / R(i, i);
      schur_diag(j) -= R(j, i) * R(j, i);
      schur_diag(j)  = std::max(schur_diag(j), 0.0);
    }
  }

  arma::mat R_final = R.cols(0, rank - 1);
  arma::mat PR(n, rank);
  for (int i = 0; i < n; ++i)
    PR.row(piv(i)) = R_final.row(i);

  if (!m.isNull() && rank < r && verbose)
    Rcpp::Rcout << "Note: m truncated at " << rank << " due to early termination\n";

  return Rcpp::List::create(
    Rcpp::Named("PR")     = PR,
    Rcpp::Named("rank")   = rank,
    Rcpp::Named("eps")    = eps,
    Rcpp::Named("kernel") = kernel
  );
}

// [[Rcpp::export]]
Rcpp::List pchol(const arma::mat& X,
                 const arma::vec& y,
                 double lambda,
                 double rho = 1.0,
                 std::string kernel = "gaussian",
                 Rcpp::Nullable<int> m = R_NilValue,
                 double eps = 1e-6,
                 bool verbose = true) {
  int n = X.n_rows;
  Rcpp::List kernel_res = pchol_kernel(X, rho, kernel, m, eps, verbose);
  arma::mat PR = kernel_res["PR"];
  int rank     = kernel_res["rank"];

  arma::mat PRtPR = PR.t() * PR;
  PRtPR.diag() += lambda * n;
  arma::vec tmp = arma::solve(PRtPR, PR.t() * y, solve_opts::likely_sympd);
  arma::vec coef_hat = y / (n * lambda) - PR * tmp / (n * lambda);

  return Rcpp::List::create(
    Rcpp::Named("PR")            = PR,
    Rcpp::Named("m")             = rank,
    Rcpp::Named("coefficients")  = coef_hat,
    Rcpp::Named("fitted.values") = PR * PR.t() * coef_hat,
    Rcpp::Named("eps")           = eps,
    Rcpp::Named("kernel")        = kernel
  );
}
