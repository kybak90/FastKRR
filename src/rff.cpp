// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List rff(const arma::mat& X,
         const arma::vec& y,
         const arma::mat& W,
         const arma::vec& b,
         double lambda,
         int n_threads = 4) {

  int max_threads = 1;
  #ifdef _OPENMP
    max_threads = omp_get_num_procs();
  #endif
  if (max_threads <= 3)
    n_threads = 1;
  else
    n_threads = std::min(n_threads, max_threads - 1);


  int n = X.n_rows;
  int m = W.n_rows;


  arma::mat Z(n, m, arma::fill::zeros);
  double sqrt_2m = sqrt(2.0 / m);
#ifdef _OPENMP
#pragma omp parallel for num_threads(n_threads)
#endif
  for (int k = 0; k < m; ++k) {
    Z.col(k) = sqrt_2m * cos(X * W.row(k).t() + b(k));
  }

  arma::mat ZtZ = Z.t() * Z;
  ZtZ.diag() += lambda * n;
  arma::mat R = chol(ZtZ);

  arma::vec tmp = solve(trimatl(R.t()), Z.t() * y);
  arma::vec beta_hat = solve(trimatu(R), tmp);

  return List::create(
    Named("coefficients") = beta_hat,
    Named("Z") = Z,
    Named("m") = m,
    Named("W") = W,
    Named("b") = b,
    Named("m") = m,
    Named("n_threads") = n_threads
  );
}


// [[Rcpp::export]]
arma::mat predict_rff(List model, arma::mat X_new) {
  arma::mat W = as<arma::mat>(model["W"]);
  arma::vec b = as<arma::vec>(model["b"]);
  arma::mat coef_hat = as<arma::mat>(model["coefficients"]);
  int n_new = X_new.n_rows;
  int m = W.n_rows;

  arma::mat Z_new = arma::cos(X_new * W.t() + arma::repmat(b.t(), n_new, 1)) * std::sqrt(2.0/m);

  return Z_new * coef_hat;
}


// [[Rcpp::export]]
arma::mat make_Z(const arma::mat& X,
                 const arma::mat& W,
                 const arma::vec& b){
  const int m = W.n_rows;

  arma::mat Z = X * W.t();
  Z.each_row() += b.t();

  Z = arma::cos(Z);
  Z *= std::sqrt(2.0 / double(m));

  return Z;
}
