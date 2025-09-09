// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
Rcpp::List nystrom(const arma::mat& K, const arma::vec& y,
                   int m, double lambda, int n_threads = 4)
{
  int max_threads = 1;
  #ifdef _OPENMP
    max_threads = omp_get_num_procs();
  #endif
  if (max_threads <= 3)
    n_threads = 1;
  else
    n_threads = std::min(n_threads, max_threads - 1);

  int n = K.n_rows;
  if (m <= 0 || m > n) {
    Rcpp::Rcout << "Reset m = n since argument 'm' must be in the range [1, nrow(K)]\n";
    m = n;
  }

  arma::mat U;
  arma::mat V;
  arma::vec s;

  bool success = arma::svd(U, s, V, K(span(0, m - 1), span(0, m - 1)));
  if(!success){
    Rcpp::stop("Singular value decomposition failed");
  }

  arma::mat R(n, m, arma::fill::zeros);

#ifdef _OPENMP
#pragma omp parallel for num_threads(n_threads)
#endif
  for(int j = 0; j < n; j++){
    for(int i = 0; i < m; i ++){
      R(j, i) = arma::dot(K.row(j).subvec(0, m - 1), U.col(i)) / sqrt(s(i));
    }
  }


  arma::mat RtR = R.t() * R;
  arma::mat identity_m = arma::eye(m, m);


  arma::vec tmp = arma::solve(RtR + (n * lambda) * identity_m, R.t() * y, solve_opts::likely_sympd);
  arma::vec coef_hat = y / (n * lambda) - R * tmp / (n * lambda);


  return Rcpp::List::create(
    Rcpp::Named("R") = R,
    Rcpp::Named("m") = m,
    Rcpp::Named("coefficients") = coef_hat,
    Rcpp::Named("n_threads") = n_threads
  );
}



// [[Rcpp::export]]
Rcpp::List nystrom_kernel(const arma::mat& K,
                          Rcpp::Nullable<int> m_in = R_NilValue,
                          int n_threads = 4)
{
  int max_threads = 1;
  #ifdef _OPENMP
    max_threads = omp_get_num_procs();
  #endif
  if (max_threads <= 3)
    n_threads = 1;
  else
    n_threads = std::min(n_threads, max_threads - 1);


  int n = K.n_rows;
  int m;
  if (m_in.isNull()) {
    m = std::max(1, n / 10); // 최소 1 이상
  } else {
    m = Rcpp::as<int>(m_in);
  }

  if (m <= 0 || m > n) {
    Rcpp::Rcout << "Reset m = n since argument 'm' must be in the range [1, nrow(K)]\n";
    m = n;
  }

  arma::mat U;
  arma::mat V;
  arma::vec s;

  bool success = arma::svd(U, s, V, K(span(0, m - 1), span(0, m - 1)));
  if(!success){
    Rcpp::stop("Singular value decomposition failed");
  }

  arma::mat R(n, m, arma::fill::zeros);
#ifdef _OPENMP
#pragma omp parallel for num_threads(n_threads)
#endif
  for(int j = 0; j < n; j++){
    for(int i = 0; i < m; i ++){
      R(j, i) = arma::dot(K.row(j).subvec(0, m - 1), U.col(i)) / sqrt(s(i));
    }
  }

  arma::mat RtR = R * R.t();

  return Rcpp::List::create(
    Rcpp::Named("K_approx") = RtR,
    Rcpp::Named("m") = m
  );
}

// [[Rcpp::export]]
Rcpp::List nystrom_vec(const arma::mat& K, const arma::vec& y,
                       int m, double lambda)
{
  int n = K.n_rows;

  if (m <= 0 || m > n) {
    Rcpp::Rcout << "Reset m = n since argument 'm' must be in the range [1, nrow(K)]\n";
    m = n;
  }

  arma::mat U;
  arma::mat V;
  arma::vec s;

  bool success = arma::svd(U, s, V, K(span(0, m - 1), span(0, m - 1)));
  if(!success){
    Rcpp::stop("Singular value decomposition failed");
  }

  arma::mat S = arma::repmat(arma::sqrt(s.t()), n, 1);
  arma::mat R = K.cols(0, m - 1) * U / S;

  arma::mat RtR = R.t() * R;
  arma::mat identity_m = arma::eye(m, m);

  arma::vec tmp = arma::solve(RtR + (n * lambda) * identity_m, R.t() * y, solve_opts::likely_sympd);
  arma::vec coef_hat = y / (n * lambda) - R * tmp / (n * lambda);


  return Rcpp::List::create(
    Rcpp::Named("R") = R,
    Rcpp::Named("m") = m,
    Rcpp::Named("coefficients") = coef_hat
  );
}
