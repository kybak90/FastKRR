// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List pchol(const arma::mat& A,
           const arma::vec& y,
           double lambda,
           Rcpp::Nullable<int> m = R_NilValue,
           double eps = 1e-6,
           bool verbose = true){

  int n = A.n_rows;
  int r = m.isNull() ? n : Rcpp::as<int>(m);
  if (!m.isNull() && (r <= 0 || r > static_cast<int>(n))) {
    Rcpp::Rcout << "Reset m=n, since argument 'm' must be in the range [1, nrow(A)]\n";
    r = n;
  }



  arma::mat R(n, r, arma::fill::zeros);
  arma::mat B = A;
  arma::mat P = arma::eye<arma::mat>(n, n);

  arma::vec diag_B = B.diag();
  arma::vec schur_diags = arma::sqrt(diag_B);

  // Select pivot for the first column
  arma::uword max_idx = diag_B.index_max();
  double max_diag = diag_B(max_idx);

  if(max_idx != 0){
    B.swap_cols(0, max_idx);
    B.swap_rows(0, max_idx);
    P.swap_cols(0, max_idx);
  }

  R(0, 0) = std::sqrt(max_diag);
  R.col(0).subvec(1, n - 1) = B.col(0).subvec(1, n - 1) / R(0, 0);
  int rank = 1;

  //Select pivot for the i-th column
  for(int i = 1; i < r; ++i){
    for(int j = i; j < n; ++j){

      double sumsq = 0.0;
      for(int k = 0; k < i; ++k) {
        double v = R(j, k);
        sumsq += v * v;
      }

      double tmp = B(j, j) - sumsq;
      schur_diags[j] = std::sqrt(std::max(tmp, 0.0));
    }

    double sqsum = 0.0;
    for (int j = i; j < n; ++j) {
      sqsum += schur_diags[j] * schur_diags[j];
    }

    if (sqsum <= eps) {
      break;
    }

    max_idx = i + schur_diags.subvec(i, n - 1).index_max();
    max_diag = schur_diags[max_idx];

    R(i, i) = max_diag;

    rank += 1;

    if(max_idx != i){
      for (int k = 0; k < i; ++k) {
        std::swap(R(i, k), R(max_idx, k));
      }


      B.swap_cols(i, max_idx);
      B.swap_rows(i, max_idx);
      P.swap_cols(i, max_idx);
    }

    if (i < n) {
      arma::rowvec Ri_prev = R.row(i).cols(0, i - 1);
      for (int j = i + 1; j < n; ++j) {
        arma::rowvec Rj_prev = R.row(j).cols(0, i - 1);
        R(j, i) = (B(j, i) - arma::accu(Ri_prev % Rj_prev)) / R(i, i);
      }
    }
  }

  arma::mat R_final = R.cols(0, rank - 1);
  arma::mat PR = P * R_final;

  // fit coef vector
  arma::mat PRtPR = PR.t() * PR;
  PRtPR.diag() += lambda * n;
  arma::mat L = chol(PRtPR);

  arma::vec tmp = arma::solve(PRtPR, PR.t() * y, solve_opts::likely_sympd);
  arma::vec coef_hat = y / (n * lambda) - PR * tmp / (n * lambda);



  std::string note = "";
  if (!m.isNull() && rank < r) {
    note = "m truncated at " + std::to_string(rank) + " due to early termination";
    if (verbose)
      Rcpp::Rcout << "Note: " << note << std::endl;
  }

  return Rcpp::List::create(
    Rcpp::Named("P") = P,
    Rcpp::Named("R") = R_final,
    Rcpp::Named("PR") = PR,
    Rcpp::Named("m") = rank,
    Rcpp::Named("coefficients") = coef_hat
  );
}

// [[Rcpp::export]]
Rcpp::List pchol_kernel(const arma::mat& A,
                        Rcpp::Nullable<int> m = R_NilValue,
                        double eps = 1e-6,
                        bool verbose = true){

  int n = A.n_rows;
  int r = m.isNull() ? n : Rcpp::as<int>(m);
  if (!m.isNull() && (r <= 0 || r > static_cast<int>(n))) {
    Rcpp::Rcout << "Reset m=n, since argument 'm' must be in the range [1, nrow(A)]\n";
    r = n;
  }



  arma::mat R(n, r, arma::fill::zeros);
  arma::mat B = A;
  arma::mat P = arma::eye<arma::mat>(n, n);

  arma::vec diag_B = B.diag();
  arma::vec schur_diags = arma::sqrt(diag_B);

  // 1st column pivot
  arma::uword max_idx = diag_B.index_max();
  double max_diag = diag_B(max_idx);

  if(max_idx != 0){
    B.swap_cols(0, max_idx);
    B.swap_rows(0, max_idx);
    P.swap_cols(0, max_idx);
  }

  R(0, 0) = std::sqrt(max_diag);
  R.col(0).subvec(1, n - 1) = B.col(0).subvec(1, n - 1) / R(0, 0);
  int rank = 1;

  for(int i = 1; i < r; ++i){
    for(int j = i; j < n; ++j){

      double sumsq = 0.0;
      for(int k = 0; k < i; ++k) {
        double v = R(j, k);
        sumsq += v * v;
      }

      double tmp = B(j, j) - sumsq;
      schur_diags[j] = std::sqrt(std::max(tmp, 0.0));
    }

    double sqsum = 0.0;
    for (int j = i; j < n; ++j) {
      sqsum += schur_diags[j] * schur_diags[j];
    }

    if (sqsum <= eps) {
      //Rcout << "Schur complement trace-norm <= eps at rank = " << rank << "\n";
      break;
    }

    max_idx = i + schur_diags.subvec(i, n - 1).index_max();
    max_diag = schur_diags[max_idx];

    R(i, i) = max_diag;

    rank += 1;

    if(max_idx != 1){
      for (int k = 0; k < i; ++k) {
        std::swap(R(i, k), R(max_idx, k));
      }


      B.swap_cols(i, max_idx);
      B.swap_rows(i, max_idx);
      P.swap_cols(i, max_idx);
    }


    if (i < n) {
      arma::rowvec Ri_prev = R.row(i).cols(0, i - 1);
      for (int j = i + 1; j < n; ++j) {
        arma::rowvec Rj_prev = R.row(j).cols(0, i - 1);
        R(j, i) = (B(j, i) - arma::accu(Ri_prev % Rj_prev)) / R(i, i);
      }
    }
  }

  arma::mat R_final = R.cols(0, rank - 1);
  arma::mat PR = P * R_final;

  // fit coef vector
  arma::mat PRtPR = PR * PR.t();

  std::string note = "";
  if (!m.isNull() && rank < r) {
    note = "m truncated at " + std::to_string(rank) + " due to early termination";
    if (verbose)
      Rcpp::Rcout << "Note: " << note << std::endl;
  }

  return Rcpp::List::create(
    Rcpp::Named("K_approx") = PRtPR,
    Rcpp::Named("rank") = rank,
    Rcpp::Named("eps") = eps
  );
}
