// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::mat make_kernel(const arma::mat& X,
                      Nullable<NumericMatrix> X_new = R_NilValue,
                      std::string kernel = "gaussian",
                      double rho = 0,
                      int n_threads = 4){

  int max_threads = 1;
  #ifdef _OPENMP
    max_threads = omp_get_num_procs();
  #endif
  if (max_threads <= 3)
    n_threads = 1;
  else
    n_threads = std::min(n_threads, max_threads - 1);


  if (X.n_rows == 0 || X.n_cols == 0) {
    stop("Input matrix X is empty");
  }

  if(rho < 0){
    stop("rho must be positive");
  }else if(rho == 0){
    rho = 1.0;
  }


  arma::mat U;
  if(X_new.isNull()){
    U = X;
  }else{
    U = as<arma::mat>(X_new);
  }


  int n = X.n_rows;
  int p = U.n_rows;
  arma::mat K(p, n, arma::fill::zeros);


  if(kernel == "gaussian")
  {
    if(X_new.isNull()){ // symmetry matrix
      arma::vec sq_norm = arma::sum(arma::square(X), 1);
      arma::mat dist_mat = arma::repmat(sq_norm, 1, X.n_rows) +
        arma::repmat(sq_norm.t(), X.n_rows, 1) -
        2.0 * (X * X.t());
      K = arma::exp(- rho * dist_mat);

    }else{ // for prediction
      arma::vec sq_norm1 = arma::sum(arma::square(U), 1);
      arma::vec sq_norm2 = arma::sum(arma::square(X), 1);
      arma::mat dist_mat = arma::repmat(sq_norm1, 1, n) +
        arma::repmat(sq_norm2.t(), p, 1) -
        2.0 * (U * X.t());
      K = arma::exp(- rho * dist_mat);
    }

  }else if(kernel == "laplace")
  {
    if(X_new.isNull()){ // symmetry matrix

#ifdef _OPENMP
#pragma omp parallel for num_threads(n_threads) schedule(static)
#endif
      for(int i = 0; i < n; ++i){
        arma::rowvec xi = X.row(i);
        for(int j = i; j < n; ++j){
          double dist = arma::accu(arma::abs(xi - X.row(j)));
          double val = std::exp(- rho * dist);

          K(i, j) = val;
          K(j, i) = val;
        }
      }
    }else{ // for prediction
#ifdef _OPENMP
#pragma omp parallel for num_threads(n_threads) schedule(static)
#endif
      for(int i = 0; i < n; ++i){
        arma::rowvec xi = X.row(i);
        for(int j = 0; j < p; ++j){
          double dist = arma::accu(arma::abs(xi - U.row(j)));
          double val = std::exp(- rho * dist);

          K(j, i) = val;
        }
      }
    }
  }else
  {
    Rcpp::stop("Unsupported kernel: Please choose either 'gaussian' or 'laplace'");
  }

  return K;
}
