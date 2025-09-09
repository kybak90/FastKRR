// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::vec solve_chol(const arma::mat& A,
                     const arma::vec& b)
{
  arma::mat L = arma::chol(A, "lower");

  arma::vec tmp = arma::solve(arma::trimatl(L), b); //Lz = b

  return arma::solve(arma::trimatu(L.t()), tmp); // arma L^Tx = z
}


// [[Rcpp::export]]
arma::vec SOLVE_sympd(const arma::mat& A,
                      const arma::vec& b){
  return arma::solve(A, b, solve_opts::likely_sympd); // LL
}
