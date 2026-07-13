// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <functional>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace Rcpp;
using namespace arma;

// ================================================================
// REML criterion, no mean/intercept term (matches exactCV/fastCV,
// which fit KRR directly on y with no centering):
//
//   h(lambda) = n/2*log(y'W^{-1}y/n) + 1/2*log|W|,   W = K + n*lambda*I
//
// h, dh/dt, d2h/dt2 (t = log(lambda)) returned together so Newton-Raphson
// can be run directly on the log scale.
// ================================================================
struct REMLEval {
  double h;
  double dh_dt;
  double d2h_dt;
};

// ----------------------------------------------------------------
// EXACT: K = V D V' eigendecomposition reused across lambda.
// v = V'y.
// ----------------------------------------------------------------
static inline REMLEval eval_exact(const double* ev,
                                   const double* vp,
                                   int n,
                                   double lambda) {
  double nl = (double)n * lambda;

  // quad=y'W^{-1}y, R2=y'W^{-2}y, R3=y'W^{-3}y, T1=tr(W^{-1}), T2=tr(W^{-2}), ldet=log|W|
  double quad = 0, R2 = 0, R3 = 0, T1 = 0, T2 = 0, ldet = 0;
  for (int j = 0; j < n; j++) {
    double d   = ev[j] + nl;
    double id  = 1.0 / d;
    double id2 = id * id;
    double vj2 = vp[j] * vp[j];

    quad += vj2 * id;
    R2   += vj2 * id2;
    R3   += vj2 * id2 * id;
    T1   += id;
    T2   += id2;
    ldet += std::log(d);
  }
  if (quad < 1e-14) quad = 1e-14;

  // dW/dlambda = n*I  =>  d/dlambda[a'W^{-1}b] = -n*a'W^{-2}b,
  //              d2/dlambda2[a'W^{-1}b] = 2n^2*a'W^{-3}b
  double dquad_dl  = -(double)n * R2;
  double d2quad_dl = 2.0 * n * n * R3;

  double dldet_dl  = (double)n * T1;
  double d2ldet_dl = -(double)n * n * T2;

  double h = 0.5*n*std::log(quad/n) + 0.5*ldet;

  double dh_dl = 0.5*n*dquad_dl/quad
               + 0.5*dldet_dl;

  double d2h_dl = 0.5*n*((quad*d2quad_dl - dquad_dl*dquad_dl)/(quad*quad))
                + 0.5*d2ldet_dl;

  double dh_dt  = lambda * dh_dl;
  double d2h_dt = dh_dt + lambda * lambda * d2h_dl;

  return {h, dh_dt, d2h_dt};
}

// ----------------------------------------------------------------
// LOW-RANK: K ~= FF', F'F = Q D Q' eigendecomposition reused across lambda.
// hy = Q'F'y. All bilinear forms a'W^{-1}b (W=FF'+nl*I) via
// Woodbury: a'W^{-1}b = (a'b - sum_j ha_j*hb_j/(mu_j+nl)) / nl.
// ----------------------------------------------------------------
static inline REMLEval eval_lowrank(const double* ev,
                                     const double* hyp,
                                     double yy,
                                     int n,
                                     int m,
                                     double lambda) {
  double nl = (double)n * lambda;

  double P1_yy = 0, P2_yy = 0, P3_yy = 0;
  double R1 = 0, mR2 = 0, ldet = 0;

  for (int j = 0; j < m; j++) {
    double d   = ev[j] + nl;
    double id  = 1.0 / d;
    double id2 = id * id;
    double hy2 = hyp[j] * hyp[j];

    P1_yy += hy2 * id;
    P2_yy += hy2 * id2;
    P3_yy += hy2 * id2 * id;
    R1    += id;
    mR2   += id2;
    ldet  += std::log(d);
  }

  double quad = (yy - P1_yy) / nl;
  if (quad < 1e-14 * yy) quad = 1e-14 * yy;

  // Woodbury-form derivative: d/dlambda[a'W^{-1}b] = (P2(a,b) - a'W^{-1}b)/lambda
  double dquad_dl  = (P2_yy - quad) / lambda;
  double d2quad_dl = -2.0*n*P3_yy/lambda - 2.0*dquad_dl/lambda;

  double log_det   = ldet + (double)(n - m) * std::log(nl);
  double dldet_dl  = (double)n * R1 + (double)(n - m) / lambda;
  double d2ldet_dl = -(double)n * n * mR2 - (double)(n - m) / (lambda*lambda);

  double h = 0.5*n*std::log(quad/n) + 0.5*log_det;

  double dh_dl = 0.5*n*dquad_dl/quad
               + 0.5*dldet_dl;

  double d2h_dl = 0.5*n*((quad*d2quad_dl - dquad_dl*dquad_dl)/(quad*quad))
                + 0.5*d2ldet_dl;

  double dh_dt  = lambda * dh_dl;
  double d2h_dt = dh_dt + lambda * lambda * d2h_dl;

  return {h, dh_dt, d2h_dt};
}

// ================================================================
// Newton-Raphson on t = log(lambda), with backtracking line search.
// If no backtracked step improves h within 10 halvings, the iteration
// stops at the current point rather than accepting an unverified step.
// ================================================================
static double newton_loop(std::function<REMLEval(double)> eval_fn,
                           double t_init,
                           double t_lo,
                           double t_hi,
                           int    max_iter,
                           double tol) {
  double t = std::max(t_lo, std::min(t_hi, t_init));

  for (int iter = 0; iter < max_iter; iter++) {
    REMLEval rv = eval_fn(std::exp(t));

    if (std::abs(rv.dh_dt) < tol) break;

    double step;
    if (rv.d2h_dt > 1e-10)
      step = rv.dh_dt / rv.d2h_dt;
    else
      step = rv.dh_dt * 0.1;  // Hessian non-positive/degenerate: gradient-descent fallback

    step = std::max(-5.0, std::min(5.0, step));

    double t_new = t - step;
    t_new = std::max(t_lo, std::min(t_hi, t_new));

    double h_old = rv.h;
    double alpha = 1.0;
    bool improved = false;
    double t_accept = t;
    for (int ls = 0; ls < 10; ls++) {
      double t_try = t + alpha * (t_new - t);
      t_try = std::max(t_lo, std::min(t_hi, t_try));
      REMLEval rv_try = eval_fn(std::exp(t_try));
      if (rv_try.h < h_old) { t_accept = t_try; improved = true; break; }
      alpha *= 0.5;
    }

    if (!improved) break;  // no verified decrease found: stop rather than guess
    t = t_accept;
  }

  return std::exp(t);
}

// Newton-Raphson REML for exact KRR
double reml_exact(const arma::mat& K,
                   const arma::vec& y,
                   const arma::vec& lambda_vec,
                   int    max_iter  = 50,
                   double tol       = 1e-8,
                   int    n_threads = 1) {
  int n = (int)K.n_rows;

  arma::vec eigval;
  arma::mat eigvec;
  arma::eig_sym(eigval, eigvec, K);
  eigval.clamp(0.0, eigval.max());

  arma::vec v = eigvec.t() * y;

  const double* ev = eigval.memptr();
  const double* vp = v.memptr();

  double t_lo   = std::log(lambda_vec.min());
  double t_hi   = std::log(lambda_vec.max());
  double t_init = 0.5 * (t_lo + t_hi);

  auto eval_fn = [&](double lam) -> REMLEval {
    return eval_exact(ev, vp, n, lam);
  };

  return newton_loop(eval_fn, t_init, t_lo, t_hi, max_iter, tol);
}

// Newton-Raphson REML (no intercept) for low-rank KRR (K ~= FF')
double reml_lowrank(const arma::mat& F,
                     const arma::vec& y,
                     const arma::vec& lambda_vec,
                     int    max_iter  = 50,
                     double tol       = 1e-8,
                     int    n_threads = 1) {
  int n = (int)F.n_rows;
  int m = (int)F.n_cols;

  arma::vec eigval;
  arma::mat eigvec;
  arma::eig_sym(eigval, eigvec, arma::mat(F.t() * F));
  eigval.clamp(0.0, eigval.max());

  arma::vec hy = eigvec.t() * (F.t() * y);

  double yy = arma::dot(y, y);

  const double* ev  = eigval.memptr();
  const double* hyp = hy.memptr();

  double t_lo   = std::log(lambda_vec.min());
  double t_hi   = std::log(lambda_vec.max());
  double t_init = 0.5 * (t_lo + t_hi);

  auto eval_fn = [&](double lam) -> REMLEval {
    return eval_lowrank(ev, hyp, yy, n, m, lam);
  };

  return newton_loop(eval_fn, t_init, t_lo, t_hi, max_iter, tol);
}
