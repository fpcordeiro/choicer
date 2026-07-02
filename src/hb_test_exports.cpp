// hb_test_exports.cpp — Thin Rcpp wrappers exposing hb_internal.h internals
// for unit testing (tests/testthat/test-hb-internal.R). These functions are
// NOT user-facing API; they are @noRd and only exported so the Phase-0
// correctness gates (hand-rolled Cholesky/trisolves vs LAPACK references,
// LSE with the implicit outside term, the sigma_d scale-mixture Gibbs) are
// callable from R.
//
// DO NOT add any of these to the public documentation or NAMESPACE.
// They are kept internal through @noRd roxygen tags. All are rng = false:
// none of them touches R's RNG (the sigma_d2 draw uses the package's own
// Xoshiro256pp streams).

// [[Rcpp::depends(RcppArmadillo)]]
#include "choicer.h"       // brings in RcppArmadillo.h; arma:: available
#include "hb_internal.h"   // included AFTER choicer.h (per contract)

//' Hand-rolled lower Cholesky for testing hb_internal.h
//'
//' @param A Symmetric matrix.
//' @return List with `ok` (FALSE on a non-positive pivot) and `L` (lower
//'   Cholesky factor; zero-filled when `ok` is FALSE).
//' @noRd
// [[Rcpp::export(rng = false)]]
Rcpp::List hb_test_chol(const arma::mat& A) {
  arma::mat L(A.n_rows, A.n_cols, arma::fill::zeros);
  const bool ok = hb_chol_lower(A, L);
  return Rcpp::List::create(Rcpp::Named("ok") = ok, Rcpp::Named("L") = L);
}

//' Hand-rolled triangular solve for testing hb_internal.h
//'
//' Solves L x = b (`transpose = FALSE`, forward substitution) or L' x = b
//' (`transpose = TRUE`, back substitution) for a lower-triangular L.
//'
//' @param L Lower-triangular matrix.
//' @param b Right-hand side vector.
//' @param transpose Solve against L' instead of L.
//' @return Solution vector.
//' @noRd
// [[Rcpp::export(rng = false)]]
arma::vec hb_test_trisolve(const arma::mat& L, const arma::vec& b,
                           const bool transpose) {
  arma::vec x;
  if (transpose) {
    hb_back_solve(L, b, x);
  } else {
    hb_forward_solve(L, b, x);
  }
  return x;
}

//' SPD solve via hand-rolled Cholesky + two trisolves for testing
//'
//' @param A Symmetric positive-definite matrix.
//' @param b Right-hand side vector.
//' @return List with `ok` (FALSE when `A` is not SPD) and `x` (solution;
//'   zero-filled when `ok` is FALSE).
//' @noRd
// [[Rcpp::export(rng = false)]]
Rcpp::List hb_test_spd_solve(const arma::mat& A, const arma::vec& b) {
  arma::vec x(b.n_elem, arma::fill::zeros);
  const bool ok = hb_spd_solve(A, b, x);
  return Rcpp::List::create(Rcpp::Named("ok") = ok, Rcpp::Named("x") = x);
}

//' Fixed-order log-sum-exp with optional implicit outside term for testing
//'
//' @param v Vector of utilities.
//' @param include_outside Add the implicit `exp(0)` outside-option term.
//' @return `log(sum(exp(v)))`, plus 1 inside the sum when
//'   `include_outside = TRUE`.
//' @noRd
// [[Rcpp::export(rng = false)]]
double hb_test_logsumexp(const arma::vec& v, const bool include_outside) {
  return hb_logsumexp(v.memptr(), static_cast<int>(v.n_elem),
                      include_outside);
}

//' Two-block Gibbs for sigma_d2 under the half-Cauchy scale mixture
//'
//' Runs `n_iter` sweeps of draw_sigma_d2_conditional() on fixed residuals
//' `xi`, starting from sigma_d2 = a_d = 1, one Xoshiro stream per iteration
//' (tag 0). Exercises both the half-Cauchy Makalic-Schmidt mixture
//' (`half_cauchy = TRUE`, scale `s_d`) and the plain IG(c0, d0) fallback.
//'
//' @param xi Vector of alternative-level residuals.
//' @param n_iter Number of Gibbs sweeps.
//' @param seed Master seed (coerced to uint64_t).
//' @param half_cauchy Use the half-Cauchy scale mixture (else IG fallback).
//' @param s_d Half-Cauchy scale.
//' @param c0 IG fallback shape.
//' @param d0 IG fallback scale.
//' @return Numeric vector of `n_iter` sigma_d2 draws.
//' @noRd
// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector hb_test_sigma_d2_gibbs(const arma::vec& xi,
                                           const int n_iter,
                                           const double seed,
                                           const bool half_cauchy,
                                           const double s_d,
                                           const double c0,
                                           const double d0) {
  if (n_iter < 1) Rcpp::stop("n_iter must be >= 1.");
  HbSigmaDPrior prior;
  prior.half_cauchy = half_cauchy;
  prior.s_d = s_d;
  prior.c0 = c0;
  prior.d0 = d0;
  double sigma_d2 = 1.0;
  double a_d = 1.0;
  const uint64_t useed = static_cast<uint64_t>(seed);
  Rcpp::NumericVector out(n_iter);
  for (int r = 0; r < n_iter; ++r) {
    Xoshiro256pp rng = make_stream(useed, static_cast<uint64_t>(r), 0);
    draw_sigma_d2_conditional(rng, xi, prior, sigma_d2, a_d);
    out[r] = sigma_d2;
  }
  return out;
}
