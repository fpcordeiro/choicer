// kernel_test_exports.cpp — Thin Rcpp wrappers exposing choicer_internal.h
// helpers for unit testing. These functions are NOT user-facing API; they are
// @noRd and only exported so that tests can pin the raw-array softmax and
// log-sum-exp to their Armadillo counterparts bit for bit.
//
// DO NOT add any of these to the public documentation or NAMESPACE.

// [[Rcpp::depends(RcppArmadillo)]]
#include "choicer.h"
#include "choicer_internal.h"

//' stable_softmax_n() and log_sum_exp_n() next to stable_softmax() and
//' logSumExp(), for a bitwise comparison in tests
//'
//' @param v Numeric vector of utilities (length >= 1).
//' @return List with the shifted utilities, probabilities and log-denominator
//'   from both softmax versions, and both log-sum-exps of v.
//' @noRd
// [[Rcpp::export]]
Rcpp::List test_softmax_n(const arma::vec& v) {
  arma::vec v_arma = v, p_arma;
  const double ld_arma = stable_softmax(v_arma, p_arma);
  arma::vec v_raw = v, p_raw(v.n_elem);
  const int n = static_cast<int>(v.n_elem);
  const double ld_raw = stable_softmax_n(v_raw.memptr(), p_raw.memptr(), n);
  return Rcpp::List::create(
      Rcpp::Named("v_arma") = v_arma, Rcpp::Named("v_raw") = v_raw,
      Rcpp::Named("p_arma") = p_arma, Rcpp::Named("p_raw") = p_raw,
      Rcpp::Named("log_denom_arma") = ld_arma,
      Rcpp::Named("log_denom_raw") = ld_raw,
      Rcpp::Named("lse_arma") = logSumExp(v),
      Rcpp::Named("lse_raw") = log_sum_exp_n(v.memptr(), n));
}
