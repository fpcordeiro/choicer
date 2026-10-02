// kernel_test_exports.cpp — Thin Rcpp wrappers exposing choicer_internal.h
// helpers for unit testing. These functions are NOT user-facing API; they are
// @noRd and only exported so that tests can pin the raw-array softmax,
// log-sum-exp and nested logit probabilities to their Armadillo counterparts
// bit for bit.
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

// nl_individual_probs() as it was written on Armadillo expressions (src/choicer.h
// before the raw-array rewrite), kept verbatim as the reference the rewrite
// must reproduce bit for bit.
static void nl_individual_probs_ref(
    const arma::vec& V_inside,
    const arma::uvec& nest_idx0_i,
    const arma::vec& lambda,
    const int n_nests,
    const bool include_outside_option,
    arma::vec& P_i,
    arma::vec& P_j_given_k,
    arma::vec& P_k,
    arma::vec& log_I_k,
    arma::vec& log_P_i,
    double& log_P_outside
) {
  const int m_i = V_inside.n_elem;

  // V_ij / lambda_k  (lambda_k = 1 for singletons)
  arma::vec V_over_lambda = V_inside / lambda.elem(nest_idx0_i);

  // --- log_I_k (inclusive value) via log-sum-exp within each nest ---
  arma::vec max_V_k = arma::vec(n_nests).fill(-arma::datum::inf);
  for (int j = 0; j < m_i; ++j) {
    const int k = nest_idx0_i[j];
    if (V_over_lambda[j] > max_V_k[k]) {
      max_V_k[k] = V_over_lambda[j];
    }
  }

  arma::vec I_k_unscaled = arma::zeros(n_nests);
  for (int j = 0; j < m_i; ++j) {
    const int k = nest_idx0_i[j];
    if (std::isfinite(max_V_k[k])) {
      I_k_unscaled[k] += std::exp(V_over_lambda[j] - max_V_k[k]);
    }
  }

  log_I_k = arma::vec(n_nests).fill(-arma::datum::inf);
  for (int k = 0; k < n_nests; ++k) {
    if (I_k_unscaled[k] > 0) {
      log_I_k[k] = max_V_k[k] + std::log(I_k_unscaled[k]);
    }
  }

  // --- log(P_k) (nest probability) ---
  arma::vec nest_terms = lambda % log_I_k;

  double max_nest_term = nest_terms.max();
  if (!std::isfinite(max_nest_term)) {
    max_nest_term = 0;
  }

  double sum_exp_nest_terms = arma::accu(arma::exp(nest_terms - max_nest_term));
  if (include_outside_option) {
    // Outside option: V=0, lambda=1 -> term = 0
    sum_exp_nest_terms += std::exp(0.0 - max_nest_term);
  }

  const double log_denom_P_nest = max_nest_term + std::log(sum_exp_nest_terms);

  arma::vec log_P_k = nest_terms - log_denom_P_nest;
  log_P_outside = include_outside_option ? (0.0 - log_denom_P_nest)
                                         : -arma::datum::inf;

  // --- log(P_j|k) and log(P_ij) ---
  arma::vec log_P_j_given_k = V_over_lambda - log_I_k.elem(nest_idx0_i);
  log_P_i = log_P_j_given_k + log_P_k.elem(nest_idx0_i);

  P_i = arma::exp(log_P_i);
  P_j_given_k = arma::exp(log_P_j_given_k);
  P_k = arma::exp(log_P_k);
}

//' nl_individual_probs() next to its Armadillo-expression reference, for a
//' bitwise comparison in tests
//'
//' @param V Inside utilities of one individual.
//' @param nest0 0-based nest of each inside alternative (same length as V).
//' @param lambda Full lambda vector, one entry per nest.
//' @param include_outside_option Whether the outside option is in the set.
//' @return List with each output of both versions (`*_ref` the reference).
//' @noRd
// [[Rcpp::export]]
Rcpp::List test_nl_individual_probs(const arma::vec& V,
                                    const Rcpp::IntegerVector& nest0,
                                    const arma::vec& lambda,
                                    const bool include_outside_option) {
  const int m = V.n_elem;
  const int n_nests = lambda.n_elem;
  if (nest0.size() != m) Rcpp::stop("nest0 must have one entry per utility.");
  for (int j = 0; j < m; ++j) {
    if (nest0[j] < 0 || nest0[j] >= n_nests) Rcpp::stop("nest0 out of range.");
  }
  arma::uvec nest_u(m);
  for (int j = 0; j < m; ++j) nest_u[j] = static_cast<arma::uword>(nest0[j]);
  arma::vec P_i, P_j_given_k, P_k, log_I_k, log_P_i;
  double log_P_outside;
  nl_individual_probs_ref(V, nest_u, lambda, n_nests, include_outside_option,
                          P_i, P_j_given_k, P_k, log_I_k, log_P_i,
                          log_P_outside);
  // Buffers larger than this individual needs and filled with R's NA (a NaN
  // with a payload no arithmetic produces), as a thread's are after earlier,
  // larger individuals: an element read or left before this call writes it
  // would show.
  NlProbs pr(m + 3, n_nests);
  for (std::vector<double>* b :
       {&pr.P_i, &pr.P_j_given_k, &pr.log_P_i, &pr.P_k, &pr.log_I_k,
        &pr.V_over_lambda, &pr.log_P_j_given_k, &pr.max_V_k,
        &pr.I_k_unscaled, &pr.nest_terms, &pr.log_P_k}) {
    std::fill(b->begin(), b->end(), NA_REAL);
  }
  pr.log_P_outside = NA_REAL;
  nl_individual_probs(V.memptr(), nest0.begin(), m, lambda, n_nests,
                      include_outside_option, pr);
  auto num = [](const std::vector<double>& x, const int n) {
    return Rcpp::NumericVector(x.begin(), x.begin() + n);
  };
  return Rcpp::List::create(
      Rcpp::Named("P_i_ref") = Rcpp::NumericVector(P_i.begin(), P_i.end()),
      Rcpp::Named("P_i") = num(pr.P_i, m),
      Rcpp::Named("P_j_given_k_ref") =
          Rcpp::NumericVector(P_j_given_k.begin(), P_j_given_k.end()),
      Rcpp::Named("P_j_given_k") = num(pr.P_j_given_k, m),
      Rcpp::Named("P_k_ref") = Rcpp::NumericVector(P_k.begin(), P_k.end()),
      Rcpp::Named("P_k") = num(pr.P_k, n_nests),
      Rcpp::Named("log_I_k_ref") =
          Rcpp::NumericVector(log_I_k.begin(), log_I_k.end()),
      Rcpp::Named("log_I_k") = num(pr.log_I_k, n_nests),
      Rcpp::Named("log_P_i_ref") =
          Rcpp::NumericVector(log_P_i.begin(), log_P_i.end()),
      Rcpp::Named("log_P_i") = num(pr.log_P_i, m),
      Rcpp::Named("log_P_outside_ref") = log_P_outside,
      Rcpp::Named("log_P_outside") = pr.log_P_outside);
}
