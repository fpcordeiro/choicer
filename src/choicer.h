#ifndef CHOICER_HPP
#define CHOICER_HPP

#include <RcppArmadillo.h>
#include <cmath>
#include <numeric>
#ifdef _OPENMP
#include <omp.h>
#endif

// Function Declarations ------------------------------------------------------
Rcpp::IntegerVector compute_prefix_sum(const Rcpp::IntegerVector& M);

// Inline Function Definitions ------------------------------------------------
inline double logSumExp(const arma::vec& x) {
  if (x.n_elem == 0) {
    return -arma::datum::inf;                   // log(sum(empty)) = log(0)
  }

  const double a = x.max();

  if (std::isnan(a))            return arma::datum::nan;  // propagate NaN
  if (a ==  arma::datum::inf)   return  arma::datum::inf; // any +inf -> +inf
  if (a == -arma::datum::inf)   return -arma::datum::inf; // all -inf -> -inf

  const double s = arma::accu(arma::exp(x - a)); // s >= 1 because at least one term is exp(0)=1
  return a + std::log(s);
}

// ---------------------------------------------------------------------------
// Nested-Logit per-individual probability helper
//
// Given one individual's inside utilities V_inside (length m_i), the 0-based
// nest index of each inside alternative (nest_idx0_i, length m_i), the *full*
// lambda vector (length n_nests, singletons fixed to 1), n_nests, and the
// outside-option flag, this fills:
//   P_i          (m_i)      joint choice probability  P_ij = P(j|k) * P_k
//   P_j_given_k  (m_i)      within-nest conditional probability P(j|k)
//   P_k          (n_nests)  marginal nest probability P_k
//   log_I_k      (n_nests)  log inclusive value of each nest (-inf if empty)
//   log_P_i      (m_i)      log joint choice probability (stabilised)
//   log_P_outside (scalar)  log probability of the outside option (-inf if absent)
//
// Uses the same two-level log-sum-exp stabilisation as the likelihood kernel.
// Outside option: V = 0, lambda = 1 -> nest term = 0.
inline void nl_individual_probs(
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

#endif