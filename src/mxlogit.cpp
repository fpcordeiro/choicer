// [[Rcpp::depends(RcppArmadillo)]]
#include "choicer.h"

// Scalar logSumExp for two values to avoid vector allocation
inline double logSumExp2(double a, double b) {
  if (a == -arma::datum::inf)
    return b;
  if (b == -arma::datum::inf)
    return a;
  const double max_val = std::max(a, b);
  return max_val + std::log(std::exp(a - max_val) + std::exp(b - max_val));
}

// Reconstruct lower-triangular choleski factor L from L_params
// [[Rcpp::export]]
arma::mat build_L_mat(const arma::vec &L_params, const int K_w,
                      const bool rc_correlation) {
  arma::mat L = arma::zeros(K_w, K_w);
  int idx = 0;
  if (rc_correlation) { // full (lower-triangular) factor
    for (int i = 0; i < K_w; ++i) {
      for (int j = 0; j <= i; ++j, ++idx) {
        double val = L_params(idx);
        if (i == j) { // diagonal - exp()
          L(i, j) = std::exp(val);
        } else {
          L(i, j) = val; // off-diagonal stays unconstrained
        }
      }
    }
  } else { // diagonal-only (Sigma is diagonal)
    for (int k = 0; k < K_w; ++k) {
      L(k, k) = std::exp(L_params(k));
    }
  }
  return L;
}

//' Reconstruct variance matrix L from L_params
//'
//' @param L_params flattened choleski decomposition version of the random coefficient parameters matrix
//' @param K_w dimension of the random coefficient parameter (symmetric) matrix
//' @param rc_correlation whether random coefficients are correlated
//' @return matrix equal to LL', where L is the choleski decomposition of random coefficient matrix
//' @export
// [[Rcpp::export]]
arma::mat build_var_mat(const arma::vec &L_params, const int K_w,
                        const bool rc_correlation) {
  arma::mat L = build_L_mat(L_params, K_w, rc_correlation);
  // Return the variance matrix
  return L * L.t();
}

//' Log-likelihood and gradient for Mixed Logit
//'
//' Computes the log-likelihood and its gradient for the Mixed Logit model using
//' OpenMP for parallelization. Allows for inclusion of alternative-specific
//' constants, outside option, observation weights, correlated random coefficients.
//'
//' @param theta vector collecting model parameters (beta, mu, L, delta (ASCs))
//' @param X design matrix for covariates with fixed coefficients; sum(M_i) x K_x
//' @param W design matrix for covariates with random coefficients; sum(M_i) x K_w or J x K_w
//' @param alt_idx sum(M) x 1 vector with indices of alternatives within each choice set; 1-based indexing
//' @param choice_idx N x 1 vector with indices of chosen alternatives; 1-based indexing relative to X; 0 is used if include_outside_option=True
//' @param M N x 1 vector with number of alternatives for each individual
//' @param weights N x 1 vector with weights for each observation
//' @param eta_draws Array with choice situation draws; K_w x S x N
//' @param rc_dist K_w x 1 integer vector indicating distribution of random coefficients: 0 = normal, 1 = log-normal
//' @param rc_correlation whether random coefficients should be correlated
//' @param rc_mean whether to estimate means for random coefficients. If so, mean parameters (mu) should be included in theta after beta parameters.
//' @param use_asc whether to use alternative-specific constants. If so, parameters should be included in theta after beta and L (and mu, if applicable).
//' @param include_outside_option whether to include outside option normalized to 0 (if so, the outside option is not included in the data)
//' @return List with loglikelihood and gradient evaluated at input arguments
//' @note For log-normal random coefficients (rc_dist=1) with rc_mean=TRUE,
//'   the distribution is a shifted log-normal: beta_k = exp(mu_k) + exp(L_k * eta),
//'   where exp(mu_k) shifts the location and exp(L_k * eta) ~ LogNormal(0, sigma_k^2).
//'   This differs from the textbook parameterization exp(mu_k + L_k * eta).
//' @export
// [[Rcpp::export]]
Rcpp::List mxl_loglik_gradient_parallel(
    const arma::vec &theta, const arma::mat &X, const arma::mat &W,
    const arma::uvec &alt_idx, const arma::uvec &choice_idx,
    const Rcpp::IntegerVector &M, const arma::vec &weights,
    const arma::cube &eta_draws, const arma::uvec &rc_dist,
    const bool rc_correlation = true, const bool rc_mean = false,
    const bool use_asc = true, const bool include_outside_option = false) {

  // Basic dimensions
  const int N = M.size();
  const int K_x = X.n_cols;
  const int K_w = W.n_cols;
  const int Sdraw = eta_draws.n_cols; // # simulations per individual
  const int n_params = theta.n_elem;
  const int L_size =
      rc_correlation ? (K_w * (K_w + 1)) / 2 : K_w; // Size of L block

  // Check rc_dist input
  if (rc_dist.n_elem != K_w) {
    Rcpp::stop("rc_dist must be a vector of length K_w (%d)", K_w);
  }
  if (eta_draws.n_slices != N) {
    Rcpp::stop("eta_draws 3rd dimension (%d) does not match N (%d)",
               eta_draws.n_slices, N);
  }
  if (eta_draws.n_rows != K_w) {
    Rcpp::stop("eta_draws 1st dimension (%d) does not match K_w (%d)",
               eta_draws.n_rows, K_w);
  }

  // Define parameter block start indices
  const int idx_beta_start = 0;
  const int idx_mu_start = K_x;
  const int idx_L_start = rc_mean ? K_x + K_w : K_x;
  const int idx_delta_start = idx_L_start + L_size;

  // Validate parameter indices
  if (K_x <= 0) {
    Rcpp::stop("K_x must be positive, got %d", K_x);
  }
  if (idx_mu_start > n_params) {
    Rcpp::stop("idx_mu_start (%d) exceeds n_params (%d)", idx_mu_start, n_params);
  }
  if (idx_L_start > n_params) {
    Rcpp::stop("idx_L_start (%d) exceeds n_params (%d)", idx_L_start, n_params);
  }
  if (idx_delta_start > n_params + 1) {
    Rcpp::stop("idx_delta_start (%d) exceeds n_params + 1 (%d)", idx_delta_start, n_params + 1);
  }
  // beta: coefficients for design matrix X
  arma::vec beta = theta.subvec(idx_beta_start, idx_mu_start - 1);

  // mu: means of random coefficients
  arma::vec mu;
  if (rc_mean) {
    if (idx_L_start > n_params)
      Rcpp::stop("Theta vector too short: missing mu parameters.");
    mu = theta.subvec(idx_mu_start, idx_L_start - 1);
  } else {
    mu = arma::zeros(K_w); // Fix means to zero if not estimated
  }

  // L: choleski decomposition of random coefficients matrix
  if (idx_delta_start > n_params)
    Rcpp::stop("Theta vector too short: missing L parameters.");
  arma::vec L_params = theta.subvec(idx_L_start, idx_delta_start - 1);
  arma::mat L = build_L_mat(L_params, K_w, rc_correlation); // K_w x K_w

  // delta (ASC)
  arma::vec delta;
  if (use_asc) {
    const int delta_free_len = n_params - idx_delta_start;
    if (delta_free_len <= 0) {
      Rcpp::stop("Theta vector too short: missing delta parameters.");
    }
    if (include_outside_option) {
      // all inside alternatives are free
      delta = theta.subvec(idx_delta_start, n_params - 1);
    } else {
      // first delta is fixed to 0 -> it's not in theta
      delta = arma::zeros(delta_free_len + 1);
      delta.subvec(1, delta_free_len) =
          theta.subvec(idx_delta_start, n_params - 1);
    }
  } else {
    delta.set_size(0); // empty
  }

  // Convenience objects shared by all threads
  arma::uvec alt_idx0 = alt_idx - 1; // 0-based
  Rcpp::IntegerVector S_prefix = compute_prefix_sum(M);

  // mu transformations for distributions
  arma::vec mu_final = mu;
  arma::vec dmu_final_dmu =
      arma::ones(K_w); // gradient without the exp transform
  if (rc_mean) {
    for (int k = 0; k < K_w; ++k) {
      if (rc_dist(k) == 1) { // 1 == log-normal
        mu_final(k) = std::exp(mu(k));
        dmu_final_dmu(k) = mu_final(k); // d(exp(mu))/dmu = exp(mu)
      }
    }
  }

  // Pre-compute base utility for all individuals (single BLAS call)
  // base_util = X*beta + W*mu_final + delta, computed once for all rows
  arma::vec base_util = X * beta;
  if (static_cast<int>(W.n_rows) == static_cast<int>(X.n_rows)) {
    base_util += W * mu_final;
  } else {
    arma::vec W_mu = W * mu_final;
    base_util += W_mu.elem(alt_idx0);
  }
  if (use_asc) {
    base_util += delta.elem(alt_idx0);
  }

  // Prepare global accumulators
  double global_loglik = 0.0;
  arma::vec global_grad = arma::zeros(n_params);

#ifdef _OPENMP
#pragma omp parallel
#endif
  {
    // Thread-local accumulators
    double local_loglik = 0.0;
    arma::vec local_grad = arma::zeros(n_params);

    // --- Pre-allocate working memory for the thread ---
    arma::vec eta_i_s(K_w);
    arma::vec gamma_i_s(K_w);
    arma::vec gamma_i_s_final(K_w);
    arma::vec dgamma_final_dgamma(K_w);

    arma::vec V_s;          // Re-sized per individual
    arma::vec inside_utils; // Re-sized per individual
    arma::vec P_s;          // Re-sized per individual

    arma::vec g_s(n_params);
    arma::vec w_ap_vec(K_w);

// Loop over individuals in parallel
#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
    for (int i = 0; i < N; ++i) {
      //  Slice data for individual i
      const int m_i = M[i]; // # inside alternatives
      const int num_choices = include_outside_option ? m_i + 1 : m_i;
      const int start_idx = S_prefix[i];
      const int end_idx = start_idx + m_i - 1;
      const double w_i = weights[i];
      const auto X_i = X.rows(start_idx, end_idx); // m_i x K_x
      const auto alt_idx0_i = alt_idx0.subvec(start_idx, end_idx);
      arma::mat W_i;
      if (W.n_rows == X.n_rows)           // row-aligned with X
        W_i = W.rows(start_idx, end_idx); // m_i x K_w
      else                                // global alt-level W
        W_i = W.rows(alt_idx0_i);

      // chosen alternative index
      int chosen_alt = choice_idx[i];
      if (!include_outside_option)
        chosen_alt -= 1; // to 0-based inside-only

      // Pre-computed base utility for this individual (includes X*beta + W*mu_final + delta)
      const arma::vec base_util_i = base_util.subvec(start_idx, end_idx);

      //  Per-draw accumulators
      double log_P_avg = -std::numeric_limits<double>::infinity();
      arma::vec grad_num = arma::zeros(n_params); // Sigma_s P_s * g_s

      // Size variable vectors for this individual
      V_s.set_size(num_choices);
      inside_utils.set_size(m_i);
      // P_s set_size happens implicitly or can be reserved

      // Loop over simulations
      for (int s = 0; s < Sdraw; ++s) {
        // Get eta_i^s and compute gamma_i^s
        eta_i_s = eta_draws.slice(i).col(s); // Length K_w
        gamma_i_s = L * eta_i_s;             // Length K_w

        // gamma transformations for distributions
        gamma_i_s_final = gamma_i_s;
        dgamma_final_dgamma.fill(1.0); // Reset

        for (int k = 0; k < K_w; ++k) {
          if (rc_dist(k) == 1) { // 1 == log-normal
            gamma_i_s_final(k) = std::exp(gamma_i_s(k));
            dgamma_final_dgamma(k) =
                gamma_i_s_final(k); // d(exp(g))/dg = exp(g)
          }
        }

        // Build utility vector V_i_s for individual i and simulation s
        V_s.zeros(); // Reset V_s

        inside_utils = base_util_i + W_i * gamma_i_s_final; // Length m_i
        if (include_outside_option)
          V_s.subvec(1, num_choices - 1) = inside_utils;
        else
          V_s = inside_utils;

        // Probabilities and log-probability for chosen alternative
        V_s -= V_s.max(); // for numerical stability
        const double log_denom = std::log(arma::accu(arma::exp(V_s)));
        P_s = arma::exp(V_s - log_denom);
        double P_choice = P_s(chosen_alt);
        double log_P = std::log(P_choice);

        // log-sum-exp over draws
        if (s == 0) {
          log_P_avg = log_P;
        } else {
          // Optimized scalar version
          log_P_avg = logSumExp2(log_P_avg, log_P);
        }

        // Gradient g_s = d(log P_choice) / d(theta)
        // ------------------------------
        g_s.zeros(); // Reset g_s

        for (int a = 0; a < num_choices; ++a) {
          const double diff = (a == chosen_alt ? 1.0 : 0.0) - P_s[a];
          // beta
          if (include_outside_option) {
            if (a > 0)
              g_s.subvec(idx_beta_start, idx_mu_start - 1) +=
                  diff * X_i.row(a - 1).t();
          } else {
            g_s.subvec(idx_beta_start, idx_mu_start - 1) +=
                diff * X_i.row(a).t();
          }
          // delta
          if (use_asc) {
            if (include_outside_option) {
              if (a > 0) {
                const int id = alt_idx0_i[a - 1];
                g_s[idx_delta_start + id] += diff;
              }
            } else { // alt 0 normalised
              const int id = alt_idx0_i[a];
              if (id > 0)
                g_s[idx_delta_start + (id - 1)] += diff;
            }
          }

          // mu parameters
          w_ap_vec.zeros(); // technically not needed if we overwrite, but safer
          if (include_outside_option) {
            if (a > 0)
              w_ap_vec = W_i.row(a - 1).t();
          } else {
            w_ap_vec = W_i.row(a).t();
          }
          if (rc_mean) {
            for (int p = 0; p < K_w; ++p) {
              // d(Ua) / d(mu_p) = W_ap * dmu_final_dmu(p)
              double dUa_dmu_p = w_ap_vec(p) * dmu_final_dmu(p);
              g_s[idx_mu_start + p] += diff * dUa_dmu_p;
            }
          }

          // L parameters
          if (rc_correlation) {
            // Full L matrix
            int lp_idx = 0;
            for (int p = 0; p < K_w; ++p) {            // loop over rows of L
              for (int q = 0; q <= p; ++q, ++lp_idx) { // loop over cols of L
                // d(L_pq) / d(theta_r)
                double dLpq_dparam = (p == q ? L(p, p) : 1.0);
                // d(gamma_p) / d(theta_r) = [d(L_pq)/d(theta_r)] * eta_q
                double dgamma_p_dtheta_r = dLpq_dparam * eta_i_s(q);
                // d(gamma_p_final) / d(gamma_p)
                double dgamma_p_final_dgamma_p = dgamma_final_dgamma(p);
                // d(Utility_a) / d(theta_r) = W_ap *
                // [d(gamma_p_final)/d(gamma_p)] * [d(gamma_p)/d(theta_r)]
                double dUa_dtheta_r =
                    w_ap_vec(p) * dgamma_p_final_dgamma_p * dgamma_p_dtheta_r;
                // d(logP) / d(theta_r) = [d(logP)/d(U_a)] * [d(U_a)/d(theta_r)]
                g_s[idx_L_start + lp_idx] += diff * dUa_dtheta_r;
              }
            }
          } else {
            // Diagonal L matrix
            for (int p = 0; p < K_w; ++p) {
              double dLpp_dparam = L(p, p);
              double dgamma_p_dtheta_r = dLpp_dparam * eta_i_s(p);
              double dgamma_p_final_dgamma_p = dgamma_final_dgamma(p);
              double dUa_dtheta_r =
                  w_ap_vec(p) * dgamma_p_final_dgamma_p * dgamma_p_dtheta_r;
              g_s[idx_L_start + p] += diff * dUa_dtheta_r;
            }
          }
        } // end alt loop

        // accumulate numerator     Sigma_s P_s * g_s
        grad_num += P_choice * g_s;

      } // end S loop
      // Compute gradient contribution for individual i
      local_grad += w_i * grad_num * std::exp(-log_P_avg);
      // Finish log-probability: divide by S
      log_P_avg -= std::log((double)Sdraw);
      //  Add weighted contribution to thread totals
      local_loglik += w_i * log_P_avg;
    } // end N loop

// Combine thread-local results
#ifdef _OPENMP
#pragma omp critical
#endif
    {
      global_loglik += local_loglik;
      global_grad += local_grad;
    }
  } // end parallel region

  // Return negative log-likelihood and gradient
  return Rcpp::List::create(Rcpp::Named("objective") = -global_loglik,
                            Rcpp::Named("gradient") = -global_grad);
}

//' Numerical Hessian of the log-likelihood via finite differences for mixed logit
//'
//' @param theta vector collecting model parameters (beta, mu, L, delta (ASCs))
//' @param X design matrix for covariates with fixed coefficients; sum(M_i) x K_x
//' @param W design matrix for covariates with random coefficients; sum(M_i) x K_w or J x K_w
//' @param alt_idx sum(M) x 1 vector with indices of alternatives within each choice set; 1-based indexing
//' @param choice_idx N x 1 vector with indices of chosen alternatives; 1-based indexing relative to X; 0 is used if include_outside_option=True
//' @param M N x 1 vector with number of alternatives for each individual
//' @param weights N x 1 vector with weights for each observation
//' @param eta_draws Array with choice situation draws; K_w x S x N
//' @param rc_dist K_w x 1 integer vector indicating distribution of random coefficients: 0 = normal, 1 = log-normal
//' @param rc_correlation whether random coefficients should be correlated
//' @param rc_mean whether to estimate means for random coefficients. If so, mean parameters (mu) should be included in theta after beta parameters.
//' @param use_asc whether to use alternative-specific constants. If so, parameters should be included in theta after beta and L (and mu, if applicable).
//' @param include_outside_option whether to include outside option normalized to 0 (if so, the outside option is not included in the data)
//' @param eps numerical tolerance
//' @return List with loglikelihood and gradient evaluated at input arguments
//' @export
// [[Rcpp::export]]
arma::mat mxl_loglik_numeric_hessian(
    const arma::vec &theta, const arma::mat &X, const arma::mat &W,
    const arma::uvec &alt_idx, const arma::uvec &choice_idx,
    const Rcpp::IntegerVector &M, const arma::vec &weights,
    const arma::cube &eta_draws, const arma::uvec &rc_dist,
    const bool rc_correlation = true, const bool rc_mean = false,
    const bool use_asc = true, const bool include_outside_option = false,
    double eps = 1e-6) {
  int p = theta.n_elem;
  arma::mat Hess(p, p, arma::fill::zeros);

  // For each dimension i, do a +/- eps perturbation
  for (int i = 0; i < p; i++) {
    // Create a copy of theta
    arma::vec theta_pos = theta;
    arma::vec theta_neg = theta;

    double eps_scaled = eps * std::max(std::abs(theta(i)), 1.0);

    theta_pos(i) += eps_scaled;
    theta_neg(i) -= eps_scaled;

    // Evaluate gradient at theta_pos
    Rcpp::List pos_eval = mxl_loglik_gradient_parallel(
        theta_pos, X, W, alt_idx, choice_idx, M, weights, eta_draws, rc_dist,
        rc_correlation, rc_mean, use_asc, include_outside_option);
    arma::vec grad_pos = Rcpp::as<arma::vec>(pos_eval["gradient"]);

    // Evaluate gradient at theta_neg
    Rcpp::List neg_eval = mxl_loglik_gradient_parallel(
        theta_neg, X, W, alt_idx, choice_idx, M, weights, eta_draws, rc_dist,
        rc_correlation, rc_mean, use_asc, include_outside_option);
    arma::vec grad_neg = Rcpp::as<arma::vec>(neg_eval["gradient"]);

    // Check for NaN or Inf in grad_pos and grad_neg
    if (!grad_pos.is_finite()) {
      Rcpp::Rcout << "Warning: NaN or Inf in grad_pos at index " << i
                  << std::endl;
      grad_pos.fill(0.0); // fallback to zero
    }
    if (!grad_neg.is_finite()) {
      Rcpp::Rcout << "Warning: NaN or Inf in grad_neg at index " << i
                  << std::endl;
      grad_neg.fill(0.0); // fallback to zero
    }

    // central difference for gradient
    arma::vec diff_grad_i = (grad_pos - grad_neg) / (2.0 * eps_scaled);

    // Check for NaN or Inf in diff_grad_i
    if (!diff_grad_i.is_finite()) {
      Rcpp::Rcout << "Warning: NaN or Inf in diff_grad_i at index " << i
                  << std::endl;
      diff_grad_i.fill(0.0); // fallback to zero
    }

    // Put diff_grad_i into column i of Hess
    for (int j = 0; j < p; j++) {
      Hess(j, i) = diff_grad_i(j);
    }
  }
  return Hess;
}

// vech(): lower-triangular vectorisation (including the diagonal)
inline arma::vec vech(const arma::mat &M) {
  arma::uword K = M.n_rows;
  arma::vec out(K * (K + 1) / 2);
  arma::uword idx = 0;
  for (arma::uword j = 0; j < K; ++j)
    for (arma::uword i = j; i < K; ++i)
      out(idx++) = M(i, j);
  return out;
}

//' Utility to compute analytical Jacobian of random coefficient matrix transformed by vech (dVech(Sigma) / dTheta)
//'
//' @param L_params flattened choleski decomposition version of the random coefficient parameters matrix
//' @param K_w dimension of the random coefficient parameter (symmetric) matrix
//' @param rc_correlation whether random coefficients are correlated
//' @return Jacobian (dVech(Sigma) / dTheta)
//' @export
// [[Rcpp::export]]
arma::mat jacobian_vech_Sigma(const arma::vec &L_params, const int K_w,
                              const bool rc_correlation = true) {
  // dimensions
  const int L_size = rc_correlation ? K_w * (K_w + 1) / 2 : K_w;

  arma::mat L = build_L_mat(L_params, K_w, rc_correlation);
  arma::mat J(L_size, L_size, arma::fill::zeros);

  // loop over parameters
  arma::mat E(K_w, K_w, arma::fill::zeros); // holds dL / dtheta_m
  std::size_t idx_param = 0;

  if (rc_correlation) {
    for (int i = 0; i < K_w; ++i) {
      for (int j = 0; j <= i; ++j, ++idx_param) {
        // reset E
        E.zeros();
        if (i == j) {        // diagonal: L_ii = exp(z_i)
          E(i, j) = L(i, i); // dL_ii / dz_i = exp(z_i)
        } else {             // off-diagonal parameter
          E(i, j) = 1.0;
        }
        arma::mat dSigma = E * L.t() + L * E.t(); // product rule
        J.col(idx_param) = vech(dSigma);
      }
    }

  } else {
    // diagonal Sigma only (no correlations)
    // Jacobian of Sigma wrt L_params  (diagonal-only case)
    for (int k = 0; k < K_w; ++k) {
      // Sigma_kk = L_kk^2  ,  L_kk = exp(z_k)
      // dSigma_kk/dz_k = 2 * exp(2 z_k) = 2 * L_kk^2
      double deriv = 2.0 * L(k, k) * L(k, k);
      J(k, k) = deriv;
    }
  }
  return J;
}

//' Analytical Hessian of the log-likelihood v2
//'
//' Computes the Hessian of the log-likelihood for the Mixed Logit model using
//' OpenMP for parallelization. Mirrors the parameters of mxl_loglik_gradient_parallel.
//'
//' @param theta vector collecting model parameters (beta, mu, L, delta (ASCs))
//' @param X design matrix for covariates with fixed coefficients; sum(M_i) x K_x
//' @param W design matrix for covariates with random coefficients; sum(M_i) x K_w or J x K_w
//' @param alt_idx sum(M) x 1 vector with indices of alternatives within each choice set; 1-based indexing
//' @param choice_idx N x 1 vector with indices of chosen alternatives; 1-based indexing relative to X; 0 is used if include_outside_option=True
//' @param M N x 1 vector with number of alternatives for each individual
//' @param weights N x 1 vector with weights for each observation
//' @param eta_draws Array with choice situation draws; K_w x S x N
//' @param rc_dist K_w x 1 integer vector indicating distribution of random coefficients: 0 = normal, 1 = log-normal
//' @param rc_correlation whether random coefficients should be correlated
//' @param rc_mean whether to estimate means for random coefficients.
//' @param use_asc whether to use alternative-specific constants.
//' @param include_outside_option whether to include outside option normalized to 0 (if so, the outside option is not included in the data)
//' @return Hessian evaluated at input arguments
//' @note For log-normal random coefficients (rc_dist=1) with rc_mean=TRUE,
//'   the distribution is a shifted log-normal: beta_k = exp(mu_k) + exp(L_k * eta),
//'   where exp(mu_k) shifts the location and exp(L_k * eta) ~ LogNormal(0, sigma_k^2).
//'   This differs from the textbook parameterization exp(mu_k + L_k * eta).
//' @export
// [[Rcpp::export]]
arma::mat mxl_hessian_parallel(
    const arma::vec &theta, const arma::mat &X, const arma::mat &W,
    const arma::uvec &alt_idx, const arma::uvec &choice_idx,
    const Rcpp::IntegerVector &M, const arma::vec &weights,
    const arma::cube &eta_draws, const arma::uvec &rc_dist,
    const bool rc_correlation = true, const bool rc_mean = false,
    const bool use_asc = true, const bool include_outside_option = false) {
  // Basic dimensions
  const int N = M.size();
  const int K_x = X.n_cols;
  const int K_w = W.n_cols;
  const int Sdraw = eta_draws.n_cols;
  const int n_params = theta.n_elem;
  const int L_size = rc_correlation ? (K_w * (K_w + 1)) / 2 : K_w;

  if (rc_dist.n_elem != K_w) {
    Rcpp::stop("rc_dist must be a vector of length K_w (%d)", K_w);
  }

  // === 1. Parameter Parsing ===
  const int idx_beta_start = 0;
  const int idx_mu_start = K_x;
  const int idx_L_start = rc_mean ? K_x + K_w : K_x;
  const int idx_delta_start = idx_L_start + L_size;

  arma::vec beta = theta.subvec(idx_beta_start, idx_mu_start - 1);
  arma::vec mu;
  if (rc_mean) {
    if (idx_L_start > n_params)
      Rcpp::stop("Theta vector too short: missing mu parameters.");
    mu = theta.subvec(idx_mu_start, idx_L_start - 1);
  } else {
    mu = arma::zeros(K_w);
  }

  if (idx_delta_start > n_params)
    Rcpp::stop("Theta vector too short: missing L parameters.");
  arma::vec L_params = theta.subvec(idx_L_start, idx_delta_start - 1);
  arma::mat L = build_L_mat(L_params, K_w, rc_correlation);

  arma::vec delta;
  if (use_asc) {
    const int delta_free_len = n_params - idx_delta_start;
    if (delta_free_len <= 0) {
      Rcpp::stop("Theta vector too short: missing delta parameters.");
    }
    if (include_outside_option) {
      delta = theta.subvec(idx_delta_start, n_params - 1);
    } else {
      delta = arma::zeros(delta_free_len + 1);
      delta.subvec(1, delta_free_len) =
          theta.subvec(idx_delta_start, n_params - 1);
    }
  } else {
    delta.set_size(0);
  }

  arma::uvec alt_idx0 = alt_idx - 1;
  Rcpp::IntegerVector S_prefix = compute_prefix_sum(M);

  // mu transformations for distributions
  arma::vec mu_final = arma::zeros(K_w);
  arma::vec dmu_final_dmu = arma::zeros(K_w);
  arma::vec dmu2_final_dmu2 = arma::zeros(K_w);
  if (rc_mean) {
    mu_final = mu;
    dmu_final_dmu = arma::ones(K_w); // gradient of transform
    for (int k = 0; k < K_w; ++k) {
      if (rc_dist(k) == 1) { // 1 == log-normal
        mu_final(k) = std::exp(mu(k));
        dmu_final_dmu(k) = mu_final(k);   // d(exp(mu))/dmu = exp(mu)
        dmu2_final_dmu2(k) = mu_final(k); // d^2(exp(mu))/dmu^2 = exp(mu)
      }
    }
  }

  // Pre-compute base utility for all individuals (single BLAS call)
  arma::vec base_util_h = X * beta;
  if (static_cast<int>(W.n_rows) == static_cast<int>(X.n_rows)) {
    base_util_h += W * mu_final;
  } else {
    arma::vec W_mu_h = W * mu_final;
    base_util_h += W_mu_h.elem(alt_idx0);
  }
  if (use_asc) {
    base_util_h += delta.elem(alt_idx0);
  }

  // Global accumulator
  arma::mat global_hess = arma::zeros(n_params, n_params);

#ifdef _OPENMP
#pragma omp parallel
#endif
  {
    arma::mat local_hess = arma::zeros(n_params, n_params);

#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
    for (int i = 0; i < N; ++i) {
      // Slice data for individual i
      const int m_i = M[i];
      const int num_choices = include_outside_option ? m_i + 1 : m_i;
      const int start_idx = S_prefix[i];
      const int end_idx = start_idx + m_i - 1;
      const double w_i = weights[i];
      const auto X_i = X.rows(start_idx, end_idx);
      const auto alt_idx0_i = alt_idx0.subvec(start_idx, end_idx);
      arma::mat W_i;
      if (W.n_rows == X.n_rows)
        W_i = W.rows(start_idx, end_idx);
      else
        W_i = W.rows(alt_idx0_i);

      int chosen_alt = choice_idx[i];
      if (!include_outside_option)
        chosen_alt -= 1;

      // Pre-computed base utility for this individual
      const arma::vec base_util_i = base_util_h.subvec(start_idx, end_idx);

      // Per-individual accumulators
      arma::vec P_s_vec = arma::zeros(Sdraw);
      arma::vec grad_numerator = arma::zeros(n_params);
      arma::mat hess_term1_numerator = arma::zeros(n_params, n_params);

      // Loop over simulations (draws) s
      for (int s = 0; s < Sdraw; ++s) {
        const arma::vec eta_i_s = eta_draws.slice(i).col(s);
        const arma::vec gamma_i_s = L * eta_i_s;

        // === 2. Transformations ===
        arma::vec gamma_i_s_final = gamma_i_s;
        arma::vec dgamma_final_dgamma = arma::ones(K_w);
        arma::vec d2gamma_final_dgamma2 = arma::zeros(K_w);

        for (int k = 0; k < K_w; ++k) {
          if (rc_dist(k) == 1) { // 1 == log-normal
            gamma_i_s_final(k) = std::exp(gamma_i_s(k));
            dgamma_final_dgamma(k) = gamma_i_s_final(k);
            d2gamma_final_dgamma2(k) = gamma_i_s_final(k);
          }
        }

        // === 3. Utility ===
        arma::vec V_s = arma::zeros(num_choices);
        arma::vec inside_utils = base_util_i + W_i * gamma_i_s_final;

        if (include_outside_option)
          V_s.subvec(1, num_choices - 1) = inside_utils;
        else
          V_s = inside_utils;

        V_s -= V_s.max();
        arma::vec P_s = arma::exp(V_s - std::log(arma::accu(arma::exp(V_s))));
        double P_choice_s = P_s(chosen_alt);
        P_s_vec(s) = P_choice_s;

        // === 4. Per-Draw Gradient (g_is) and Hessian (H_is) ===
        arma::vec g_is = arma::zeros(n_params);
        arma::vec sum_Pz = arma::zeros(n_params);
        arma::mat sum_Pzz = arma::zeros(n_params, n_params);
        arma::mat sum_diff_H_V = arma::zeros(n_params, n_params);

        for (int a = 0; a < num_choices; ++a) {
          arma::vec z_as = arma::zeros(n_params);
          arma::mat H_V_as = arma::zeros(n_params, n_params);

          if (!include_outside_option || a > 0) {
            int current_a_idx = include_outside_option ? a - 1 : a;
            arma::vec w_ap_vec = W_i.row(current_a_idx).t();

            // --- beta block ---
            z_as.subvec(idx_beta_start, idx_mu_start - 1) =
                X_i.row(current_a_idx).t();

            // --- delta (ASC) block ---
            if (use_asc) {
              const int id = alt_idx0_i[current_a_idx];
              if (include_outside_option) {
                z_as[idx_delta_start + id] = 1.0;
              } else if (id > 0) {
                z_as[idx_delta_start + (id - 1)] = 1.0;
              }
            }

            // --- mu block ---
            // dV/dmu = W * dmu_final/dmu
            if (rc_mean) {
              for (int p = 0; p < K_w; ++p) {
                // z_as: dV/d(mu_p)
                z_as[idx_mu_start + p] = w_ap_vec(p) * dmu_final_dmu(p);
                // H_V_as: d^2V/d(mu_p)^2
                H_V_as(idx_mu_start + p, idx_mu_start + p) =
                    w_ap_vec(p) * dmu2_final_dmu2(p);
              }
            }

            // --- L block ---
            if (rc_correlation) {
              int lp_idx_r = 0;
              for (int p = 0; p < K_w; ++p) {
                double f_p_prime = dgamma_final_dgamma(p);
                double f_p_double_prime = d2gamma_final_dgamma2(p);

                for (int q = 0; q <= p; ++q, ++lp_idx_r) {
                  int r = idx_L_start + lp_idx_r;

                  double dLpq_dparam_r = (p == q) ? L(p, p) : 1.0;
                  double d2Lpq_dparam2_r = (p == q) ? L(p, p) : 0.0;
                  double dgamma_p_dparam_r = dLpq_dparam_r * eta_i_s(q);
                  double d2gamma_p_dparam2_r = d2Lpq_dparam2_r * eta_i_s(q);

                  // z_as
                  z_as[r] = w_ap_vec(p) * f_p_prime * dgamma_p_dparam_r;

                  // H_V_as: L-L Diagonal
                  H_V_as(r, r) =
                      w_ap_vec(p) *
                      (f_p_double_prime * std::pow(dgamma_p_dparam_r, 2) +
                       f_p_prime * d2gamma_p_dparam2_r);

                  // H_V_as: mu-L Cross-term is zero (not included)

                  // H_V_as: L-L Off-Diagonal
                  for (int q_s = 0; q_s < q; ++q_s) {
                    int param_s = idx_L_start + lp_idx_r - (q - q_s);
                    double dLpq_s_dparam = (p == q_s) ? L(p, p) : 1.0;
                    double dgamma_p_dparam_s = dLpq_s_dparam * eta_i_s(q_s);
                    double d2V_dLr_dLs = w_ap_vec(p) * f_p_double_prime *
                                         dgamma_p_dparam_r * dgamma_p_dparam_s;
                    H_V_as(r, param_s) = d2V_dLr_dLs;
                    H_V_as(param_s, r) = d2V_dLr_dLs;
                  }
                }
              }
            } else { // Diagonal L matrix
              for (int p = 0; p < K_w; ++p) {
                int r = idx_L_start + p;

                double f_p_prime = dgamma_final_dgamma(p);
                double f_p_double_prime = d2gamma_final_dgamma2(p);

                double dLpp_dparam = L(p, p);
                double d2Lpp_dparam2 = L(p, p);
                double dgamma_p_dparam = dLpp_dparam * eta_i_s(p);
                double d2gamma_p_dparam2 = d2Lpp_dparam2 * eta_i_s(p);

                // z_as
                z_as[r] = w_ap_vec(p) * f_p_prime * dgamma_p_dparam;
                // H_V_as: L-L Diagonal
                H_V_as(r, r) = w_ap_vec(p) * (f_p_double_prime *
                                                  std::pow(dgamma_p_dparam, 2) +
                                              f_p_prime * d2gamma_p_dparam2);

                // H_V_as: mu-L Cross-term removed (is zero)
              }
            }
          } // end if not outside option

          // === 5. Accumulate Hessian Components ===
          double diff = (a == chosen_alt ? 1.0 : 0.0) - P_s(a);
          g_is += diff * z_as;
          sum_Pz += P_s(a) * z_as;
          sum_Pzz += P_s(a) * z_as * z_as.t();
          sum_diff_H_V += diff * H_V_as;
        } // end alt loop

        // H_is = d^2(log P_is) / d(theta)^2
        arma::mat H_is = -sum_Pzz + sum_Pz * sum_Pz.t() + sum_diff_H_V;

        // Accumulate for individual i
        grad_numerator += P_choice_s * g_is;
        hess_term1_numerator += P_choice_s * (g_is * g_is.t() + H_is);
      } // end S loop

      // === 6. Finalize Hessian for Individual i ===
      double P_i_hat = arma::sum(P_s_vec);
      if (P_i_hat > 1e-12) {
        arma::vec g_i = grad_numerator / P_i_hat;
        arma::mat H_i_term1 = hess_term1_numerator / P_i_hat;
        arma::mat H_i = H_i_term1 - g_i * g_i.t();
        local_hess += w_i * H_i;
      }
    } // end N loop

#ifdef _OPENMP
#pragma omp critical
#endif
    {
      global_hess += local_hess;
    }
  } // end parallel region

  return -global_hess;
}

// ============================================================================
// Mixed Logit: Share Prediction and BLP Contraction
// ============================================================================

// Internal function for computing simulated market shares
arma::vec mxl_predict_shares_internal(
    const arma::mat& X,
    const arma::mat& W,
    const arma::vec& beta,
    const arma::vec& mu_final,         // Transformed mu (exp(mu) for log-normal)
    const arma::mat& L,                // Cholesky factor
    const arma::uvec& alt_idx0,        // 0-based indexing
    const Rcpp::IntegerVector& M,
    const Rcpp::IntegerVector& S_prefix,
    const arma::vec& weights,
    const arma::vec& delta,            // Full J-element delta
    const arma::cube& eta_draws,
    const arma::uvec& rc_dist,
    const int num_alts,
    const bool use_asc,
    const bool include_outside_option
) {
  const int N = M.size();
  const int K_w = W.n_cols;
  const int Sdraw = eta_draws.n_cols;
  const double weight_sum = arma::accu(weights);

  if (weight_sum <= 0) {
    Rcpp::stop("Error: Sum of weights must be positive.");
  }

  // Pre-compute base utility for all individuals (single BLAS call)
  arma::vec base_util_s = X * beta;
  if (static_cast<int>(W.n_rows) == static_cast<int>(X.n_rows)) {
    base_util_s += W * mu_final;
  } else {
    arma::vec W_mu_s = W * mu_final;
    base_util_s += W_mu_s.elem(alt_idx0);
  }
  if (use_asc) {
    base_util_s += delta.elem(alt_idx0);
  }

  // Initialize global accumulator for predicted shares
  arma::vec global_shares = arma::zeros(num_alts);

#ifdef _OPENMP
#pragma omp parallel
#endif
  {
    // Thread-local accumulator
    arma::vec local_shares = arma::zeros(num_alts);

    // Thread-local working vectors
    arma::vec eta_i_s(K_w);
    arma::vec gamma_i_s(K_w);
    arma::vec gamma_i_s_final(K_w);

#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
    for (int i = 0; i < N; ++i) {
      const int m_i = M[i];
      const int num_choices = include_outside_option ? m_i + 1 : m_i;
      const int start_idx = S_prefix[i];
      const int end_idx = start_idx + m_i - 1;
      const double w_i = weights[i];
      const arma::uvec alt_idx0_i = alt_idx0.subvec(start_idx, end_idx);

      arma::mat W_i;
      if (W.n_rows == X.n_rows)
        W_i = W.rows(start_idx, end_idx);
      else
        W_i = W.rows(alt_idx0_i);

      // Pre-computed base utility for this individual
      const arma::vec base_util_i = base_util_s.subvec(start_idx, end_idx);

      // Accumulate probabilities over draws
      arma::vec P_bar_i = arma::zeros(num_choices);

      for (int s = 0; s < Sdraw; ++s) {
        // Get eta and compute gamma
        eta_i_s = eta_draws.slice(i).col(s);
        gamma_i_s = L * eta_i_s;

        // Apply distribution transforms
        gamma_i_s_final = gamma_i_s;
        for (int k = 0; k < K_w; ++k) {
          if (rc_dist(k) == 1) {  // log-normal
            gamma_i_s_final(k) = std::exp(gamma_i_s(k));
          }
        }

        // Build utility vector
        arma::vec V_s = arma::zeros(num_choices);
        arma::vec inside_utils = base_util_i + W_i * gamma_i_s_final;

        if (include_outside_option) {
          V_s.subvec(1, num_choices - 1) = inside_utils;
        } else {
          V_s = inside_utils;
        }

        // Compute probabilities with numerical stability
        V_s -= V_s.max();
        double log_denom = std::log(arma::accu(arma::exp(V_s)));
        arma::vec P_s = arma::exp(V_s - log_denom);

        P_bar_i += P_s;
      }  // end s loop

      P_bar_i /= static_cast<double>(Sdraw);

      // Accumulate shares by alternative
      if (include_outside_option) {
        local_shares(0) += w_i * P_bar_i(0);
      }
      for (int a = 0; a < m_i; ++a) {
        if (include_outside_option) {
          local_shares(alt_idx0_i(a) + 1) += w_i * P_bar_i(a + 1);
        } else {
          local_shares(alt_idx0_i(a)) += w_i * P_bar_i(a);
        }
      }
    }  // end i loop

#ifdef _OPENMP
#pragma omp critical
#endif
    {
      global_shares += local_shares;
    }
  }  // end parallel region

  return global_shares / weight_sum;
}

//' BLP contraction mapping for mixed logit
//'
//' Finds the ASC (delta) parameters such that predicted market shares
//' match target shares, using the contraction mapping of Berry, Levinsohn,
//' and Pakes (1995).
//'
//' @param delta J-1 or J vector with initial guess for deltas (ASCs)
//' @param target_shares J vector with target market shares
//' @param X design matrix for fixed coefficients; sum(M_i) x K_x
//' @param W design matrix for random coefficients; sum(M_i) x K_w or J x K_w
//' @param beta K_x vector with fixed coefficients
//' @param mu K_w vector with mean parameters (raw, will be transformed if log-normal)
//' @param L_params Cholesky parameters vector
//' @param alt_idx sum(M) x 1 vector with indices of alternatives; 1-based indexing
//' @param M N x 1 vector with number of alternatives for each individual
//' @param weights N x 1 vector with weights for each observation
//' @param eta_draws Array with draws; K_w x S x N
//' @param rc_dist K_w vector indicating distribution (0=normal, 1=log-normal)
//' @param rc_correlation whether random coefficients are correlated
//' @param rc_mean whether mu parameters represent means (TRUE) or are zero (FALSE)
//' @param include_outside_option whether outside option is included
//' @param tol convergence tolerance (default 1e-8)
//' @param max_iter maximum iterations (default 1000)
//' @return vector with converged delta (ASC) values
//' @export
// [[Rcpp::export]]
arma::vec mxl_blp_contraction(
    const arma::vec& delta,
    const arma::vec& target_shares,
    const arma::mat& X,
    const arma::mat& W,
    const arma::vec& beta,
    const arma::vec& mu,
    const arma::vec& L_params,
    const arma::uvec& alt_idx,
    const Rcpp::IntegerVector& M,
    const arma::vec& weights,
    const arma::cube& eta_draws,
    const arma::uvec& rc_dist,
    const bool rc_correlation = true,
    const bool rc_mean = false,
    const bool include_outside_option = false,
    const double tol = 1e-8,
    const int max_iter = 1000
) {
  const int K_w = W.n_cols;
  const bool use_asc = true;

  // Build L matrix
  arma::mat L = build_L_mat(L_params, K_w, rc_correlation);

  // Transform mu for log-normal coefficients
  arma::vec mu_final = mu;
  if (rc_mean) {
    for (int k = 0; k < K_w; ++k) {
      if (rc_dist(k) == 1) {  // log-normal
        mu_final(k) = std::exp(mu(k));
      }
    }
  }

  // Convert to 0-based indexing
  arma::uvec alt_idx0 = alt_idx - 1;

  // Deduce number of alternatives from data
  const int J_inside = static_cast<int>(arma::max(alt_idx0)) + 1;  // inside options
  const int num_alts = include_outside_option ? (J_inside + 1) : J_inside;  // total options incl. outside

  // Validate target shares length
  if (static_cast<int>(target_shares.n_elem) != num_alts) {
    Rcpp::stop("Error: target_shares must have length %d (total alternatives, incl. outside if present).", num_alts);
  }

  // Compute prefix sums
  Rcpp::IntegerVector S_prefix = compute_prefix_sum(M);

  // ---------------------------------------------------------------------------
  // Harmonize delta input:
  // - include_outside_option = TRUE: delta must have length J_inside (inside
  //   alternatives only). The outside option has no ASC (utility normalized to 0).
  // - include_outside_option = FALSE: allow J-1 (free) or J (full) length. Pad
  //   a leading zero if only free deltas are provided; keep baseline anchored.
  // ---------------------------------------------------------------------------
  arma::vec delta_current;  // length J_inside
  if (include_outside_option) {
    if (static_cast<int>(delta.n_elem) == J_inside) {
      delta_current = delta;
    } else {
      Rcpp::stop("Error: delta must have length %d (inside alternatives only).", J_inside);
    }
  } else {
    if (static_cast<int>(delta.n_elem) == J_inside) {
      delta_current = delta;
    } else if (static_cast<int>(delta.n_elem) == J_inside - 1) {
      delta_current = arma::zeros(J_inside);
      delta_current.subvec(1, J_inside - 1) = delta;  // pad baseline with 0
    } else {
      Rcpp::stop("Error: delta must have length %d (full) or %d (free, with baseline omitted).",
                 J_inside, J_inside - 1);
    }
    // Anchor baseline at zero for identification
    delta_current -= delta_current(0);
  }

  // Compute initial predicted shares
  arma::vec shares_pred = mxl_predict_shares_internal(
    X, W, beta, mu_final, L, alt_idx0, M, S_prefix, weights,
    delta_current, eta_draws, rc_dist, num_alts, use_asc, include_outside_option
  );

  // Work with inside shares only for the contraction step
  arma::vec shares_pred_inside = include_outside_option
                                 ? shares_pred.subvec(1, num_alts - 1)
                                 : shares_pred;
  arma::vec target_shares_inside = include_outside_option
                                   ? target_shares.subvec(1, num_alts - 1)
                                   : target_shares;

  // Guard against zeros before taking logs
  auto validate_positive = [](const arma::vec& v, const char* name) {
    if (v.min() <= 0.0) {
      Rcpp::stop("%s contains non-positive entries; cannot take logarithm.", name);
    }
  };
  validate_positive(shares_pred_inside, "Predicted shares");
  validate_positive(target_shares_inside, "Target shares");

  arma::vec log_shares_old = arma::log(shares_pred_inside);
  arma::vec log_shares_target = arma::log(target_shares_inside);

  // Iteration
  int iter = 0;
  double residual = 10.0;

  while (iter < max_iter) {
    arma::vec delta_new = delta_current + (log_shares_target - log_shares_old);

    // Re-anchor baseline each iteration when there is no outside option
    if (!include_outside_option) {
      delta_new -= delta_new(0);
    }

    residual = arma::max(arma::abs(delta_new - delta_current));

    if (residual < tol) {
      break;
    }

    delta_current = delta_new;
    shares_pred = mxl_predict_shares_internal(
      X, W, beta, mu_final, L, alt_idx0, M, S_prefix, weights,
      delta_current, eta_draws, rc_dist, num_alts, use_asc, include_outside_option
    );

    shares_pred_inside = include_outside_option
                         ? shares_pred.subvec(1, num_alts - 1)
                         : shares_pred;
    validate_positive(shares_pred_inside, "Predicted shares");
    log_shares_old = arma::log(shares_pred_inside);
    ++iter;
  }

  if (iter >= max_iter) {
    Rcpp::Rcout << "Warning: Maximum iterations reached without convergence. Residual: "
                << residual << std::endl;
  }

  return delta_current;
}

//' Compute aggregate elasticities for mixed logit model
//'
//' Computes the aggregate elasticity matrix (weighted average of individual
//' elasticities) for the Mixed Logit model. The elasticity E(i,j) represents
//' the percentage change in the probability of choosing alternative i when
//' the attribute of alternative j changes by 1%.
//'
//' @param theta parameter vector (beta, \[mu\], L, delta)
//' @param X design matrix for fixed coefficients; sum(M_i) x K_x
//' @param W design matrix for random coefficients; sum(M_i) x K_w or J x K_w
//' @param alt_idx sum(M) x 1 vector with indices of alternatives; 1-based indexing
//' @param choice_idx N x 1 vector (kept for API consistency, not used)
//' @param M N x 1 vector with number of alternatives for each individual
//' @param weights N x 1 vector with weights for each observation
//' @param eta_draws Array with draws; K_w x S x N
//' @param rc_dist K_w vector indicating distribution (0=normal, 1=log-normal)
//' @param elast_var_idx 1-based index of the variable for elasticity computation
//' @param is_random_coef TRUE if variable is in W (random coef), FALSE if in X (fixed coef)
//' @param rc_correlation whether random coefficients are correlated
//' @param rc_mean whether mu parameters are estimated
//' @param use_asc whether ASCs are included
//' @param include_outside_option whether outside option is included
//' @return J x J matrix of aggregate elasticities
//' @export
// [[Rcpp::export]]
arma::mat mxl_elasticities_parallel(
    const arma::vec& theta,
    const arma::mat& X,
    const arma::mat& W,
    const arma::uvec& alt_idx,
    const arma::uvec& choice_idx,
    const Rcpp::IntegerVector& M,
    const arma::vec& weights,
    const arma::cube& eta_draws,
    const arma::uvec& rc_dist,
    const int elast_var_idx,
    const bool is_random_coef,
    const bool rc_correlation = true,
    const bool rc_mean = false,
    const bool use_asc = true,
    const bool include_outside_option = false
) {
  (void)choice_idx;  // unused, kept for API consistency

  // Basic dimensions
  const int N = M.size();
  const int K_x = X.n_cols;
  const int K_w = W.n_cols;
  const int Sdraw = eta_draws.n_cols;
  const int n_params = theta.n_elem;
  const int L_size = rc_correlation ? (K_w * (K_w + 1)) / 2 : K_w;

  // Convert 1-based R index to 0-based C++ index
  const int var_idx = elast_var_idx - 1;

  // Validate variable index
  if (is_random_coef) {
    if (var_idx < 0 || var_idx >= K_w) {
      Rcpp::stop("elast_var_idx (%d) is out of bounds for W matrix (K_w=%d).",
                 elast_var_idx, K_w);
    }
  } else {
    if (var_idx < 0 || var_idx >= K_x) {
      Rcpp::stop("elast_var_idx (%d) is out of bounds for X matrix (K_x=%d).",
                 elast_var_idx, K_x);
    }
  }

  // Parameter block indices
  const int idx_beta_start = 0;
  const int idx_mu_start = K_x;
  const int idx_L_start = rc_mean ? K_x + K_w : K_x;
  const int idx_delta_start = idx_L_start + L_size;

  // Extract beta
  arma::vec beta = theta.subvec(idx_beta_start, idx_mu_start - 1);
  const double beta_k = is_random_coef ? 0.0 : beta(var_idx);

  // Extract and transform mu
  arma::vec mu_final = arma::zeros(K_w);
  if (rc_mean) {
    arma::vec mu = theta.subvec(idx_mu_start, idx_L_start - 1);
    mu_final = mu;
    for (int k = 0; k < K_w; ++k) {
      if (rc_dist(k) == 1) {  // log-normal
        mu_final(k) = std::exp(mu(k));
      }
    }
  }

  // Extract L parameters and build L matrix
  arma::vec L_params = theta.subvec(idx_L_start, idx_delta_start - 1);
  arma::mat L = build_L_mat(L_params, K_w, rc_correlation);

  // Extract delta
  arma::vec delta;
  if (use_asc) {
    const int delta_free_len = n_params - idx_delta_start;
    if (delta_free_len <= 0) {
      Rcpp::stop("Theta vector too short: missing delta parameters.");
    }
    if (include_outside_option) {
      delta = theta.subvec(idx_delta_start, n_params - 1);
    } else {
      delta = arma::zeros(delta_free_len + 1);
      delta.subvec(1, delta_free_len) = theta.subvec(idx_delta_start, n_params - 1);
    }
  } else {
    delta.set_size(0);
  }

  // Convert to 0-based indexing
  arma::uvec alt_idx0 = alt_idx - 1;

  // Determine total number of alternatives
  const int J_inside = use_asc ? static_cast<int>(delta.n_elem) : (arma::max(alt_idx0) + 1);
  const int J_total = include_outside_option ? J_inside + 1 : J_inside;

  // Compute prefix sums
  const Rcpp::IntegerVector S_prefix = compute_prefix_sum(M);

  // Pre-compute base utility for all individuals (single BLAS call)
  arma::vec base_util_e = X * beta;
  if (static_cast<int>(W.n_rows) == static_cast<int>(X.n_rows)) {
    base_util_e += W * mu_final;
  } else {
    arma::vec W_mu_e = W * mu_final;
    base_util_e += W_mu_e.elem(alt_idx0);
  }
  if (use_asc) {
    base_util_e += delta.elem(alt_idx0);
  }

  // Global accumulators
  arma::mat global_elas_matrix = arma::zeros(J_total, J_total);
  double global_total_weight = 0.0;

#ifdef _OPENMP
#pragma omp parallel
#endif
  {
    // Thread-local accumulators
    arma::mat local_elas_matrix = arma::zeros(J_total, J_total);
    double local_total_weight = 0.0;

    // Thread-local working vectors
    arma::vec eta_i_s(K_w);
    arma::vec gamma_i_s(K_w);
    arma::vec gamma_i_s_final(K_w);

#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
    for (int i = 0; i < N; ++i) {
      const int m_i = M[i];
      const int num_choices = include_outside_option ? m_i + 1 : m_i;
      const int start_idx = S_prefix[i];
      const int end_idx = start_idx + m_i - 1;
      const double w_i = weights[i];
      const auto X_i = X.rows(start_idx, end_idx);
      const arma::uvec alt_idx0_i = alt_idx0.subvec(start_idx, end_idx);

      arma::mat W_i;
      if (W.n_rows == X.n_rows)
        W_i = W.rows(start_idx, end_idx);
      else
        W_i = W.rows(alt_idx0_i);

      // Pre-computed base utility for this individual
      const arma::vec base_util_i = base_util_e.subvec(start_idx, end_idx);

      // Map local indices to global alternative indices
      arma::uvec global_j_map(num_choices);
      if (include_outside_option) {
        global_j_map(0) = 0;  // Outside option is global index 0
        global_j_map.subvec(1, m_i) = alt_idx0_i + 1;
      } else {
        global_j_map = alt_idx0_i;
      }

      // Get attribute values for the elasticity variable
      arma::vec x_k_i = arma::zeros(num_choices);
      if (is_random_coef) {
        if (include_outside_option) {
          x_k_i.subvec(1, num_choices - 1) = W_i.col(var_idx);
        } else {
          x_k_i = W_i.col(var_idx);
        }
      } else {
        if (include_outside_option) {
          x_k_i.subvec(1, num_choices - 1) = X_i.col(var_idx);
        } else {
          x_k_i = X_i.col(var_idx);
        }
      }

      // Compute P_bar (average probabilities) and accumulate elasticity terms
      arma::vec P_bar_i = arma::zeros(num_choices);
      arma::mat elas_accum = arma::zeros(num_choices, num_choices);

      for (int s = 0; s < Sdraw; ++s) {
        // Get eta and compute gamma
        eta_i_s = eta_draws.slice(i).col(s);
        gamma_i_s = L * eta_i_s;

        // Apply distribution transforms
        gamma_i_s_final = gamma_i_s;
        for (int k = 0; k < K_w; ++k) {
          if (rc_dist(k) == 1) {  // log-normal
            gamma_i_s_final(k) = std::exp(gamma_i_s(k));
          }
        }

        // Get effective coefficient for this draw
        double beta_k_eff;
        if (is_random_coef) {
          beta_k_eff = mu_final(var_idx) + gamma_i_s_final(var_idx);
        } else {
          beta_k_eff = beta_k;
        }

        // Build utility vector
        arma::vec V_s = arma::zeros(num_choices);
        arma::vec inside_utils = base_util_i + W_i * gamma_i_s_final;

        if (include_outside_option) {
          V_s.subvec(1, num_choices - 1) = inside_utils;
        } else {
          V_s = inside_utils;
        }

        // Compute probabilities
        V_s -= V_s.max();
        double log_denom = std::log(arma::accu(arma::exp(V_s)));
        arma::vec P_s = arma::exp(V_s - log_denom);

        P_bar_i += P_s;

        // Accumulate elasticity terms for this draw
        for (int j_local = 0; j_local < num_choices; ++j_local) {
          const double P_j = P_s(j_local);

          for (int m_local = 0; m_local < num_choices; ++m_local) {
            const double P_m = P_s(m_local);
            const double x_km = x_k_i(m_local);

            double elas_term;
            if (j_local == m_local) {
              // Own-elasticity: beta_k * x_k * P_j * (1 - P_j)
              elas_term = beta_k_eff * x_km * P_j * (1.0 - P_j);
            } else {
              // Cross-elasticity: -beta_k * x_km * P_j * P_m
              elas_term = -beta_k_eff * x_km * P_j * P_m;
            }

            elas_accum(j_local, m_local) += elas_term;
          }
        }
      }  // end s loop

      P_bar_i /= static_cast<double>(Sdraw);
      elas_accum /= static_cast<double>(Sdraw);

      // Compute final elasticities: E = elas_accum / P_bar
      // and map to global indices
      for (int j_local = 0; j_local < num_choices; ++j_local) {
        const int global_j = global_j_map(j_local);
        const double P_bar_j = P_bar_i(j_local);

        if (P_bar_j > 1e-12) {  // Avoid division by zero
          for (int m_local = 0; m_local < num_choices; ++m_local) {
            const int global_m = global_j_map(m_local);
            double elasticity = elas_accum(j_local, m_local) / P_bar_j;
            local_elas_matrix(global_j, global_m) += w_i * elasticity;
          }
        }
      }

      local_total_weight += w_i;
    }  // end i loop

#ifdef _OPENMP
#pragma omp critical
#endif
    {
      global_elas_matrix += local_elas_matrix;
      global_total_weight += local_total_weight;
    }
  }  // end parallel region

  // Compute weighted average
  if (global_total_weight > 1e-10) {
    global_elas_matrix /= global_total_weight;
  }

  return global_elas_matrix;
}
