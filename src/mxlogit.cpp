// [[Rcpp::depends(RcppArmadillo)]]
#include "choicer.hpp"

// Reconstruct lower-triangular choleski factor L from L_params
// [[Rcpp::export]]
arma::mat build_L_mat(const arma::vec& L_params, const int K_w, const bool rc_correlation) {
  arma::mat L = arma::zeros(K_w, K_w);
  int idx = 0;
  if (rc_correlation) {                   // full (lower‑triangular) factor
    for (int i = 0; i < K_w; ++i) {
      for (int j = 0; j <= i; ++j, ++idx) {
        double val = L_params(idx);
        if (i == j) {                     // diagonal → exp()
          L(i, j) = std::exp(val);
        } else {
          L(i, j) = val;                 // off‑diagonal stays unconstrained
        }
      }
    }
  } else {                                // diagonal‑only (Σ is diagonal)
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
arma::mat build_var_mat(const arma::vec& L_params, const int K_w, const bool rc_correlation) {
  arma::mat L = build_L_mat(L_params, K_w, rc_correlation);
  // Return the variance matrix
  return L * L.t();
}

//' Log-likelihood and gradient for Mixed Logit
//'
//' @param theta vector collecting model parameters (beta, L, delta (ASCs))
//' @param X design matrix for covariates with fixed coefficients; sum(M_i) × K_x
//' @param W design matrix for covariates with random coefficients; sum(M_i) × K_w or J x K_w
//' @param alt_idx sum(M) x 1 vector with indices of alternatives within each choice set; 1-based indexing
//' @param choice_idx N x 1 vector with indices of chosen alternatives; 1-based indexing relative to X; 0 is used if include_outside_option=True
//' @param M N x 1 vector with number of alternatives for each individual
//' @param weights N x 1 vector with weights for each observation
//' @param eta_draws Array with choice situation draws; K_w × S × N 
//' @param rc_correlation whether random coefficients should be correlated
//' @param use_asc whether to use alternative-specific constants
//' @param include_outside_option whether to include outside option normalized to 0 (if so, the outside option is not included in the data)
//' @return List with loglikelihood and gradient evaluated at input arguments
//' @export
// [[Rcpp::export]]
Rcpp::List mxl_loglik_gradient_parallel(
    const arma::vec&            theta,
    const arma::mat&            X,
    const arma::mat&            W,
    const arma::uvec&           alt_idx,
    const arma::uvec&           choice_idx,
    const Rcpp::IntegerVector&  M,
    const arma::vec&            weights,
    const arma::cube&           eta_draws,
    const bool                  rc_correlation  = true,
    const bool                  use_asc         = true,
    const bool                  include_outside_option = false
) {
  // Basic dimensions
  const int N        = M.size();
  const int K_x      = X.n_cols;
  const int K_w      = W.n_cols;
  const int Sdraw    = eta_draws.n_cols;        // # simulations per individual
  const int n_params = theta.n_elem;
  const int L_size   = rc_correlation ? (K_w * (K_w + 1)) / 2 : K_w; // Size of L block

  // beta: coefficients for design matrix X
  arma::vec beta = theta.subvec(0, K_x - 1);

  // L: choleski decomposition of random coefficients matrix
  arma::vec L_params = theta.subvec(K_x, K_x + L_size - 1);
  arma::mat L        = build_L_mat(L_params, K_w, rc_correlation);   // K_w × K_w

  // delta (ASC)
  arma::vec delta;
  if (use_asc) {
    const int delta_free_len = n_params - K_x - L_size;
    if (delta_free_len <= 0) {
      Rcpp::stop("Theta vector too short: missing delta parameters.");
    }
    if (include_outside_option) {
      // all inside alternatives are free
      delta = theta.subvec(K_x + L_size, n_params - 1);
    } else {
      // first delta is fixed to 0 ⇒ it's not in θ
      delta = arma::zeros(delta_free_len + 1);
      delta.subvec(1, delta_free_len) = theta.subvec(K_x + L_size, n_params - 1);
    }
  } else {
    delta.set_size(0);   // empty
  }

  // Convenience objects shared by all threads
  arma::uvec  alt_idx0 = alt_idx - 1; // 0-based
  Rcpp::IntegerVector S_prefix = compute_prefix_sum(M);

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

    // Loop over individuals in parallel
    #ifdef _OPENMP
    #pragma omp for schedule(dynamic)
    #endif
      for (int i = 0; i < N; ++i) {
        //  Slice data for individual i
        const int m_i               = M[i];                       // # inside alternatives
        const int num_choices       = include_outside_option ? m_i + 1 : m_i;
        const int start_idx         = S_prefix[i];
        const int end_idx           = start_idx + m_i - 1;
        const double w_i            = weights[i];
        const auto X_i              = X.rows(start_idx, end_idx);         // m_i × K_x
        const auto alt_idx0_i       = alt_idx0.subvec(start_idx, end_idx);
        arma::mat W_i;
        if (W.n_rows == X.n_rows)                           // row-aligned with X
          W_i = W.rows(start_idx, end_idx);                 // m_i × K_w
        else                                                // global alt-level W
          W_i = W.rows(alt_idx0_i);

        // chosen alternative index
        int chosen_alt = choice_idx[i];
        if (!include_outside_option) chosen_alt -= 1;      // to 0-based inside-only
        if (chosen_alt < 0 || chosen_alt >= num_choices)
          Rcpp::stop("Invalid choice index for individual %d", i);
        
        //  Per-draw accumulators
        double    log_P_avg = -std::numeric_limits<double>::infinity();
        arma::vec grad_num  = arma::zeros(n_params);           // Σ_s P_s * g_s

        // Loop over simulations
        for (int s = 0; s < Sdraw; ++s) {
          // Get eta_i^s and compute gamma_i^s
          const arma::vec eta_i_s = eta_draws.slice(i).col(s); // Size K_w
          const arma::vec gamma_i_s = L * eta_i_s;             // Size K_w

          // Build utility vector V_i_s for individual i and simulation s
          arma::vec V_s = arma::zeros(num_choices);
          arma::vec inside_utils = X_i * beta + W_i * gamma_i_s; // m_i
          if (use_asc)  inside_utils += delta.elem(alt_idx0_i);
          if (include_outside_option)
            V_s.subvec(1, num_choices - 1) = inside_utils;
          else
            V_s = inside_utils;

          // Probabilities
          V_s -= V_s.max();
          const double log_denom = std::log(arma::sum(arma::exp(V_s)));
          arma::vec P_s          = arma::exp(V_s - log_denom);
          double P_choice        = P_s(chosen_alt);
          double log_P           = std::log(P_choice);

          // log-sum-exp over draws
          if (s == 0) {
            log_P_avg = log_P;
          } else {
            log_P_avg = logSumExp( arma::vec({log_P_avg, log_P}) );
          }
          
          // Gradient g_s = ∂ log P_choice / ∂ θ ------------------------------
          arma::vec g_s(n_params, arma::fill::zeros);

          for (int a = 0; a < num_choices; ++a) {
            const double diff = (a == chosen_alt ? 1.0 : 0.0) - P_s[a];
            // β
            if (include_outside_option) {
              if (a > 0)
                g_s.subvec(0, K_x - 1) += diff * X_i.row(a - 1).t();
            } else {
              g_s.subvec(0, K_x - 1) += diff * X_i.row(a).t();
            }
            // delta
            if (use_asc) {
              if (include_outside_option) {
                if (a > 0) {
                  const int id = alt_idx0_i[a - 1];
                  g_s[K_x + L_size + id] += diff;
                }
              } else {                                        // alt 0 normalised
                const int id = alt_idx0_i[a];
                if (id > 0)
                  g_s[K_x + L_size + (id - 1)] += diff;
              }
            }
            // L parameters
            int lp_idx = 0;
            for (int p = 0; p < K_w; ++p) {
              for (int q = 0; q <= p; ++q, ++lp_idx) {
                double dLpq_dparam = (p == q ? L(p,p) : 1.0);   // ∂L/∂θ_diagonal=exp(val)
                double dgamma_p = dLpq_dparam * eta_i_s(q);
                double w_ap = 0.0;
                if (include_outside_option) {
                  if (a > 0)  w_ap = W_i(a - 1, p);
                } else {
                  w_ap = W_i(a, p);
                }
                g_s[K_x + lp_idx] += diff * w_ap * dgamma_p;
              }
            }
          } // end alt loop
          // accumulate numerator     Σ_s P_s * g_s
          grad_num += P_choice * g_s;
        } // end S loop
        // Compute gradient contribution for individual i
        local_grad += w_i * grad_num * std::exp(-log_P_avg);
        // Finish log-probability: divide by S
        log_P_avg -= std::log( (double) Sdraw );
        //  Add weighted contribution to thread totals
        local_loglik += w_i * log_P_avg;
      } // end N loop

    // Combine thread-local results
    #ifdef _OPENMP
    #pragma omp critical
    #endif
    {
      global_loglik += local_loglik;
      global_grad   += local_grad;
    }
  } // end parallel region

  // Return negative log-likelihood and gradient
  return Rcpp::List::create(
    Rcpp::Named("objective") = -global_loglik,
    Rcpp::Named("gradient")  = -global_grad
  );
}

//' Numerical Hessian of the log-likelihood via finite differences for mixed logit
//'
//' @param theta vector collecting model parameters (beta, L, delta (ASCs))
//' @param X design matrix for covariates with fixed coefficients; sum(M_i) × K_x
//' @param W design matrix for covariates with random coefficients; sum(M_i) × K_w or J x K_w
//' @param alt_idx sum(M) x 1 vector with indices of alternatives within each choice set; 1-based indexing
//' @param choice_idx N x 1 vector with indices of chosen alternatives; 1-based indexing relative to X; 0 is used if include_outside_option=True
//' @param M N x 1 vector with number of alternatives for each individual
//' @param weights N x 1 vector with weights for each observation
//' @param eta_draws Array with choice situation draws; K_w × S × N 
//' @param rc_correlation whether random coefficients should be correlated
//' @param use_asc whether to use alternative-specific constants
//' @param include_outside_option whether to include outside option normalized to 0 (if so, the outside option is not included in the data)
//' @param eps numerical tolerance
//' @return List with loglikelihood and gradient evaluated at input arguments
//' @export
// [[Rcpp::export]]
arma::mat mxl_loglik_numeric_hessian(
  const arma::vec&            theta,
  const arma::mat&            X,
  const arma::mat&            W,
  const arma::uvec&           alt_idx,
  const arma::uvec&           choice_idx,
  const Rcpp::IntegerVector&  M,
  const arma::vec&            weights,
  const arma::cube&           eta_draws,
  const bool                  rc_correlation  = true,
  const bool                  use_asc         = true,
  const bool                  include_outside_option = false,
  double eps = 1e-6
) {
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
    Rcpp::List pos_eval = mxl_loglik_gradient_parallel(theta_pos, X, W, alt_idx, choice_idx, M, weights, eta_draws, rc_correlation, use_asc, include_outside_option);
    arma::vec grad_pos  = Rcpp::as<arma::vec>(pos_eval["gradient"]);

    // Evaluate gradient at theta_neg
    Rcpp::List neg_eval = mxl_loglik_gradient_parallel(theta_neg, X, W, alt_idx, choice_idx, M, weights, eta_draws, rc_correlation, use_asc, include_outside_option);
    arma::vec grad_neg  = Rcpp::as<arma::vec>(neg_eval["gradient"]);

    // Check for NaN or Inf in grad_pos and grad_neg
    if (!grad_pos.is_finite()) {
      Rcpp::Rcout << "Warning: NaN or Inf in grad_pos at index " << i << std::endl;
      grad_pos.fill(0.0); // fallback to zero
    }
    if (!grad_neg.is_finite()) {
      Rcpp::Rcout << "Warning: NaN or Inf in grad_neg at index " << i << std::endl;
      grad_neg.fill(0.0); // fallback to zero
    }

    // central difference for gradient
    arma::vec diff_grad_i = (grad_pos - grad_neg) / (2.0 * eps_scaled);

    // Check for NaN or Inf in diff_grad_i
    if (!diff_grad_i.is_finite()) {
      Rcpp::Rcout << "Warning: NaN or Inf in diff_grad_i at index " << i << std::endl;
      diff_grad_i.fill(0.0); // fallback to zero
    }

    // Put diff_grad_i into column i of Hess
    for (int j = 0; j < p; j++) {
      Hess(j, i) = diff_grad_i(j);
    }
  }
  return Hess;
}

// vech(): lower‑triangular vectorisation (including the diagonal)
inline arma::vec vech(const arma::mat& M)
{
  arma::uword K = M.n_rows;
  arma::vec out(K * (K + 1) / 2);
  arma::uword idx = 0;
  for (arma::uword j = 0; j < K; ++j)
    for (arma::uword i = j; i < K; ++i)
      out(idx++) = M(i, j);
  return out;
}

//' Utility to compute analytical Jacobian of random coefficient matrix transformed by vech (d[vech(Σ)] / dθ)
//'
//' @param L_params flattened choleski decomposition version of the random coefficient parameters matrix
//' @param K_w dimension of the random coefficient parameter (symmetric) matrix
//' @param rc_correlation whether random coefficients are correlated
//' @return Jacobian (d[vech(Σ)] / dθ)
//' @export
// [[Rcpp::export]]
arma::mat jacobian_vech_Sigma(
  const arma::vec& L_params,
  const int K_w,
  const bool rc_correlation = true
) {
  // dimensions
  const int L_size  = rc_correlation ? K_w * (K_w + 1) / 2 : K_w;

  arma::mat L = build_L_mat(L_params, K_w, rc_correlation);
  arma::mat J(L_size, L_size, arma::fill::zeros);

  // loop over parameters
  arma::mat E(K_w, K_w, arma::fill::zeros);   // holds dL / dθ_m
  std::size_t idx_param = 0;

  if (rc_correlation) {
    for (int i = 0; i < K_w; ++i) {
      for (int j = 0; j <= i; ++j, ++idx_param) {
        // reset E
        E.zeros();
        if (i == j) {                      // diagonal: L_ii = exp(z_i)
          E(i, j) = L(i, i);               // dL_ii / dz_i = exp(z_i)
        } else {                           // off‑diagonal parameter
          E(i, j) = 1.0;
        }
        arma::mat dSigma = E * L.t() + L * E.t();   // product rule
        J.col(idx_param) = vech(dSigma);
      }
    }

  } else {
    // diagonal Σ only (no correlations)
    // Jacobian of Σ wrt L_params  (diagonal-only case)
    for (int k = 0; k < K_w; ++k) {
        // Σ_kk = L_kk^2  ,  L_kk = exp(z_k)
        // dΣ_kk/dz_k = 2 * exp(2 z_k) = 2 * L_kk^2
        double deriv = 2.0 * L(k,k) * L(k,k);
        J(k, k) = deriv;
    }
  }
  return J;
}

//' Analytical Hessian of the log-likelihood
//'
//' @param theta vector collecting model parameters (beta, L, delta (ASCs))
//' @param X design matrix for covariates with fixed coefficients; sum(M_i) × K_x
//' @param W design matrix for covariates with random coefficients; sum(M_i) × K_w or J x K_w
//' @param alt_idx sum(M) x 1 vector with indices of alternatives within each choice set; 1-based indexing
//' @param choice_idx N x 1 vector with indices of chosen alternatives; 1-based indexing relative to X; 0 is used if include_outside_option=True
//' @param M N x 1 vector with number of alternatives for each individual
//' @param weights N x 1 vector with weights for each observation
//' @param eta_draws Array with choice situation draws; K_w × S × N 
//' @param rc_correlation whether random coefficients should be correlated
//' @param use_asc whether to use alternative-specific constants
//' @param include_outside_option whether to include outside option normalized to 0 (if so, the outside option is not included in the data)
//' @return Hessian evaluated at input arguments
//' @export
// [[Rcpp::export]]
arma::mat mxl_hessian_parallel(
    const arma::vec&            theta,
    const arma::mat&            X,
    const arma::mat&            W,
    const arma::uvec&           alt_idx,
    const arma::uvec&           choice_idx,
    const Rcpp::IntegerVector&  M,
    const arma::vec&            weights,
    const arma::cube&           eta_draws,
    const bool                  rc_correlation = true,
    const bool                  use_asc = true,
    const bool                  include_outside_option = false
) {
  const int N = M.size();
  const int K_x = X.n_cols;
  const int K_w = W.n_cols;
  const int Sdraw = eta_draws.n_cols;
  const int n_params = theta.n_elem;
  const int L_size = rc_correlation ? (K_w * (K_w + 1)) / 2 : K_w;

  arma::vec beta = theta.subvec(0, K_x - 1);
  arma::vec L_params = theta.subvec(K_x, K_x + L_size - 1);
  arma::mat L = build_L_mat(L_params, K_w, rc_correlation);

  arma::vec delta;
  int delta_free_len = 0;
  if (use_asc) {
    delta_free_len = n_params - K_x - L_size;
    if (delta_free_len < 0) Rcpp::stop("Theta vector too short for ASCs.");
    if (include_outside_option) {
      delta = theta.subvec(K_x + L_size, n_params - 1);
    } else {
      delta = arma::zeros(delta_free_len + 1);
      delta.subvec(1, delta_free_len) = theta.subvec(K_x + L_size, n_params - 1);
    }
  }

  arma::uvec alt_idx0 = alt_idx - 1;
  Rcpp::IntegerVector S_prefix = compute_prefix_sum(M);

  // Global accumulator for the Hessian
  arma::mat global_hess = arma::zeros(n_params, n_params);

  #ifdef _OPENMP
  #pragma omp parallel
  #endif
  {
    // Thread-local accumulator
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
      if (W.n_rows == X.n_rows) W_i = W.rows(start_idx, end_idx);
      else W_i = W.rows(alt_idx0_i);

      int chosen_alt = choice_idx[i];
      if (!include_outside_option) chosen_alt -= 1;

      // Per-individual accumulators
      arma::vec P_s_vec = arma::zeros(Sdraw);
      arma::vec grad_numerator = arma::zeros(n_params);
      arma::mat hess_term1_numerator = arma::zeros(n_params, n_params);

      // Loop over simulations (draws) s
      for (int s = 0; s < Sdraw; ++s) {
        const arma::vec eta_i_s = eta_draws.slice(i).col(s);
        const arma::vec gamma_i_s = L * eta_i_s;

        arma::vec V_s = arma::zeros(num_choices);
        arma::vec inside_utils = X_i * beta + W_i * gamma_i_s;
        if (use_asc) inside_utils += delta.elem(alt_idx0_i);
        if (include_outside_option) V_s.subvec(1, num_choices - 1) = inside_utils;
        else V_s = inside_utils;

        V_s -= V_s.max();
        arma::vec P_s = arma::exp(V_s) / arma::sum(arma::exp(V_s));
        double P_choice_s = P_s(chosen_alt);
        P_s_vec(s) = P_choice_s;

        // Calculate g_is and H_is
        arma::vec g_is = arma::zeros(n_params);
        arma::vec sum_Pz = arma::zeros(n_params);
        arma::mat sum_Pzz = arma::zeros(n_params, n_params);
        arma::mat sum_diff_H_V = arma::zeros(n_params, n_params);

        for (int a = 0; a < num_choices; ++a) {
          arma::vec z_as = arma::zeros(n_params); // dV_as/d(theta)
          arma::mat H_V_as = arma::zeros(n_params, n_params); // d^2V_as/d(theta)^2

          int current_a_idx = include_outside_option ? a - 1 : a;
          if (!include_outside_option || a > 0) { // If not the outside option
            // beta part
            z_as.subvec(0, K_x - 1) = X_i.row(current_a_idx).t();
            // delta part
            if (use_asc) {
              const int id = alt_idx0_i[current_a_idx];
              if (include_outside_option) { // delta_0 fixed at 0
                 z_as[K_x + L_size + id] = 1.0;
              } else if (id > 0) { // delta_0 fixed at 0
                 z_as[K_x + L_size + id - 1] = 1.0;
              }
            }
            // L_params part
            const arma::rowvec W_ia = W_i.row(current_a_idx);
            int lp_idx = 0;
            if (rc_correlation) {
              for (int p = 0; p < K_w; ++p) {
                for (int q = 0; q <= p; ++q, ++lp_idx) {
                  double dLpq_dparam = (p == q) ? L(p, p) : 1.0;
                  z_as(K_x + lp_idx) = W_ia(p) * dLpq_dparam * eta_i_s(q);
                  if (p == q) { // Second derivative part (only for diagonal exp() terms)
                    H_V_as(K_x + lp_idx, K_x + lp_idx) = z_as(K_x + lp_idx);
                  }
                }
              }
            } else {
              for (int k=0; k < K_w; ++k) {
                z_as(K_x + k) = W_ia(k) * L(k, k) * eta_i_s(k);
                H_V_as(K_x + k, K_x + k) = z_as(K_x + k);
              }
            }
          }

          double diff = (a == chosen_alt ? 1.0 : 0.0) - P_s(a);
          g_is += diff * z_as;
          sum_Pz += P_s(a) * z_as;
          sum_Pzz += P_s(a) * z_as * z_as.t();
          sum_diff_H_V += diff * H_V_as;
        }
        arma::mat H_is = -sum_Pzz + sum_Pz * sum_Pz.t() + sum_diff_H_V;

        // Accumulate for individual i
        grad_numerator += P_choice_s * g_is;
        hess_term1_numerator += P_choice_s * (g_is * g_is.t() + H_is);
      } // end S loop

      // Finalize Hessian for individual i
      double P_i_hat = arma::sum(P_s_vec); // Denominator: sum over s of P_is
      if (P_i_hat > 1e-12) {
        arma::vec g_i = grad_numerator / P_i_hat;
        arma::mat H_i_term1 = hess_term1_numerator / P_i_hat;
        arma::mat H_i = H_i_term1 - g_i * g_i.t();
        local_hess += w_i * H_i;
      }
    } // end N loop

    // Combine thread-local results
    #ifdef _OPENMP
    #pragma omp critical
    #endif
    {
      global_hess += local_hess;
    }
  } // end parallel region

  return -global_hess; // Return Hessian of negative log-likelihood
}