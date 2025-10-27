// [[Rcpp::depends(RcppArmadillo)]]
#include "choicer.h"

//' Log-likelihood and gradient for multinomial logit model
//'
//' @param theta K + J - 1 or K + J vector with model parameters
//' @param X sum(M) x K design matrix with covariates. Stacks M\[i] x K matrices for individual i.
//' @param alt_idx sum(M) x 1 vector with indices of alternatives within each choice set; 1-based indexing
//' @param choice_idx N x 1 vector with indices of chosen alternatives; 1-based indexing relative to X; 0 is used if include_outside_option=True
//' @param M N x 1 vector with number of alternatives for each individual
//' @param weights N x 1 vector with weights for each observation
//' @param use_asc whether to use alternative-specific constants
//' @param include_outside_option whether to include outside option normalized to 0 (if so, the outside option is not included in the data)
//' @return List with loglikelihood and gradient evaluated at input arguments
//' @export
// [[Rcpp::export]]
Rcpp::List mnl_loglik_gradient_parallel(
    const arma::vec& theta,
    const arma::mat& X,
    const arma::uvec& alt_idx,
    const arma::uvec& choice_idx,
    const Rcpp::IntegerVector& M,
    const arma::vec& weights,
    const bool use_asc = true,
    const bool include_outside_option = false
) {
  // Extract beta and delta from theta
  const int N = M.size();
  const int K = X.n_cols;
  const int n_params = theta.n_elem;

  arma::vec beta = theta.subvec(0, K - 1);
  arma::vec delta;
  int delta_length = 0;

  if (use_asc) {
    delta_length = n_params - K;
    if (delta_length <= 0) {
      Rcpp::stop("Error: ASC parameters expected but not provided.");
    }
    if (include_outside_option) {
      // delta covers all J inside alternatives
      delta = theta.subvec(K, n_params - 1);
    } else {
      // delta_1 = 0 fixed
      delta = arma::zeros(delta_length + 1);
      delta.subvec(1, delta_length) = theta.subvec(K, n_params - 1);
    }
  } else {
    delta_length = 0;
    delta = arma::zeros(delta_length);
  }
  // alt_idx is 1-based indexing => shift to 0-based indexing
  arma::uvec alt_idx0 = alt_idx - 1;

  // Compute prefix sums for indexing
  const Rcpp::IntegerVector S = compute_prefix_sum(M);

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

#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
    for (int i = 0; i < N; ++i) {
      const int m_i         = M[i];
      const int num_choices = include_outside_option ? m_i + 1 : m_i;
      const int start_idx   = S[i];
      const int end_idx     = start_idx + m_i - 1;
      const double w_i      = weights[i];
      const auto X_i        = X.rows(start_idx, end_idx); // M[i] x K
      arma::uvec alt_idx0_i = alt_idx0.subvec(start_idx, end_idx); // M[i]

      // Build utility vector V_i
      arma::vec V_i = arma::zeros(num_choices);
      arma::vec inside_utils = X_i * beta;

      // Add delta to inside utilities
      if (use_asc) inside_utils += delta.elem(alt_idx0_i);

      if (include_outside_option) {
        // outside option at index 0
        V_i.subvec(1, num_choices - 1) = inside_utils;
      } else {
        V_i = inside_utils;
      }

      // log-likelihood -------------------------------------------------------
      
      // Vector of choice probabilities
      V_i -= V_i.max(); // for numerical stability
      double log_denom = std::log(arma::accu(arma::exp(V_i)));
      arma::vec P_i = arma::exp(V_i - log_denom);
      
      // Identify chosen alternative
      int chosen_alt = choice_idx[i];
      if (!include_outside_option) {
        chosen_alt -= 1; // shift by 1 for inside-only indexing
      }
      if (chosen_alt < 0 || chosen_alt >= num_choices) {
        Rcpp::stop("Invalid chosen alternative index for individual %d", i);
      }

      // Probability of chosen alternative
      double V_choice = V_i(chosen_alt);
      double log_P_choice = V_choice - log_denom;
      if (!std::isfinite(log_P_choice)) {
        Rcpp::Rcout << "Warning: log_P_choice is not finite at individual " << i << std::endl;
        log_P_choice = -1e10; // fallback to a large negative value
      }

      // Accumulate local weighted log-likelihood
      local_loglik += w_i * log_P_choice;

      // Gradient -------------------------------------------------------------
      for (int a = 0; a < num_choices; ++a) {
        const double P_ia  = P_i[a];
        const double val   = w_i * ( (a == chosen_alt ? 1.0 : 0.0) - P_ia );

        // beta gradient
        if (include_outside_option) {
            if (a > 0) {                                      // skip outside
                local_grad.subvec(0, K - 1) += val * X_i.row(a - 1).t();
            }
        } else {
            local_grad.subvec(0, K - 1) += val * X_i.row(a).t();
        }

        // delta gradient (if any)
        if (use_asc) {
            if (include_outside_option) {
                if (a > 0) {                                    // only inside alt's
                    const int a_id = alt_idx0_i[a - 1];
                    local_grad[K + a_id] += val;
                }
            } else {
                const int a_id = alt_idx0_i[a];                  // 1 ... J
                if (a_id > 0) local_grad[K + (a_id - 1)] += val; // delta_1 is normalised 0
            }
        }
      } // end of alt loop
    } // end of i loop

    // Combine partial accumulators into global accumulators
#ifdef _OPENMP
#pragma omp critical
#endif
    {
      global_loglik += local_loglik;
      global_grad += local_grad;
    }
  } // end parallel region

  // Return negative log-likelihood and gradient
  return Rcpp::List::create(
    Rcpp::Named("objective") = -global_loglik,
    Rcpp::Named("gradient")  = -global_grad
  );
}

//' Numerical Hessian of the log-likelihood via finite differences
//'
//' @param theta K + J - 1 or K + J vector with model parameters
//' @param X sum(M) x K design matrix with covariates. Stacks M\[i] x K matrices for individual i.
//' @param alt_idx sum(M) x 1 vector with indices of alternatives within each choice set; 1-based indexing
//' @param choice_idx N x 1 vector with indices of chosen alternatives; 1-based indexing relative to X; 0 is used if include_outside_option=True
//' @param M N x 1 vector with number of alternatives for each individual
//' @param weights N x 1 vector with weights for each observation
//' @param use_asc whether to use alternative-specific constants
//' @param include_outside_option whether to include outside option normalized to 0 (if so, the outside option is not included in the data)
//' @param eps finite difference step size
//' @return Hessian evaluated at input arguments
//' @export
// [[Rcpp::export]]
arma::mat mnl_loglik_numeric_hessian(
  const arma::vec& theta,
  const arma::mat& X,
  const arma::uvec& alt_idx,
  const arma::uvec& choice_idx,
  const Rcpp::IntegerVector& M,
  const arma::vec& weights,
  bool use_asc = true,
  bool include_outside_option = true,
  double eps = 1e-6
) {
  int p = theta.n_elem;
  arma::mat Hess(p, p, arma::fill::zeros);

  // For each dimension i, we do a +/- eps perturbation
  // Then measure the difference in the gradient
  for (int i = 0; i < p; i++) {
    // Create a copy of theta
    arma::vec theta_pos = theta;
    arma::vec theta_neg = theta;

    double eps_scaled = eps * std::max(std::abs(theta(i)), 1.0);

    theta_pos(i) += eps_scaled;
    theta_neg(i) -= eps_scaled;

    // Evaluate gradient at theta_pos
    Rcpp::List pos_eval = mnl_loglik_gradient_parallel(theta_pos, X, alt_idx, choice_idx, M, weights, use_asc, include_outside_option);
    arma::vec grad_pos  = Rcpp::as<arma::vec>(pos_eval["gradient"]);

    // Evaluate gradient at theta_neg
    Rcpp::List neg_eval = mnl_loglik_gradient_parallel(theta_neg, X, alt_idx, choice_idx, M, weights, use_asc, include_outside_option);
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
    // Because Hess[i, j] = derivative of gradient[j] wrt theta[i]
    // i.e. row = j, col = i
    for (int j = 0; j < p; j++) {
      Hess(j, i) = diff_grad_i(j);
    }
  }
  return Hess;
}

//' Prediction of choice probabilities and utilities based on fitted model
//'
//' @param theta K + J - 1 or K + J vector with model parameters
//' @param X sum(M) x K design matrix with covariates. Stacks M\[i] x K matrices for individual i.
//' @param alt_idx sum(M) x 1 vector with indices of alternatives within each choice set; 1-based indexing
//' @param M N x 1 vector with number of alternatives for each individual
//' @param use_asc whether to use alternative-specific constants
//' @param include_outside_option whether to include outside option normalized to 0 (if so, the outside option is not included in the data)
//' @return List with choice probability and utility for each choice situation evaluated at input arguments
//' @export
// [[Rcpp::export]]
Rcpp::List mnl_predict(
    const arma::vec& theta,
    const arma::mat& X,
    const arma::uvec& alt_idx,
    const Rcpp::IntegerVector& M,
    const bool use_asc = true,
    const bool include_outside_option = false
) {
  // Extract beta and delta from theta
  const int N = M.size();
  const int K = X.n_cols;
  const int n_params = theta.n_elem;

  arma::vec beta = theta.subvec(0, K - 1);
  arma::vec delta;
  int delta_length = 0;

  if (use_asc) {
    delta_length = n_params - K;
    if (delta_length <= 0) {
      Rcpp::stop("Error: ASC parameters expected but not provided.");
    }
    if (include_outside_option) {
      // delta covers all J inside alternatives
      delta = theta.subvec(K, n_params - 1);
    } else {
      // delta_1 = 0 fixed
      delta = arma::zeros(delta_length + 1);
      delta.subvec(1, delta_length) = theta.subvec(K, n_params - 1);
    }
  } else {
    delta_length = 0;
    delta = arma::zeros(delta_length);
  }
  // alt_idx is 1-based indexing => shift to 0-based indexing
  arma::uvec alt_idx0 = alt_idx - 1;

  // Compute prefix sums for indexing
  Rcpp::IntegerVector S = compute_prefix_sum(M);

  // Preallocate predicted values
  arma::vec V_all = arma::zeros(X.n_rows);
  arma::vec P_all = arma::zeros(X.n_rows);

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
  for (int i = 0; i < N; ++i) {
    const int m_i         = M[i];
    const int num_choices = include_outside_option ? m_i + 1 : m_i;
    const int start_idx   = S[i];
    const int end_idx     = start_idx + m_i - 1;
    const auto X_i        = X.rows(start_idx, end_idx); // M[i] x K
    arma::uvec alt_idx0_i = alt_idx0.subvec(start_idx, end_idx); // M[i]

    // Build utility vector V_i
    arma::vec V_i = arma::zeros(num_choices);
    arma::vec inside_utils = X_i * beta;

    // Add delta to inside utilities
    if (use_asc) inside_utils += delta.elem(alt_idx0_i);

    if (include_outside_option) {
      // outside option at index 0
      V_i.subvec(1, num_choices - 1) = inside_utils;
    } else {
      V_i = inside_utils;
    }

    V_all.subvec(start_idx, end_idx) = inside_utils;
    
    // Vector of choice probabilities
    V_i -= V_i.max(); // for numerical stability
    double log_denom = std::log(arma::sum(arma::exp(V_i)));
    const arma::vec P_i = arma::exp(V_i - log_denom);

    if (include_outside_option) {
      P_all.subvec(start_idx, end_idx) = P_i.subvec(1, num_choices - 1);
    } else {
      P_all.subvec(start_idx, end_idx) = P_i;
    }

  } // end of i loop

  // Return predicted choice probabilities and utilities
  return Rcpp::List::create(
    Rcpp::Named("choice_prob") = P_all,
    Rcpp::Named("utility")  = V_all
  );
}

// Prediciton of shares (internal function)
arma::vec mnl_predict_shares_internal(
    const arma::mat& X,
    const arma::vec& beta,
    const arma::uvec& alt_idx0,        // 0-based indexing of alternatives
    const Rcpp::IntegerVector& M,      // N x 1 vector with number of alternatives for each individual
    const Rcpp::IntegerVector& S,      // N x 1 vector with prefix sums for alternative indices
    const arma::vec& weights,          // N x 1 vector with weights for each observation
    const arma::vec& delta,            // J x 1 vector with alternative-specific constants
    const int num_alts,                // total number of distinct alternatives
    const bool use_asc = true,         // whether to use alternative-specific constants
    const bool include_outside_option = false  // whether to include
) {
  const int N = M.size();
  const double denominator = arma::sum(weights);
  if (denominator <= 0) {
    Rcpp::stop("Error: Sum of weights must be positive.");
  }
  // Initialize global accumulator for predicted shares
  arma::vec global_shares = arma::zeros(num_alts);

#ifdef _OPENMP
#pragma omp parallel
#endif
{
  // Thread-local accumulator for predicted shares
  arma::vec local_shares = arma::zeros(num_alts);

  // Loop over individuals in parallel
#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
  for (int i = 0; i < N; ++i) {
    const int m_i         = M[i];
    const int num_choices = include_outside_option ? m_i + 1 : m_i;
    const int start_idx   = S[i];
    const int end_idx     = start_idx + m_i - 1;
    const auto X_i        = X.rows(start_idx, end_idx); // M[i] x K
    arma::uvec alt_idx0_i = alt_idx0.subvec(start_idx, end_idx); // M[i]

    // Build utility vector V_i
    arma::vec V_i = arma::zeros(num_choices);
    arma::vec inside_utils = X_i * beta;

    // Add delta to inside utilities
    if (use_asc) inside_utils += delta.elem(alt_idx0_i);

    if (include_outside_option) {
      // outside option at index 0
      V_i.subvec(1, num_choices - 1) = inside_utils;
    } else {
      V_i = inside_utils;
    }
   
    // Vector of choice probabilities
    V_i -= V_i.max(); // for numerical stability
    double log_denom = std::log(arma::sum(arma::exp(V_i)));
    const arma::vec P_i = arma::exp(V_i - log_denom);

    // Sum probabilities by alternative
    if (include_outside_option) {
      // outside option at index 0
      local_shares(0) += weights[i] * P_i(0); // outside option
    }
    for (int a = 0; a < m_i; ++a) {
      if (include_outside_option) {
        // outside option at index 0
        local_shares(alt_idx0_i(a) + 1) += weights[i] * P_i(a + 1);
      } else {
        local_shares(alt_idx0_i(a)) += weights[i] * P_i(a);
      }
    }

  } // end of i loop

#ifdef _OPENMP
#pragma omp critical
#endif
  {
    // Combine local sums into global sum
    global_shares += local_shares;
  }
}
  // Return predicted choice probabilities and utilities
  return global_shares / denominator;
}

//' Prediction of market shares based on fitted model
//'
//' @param theta K + J - 1 or K + J vector with model parameters
//' @param X sum(M) x K design matrix with covariates. Stacks M\[i] x K matrices for individual i.
//' @param alt_idx sum(M) x 1 vector with indices of alternatives within each choice set; 1-based indexing
//' @param M N x 1 vector with number of alternatives for each individual
//' @param weights N x 1 vector with weights for each observation
//' @param use_asc whether to use alternative-specific constants
//' @param include_outside_option whether to include outside option normalized to 0 (if so, the outside option is not included in the data)
//' @return vector with predicted market shares for each alternative
//' @export
// [[Rcpp::export]]
arma::vec mnl_predict_shares(
    const arma::vec& theta,            // K + J - 1 or K + J vector with model parameters
    const arma::mat& X,                // sum(M) x K design matrix with covariates. M\[i] x K matrix for individual i
    const arma::uvec& alt_idx,         // sum(M) x 1 vector with indices of alternatives within each choice set; 1-based indexing
    const Rcpp::IntegerVector& M,      // N x 1 vector with number of alternatives for each individual
    const arma::vec& weights,          // N x 1 vector with weights for each observation
    const bool use_asc = true,         // whether to use alternative-specific constants
    const bool include_outside_option = false  // whether to include outside option normalized to 0 (if so, the outside option is not included in the data)
) {
  // Extract beta and delta from theta
  const int K = X.n_cols;
  const int n_params = theta.n_elem;

  arma::vec beta = theta.subvec(0, K - 1);
  arma::vec delta;
  int delta_length = 0;

  if (use_asc) {
    delta_length = n_params - K;
    if (delta_length <= 0) {
      Rcpp::stop("Error: ASC parameters expected but not provided.");
    }
    if (include_outside_option) {
      // delta covers all J inside alternatives
      delta = theta.subvec(K, n_params - 1);
    } else {
      // delta_1 = 0 fixed
      delta = arma::zeros(delta_length + 1);
      delta.subvec(1, delta_length) = theta.subvec(K, n_params - 1);
    }
  } else {
    delta_length = 0;
    delta = arma::zeros(delta_length);
  }
  // alt_idx is 1-based indexing => shift to 0-based indexing
  arma::uvec alt_idx0 = alt_idx - 1;

  // Compute prefix sums for indexing
  Rcpp::IntegerVector S = compute_prefix_sum(M);

  // total number of distinct alternatives
  int num_alts = include_outside_option ? (alt_idx.max() + 1) :  alt_idx.max();

  arma::vec global_shares = mnl_predict_shares_internal(
    X, beta, alt_idx0, M, S, weights, delta, num_alts, use_asc, include_outside_option
  );

  return global_shares;
}

//' BLP95 contraction mapping to find delta given target shares
//'
//' @param delta J x 1 vector with initial guess for deltas (ASCs)
//' @param target_shares J x 1 vector with target shares for each alternative
//' @param X sum(M) x K design matrix with covariates. M\[i] x K matrix for individual i
//' @param beta K x 1 vector with model parameters
//' @param alt_idx sum(M) x 1 vector with indices of alternatives within each choice set; 1-based indexing
//' @param M N x 1 vector with number of alternatives for each individual
//' @param weights N x 1 vector with weights for each observation
//' @param include_outside_option whether to include outside option normalized to 0 (if so, the outside option is not included in the data)
//' @param tol convergence tolerance
//' @param max_iter maximum number of iterations
//' @return vector with contraction's delta (ASCs) output
//' @export
// [[Rcpp::export]]
arma::vec blp_contraction(
  const arma::vec& delta,
  const arma::vec& target_shares,
  const arma::mat& X,
  const arma::vec& beta,
  const arma::uvec& alt_idx,
  const Rcpp::IntegerVector& M,
  const arma::vec& weights,
  const bool include_outside_option = false,
  const double tol      = 1e-8,
  const int    max_iter = 1000
) {
  // Extract beta and delta from theta
  const bool use_asc = true;

  // total number of distinct alternatives
  int num_alts = include_outside_option ? (delta.n_elem + 1) :  delta.n_elem;
  if (target_shares.n_elem != num_alts) {
    Rcpp::stop("Error: target_shares must have the same length as the total number of alternatives.");
  }

  // alt_idx is 1-based indexing => shift to 0-based indexing
  arma::uvec alt_idx0 = alt_idx - 1;

  // Compute prefix sums for indexing
  Rcpp::IntegerVector S = compute_prefix_sum(M);

  // Initialize delta_old
  arma::vec delta_old = arma::zeros(num_alts);
  if (include_outside_option) {
    delta_old.subvec(1, num_alts - 1) = delta; // outside option at index 0
  } else {
    delta_old = delta; // no outside option, delta already has J elements
    delta_old -= delta_old[0];
  }

  // Initialize shares and delta_new
  arma::vec log_shares_old = mnl_predict_shares_internal(
    X, beta, alt_idx0, M, S, weights, delta_old, num_alts, use_asc, include_outside_option
  );
  log_shares_old = arma::log(log_shares_old);
  arma::vec log_shares_target = arma::log(target_shares);
  arma::vec delta_new = delta_old;

  // Initialize convergence criteria
  arma::vec delta_diff(num_alts, arma::fill::ones);
  double residual = 10.0;
  int iter = 0;

  // iterate until convergence or until max_iter reached
  while (iter < max_iter) {
      delta_new = delta_old + (log_shares_target - log_shares_old);
      delta_diff = arma::abs(delta_new - delta_old);
      residual = arma::max(delta_diff);
      if (residual < tol) {
          break;
      }
      delta_old = delta_new;
      log_shares_old = mnl_predict_shares_internal(
        X, beta, alt_idx0, M, S, weights, delta_old, num_alts, use_asc, include_outside_option
      );
      log_shares_old = arma::log(log_shares_old);
      ++iter;
  }

  if (iter >= max_iter) {
    Rcpp::Rcout << "Warning: Maximum iterations reached without convergence." << std::endl;
  }

  // Return the new delta vector, normalized by subtracting the first element
  delta_new -= delta_new[0];

  // Return appropriate delta vector based on include_outside_option
  if (include_outside_option) {
      return delta_new.subvec(1, num_alts - 1);
  } else {
      return delta_new;
  }
}

//' Hessian matrix for multinomial logit model
//'
//' @param theta K + J - 1 or K + J vector with model parameters
//' @param X sum(M) x K design matrix with covariates. Stacks M\[i] x K matrices for individual i.
//' @param alt_idx sum(M) x 1 vector with indices of alternatives within each choice set; 1-based indexing
//' @param choice_idx N x 1 vector with indices of chosen alternatives; 1-based indexing relative to X; 0 is used if include_outside_option=True
//' @param M N x 1 vector with number of alternatives for each individual
//' @param weights N x 1 vector with weights for each observation
//' @param use_asc whether to use alternative-specific constants
//' @param include_outside_option whether to include outside option normalized to 0 (if so, the outside option is not included in the data)
//' @return Hessian matrix of the negative log-likelihood
//' @export
// [[Rcpp::export]]
arma::mat mnl_loglik_hessian_parallel(
    const arma::vec& theta,
    const arma::mat& X,
    const arma::uvec& alt_idx,
    const arma::uvec& choice_idx, // Note: choice_idx is not needed for Hessian
    const Rcpp::IntegerVector& M,
    const arma::vec& weights,
    const bool use_asc = true,
    const bool include_outside_option = false
) {
  (void)choice_idx; // silence unused-parameter warning

  // Extract beta and delta from theta
  const int N = M.size();
  const int K = X.n_cols;
  const int n_params = theta.n_elem;

  // Split theta into beta and delta (ASCs) as in your gradient function
  arma::vec beta = theta.subvec(0, K - 1);
  arma::vec delta;
  int delta_length = 0;
  
  if (use_asc) {
    delta_length = n_params - K;
    if (delta_length <= 0) {
      Rcpp::stop("Error: ASC parameters expected but not provided.");
    }
    if (include_outside_option) {
      delta = theta.subvec(K, n_params - 1);
    } else {
      delta = arma::zeros(delta_length + 1);
      delta.subvec(1, delta_length) = theta.subvec(K, n_params - 1);
    }
  } else {
    delta_length = 0;
    delta = arma::zeros(delta_length);
  }

  // alt_idx is 1-based indexing => shift to 0-based indexing
  arma::uvec alt_idx0 = alt_idx - 1;
  
  // Compute prefix sums for indexing each individual's block in X / alt_idx
  const Rcpp::IntegerVector S = compute_prefix_sum(M);
  
  // Prepare global accumulator
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
      const int m_i         = M[i];
      const int num_choices = include_outside_option ? m_i + 1 : m_i;
      const int start_idx   = S[i];
      const int end_idx     = start_idx + m_i - 1;
      const double w_i      = weights[i];

      // M[i] x K view of design rows for this individual (inside alts only)
      const auto X_i        = X.rows(start_idx, end_idx); // M[i] x K
      arma::uvec alt_idx0_i = alt_idx0.subvec(start_idx, end_idx); // inside alt IDs (0-based), length M[i]

      // Build utility vector V_i
      arma::vec V_i = arma::zeros(num_choices);
      arma::vec inside_utils = X_i * beta;
      
      if (use_asc) inside_utils += delta.elem(alt_idx0_i);

      if (include_outside_option) {
        V_i.subvec(1, num_choices - 1) = inside_utils;
      } else {
        V_i = inside_utils;
      }
      
      // Calculate Probabilities ---------------------------------------------
      V_i -= V_i.max(); // for numerical stability
      double log_denom = std::log(arma::sum(arma::exp(V_i)));
      arma::vec P_i = arma::exp(V_i - log_denom);
      
      // Calculate Hessian components for individual i -----------------------
      arma::vec sum_P_Z    = arma::zeros(n_params);
      arma::mat sum_P_Z_Zt = arma::zeros(n_params, n_params);
      arma::vec Z_a        = arma::zeros(n_params);

      for (int a = 0; a < num_choices; ++a) {
        const double P_ia  = P_i[a];
        
        // Reset Z_a to zero
        Z_a.zeros(); 
        
        // Construct Z_a (covariate vector for alternative a)
        if (include_outside_option) {
            if (a > 0) { // skip outside option (a=0)
                // beta part
                Z_a.subvec(0, K - 1) = X_i.row(a - 1).t();
                // delta part
                if (use_asc) {
                    const int a_id = alt_idx0_i[a - 1]; // 0 ... J-1
                    Z_a[K + a_id] = 1.0;
                }
            }
            // if a == 0, Z_a remains all zeros
        } else {
            // beta part
            Z_a.subvec(0, K - 1) = X_i.row(a).t();
            // delta part
            if (use_asc) {
                const int a_id = alt_idx0_i[a]; // 0 ... J-1
                if (a_id > 0) { // delta_0 is normalized to 0
                    Z_a[K + (a_id - 1)] = 1.0;
                }
            }
        }
        
        // Accumulate H_i components
        sum_P_Z    += P_ia * Z_a;
        sum_P_Z_Zt += P_ia * (Z_a * Z_a.t()); // outer product
      } // end of alt loop

      // Compute H_i and add to local accumulator
      // H_i = sum_P_Z_Zt - (sum_P_Z * sum_P_Z.t())
      local_hess += w_i * ( (sum_P_Z * sum_P_Z.t()) - sum_P_Z_Zt );

    } // end of i loop

    // Combine partial accumulators into global accumulator
#ifdef _OPENMP
#pragma omp critical
#endif
    {
      global_hess += local_hess;
    }
  } // end parallel region

  // Enforce exact symmetry (numerical tidy-up)
  global_hess = 0.5 * (global_hess + global_hess.t());

  // Return Hessian of the *negative* log-likelihood
  return -global_hess;
}