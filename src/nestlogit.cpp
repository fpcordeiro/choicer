// [[Rcpp::depends(RcppArmadillo)]]
#include "choicer.h"

//' Log-likelihood and gradient for Nested Logit model
//'
//' @param theta (K + n_non_singleton_nests + n_delta) vector with model parameters.
//'        Order: [beta (K), lambda (n_non_singleton_nests), delta (n_delta)]
//' @param X sum(M) x K design matrix with covariates.
//' @param alt_idx sum(M) x 1 vector with indices of alternatives; 1-based indexing.
//' @param choice_idx N x 1 vector with indices of chosen alternatives; 0 for outside option,
//'        1-based index relative to rows in X_i otherwise.
//' @param nest_idx J x 1 vector with indices of nests for each alternative; 1-based indexing (1 to n_nests).
//' @param M N x 1 vector with number of alternatives for each individual.
//' @param weights N x 1 vector with weights for each observation.
//' @param use_asc whether to use alternative-specific constants.
//' @param include_outside_option whether to include outside option normalized to V=0, lambda=1.
//' @return List with loglikelihood and gradient evaluated at input arguments
//' @export
// [[Rcpp::export]]
Rcpp::List nl_loglik_gradient_parallel(
    const arma::vec& theta,
    const arma::mat& X,
    const arma::uvec& alt_idx,
    const arma::uvec& choice_idx,
    const arma::uvec& nest_idx,
    const Rcpp::IntegerVector& M,
    const arma::vec& weights,
    const bool use_asc = true,
    const bool include_outside_option = false
) {
  // Extract dimensions
  const int N = M.size();
  const int K = X.n_cols;
  const int n_params = theta.n_elem;
  const int n_nests = arma::max(nest_idx); // assuming nest_idx uses 1-based indexing

  // --- 1. Parameter Parsing ---
  arma::vec beta = theta.subvec(0, K - 1);
  int delta_length = 0;
  
  // Identify singleton nests (nests with only 1 alternative)
  arma::uvec nest_counts = arma::zeros<arma::uvec>(n_nests);
  for (unsigned int j = 0; j < nest_idx.n_elem; ++j) {
      if (nest_idx[j] > 0 && nest_idx[j] <= n_nests) {
          nest_counts[nest_idx[j] - 1]++; // nest_idx is 1-based
      } else {
          Rcpp::stop("Invalid nest index found in nest_idx.");
      }
  }
  arma::uvec is_singleton = (nest_counts == 1);
  const int n_non_singleton_nests = arma::accu(is_singleton == 0);

  // Lambda parameters (inclusive value coefficients)
  const int lambda_start_idx = K;
  
  // Create a *full* lambda vector (size n_nests), initialized to 1 (ingletons will keep this value)
  arma::vec lambda = arma::ones(n_nests); 
  
  // Create a map to link full nest index 'k' to its gradient position in 'theta'
  // -1 indicates a singleton nest (lambda=1, not a parameter)
  arma::ivec nest_k_to_theta_idx = arma::ivec(n_nests).fill(-1);
  
  if (n_non_singleton_nests > 0) {
    // Check if theta is long enough
    if (theta.n_elem < K + n_non_singleton_nests) {
        Rcpp::stop("Error: theta vector is too short for K + n_non_singleton_nests.");
    }
    // Extract only the non-singleton lambdas from theta
    arma::vec non_singleton_lambdas = theta.subvec(lambda_start_idx, lambda_start_idx + n_non_singleton_nests - 1);
    if (arma::any(non_singleton_lambdas <= 0)) {
        Rcpp::stop("Error: All non-singleton lambda (nest) parameters must be > 0.");
    }
    // "Scatter" non-singleton lambdas into the full lambda vector; build the k -> theta_index map
    int current_lambda_idx = 0;
    for (int k = 0; k < n_nests; ++k) {
        if (is_singleton[k] == 0) { // if NOT a singleton
            lambda[k] = non_singleton_lambdas[current_lambda_idx];
            nest_k_to_theta_idx[k] = lambda_start_idx + current_lambda_idx;
            current_lambda_idx++;
        }
    }
  } else {
    Rcpp::stop("Error: No non-singleton nests found. At least one nest must have multiple alternatives.");
  }
  
  // ASC parameters (delta)
  const int delta_start_idx = K + n_non_singleton_nests;
  arma::vec delta;
  
  if (use_asc) {
    // The length of delta is also calculated from the new start index
    delta_length = n_params - K - n_non_singleton_nests;
    
    if (delta_length < 0) {
      Rcpp::stop("Error: Not enough parameters for K + n_non_singleton_nests + n_delta.");
    }
    
    if (delta_length > 0) {
        if (include_outside_option) {
          delta = theta.subvec(delta_start_idx, delta_start_idx + delta_length - 1);
        } else {
          delta = arma::zeros(delta_length + 1);
          delta.subvec(1, delta_length) = theta.subvec(delta_start_idx, delta_start_idx + delta_length - 1);
        }
    } else {
        // Check against K + n_non_singleton_nests
        if (n_params != K + n_non_singleton_nests) {
             Rcpp::stop("Error: Mismatch in parameters. Expected K + n_non_singleton_nests.");
        }
        delta = arma::zeros((include_outside_option) ? 0 : 1);
    }
  } else {
    delta_length = 0;
    delta = arma::zeros(delta_length);
  }

  // 0-based indexing for inputs
  arma::uvec alt_idx0 = alt_idx - 1;
  arma::uvec nest_idx0 = nest_idx - 1;
  
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
    
    // Make the index map thread-private
    const arma::ivec thread_nest_k_to_theta_idx = nest_k_to_theta_idx;

#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
    for (int i = 0; i < N; ++i) {
      const int m_i         = M[i];
      const int start_idx   = S[i];
      const int end_idx     = start_idx + m_i - 1;
      const double w_i      = weights[i];
      
      // Get individual-specific data
      const auto X_i        = X.rows(start_idx, end_idx); // m_i x K
      arma::uvec alt_idx0_i = alt_idx0.subvec(start_idx, end_idx); // m_i
      arma::uvec nest_idx0_i= nest_idx0.elem(alt_idx0_i); // m_i

      // --- 2. Calculate Utilities and Probabilities ---
      
      // V_ij = X_ij * beta + delta_j
      arma::vec V_inside = X_i * beta;
      if (use_asc) V_inside += delta.elem(alt_idx0_i);

      // V_ij / lambda_k
      // This now correctly uses lambda_k=1 for singleton nests
      arma::vec V_over_lambda = V_inside / lambda.elem(nest_idx0_i);

      // --- 2a. Calculate log_I_k (Inclusive Value) ---
      // log_I_k = log( sum_{j in B_k} exp(V_ij / lambda_k) )
      // Use log-sum-exp trick for numerical stability
      arma::vec max_V_k = arma::vec(n_nests).fill(-arma::datum::inf);
      for(int j = 0; j < m_i; ++j) {
        const int k = nest_idx0_i[j];
        if (V_over_lambda[j] > max_V_k[k]) {
          max_V_k[k] = V_over_lambda[j];
        }
      }

      arma::vec I_k_unscaled = arma::zeros(n_nests);
      for(int j = 0; j < m_i; ++j) {
        const int k = nest_idx0_i[j];
        if (std::isfinite(max_V_k[k])) { // Avoid exp(-inf - (-inf)) -> NaN
            I_k_unscaled[k] += std::exp(V_over_lambda[j] - max_V_k[k]);
        }
      }
      
      arma::vec log_I_k = arma::vec(n_nests).fill(-arma::datum::inf);
      for(int k = 0; k < n_nests; ++k) {
          if (I_k_unscaled[k] > 0) {
              log_I_k[k] = max_V_k[k] + std::log(I_k_unscaled[k]);
          }
      }

      // --- 2b. Calculate log(P_k) (Nest Probability) ---
      // log_P_k = lambda_k * log_I_k - log( sum_l(exp(lambda_l * log_I_l)) )
      arma::vec nest_terms = lambda % log_I_k;
      
      double max_nest_term = nest_terms.max();
      if (!std::isfinite(max_nest_term)) {
        max_nest_term = 0;
      }

      double sum_exp_nest_terms = arma::accu(arma::exp(nest_terms - max_nest_term));
      double log_denom_P_nest;

      if (include_outside_option) {
        // Add outside option: V=0, lambda=1 -> term = 1*log(exp(0/1)) = 0
        sum_exp_nest_terms += std::exp(0.0 - max_nest_term);
      }
      
      log_denom_P_nest = max_nest_term + std::log(sum_exp_nest_terms);
      
      arma::vec log_P_k = nest_terms - log_denom_P_nest;
      double log_P_outside = (include_outside_option) ? (0.0 - log_denom_P_nest) : -arma::datum::inf;

      // --- 2c. Calculate log(P_j|k) and log(P_ij) ---
      // log_P_j_given_k = (V_ij / lambda_k) - log_I_k
      arma::vec log_P_j_given_k = V_over_lambda - log_I_k.elem(nest_idx0_i);
      
      // log_P_i = log_P_j_given_k + log_P_k
      arma::vec log_P_i = log_P_j_given_k + log_P_k.elem(nest_idx0_i);

      // --- 3. Log-Likelihood Calculation ---
      
      const int chosen_alt_idx = choice_idx[i];
      double log_P_choice;
      
      int chosen_nest_k = -1;
      int chosen_inside_idx = -1;
      
      if (include_outside_option && chosen_alt_idx == 0) {
        log_P_choice = log_P_outside;
      } else {
        chosen_inside_idx = chosen_alt_idx - 1; 

        if (chosen_inside_idx < 0 || chosen_inside_idx >= m_i) {
          Rcpp::stop("Invalid chosen alternative index for individual %d", i);
        }
        
        log_P_choice = log_P_i[chosen_inside_idx];
        chosen_nest_k = nest_idx0_i[chosen_inside_idx];
      }

      if (!std::isfinite(log_P_choice)) {
        log_P_choice = -1e10; 
      }
      local_loglik += w_i * log_P_choice;

      // --- 4. Gradient Calculation ---
      arma::vec P_i = arma::exp(log_P_i);
      arma::vec P_j_given_k = arma::exp(log_P_j_given_k);
      arma::vec P_k = arma::exp(log_P_k);

      // 4.1: Pre-calculate sum_{j in B_k} P(j|B_k) * V_ij
      arma::vec sum_P_V_k = arma::zeros(n_nests);
      for (int j = 0; j < m_i; ++j) {
        sum_P_V_k[nest_idx0_i[j]] += P_j_given_k[j] * V_inside[j];
      }

      // 4.2: Gradient w.r.t. beta and delta
      for (int j = 0; j < m_i; ++j) {
        const int k_j = nest_idx0_i[j];
        const double lambda_k_j = lambda[k_j];
        
        double grad_term_j = -P_i[j];
        
        if (k_j == chosen_nest_k) { 
          grad_term_j += P_j_given_k[j] * (1.0 - 1.0 / lambda_k_j);
          if (j == chosen_inside_idx) { 
            grad_term_j += 1.0 / lambda_k_j;
          }
        }
        
        grad_term_j *= w_i;
        
        // beta gradient
        local_grad.subvec(0, K - 1) += grad_term_j * X_i.row(j).t();
        
        // delta gradient
        if (use_asc && delta_length > 0) {
          const int a_id = alt_idx0_i[j];
          if (include_outside_option) {
            local_grad[delta_start_idx + a_id] += grad_term_j;
          } else {
            if (a_id > 0) {
                local_grad[delta_start_idx + (a_id - 1)] += grad_term_j;
            }
          }
        }
      } // end gradient loop for beta/delta 

      // 4.3: Gradient w.r.t. lambda_k
      for (int k = 0; k < n_nests; ++k) {
        
        // Check if this nest 'k' corresponds to an estimated parameter
        const int theta_idx_k = thread_nest_k_to_theta_idx[k];
        if (theta_idx_k == -1) continue;
 
        if (!std::isfinite(log_I_k[k])) continue;

        const double lambda_k = lambda[k];
        double grad_lambda_k = 0.0;
        
        double term_in_brackets = log_I_k[k] - (1.0 / lambda_k) * sum_P_V_k[k];
        
        double P_k_k = P_k[k]; // P(B_k)
        
        if (k == chosen_nest_k) {
          grad_lambda_k += (1.0 - P_k_k) * term_in_brackets;
          
          double V_chosen = V_inside[chosen_inside_idx];
          grad_lambda_k += (1.0 / (lambda_k * lambda_k)) * (sum_P_V_k[k] - V_chosen);
          
        } else {
          grad_lambda_k += -P_k_k * term_in_brackets;
        }
        
        local_grad[theta_idx_k] += w_i * grad_lambda_k;
        
      } // end gradient loop for lambda
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

  // Return negative log-likelihood and gradient for minimization
  return Rcpp::List::create(
    Rcpp::Named("objective") = -global_loglik,
    Rcpp::Named("gradient")  = -global_grad
  );
}

//' Numerical Hessian of the log-likelihood via finite differences
//'
//' @param theta (K + n_delta + n_nests) vector with model parameters.
//'        Order: [beta (K), delta (n_delta), lambda (n_nests)]
//' @param X sum(M) x K design matrix with covariates.
//' @param alt_idx sum(M) x 1 vector with indices of alternatives; 1-based indexing.
//' @param choice_idx N x 1 vector with indices of chosen alternatives; 0 for outside option,
//'        1-based index relative to rows in X_i otherwise.
//' @param nest_idx J x 1 vector with indices of nests for each alternative; 1-based indexing (1 to n_nests).
//' @param M N x 1 vector with number of alternatives for each individual.
//' @param weights N x 1 vector with weights for each observation.
//' @param use_asc whether to use alternative-specific constants.
//' @param include_outside_option whether to include outside option normalized to V=0, lambda=1.
//' @param eps finite difference step size
//' @return Hessian evaluated at input arguments
//' @export
// [[Rcpp::export]]
arma::mat nl_loglik_numeric_hessian(
  const arma::vec& theta,
  const arma::mat& X,
  const arma::uvec& alt_idx,
  const arma::uvec& choice_idx,
  const arma::uvec& nest_idx,
  const Rcpp::IntegerVector& M,
  const arma::vec& weights,
  const bool use_asc = true,
  const bool include_outside_option = false,
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
    Rcpp::List pos_eval = nl_loglik_gradient_parallel(theta_pos, X, alt_idx, choice_idx, nest_idx, M, weights, use_asc, include_outside_option);
    arma::vec grad_pos  = Rcpp::as<arma::vec>(pos_eval["gradient"]);

    // Evaluate gradient at theta_neg
    Rcpp::List neg_eval = nl_loglik_gradient_parallel(theta_neg, X, alt_idx, choice_idx, nest_idx, M, weights, use_asc, include_outside_option);
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
