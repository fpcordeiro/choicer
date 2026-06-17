// [[Rcpp::depends(RcppArmadillo)]]
#include "choicer.h"
#include "choicer_internal.h"

// Forward declaration: shared NL parameter parsing (defined further below,
// after the likelihood kernel; see the full doc comment at the definition).
static void nl_parse_theta(
    const arma::vec& theta,
    const arma::uvec& nest_idx,
    const int K,
    const bool use_asc,
    const bool include_outside_option,
    arma::vec& beta,
    arma::vec& lambda,
    arma::vec& delta,
    int& delta_start_idx,
    int& delta_length,
    arma::ivec* nest_k_to_theta_idx = nullptr);

//' Log-likelihood and gradient for Nested Logit model
//'
//' Computes the log-likelihood and its gradient for the Nested Logit model using OpenMP for parallelization.
//' Especially handles singleton nests by fixing their lambda parameters to 1. Only non-singleton nests have a inclusive value coefficient estimated in theta.
//'
//' @param theta (K + n_non_singleton_nests + n_delta) vector with model parameters.
//'        Order: `[beta (K), lambda (n_non_singleton_nests), delta (n_delta)]`
//' @param X sum(M) x K design matrix with covariates.
//' @param alt_idx sum(M) x 1 vector with indices of alternatives; 1-based indexing.
//' @param choice_idx N x 1 vector with indices of chosen alternatives; 0 for outside option,
//'        1-based index relative to rows in X_i otherwise.
//' @param nest_idx J x 1 vector with indices of nests for each alternative; 1-based indexing (1 to n_nests).
//' @param M N x 1 vector with number of alternatives for each individual.
//' @param weights N x 1 vector with weights for each observation.
//' @param use_asc whether to use alternative-specific constants.
//' @param include_outside_option whether to include outside option normalized to V=0, lambda=1.
//' @returns List with loglikelihood and gradient evaluated at input arguments
//' @examples
//' \donttest{
//' library(data.table)
//' set.seed(42)
//' N <- 50; J <- 4
//' dt <- data.table(id = rep(1:N, each = J), alt = rep(1:J, N))
//' dt[, `:=`(x1 = rnorm(.N), x2 = rnorm(.N))]
//' dt[, nest := ifelse(alt <= 2, "A", "B")]
//' dt[, choice := 0L]
//' dt[, choice := sample(c(1L, rep(0L, J - 1))), by = id]
//' d <- prepare_nl_data(dt, "id", "alt", "choice", c("x1", "x2"), "nest")
//' K_x <- ncol(d$X); K_l <- length(unique(d$nest_idx))
//' theta <- c(rep(0, K_x), rep(0.5, K_l), rep(0, J - 1))
//' result <- nl_loglik_gradient_parallel(theta, d$X, d$alt_idx,
//'   d$choice_idx, d$nest_idx, d$M, d$weights)
//' result$objective
//' }
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

  // --- 1. Parameter Parsing (shared helper; also returns the map from full
  // nest index k to lambda_k's position in theta, -1 for singleton nests) ---
  arma::vec beta, lambda, delta;
  int delta_start_idx, delta_length;
  arma::ivec nest_k_to_theta_idx;
  nl_parse_theta(theta, nest_idx, K, use_asc, include_outside_option,
                 beta, lambda, delta, delta_start_idx, delta_length,
                 &nest_k_to_theta_idx);
  validate_nl_inputs(X, alt_idx, nest_idx, M, use_asc, delta,
                     &weights, &choice_idx);

  // 0-based indexing for inputs
  arma::uvec alt_idx0 = alt_idx - 1;
  arma::uvec nest_idx0 = nest_idx - 1;
  
  // Compute prefix sums for indexing
  const Rcpp::IntegerVector S = compute_prefix_sum(M);

  // Pre-compute base utility for all individuals (single BLAS call)
  arma::vec base_util = compute_base_util(X, beta, alt_idx0, use_asc, delta);

  // --- H2: Serial pre-loop validation of chosen-alternative indices ---
  // Rcpp::stop() is only safe outside parallel regions.
  for (int i = 0; i < N; ++i) {
    const int chosen_alt_idx_check = choice_idx[i];
    if (include_outside_option && chosen_alt_idx_check == 0) continue; // outside option is valid
    const int chosen_inside = chosen_alt_idx_check - 1;
    if (chosen_inside < 0 || chosen_inside >= M[i]) {
      Rcpp::stop("Invalid chosen alternative index for individual %d", i);
    }
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
    arma::vec grad_vec; // pre-allocated per-thread, resized per individual

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

      // V_ij = X_ij * beta + delta_j (pre-computed)
      arma::vec V_inside = base_util.subvec(start_idx, end_idx);

      // Per-individual probability block (shared helper)
      arma::vec P_i, P_j_given_k, P_k, log_I_k, log_P_i;
      double log_P_outside;
      nl_individual_probs(V_inside, nest_idx0_i, lambda, n_nests,
                          include_outside_option,
                          P_i, P_j_given_k, P_k, log_I_k, log_P_i, log_P_outside);

      // --- 3. Log-Likelihood Calculation ---

      const int chosen_alt_idx = choice_idx[i];
      double log_P_choice;

      int chosen_nest_k = -1;
      int chosen_inside_idx = -1;

      if (include_outside_option && chosen_alt_idx == 0) {
        log_P_choice = log_P_outside;
      } else {
        // validated serially above; no in-loop Rcpp::stop needed
        chosen_inside_idx = chosen_alt_idx - 1;
        log_P_choice = log_P_i[chosen_inside_idx];
        chosen_nest_k = nest_idx0_i[chosen_inside_idx];
      }

      if (!std::isfinite(log_P_choice)) {
        log_P_choice = -1e10; // H1: clamp without in-loop R-API call
      }
      local_loglik += w_i * log_P_choice;

      // --- 4. Gradient Calculation ---

      // 4.1: Pre-calculate sum_{j in B_k} P(j|B_k) * V_ij
      arma::vec sum_P_V_k = arma::zeros(n_nests);
      for (int j = 0; j < m_i; ++j) {
        sum_P_V_k[nest_idx0_i[j]] += P_j_given_k[j] * V_inside[j];
      }

      // 4.2: Gradient w.r.t. beta and delta
      // Build grad_vec (NL analogue of diff_vec in MNL):
      //   grad_vec[j] = -P_i[j]
      //               + P_j_given_k[j] * (1 - 1/lambda_k)  if j in chosen nest
      //               + 1/lambda_k                          if j == chosen alt
      grad_vec.set_size(m_i);
      grad_vec = -P_i;
      if (chosen_nest_k >= 0) { // not outside option
        const double lambda_chosen = lambda[chosen_nest_k];
        const double term_scale = 1.0 - 1.0 / lambda_chosen;
        for (int j = 0; j < m_i; ++j) {
          if ((int)nest_idx0_i[j] == chosen_nest_k) {
            grad_vec[j] += P_j_given_k[j] * term_scale;
          }
        }
        grad_vec[chosen_inside_idx] += 1.0 / lambda_chosen;
      }

      // Beta block: single BLAS dgemv (analogous to MNL)
      local_grad.subvec(0, K - 1) += w_i * X_i.t() * grad_vec;

      // Delta block: scatter loop (irregular alt-index mapping)
      if (use_asc && delta_length > 0) {
        scatter_delta_grad(local_grad, delta_start_idx, grad_vec, alt_idx0_i,
                           m_i, include_outside_option, w_i);
      } // end gradient block for beta/delta

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

  // H3: Sanitize NaN/Inf (matching mxl_loglik_gradient_parallel)
  double obj = -global_loglik;
  arma::vec grad = -global_grad;
  if (!std::isfinite(obj)) {
    obj = 1e10;
    grad.zeros();
  } else {
    grad.elem(arma::find_nonfinite(grad)).zeros();
  }
  return Rcpp::List::create(
    Rcpp::Named("objective") = obj,
    Rcpp::Named("gradient")  = grad
  );
}

//' Numerical Hessian of the log-likelihood via finite differences
//'
//' @param theta (K + n_delta + n_nests) vector with model parameters.
//'        Order: `[beta (K), delta (n_delta), lambda (n_nests)]`
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
//' @returns Hessian evaluated at input arguments
//' @examples
//' \donttest{
//' library(data.table)
//' set.seed(42)
//' N <- 50; J <- 4
//' dt <- data.table(id = rep(1:N, each = J), alt = rep(1:J, N))
//' dt[, `:=`(x1 = rnorm(.N), x2 = rnorm(.N))]
//' dt[, nest := ifelse(alt <= 2, "A", "B")]
//' dt[, choice := 0L]
//' dt[, choice := sample(c(1L, rep(0L, J - 1))), by = id]
//' d <- prepare_nl_data(dt, "id", "alt", "choice", c("x1", "x2"), "nest")
//' K_x <- ncol(d$X); K_l <- length(unique(d$nest_idx))
//' theta <- c(rep(0, K_x), rep(0.5, K_l), rep(0, J - 1))
//' H <- nl_loglik_numeric_hessian(theta, d$X, d$alt_idx, d$choice_idx,
//'   d$nest_idx, d$M, d$weights)
//' dim(H)
//' }
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

// =============================================================================
// Shared NL parameter parsing (file-local)
//
// Given theta in native layout [beta(K), lambda(non-singleton), delta], the
// nest index of each alternative, and the ASC/outside-option flags, this
// resolves:
//   beta        (K)         fixed coefficients
//   lambda      (n_nests)   full lambda vector (singletons fixed to 1)
//   delta       (J or J-1)  alternative-specific constants (index 0 = first
//                           inside alt, fixed to 0 when no outside option)
//   delta_start_idx         offset of delta block in theta
//   delta_length            number of *estimated* delta parameters
//   nest_k_to_theta_idx     (optional) map from full nest index k to
//                           lambda_k's position in theta; -1 for singletons
// Used by every entry point, including nl_loglik_gradient_parallel.
static void nl_parse_theta(
    const arma::vec& theta,
    const arma::uvec& nest_idx,
    const int K,
    const bool use_asc,
    const bool include_outside_option,
    arma::vec& beta,
    arma::vec& lambda,
    arma::vec& delta,
    int& delta_start_idx,
    int& delta_length,
    arma::ivec* nest_k_to_theta_idx
) {
  const int n_params = theta.n_elem;
  if (K <= 0) {
    Rcpp::stop("K must be positive, got %d", K);
  }
  if (n_params < K) {
    Rcpp::stop("Theta vector too short: missing beta parameters "
               "(expected at least %d, got %d).", K, n_params);
  }
  if (nest_idx.n_elem == 0) {
    Rcpp::stop("nest_idx must be non-empty.");
  }
  const int n_nests = arma::max(nest_idx); // 1-based nest_idx

  beta = theta.subvec(0, K - 1);

  // Identify singleton nests
  arma::uvec nest_counts = arma::zeros<arma::uvec>(n_nests);
  for (unsigned int j = 0; j < nest_idx.n_elem; ++j) {
    if (nest_idx[j] > 0 && (int)nest_idx[j] <= n_nests) {
      nest_counts[nest_idx[j] - 1]++;
    } else {
      Rcpp::stop("Invalid nest index found in nest_idx.");
    }
  }
  arma::uvec is_singleton = (nest_counts == 1);
  const int n_non_singleton_nests = arma::accu(is_singleton == 0);

  const int lambda_start_idx = K;
  lambda = arma::ones(n_nests);

  // Optional output: map from full nest index k to lambda_k's position in
  // theta; -1 marks singleton nests (lambda fixed to 1, not a parameter).
  if (nest_k_to_theta_idx) {
    nest_k_to_theta_idx->set_size(n_nests);
    nest_k_to_theta_idx->fill(-1);
  }

  if (n_non_singleton_nests > 0) {
    if (theta.n_elem < (unsigned)(K + n_non_singleton_nests)) {
      Rcpp::stop("Error: theta vector is too short for K + n_non_singleton_nests.");
    }
    arma::vec non_singleton_lambdas =
      theta.subvec(lambda_start_idx, lambda_start_idx + n_non_singleton_nests - 1);
    if (arma::any(non_singleton_lambdas <= 0)) {
      Rcpp::stop("Error: All non-singleton lambda (nest) parameters must be > 0.");
    }
    int current_lambda_idx = 0;
    for (int k = 0; k < n_nests; ++k) {
      if (is_singleton[k] == 0) {
        lambda[k] = non_singleton_lambdas[current_lambda_idx];
        if (nest_k_to_theta_idx) {
          (*nest_k_to_theta_idx)[k] = lambda_start_idx + current_lambda_idx;
        }
        current_lambda_idx++;
      }
    }
  } else {
    Rcpp::stop("Error: No non-singleton nests found. At least one nest must have multiple alternatives.");
  }

  delta_start_idx = K + n_non_singleton_nests;

  if (use_asc) {
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
      delta = arma::zeros((include_outside_option) ? 0 : 1);
    }
  } else {
    if (n_params != delta_start_idx) {
      Rcpp::stop("Theta vector too long: %d parameters given but the model "
                 "expects %d. Did you mean use_asc = TRUE?",
                 n_params, delta_start_idx);
    }
    delta_length = 0;
    delta = arma::zeros(0);
  }
}

//' Prediction of choice probabilities and utilities for the Nested Logit model
//'
//' @param theta (K + n_non_singleton_nests + n_delta) vector with model
//'        parameters. Order: `[beta (K), lambda (non-singleton), delta]`.
//' @param X sum(M) x K design matrix with covariates.
//' @param alt_idx sum(M) x 1 vector with indices of alternatives; 1-based indexing.
//' @param M N x 1 vector with number of alternatives for each individual.
//' @param nest_idx J x 1 vector with nest indices for each alternative; 1-based indexing.
//' @param use_asc whether to use alternative-specific constants.
//' @param include_outside_option whether to include outside option normalized to V=0, lambda=1.
//' @returns List with `choice_prob` (joint P_ij per stacked row) and `utility` (V_ij).
//' @examples
//' \donttest{
//' library(data.table)
//' set.seed(42)
//' N <- 50; J <- 4
//' dt <- data.table(id = rep(1:N, each = J), alt = rep(1:J, N))
//' dt[, `:=`(x1 = rnorm(.N), x2 = rnorm(.N))]
//' dt[, nest := ifelse(alt <= 2, "A", "B")]
//' dt[, choice := 0L]
//' dt[, choice := sample(c(1L, rep(0L, J - 1))), by = id]
//' fit <- run_nestlogit(dt, "id", "alt", "choice", c("x1", "x2"), "nest")
//' pred <- nl_predict(coef(fit), fit$data$X, fit$data$alt_idx, fit$data$M,
//'   fit$data$nest_idx, use_asc = TRUE)
//' head(pred$choice_prob)
//' }
//' @export
// [[Rcpp::export]]
Rcpp::List nl_predict(
    const arma::vec& theta,
    const arma::mat& X,
    const arma::uvec& alt_idx,
    const Rcpp::IntegerVector& M,
    const arma::uvec& nest_idx,
    const bool use_asc = true,
    const bool include_outside_option = false
) {
  const int N = M.size();
  const int K = X.n_cols;
  const int n_nests = arma::max(nest_idx);

  arma::vec beta, lambda, delta;
  int delta_start_idx, delta_length;
  nl_parse_theta(theta, nest_idx, K, use_asc, include_outside_option,
                 beta, lambda, delta, delta_start_idx, delta_length);
  validate_nl_inputs(X, alt_idx, nest_idx, M, use_asc, delta);

  arma::uvec alt_idx0 = alt_idx - 1;
  arma::uvec nest_idx0 = nest_idx - 1;
  const Rcpp::IntegerVector S = compute_prefix_sum(M);

  arma::vec base_util = compute_base_util(X, beta, alt_idx0, use_asc, delta);

  arma::vec V_all = arma::zeros(X.n_rows);
  arma::vec P_all = arma::zeros(X.n_rows);

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
  for (int i = 0; i < N; ++i) {
    const int m_i       = M[i];
    const int start_idx = S[i];
    const int end_idx   = start_idx + m_i - 1;

    arma::uvec alt_idx0_i = alt_idx0.subvec(start_idx, end_idx);
    arma::uvec nest_idx0_i = nest_idx0.elem(alt_idx0_i);
    arma::vec V_inside = base_util.subvec(start_idx, end_idx);

    arma::vec P_i, P_j_given_k, P_k, log_I_k, log_P_i;
    double log_P_outside;
    nl_individual_probs(V_inside, nest_idx0_i, lambda, n_nests,
                        include_outside_option,
                        P_i, P_j_given_k, P_k, log_I_k, log_P_i, log_P_outside);

    V_all.subvec(start_idx, end_idx) = V_inside;
    P_all.subvec(start_idx, end_idx) = P_i;
  }

  return Rcpp::List::create(
    Rcpp::Named("choice_prob") = P_all,
    Rcpp::Named("utility")     = V_all
  );
}

// Predicted NL market shares (file-local internal).
// alt_idx0/nest_idx0 are 0-based. delta is the *full* delta vector
// (index 0 = first inside alt). Returns weighted shares of length num_alts
// (index 0 = outside option when present, then inside alts in order).
static arma::vec nl_predict_shares_internal(
    const arma::mat& X,
    const arma::vec& beta,
    const arma::vec& lambda,
    const arma::uvec& alt_idx0,
    const arma::uvec& nest_idx0,
    const Rcpp::IntegerVector& M,
    const Rcpp::IntegerVector& S,
    const arma::vec& weights,
    const arma::vec& delta,
    const int n_nests,
    const int num_alts,
    const bool use_asc,
    const bool include_outside_option
) {
  const int N = M.size();
  const double denominator = arma::sum(weights);
  if (denominator <= 0) {
    Rcpp::stop("Error: Sum of weights must be positive.");
  }

  arma::vec base_util = compute_base_util(X, beta, alt_idx0, use_asc, delta);

  arma::vec global_shares = arma::zeros(num_alts);

#ifdef _OPENMP
#pragma omp parallel
#endif
  {
    arma::vec local_shares = arma::zeros(num_alts);

#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
    for (int i = 0; i < N; ++i) {
      const int m_i       = M[i];
      const int start_idx = S[i];
      const int end_idx   = start_idx + m_i - 1;
      const double w_i    = weights[i];

      arma::uvec alt_idx0_i = alt_idx0.subvec(start_idx, end_idx);
      arma::uvec nest_idx0_i = nest_idx0.elem(alt_idx0_i);
      arma::vec V_inside = base_util.subvec(start_idx, end_idx);

      arma::vec P_i, P_j_given_k, P_k, log_I_k, log_P_i;
      double log_P_outside;
      nl_individual_probs(V_inside, nest_idx0_i, lambda, n_nests,
                          include_outside_option,
                          P_i, P_j_given_k, P_k, log_I_k, log_P_i, log_P_outside);

      if (include_outside_option) {
        local_shares(0) += w_i * std::exp(log_P_outside);
      }
      for (int a = 0; a < m_i; ++a) {
        if (include_outside_option) {
          local_shares(alt_idx0_i(a) + 1) += w_i * P_i(a);
        } else {
          local_shares(alt_idx0_i(a)) += w_i * P_i(a);
        }
      }
    }

#ifdef _OPENMP
#pragma omp critical
#endif
    {
      global_shares += local_shares;
    }
  }

  return global_shares / denominator;
}

//' Prediction of market shares for the Nested Logit model
//'
//' @param theta (K + n_non_singleton_nests + n_delta) vector with model
//'        parameters. Order: `[beta (K), lambda (non-singleton), delta]`.
//' @param X sum(M) x K design matrix with covariates.
//' @param alt_idx sum(M) x 1 vector with indices of alternatives; 1-based indexing.
//' @param M N x 1 vector with number of alternatives for each individual.
//' @param weights N x 1 vector with weights for each observation.
//' @param nest_idx J x 1 vector with nest indices for each alternative; 1-based indexing.
//' @param use_asc whether to use alternative-specific constants.
//' @param include_outside_option whether to include outside option normalized to V=0, lambda=1.
//' @returns vector with predicted market shares (outside-option share first when present).
//' @examples
//' \donttest{
//' library(data.table)
//' set.seed(42)
//' N <- 50; J <- 4
//' dt <- data.table(id = rep(1:N, each = J), alt = rep(1:J, N))
//' dt[, `:=`(x1 = rnorm(.N), x2 = rnorm(.N))]
//' dt[, nest := ifelse(alt <= 2, "A", "B")]
//' dt[, choice := 0L]
//' dt[, choice := sample(c(1L, rep(0L, J - 1))), by = id]
//' fit <- run_nestlogit(dt, "id", "alt", "choice", c("x1", "x2"), "nest")
//' shares <- nl_predict_shares(coef(fit), fit$data$X, fit$data$alt_idx,
//'   fit$data$M, fit$data$weights, fit$data$nest_idx, use_asc = TRUE)
//' shares
//' }
//' @export
// [[Rcpp::export]]
arma::vec nl_predict_shares(
    const arma::vec& theta,
    const arma::mat& X,
    const arma::uvec& alt_idx,
    const Rcpp::IntegerVector& M,
    const arma::vec& weights,
    const arma::uvec& nest_idx,
    const bool use_asc = true,
    const bool include_outside_option = false
) {
  const int K = X.n_cols;
  const int n_nests = arma::max(nest_idx);

  arma::vec beta, lambda, delta;
  int delta_start_idx, delta_length;
  nl_parse_theta(theta, nest_idx, K, use_asc, include_outside_option,
                 beta, lambda, delta, delta_start_idx, delta_length);
  validate_nl_inputs(X, alt_idx, nest_idx, M, use_asc, delta, &weights);

  arma::uvec alt_idx0 = alt_idx - 1;
  arma::uvec nest_idx0 = nest_idx - 1;
  const Rcpp::IntegerVector S = compute_prefix_sum(M);

  int num_alts = include_outside_option ? (alt_idx.max() + 1) : alt_idx.max();

  return nl_predict_shares_internal(
    X, beta, lambda, alt_idx0, nest_idx0, M, S, weights, delta,
    n_nests, num_alts, use_asc, include_outside_option
  );
}

//' Compute aggregate elasticities for the Nested Logit model
//'
//' Computes the aggregate (weighted-average) elasticity matrix for the Nested
//' Logit model. Reduces to the MNL elasticities when all lambda = 1.
//'
//' @param theta (K + n_non_singleton_nests + n_delta) vector with model
//'        parameters. Order: `[beta (K), lambda (non-singleton), delta]`.
//' @param X sum(M) x K design matrix with covariates.
//' @param alt_idx sum(M) x 1 vector with indices of alternatives; 1-based indexing.
//' @param choice_idx N x 1 vector (kept for API consistency, not used).
//' @param nest_idx J x 1 vector with nest indices for each alternative; 1-based indexing.
//' @param M N x 1 vector with number of alternatives for each individual.
//' @param weights N x 1 vector with weights for each observation.
//' @param elast_var_idx 1-based index of the column in X for which to compute the elasticity.
//' @param use_asc whether to use alternative-specific constants.
//' @param include_outside_option whether to include outside option normalized to V=0, lambda=1.
//' @returns J x J matrix of aggregate elasticities (row = responding alt, col = perturbed alt).
//' @examples
//' \donttest{
//' library(data.table)
//' set.seed(42)
//' N <- 50; J <- 4
//' dt <- data.table(id = rep(1:N, each = J), alt = rep(1:J, N))
//' dt[, `:=`(x1 = rnorm(.N), x2 = rnorm(.N))]
//' dt[, nest := ifelse(alt <= 2, "A", "B")]
//' dt[, choice := 0L]
//' dt[, choice := sample(c(1L, rep(0L, J - 1))), by = id]
//' fit <- run_nestlogit(dt, "id", "alt", "choice", c("x1", "x2"), "nest")
//' elas <- nl_elasticities_parallel(coef(fit), fit$data$X, fit$data$alt_idx,
//'   fit$data$choice_idx, fit$data$nest_idx, fit$data$M, fit$data$weights,
//'   elast_var_idx = 1L)
//' elas
//' }
//' @export
// [[Rcpp::export]]
arma::mat nl_elasticities_parallel(
    const arma::vec& theta,
    const arma::mat& X,
    const arma::uvec& alt_idx,
    const arma::uvec& choice_idx, // kept for consistency, not used
    const arma::uvec& nest_idx,
    const Rcpp::IntegerVector& M,
    const arma::vec& weights,
    const int elast_var_idx,
    const bool use_asc = true,
    const bool include_outside_option = false
) {
  (void)choice_idx;
  const int N = M.size();
  const int K = X.n_cols;
  const int n_nests = arma::max(nest_idx);

  const int var_idx = elast_var_idx - 1;
  if (var_idx < 0 || var_idx >= K) {
    Rcpp::stop("elast_var_idx is out of bounds for design matrix X.");
  }

  arma::vec beta, lambda, delta;
  int delta_start_idx, delta_length;
  nl_parse_theta(theta, nest_idx, K, use_asc, include_outside_option,
                 beta, lambda, delta, delta_start_idx, delta_length);
  validate_nl_inputs(X, alt_idx, nest_idx, M, use_asc, delta, &weights);
  const double beta_k = beta(var_idx);

  arma::uvec alt_idx0 = alt_idx - 1;
  arma::uvec nest_idx0 = nest_idx - 1;

  const int J_inside = compute_J_inside(use_asc, delta, alt_idx0);
  const int J_total = compute_J_total(J_inside, include_outside_option);

  const Rcpp::IntegerVector S = compute_prefix_sum(M);

  arma::vec base_util = compute_base_util(X, beta, alt_idx0, use_asc, delta);

  arma::mat global_elas_matrix = arma::zeros(J_total, J_total);
  double global_total_weight = 0.0;

#ifdef _OPENMP
#pragma omp parallel
#endif
  {
    arma::mat local_elas_matrix = arma::zeros(J_total, J_total);
    double local_total_weight = 0.0;

#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
    for (int i = 0; i < N; ++i) {
      const int m_i       = M[i];
      const int start_idx = S[i];
      const int end_idx   = start_idx + m_i - 1;
      const double w_i    = weights[i];
      const auto X_i      = X.rows(start_idx, end_idx);

      arma::uvec alt_idx0_i = alt_idx0.subvec(start_idx, end_idx);
      arma::uvec nest_idx0_i = nest_idx0.elem(alt_idx0_i);
      arma::vec V_inside = base_util.subvec(start_idx, end_idx);

      arma::vec P_i, P_j_given_k, P_k, log_I_k, log_P_i;
      double log_P_outside;
      nl_individual_probs(V_inside, nest_idx0_i, lambda, n_nests,
                          include_outside_option,
                          P_i, P_j_given_k, P_k, log_I_k, log_P_i, log_P_outside);

      // x_k for each inside alternative
      arma::vec x_k_i = X_i.col(var_idx);

      // Global alternative index for each inside alt (outside option = 0)
      arma::uvec global_map =
          build_global_alt_map_inside(alt_idx0_i, include_outside_option);

      // Elasticity of P_ij (row j) w.r.t. attribute x of alt a (col a):
      //   E_ja = beta_k * x_{ia} * d log P_ij / d V_ia
      // d log P_ij/dV_ia = 1{a=j}/lambda_r + (1-1/lambda_r) 1{s=r} P(a|r) - P_ia
      for (int j = 0; j < m_i; ++j) {
        const int global_j = global_map[j];
        const int r = nest_idx0_i[j];   // nest of responding alt j
        const double lam_r = lambda[r];

        for (int a = 0; a < m_i; ++a) {
          const int global_a = global_map[a];
          const int s = nest_idx0_i[a]; // nest of perturbed alt a
          const double x_ia = x_k_i[a];
          const double P_ia = P_i[a];

          double dlogP;
          if (a == j) {
            // own
            dlogP = 1.0 / lam_r + (1.0 - 1.0 / lam_r) * P_j_given_k[a] - P_ia;
          } else if (s == r) {
            // cross, same nest
            dlogP = (1.0 - 1.0 / lam_r) * P_j_given_k[a] - P_ia;
          } else {
            // cross, different nest
            dlogP = -P_ia;
          }

          const double elasticity = beta_k * x_ia * dlogP;
          local_elas_matrix(global_j, global_a) += w_i * elasticity;
        }
      }

      // Outside-option row (responding alt = outside option, global index 0).
      // Outside has V=0, its own singleton nest -> for any inside perturbed
      // alt a: d log P_out / d V_ia = -P_ia (different nest, a != out).
      // (The outside-option column is 0 because its own x_k = 0.)
      if (include_outside_option) {
        for (int a = 0; a < m_i; ++a) {
          const int global_a = global_map[a];
          const double x_ia = x_k_i[a];
          const double P_ia = P_i[a];
          const double elasticity = beta_k * x_ia * (-P_ia);
          local_elas_matrix(0, global_a) += w_i * elasticity;
        }
      }

      local_total_weight += w_i;
    }

#ifdef _OPENMP
#pragma omp critical
#endif
    {
      global_elas_matrix += local_elas_matrix;
      global_total_weight += local_total_weight;
    }
  }

  if (global_total_weight > 1e-10) {
    global_elas_matrix /= global_total_weight;
  }

  return global_elas_matrix;
}

//' Compute Nested Logit diversion ratios (parallelized over individuals)
//'
//' Computes the diversion ratio matrix DR(j->k) for the Nested Logit model.
//' Entry (k, j) = fraction of demand lost by alternative j captured by k.
//' Reduces to the MNL diversion ratios when all lambda = 1.
//'
//' @param theta (K + n_non_singleton_nests + n_delta) vector with model
//'        parameters. Order: `[beta (K), lambda (non-singleton), delta]`.
//' @param X sum(M) x K design matrix with covariates.
//' @param alt_idx sum(M) x 1 vector with indices of alternatives; 1-based indexing.
//' @param nest_idx J x 1 vector with nest indices for each alternative; 1-based indexing.
//' @param M N x 1 vector with number of alternatives for each individual.
//' @param weights N x 1 vector with weights for each observation.
//' @param use_asc whether to use alternative-specific constants.
//' @param include_outside_option whether to include outside option normalized to V=0, lambda=1.
//' @returns J x J matrix where entry (k, j) = DR(j->k). Diagonal is 0.
//' @examples
//' \donttest{
//' library(data.table)
//' set.seed(42)
//' N <- 50; J <- 4
//' dt <- data.table(id = rep(1:N, each = J), alt = rep(1:J, N))
//' dt[, `:=`(x1 = rnorm(.N), x2 = rnorm(.N))]
//' dt[, nest := ifelse(alt <= 2, "A", "B")]
//' dt[, choice := 0L]
//' dt[, choice := sample(c(1L, rep(0L, J - 1))), by = id]
//' fit <- run_nestlogit(dt, "id", "alt", "choice", c("x1", "x2"), "nest")
//' dr <- nl_diversion_ratios_parallel(coef(fit), fit$data$X, fit$data$alt_idx,
//'   fit$data$nest_idx, fit$data$M, fit$data$weights)
//' dr
//' }
//' @export
// [[Rcpp::export]]
arma::mat nl_diversion_ratios_parallel(
    const arma::vec& theta,
    const arma::mat& X,
    const arma::uvec& alt_idx,
    const arma::uvec& nest_idx,
    const Rcpp::IntegerVector& M,
    const arma::vec& weights,
    const bool use_asc = true,
    const bool include_outside_option = false
) {
  const int N = M.size();
  const int K = X.n_cols;
  const int n_nests = arma::max(nest_idx);

  arma::vec beta, lambda, delta;
  int delta_start_idx, delta_length;
  nl_parse_theta(theta, nest_idx, K, use_asc, include_outside_option,
                 beta, lambda, delta, delta_start_idx, delta_length);
  validate_nl_inputs(X, alt_idx, nest_idx, M, use_asc, delta, &weights);

  arma::uvec alt_idx0 = alt_idx - 1;
  arma::uvec nest_idx0 = nest_idx - 1;

  const int J_inside = compute_J_inside(use_asc, delta, alt_idx0);
  const int J_total = compute_J_total(J_inside, include_outside_option);

  const Rcpp::IntegerVector S = compute_prefix_sum(M);

  arma::vec base_util = compute_base_util(X, beta, alt_idx0, use_asc, delta);

  // numerator(k, j) = sum_i w_i * (-dP_ik/dV_ij), k != j
  // denominator(j)  = sum_i w_i * (dP_ij/dV_ij)
  arma::mat global_numerator = arma::zeros(J_total, J_total);
  arma::vec global_denominator = arma::zeros(J_total);

#ifdef _OPENMP
#pragma omp parallel
#endif
  {
    arma::mat local_numerator = arma::zeros(J_total, J_total);
    arma::vec local_denominator = arma::zeros(J_total);

#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
    for (int i = 0; i < N; ++i) {
      const int m_i       = M[i];
      const int start_idx = S[i];
      const int end_idx   = start_idx + m_i - 1;
      const double w_i    = weights[i];

      arma::uvec alt_idx0_i = alt_idx0.subvec(start_idx, end_idx);
      arma::uvec nest_idx0_i = nest_idx0.elem(alt_idx0_i);
      arma::vec V_inside = base_util.subvec(start_idx, end_idx);

      arma::vec P_i, P_j_given_k, P_k, log_I_k, log_P_i;
      double log_P_outside;
      nl_individual_probs(V_inside, nest_idx0_i, lambda, n_nests,
                          include_outside_option,
                          P_i, P_j_given_k, P_k, log_I_k, log_P_i, log_P_outside);

      const double P_out = include_outside_option ? std::exp(log_P_outside) : 0.0;

      // Global alternative index for each inside alt (outside option = 0)
      arma::uvec global_map =
          build_global_alt_map_inside(alt_idx0_i, include_outside_option);

      // For each perturbed alt j (column), accumulate response of every alt.
      // dP_im/dV_ij = P_im * (d log P_im/d V_ij), with r = nest of j:
      //   own (m=j):            P_ij [ 1/lam_r + (1-1/lam_r) P(j|r) - P_ij ]
      //   cross same nest:      P_im [ (1-1/lam_r) P(j|r) - P_ij ]
      //   cross diff nest:      -P_im P_ij
      for (int j = 0; j < m_i; ++j) {
        const int global_j = global_map[j];
        const int r = nest_idx0_i[j];
        const double lam_r = lambda[r];
        const double P_ij = P_i[j];
        const double Pj_given_r = P_j_given_k[j];

        // denominator (own derivative)
        const double dP_own =
          P_ij * (1.0 / lam_r + (1.0 - 1.0 / lam_r) * Pj_given_r - P_ij);
        local_denominator(global_j) += w_i * dP_own;

        // numerator: inside alts m != j
        for (int m = 0; m < m_i; ++m) {
          if (m == j) continue;
          const int global_m = global_map[m];
          const int s = nest_idx0_i[m];
          const double P_im = P_i[m];

          double dP_cross;
          if (s == r) {
            dP_cross = P_im * ((1.0 - 1.0 / lam_r) * Pj_given_r - P_ij);
          } else {
            dP_cross = -P_im * P_ij;
          }
          local_numerator(global_m, global_j) += w_i * (-dP_cross);
        }

        // outside option as a *destination* (always a different nest from j)
        if (include_outside_option) {
          const double dP_out_cross = -P_out * P_ij;
          local_numerator(0, global_j) += w_i * (-dP_out_cross);
        }
      }

      // outside option as a *source* (its own singleton nest, lambda = 1):
      //   own:   dP_out/dV_out = P_out (1 - P_out)
      //   cross: dP_im/dV_out  = -P_im P_out  (different nest)
      if (include_outside_option) {
        local_denominator(0) += w_i * P_out * (1.0 - P_out);
        for (int m = 0; m < m_i; ++m) {
          const int global_m = global_map[m];
          const double P_im = P_i[m];
          const double dP_cross = -P_im * P_out;
          local_numerator(global_m, 0) += w_i * (-dP_cross);
        }
      }
    }

#ifdef _OPENMP
#pragma omp critical
#endif
    {
      global_numerator += local_numerator;
      global_denominator += local_denominator;
    }
  }

  arma::mat diversion_matrix = arma::zeros(J_total, J_total);
  for (int j = 0; j < J_total; ++j) {
    if (global_denominator(j) > 1e-15) {
      for (int k = 0; k < J_total; ++k) {
        if (k != j) {
          diversion_matrix(k, j) = global_numerator(k, j) / global_denominator(j);
        }
      }
    }
  }

  return diversion_matrix;
}

//' BLP95 contraction mapping for the Nested Logit model
//'
//' Damped iterative fixed point recovering delta given target shares, using the
//' NL probability structure. `damping = 1` reproduces the plain BLP update.
//'
//' @param delta J x 1 vector with initial guess for deltas (ASCs).
//' @param target_shares vector with target shares (outside-option share first when present).
//' @param X sum(M) x K design matrix with covariates.
//' @param beta K x 1 vector with fixed coefficients.
//' @param lambda full nest dissimilarity vector of length n_nests (singletons = 1).
//' @param alt_idx sum(M) x 1 vector with indices of alternatives; 1-based indexing.
//' @param nest_idx J x 1 vector with nest indices for each alternative; 1-based indexing.
//' @param M N x 1 vector with number of alternatives for each individual.
//' @param weights N x 1 vector with weights for each observation.
//' @param include_outside_option whether to include outside option normalized to V=0, lambda=1.
//' @param damping damping factor for the update (default 1.0 = plain BLP).
//' @param tol convergence tolerance.
//' @param max_iter maximum number of iterations.
//' @returns vector with contraction's delta (ASCs) output.
//' @examples
//' \donttest{
//' library(data.table)
//' set.seed(42)
//' N <- 50; J <- 4
//' dt <- data.table(id = rep(1:N, each = J), alt = rep(1:J, N))
//' dt[, `:=`(x1 = rnorm(.N), x2 = rnorm(.N))]
//' dt[, nest := ifelse(alt <= 2, "A", "B")]
//' dt[, choice := 0L]
//' dt[, choice := sample(c(1L, rep(0L, J - 1))), by = id]
//' fit <- run_nestlogit(dt, "id", "alt", "choice", c("x1", "x2"), "nest")
//' beta <- coef(fit)[fit$param_map$beta]
//' lambda <- rep(1, length(unique(fit$data$nest_idx)))
//' lambda[as.integer(names(which(table(fit$data$nest_idx) > 1)))] <-
//'   coef(fit)[fit$param_map$lambda]
//' delta <- nl_blp_contraction(rep(0, J), rep(1/J, J), fit$data$X, beta, lambda,
//'   fit$data$alt_idx, fit$data$nest_idx, fit$data$M, fit$data$weights)
//' delta
//' }
//' @export
// [[Rcpp::export]]
arma::vec nl_blp_contraction(
    const arma::vec& delta,
    const arma::vec& target_shares,
    const arma::mat& X,
    const arma::vec& beta,
    const arma::vec& lambda,
    const arma::uvec& alt_idx,
    const arma::uvec& nest_idx,
    const Rcpp::IntegerVector& M,
    const arma::vec& weights,
    const bool include_outside_option = false,
    const double damping  = 1.0,
    const double tol      = 1e-8,
    const int    max_iter = 1000
) {
  const bool use_asc = true;
  const int n_nests = arma::max(nest_idx);

  int num_alts = include_outside_option ? (delta.n_elem + 1) : delta.n_elem;
  if ((int)target_shares.n_elem != num_alts) {
    Rcpp::stop("Error: target_shares must have the same length as the total number of alternatives.");
  }
  if (arma::any(target_shares <= 0)) {
    Rcpp::stop("Error: all target_shares must be strictly positive (log(share) is undefined otherwise).");
  }
  validate_nl_inputs(X, alt_idx, nest_idx, M, use_asc, delta, &weights);

  arma::uvec alt_idx0 = alt_idx - 1;
  arma::uvec nest_idx0 = nest_idx - 1;
  const Rcpp::IntegerVector S = compute_prefix_sum(M);

  // The iteration bookkeeping (delta_old/delta_new, target/predicted log-shares,
  // residual) lives in the outside-inclusive share space of length num_alts:
  //   index 0 = outside option (when present), indices 1..J = inside alts.
  // But nl_predict_shares_internal indexes delta by INSIDE-alt index
  // (alt_idx0 in {0..J-1}) and expects a length-J inside-delta vector (the
  // outside option is handled separately via include_outside_option). We
  // therefore feed it delta_old.subvec(1, num_alts - 1) when an outside option
  // is present, and pin the outside slot delta_old[0] = 0 throughout (the
  // outside option's utility is the fixed normalization).
  arma::vec delta_old = arma::zeros(num_alts);
  if (include_outside_option) {
    delta_old.subvec(1, num_alts - 1) = delta;
    delta_old[0] = 0.0;
  } else {
    delta_old = delta;
    delta_old -= delta_old[0];
  }

  arma::vec inside_delta_old = include_outside_option
    ? arma::vec(delta_old.subvec(1, num_alts - 1))
    : delta_old;

  arma::vec log_shares_old = nl_predict_shares_internal(
    X, beta, lambda, alt_idx0, nest_idx0, M, S, weights, inside_delta_old,
    n_nests, num_alts, use_asc, include_outside_option
  );
  log_shares_old = arma::log(log_shares_old);
  arma::vec log_shares_target = arma::log(target_shares);
  arma::vec delta_new = delta_old;

  arma::vec delta_diff(num_alts, arma::fill::ones);
  double residual = 10.0;
  int iter = 0;

  while (iter < max_iter) {
    Rcpp::checkUserInterrupt(); // H4: allow user to interrupt long-running contraction
    delta_new = delta_old + damping * (log_shares_target - log_shares_old);
    if (include_outside_option) {
      // Outside option's delta is the fixed normalization; pin it at 0.
      delta_new[0] = 0.0;
    }
    delta_diff = arma::abs(delta_new - delta_old);
    residual = arma::max(delta_diff);
    if (residual < tol) {
      break;
    }
    delta_old = delta_new;
    inside_delta_old = include_outside_option
      ? arma::vec(delta_old.subvec(1, num_alts - 1))
      : delta_old;
    log_shares_old = nl_predict_shares_internal(
      X, beta, lambda, alt_idx0, nest_idx0, M, S, weights, inside_delta_old,
      n_nests, num_alts, use_asc, include_outside_option
    );
    log_shares_old = arma::log(log_shares_old);
    ++iter;
  }

  if (iter >= max_iter) {
    Rcpp::Rcout << "Warning: Maximum iterations reached without convergence." << std::endl;
  }

  delta_new -= delta_new[0];

  if (include_outside_option) {
    return delta_new.subvec(1, num_alts - 1);
  } else {
    return delta_new;
  }
}
