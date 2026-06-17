// [[Rcpp::depends(RcppArmadillo)]]
#include "choicer.h"
#include "choicer_internal.h"

//' Log-likelihood and gradient for multinomial logit model
//'
//' Computes the log-likelihood and its gradient for the Multinomial Logit model using OpenMP for parallelization.
//' Allows for inclusion of alternative-specific constants, outside option, and observation weights.
//'
//' @param theta K + J - 1 or K + J vector with model parameters
//' @param X sum(M) x K design matrix with covariates. Stacks M\[i] x K matrices for individual i.
//' @param alt_idx sum(M) x 1 vector with indices of alternatives within each choice set; 1-based indexing
//' @param choice_idx N x 1 vector with indices of chosen alternatives; 1-based indexing relative to X; 0 is used if include_outside_option=True
//' @param M N x 1 vector with number of alternatives for each individual
//' @param weights N x 1 vector with weights for each observation
//' @param use_asc whether to use alternative-specific constants
//' @param include_outside_option whether to include outside option normalized to 0 (if so, the outside option is not included in the data)
//' @returns List with loglikelihood and gradient evaluated at input arguments
//' @examples
//' \donttest{
//' library(data.table)
//' set.seed(42)
//' N <- 50; J <- 3
//' dt <- data.table(id = rep(1:N, each = J), alt = rep(1:J, N))
//' dt[, `:=`(x1 = rnorm(.N), x2 = rnorm(.N))]
//' dt[, choice := 0L]
//' dt[, choice := sample(c(1L, rep(0L, J - 1))), by = id]
//' d <- prepare_mnl_data(dt, "id", "alt", "choice", c("x1", "x2"))
//' theta <- rep(0, ncol(d$X) + nrow(d$alt_mapping) - 1)
//' result <- mnl_loglik_gradient_parallel(theta, d$X, d$alt_idx,
//'   d$choice_idx, d$M, d$weights)
//' result$objective  # negative log-likelihood
//' }
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

  const MnlParams par = parse_mnl_theta(theta, K, use_asc, include_outside_option);
  validate_choice_data(X, alt_idx, M, use_asc, par.delta, &weights, &choice_idx);

  // alt_idx is 1-based indexing => shift to 0-based indexing
  arma::uvec alt_idx0 = alt_idx - 1;

  // Compute prefix sums for indexing
  const Rcpp::IntegerVector S = compute_prefix_sum(M);

  // Pre-compute base utility for all individuals (single BLAS call)
  arma::vec base_util = compute_base_util(X, par.beta, alt_idx0, use_asc, par.delta);

  // --- H2: Serial pre-loop validation of chosen-alternative indices ---
  // Rcpp::stop() is only safe outside parallel regions.
  for (int i = 0; i < N; ++i) {
    int chosen = choice_idx[i];
    if (!include_outside_option) chosen -= 1;
    const int num_choices_i = include_outside_option ? M[i] + 1 : M[i];
    if (chosen < 0 || chosen >= num_choices_i) {
      Rcpp::stop("Invalid chosen alternative index for individual %d", i);
    }
  }

  // Prepare global accumulators
  double global_loglik = 0.0;
  arma::vec global_grad = arma::zeros(n_params);
  bool nonfinite_seen = false; // H1: flag set inside parallel, warning emitted after

#ifdef _OPENMP
#pragma omp parallel
#endif
  {
    // Thread-local accumulators
    double local_loglik = 0.0;
    arma::vec local_grad = arma::zeros(n_params);
    arma::vec diff_vec; // pre-allocated per-thread, resized per individual
    bool local_nonfinite = false;

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
      arma::vec V_i(num_choices);
      arma::vec inside_utils = base_util.subvec(start_idx, end_idx);
      fill_choice_utilities(V_i, inside_utils, num_choices, include_outside_option);

      // log-likelihood -------------------------------------------------------

      // Vector of choice probabilities
      arma::vec P_i;
      const double log_denom = stable_softmax(V_i, P_i);

      // Identify chosen alternative (validated serially above)
      int chosen_alt = choice_idx[i];
      if (!include_outside_option) {
        chosen_alt -= 1; // shift by 1 for inside-only indexing
      }

      // Probability of chosen alternative
      double V_choice = V_i(chosen_alt);
      double log_P_choice = V_choice - log_denom;
      if (!std::isfinite(log_P_choice)) {
        local_nonfinite = true; // H1: note without printing (R-API unsafe here)
        log_P_choice = -1e10;  // fallback to a large negative value
      }

      // Accumulate local weighted log-likelihood
      local_loglik += w_i * log_P_choice;

      // Gradient -------------------------------------------------------------
      // diff_vec[a] = 1{a == chosen_alt} - P_i[a]
      diff_vec.set_size(num_choices);
      diff_vec = -P_i;
      diff_vec(chosen_alt) += 1.0;

      // ---- Beta block: single BLAS dgemv ----
      if (include_outside_option) {
        const auto diff_inside = diff_vec.subvec(1, m_i); // subview, no copy; m_i elements
        local_grad.subvec(0, K - 1) += w_i * X_i.t() * diff_inside;
      } else {
        local_grad.subvec(0, K - 1) += w_i * X_i.t() * diff_vec;
      }

      // ---- Delta block: scatter (irregular alt-index mapping) ----
      if (use_asc) {
        if (include_outside_option) {
          scatter_delta_grad(local_grad, K, diff_vec.subvec(1, m_i), alt_idx0_i,
                             m_i, include_outside_option, w_i);
        } else {
          scatter_delta_grad(local_grad, K, diff_vec, alt_idx0_i,
                             m_i, include_outside_option, w_i);
        }
      }
    } // end of i loop

    // Combine partial accumulators into global accumulators
#ifdef _OPENMP
#pragma omp critical
#endif
    {
      global_loglik += local_loglik;
      global_grad += local_grad;
      if (local_nonfinite) nonfinite_seen = true;
    }
  } // end parallel region

  // H1: emit warning once, outside the parallel region (R-API safe here)
  if (nonfinite_seen) {
    Rcpp::warning("Non-finite log-probability encountered; clamped to -1e10. Check for extreme utility values.");
  }

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

//' Prediction of choice probabilities and utilities based on fitted model
//'
//' @param theta K + J - 1 or K + J vector with model parameters
//' @param X sum(M) x K design matrix with covariates. Stacks M\[i] x K matrices for individual i.
//' @param alt_idx sum(M) x 1 vector with indices of alternatives within each choice set; 1-based indexing
//' @param M N x 1 vector with number of alternatives for each individual
//' @param use_asc whether to use alternative-specific constants
//' @param include_outside_option whether to include outside option normalized to 0 (if so, the outside option is not included in the data)
//' @returns List with choice probability and utility for each choice situation evaluated at input arguments
//' @examples
//' \donttest{
//' library(data.table)
//' set.seed(42)
//' N <- 50; J <- 3
//' dt <- data.table(id = rep(1:N, each = J), alt = rep(1:J, N))
//' dt[, `:=`(x1 = rnorm(.N), x2 = rnorm(.N))]
//' dt[, choice := 0L]
//' dt[, choice := sample(c(1L, rep(0L, J - 1))), by = id]
//' fit <- run_mnlogit(dt, "id", "alt", "choice", c("x1", "x2"))
//' pred <- mnl_predict(coef(fit), fit$data$X, fit$data$alt_idx,
//'   fit$data$M, use_asc = TRUE)
//' head(pred$choice_prob)
//' }
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

  const MnlParams par = parse_mnl_theta(theta, K, use_asc, include_outside_option);
  validate_choice_data(X, alt_idx, M, use_asc, par.delta);

  // alt_idx is 1-based indexing => shift to 0-based indexing
  arma::uvec alt_idx0 = alt_idx - 1;

  // Compute prefix sums for indexing
  Rcpp::IntegerVector S = compute_prefix_sum(M);

  // Pre-compute base utility for all individuals (single BLAS call)
  arma::vec base_util = compute_base_util(X, par.beta, alt_idx0, use_asc, par.delta);

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

    // Build utility vector V_i
    arma::vec V_i(num_choices);
    arma::vec inside_utils = base_util.subvec(start_idx, end_idx);
    fill_choice_utilities(V_i, inside_utils, num_choices, include_outside_option);

    V_all.subvec(start_idx, end_idx) = inside_utils;

    // Vector of choice probabilities
    arma::vec P_i;
    stable_softmax(V_i, P_i);

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
  // Pre-compute base utility for all individuals (single BLAS call)
  arma::vec base_util = compute_base_util(X, beta, alt_idx0, use_asc, delta);

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
    arma::uvec alt_idx0_i = alt_idx0.subvec(start_idx, end_idx); // M[i]

    // Build utility vector V_i
    arma::vec V_i(num_choices);
    arma::vec inside_utils = base_util.subvec(start_idx, end_idx);
    fill_choice_utilities(V_i, inside_utils, num_choices, include_outside_option);

    // Vector of choice probabilities
    arma::vec P_i;
    stable_softmax(V_i, P_i);

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
//' @returns vector with predicted market shares for each alternative
//' @examples
//' \donttest{
//' library(data.table)
//' set.seed(42)
//' N <- 50; J <- 3
//' dt <- data.table(id = rep(1:N, each = J), alt = rep(1:J, N))
//' dt[, `:=`(x1 = rnorm(.N), x2 = rnorm(.N))]
//' dt[, choice := 0L]
//' dt[, choice := sample(c(1L, rep(0L, J - 1))), by = id]
//' fit <- run_mnlogit(dt, "id", "alt", "choice", c("x1", "x2"))
//' shares <- mnl_predict_shares(coef(fit), fit$data$X, fit$data$alt_idx,
//'   fit$data$M, fit$data$weights, use_asc = TRUE)
//' shares
//' }
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

  const MnlParams par = parse_mnl_theta(theta, K, use_asc, include_outside_option);
  validate_choice_data(X, alt_idx, M, use_asc, par.delta, &weights);

  // alt_idx is 1-based indexing => shift to 0-based indexing
  arma::uvec alt_idx0 = alt_idx - 1;

  // Compute prefix sums for indexing
  Rcpp::IntegerVector S = compute_prefix_sum(M);

  // total number of distinct alternatives
  int num_alts = include_outside_option ? (alt_idx.max() + 1) :  alt_idx.max();

  arma::vec global_shares = mnl_predict_shares_internal(
    X, par.beta, alt_idx0, M, S, weights, par.delta, num_alts, use_asc, include_outside_option
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
//' @returns vector with contraction's delta (ASCs) output
//' @examples
//' \donttest{
//' library(data.table)
//' set.seed(42)
//' N <- 50; J <- 3
//' dt <- data.table(id = rep(1:N, each = J), alt = rep(1:J, N))
//' dt[, `:=`(x1 = rnorm(.N), x2 = rnorm(.N))]
//' dt[, choice := 0L]
//' dt[, choice := sample(c(1L, rep(0L, J - 1))), by = id]
//' fit <- run_mnlogit(dt, "id", "alt", "choice", c("x1", "x2"))
//' beta <- coef(fit)[fit$param_map$beta]
//' delta <- blp_contraction(rep(0, J), rep(1/J, J), fit$data$X,
//'   beta, fit$data$alt_idx, fit$data$M, fit$data$weights)
//' delta
//' }
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
  if (arma::any(target_shares <= 0)) {
    Rcpp::stop("Error: all target_shares must be strictly positive (log(share) is undefined otherwise).");
  }
  validate_choice_data(X, alt_idx, M, use_asc, delta, &weights);

  // alt_idx is 1-based indexing => shift to 0-based indexing
  arma::uvec alt_idx0 = alt_idx - 1;

  // Compute prefix sums for indexing
  Rcpp::IntegerVector S = compute_prefix_sum(M);

  // The iteration bookkeeping (delta_old/delta_new, target/predicted log-shares,
  // residual) lives in the outside-inclusive share space of length num_alts:
  //   index 0 = outside option (when present), indices 1..J = inside alts.
  // But mnl_predict_shares_internal indexes delta by INSIDE-alt index
  // (alt_idx0 in {0..J-1}) and expects a length-J inside-delta vector (the
  // outside option is handled separately via include_outside_option). We
  // therefore feed it delta_old.subvec(1, num_alts - 1) when an outside option
  // is present, and pin the outside slot delta_old[0] = 0 throughout (the
  // outside option's utility is the fixed normalization).
  arma::vec delta_old = arma::zeros(num_alts);
  if (include_outside_option) {
    delta_old.subvec(1, num_alts - 1) = delta; // outside option at index 0
    delta_old[0] = 0.0;
  } else {
    delta_old = delta; // no outside option, delta already has J elements
    delta_old -= delta_old[0];
  }

  arma::vec inside_delta_old = include_outside_option
    ? arma::vec(delta_old.subvec(1, num_alts - 1))
    : delta_old;

  // Initialize shares and delta_new
  arma::vec log_shares_old = mnl_predict_shares_internal(
    X, beta, alt_idx0, M, S, weights, inside_delta_old, num_alts, use_asc, include_outside_option
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
      Rcpp::checkUserInterrupt(); // H4: allow user to interrupt long-running contraction
      delta_new = delta_old + (log_shares_target - log_shares_old);
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
      log_shares_old = mnl_predict_shares_internal(
        X, beta, alt_idx0, M, S, weights, inside_delta_old, num_alts, use_asc, include_outside_option
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
//' @returns Hessian matrix of the negative log-likelihood
//' @examples
//' \donttest{
//' library(data.table)
//' set.seed(42)
//' N <- 50; J <- 3
//' dt <- data.table(id = rep(1:N, each = J), alt = rep(1:J, N))
//' dt[, `:=`(x1 = rnorm(.N), x2 = rnorm(.N))]
//' dt[, choice := 0L]
//' dt[, choice := sample(c(1L, rep(0L, J - 1))), by = id]
//' fit <- run_mnlogit(dt, "id", "alt", "choice", c("x1", "x2"))
//' H <- mnl_loglik_hessian_parallel(coef(fit), fit$data$X, fit$data$alt_idx,
//'   fit$data$choice_idx, fit$data$M, fit$data$weights)
//' dim(H)
//' }
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
  const MnlParams par = parse_mnl_theta(theta, K, use_asc, include_outside_option);
  validate_choice_data(X, alt_idx, M, use_asc, par.delta, &weights);

  // alt_idx is 1-based indexing => shift to 0-based indexing
  arma::uvec alt_idx0 = alt_idx - 1;

  // Compute prefix sums for indexing each individual's block in X / alt_idx
  const Rcpp::IntegerVector S = compute_prefix_sum(M);

  // Pre-compute base utility for all individuals (single BLAS call)
  arma::vec base_util = compute_base_util(X, par.beta, alt_idx0, use_asc, par.delta);

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
      arma::vec V_i(num_choices);
      arma::vec inside_utils = base_util.subvec(start_idx, end_idx);
      fill_choice_utilities(V_i, inside_utils, num_choices, include_outside_option);

      // Calculate Probabilities ---------------------------------------------
      arma::vec P_i;
      stable_softmax(V_i, P_i);

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

//' Compute aggregate elasticities for MNL model
//'
//' Computes the aggregate elasticity matrix (weighted average of individual
//' elasticities) for the Multinomial Logit model.
//'
//' @param theta K + J - 1 or K + J vector with model parameters
//' @param X sum(M) x K design matrix with covariates.
//' @param alt_idx sum(M) x 1 vector with indices of alternatives; 1-based indexing
//' @param choice_idx N x 1 vector (kept for API consistency, but not used)
//' @param M N x 1 vector with number of alternatives for each individual
//' @param weights N x 1 vector with weights for each observation
//' @param elast_var_idx 1-based index of the column in X for which to compute the elasticity
//' @param use_asc whether to use alternative-specific constants
//' @param include_outside_option whether to include outside option
//' @returns J x J matrix of aggregate elasticities
//' @examples
//' \donttest{
//' library(data.table)
//' set.seed(42)
//' N <- 50; J <- 3
//' dt <- data.table(id = rep(1:N, each = J), alt = rep(1:J, N))
//' dt[, `:=`(x1 = rnorm(.N), x2 = rnorm(.N))]
//' dt[, choice := 0L]
//' dt[, choice := sample(c(1L, rep(0L, J - 1))), by = id]
//' fit <- run_mnlogit(dt, "id", "alt", "choice", c("x1", "x2"))
//' elas <- mnl_elasticities_parallel(coef(fit), fit$data$X, fit$data$alt_idx,
//'   fit$data$choice_idx, fit$data$M, fit$data$weights, elast_var_idx = 1L)
//' elas
//' }
//' @export
// [[Rcpp::export]]
arma::mat mnl_elasticities_parallel(
    const arma::vec& theta,
    const arma::mat& X,
    const arma::uvec& alt_idx,
    const arma::uvec& choice_idx, // Kept for consistency, but not used
    const Rcpp::IntegerVector& M,
    const arma::vec& weights,
    const int elast_var_idx,
    const bool use_asc = true,
    const bool include_outside_option = false
) {
  // --- 1. Parameter and Variable Setup ---
  const int N = M.size();
  const int K = X.n_cols;

  // Convert 1-based R index to 0-based C++ index
  const int var_idx = elast_var_idx - 1;
  if (var_idx < 0 || var_idx >= K) {
    Rcpp::stop("elast_var_idx is out of bounds for design matrix X.");
  }

  // Extract beta and delta (same logic as loglik function)
  const MnlParams par = parse_mnl_theta(theta, K, use_asc, include_outside_option);
  validate_choice_data(X, alt_idx, M, use_asc, par.delta, &weights);
  const double beta_k = par.beta(var_idx); // The coefficient for our variable

  // alt_idx is 1-based indexing => shift to 0-based indexing
  arma::uvec alt_idx0 = alt_idx - 1;

  // Determine total number of alternatives for the output matrix
  // J_inside = number of alternatives (with or without normalization)
  // J_total = total size of matrix (J_inside + 1 if outside option)
  const int J_inside = compute_J_inside(use_asc, par.delta, alt_idx0);
  const int J_total = compute_J_total(J_inside, include_outside_option);

  // Compute prefix sums for indexing
  const Rcpp::IntegerVector S = compute_prefix_sum(M);

  // Pre-compute base utility for all individuals (single BLAS call)
  arma::vec base_util = compute_base_util(X, par.beta, alt_idx0, use_asc, par.delta);

  // Prepare global accumulators
  arma::mat global_elas_matrix = arma::zeros(J_total, J_total);
  double global_total_weight = 0.0;

#ifdef _OPENMP
#pragma omp parallel
#endif
  {
    // Thread-local accumulators
    arma::mat local_elas_matrix = arma::zeros(J_total, J_total);
    double local_total_weight = 0.0;

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

      // --- 2. Compute Probabilities (same as loglik) ---
      arma::vec V_i(num_choices);
      arma::vec inside_utils = base_util.subvec(start_idx, end_idx);
      fill_choice_utilities(V_i, inside_utils, num_choices, include_outside_option);

      arma::vec P_i;
      stable_softmax(V_i, P_i);

      // --- 3. Prepare Data for Elasticity Calculation ---

      // Get the values of the variable k for this individual's choice set
      arma::vec x_k_i = arma::zeros(num_choices);
      if (include_outside_option) {
        // x_k for outside option (index 0) is 0
        x_k_i.subvec(1, num_choices - 1) = X_i.col(var_idx);
      } else {
        x_k_i = X_i.col(var_idx);
      }

      // Map local choice set indices (0...num_choices-1) to global
      // alternative indices (0...J_total-1)
      arma::uvec global_j_map =
          build_global_alt_map(alt_idx0_i, m_i, include_outside_option);

      // --- 4. Compute Individual Elasticity Matrix ---
      // E_n(i, j) = elasticity of P_ni w.r.t attribute x_njk
      for (int i_local = 0; i_local < num_choices; ++i_local) {
        const int global_i = global_j_map[i_local]; // Global row index
        const double P_ni = P_i[i_local];

        for (int j_local = 0; j_local < num_choices; ++j_local) {
          const int global_j = global_j_map[j_local]; // Global col index
          const double P_nj = P_i[j_local];
          const double x_njk = x_k_i[j_local];
          
          double elasticity;
          if (global_i == global_j) {
            // Own-elasticity: beta_k * x_nik * (1 - P_ni)
            elasticity = beta_k * x_njk * (1.0 - P_ni);
          } else {
            // Cross-elasticity: -beta_k * x_njk * P_nj
            elasticity = -beta_k * x_njk * P_nj;
          }

          // Accumulate the weighted elasticity
          local_elas_matrix(global_i, global_j) += w_i * elasticity;
        }
      }
      local_total_weight += w_i;
    } // end of i loop

    // Combine partial accumulators into global accumulators
#ifdef _OPENMP
#pragma omp critical
#endif
    {
      global_elas_matrix += local_elas_matrix;
      global_total_weight += local_total_weight;
    }
  } // end parallel region

  // --- 5. Finalize ---
  // Compute the average by dividing by the total weight
  if (global_total_weight > 1e-10) {
    global_elas_matrix /= global_total_weight;
  }

  return global_elas_matrix;
}


// =============================================================================
// Diversion Ratios
// =============================================================================

//' Compute MNL diversion ratios (parallelized over individuals)
//'
//' Computes the diversion ratio matrix DR(j->k), which measures the fraction
//' of demand lost by alternative j that is captured by alternative k.
//' For MNL: DR(j->k) = sum_n(w_n * P_nj * P_nk) / sum_n(w_n * P_nj * (1 - P_nj))
//'
//' @param theta K + J - 1 or K + J vector with model parameters
//' @param X sum(M) x K design matrix with covariates.
//' @param alt_idx sum(M) x 1 vector with indices of alternatives; 1-based indexing
//' @param M N x 1 vector with number of alternatives for each individual
//' @param weights N x 1 vector with weights for each observation
//' @param use_asc whether to use alternative-specific constants
//' @param include_outside_option whether to include outside option
//' @returns J x J matrix where entry (k, j) = DR(j->k). Diagonal is 0.
//' @examples
//' \donttest{
//' library(data.table)
//' set.seed(42)
//' N <- 50; J <- 3
//' dt <- data.table(id = rep(1:N, each = J), alt = rep(1:J, N))
//' dt[, `:=`(x1 = rnorm(.N), x2 = rnorm(.N))]
//' dt[, choice := 0L]
//' dt[, choice := sample(c(1L, rep(0L, J - 1))), by = id]
//' fit <- run_mnlogit(dt, "id", "alt", "choice", c("x1", "x2"))
//' dr <- mnl_diversion_ratios_parallel(coef(fit), fit$data$X, fit$data$alt_idx,
//'   fit$data$M, fit$data$weights)
//' dr
//' }
//' @export
// [[Rcpp::export]]
arma::mat mnl_diversion_ratios_parallel(
    const arma::vec& theta,
    const arma::mat& X,
    const arma::uvec& alt_idx,
    const Rcpp::IntegerVector& M,
    const arma::vec& weights,
    const bool use_asc = true,
    const bool include_outside_option = false
) {
  // --- 1. Parameter and Variable Setup ---
  const int N = M.size();
  const int K = X.n_cols;

  // Extract beta and delta (same logic as elasticities/loglik)
  const MnlParams par = parse_mnl_theta(theta, K, use_asc, include_outside_option);
  validate_choice_data(X, alt_idx, M, use_asc, par.delta, &weights);

  // alt_idx is 1-based indexing => shift to 0-based indexing
  arma::uvec alt_idx0 = alt_idx - 1;

  // Determine total number of alternatives for the output matrix
  const int J_inside = compute_J_inside(use_asc, par.delta, alt_idx0);
  const int J_total = compute_J_total(J_inside, include_outside_option);

  // Compute prefix sums for indexing
  const Rcpp::IntegerVector S = compute_prefix_sum(M);

  // Pre-compute base utility for all individuals (single BLAS call)
  arma::vec base_util = compute_base_util(X, par.beta, alt_idx0, use_asc, par.delta);

  // Prepare global accumulators
  // numerator(k, j) = sum_n w_n * P_nj * P_nk  (for k != j)
  // denominator(j) = sum_n w_n * P_nj * (1 - P_nj)
  arma::mat global_numerator = arma::zeros(J_total, J_total);
  arma::vec global_denominator = arma::zeros(J_total);

#ifdef _OPENMP
#pragma omp parallel
#endif
  {
    // Thread-local accumulators
    arma::mat local_numerator = arma::zeros(J_total, J_total);
    arma::vec local_denominator = arma::zeros(J_total);

#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
    for (int i = 0; i < N; ++i) {
      const int m_i         = M[i];
      const int num_choices = include_outside_option ? m_i + 1 : m_i;
      const int start_idx   = S[i];
      const int end_idx     = start_idx + m_i - 1;
      const double w_i      = weights[i];

      arma::uvec alt_idx0_i = alt_idx0.subvec(start_idx, end_idx); // M[i]

      // --- 2. Compute Probabilities ---
      arma::vec V_i(num_choices);
      arma::vec inside_utils = base_util.subvec(start_idx, end_idx);
      fill_choice_utilities(V_i, inside_utils, num_choices, include_outside_option);

      arma::vec P_i;
      stable_softmax(V_i, P_i);

      // --- 3. Map local indices to global alternative indices ---
      arma::uvec global_j_map =
          build_global_alt_map(alt_idx0_i, m_i, include_outside_option);

      // --- 4. Accumulate numerator and denominator ---
      for (int j_local = 0; j_local < num_choices; ++j_local) {
        const int global_j = global_j_map[j_local];
        const double P_nj = P_i[j_local];

        // Denominator: w_i * P_nj * (1 - P_nj)
        local_denominator(global_j) += w_i * P_nj * (1.0 - P_nj);

        // Numerator: w_i * P_nj * P_nk for k != j
        for (int k_local = 0; k_local < num_choices; ++k_local) {
          if (k_local == j_local) continue;
          const int global_k = global_j_map[k_local];
          const double P_nk = P_i[k_local];
          local_numerator(global_k, global_j) += w_i * P_nj * P_nk;
        }
      }
    } // end of i loop

    // Combine partial accumulators into global accumulators
#ifdef _OPENMP
#pragma omp critical
#endif
    {
      global_numerator += local_numerator;
      global_denominator += local_denominator;
    }
  } // end parallel region

  // --- 5. Compute diversion ratios ---
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