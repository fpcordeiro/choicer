// [[Rcpp::depends(RcppArmadillo)]]
#include "choicer.h"
#include "choicer_internal.h"
#include "choicer_derivs.h"

// Memory of the MNL kernels' loop over the choice situations, an upper bound
// per thread: the utilities, probabilities and score differences (V_i,
// inside_utils, P_i, diff_vec; max_m + 1 doubles each), the copy of X_i that
// Armadillo makes for X_i.t() * v (max_m x K) and its product (K), and n_vecs
// vectors of the n parameters (the thread's sums or scores).
inline double mnl_thread_bytes(const ChoiceLayout& lay, const int K,
                               const int n_params, const double n_vecs) {
  const double m = static_cast<double>(lay.max_m);
  return sizeof(double) * (4.0 * (m + 1.0) + m * K + K +
                           n_vecs * static_cast<double>(n_params));
}

// ============================================================================
// Compact blocks of the MNL derivative kernels
//
// A situation's score and its Hessian live on its compact block: the K
// coefficients, then the free constants of the alternatives it offers, in
// increasing global order. The normalized reference and the outside option
// have none; a situation that lists an alternative twice adds both rows into
// one slot, in row order, as scatter_delta_grad() adds them into one
// constant. MnlSituationMap holds the slots, one per thread, sized on the
// primary thread for the largest situation, so that building them allocates
// nothing; slot_of is -1 except while a situation is mapped.
// ============================================================================
struct MnlSituationMap {
  std::vector<int> slot_of;   // free constant j -> slot, or -1 (J_d entries)
  std::vector<int> delta;     // slot -> free constant, increasing
  std::vector<int> row_slot;  // row of the situation -> slot, or -1
  int a = 0;                  // free constants of the situation
  char pad[128];              // keeps the next thread's map off these lines

  // Each array is reserved 32 entries (128 bytes) beyond its largest size,
  // so that no two threads' arrays share a cache line.
  MnlSituationMap(const int J_d, const choicer_off max_m)
      : slot_of(static_cast<std::size_t>(J_d) + 32, -1) {
    delta.reserve(static_cast<std::size_t>(std::min<choicer_off>(J_d, max_m)) + 32);
    row_slot.reserve(static_cast<std::size_t>(max_m) + 32);
  }
  // Global parameter of local index i (the K coefficients first).
  int global(const int i, const int K) const {
    return i < K ? i : K + delta[i - K];
  }
  // Local index of global parameter g, or -1 outside the situation's block.
  int local(const int g, const int K) const {
    if (g < K) return g;
    const int k = slot_of[g - K];
    return k >= 0 ? K + k : -1;
  }
};

// Bytes of one thread's MnlSituationMap.
inline double mnl_map_bytes(const int J_d, const choicer_off max_m) {
  const double m = static_cast<double>(max_m);
  return sizeof(int) * (J_d + std::min<double>(J_d, m) + m + 96.0) + 256.0;
}

// Every thread's map, each in its own allocation, on the primary thread.
using MnlSituationMaps = std::vector<std::unique_ptr<MnlSituationMap>>;
inline MnlSituationMaps mnl_situation_maps(const int J_d, const choicer_off max_m,
                                           const int n_threads) {
  MnlSituationMaps maps;
  maps.reserve(n_threads);
  for (int i = 0; i < n_threads; ++i) {
    maps.emplace_back(new MnlSituationMap(J_d, max_m));
  }
  return maps;
}

// Free constant of a row with 0-based alternative alt0 (as
// scatter_delta_grad() places it), or -1.
inline int mnl_free_delta(const int alt0, const bool use_asc,
                          const bool include_outside_option) {
  if (!use_asc) return -1;
  return include_outside_option ? alt0 : alt0 - 1;
}

// The slots of situation i: its distinct free constants, in increasing
// order, and the slot of each of its rows. mnl_situation_map_clear() resets
// the markers.
inline void mnl_situation_map(const ChoiceLayout& lay, const choicer_off i,
                              const bool use_asc,
                              const bool include_outside_option,
                              MnlSituationMap& mp) {
  const int* alt = lay.alt + lay.row_off[i];
  const int m = lay.m(i);
  mp.delta.clear();
  for (int r = 0; r < m; ++r) {
    const int j = mnl_free_delta(alt[r] - 1, use_asc, include_outside_option);
    if (j >= 0 && mp.slot_of[j] < 0) {
      mp.slot_of[j] = 0;  // seen
      mp.delta.push_back(j);
    }
  }
  std::sort(mp.delta.begin(), mp.delta.end());
  mp.a = static_cast<int>(mp.delta.size());
  for (int k = 0; k < mp.a; ++k) mp.slot_of[mp.delta[k]] = k;
  mp.row_slot.resize(static_cast<std::size_t>(m));
  for (int r = 0; r < m; ++r) {
    const int j = mnl_free_delta(alt[r] - 1, use_asc, include_outside_option);
    mp.row_slot[r] = j >= 0 ? mp.slot_of[j] : -1;
  }
}

inline void mnl_situation_map_clear(MnlSituationMap& mp) {
  for (int k = 0; k < mp.a; ++k) mp.slot_of[mp.delta[k]] = -1;
}

// The weight-free score of situation i on its compact block, into s_c (K + a
// entries, the map built for i), and the situation's 1{chosen} - P into
// diff_vec: the coefficients' entries as X_i' diff, the expression of the
// gradient kernel's OPG diagonal, then each slot's differences in row order
// at scale 1.0, as scatter_delta_grad() adds them into the n_params-vector.
inline void mnl_situation_score(const arma::mat& X, const ChoiceLayout& lay,
                                const arma::vec& base_util, const choicer_off i,
                                const bool include_outside_option,
                                const MnlSituationMap& mp, arma::vec& diff_vec,
                                arma::vec& s_c) {
  const int m_i         = lay.m(i);
  const int num_choices = include_outside_option ? m_i + 1 : m_i;
  const arma::uword start_idx = static_cast<arma::uword>(lay.row_off[i]);
  const arma::uword end_idx   = start_idx + m_i - 1;
  const auto X_i        = X.rows(start_idx, end_idx); // M[i] x K
  const int K = X.n_cols;
  const int o = include_outside_option ? 1 : 0;  // slot of the first inside alt

  // Build utility vector V_i
  arma::vec V_i(num_choices);
  arma::vec inside_utils = base_util.subvec(start_idx, end_idx);
  fill_choice_utilities(V_i, inside_utils, num_choices, include_outside_option);

  // Probabilities
  arma::vec P_i;
  stable_softmax(V_i, P_i);

  // Chosen alternative (validated serially by the caller)
  int chosen_alt = lay.choice[i];
  if (!include_outside_option) chosen_alt -= 1;

  // diff_vec[a] = 1{a == chosen_alt} - P_i[a]  (same as gradient kernel)
  diff_vec.set_size(num_choices);
  diff_vec = -P_i;
  diff_vec(chosen_alt) += 1.0;

  s_c.zeros();
  if (include_outside_option) {
    s_c.subvec(0, K - 1) = X_i.t() * diff_vec.subvec(1, m_i);
  } else {
    s_c.subvec(0, K - 1) = X_i.t() * diff_vec;
  }
  for (int r = 0; r < m_i; ++r) {
    const int k = mp.row_slot[r];
    if (k >= 0) s_c[K + k] += 1.0 * diff_vec[r + o];
  }
}

// The MNL gradient kernel's pass over the choice situations: the weighted
// sums of their log-probabilities and scores and, with OpgDiag, the diagonal
// of the BHHH matrix, sum_i w_i s_i o s_i, its terms formed as
// mnl_bhhh_parallel() forms them (the beta block from the unweighted product
// X_i' diff, the constants scattered at scale 1.0 and squared once each) and
// rounded before they are added (OpgDiagTerms). Without OpgDiag it is the
// kernel's loop as it was.
template <bool OpgDiag>
static void mnl_gradient_pass(const arma::mat& X, const ChoiceLayout& lay,
                              const arma::vec& base_util,
                              const arma::vec& weights, const int K,
                              const int n_params, const bool use_asc,
                              const bool include_outside_option,
                              double& global_loglik, arma::vec& global_grad,
                              arma::vec& global_diag, bool& nonfinite_seen,
                              const int n_threads) {
  const choicer_off N = lay.N;

#ifdef _OPENMP
#pragma omp parallel num_threads(n_threads)
#endif
  {
    // Thread-local accumulators
    double local_loglik = 0.0;
    arma::vec local_grad = arma::zeros(n_params);
    arma::vec diff_vec; // pre-allocated per-thread, resized per individual
    bool local_nonfinite = false;
    // The diagonal, the constants' scores (zero between situations), the
    // beta block's scores and a situation's terms
    arma::vec local_diag, s_d, s_b;
    OpgDiagTerms terms(OpgDiag ? K + lay.max_m : 0);
    if (OpgDiag) {
      local_diag.zeros(n_params);
      if (use_asc) s_d.zeros(n_params);
    }

#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
    for (choicer_off i = 0; i < N; ++i) {
      const int m_i         = lay.m(i);
      const int num_choices = include_outside_option ? m_i + 1 : m_i;
      const arma::uword start_idx = static_cast<arma::uword>(lay.row_off[i]);
      const arma::uword end_idx   = start_idx + m_i - 1;
      const double w_i      = weights[i];
      const auto X_i        = X.rows(start_idx, end_idx); // M[i] x K
      const AltCodes0 alt_idx0_i{lay.alt + lay.row_off[i]}; // M[i]

      // Build utility vector V_i
      arma::vec V_i(num_choices);
      arma::vec inside_utils = base_util.subvec(start_idx, end_idx);
      fill_choice_utilities(V_i, inside_utils, num_choices, include_outside_option);

      // log-likelihood -------------------------------------------------------

      // Vector of choice probabilities
      arma::vec P_i;
      const double log_denom = stable_softmax(V_i, P_i);

      // Identify chosen alternative (validated serially above)
      int chosen_alt = lay.choice[i];
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

      // ---- OPG diagonal (opg_diag): the weight-free score's squares ----
      if (OpgDiag) {
        if (include_outside_option) {
          const auto diff_inside = diff_vec.subvec(1, m_i);
          s_b = X_i.t() * diff_inside;
          if (use_asc) {
            scatter_delta_grad(s_d, K, diff_inside, alt_idx0_i, m_i,
                               include_outside_option, 1.0);
          }
        } else {
          s_b = X_i.t() * diff_vec;
          if (use_asc) {
            scatter_delta_grad(s_d, K, diff_vec, alt_idx0_i, m_i,
                               include_outside_option, 1.0);
          }
        }
        for (int k = 0; k < K; ++k) terms.push(k, s_b[k], w_i);
        if (use_asc) {
          opg_diag_push_delta(terms, s_d, K, alt_idx0_i, m_i,
                              include_outside_option, w_i);
        }
        terms.flush(local_diag);
      }
    } // end of i loop

    // Combine partial accumulators into global accumulators
#ifdef _OPENMP
#pragma omp critical
#endif
    {
      global_loglik += local_loglik;
      global_grad += local_grad;
      if (OpgDiag) global_diag += local_diag;
      if (local_nonfinite) nonfinite_seen = true;
    }
  } // end parallel region
}

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
//' @param opg_diag Logical; \code{TRUE} also returns the diagonal of the
//'   BHHH (outer product of gradients) matrix with the same weights,
//'   \eqn{\sum_i w_i s_i \circ s_i}, accumulated in the same pass.
//' @returns List with loglikelihood and gradient evaluated at input
//'   arguments; with \code{opg_diag}, also \code{opg_diag}, as computed
//'   (not sanitized)
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
//' result <- choicer:::mnl_loglik_gradient_parallel(theta, d$X, d$alt_idx,
//'   d$choice_idx, d$M, d$weights)
//' result$objective  # negative log-likelihood
//' }
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List mnl_loglik_gradient_parallel(
    const arma::vec& theta,
    const arma::mat& X,
    const Rcpp::IntegerVector& alt_idx,
    const Rcpp::IntegerVector& choice_idx,
    const Rcpp::IntegerVector& M,
    const arma::vec& weights,
    const bool use_asc = true,
    const bool include_outside_option = false,
    const bool opg_diag = false
) {
  // Extract beta and delta from theta
  const int K = X.n_cols;
  const int n_params = theta.n_elem;

  const MnlParams par = parse_mnl_theta(theta, K, use_asc, include_outside_option);
  // Layout of the stacked design: 64-bit row offsets, and the alternative
  // codes and choices read in place
  const ChoiceLayout lay = choice_layout_build(X, alt_idx, M, use_asc,
                                               par.delta, &weights, &choice_idx);

  // --- H2: Serial pre-loop validation of chosen-alternative indices ---
  // Rcpp::stop() is only safe outside parallel regions.
  validate_choices(lay, include_outside_option, nullptr);

  // Memory: each thread's sums and situation buffers (with opg_diag also its
  // diagonal, the constants' scores and the terms), X beta, the sums, the
  // negated gradient and R's copies of it and of the diagonal. Each thread
  // allocates its sums before its first situation, so the pass runs on the
  // team checked here.
  const int n_threads = kernel_team_threads(lay.N);
  kernel_check_memory(
      "gradient", n_threads,
      mnl_thread_bytes(lay, K, n_params, opg_diag ? 3.0 : 1.0) +
          (opg_diag ? 16.0 * (K + static_cast<double>(lay.max_m)) + 8.0 * K : 0.0),
      [&](const int) {
        return 8.0 * static_cast<double>(lay.n_rows) +
               (opg_diag ? 40.0 : 24.0) * n_params;
      });

  // Pre-compute base utility for all individuals (single BLAS call)
  arma::vec base_util = compute_base_util(X, par.beta, lay, use_asc, par.delta);

  // Prepare global accumulators
  double global_loglik = 0.0;
  arma::vec global_grad = arma::zeros(n_params);
  arma::vec global_diag;
  bool nonfinite_seen = false; // H1: flag set inside parallel, warning emitted after
  if (opg_diag) {
    global_diag.zeros(n_params);
    mnl_gradient_pass<true>(X, lay, base_util, weights, K, n_params, use_asc,
                            include_outside_option, global_loglik, global_grad,
                            global_diag, nonfinite_seen, n_threads);
  } else {
    mnl_gradient_pass<false>(X, lay, base_util, weights, K, n_params, use_asc,
                             include_outside_option, global_loglik, global_grad,
                             global_diag, nonfinite_seen, n_threads);
  }

  // H1: emit warning once, outside the parallel region (R-API safe here)
  if (nonfinite_seen) {
    Rcpp::warning("Non-finite log-probability encountered; clamped to -1e10. Check for extreme utility values.");
  }

  // H3: Sanitize NaN/Inf (matching mxl_loglik_gradient_parallel); the
  // diagonal is returned as computed
  double obj = -global_loglik;
  arma::vec grad = -global_grad;
  if (!std::isfinite(obj)) {
    obj = 1e10;
    grad.zeros();
  } else {
    grad.elem(arma::find_nonfinite(grad)).zeros();
  }
  if (opg_diag) {
    return Rcpp::List::create(
      Rcpp::Named("objective") = obj,
      Rcpp::Named("gradient")  = grad,
      Rcpp::Named("opg_diag")  = global_diag
    );
  }
  return Rcpp::List::create(
    Rcpp::Named("objective") = obj,
    Rcpp::Named("gradient")  = grad
  );
}

//' BHHH/OPG information matrix for multinomial logit model
//'
//' Computes the weighted outer product of per-situation scores
//' \eqn{\sum_i w_i\, s_i s_i^\top} for the Multinomial Logit model. The
//' per-situation score \eqn{s_i} is the (positive) gradient of situation
//' \eqn{i}'s log-likelihood contribution and is weight-free; the supplied
//' \code{weights} enter only as the leading multiplier. Passing
//' \code{weights = w} yields the ordinary weighted BHHH/OPG information; passing
//' \code{weights = w^2} yields the sandwich \emph{meat}
//' \eqn{B = \sum_i w_i^2 s_i s_i^\top} used for robust (WESML) inference.
//' Each situation's outer product is formed only over the parameters its
//' likelihood involves (the coefficients and the constants of the
//' alternatives it offers), so the work per situation does not grow with the
//' number of alternatives; the n_params x n_params result and the threads'
//' accumulators do (see acc_bytes).
//'
//' @param theta K + J - 1 or K + J vector with model parameters
//' @param X sum(M) x K design matrix with covariates. Stacks M\[i] x K matrices for individual i.
//' @param alt_idx sum(M) x 1 vector with indices of alternatives within each choice set; 1-based indexing
//' @param choice_idx N x 1 vector with indices of chosen alternatives; 1-based indexing relative to X; 0 is used if include_outside_option=True
//' @param M N x 1 vector with number of alternatives for each individual
//' @param weights N x 1 vector with weights for each observation
//' @param use_asc whether to use alternative-specific constants
//' @param include_outside_option whether to include outside option normalized to 0 (if so, the outside option is not included in the data)
//' @param acc_bytes Bytes for the threads' private upper triangles (2 GiB by
//'   default). While the triangles of the T threads in use, T n_params
//'   (n_params + 1) / 2 doubles, fit, each thread adds into its own; above,
//'   or for a NaN or negative value, the threads add into the result (the
//'   coefficients' rows through smaller private buffers, the constants' block
//'   atomically). One thread adds into the result whatever the value. The two
//'   ways differ by rounding only; \code{0} forces the shared result (for
//'   tests).
//' @returns A symmetric positive-semidefinite information matrix
//'   \eqn{\sum_i w_i\, s_i s_i^\top} (same sign convention as the negated Hessian).
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
//' B <- choicer:::mnl_bhhh_parallel(coef(fit), fit$data$X, fit$data$alt_idx,
//'   fit$data$choice_idx, fit$data$M, fit$data$weights)
//' dim(B)
//' }
//' @keywords internal
// [[Rcpp::export]]
Rcpp::NumericMatrix mnl_bhhh_parallel(
    const arma::vec& theta,
    const arma::mat& X,
    const Rcpp::IntegerVector& alt_idx,
    const Rcpp::IntegerVector& choice_idx,
    const Rcpp::IntegerVector& M,
    const arma::vec& weights,
    const bool use_asc = true,
    const bool include_outside_option = false,
    const double acc_bytes = 2147483648.0
) {
  const int K = X.n_cols;
  const int n_params = theta.n_elem;

  const MnlParams par = parse_mnl_theta(theta, K, use_asc, include_outside_option);
  const ChoiceLayout lay = choice_layout_build(X, alt_idx, M, use_asc,
                                               par.delta, &weights, &choice_idx);
  const choicer_off N = lay.N;

  // --- Serial pre-loop validation of chosen-alternative indices ---
  validate_choices(lay, include_outside_option, "mnl_bhhh_parallel");

  // Memory: each thread's slot map, compact score, a column of the update
  // (every parameter, for the dense fallback) and situation buffers; X beta,
  // the result and the threads' accumulators
  const int J_d = n_params - K;  // free constants
  const double a_max = std::min<double>(J_d, static_cast<double>(lay.max_m));
  const int n_threads = kernel_team_threads(N);
  kernel_check_memory("BHHH matrix", n_threads,
                      mnl_map_bytes(J_d, lay.max_m) + 8.0 * (K + a_max + n_params + 32.0) +
                          mnl_thread_bytes(lay, K, n_params, 0.0),
                      [&](const int t) {
                        return 8.0 * static_cast<double>(lay.n_rows) +
                               sym_acc_need(n_params, K, t, acc_bytes);
                      });

  // Pre-compute base utility for all individuals (single BLAS call)
  arma::vec base_util = compute_base_util(X, par.beta, lay, use_asc, par.delta);

  // The result and the threads' accumulators, slot maps and buffers, all
  // allocated here on the primary thread
  SymAcc acc(n_params, n_threads, K, "BHHH", acc_bytes);
  MnlSituationMaps maps = mnl_situation_maps(J_d, lay.max_m, n_threads);
  std::vector<std::vector<double>> scores =
      thread_vectors<double>(static_cast<std::size_t>(K + a_max), n_threads);
  std::vector<std::vector<double>> cols =
      thread_vectors<double>(static_cast<std::size_t>(n_params), n_threads);

#ifdef _OPENMP
#pragma omp parallel num_threads(n_threads)
#endif
  {
    const int tid = kernel_thread_num();
    MnlSituationMap& mp = *maps[tid];
    double* acc_t = acc.base(tid);
    double* s_buf = scores[tid].data();
    double* col = cols[tid].data();
    arma::vec diff_vec; // resized per individual

#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
    for (choicer_off i = 0; i < N; ++i) {
      const double w_i = weights[i];
      mnl_situation_map(lay, i, use_asc, include_outside_option, mp);
      arma::vec s_c(s_buf, K + mp.a, false, true);
      mnl_situation_score(X, lay, base_util, i, include_outside_option, mp,
                          diff_vec, s_c);

      // w_i s_c s_c' over the compact block, column by column, each entry
      // formed as the dense rank-one update formed it (`B += w_i * s_i *
      // s_i.t()` is Armadillo's syrk_vec: the product of the two scores,
      // times the weight, then added). A column's terms are stored before
      // they are added, so that the weight's product is never fused into the
      // addition (bit for bit under clang and under GCC at R's default -O2;
      // GCC's -O3 loop vectorizer rounds the dense update differently). The
      // dense update added w_i * (0 * s) outside the block, +-0, which
      // changes nothing: no accumulator, which starts at +0, is ever -0.
      if (std::isfinite(w_i) && s_c.is_finite()) {
        if (acc.shared) {
          bhhh_add<true>(acc, acc_t, mp, K, s_c, w_i, col);
        } else {
          bhhh_add<false>(acc, acc_t, mp, K, s_c, w_i, col);
        }
      } else {
        // A non-finite weight, or a non-finite score (utilities that
        // overflowed): the dense update over every entry, whose 0 * NaN and
        // 0 * Inf spread it over whole rows and columns. The scores come
        // from the slot map, zero outside the block.
        if (acc.shared) {
          bhhh_add_dense<true>(acc, acc_t, mp, K, s_c, w_i, col);
        } else {
          bhhh_add_dense<false>(acc, acc_t, mp, K, s_c, w_i, col);
        }
      }
      mnl_situation_map_clear(mp);
    } // end of i loop
  } // end parallel region

  // Return PSD information matrix (same sign convention as negated Hessian).
  acc.finish(false);
  return acc.result;
}

// Cluster-robust meat for the multinomial logit model (internal).
//
// Returns the n_params x n_params matrix sum_g s_g s_g' with
// s_g = sum_{i in g} w_i s_i the weighted score of cluster g, s_i the
// weight-free score of choice situation i and w_i its weight: what
// crossprod(rowsum(w * S, cluster)) forms from the N x n_params score matrix
// of mnl_scores_parallel(), without that matrix. `cluster` holds one code per
// choice situation, 1..G (R passes match(as.character(labels), unique(...)),
// the grouping rowsum() uses); codes need not all occur, but the sort's
// arrays have G + 1 entries, so G should not be far above the number of
// situations.
//
// As mxl_cluster_meat_parallel(), with choice situations as the units: the
// situations sorted by cluster (a stable counting sort, so a cluster's
// situations keep their order, as rowsum() adds them), the sorted sequence
// cut into chunks of similar numbers of rows processed in parallel, each
// thread adding the weighted compact scores of a cluster's situations into a
// dense n-vector over the parameters they touch (the cluster's support) and
// updating its accumulator once the cluster ends. A cluster cut by a chunk
// boundary leaves its partial sums in the chunks; they are added in chunk
// order after the parallel region and the cluster updated then. A cluster with
// a non-finite weight or score is updated over every parameter, as
// crossprod() spread its NaN or Inf (0 * NaN, 0 * Inf). acc_bytes chooses
// between the threads' private triangles and a shared result, as for
// mnl_bhhh_parallel().
// [[Rcpp::export]]
Rcpp::NumericMatrix mnl_cluster_meat_parallel(
    const arma::vec& theta,
    const arma::mat& X,
    const Rcpp::IntegerVector& alt_idx,
    const Rcpp::IntegerVector& choice_idx,
    const Rcpp::IntegerVector& M,
    const arma::vec& weights,
    const Rcpp::IntegerVector& cluster,
    const bool use_asc = true,
    const bool include_outside_option = false,
    const double acc_bytes = 2147483648.0
) {
  const int K = X.n_cols;
  const int n_params = theta.n_elem;

  // Inputs, layout and choices validated on the primary thread as in the
  // other weighted kernels; then the cluster codes.
  const MnlParams par = parse_mnl_theta(theta, K, use_asc, include_outside_option);
  const ChoiceLayout lay = choice_layout_build(X, alt_idx, M, use_asc,
                                               par.delta, &weights, &choice_idx);
  const choicer_off N = lay.N;
  validate_choices(lay, include_outside_option, "mnl_cluster_meat_parallel");
  if (static_cast<choicer_off>(cluster.size()) != N) {
    Rcpp::stop("cluster length (%d) does not match the number of choice "
               "situations (%d).", cluster.size(), N);
  }
  const int* cl = cluster.begin();
  int G = 0;
  for (choicer_off i = 0; i < N; ++i) {
    if (cl[i] < 1) {  // NA_INTEGER is negative
      Rcpp::stop("cluster codes must be positive integers (situation %d has %s).",
                 i + 1, cl[i] == NA_INTEGER ? "NA" : std::to_string(cl[i]));
    }
    G = std::max(G, cl[i]);
  }

  // Memory: per thread its slot map, compact score, situation buffers and
  // cluster sum, and its share of the chunks' cut parts (two per chunk, four
  // chunks per thread); shared X beta, the result and the accumulators, and
  // the sort's and the chunks' arrays.
  const int J_d = n_params - K;  // free constants
  const double a_max = std::min<double>(J_d, static_cast<double>(lay.max_m));
  const int n_threads = kernel_team_threads(N);
  kernel_check_memory("cluster meat", n_threads,
                      mnl_map_bytes(J_d, lay.max_m) + 8.0 * (K + a_max + 32.0) +
                          mnl_thread_bytes(lay, K, n_params, 0.0) +
                          9.0 * ClusterSum::bytes(n_params),
                      [&](const int t) {
                        return 8.0 * static_cast<double>(lay.n_rows) +
                               sym_acc_need(n_params, K, t, acc_bytes) +
                               8.0 * (static_cast<double>(N) + 2.0 * G +
                                      4.0 * t + 2.0);
                      });

  // Situations sorted by cluster: order[off[g - 1], off[g]) holds cluster g's
  // situations in increasing order.
  std::vector<choicer_off> off(static_cast<std::size_t>(G) + 1, 0);
  for (choicer_off i = 0; i < N; ++i) ++off[cl[i]];
  for (choicer_off g = 1; g <= G; ++g) off[g] += off[g - 1];
  std::vector<choicer_off> order(static_cast<std::size_t>(N));
  {
    std::vector<choicer_off> next(off.begin(), off.end() - 1);
    for (choicer_off i = 0; i < N; ++i) order[next[cl[i] - 1]++] = i;
  }

  // Chunks of the sorted sequence with similar numbers of rows (a
  // situation's cost grows with its rows): about four per thread.
  const choicer_off n_chunks = std::max<choicer_off>(
      1, std::min<choicer_off>(N, 4 * static_cast<choicer_off>(n_threads)));
  std::vector<choicer_off> chunk(static_cast<std::size_t>(n_chunks) + 1, N);
  {
    const double total = static_cast<double>(lay.n_rows);
    double cum = 0.0;
    choicer_off c = 1;
    chunk[0] = 0;
    for (choicer_off q = 0; q < N && c < n_chunks; ++q) {
      cum += static_cast<double>(lay.m(order[q]));
      if (cum >= total * static_cast<double>(c) / static_cast<double>(n_chunks)) {
        chunk[c++] = q + 1;
      }
    }
    for (; c < n_chunks; ++c) chunk[c] = N;
  }

  // Pre-compute base utility for all individuals (single BLAS call)
  arma::vec base_util = compute_base_util(X, par.beta, lay, use_asc, par.delta);

  // The result and the threads' accumulators, slot maps, score buffers and
  // cluster sums, and the chunks' partial sums of their cut clusters (two per
  // chunk: the cluster it starts in and the one it ends in), all allocated
  // here on the primary thread.
  SymAcc acc(n_params, n_threads, K, "cluster meat", acc_bytes);
  MnlSituationMaps maps = mnl_situation_maps(J_d, lay.max_m, n_threads);
  std::vector<std::vector<double>> scores =
      thread_vectors<double>(static_cast<std::size_t>(K + a_max), n_threads);
  std::vector<ClusterSum> sums;
  sums.reserve(n_threads);
  for (int t = 0; t < n_threads; ++t) sums.emplace_back(n_params);
  std::vector<ClusterSum> cut;  // chunk c: 2 c (head), 2 c + 1 (tail)
  cut.reserve(static_cast<std::size_t>(2 * n_chunks));
  for (choicer_off t = 0; t < 2 * n_chunks; ++t) cut.emplace_back(n_params);

#ifdef _OPENMP
#pragma omp parallel num_threads(n_threads)
#endif
  {
    const int tid = kernel_thread_num();
    MnlSituationMap& mp = *maps[tid];
    ClusterSum& cs = sums[tid];
    double* acc_t = acc.base(tid);
    double* s_buf = scores[tid].data();
    arma::vec diff_vec; // resized per situation

#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
    for (choicer_off c = 0; c < n_chunks; ++c) {
      const choicer_off q0 = chunk[c], q1 = chunk[c + 1];
      for (choicer_off q = q0; q < q1; ++q) {
        const choicer_off i = order[q];
        const int g = cl[i];
        if (q == q0 || cl[order[q - 1]] != g) cs.start(g, K);
        mnl_situation_map(lay, i, use_asc, include_outside_option, mp);
        arma::vec s_c(s_buf, K + mp.a, false, true);
        mnl_situation_score(X, lay, base_util, i, include_outside_option, mp,
                            diff_vec, s_c);
        cs.add(mp, K, s_c, weights[i]);
        mnl_situation_map_clear(mp);
        if (q + 1 == q1 || cl[order[q + 1]] != g) {
          // The cluster ends here or at the chunk's end: update it, or keep
          // the part of a cluster the chunk boundary cuts.
          const bool head = off[g - 1] < q0;   // began before this chunk
          const bool tail = off[g] > q1;       // continues after it
          if (!head && !tail) {
            cs.update(acc, acc_t);
            cs.clear();
          } else {
            cs.hand_to(cut[2 * c + (head ? 0 : 1)]);
          }
        }
      }
    }
  } // end parallel region

  // Clusters cut by chunk boundaries: their partial sums in chunk order, then
  // the rank-one update, through the first thread's buffer and with no other
  // thread running (no atomic additions into a shared result).
  ClusterSum& cs = sums[0];
  int g_open = 0;
  for (choicer_off c = 0; c < n_chunks; ++c) {
    for (int k = 0; k < 2; ++k) {
      ClusterSum& part = cut[2 * c + k];
      if (!part.used) continue;
      if (part.g != g_open) {
        if (g_open > 0) {
          cs.update(acc, acc.base(0), true);
          cs.clear();
        }
        cs.start(part.g, K);
        g_open = part.g;
      }
      cs.merge(part);
    }
  }
  if (g_open > 0) cs.update(acc, acc.base(0), true);

  acc.finish(false);
  return acc.result;
}

// Per-situation score matrix for the multinomial logit model (internal).
//
// Returns the N x n_params matrix whose row i is the weight-free score
// s_i = d log P_i / d theta (positive gradient of situation i's
// log-likelihood contribution). Same loop body as mnl_bhhh_parallel with the
// outer-product accumulator replaced by a row write; weights are applied on
// the R side (see .assemble_score_vcov in R/classes.R).
// [[Rcpp::export]]
arma::mat mnl_scores_parallel(
    const arma::vec& theta,
    const arma::mat& X,
    const Rcpp::IntegerVector& alt_idx,
    const Rcpp::IntegerVector& choice_idx,
    const Rcpp::IntegerVector& M,
    const bool use_asc = true,
    const bool include_outside_option = false
) {
  const int K = X.n_cols;
  const int n_params = theta.n_elem;

  const MnlParams par = parse_mnl_theta(theta, K, use_asc, include_outside_option);
  const ChoiceLayout lay = choice_layout_build(X, alt_idx, M, use_asc,
                                               par.delta, nullptr, &choice_idx);
  const choicer_off N = lay.N;

  // --- Serial pre-loop validation of chosen-alternative indices ---
  validate_choices(lay, include_outside_option, "mnl_scores_parallel");

  // Memory: each thread's score and situation buffers; X beta, the N x n
  // matrix and its copy into R
  kernel_check_memory("score matrix", kernel_team_threads(N),
                      mnl_thread_bytes(lay, K, n_params, 1.0), [&](const int) {
                        return 8.0 * static_cast<double>(lay.n_rows) +
                               16.0 * static_cast<double>(N) * n_params;
                      });

  // Pre-compute base utility for all individuals (single BLAS call)
  arma::vec base_util = compute_base_util(X, par.beta, lay, use_asc, par.delta);

  // Output: one row per choice situation (each written by exactly one
  // iteration, so no accumulator or critical section is needed).
  arma::mat scores(static_cast<arma::uword>(N), n_params);

#ifdef _OPENMP
#pragma omp parallel
#endif
  {
    arma::vec diff_vec; // resized per individual
    arma::vec s_i;      // per-individual score

#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
    for (choicer_off i = 0; i < N; ++i) {
      const int m_i         = lay.m(i);
      const int num_choices = include_outside_option ? m_i + 1 : m_i;
      const arma::uword start_idx = static_cast<arma::uword>(lay.row_off[i]);
      const arma::uword end_idx   = start_idx + m_i - 1;
      const auto X_i        = X.rows(start_idx, end_idx); // M[i] x K
      const AltCodes0 alt_idx0_i{lay.alt + lay.row_off[i]}; // M[i]

      // Build utility vector V_i
      arma::vec V_i(num_choices);
      arma::vec inside_utils = base_util.subvec(start_idx, end_idx);
      fill_choice_utilities(V_i, inside_utils, num_choices, include_outside_option);

      // Probabilities
      arma::vec P_i;
      stable_softmax(V_i, P_i);

      // Chosen alternative (validated serially above)
      int chosen_alt = lay.choice[i];
      if (!include_outside_option) chosen_alt -= 1;

      // diff_vec[a] = 1{a == chosen_alt} - P_i[a]  (same as gradient kernel)
      diff_vec.set_size(num_choices);
      diff_vec = -P_i;
      diff_vec(chosen_alt) += 1.0;

      // Assemble the (weight-free) per-individual score s_i.
      s_i.zeros(n_params);

      // Beta block
      if (include_outside_option) {
        const auto diff_inside = diff_vec.subvec(1, m_i); // m_i elements
        s_i.subvec(0, K - 1) = X_i.t() * diff_inside;
      } else {
        s_i.subvec(0, K - 1) = X_i.t() * diff_vec;
      }

      // Delta block (scatter; scale = 1)
      if (use_asc) {
        if (include_outside_option) {
          scatter_delta_grad(s_i, K, diff_vec.subvec(1, m_i), alt_idx0_i,
                             m_i, include_outside_option, 1.0);
        } else {
          scatter_delta_grad(s_i, K, diff_vec, alt_idx0_i,
                             m_i, include_outside_option, 1.0);
        }
      }

      scores.row(i) = s_i.t();
    } // end of i loop
  } // end parallel region

  return scores;
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
//' pred <- choicer:::mnl_predict(coef(fit), fit$data$X, fit$data$alt_idx,
//'   fit$data$M, use_asc = TRUE)
//' head(pred$choice_prob)
//' }
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List mnl_predict(
    const arma::vec& theta,
    const arma::mat& X,
    const Rcpp::IntegerVector& alt_idx,
    const Rcpp::IntegerVector& M,
    const bool use_asc = true,
    const bool include_outside_option = false
) {
  // Extract beta and delta from theta
  const int K = X.n_cols;

  const MnlParams par = parse_mnl_theta(theta, K, use_asc, include_outside_option);
  const ChoiceLayout lay = choice_layout_build(X, alt_idx, M, use_asc, par.delta);
  const choicer_off N = lay.N;

  // Pre-compute base utility for all individuals (single BLAS call)
  arma::vec base_util = compute_base_util(X, par.beta, lay, use_asc, par.delta);

  // Preallocate predicted values
  arma::vec V_all = arma::zeros(X.n_rows);
  arma::vec P_all = arma::zeros(X.n_rows);

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
  for (choicer_off i = 0; i < N; ++i) {
    const int m_i         = lay.m(i);
    const int num_choices = include_outside_option ? m_i + 1 : m_i;
    const arma::uword start_idx = static_cast<arma::uword>(lay.row_off[i]);
    const arma::uword end_idx   = start_idx + m_i - 1;

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

// Prediction of shares (internal function), from the base utilities of all
// rows (X * beta plus the rows' ASCs). denominator is
// shares_denominator(weights), which the callers check before they form
// the utilities.
static arma::vec mnl_predict_shares_internal(
    const arma::vec& base_util,        // sum(M) x 1 vector of base utilities
    const ChoiceLayout& lay,           // layout of the stacked design
    const arma::vec& weights,          // N x 1 vector with weights for each observation
    const double denominator,          // sum of the weights
    const int num_alts,                // total number of distinct alternatives
    const bool include_outside_option  // whether to include an outside option
) {
  const choicer_off N = lay.N;

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
  for (choicer_off i = 0; i < N; ++i) {
    const int m_i         = lay.m(i);
    const int num_choices = include_outside_option ? m_i + 1 : m_i;
    const arma::uword start_idx = static_cast<arma::uword>(lay.row_off[i]);
    const arma::uword end_idx   = start_idx + m_i - 1;
    const AltCodes0 alt_idx0_i{lay.alt + lay.row_off[i]}; // M[i]

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
        local_shares(alt_idx0_i[a] + 1) += weights[i] * P_i(a + 1);
      } else {
        local_shares(alt_idx0_i[a]) += weights[i] * P_i(a);
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
//' shares <- choicer:::mnl_predict_shares(coef(fit), fit$data$X, fit$data$alt_idx,
//'   fit$data$M, fit$data$weights, use_asc = TRUE)
//' shares
//' }
//' @keywords internal
// [[Rcpp::export]]
arma::vec mnl_predict_shares(
    const arma::vec& theta,            // K + J - 1 or K + J vector with model parameters
    const arma::mat& X,                // sum(M) x K design matrix with covariates. M\[i] x K matrix for individual i
    const Rcpp::IntegerVector& alt_idx, // sum(M) x 1 vector with indices of alternatives within each choice set; 1-based indexing
    const Rcpp::IntegerVector& M,      // N x 1 vector with number of alternatives for each individual
    const arma::vec& weights,          // N x 1 vector with weights for each observation
    const bool use_asc = true,         // whether to use alternative-specific constants
    const bool include_outside_option = false  // whether to include outside option normalized to 0 (if so, the outside option is not included in the data)
) {
  // Extract beta and delta from theta
  const int K = X.n_cols;

  const MnlParams par = parse_mnl_theta(theta, K, use_asc, include_outside_option);
  const ChoiceLayout lay = choice_layout_build(X, alt_idx, M, use_asc, par.delta,
                                               &weights);

  // total number of distinct alternatives (alt_idx.max() on the arma::uvec
  // argument this replaces threw on an empty alt_idx; kept)
  if (lay.n_rows == 0) {
    throw std::logic_error("Mat::max(): object has no elements");
  }
  const int num_alts = compute_J_total(lay.J, include_outside_option);

  const double denominator = shares_denominator(weights);
  const arma::vec base_util =
      compute_base_util(X, par.beta, lay, use_asc, par.delta);
  arma::vec global_shares = mnl_predict_shares_internal(
    base_util, lay, weights, denominator, num_alts, include_outside_option
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
  const Rcpp::IntegerVector& alt_idx,
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
  if (target_shares.n_elem != static_cast<arma::uword>(num_alts)) {
    Rcpp::stop("Error: target_shares must have the same length as the total number of alternatives.");
  }
  if (arma::any(target_shares <= 0)) {
    Rcpp::stop("Error: all target_shares must be strictly positive (log(share) is undefined otherwise).");
  }
  const ChoiceLayout lay = choice_layout_build(X, alt_idx, M, use_asc, delta,
                                               &weights);

  // X * beta, once: only delta changes across the iterations. Each one adds
  // its ASCs to a copy of it, in one buffer kept across them.
  const double denominator = shares_denominator(weights);
  const arma::vec Xb = X * beta;
  arma::vec base_util;
  auto predict_shares = [&](const arma::vec& inside_delta) {
    base_util = Xb;
    add_row_asc(base_util, lay, inside_delta);
    return mnl_predict_shares_internal(base_util, lay, weights, denominator,
                                       num_alts, include_outside_option);
  };

  // The iteration bookkeeping (delta_old/delta_new, target/predicted log-shares,
  // residual) lives in the outside-inclusive share space of length num_alts:
  //   index 0 = outside option (when present), indices 1..J = inside alts.
  // But predict_shares() adds the ASCs by inside alternative (add_row_asc()
  // reads delta at each row's code minus 1, in 0..J-1), so it takes a length-J
  // inside-delta vector (the outside option is handled separately via
  // include_outside_option). We therefore feed it
  // delta_old.subvec(1, num_alts - 1) when an outside option is present, and
  // pin the outside slot delta_old[0] = 0 throughout (the outside option's
  // utility is the fixed normalization).
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
  arma::vec log_shares_old = predict_shares(inside_delta_old);
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
      log_shares_old = predict_shares(inside_delta_old);
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

// Add situation i's terms of the Hessian into the accumulator's upper
// triangle, each entry once, as the dense code added them: it subtracted w_i
// times each block from its accumulator and negated the sum at the end, so
// each term here is -(entry * w_i) and finish() negates (entries no situation
// touches are -0, as before). Inside the block (the K coefficients, then the
// constants in slot order): BB (its upper triangle), BD (K x a, column s the
// constant in slot s) and the constants' block, formed here from the slots'
// probabilities mu_d: mu_d[s] - mu_d[s]^2 on the diagonal, 0 - mu_d[s]
// mu_d[t] off it. Outside the block the dense code's full-width subtractions
// left 0 - mean * 0 in the row of each mean: +-0 for a finite mean, which
// changes nothing, but NaN, or NA for an NA, for a non-finite one
// (overflowing utilities, a covariate past the largest double), added here
// to that row's entries outside the block. The entries inside keep their
// values (an infinite mean leaves infinities there, as before). A NaN mean
// makes its whole row NaN in this thread's sums, inside the block through
// the products and outside through these terms, so nan_row marks the row
// and its later spreads, which change nothing, are skipped.
template <bool Shared>
inline void mnl_hessian_add(const SymAcc& acc, double* buf,
                            const MnlSituationMap& mp, const int K,
                            const int J_d, const double* BB,
                            const double* mu_b, const double* BD,
                            const double* mu_d, const double w,
                            char* nan_row) {
  const std::size_t Ks = static_cast<std::size_t>(K);
  for (int c = 0; c < K; ++c) {
    for (int r = 0; r <= c; ++r) {
      double v = BB[r + c * Ks];
      v *= w;
      acc.add<Shared>(buf, r, c, -v);
    }
  }
  for (int s = 0; s < mp.a; ++s) {
    const int g = K + mp.delta[s];
    for (int r = 0; r < K; ++r) {
      double v = BD[r + s * Ks];
      v *= w;
      acc.add<Shared>(buf, r, g, -v);
    }
    const double ms = mu_d[s];
    for (int t = 0; t < s; ++t) {
      double v = 0.0 - mu_d[t] * ms;
      v *= w;
      acc.add<Shared>(buf, K + mp.delta[t], g, -v);
    }
    double v = ms - ms * ms;
    v *= w;
    acc.add<Shared>(buf, g, g, -v);
  }
  for (int l = 0; l < K + mp.a; ++l) {
    const double mean = l < K ? mu_b[l] : mu_d[l - K];
    if (std::isfinite(mean)) continue;
    const int g = mp.global(l, K);
    if (std::isnan(mean)) {
      if (nan_row[g]) continue;
      nan_row[g] = 1;
    }
    double v = 0.0 - mean * 0.0;
    v *= w;
    for (int j = 0; j < J_d; ++j) {
      if (mp.slot_of[j] >= 0) continue;  // inside the block
      const int h = K + j;
      acc.add<Shared>(buf, std::min(g, h), std::max(g, h), -v);
    }
  }
}

//' Hessian matrix for multinomial logit model
//'
//' The Hessian of the negated log-likelihood: the weighted sum over choice
//' situations of the probability-weighted covariance of each alternative's
//' covariates and constant indicator. Each situation's block is formed only
//' over the parameters its likelihood involves (the coefficients and the
//' constants of the alternatives it offers), so the work per situation does
//' not grow with the number of alternatives; the n_params x n_params result
//' and the threads' accumulators do (see acc_bytes).
//'
//' @param theta K + J - 1 or K + J vector with model parameters
//' @param X sum(M) x K design matrix with covariates. Stacks M\[i] x K matrices for individual i.
//' @param alt_idx sum(M) x 1 vector with indices of alternatives within each choice set; 1-based indexing
//' @param choice_idx N x 1 vector with indices of chosen alternatives; 1-based indexing relative to X; 0 is used if include_outside_option=True
//' @param M N x 1 vector with number of alternatives for each individual
//' @param weights N x 1 vector with weights for each observation
//' @param use_asc whether to use alternative-specific constants
//' @param include_outside_option whether to include outside option normalized to 0 (if so, the outside option is not included in the data)
//' @param acc_bytes Bytes for the threads' private upper triangles (2 GiB by
//'   default), as in \code{mnl_bhhh_parallel()}; \code{0} forces the shared
//'   result (for tests).
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
//' H <- choicer:::mnl_loglik_hessian_parallel(coef(fit), fit$data$X, fit$data$alt_idx,
//'   fit$data$choice_idx, fit$data$M, fit$data$weights)
//' dim(H)
//' }
//' @keywords internal
// [[Rcpp::export]]
Rcpp::NumericMatrix mnl_loglik_hessian_parallel(
    const arma::vec& theta,
    const arma::mat& X,
    const Rcpp::IntegerVector& alt_idx,
    const Rcpp::IntegerVector& choice_idx, // Note: choice_idx is not needed for Hessian
    const Rcpp::IntegerVector& M,
    const arma::vec& weights,
    const bool use_asc = true,
    const bool include_outside_option = false,
    const double acc_bytes = 2147483648.0
) {
  (void)choice_idx; // silence unused-parameter warning

  // Extract beta and delta from theta
  const int K = X.n_cols;
  const int n_params = theta.n_elem;

  // Split theta into beta and delta (ASCs) as in your gradient function
  const MnlParams par = parse_mnl_theta(theta, K, use_asc, include_outside_option);
  // Layout of the stacked design (each individual's block in X / alt_idx)
  const ChoiceLayout lay = choice_layout_build(X, alt_idx, M, use_asc, par.delta,
                                               &weights);
  const choicer_off N = lay.N;

  // Memory: each thread's slot map, a situation's blocks (BB, BD, the means
  // and K products), the non-finite-weight records, the NaN rows' marks and
  // situation buffers; X beta, the result and the threads' accumulators
  const int J_d = n_params - K;  // free constants
  const double a_max = std::min<double>(J_d, static_cast<double>(lay.max_m));
  const int n_threads = kernel_team_threads(N);
  kernel_check_memory("Hessian", n_threads,
                      mnl_map_bytes(J_d, lay.max_m) +
                          8.0 * (K * (K + 2.0) + (K + 1.0) * a_max + 64.0) +
                          7.0 * n_params + 128.0 + mnl_thread_bytes(lay, K, n_params, 0.0),
                      [&](const int t) {
                        return 8.0 * static_cast<double>(lay.n_rows) +
                               sym_acc_need(n_params, K, t, acc_bytes);
                      });

  // Pre-compute base utility for all individuals (single BLAS call)
  arma::vec base_util = compute_base_util(X, par.beta, lay, use_asc, par.delta);

  // The result and the threads' accumulators, slot maps, block buffers and
  // non-finite-weight records, all allocated here on the primary thread
  SymAcc acc(n_params, n_threads, K, "Hessian", acc_bytes);
  MnlSituationMaps maps = mnl_situation_maps(J_d, lay.max_m, n_threads);
  const std::size_t blk = static_cast<std::size_t>(K * (K + 2.0) + (K + 1.0) * a_max);
  std::vector<std::vector<double>> blocks = thread_vectors<double>(blk, n_threads);
  std::vector<NonFiniteRecord> nonfinite = nonfinite_records(n_params, n_threads);
  std::vector<std::vector<char>> nan_rows =
      thread_vectors<char>(static_cast<std::size_t>(n_params), n_threads, 0);
  const int o = include_outside_option ? 1 : 0;  // slot of the first inside alt
  const double* Xm = X.memptr();
  const arma::uword n_rows = X.n_rows;

#ifdef _OPENMP
#pragma omp parallel num_threads(n_threads)
#endif
  {
    const int tid = kernel_thread_num();
    MnlSituationMap& mp = *maps[tid];
    NonFiniteRecord& nf = nonfinite[tid];
    char* nan_row = nan_rows[tid].data();
    double* acc_t = acc.base(tid);
    // BB (K x K), mu_b (K), BD (K x a), mu_d (a) and K products, in one
    // buffer
    const std::size_t Ks = static_cast<std::size_t>(K);
    double* BB = blocks[tid].data();
    double* mu_b = BB + Ks * Ks;
    double* BD = mu_b + K;

#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
    for (choicer_off i = 0; i < N; ++i) {
      const int m_i         = lay.m(i);
      const int num_choices = include_outside_option ? m_i + 1 : m_i;
      const arma::uword start_idx = static_cast<arma::uword>(lay.row_off[i]);
      const arma::uword end_idx   = start_idx + m_i - 1;
      const double w_i      = weights[i];

      // Build utility vector V_i
      arma::vec V_i(num_choices);
      arma::vec inside_utils = base_util.subvec(start_idx, end_idx);
      fill_choice_utilities(V_i, inside_utils, num_choices, include_outside_option);

      // Calculate Probabilities ---------------------------------------------
      arma::vec P_i;
      stable_softmax(V_i, P_i);

      // The situation's slots and blocks: BB = sum_a p_a x_a x_a' - mu_b mu_b'
      // (upper triangle), mu_b = sum_a p_a x_a, BD = sum_a p_a x_a e_s(a)' -
      // mu_b mu_d', mu_d[s] = sum of p_a over the rows of slot s
      mnl_situation_map(lay, i, use_asc, include_outside_option, mp);
      const int a = mp.a;
      double* mu_d = BD + Ks * a;
      std::fill(BB, mu_d + a, 0.0);
      for (int r = 0; r < m_i; ++r) {
        const double p_r = P_i[r + o];
        const double* x = Xm + start_idx + r;  // x_r[k] = x[k * n_rows]
        const int k = mp.row_slot[r];
        // Row c of BB's upper triangle from (p x_c) x_q, q >= c, as the
        // dense code formed it
        for (int c = 0; c < K; ++c) {
          const double p_x_c = p_r * x[c * n_rows];
          for (int q = c; q < K; ++q) BB[c + q * Ks] += p_x_c * x[q * n_rows];
          mu_b[c] += p_x_c;
          if (k >= 0) BD[c + k * Ks] += p_x_c;
        }
        if (k >= 0) mu_d[k] += p_r;
      }
      // Each product of means rounded before it is subtracted, as the dense
      // code rounded them (Armadillo's syrk_vec for mu_b mu_b', BLAS for
      // mu_b mu_d'): a column's products are stored first, so that no
      // compiler fuses them into the subtraction (a fused product past the
      // largest double, from p x x' past it, would leave Inf where the dense
      // code left Inf - Inf = NaN)
      double* prod = mu_d + a;  // K entries
      for (int c = 0; c < K; ++c) {
        for (int q = 0; q <= c; ++q) prod[q] = mu_b[q] * mu_b[c];
        for (int q = 0; q <= c; ++q) BB[q + c * Ks] -= prod[q];
      }
      for (int s = 0; s < a; ++s) {
        for (int c = 0; c < K; ++c) prod[c] = mu_b[c] * mu_d[s];
        for (int c = 0; c < K; ++c) BD[c + s * Ks] -= prod[c];
      }

      // A non-finite weight multiplied the zeros outside the block too:
      // recorded, and spread after the parallel region
      if (!std::isfinite(w_i)) {
        nonfinite_weight(w_i, K + a, [&](const int l) { return mp.global(l, K); }, nf);
      }

      if (acc.shared) {
        mnl_hessian_add<true>(acc, acc_t, mp, K, J_d, BB, mu_b, BD, mu_d, w_i,
                              nan_row);
      } else {
        mnl_hessian_add<false>(acc, acc_t, mp, K, J_d, BB, mu_b, BD, mu_d, w_i,
                               nan_row);
      }
      mnl_situation_map_clear(mp);
    } // end of i loop
  } // end parallel region

  // The information (Hessian of the negated log-likelihood), with the dense
  // code's spread of non-finite weights
  acc.finish(true);
  apply_nonfinite(acc, nonfinite);
  return acc.result;
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
//' elas <- choicer:::mnl_elasticities_parallel(coef(fit), fit$data$X, fit$data$alt_idx,
//'   fit$data$choice_idx, fit$data$M, fit$data$weights, elast_var_idx = 1L)
//' elas
//' }
//' @keywords internal
// [[Rcpp::export]]
arma::mat mnl_elasticities_parallel(
    const arma::vec& theta,
    const arma::mat& X,
    const Rcpp::IntegerVector& alt_idx,
    const Rcpp::IntegerVector& choice_idx, // Kept for consistency, but not used
    const Rcpp::IntegerVector& M,
    const arma::vec& weights,
    const int elast_var_idx,
    const bool use_asc = true,
    const bool include_outside_option = false
) {
  (void)choice_idx;

  // --- 1. Parameter and Variable Setup ---
  const int K = X.n_cols;

  // Convert 1-based R index to 0-based C++ index
  const int var_idx = elast_var_idx - 1;
  if (var_idx < 0 || var_idx >= K) {
    Rcpp::stop("elast_var_idx is out of bounds for design matrix X.");
  }

  // Extract beta and delta (same logic as loglik function)
  const MnlParams par = parse_mnl_theta(theta, K, use_asc, include_outside_option);
  const ChoiceLayout lay = choice_layout_build(X, alt_idx, M, use_asc, par.delta,
                                               &weights);
  const choicer_off N = lay.N;
  const double beta_k = par.beta(var_idx); // The coefficient for our variable

  // Determine total number of alternatives for the output matrix
  // J_inside = number of alternatives (with or without normalization)
  // J_total = total size of matrix (J_inside + 1 if outside option)
  const int J_inside = compute_J_inside(use_asc, par.delta, lay);
  const int J_total = compute_J_total(J_inside, include_outside_option);

  // Pre-compute base utility for all individuals (single BLAS call)
  arma::vec base_util = compute_base_util(X, par.beta, lay, use_asc, par.delta);

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
    // Local-to-global alternative map, sized once for the largest choice set
    std::vector<int> global_j_map(static_cast<std::size_t>(lay.max_m) + 1);

#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
    for (choicer_off i = 0; i < N; ++i) {
      const int m_i         = lay.m(i);
      const int num_choices = include_outside_option ? m_i + 1 : m_i;
      const arma::uword start_idx = static_cast<arma::uword>(lay.row_off[i]);
      const arma::uword end_idx   = start_idx + m_i - 1;
      const double w_i      = weights[i];
      const auto X_i        = X.rows(start_idx, end_idx); // M[i] x K

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
      fill_global_alt_map(global_j_map.data(), lay.alt + lay.row_off[i], m_i,
                          include_outside_option);

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
//' dr <- choicer:::mnl_diversion_ratios_parallel(coef(fit), fit$data$X, fit$data$alt_idx,
//'   fit$data$M, fit$data$weights)
//' dr
//' }
//' @keywords internal
// [[Rcpp::export]]
arma::mat mnl_diversion_ratios_parallel(
    const arma::vec& theta,
    const arma::mat& X,
    const Rcpp::IntegerVector& alt_idx,
    const Rcpp::IntegerVector& M,
    const arma::vec& weights,
    const bool use_asc = true,
    const bool include_outside_option = false
) {
  // --- 1. Parameter and Variable Setup ---
  const int K = X.n_cols;

  // Extract beta and delta (same logic as elasticities/loglik)
  const MnlParams par = parse_mnl_theta(theta, K, use_asc, include_outside_option);
  const ChoiceLayout lay = choice_layout_build(X, alt_idx, M, use_asc, par.delta,
                                               &weights);
  const choicer_off N = lay.N;

  // Determine total number of alternatives for the output matrix
  const int J_inside = compute_J_inside(use_asc, par.delta, lay);
  const int J_total = compute_J_total(J_inside, include_outside_option);

  // Pre-compute base utility for all individuals (single BLAS call)
  arma::vec base_util = compute_base_util(X, par.beta, lay, use_asc, par.delta);

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
    // Local-to-global alternative map, sized once for the largest choice set
    std::vector<int> global_j_map(static_cast<std::size_t>(lay.max_m) + 1);

#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
    for (choicer_off i = 0; i < N; ++i) {
      const int m_i         = lay.m(i);
      const int num_choices = include_outside_option ? m_i + 1 : m_i;
      const arma::uword start_idx = static_cast<arma::uword>(lay.row_off[i]);
      const arma::uword end_idx   = start_idx + m_i - 1;
      const double w_i      = weights[i];

      // --- 2. Compute Probabilities ---
      arma::vec V_i(num_choices);
      arma::vec inside_utils = base_util.subvec(start_idx, end_idx);
      fill_choice_utilities(V_i, inside_utils, num_choices, include_outside_option);

      arma::vec P_i;
      stable_softmax(V_i, P_i);

      // --- 3. Map local indices to global alternative indices ---
      fill_global_alt_map(global_j_map.data(), lay.alt + lay.row_off[i], m_i,
                          include_outside_option);

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
