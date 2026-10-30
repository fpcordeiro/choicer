// [[Rcpp::depends(RcppArmadillo)]]
#include "choicer.h"
#include "choicer_internal.h"
#include "choicer_derivs.h"
#include "blp_engine.h"

// ============================================================================
// The MNL kernels' working arrays
//
// One thread's arrays for the choice situation in hand, allocated on the
// primary thread for the largest one (max_m rows, K covariates) and first
// touched by the thread, so that no loop over situations allocates: the
// utilities and probabilities of its choices (the outside option first),
// 1{chosen} - P, the covariate of an elasticity, a copy of its rows of X for
// the products X_i' v that do not reach BLAS (xt_times(); Armadillo copied
// X.rows() for every product), and a K-vector for a product added into a sum
// afterwards. Each array starts a 64-byte line.
// ============================================================================
struct MnlWork {
  AlignedBuffer buf;
  double *v = nullptr, *p = nullptr, *diff = nullptr, *xk = nullptr;
  double *xi = nullptr, *kv = nullptr;

  MnlWork(const choicer_off max_m, const int K) {
    const std::size_t m1 = aligned_len(static_cast<std::size_t>(max_m) + 1);
    const std::size_t mk = aligned_len(static_cast<std::size_t>(max_m) *
                                       static_cast<std::size_t>(K));
    buf = AlignedBuffer(4 * m1 + mk + aligned_len(static_cast<std::size_t>(K)));
    v = buf.data;
    p = v + m1;
    diff = p + m1;
    xk = diff + m1;
    xi = xk + m1;
    kv = xi + mk;
  }
};

// Every thread's working arrays, one allocation each, on the primary thread.
inline std::vector<MnlWork> mnl_work(const ChoiceLayout& lay, const int K,
                                     const int n_threads) {
  std::vector<MnlWork> work;
  work.reserve(static_cast<std::size_t>(n_threads));
  for (int t = 0; t < n_threads; ++t) work.emplace_back(lay.max_m, K);
  return work;
}

// Memory of the MNL kernels' loop over the choice situations, an upper bound
// per thread: MnlWork's arrays (each rounded up to a 64-byte line) and n_vecs
// padded vectors of the n parameters (the thread's sums or scores).
inline double mnl_thread_bytes(const ChoiceLayout& lay, const int K,
                               const int n_params, const double n_vecs) {
  const double m = static_cast<double>(lay.max_m);
  return sizeof(MnlWork) +
         AlignedBuffer::bytes(4.0 * (m + 8.0) + (m * K + 7.0) + (K + 7.0)) +
         sizeof(double) * n_vecs * (static_cast<double>(n_params) + 16.0);
}

// Situation i's choice probabilities into w.p (the outside option first), its
// utilities into w.v shifted by their maximum, and the log of the denominator
// returned: fill_choice_utilities() and stable_softmax(), on the thread's
// arrays (stable_softmax_n() applies Armadillo's operations in its order).
inline double mnl_softmax(const arma::vec& base_util, const ChoiceLayout& lay,
                          const choicer_off i, const bool include_outside_option,
                          MnlWork& w) {
  const int m_i = lay.m(i);
  const double* bu = base_util.memptr() + lay.row_off[i];
  double* v = w.v;
  if (include_outside_option) {
    v[0] = 0.0;
    std::copy(bu, bu + m_i, v + 1);
  } else {
    std::copy(bu, bu + m_i, v);
  }
  return stable_softmax_n(v, w.p, include_outside_option ? m_i + 1 : m_i);
}

// 1{a == chosen} - P[a] into w.diff, as diff_vec = -P; diff_vec(chosen) += 1.
inline void mnl_diff(MnlWork& w, const int num_choices, const int chosen) {
  for (int a = 0; a < num_choices; ++a) w.diff[a] = -w.p[a];
  w.diff[chosen] += 1.0;
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
// nothing; slot_of holds the situation the map was built for (-1 elsewhere),
// kept while the next situations have its codes.
// ============================================================================
struct MnlSituationMap {
  std::vector<int> slot_of;   // free constant j -> slot, or -1 (J_d entries)
  std::vector<int> delta;     // slot -> free constant, increasing (a of them)
  std::vector<int> row_slot;  // row of the situation -> slot, or -1
  std::vector<int> glob;      // local index -> global parameter (the K coefficients first)
  int a = 0;                  // free constants of the situation
  int K = 0;                  // coefficients
  choicer_off built = -1;     // first row of the situation the map was built for, or -1
  int built_m = 0;            // and its rows
  char pad[128];              // keeps the next thread's map off these lines

  // Each array is 32 entries (128 bytes) beyond its largest size, so that no
  // two threads' arrays share a cache line.
  MnlSituationMap(const int J_d, const choicer_off max_m, const int K_)
      : slot_of(static_cast<std::size_t>(J_d) + 32, -1),
        delta(static_cast<std::size_t>(std::min<choicer_off>(J_d, max_m)) + 32),
        row_slot(static_cast<std::size_t>(max_m) + 32),
        glob(static_cast<std::size_t>(K_) +
             static_cast<std::size_t>(std::min<choicer_off>(J_d, max_m)) + 32),
        K(K_) {
    for (int k = 0; k < K; ++k) glob[k] = k;
  }
  // Global parameter of local index i (the K coefficients first).
  int global(const int i, const int) const { return glob[i]; }
  // Local index of global parameter g, or -1 outside the situation's block.
  int local(const int g, const int K_) const {
    if (g < K_) return g;
    const int k = slot_of[g - K_];
    return k >= 0 ? K_ + k : -1;
  }
};

// Bytes of one thread's MnlSituationMap.
inline double mnl_map_bytes(const int J_d, const choicer_off max_m, const int K) {
  const double m = static_cast<double>(max_m);
  return sizeof(int) * (J_d + 2.0 * std::min<double>(J_d, m) + m + K + 128.0) +
         sizeof(MnlSituationMap);
}

// Every thread's map, each in its own allocation, on the primary thread.
using MnlSituationMaps = std::vector<std::unique_ptr<MnlSituationMap>>;
inline MnlSituationMaps mnl_situation_maps(const int J_d, const choicer_off max_m,
                                           const int K, const int n_threads) {
  MnlSituationMaps maps;
  maps.reserve(n_threads);
  for (int i = 0; i < n_threads; ++i) {
    maps.emplace_back(new MnlSituationMap(J_d, max_m, K));
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
// order, the slot of each of its rows and the global parameter of each
// local index. The map depends on the rows' codes only, so it is kept while
// the situations a thread meets list the codes of the one it was built for.
// Constants met in increasing order (prepared data list a situation's
// alternatives in order) need no sort; otherwise they are sorted (distinct
// keys, so any sort gives the same order) and the rows' slots set again.
CHOICER_ALWAYS_INLINE void mnl_situation_map(const ChoiceLayout& lay,
                                             const choicer_off i,
                                             const bool use_asc,
                                             const bool include_outside_option,
                                             MnlSituationMap& mp) {
  const choicer_off r0 = lay.row_off[i];
  const int* alt = lay.alt + r0;
  const int m = lay.m(i);
  if (mp.built >= 0 && mp.built_m == m &&
      std::equal(alt, alt + m, lay.alt + mp.built)) {
    return;
  }
  int* slot_of = mp.slot_of.data();
  int* delta = mp.delta.data();
  int* row_slot = mp.row_slot.data();
  for (int k = 0; k < mp.a; ++k) slot_of[delta[k]] = -1;  // the previous map's
  int a = 0;
  int last = -1;
  bool increasing = true;
  for (int r = 0; r < m; ++r) {
    const int j = mnl_free_delta(alt[r] - 1, use_asc, include_outside_option);
    if (j < 0) {
      row_slot[r] = -1;
      continue;
    }
    if (slot_of[j] < 0) {
      if (j < last) increasing = false;
      last = j;
      slot_of[j] = a;
      delta[a++] = j;
    }
    row_slot[r] = slot_of[j];
  }
  if (!increasing) {
    std::sort(delta, delta + a);
    for (int k = 0; k < a; ++k) slot_of[delta[k]] = k;
    for (int r = 0; r < m; ++r) {
      const int j = mnl_free_delta(alt[r] - 1, use_asc, include_outside_option);
      row_slot[r] = j >= 0 ? slot_of[j] : -1;
    }
  }
  const int K = mp.K;
  for (int k = 0; k < a; ++k) mp.glob[K + k] = K + delta[k];
  mp.a = a;
  mp.built = r0;
  mp.built_m = m;
}

// The weight-free score of situation i on its compact block, into s_c (K + a
// entries, the map built for i), and the situation's 1{chosen} - P into
// w.diff: the coefficients' entries as X_i' diff (xt_times(), the expression
// of the gradient kernel's OPG diagonal, written into s_c's first K entries,
// where the dense code copied the product), then each slot's differences in
// row order at scale 1.0, as scatter_delta_grad() adds them into the
// n_params-vector.
inline void mnl_situation_score(const arma::mat& X, const ChoiceLayout& lay,
                                const arma::vec& base_util, const choicer_off i,
                                const bool include_outside_option,
                                const MnlSituationMap& mp, MnlWork& w,
                                arma::vec& s_c) {
  const int m_i         = lay.m(i);
  const int num_choices = include_outside_option ? m_i + 1 : m_i;
  const int K = X.n_cols;
  const int o = include_outside_option ? 1 : 0;  // slot of the first inside alt

  // Probabilities, and diff[a] = 1{a == chosen_alt} - P[a] (as the gradient
  // kernel forms it; the choice was validated serially by the caller)
  mnl_softmax(base_util, lay, i, include_outside_option, w);
  int chosen_alt = lay.choice[i];
  if (!include_outside_option) chosen_alt -= 1;
  mnl_diff(w, num_choices, chosen_alt);
  const double* diff = w.diff;

  s_c.zeros();
  xt_times(X, static_cast<arma::uword>(lay.row_off[i]), m_i, diff + o, false,
           0.0, s_c.memptr(), w.xi);
  for (int r = 0; r < m_i; ++r) {
    const int k = mp.row_slot[r];
    if (k >= 0) s_c[K + k] += 1.0 * diff[r + o];
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
  const int o = include_outside_option ? 1 : 0;  // slot of the first inside alt

  // Every thread's working arrays and sums, allocated here on the primary
  // thread: the gradient and, with OpgDiag, the diagonal, the constants'
  // scores (zero between situations) and a situation's terms
  std::vector<MnlWork> work = mnl_work(lay, K, n_threads);
  std::vector<std::vector<double>> grads =
      thread_vectors<double>(static_cast<std::size_t>(n_params), n_threads, 0.0);
  std::vector<std::vector<double>> diags = thread_vectors<double>(
      OpgDiag ? static_cast<std::size_t>(n_params) : 0, n_threads, 0.0);
  std::vector<std::vector<double>> sds = thread_vectors<double>(
      OpgDiag && use_asc ? static_cast<std::size_t>(n_params) : 0, n_threads, 0.0);
  std::vector<OpgDiagTerms> terms_t;
  terms_t.reserve(static_cast<std::size_t>(n_threads));
  for (int t = 0; t < n_threads; ++t) {
    terms_t.emplace_back(OpgDiag ? K + lay.max_m : 0);
  }

#ifdef _OPENMP
  const int chunk = situation_chunk(lay);  // situations per dispatch
#pragma omp parallel num_threads(n_threads)
#endif
  {
    const int tid = kernel_thread_num();
    MnlWork& w = work[tid];
    // Thread-local accumulators
    double local_loglik = 0.0;
    arma::vec local_grad(grads[tid].data(), n_params, false, true);
    bool local_nonfinite = false;
    arma::vec local_diag(diags[tid].data(), OpgDiag ? n_params : 0, false, true);
    arma::vec s_d(sds[tid].data(), OpgDiag && use_asc ? n_params : 0, false, true);
    OpgDiagTerms& terms = terms_t[tid];

#ifdef _OPENMP
#pragma omp for schedule(dynamic, chunk)
#endif
    for (choicer_off i = 0; i < N; ++i) {
      const int m_i         = lay.m(i);
      const int num_choices = include_outside_option ? m_i + 1 : m_i;
      const arma::uword start_idx = static_cast<arma::uword>(lay.row_off[i]);
      const double w_i      = weights[i];
      const AltCodes0 alt_idx0_i{lay.alt + lay.row_off[i]}; // M[i]

      // log-likelihood -------------------------------------------------------

      // Choice probabilities (w.v holds the utilities, shifted)
      const double log_denom = mnl_softmax(base_util, lay, i, include_outside_option, w);

      // Identify chosen alternative (validated serially above)
      int chosen_alt = lay.choice[i];
      if (!include_outside_option) {
        chosen_alt -= 1; // shift by 1 for inside-only indexing
      }

      // Probability of chosen alternative
      double V_choice = w.v[chosen_alt];
      double log_P_choice = V_choice - log_denom;
      if (!std::isfinite(log_P_choice)) {
        local_nonfinite = true; // H1: note without printing (R-API unsafe here)
        log_P_choice = -1e10;  // fallback to a large negative value
      }

      // Accumulate local weighted log-likelihood
      local_loglik += w_i * log_P_choice;

      // Gradient -------------------------------------------------------------
      // diff_vec[a] = 1{a == chosen_alt} - P_i[a]
      mnl_diff(w, num_choices, chosen_alt);
      const arma::vec diff_vec(w.diff, num_choices, false, true);

      // ---- Beta block: single BLAS dgemv (w_i X_i' diff, the inside
      // alternatives' differences), then added ----
      xt_times(X, start_idx, m_i, w.diff + o, true, w_i, w.kv, w.xi);
      const arma::vec x_d(w.kv, K, false, true);
      local_grad.subvec(0, K - 1) += x_d;

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
        // the beta block's scores, X_i' diff
        xt_times(X, start_idx, m_i, w.diff + o, false, 0.0, w.kv, w.xi);
        const double* s_b = w.kv;
        if (use_asc) {
          if (include_outside_option) {
            scatter_delta_grad(s_d, K, diff_vec.subvec(1, m_i), alt_idx0_i, m_i,
                               include_outside_option, 1.0);
          } else {
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
//'   arguments, and \code{overflow}: whether the objective or the gradient
//'   is not the negated log-likelihood's at \code{theta}, because a
//'   situation's log-probability, the total or a gradient entry was not
//'   finite and was replaced (the log-probability by -1e10, the total by
//'   1e10 with a zero gradient, the entry by 0; the mixed logit kernel's flag
//'   marks the 1e10 sentinel only); with \code{opg_diag}, also
//'   \code{opg_diag}, as computed (not sanitized)
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
  // negated gradient and R's copies of it and of the diagonal. The pass
  // allocates every thread's sums and buffers on the primary thread, for the
  // team checked here.
  const int n_threads = situation_team(lay);
  kernel_check_memory(
      "gradient", n_threads,
      mnl_thread_bytes(lay, K, n_params, opg_diag ? 3.0 : 1.0) +
          sizeof(OpgDiagTerms) + 3.0 * sizeof(std::vector<double>) +
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

  // H3: Sanitize NaN/Inf (logit_sanitize(): the sentinel, zeroed entries,
  // and the `overflow` flag by which run_polish() fails a trial point; the
  // optimizers' path drops the flag, so they see the list as before); the
  // diagonal is returned as computed
  double obj = -global_loglik;
  arma::vec grad = -global_grad;
  const bool overflow = logit_sanitize(obj, grad, nonfinite_seen);
  if (opg_diag) {
    return Rcpp::List::create(
      Rcpp::Named("objective") = obj,
      Rcpp::Named("gradient")  = grad,
      Rcpp::Named("overflow")  = overflow,
      Rcpp::Named("opg_diag")  = global_diag
    );
  }
  return Rcpp::List::create(
    Rcpp::Named("objective") = obj,
    Rcpp::Named("gradient")  = grad,
    Rcpp::Named("overflow")  = overflow
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
  const int n_threads = situation_team(lay);
  kernel_check_memory("BHHH matrix", n_threads,
                      mnl_map_bytes(J_d, lay.max_m, K) + 8.0 * (K + a_max + n_params + 32.0) +
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
  MnlSituationMaps maps = mnl_situation_maps(J_d, lay.max_m, K, n_threads);
  std::vector<MnlWork> work = mnl_work(lay, K, n_threads);
  std::vector<std::vector<double>> scores =
      thread_vectors<double>(static_cast<std::size_t>(K + a_max), n_threads);
  std::vector<std::vector<double>> cols =
      thread_vectors<double>(static_cast<std::size_t>(n_params), n_threads);

#ifdef _OPENMP
  const int chunk = situation_chunk(lay);  // situations per dispatch
#pragma omp parallel num_threads(n_threads)
#endif
  {
    const int tid = kernel_thread_num();
    MnlSituationMap& mp = *maps[tid];
    MnlWork& w = work[tid];
    double* acc_t = acc.base(tid);
    double* s_buf = scores[tid].data();
    double* col = cols[tid].data();

#ifdef _OPENMP
#pragma omp for schedule(dynamic, chunk)
#endif
    for (choicer_off i = 0; i < N; ++i) {
      const double w_i = weights[i];
      mnl_situation_map(lay, i, use_asc, include_outside_option, mp);
      arma::vec s_c(s_buf, K + mp.a, false, true);
      mnl_situation_score(X, lay, base_util, i, include_outside_option, mp, w,
                          s_c);

      // w_i s_c s_c' over the compact block, column by column, each entry
      // formed as the dense rank-one update formed it (`B += w_i * s_i *
      // s_i.t()` is Armadillo's syrk_vec: the product of the two scores,
      // times the weight, then added). Each term is the rounded product,
      // never fused into its addition (CHOICER_STORE_TERMS: under GCC a
      // column's terms are stored before they are added), bit for bit under
      // clang and under GCC at R's default -O2 (GCC's -O3 loop vectorizer
      // rounds the dense update differently). The dense update added w_i *
      // (0 * s) outside the block, +-0, which changes nothing: no
      // accumulator, which starts at +0, is ever -0.
      if (unit_finite(w_i, s_c.memptr(), static_cast<int>(s_c.n_elem))) {
        switch (acc.mode()) {
          case SymMode::Single:
            bhhh_add<SymMode::Single>(acc, acc_t, mp, K, s_c, w_i, col);
            break;
          case SymMode::Triangles:
            bhhh_add<SymMode::Triangles>(acc, acc_t, mp, K, s_c, w_i, col);
            break;
          default:
            bhhh_add<SymMode::Shared>(acc, acc_t, mp, K, s_c, w_i, col);
        }
      } else {
        // A non-finite weight, or a non-finite score (utilities that
        // overflowed): the dense update over every entry, whose 0 * NaN and
        // 0 * Inf spread it over whole rows and columns. The scores come
        // from the slot map, zero outside the block.
        switch (acc.mode()) {
          case SymMode::Single:
            bhhh_add_dense<SymMode::Single>(acc, acc_t, mp, K, s_c, w_i, col);
            break;
          case SymMode::Triangles:
            bhhh_add_dense<SymMode::Triangles>(acc, acc_t, mp, K, s_c, w_i, col);
            break;
          default:
            bhhh_add_dense<SymMode::Shared>(acc, acc_t, mp, K, s_c, w_i, col);
        }
      }
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
                      mnl_map_bytes(J_d, lay.max_m, K) + 8.0 * (K + a_max + 32.0) +
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
  MnlSituationMaps maps = mnl_situation_maps(J_d, lay.max_m, K, n_threads);
  std::vector<MnlWork> work = mnl_work(lay, K, n_threads);
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
    MnlWork& w = work[tid];
    ClusterSum& cs = sums[tid];
    double* acc_t = acc.base(tid);
    double* s_buf = scores[tid].data();

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
        mnl_situation_score(X, lay, base_util, i, include_outside_option, mp, w,
                            s_c);
        cs.add(mp, K, s_c, weights[i]);
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
// log-likelihood contribution). Same loop body as the dense BHHH update had,
// with the outer product replaced by a row write. The variance routes no
// longer form it (mnl_bhhh_parallel and mnl_cluster_meat_parallel sum the
// scores in C++); at census scale it would hold 200,000 x 2,611 values
// (4.2 GB), so it stops above max_bytes (2 GiB by default) before allocating.
// [[Rcpp::export]]
arma::mat mnl_scores_parallel(
    const arma::vec& theta,
    const arma::mat& X,
    const Rcpp::IntegerVector& alt_idx,
    const Rcpp::IntegerVector& choice_idx,
    const Rcpp::IntegerVector& M,
    const bool use_asc = true,
    const bool include_outside_option = false,
    const double max_bytes = 2147483648.0
) {
  const int K = X.n_cols;
  const int n_params = theta.n_elem;

  const MnlParams par = parse_mnl_theta(theta, K, use_asc, include_outside_option);
  const ChoiceLayout lay = choice_layout_build(X, alt_idx, M, use_asc,
                                               par.delta, nullptr, &choice_idx);
  const choicer_off N = lay.N;

  // --- Serial pre-loop validation of chosen-alternative indices ---
  validate_choices(lay, include_outside_option, "mnl_scores_parallel");

  const double bytes = sizeof(double) * static_cast<double>(N) *
                       static_cast<double>(n_params);
  if (!(bytes <= max_bytes)) {
    Rcpp::stop("The score matrix would take %.2f GiB (%d choice situations x %d "
               "parameters), more than max_bytes (%.2f GiB). vcov(type = ), "
               "wesml_vcov() and the fit-time variances do not form it; pass "
               "a larger max_bytes to compute it.",
               bytes / 1073741824.0, N, n_params, max_bytes / 1073741824.0);
  }

  // Memory: each thread's score and situation buffers; X beta, the N x n
  // matrix and its copy into R
  const int n_threads = situation_team(lay);
  kernel_check_memory("score matrix", n_threads,
                      mnl_thread_bytes(lay, K, n_params, 1.0), [&](const int) {
                        return 8.0 * static_cast<double>(lay.n_rows) +
                               16.0 * static_cast<double>(N) * n_params;
                      });

  // Pre-compute base utility for all individuals (single BLAS call)
  arma::vec base_util = compute_base_util(X, par.beta, lay, use_asc, par.delta);

  // Output: one row per choice situation (each written by exactly one
  // iteration, so no accumulator or critical section is needed).
  arma::mat scores(static_cast<arma::uword>(N), n_params);

  // Every thread's working arrays and score, allocated here on the primary
  // thread
  std::vector<MnlWork> work = mnl_work(lay, K, n_threads);
  std::vector<std::vector<double>> score_t =
      thread_vectors<double>(static_cast<std::size_t>(n_params), n_threads);

#ifdef _OPENMP
  const int chunk = situation_chunk(lay);  // situations per dispatch
#pragma omp parallel num_threads(n_threads)
#endif
  {
    const int tid = kernel_thread_num();
    MnlWork& w = work[tid];
    arma::vec s_i(score_t[tid].data(), n_params, false, true); // per-individual score

#ifdef _OPENMP
#pragma omp for schedule(dynamic, chunk)
#endif
    for (choicer_off i = 0; i < N; ++i) {
      const int m_i         = lay.m(i);
      const int num_choices = include_outside_option ? m_i + 1 : m_i;
      const arma::uword start_idx = static_cast<arma::uword>(lay.row_off[i]);
      const AltCodes0 alt_idx0_i{lay.alt + lay.row_off[i]}; // M[i]

      // Probabilities
      mnl_softmax(base_util, lay, i, include_outside_option, w);

      // Chosen alternative (validated serially above)
      int chosen_alt = lay.choice[i];
      if (!include_outside_option) chosen_alt -= 1;

      // diff_vec[a] = 1{a == chosen_alt} - P_i[a]  (same as gradient kernel)
      mnl_diff(w, num_choices, chosen_alt);
      const arma::vec diff_vec(w.diff, num_choices, false, true);

      // Assemble the (weight-free) per-individual score s_i.
      s_i.zeros();

      // Beta block, X_i' diff into s_i's first K entries
      xt_times(X, start_idx, m_i, w.diff + (include_outside_option ? 1 : 0),
               false, 0.0, s_i.memptr(), w.xi);

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

  // Pre-compute base utility for all individuals (single BLAS call). Every
  // row belongs to one situation, so these are the utilities returned (the
  // situations' copies of them made the same vector).
  arma::vec base_util = compute_base_util(X, par.beta, lay, use_asc, par.delta);

  // Preallocate predicted values, and every thread's working arrays
  arma::vec P_all = arma::zeros(X.n_rows);
  const int n_threads = situation_team(lay);
  std::vector<MnlWork> work = mnl_work(lay, 0, n_threads);  // no products
  const int o = include_outside_option ? 1 : 0;  // slot of the first inside alt

#ifdef _OPENMP
  const int chunk = situation_chunk(lay);  // situations per dispatch
#pragma omp parallel num_threads(n_threads)
#endif
  {
    MnlWork& w = work[kernel_thread_num()];
#ifdef _OPENMP
#pragma omp for schedule(dynamic, chunk)
#endif
    for (choicer_off i = 0; i < N; ++i) {
      const int m_i = lay.m(i);
      const arma::uword start_idx = static_cast<arma::uword>(lay.row_off[i]);

      // Vector of choice probabilities
      mnl_softmax(base_util, lay, i, include_outside_option, w);
      std::copy(w.p + o, w.p + o + m_i, P_all.memptr() + start_idx);
    } // end of i loop
  }

  // Return predicted choice probabilities and utilities
  return Rcpp::List::create(
    Rcpp::Named("choice_prob") = P_all,
    Rcpp::Named("utility")  = base_util
  );
}

// The working arrays of mnl_predict_shares_internal(): every thread's
// situation arrays and shares, allocated on the primary thread once per call.
struct MnlSharesWork {
  std::vector<MnlWork> work;
  std::vector<std::vector<double>> local;  // each thread's shares
  MnlSharesWork(const ChoiceLayout& lay, const int num_alts)
      : work(mnl_work(lay, 0, situation_team(lay))),  // no products
        local(thread_vectors<double>(static_cast<std::size_t>(num_alts),
                                     static_cast<int>(work.size()))) {}
};

// Prediction of shares (internal function), from the base utilities of all
// rows (X * beta plus the rows' ASCs), into `shares` (num_alts entries).
// denominator is shares_denominator(weights), which the callers check before
// they form the utilities.
static void mnl_predict_shares_internal(
    const arma::vec& base_util,        // sum(M) x 1 vector of base utilities
    const ChoiceLayout& lay,           // layout of the stacked design
    const arma::vec& weights,          // N x 1 vector with weights for each observation
    const double denominator,          // sum of the weights
    const int num_alts,                // total number of distinct alternatives
    const bool include_outside_option, // whether to include an outside option
    MnlSharesWork& ws,                 // the threads' working arrays and shares
    arma::vec& shares                  // the result
) {
  const choicer_off N = lay.N;
  const int n_threads = static_cast<int>(ws.work.size());

  // Initialize global accumulator for predicted shares
  shares.zeros(num_alts);

#ifdef _OPENMP
  const int chunk = situation_chunk(lay);  // situations per dispatch
#pragma omp parallel num_threads(n_threads)
#endif
{
  const int tid = kernel_thread_num();
  MnlWork& w = ws.work[tid];
  // Thread-local accumulator for predicted shares
  arma::vec local_shares(ws.local[tid].data(), num_alts, false, true);
  local_shares.zeros();

  // Loop over individuals in parallel
#ifdef _OPENMP
#pragma omp for schedule(dynamic, chunk)
#endif
  for (choicer_off i = 0; i < N; ++i) {
    const int m_i         = lay.m(i);
    const AltCodes0 alt_idx0_i{lay.alt + lay.row_off[i]}; // M[i]

    // Vector of choice probabilities
    mnl_softmax(base_util, lay, i, include_outside_option, w);
    const double* P_i = w.p;

    // Sum probabilities by alternative
    if (include_outside_option) {
      // outside option at index 0
      local_shares(0) += weights[i] * P_i[0]; // outside option
    }
    for (int a = 0; a < m_i; ++a) {
      if (include_outside_option) {
        // outside option at index 0
        local_shares(alt_idx0_i[a] + 1) += weights[i] * P_i[a + 1];
      } else {
        local_shares(alt_idx0_i[a]) += weights[i] * P_i[a];
      }
    }

  } // end of i loop

#ifdef _OPENMP
#pragma omp critical
#endif
  {
    // Combine local sums into global sum
    shares += local_shares;
  }
}
  // The predicted shares
  shares /= denominator;
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
  MnlSharesWork ws(lay, num_alts);
  arma::vec global_shares(num_alts);
  mnl_predict_shares_internal(base_util, lay, weights, denominator, num_alts,
                              include_outside_option, ws, global_shares);

  return global_shares;
}

// ============================================================================
// The share inversion's model (blp_engine.h)
//
// One pass per evaluation of the shares: each situation's utilities, X beta
// plus the constants, gathered into the thread's MnlWork (the gather and the
// softmax's maximum in separate loops; interleaved they ran slower), its
// probabilities by stable_softmax_n(), then the shares in the statements of
// mnl_predict_shares_internal(), so that they are its shares at the same
// delta (bit for bit at one thread where the compiler fuses both statements
// alike, and the sums are not compensated), and A_jj = sum w P_j (1 - P_j) and C_jb = sum w P_j P_b, b the
// outside option or the situation's component's reference, in statements of
// their own (a product shared with the shares' would change how GCC fuses
// theirs). X beta is formed once per call and the threads' sums are added in
// thread order. Where a situation lists a code twice, A and C take the
// situation's total probability of each code (the derivative's exact terms).
//
// With the exp-space pass (blp_engine.h), the pass at an origin also keeps
// each row's exponential exp(v - max) before the division (and each
// situation's outside option's), and an evaluation near the origin forms the
// probabilities as those exponentials times exp(delta_j - delta0_j), over
// their sum, in the order of stable_softmax_n()'s sum and divisions.
// ============================================================================
struct MnlBlpModel {
  const ChoiceLayout& lay;
  const arma::vec& Xb;
  const arma::vec& weights;
  const double W;
  const bool outside;
  const int J;
  const BlpComponents& comp;
  std::vector<MnlWork> work;
  std::vector<std::vector<double>> acc_s;   // J shares, then each component's outside share
  std::vector<std::vector<double>> acc_ac;  // A and C of each alternative, interleaved
  std::vector<std::vector<int>> slot;       // repeated codes: alternative -> slot, or -1
  std::vector<std::vector<int>> code;       // and each slot's alternative
  std::vector<std::vector<double>> tot;     // and its probability
  BlpExpCache xc;                           // the exp-space cache's state
  std::unique_ptr<double[]> E, E0;          // its exponentials: each row's, and
                                            // each situation's outside option's
  bool kahan = false;                       // compensated share sums
  std::vector<std::vector<double>> acc_k;   // and acc_s's compensations

  MnlBlpModel(const ChoiceLayout& lay_, const arma::vec& Xb_, const arma::vec& w_,
              const double W_, const bool outside_, const int J_,
              const BlpComponents& comp_)
      : lay(lay_), Xb(Xb_), weights(w_), W(W_), outside(outside_), J(J_), comp(comp_),
        work(mnl_work(lay_, 0, situation_team(lay_))), xc(J_, outside_) {
    const int T = static_cast<int>(work.size());
    acc_s = thread_vectors<double>(static_cast<std::size_t>(J + comp.n), T);
    acc_ac = thread_vectors<double>(2 * static_cast<std::size_t>(J), T);
    if (comp.repeated) {
      const std::size_t m = static_cast<std::size_t>(lay.max_m);
      slot = thread_vectors<int>(static_cast<std::size_t>(J), T, -1);
      code = thread_vectors<int>(m, T);
      tot = thread_vectors<double>(m, T);
    }
  }

  // Keep the exp-space cache if its bytes (a double per row, and per
  // situation with an outside option) are within budget and can be had.
  bool keep_exp(const double budget) {
    const double n = static_cast<double>(lay.n_rows) +
                     (outside ? static_cast<double>(lay.N) : 0.0);
    if (lay.n_rows <= 0 || !(sizeof(double) * n <= budget)) return false;
    try {
      E.reset(new double[static_cast<std::size_t>(lay.n_rows)]);
      if (outside) E0.reset(new double[static_cast<std::size_t>(lay.N)]);
    } catch (const std::bad_alloc&) {
      E.reset();
      E0.reset();
      return false;
    }
    return true;
  }

  void pass(const arma::vec& delta, BlpPass& out) {
    if (E) {
      blp_exp_pass(*this, xc, delta, out);
    } else {
      sweep<false, false>(delta, out);
    }
  }

  // Form the exp-space cache at delta0 (tests: control$origin)
  void form(const arma::vec& delta0) {
    BlpPass scratch(J, comp.n);
    sweep<false, true>(delta0, scratch);
    xc.formed(delta0, scratch);
  }

  // Compensate the share sums (blp_engine.h)
  void compensate() {
    kahan = true;
    acc_k = thread_vectors<double>(static_cast<std::size_t>(J + comp.n),
                                   static_cast<int>(work.size()));
  }

  // One pass: the usual one (Exp false), storing the exp-space cache (Store),
  // or from the cache (Exp, its factors in xc.factor).
  template <bool Exp, bool Store>
  void sweep(const arma::vec& delta, BlpPass& out) {
    if (kahan) {
      sweep_k<Exp, Store, true>(delta, out);
    } else {
      sweep_k<Exp, Store, false>(delta, out);
    }
  }

  // ... with the share sums compensated (Kahan) or plain
  template <bool Exp, bool Store, bool Kahan>
  void sweep_k(const arma::vec& delta, BlpPass& out) {
    const choicer_off N = lay.N;
    const int T = static_cast<int>(work.size());
    const int o = outside ? 1 : 0;
    const arma::uword n_s = static_cast<arma::uword>(J + comp.n);
    const double* d = delta.memptr();
    const double* xb = Xb.memptr();
    double* e = E.get();
    double* e0 = E0.get();
    const double* fa = xc.factor.memptr();
    const bool repeated = comp.repeated;
    // the threads that ran: a smaller team than planned (OMP_DYNAMIC) leaves
    // the other buffers holding an earlier pass's sums
    int T_used = 1;
#ifdef _OPENMP
    const int chunk = situation_chunk(lay);  // situations per dispatch
#pragma omp parallel num_threads(T)
#endif
    {
      const int tid = kernel_thread_num();
#ifdef _OPENMP
      if (tid == 0) T_used = omp_get_num_threads();
#endif
      MnlWork& wk = work[tid];
      arma::vec local_shares(acc_s[tid].data(), n_s, false, true);
      local_shares.zeros();
      double* ls = acc_s[tid].data();
      double* lk = Kahan ? acc_k[tid].data() : nullptr;
      if (Kahan) std::fill(lk, lk + n_s, 0.0);
      double* ac = acc_ac[tid].data();
      std::fill(ac, ac + 2 * static_cast<std::size_t>(J), 0.0);
      int* slot_t = repeated ? slot[tid].data() : nullptr;
      int* code_t = repeated ? code[tid].data() : nullptr;
      double* tot_t = repeated ? tot[tid].data() : nullptr;
#ifdef _OPENMP
#pragma omp for schedule(dynamic, chunk)
#endif
      for (choicer_off i = 0; i < N; ++i) {
        const int m = lay.m(i);
        const choicer_off r0 = lay.row_off[i];
        const int* alt = lay.alt + r0;
        double* v = wk.v;
        if (Exp) {
          double* p = wk.p;
          if (outside) p[0] = e0[i];
          const double* er = e + r0;
          for (int a = 0; a < m; ++a) p[o + a] = er[a] * fa[alt[a] - 1];
          softmax_normalize_n(p, m + o);
        } else {
          if (outside) v[0] = 0.0;
          for (int a = 0; a < m; ++a) v[o + a] = xb[r0 + a] + d[alt[a] - 1];
          if (Store) {
            double* ek = wk.diff;  // the exponentials, the outside option's first
            stable_softmax_keep_n(v, wk.p, m + o, ek);
            std::copy(ek + o, ek + o + m, e + r0);
            if (outside) e0[i] = ek[0];
          } else {
            stable_softmax_n(v, wk.p, m + o);
          }
        }
        const double* P_i = wk.p;
        const double w_i = weights[i];
        const int c = comp.of[alt[0] - 1];
        double pb = 0.0;  // the base's probability
        if (outside) {
          if (Kahan) {
            kahan_add(ls[J + c], lk[J + c], w_i * P_i[0]);
          } else {
            local_shares(J + c) += w_i * P_i[0];
          }
          pb = P_i[0];
        } else {
          const int b = comp.ref[c] + 1;
          for (int a = 0; a < m; ++a) {
            if (alt[a] == b) pb += P_i[a];
          }
        }
        for (int a = 0; a < m; ++a) {
          if (Kahan) {
            kahan_add(ls[alt[a] - 1], lk[alt[a] - 1], w_i * P_i[o + a]);
          } else {
            local_shares(alt[a] - 1) += w_i * P_i[o + a];
          }
        }
        if (!repeated) {
          for (int a = 0; a < m; ++a) {
            const int j = alt[a] - 1;
            const double p = P_i[o + a];
            ac[2 * j] += w_i * (p * (1.0 - p));
            ac[2 * j + 1] += w_i * (p * pb);
          }
        } else {
          int n_d = 0;
          for (int a = 0; a < m; ++a) {
            const int j = alt[a] - 1;
            if (slot_t[j] < 0) {
              slot_t[j] = n_d;
              code_t[n_d] = j;
              tot_t[n_d++] = 0.0;
            }
            tot_t[slot_t[j]] += P_i[o + a];
          }
          for (int k = 0; k < n_d; ++k) {
            const int j = code_t[k];
            const double p = tot_t[k];
            ac[2 * j] += w_i * (p * (1.0 - p));
            ac[2 * j + 1] += w_i * (p * pb);
            slot_t[j] = -1;
          }
        }
      }
    }
    for (int j = 0; j < J; ++j) {
      double s = 0.0, A = 0.0, C = 0.0, sk = 0.0;
      for (int t = 0; t < T_used; ++t) {
        if (Kahan) {
          kahan_add(s, sk, acc_s[t][j] - acc_k[t][j]);
        } else {
          s += acc_s[t][j];
        }
        A += acc_ac[t][2 * j];
        C += acc_ac[t][2 * j + 1];
      }
      if (Kahan) s -= sk;
      out.s[j] = s / W;
      out.A[j] = A / W;
      out.C[j] = C / W;
    }
    for (int k = 0; k < comp.n; ++k) {
      double s0 = 0.0, sk = 0.0;
      for (int t = 0; t < T_used; ++t) {
        if (Kahan) {
          kahan_add(s0, sk, acc_s[t][J + k] - acc_k[t][J + k]);
        } else {
          s0 += acc_s[t][J + k];
        }
      }
      if (Kahan) s0 -= sk;
      out.s0[k] = s0 / W;
    }
  }
};

//' The share inversion of the multinomial logit
//'
//' The engine behind blp_contraction() and blp() (blp_engine.h): the input
//' checks of blp_contraction(), in their order, then the share inversion's own
//' (a finite start, finite and non-negative weights with a finite sum, the
//' targets' sums and bounds), then the solver.
//'
//' @inheritParams blp_contraction
//' @param control List from .blp_control(): iteration, map, max_evals, tol,
//'   damping, patience, window, Anderson acceleration's memory,
//'   max_condition, restart_patience, restart_progress and restart_best,
//'   SQUAREM's scheme, step_min, step_max and step_factor, keep_pass, and the
//'   passes' arithmetic: pass ("log" or "exp"), cache_bytes and origin
//'   (blp_engine.h's BlpExpCache), and compensate ("off", "on" or "auto": the
//'   share sums by Kahan's compensated summation).
//' @return List: delta (n x 1, as blp_contraction() returns it), code (0
//'   converged, 1 max_iter, 2 not finite after backtracking, 3 not finite at the
//'   start), iterations, evaluations, backtracks, stalls, race (the evaluations
//'   before Berry's contraction joined; 0 if it did not), resumes (the times
//'   it went on after stopping as behind), restarts (Anderson acceleration's
//'   restarts; SQUAREM's extrapolated points that failed),
//'   berry (whether the point returned is its), residual, fallbacks; with
//'   keep_pass, pass (s, A, C, s0 of the last evaluation); with pass "exp",
//'   exp_kept (whether the cache was allocated), exp_evals (evaluations from
//'   it), exp_forms (usual passes that formed it) and exp_floors (evaluations
//'   from it that the floor sent to the usual pass); with compensate other
//'   than "off", compensated (whether the passes compensated their sums).
//' @noRd
// [[Rcpp::export]]
Rcpp::List mnl_blp_solve(
  const arma::vec& delta,
  const arma::vec& target_shares,
  const arma::mat& X,
  const arma::vec& beta,
  const Rcpp::IntegerVector& alt_idx,
  const Rcpp::IntegerVector& M,
  const arma::vec& weights,
  const bool include_outside_option,
  const Rcpp::List& control
) {
  const BlpControl ctl = blp_control(control);
  const bool use_asc = true;

  // total number of distinct alternatives
  const int J = static_cast<int>(delta.n_elem);
  const int num_alts = include_outside_option ? J + 1 : J;
  if (target_shares.n_elem != static_cast<arma::uword>(num_alts)) {
    Rcpp::stop("Error: target_shares must have the same length as the total number of alternatives.");
  }
  for (arma::uword k = 0; k < target_shares.n_elem; ++k) {
    if (!(target_shares[k] > 0)) {
      Rcpp::stop("Error: all target_shares must be strictly positive (log(share) is undefined otherwise).");
    }
  }
  const ChoiceLayout lay = choice_layout_build(X, alt_idx, M, use_asc, delta,
                                               &weights);
  const bool kahan = blp_compensate(ctl, static_cast<double>(lay.N), situation_team(lay));
  const bool kahan_w = blp_compensate_serial(ctl, static_cast<double>(lay.N));
  const double W = kahan_w ? blp_weight_sum(weights) : shares_denominator(weights);
  const arma::vec Xb = X * beta;  // once: only delta changes across passes

  blp_check_start(delta, weights);
  const BlpComponents comp = blp_components(lay, weights.memptr(), W, J,
                                            include_outside_option, kahan_w);
  const BlpProblem pb = blp_problem(comp, target_shares, include_outside_option, J);
  arma::vec d = delta;
  if (!include_outside_option && J > 0) d -= delta[0];
  MnlBlpModel model(lay, Xb, weights, W, include_outside_option, J, comp);
  if (kahan) model.compensate();
  bool kept = false;
  if (ctl.pass == BLP_PASS_EXP && !pb.free.empty() &&
      (ctl.max_evals > 1 || (ctl.origin.n_elem > 0 && ctl.max_evals > 0))) {
    // (the call holds the design and X beta)
    const double held = sizeof(double) * (static_cast<double>(X.n_elem) +
                                          static_cast<double>(Xb.n_elem));
    kept = model.keep_exp(blp_cache_budget(ctl.cache_bytes, held));
    if (kept && ctl.origin.n_elem > 0) model.form(blp_origin(ctl, J, include_outside_option));
  }
  BlpPass last(J, comp.n);
  const BlpResult res = blp_solve(model, comp, pb, d, ctl, &last);
  Rcpp::List out = blp_result(d, res, ctl, last);
  blp_exp_record(out, ctl, kept, model.xc);
  blp_comp_record(out, ctl, kahan);
  return out;
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
template <SymMode Mode>
inline void mnl_hessian_add(const SymAcc& acc, double* buf,
                            const MnlSituationMap& mp, const int K,
                            const int J_d, const double* BB,
                            const double* mu_b, const double* BD,
                            const double* mu_d, const double w,
                            char* nan_row, double* ct) {
  // Each term is the rounded product, never fused into its addition: where
  // the compiler contracts across statements a column's terms are stored in
  // ct first (CHOICER_STORE_TERMS).
  const std::size_t Ks = static_cast<std::size_t>(K);
  for (int c = 0; c < K; ++c) {
    const SymCol<Mode> col = acc.col<Mode>(buf, c);
    if (CHOICER_STORE_TERMS) {
      for (int r = 0; r <= c; ++r) {
        double v = BB[r + c * Ks];
        v *= w;
        ct[r] = v;
      }
      for (int r = 0; r <= c; ++r) col.add(static_cast<std::size_t>(r), -ct[r]);
    } else {
      for (int r = 0; r <= c; ++r) {
        double v = BB[r + c * Ks];
        v *= w;
        col.add(static_cast<std::size_t>(r), -v);
      }
    }
  }
  for (int s = 0; s < mp.a; ++s) {
    const int g = K + mp.delta[s];
    const SymCol<Mode> col = acc.col<Mode>(buf, g);
    const double ms = mu_d[s];
    if (CHOICER_STORE_TERMS) {
      for (int r = 0; r < K; ++r) {
        double v = BD[r + s * Ks];
        v *= w;
        ct[r] = v;
      }
      for (int t = 0; t < s; ++t) {
        double v = 0.0 - mu_d[t] * ms;
        v *= w;
        ct[K + t] = v;
      }
      double v = ms - ms * ms;
      v *= w;
      ct[K + s] = v;
      for (int r = 0; r < K; ++r) col.add(static_cast<std::size_t>(r), -ct[r]);
      for (int t = 0; t < s; ++t) {
        col.add(static_cast<std::size_t>(K + mp.delta[t]), -ct[K + t]);
      }
      col.add(static_cast<std::size_t>(g), -ct[K + s]);
    } else {
      for (int r = 0; r < K; ++r) {
        double v = BD[r + s * Ks];
        v *= w;
        col.add(static_cast<std::size_t>(r), -v);
      }
      for (int t = 0; t < s; ++t) {
        double v = 0.0 - mu_d[t] * ms;
        v *= w;
        col.add(static_cast<std::size_t>(K + mp.delta[t]), -v);
      }
      double v = ms - ms * ms;
      v *= w;
      col.add(static_cast<std::size_t>(g), -v);
    }
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
      acc.add<Mode>(buf, std::min(g, h), std::max(g, h), -v);
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
  const int n_threads = situation_team(lay);
  kernel_check_memory("Hessian", n_threads,
                      mnl_map_bytes(J_d, lay.max_m, K) +
                          8.0 * (K * (K + 2.0) + (K + 1.0) * a_max + (K + a_max + 1.0) + 64.0) +
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
  MnlSituationMaps maps = mnl_situation_maps(J_d, lay.max_m, K, n_threads);
  std::vector<MnlWork> work = mnl_work(lay, K, n_threads);
  const std::size_t blk = static_cast<std::size_t>(K * (K + 2.0) + (K + 1.0) * a_max);
  const std::size_t n_ct = static_cast<std::size_t>(K + a_max + 1.0);
  std::vector<std::vector<double>> blocks = thread_vectors<double>(blk + n_ct, n_threads);
  std::vector<NonFiniteRecord> nonfinite = nonfinite_records(n_params, n_threads);
  std::vector<std::vector<char>> nan_rows =
      thread_vectors<char>(static_cast<std::size_t>(n_params), n_threads, 0);
  const int o = include_outside_option ? 1 : 0;  // slot of the first inside alt
  const double* Xm = X.memptr();
  const arma::uword n_rows = X.n_rows;

#ifdef _OPENMP
  const int chunk = situation_chunk(lay);  // situations per dispatch
#pragma omp parallel num_threads(n_threads)
#endif
  {
    const int tid = kernel_thread_num();
    MnlSituationMap& mp = *maps[tid];
    MnlWork& w = work[tid];
    NonFiniteRecord& nf = nonfinite[tid];
    char* nan_row = nan_rows[tid].data();
    double* acc_t = acc.base(tid);
    // BB (K x K), mu_b (K), BD (K x a), mu_d (a), K products and a column's
    // terms (mnl_hessian_add()), in one buffer
    const std::size_t Ks = static_cast<std::size_t>(K);
    double* BB = blocks[tid].data();
    double* mu_b = BB + Ks * Ks;
    double* BD = mu_b + K;
    double* ct = blocks[tid].data() + blk;

#ifdef _OPENMP
#pragma omp for schedule(dynamic, chunk)
#endif
    for (choicer_off i = 0; i < N; ++i) {
      const int m_i         = lay.m(i);
      const arma::uword start_idx = static_cast<arma::uword>(lay.row_off[i]);
      const double w_i      = weights[i];

      // Calculate Probabilities ---------------------------------------------
      mnl_softmax(base_util, lay, i, include_outside_option, w);
      const double* P_i = w.p;

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

      switch (acc.mode()) {
        case SymMode::Single:
          mnl_hessian_add<SymMode::Single>(acc, acc_t, mp, K, J_d, BB, mu_b, BD, mu_d, w_i,
                                           nan_row, ct);
          break;
        case SymMode::Triangles:
          mnl_hessian_add<SymMode::Triangles>(acc, acc_t, mp, K, J_d, BB, mu_b, BD, mu_d, w_i,
                                           nan_row, ct);
          break;
        default:
          mnl_hessian_add<SymMode::Shared>(acc, acc_t, mp, K, J_d, BB, mu_b, BD, mu_d, w_i,
                                           nan_row, ct);
      }
    } // end of i loop
  } // end parallel region

  // The information (Hessian of the negated log-likelihood), with the dense
  // code's spread of non-finite weights
  acc.finish(true);
  apply_nonfinite(acc, nonfinite);
  return acc.result;
}

// =============================================================================
// Elasticities and diversion ratios
//
// Each kernel adds every situation's terms, by alternative, into a J x J sum
// (DenseAcc in choicer_derivs.h: the result itself with one thread, the
// threads' copies while they fit in acc_bytes, else the shared result) and
// returns that R matrix, divided in place. The passes over the situations
// are templates on whether the sum is shared, with each statement of the
// code they replace: a situation's probabilities in the thread's MnlWork, its
// alternatives' rows and columns in the thread's map, x_k read from X's
// column in place.
// =============================================================================

// Every thread's map of a situation's choices to the rows and columns of the
// J x J sums, on the primary thread.
static inline std::vector<std::vector<int>> mnl_global_maps(const ChoiceLayout& lay,
                                                            const int n_threads) {
  return thread_vectors<int>(static_cast<std::size_t>(lay.max_m) + 1, n_threads);
}

// The memory of an elasticity or diversion-ratio kernel with t threads, an
// upper bound: per thread its MnlWork, its map and n_vecs padded vectors of
// J_total (its sums); shared, X beta, the J x J sums (dense_acc_need()), a
// vector of J_total (the summed denominators) and the threads' totals.
static inline void mnl_dense_check_memory(const char* what, const ChoiceLayout& lay,
                                          const int n_threads, const int J_total,
                                          const double n_vecs, const bool outside,
                                          const double acc_bytes) {
  const std::size_t J_n = static_cast<std::size_t>(J_total);
  kernel_check_memory(
      what, n_threads,
      mnl_thread_bytes(lay, 0, J_total, n_vecs) +
          sizeof(int) * (static_cast<double>(lay.max_m) + 1.0 + 32.0),
      [&](const int t) {
        return 8.0 * (static_cast<double>(lay.n_rows) + J_total + t) +
               dense_acc_need(J_n, t, acc_bytes, outside);
      });
}

// The elasticity kernel's pass over the choice situations: each one's
// weighted elasticities added into acc, its weight into its thread's total.
template <bool Shared>
static void mnl_elasticity_pass(const arma::mat& X, const ChoiceLayout& lay,
                                const arma::vec& base_util,
                                const arma::vec& weights, const int var_idx,
                                const double beta_k,
                                const bool include_outside_option,
                                std::vector<MnlWork>& work,
                                std::vector<std::vector<int>>& maps,
                                const DenseAcc& acc,
                                std::vector<double>& total_weight) {
  const choicer_off N = lay.N;
  const int o = include_outside_option ? 1 : 0;  // slot of the first inside alt

#ifdef _OPENMP
  const int chunk = situation_chunk(lay);  // situations per dispatch
#pragma omp parallel num_threads(acc.T)
#endif
  {
    const int tid = kernel_thread_num();
    MnlWork& w = work[tid];
    int* global_j_map = maps[tid].data();
    double* buf = acc.base(tid);
    double local_total_weight = 0.0;

#ifdef _OPENMP
#pragma omp for schedule(dynamic, chunk)
#endif
    for (choicer_off i = 0; i < N; ++i) {
      const int m_i         = lay.m(i);
      const int num_choices = m_i + o;
      const double w_i      = weights[i];

      // --- 2. Compute Probabilities (same as loglik) ---
      mnl_softmax(base_util, lay, i, include_outside_option, w);
      const double* P_i = w.p;

      // --- 3. Prepare Data for Elasticity Calculation ---

      // The values of the variable k for this individual's choice set, in
      // place (x_k for the outside option is 0)
      const double* x_k_i = X.colptr(var_idx) + lay.row_off[i];

      // Map local choice set indices (0...num_choices-1) to global
      // alternative indices (0...J_total-1)
      fill_global_alt_map(global_j_map, lay.alt + lay.row_off[i], m_i,
                          include_outside_option);

      // --- 4. Compute Individual Elasticity Matrix ---
      // E_n(i, j) = elasticity of P_ni w.r.t attribute x_njk
      for (int i_local = 0; i_local < num_choices; ++i_local) {
        const int global_i = global_j_map[i_local]; // Global row index
        const double P_ni = P_i[i_local];

        for (int j_local = 0; j_local < num_choices; ++j_local) {
          const int global_j = global_j_map[j_local]; // Global col index
          const double P_nj = P_i[j_local];
          const double x_njk = j_local < o ? 0.0 : x_k_i[j_local - o];

          double elasticity;
          if (global_i == global_j) {
            // Own-elasticity: beta_k * x_nik * (1 - P_ni)
            elasticity = beta_k * x_njk * (1.0 - P_ni);
          } else {
            // Cross-elasticity: -beta_k * x_njk * P_nj
            elasticity = -beta_k * x_njk * P_nj;
          }

          // Accumulate the weighted elasticity
          acc.add<Shared>(buf, static_cast<std::size_t>(global_i),
                          static_cast<std::size_t>(global_j), w_i, elasticity);
        }
      }
      local_total_weight += w_i;
    } // end of i loop
    total_weight[tid] = local_total_weight;
  } // end parallel region
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
//' @param acc_bytes Bytes for the threads' copies of the J x J sums (2 GiB by
//'   default): with more threads than one, each adds into its own copy while
//'   the copies fit, and into the shared result otherwise (its last bits then
//'   depend on the threads' timing); \code{0} forces the shared result (for
//'   tests).
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
Rcpp::NumericMatrix mnl_elasticities_parallel(
    const arma::vec& theta,
    const arma::mat& X,
    const Rcpp::IntegerVector& alt_idx,
    const Rcpp::IntegerVector& choice_idx, // Kept for consistency, but not used
    const Rcpp::IntegerVector& M,
    const arma::vec& weights,
    const int elast_var_idx,
    const bool use_asc = true,
    const bool include_outside_option = false,
    const double acc_bytes = 2147483648.0
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
  const double beta_k = par.beta(var_idx); // The coefficient for our variable

  // Determine total number of alternatives for the output matrix
  // J_inside = number of alternatives (with or without normalization)
  // J_total = total size of matrix (J_inside + 1 if outside option)
  const int J_inside = compute_J_inside(use_asc, par.delta, lay);
  const int J_total = compute_J_total(J_inside, include_outside_option);
  check_alternative_count(J_total);

  // Memory: each thread's working arrays and map; X beta and the J x J sums
  const int n_threads = situation_team(lay);
  mnl_dense_check_memory("elasticity matrix", lay, n_threads, J_total, 0.0,
                         include_outside_option, acc_bytes);

  // The J x J sums first (an R allocation), then X beta
  DenseAcc acc(J_total, n_threads, include_outside_option, "elasticity",
               acc_bytes);

  // Pre-compute base utility for all individuals (single BLAS call)
  arma::vec base_util = compute_base_util(X, par.beta, lay, use_asc, par.delta);

  // Every thread's working arrays, map and total weight
  std::vector<MnlWork> work = mnl_work(lay, 0, n_threads);
  std::vector<std::vector<int>> maps = mnl_global_maps(lay, n_threads);
  std::vector<double> total_weight(static_cast<std::size_t>(n_threads), 0.0);
  if (acc.shared()) {
    mnl_elasticity_pass<true>(X, lay, base_util, weights, var_idx, beta_k,
                              include_outside_option, work, maps, acc,
                              total_weight);
  } else {
    mnl_elasticity_pass<false>(X, lay, base_util, weights, var_idx, beta_k,
                               include_outside_option, work, maps, acc,
                               total_weight);
  }
  acc.finish();
  double global_total_weight = 0.0;
  for (int t = 0; t < n_threads; ++t) global_total_weight += total_weight[t];

  // --- 5. Finalize ---
  // Compute the average by dividing by the total weight
  if (global_total_weight > 1e-10) {
    const std::size_t nn = acc.n * acc.n;
    for (std::size_t e = 0; e < nn; ++e) acc.out[e] /= global_total_weight;
  }

  return acc.result;
}

// The diversion-ratio kernel's pass over the choice situations: each one's
// weighted numerators added into acc, its denominators into its thread's.
template <bool Shared>
static void mnl_diversion_pass(const ChoiceLayout& lay, const arma::vec& base_util,
                               const arma::vec& weights,
                               const bool include_outside_option,
                               std::vector<MnlWork>& work,
                               std::vector<std::vector<int>>& maps,
                               const DenseAcc& acc,
                               std::vector<std::vector<double>>& dens) {
  const choicer_off N = lay.N;

#ifdef _OPENMP
  const int chunk = situation_chunk(lay);  // situations per dispatch
#pragma omp parallel num_threads(acc.T)
#endif
  {
    const int tid = kernel_thread_num();
    MnlWork& w = work[tid];
    int* global_j_map = maps[tid].data();
    double* buf = acc.base(tid);
    double* local_denominator = dens[tid].data();

#ifdef _OPENMP
#pragma omp for schedule(dynamic, chunk)
#endif
    for (choicer_off i = 0; i < N; ++i) {
      const int m_i         = lay.m(i);
      const int num_choices = include_outside_option ? m_i + 1 : m_i;
      const double w_i      = weights[i];

      // --- 2. Compute Probabilities ---
      mnl_softmax(base_util, lay, i, include_outside_option, w);
      const double* P_i = w.p;

      // --- 3. Map local indices to global alternative indices ---
      fill_global_alt_map(global_j_map, lay.alt + lay.row_off[i], m_i,
                          include_outside_option);

      // --- 4. Accumulate numerator and denominator ---
      for (int j_local = 0; j_local < num_choices; ++j_local) {
        const int global_j = global_j_map[j_local];
        const double P_nj = P_i[j_local];

        // Denominator: w_i * P_nj * (1 - P_nj)
        local_denominator[global_j] += w_i * P_nj * (1.0 - P_nj);

        // Numerator: w_i * P_nj * P_nk for k != j
        for (int k_local = 0; k_local < num_choices; ++k_local) {
          if (k_local == j_local) continue;
          const int global_k = global_j_map[k_local];
          const double P_nk = P_i[k_local];
          acc.add<Shared>(buf, static_cast<std::size_t>(global_k),
                          static_cast<std::size_t>(global_j), w_i * P_nj, P_nk);
        }
      }
    } // end of i loop
  } // end parallel region
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
//' @param acc_bytes Bytes for the threads' copies of the J x J numerators (2
//'   GiB by default), as in \code{mnl_elasticities_parallel()}; \code{0}
//'   forces the shared result (for tests).
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
Rcpp::NumericMatrix mnl_diversion_ratios_parallel(
    const arma::vec& theta,
    const arma::mat& X,
    const Rcpp::IntegerVector& alt_idx,
    const Rcpp::IntegerVector& M,
    const arma::vec& weights,
    const bool use_asc = true,
    const bool include_outside_option = false,
    const double acc_bytes = 2147483648.0
) {
  // --- 1. Parameter and Variable Setup ---
  const int K = X.n_cols;

  // Extract beta and delta (same logic as elasticities/loglik)
  const MnlParams par = parse_mnl_theta(theta, K, use_asc, include_outside_option);
  const ChoiceLayout lay = choice_layout_build(X, alt_idx, M, use_asc, par.delta,
                                               &weights);

  // Determine total number of alternatives for the output matrix
  const int J_inside = compute_J_inside(use_asc, par.delta, lay);
  const int J_total = compute_J_total(J_inside, include_outside_option);
  check_alternative_count(J_total);

  // Memory: each thread's working arrays, map and denominators; X beta, the
  // J x J numerators and the denominators' sum
  const int n_threads = situation_team(lay);
  mnl_dense_check_memory("diversion-ratio matrix", lay, n_threads, J_total, 1.0,
                         include_outside_option, acc_bytes);

  // The J x J sums first (an R allocation), then X beta
  DenseAcc acc(J_total, n_threads, include_outside_option, "diversion-ratio",
               acc_bytes);

  // Pre-compute base utility for all individuals (single BLAS call)
  arma::vec base_util = compute_base_util(X, par.beta, lay, use_asc, par.delta);

  // Accumulators
  // numerator(k, j) = sum_n w_n * P_nj * P_nk  (for k != j), into acc
  // denominator(j) = sum_n w_n * P_nj * (1 - P_nj), by thread
  const std::size_t J_n = static_cast<std::size_t>(J_total);
  std::vector<MnlWork> work = mnl_work(lay, 0, n_threads);
  std::vector<std::vector<int>> maps = mnl_global_maps(lay, n_threads);
  std::vector<std::vector<double>> dens = thread_vectors<double>(J_n, n_threads, 0.0);
  if (acc.shared()) {
    mnl_diversion_pass<true>(lay, base_util, weights, include_outside_option,
                             work, maps, acc, dens);
  } else {
    mnl_diversion_pass<false>(lay, base_util, weights, include_outside_option,
                              work, maps, acc, dens);
  }
  acc.finish();
  std::vector<double> global_denominator(J_n, 0.0);
  for (int t = 0; t < n_threads; ++t) {
    for (std::size_t j = 0; j < J_n; ++j) global_denominator[j] += dens[t][j];
  }

  // --- 5. Compute diversion ratios ---
  dense_diversion_ratios(acc.out, J_n, global_denominator,
                         [](const double d) { return d > 1e-15; });
  return acc.result;
}
