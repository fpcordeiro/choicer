// [[Rcpp::depends(RcppArmadillo)]]
#include "choicer.h"
#include "choicer_internal.h"
#include "choicer_derivs.h"

// Forward declaration: shared NL parameter parsing (defined further below,
// after the likelihood kernel; see the full doc comment at the definition).
static void nl_parse_theta(
    const arma::vec& theta,
    const Rcpp::IntegerVector& nest_idx,
    const int K,
    const bool use_asc,
    const bool include_outside_option,
    arma::vec& beta,
    arma::vec& lambda,
    arma::vec& delta,
    int& delta_start_idx,
    int& delta_length,
    arma::ivec* nest_k_to_theta_idx = nullptr);

// The number of nests: the largest nest code, as arma::max() gave it on the
// arma::uvec argument this replaces, which threw std::logic_error on an empty
// nest_idx (kept). Codes below 1 are reported by nl_parse_theta() and
// nl_blp_contraction().
static int nl_n_nests(const Rcpp::IntegerVector& nest_idx) {
  if (nest_idx.size() == 0) {
    throw std::logic_error("max(): object has no elements");
  }
  int n_nests = std::numeric_limits<int>::min();
  for (const int k : nest_idx) n_nests = std::max(n_nests, k);
  return n_nests;
}

// The 0-based nest of each of a situation's m alternatives, from their
// 1-based codes alt and the 1-based nest codes nest. The kernels index lambda
// and their per-nest buffers with the result unchecked, so two checks must
// come first: nl_layout_build() (every alternative code at most nest_idx's
// length) and nl_parse_theta() or nl_blp_contraction()'s own check (every
// nest code at least 1).
static inline void nl_gather_nests(int* nest0, const int* alt, const int* nest,
                                   const int m) {
  for (int j = 0; j < m; ++j) nest0[j] = nest[alt[j] - 1] - 1;
}

// Memory of the NL kernels' loop over the choice situations, an upper bound
// per thread: the probability buffers (NlProbs: five of max_m and six of
// n_nests doubles), the thread's copy of the nest map and the per-nest sums
// (two of n_nests), the score differences (max_m), each row's nest (max_m
// ints), the copy of X_i that Armadillo makes for X_i.t() * v (max_m x K) and
// its product (K), and n_vecs vectors of the n parameters (the thread's sums
// or scores).
static inline double nl_thread_bytes(const ChoiceLayout& lay, const int K,
                                     const int n_nests, const int n_params,
                                     const double n_vecs) {
  const double m = static_cast<double>(lay.max_m);
  return sizeof(double) * (6.0 * m + 8.0 * n_nests + m * K + K +
                           n_vecs * static_cast<double>(n_params)) +
         sizeof(int) * m;
}

// ============================================================================
// Compact blocks of the NL derivative kernels
//
// A situation's score and Hessian live on its compact block: the K
// coefficients; then the dissimilarities of the non-singleton nests it offers
// whose inclusive value log I_k is finite (the lambda entries the dense
// kernels touched; the others' scores are 0); then the free constants of the
// alternatives it offers. Each range is in increasing order, so the block's
// global order is increasing (theta holds beta, then the dissimilarities in
// nest order, then the constants). A situation listing an alternative twice
// adds both rows into one constant slot, in row order, as
// scatter_delta_grad() adds them into one constant. NlSituationMap holds the
// slots, one per thread, sized on the primary thread for the largest
// situation, so that building them allocates nothing; the markers are -1
// except while a situation is mapped.
// ============================================================================
struct NlSituationMap {
  std::vector<int> slot_of;      // free constant j -> slot, or -1 (J_d entries)
  std::vector<int> delta;        // constant slot -> free constant, increasing
  std::vector<int> row_slot;     // row of the situation -> constant slot, or -1
  std::vector<int> lam_slot_of;  // dissimilarity d -> slot, or -1 (K_l entries)
  std::vector<int> lam;          // lambda slot -> dissimilarity, increasing
  int a = 0;                     // constant slots
  int n_lam = 0;                 // lambda slots
  int K_l = 0;                   // dissimilarities in theta
  char pad[128];                 // keeps the next thread's map off these lines

  // Each array is reserved 32 entries (128 bytes) beyond its largest size,
  // so that no two threads' arrays share a cache line.
  NlSituationMap(const int J_d, const int K_l_, const choicer_off max_m)
      : slot_of(static_cast<std::size_t>(J_d) + 32, -1),
        lam_slot_of(static_cast<std::size_t>(K_l_) + 32, -1), K_l(K_l_) {
    delta.reserve(static_cast<std::size_t>(std::min<choicer_off>(J_d, max_m)) + 32);
    row_slot.reserve(static_cast<std::size_t>(max_m) + 32);
    lam.reserve(static_cast<std::size_t>(std::min<choicer_off>(K_l_, max_m)) + 32);
  }
  // Global parameter of local index i (the K coefficients first).
  int global(const int i, const int K) const {
    if (i < K) return i;
    if (i < K + n_lam) return K + lam[i - K];
    return K + K_l + delta[i - K - n_lam];
  }
  // Local index of global parameter g, or -1 outside the situation's block.
  int local(const int g, const int K) const {
    if (g < K) return g;
    if (g < K + K_l) {
      const int k = lam_slot_of[g - K];
      return k >= 0 ? K + k : -1;
    }
    const int k = slot_of[g - K - K_l];
    return k >= 0 ? K + n_lam + k : -1;
  }
};

// Bytes of one thread's NlSituationMap.
static inline double nl_map_bytes(const int J_d, const int K_l,
                                  const choicer_off max_m) {
  const double m = static_cast<double>(max_m);
  return sizeof(int) * (J_d + std::min<double>(J_d, m) + m + K_l +
                        std::min<double>(K_l, m) + 160.0) + 256.0;
}

// Every thread's map, each in its own allocation, on the primary thread.
using NlSituationMaps = std::vector<std::unique_ptr<NlSituationMap>>;
static inline NlSituationMaps nl_situation_maps(const int J_d, const int K_l,
                                                const choicer_off max_m,
                                                const int n_threads) {
  NlSituationMaps maps;
  maps.reserve(n_threads);
  for (int i = 0; i < n_threads; ++i) {
    maps.emplace_back(new NlSituationMap(J_d, K_l, max_m));
  }
  return maps;
}

// Every thread's probability buffers, nests of a situation's rows and per-nest
// sums, on the primary thread.
struct NlThreadBuffers {
  NlProbs pr;
  std::vector<int> nest0;        // each row's 0-based nest
  std::vector<double> sum_P_V;   // sum_{j in B_k} P(j|k) V_j, zero between situations
  NlThreadBuffers(const choicer_off max_m, const int n_nests)
      : pr(max_m, n_nests), nest0(static_cast<std::size_t>(max_m) + 32),
        sum_P_V(static_cast<std::size_t>(n_nests) + 16, 0.0) {}
};
static inline std::vector<NlThreadBuffers> nl_thread_buffers(const choicer_off max_m,
                                                             const int n_nests,
                                                             const int n_threads) {
  std::vector<NlThreadBuffers> buf;
  buf.reserve(n_threads);
  for (int i = 0; i < n_threads; ++i) buf.emplace_back(max_m, n_nests);
  return buf;
}

// Each thread's buffer of n doubles, on the primary thread.
static inline std::vector<std::vector<double>> nl_thread_vectors(const std::size_t n,
                                                                 const int n_threads) {
  return std::vector<std::vector<double>>(static_cast<std::size_t>(n_threads),
                                          std::vector<double>(n + 16));
}

// Free constant of a row with 0-based alternative alt0 (as
// scatter_delta_grad() places it), or -1.
static inline int nl_free_delta(const int alt0, const bool has_delta,
                                const bool include_outside_option) {
  if (!has_delta) return -1;
  return include_outside_option ? alt0 : alt0 - 1;
}

// The slots of situation i: its distinct free constants and the slot of each
// of its rows (has_delta: use_asc with constants in theta), and the
// dissimilarities of its rows' non-singleton nests whose log I_k is finite
// (lam_theta: each nest's position in theta, -1 for a singleton).
// nl_situation_map_clear() resets the markers.
static inline void nl_situation_map(const ChoiceLayout& lay, const choicer_off i,
                                    const int* nest0, const double* log_I_k,
                                    const arma::ivec& lam_theta, const int K,
                                    const bool has_delta,
                                    const bool include_outside_option,
                                    NlSituationMap& mp) {
  const int* alt = lay.alt + lay.row_off[i];
  const int m = lay.m(i);
  mp.delta.clear();
  mp.lam.clear();
  for (int r = 0; r < m; ++r) {
    const int j = nl_free_delta(alt[r] - 1, has_delta, include_outside_option);
    if (j >= 0 && mp.slot_of[j] < 0) {
      mp.slot_of[j] = 0;  // seen
      mp.delta.push_back(j);
    }
    const int k = nest0[r];
    const int t = static_cast<int>(lam_theta[k]);
    if (t >= 0 && std::isfinite(log_I_k[k]) && mp.lam_slot_of[t - K] < 0) {
      mp.lam_slot_of[t - K] = 0;  // seen
      mp.lam.push_back(t - K);
    }
  }
  std::sort(mp.delta.begin(), mp.delta.end());
  std::sort(mp.lam.begin(), mp.lam.end());
  mp.a = static_cast<int>(mp.delta.size());
  mp.n_lam = static_cast<int>(mp.lam.size());
  for (int k = 0; k < mp.a; ++k) mp.slot_of[mp.delta[k]] = k;
  for (int k = 0; k < mp.n_lam; ++k) mp.lam_slot_of[mp.lam[k]] = k;
  mp.row_slot.resize(static_cast<std::size_t>(m));
  for (int r = 0; r < m; ++r) {
    const int j = nl_free_delta(alt[r] - 1, has_delta, include_outside_option);
    mp.row_slot[r] = j >= 0 ? mp.slot_of[j] : -1;
  }
}

static inline void nl_situation_map_clear(NlSituationMap& mp) {
  for (int k = 0; k < mp.a; ++k) mp.slot_of[mp.delta[k]] = -1;
  for (int k = 0; k < mp.n_lam; ++k) mp.lam_slot_of[mp.lam[k]] = -1;
}

// The weight-free score of situation i on its compact block (the map built
// for i), into s_c, and its V-space gradient g into grad_vec, with the
// expressions of the dense kernels: the coefficients' entries as X_i' g, each
// slot's differences in row order at scale 1.0 (as scatter_delta_grad()), and
// each dissimilarity's score from its nest's sum of P(j|k) V_j. lam_nest maps
// a dissimilarity to its nest; sum_P_V (one entry per nest) is zero on entry
// and left zero.
static inline void nl_situation_score(const arma::mat& X,
                                      const arma::uword start_idx, const int m_i,
                                      const int* nest0, const double* V_inside,
                                      NlProbs& pr, const arma::vec& lambda,
                                      const int* lam_nest,
                                      const int chosen_nest_k,
                                      const int chosen_inside_idx,
                                      const NlSituationMap& mp, const int K,
                                      double* sum_P_V, arma::vec& grad_vec,
                                      arma::vec& s_c) {
  const auto X_i = X.rows(start_idx, start_idx + m_i - 1); // m_i x K
  const arma::vec P_i(pr.P_i.data(), m_i, false, true);
  const double* P_j_given_k = pr.P_j_given_k.data();
  const double* P_k = pr.P_k.data();
  const double* log_I_k = pr.log_I_k.data();

  // grad_vec[j] = -P_i[j]
  //             + P_j_given_k[j] * (1 - 1/lambda_k)  if j in chosen nest
  //             + 1/lambda_k                          if j == chosen alt
  grad_vec.set_size(m_i);
  grad_vec = -P_i;
  if (chosen_nest_k >= 0) { // not outside option
    const double lambda_chosen = lambda[chosen_nest_k];
    const double term_scale = 1.0 - 1.0 / lambda_chosen;
    for (int j = 0; j < m_i; ++j) {
      if (nest0[j] == chosen_nest_k) {
        grad_vec[j] += P_j_given_k[j] * term_scale;
      }
    }
    grad_vec[chosen_inside_idx] += 1.0 / lambda_chosen;
  }

  s_c.zeros();
  s_c.subvec(0, K - 1) = X_i.t() * grad_vec;
  const int d0 = K + mp.n_lam;  // the first constant slot
  for (int r = 0; r < m_i; ++r) {
    const int k = mp.row_slot[r];
    if (k >= 0) s_c[d0 + k] += 1.0 * grad_vec[r];
  }

  // The dissimilarities: sum_{j in B_k} P(j|B_k) * V_ij, then each touched
  // nest's score
  for (int j = 0; j < m_i; ++j) {
    sum_P_V[nest0[j]] += P_j_given_k[j] * V_inside[j];
  }
  for (int l = 0; l < mp.n_lam; ++l) {
    const int k = lam_nest[mp.lam[l]];
    const double lambda_k = lambda[k];
    double grad_lambda_k = 0.0;
    double term_in_brackets = log_I_k[k] - (1.0 / lambda_k) * sum_P_V[k];
    double P_k_k = P_k[k];

    if (k == chosen_nest_k) {
      grad_lambda_k += (1.0 - P_k_k) * term_in_brackets;
      double V_chosen = V_inside[chosen_inside_idx];
      grad_lambda_k += (1.0 / (lambda_k * lambda_k)) * (sum_P_V[k] - V_chosen);
    } else {
      grad_lambda_k += -P_k_k * term_in_brackets;
    }
    s_c[K + l] += grad_lambda_k;
  }
  for (int j = 0; j < m_i; ++j) sum_P_V[nest0[j]] = 0.0;
}

// The NL gradient kernel's pass over the choice situations: the weighted
// sums of their log-probabilities and scores and, with OpgDiag, the diagonal
// of the BHHH matrix, sum_i w_i s_i o s_i, its terms formed as
// nl_bhhh_parallel() forms them (the beta block from the unweighted product
// X_i' grad_vec, the constants scattered at scale 1.0 and squared once each,
// lambda's weight-free terms) and rounded before they are added
// (OpgDiagTerms). Without OpgDiag it is the kernel's loop as it was.
template <bool OpgDiag>
static void nl_gradient_pass(const arma::mat& X, const ChoiceLayout& lay,
                             const int* nest, const arma::vec& base_util,
                             const arma::vec& lambda, const int n_nests,
                             const arma::ivec& nest_k_to_theta_idx,
                             const arma::vec& weights, const int K,
                             const int n_params, const int delta_start_idx,
                             const int delta_length,
                             const bool include_outside_option,
                             const bool use_asc, double& global_loglik,
                             arma::vec& global_grad, arma::vec& global_diag,
                             const int n_threads) {
  const choicer_off N = lay.N;
  const choicer_off max_m = lay.max_m;

#ifdef _OPENMP
#pragma omp parallel num_threads(n_threads)
#endif
  {
    // Thread-local accumulators
    double local_loglik = 0.0;
    arma::vec local_grad = arma::zeros(n_params);
    arma::vec grad_vec; // pre-allocated per-thread, resized per individual
    NlProbs pr(max_m, n_nests);         // nl_individual_probs() buffers
    std::vector<int> nest_buf(max_m);   // each inside alternative's nest
    // The diagonal, the constants' scores (zero between situations), the
    // beta block's scores and a situation's terms
    arma::vec local_diag, s_d, s_b;
    OpgDiagTerms terms(OpgDiag ? K + max_m + n_nests : 0);
    if (OpgDiag) {
      local_diag.zeros(n_params);
      if (use_asc && delta_length > 0) s_d.zeros(n_params);
    }

    // Make the index map thread-private
    const arma::ivec thread_nest_k_to_theta_idx = nest_k_to_theta_idx;

#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
    for (choicer_off i = 0; i < N; ++i) {
      const int m_i         = lay.m(i);
      const arma::uword start_idx = static_cast<arma::uword>(lay.row_off[i]);
      const arma::uword end_idx   = start_idx + m_i - 1;
      const double w_i      = weights[i];

      // Get individual-specific data
      const auto X_i        = X.rows(start_idx, end_idx); // m_i x K
      const AltCodes0 alt_idx0_i{lay.alt + lay.row_off[i]}; // m_i
      int* nest_idx0_i      = nest_buf.data(); // m_i
      nl_gather_nests(nest_idx0_i, lay.alt + lay.row_off[i], nest, m_i);

      // --- 2. Calculate Utilities and Probabilities ---

      // V_ij = X_ij * beta + delta_j (pre-computed)
      const double* V_inside = base_util.memptr() + start_idx;

      // Per-individual probability block (shared helper)
      nl_individual_probs(V_inside, nest_idx0_i, m_i, lambda, n_nests,
                          include_outside_option, pr);
      const arma::vec P_i(pr.P_i.data(), m_i, false, true);
      const double* P_j_given_k = pr.P_j_given_k.data();
      const double* P_k = pr.P_k.data();
      const double* log_I_k = pr.log_I_k.data();
      const double* log_P_i = pr.log_P_i.data();
      const double log_P_outside = pr.log_P_outside;

      // --- 3. Log-Likelihood Calculation ---

      const int chosen_alt_idx = lay.choice[i];
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

      // OPG diagonal (opg_diag), beta and delta: the weight-free score's
      // squares, the beta block from the unweighted product as
      // nl_bhhh_parallel() forms it
      if (OpgDiag) {
        s_b = X_i.t() * grad_vec;
        for (int k = 0; k < K; ++k) terms.push(k, s_b[k], w_i);
        if (use_asc && delta_length > 0) {
          scatter_delta_grad(s_d, delta_start_idx, grad_vec, alt_idx0_i,
                             m_i, include_outside_option, 1.0);
          opg_diag_push_delta(terms, s_d, delta_start_idx, alt_idx0_i,
                              m_i, include_outside_option, w_i);
        }
      }

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
        if (OpgDiag) terms.push(theta_idx_k, grad_lambda_k, w_i);

      } // end gradient loop for lambda
      if (OpgDiag) terms.flush(local_diag);
    } // end of i loop

    // Combine partial accumulators into global accumulators
#ifdef _OPENMP
#pragma omp critical
#endif
    {
      global_loglik += local_loglik;
      global_grad += local_grad;
      if (OpgDiag) global_diag += local_diag;
    }
  } // end parallel region
}

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
//' N <- 50; J <- 4
//' dt <- data.table(id = rep(1:N, each = J), alt = rep(1:J, N))
//' dt[, `:=`(x1 = rnorm(.N), x2 = rnorm(.N))]
//' dt[, nest := ifelse(alt <= 2, "A", "B")]
//' dt[, choice := 0L]
//' dt[, choice := sample(c(1L, rep(0L, J - 1))), by = id]
//' d <- prepare_nl_data(dt, "id", "alt", "choice", c("x1", "x2"), "nest")
//' K_x <- ncol(d$X); K_l <- length(unique(d$nest_idx))
//' theta <- c(rep(0, K_x), rep(0.5, K_l), rep(0, J - 1))
//' result <- choicer:::nl_loglik_gradient_parallel(theta, d$X, d$alt_idx,
//'   d$choice_idx, d$nest_idx, d$M, d$weights)
//' result$objective
//' }
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List nl_loglik_gradient_parallel(
    const arma::vec& theta,
    const arma::mat& X,
    const Rcpp::IntegerVector& alt_idx,
    const Rcpp::IntegerVector& choice_idx,
    const Rcpp::IntegerVector& nest_idx,
    const Rcpp::IntegerVector& M,
    const arma::vec& weights,
    const bool use_asc = true,
    const bool include_outside_option = false,
    const bool opg_diag = false
) {
  // Extract dimensions
  const int K = X.n_cols;
  const int n_params = theta.n_elem;
  const int n_nests = nl_n_nests(nest_idx);

  // --- 1. Parameter Parsing (shared helper; also returns the map from full
  // nest index k to lambda_k's position in theta, -1 for singleton nests) ---
  arma::vec beta, lambda, delta;
  int delta_start_idx, delta_length;
  arma::ivec nest_k_to_theta_idx;
  nl_parse_theta(theta, nest_idx, K, use_asc, include_outside_option,
                 beta, lambda, delta, delta_start_idx, delta_length,
                 &nest_k_to_theta_idx);
  // Layout of the stacked design: 64-bit row offsets, and the alternative
  // codes, choices and nest codes read in place
  const ChoiceLayout lay = nl_layout_build(X, alt_idx, nest_idx, M, use_asc,
                                           delta, &weights, &choice_idx);
  const int* nest = nest_idx.begin(); // alternative a is in nest nest[a] - 1

  // --- H2: Serial pre-loop validation of chosen-alternative indices ---
  // Rcpp::stop() is only safe outside parallel regions.
  validate_choices(lay, include_outside_option, nullptr);

  // Memory: each thread's sums and situation buffers (with opg_diag also its
  // diagonal, the constants' scores and the terms), X beta, the sums and the
  // returned copies. Each thread allocates its sums and buffers before its
  // first situation, so the pass runs on the team checked here.
  const int n_threads = kernel_team_threads(lay.N);
  kernel_check_memory(
      "gradient", n_threads,
      nl_thread_bytes(lay, K, n_nests, n_params, opg_diag ? 3.0 : 1.0) +
          (opg_diag ? 16.0 * (K + static_cast<double>(lay.max_m) + n_nests) + 8.0 * K
                    : 0.0),
      [&](const int) {
        return 8.0 * static_cast<double>(lay.n_rows) +
               (opg_diag ? 32.0 : 16.0) * n_params;
      });

  // Pre-compute base utility for all individuals (single BLAS call)
  arma::vec base_util = compute_base_util(X, beta, lay, use_asc, delta);

  // Prepare global accumulators
  double global_loglik = 0.0;
  arma::vec global_grad = arma::zeros(n_params);
  arma::vec global_diag;
  if (opg_diag) {
    global_diag.zeros(n_params);
    nl_gradient_pass<true>(X, lay, nest, base_util, lambda, n_nests,
                           nest_k_to_theta_idx, weights, K, n_params,
                           delta_start_idx, delta_length,
                           include_outside_option, use_asc, global_loglik,
                           global_grad, global_diag, n_threads);
  } else {
    nl_gradient_pass<false>(X, lay, nest, base_util, lambda, n_nests,
                            nest_k_to_theta_idx, weights, K, n_params,
                            delta_start_idx, delta_length,
                            include_outside_option, use_asc, global_loglik,
                            global_grad, global_diag, n_threads);
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

//' BHHH/OPG information matrix for the Nested Logit model
//'
//' Computes the weighted outer product of per-individual scores
//' \eqn{\sum_i w_i\, s_i s_i^\top} for the Nested Logit model. The
//' per-individual score \eqn{s_i} (over the beta, lambda and delta/ASC blocks)
//' is the (positive) gradient of individual \eqn{i}'s log-likelihood
//' contribution and is weight-free; the supplied \code{weights} enter only as
//' the leading multiplier. Passing \code{weights = w} yields the ordinary
//' weighted BHHH/OPG information; passing \code{weights = w^2} yields the
//' sandwich \emph{meat} \eqn{B = \sum_i w_i^2 s_i s_i^\top} for robust (WESML)
//' inference. Singleton-nest lambdas are fixed to 1 and contribute no score
//' (mirroring the gradient kernel). Each situation's outer product is formed
//' only over the parameters its likelihood involves (the coefficients, the
//' dissimilarities of the nests it offers and the constants of its
//' alternatives), so the work per situation does not grow with the number of
//' alternatives or nests; the n_params x n_params result and the threads'
//' accumulators do (see acc_bytes).
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
//' @param acc_bytes Bytes for the threads' private upper triangles (2 GiB by
//'   default), as in \code{mnl_bhhh_parallel()}; in a shared result the
//'   coefficients' and the dissimilarities' rows stay private. \code{0}
//'   forces the shared result (for tests).
//' @returns A symmetric positive-semidefinite information matrix
//'   \eqn{\sum_i w_i\, s_i s_i^\top} (same sign convention as the negated Hessian).
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
//' B <- choicer:::nl_bhhh_parallel(theta, d$X, d$alt_idx, d$choice_idx,
//'   d$nest_idx, d$M, d$weights)
//' dim(B)
//' }
//' @keywords internal
// [[Rcpp::export]]
Rcpp::NumericMatrix nl_bhhh_parallel(
    const arma::vec& theta,
    const arma::mat& X,
    const Rcpp::IntegerVector& alt_idx,
    const Rcpp::IntegerVector& choice_idx,
    const Rcpp::IntegerVector& nest_idx,
    const Rcpp::IntegerVector& M,
    const arma::vec& weights,
    const bool use_asc = true,
    const bool include_outside_option = false,
    const double acc_bytes = 2147483648.0
) {
  const int K = X.n_cols;
  const int n_params = theta.n_elem;
  const int n_nests = nl_n_nests(nest_idx);

  // Parameter parsing (shared helper; also returns nest-k -> theta index map)
  arma::vec beta, lambda, delta;
  int delta_start_idx, delta_length;
  arma::ivec nest_k_to_theta_idx;
  nl_parse_theta(theta, nest_idx, K, use_asc, include_outside_option,
                 beta, lambda, delta, delta_start_idx, delta_length,
                 &nest_k_to_theta_idx);
  // Layout of the stacked design: 64-bit row offsets, and the alternative
  // codes, choices and nest codes read in place
  const ChoiceLayout lay = nl_layout_build(X, alt_idx, nest_idx, M, use_asc,
                                           delta, &weights, &choice_idx);
  const choicer_off N = lay.N;
  const int* nest = nest_idx.begin(); // alternative a is in nest nest[a] - 1

  // --- Serial pre-loop validation of chosen-alternative indices ---
  validate_choices(lay, include_outside_option, "nl_bhhh_parallel");

  // Memory: each thread's slot map, compact score, a column of the update
  // (every parameter, for the dense fallback) and situation buffers; X beta,
  // the result and the threads' accumulators (the coefficients' and
  // dissimilarities' rows private in a shared result)
  const int K_l = delta_start_idx - K;  // dissimilarities
  const bool has_delta = use_asc && delta_length > 0;
  const int J_d = has_delta ? delta_length : 0;  // free constants
  const double m_max = static_cast<double>(lay.max_m);
  const double blk_max = K + std::min<double>(K_l, m_max) + std::min<double>(J_d, m_max);
  const int n_threads = kernel_team_threads(N);
  kernel_check_memory("BHHH matrix", n_threads,
                      nl_map_bytes(J_d, K_l, lay.max_m) + 8.0 * (blk_max + n_params + 64.0) +
                          nl_thread_bytes(lay, K, n_nests, n_params, 0.0),
                      [&](const int t) {
                        return 8.0 * static_cast<double>(lay.n_rows) + 4.0 * K_l +
                               sym_acc_need(n_params, K + K_l, t, acc_bytes);
                      });

  // Pre-compute base utility for all individuals (single BLAS call)
  arma::vec base_util = compute_base_util(X, beta, lay, use_asc, delta);

  // Each dissimilarity's nest
  std::vector<int> lam_nest(static_cast<std::size_t>(K_l) + 1);
  for (int k = 0; k < n_nests; ++k) {
    const arma::sword t = nest_k_to_theta_idx[k];
    if (t >= 0) lam_nest[static_cast<std::size_t>(t - K)] = k;
  }

  // The result and the threads' accumulators, slot maps and buffers, all
  // allocated here on the primary thread
  SymAcc acc(n_params, n_threads, K + K_l, "BHHH", acc_bytes);
  NlSituationMaps maps = nl_situation_maps(J_d, K_l, lay.max_m, n_threads);
  std::vector<NlThreadBuffers> bufs = nl_thread_buffers(lay.max_m, n_nests, n_threads);
  std::vector<std::vector<double>> scores =
      nl_thread_vectors(static_cast<std::size_t>(blk_max), n_threads);
  std::vector<std::vector<double>> cols =
      nl_thread_vectors(static_cast<std::size_t>(n_params), n_threads);

#ifdef _OPENMP
#pragma omp parallel num_threads(n_threads)
#endif
  {
    const int tid = kernel_thread_num();
    NlSituationMap& mp = *maps[tid];
    NlThreadBuffers& b = bufs[tid];
    int* nest0 = b.nest0.data();
    double* acc_t = acc.base(tid);
    double* s_buf = scores[tid].data();
    double* col = cols[tid].data();
    arma::vec grad_vec; // resized per situation

#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
    for (choicer_off i = 0; i < N; ++i) {
      const int m_i         = lay.m(i);
      const arma::uword start_idx = static_cast<arma::uword>(lay.row_off[i]);
      const double w_i      = weights[i];
      nl_gather_nests(nest0, lay.alt + lay.row_off[i], nest, m_i);

      // V_ij = X_ij * beta + delta_j (pre-computed)
      const double* V_inside = base_util.memptr() + start_idx;

      // Per-individual probability block (shared helper)
      nl_individual_probs(V_inside, nest0, m_i, lambda, n_nests,
                          include_outside_option, b.pr);

      const int chosen_alt_idx = lay.choice[i];
      int chosen_nest_k = -1;
      int chosen_inside_idx = -1;
      if (!(include_outside_option && chosen_alt_idx == 0)) {
        chosen_inside_idx = chosen_alt_idx - 1;
        chosen_nest_k = nest0[chosen_inside_idx];
      }

      nl_situation_map(lay, i, nest0, b.pr.log_I_k.data(), nest_k_to_theta_idx,
                       K, has_delta, include_outside_option, mp);
      arma::vec s_c(s_buf, K + mp.n_lam + mp.a, false, true);
      nl_situation_score(X, start_idx, m_i, nest0, V_inside, b.pr, lambda,
                         lam_nest.data(), chosen_nest_k, chosen_inside_idx, mp,
                         K, b.sum_P_V.data(), grad_vec, s_c);

      // w_i s_c s_c' over the compact block, each entry formed as the dense
      // rank-one update formed it (mnl_bhhh_parallel()); a non-finite weight
      // or score takes the dense update over every entry.
      if (std::isfinite(w_i) && s_c.is_finite()) {
        if (acc.shared) {
          bhhh_add<true>(acc, acc_t, mp, K, s_c, w_i, col);
        } else {
          bhhh_add<false>(acc, acc_t, mp, K, s_c, w_i, col);
        }
      } else {
        if (acc.shared) {
          bhhh_add_dense<true>(acc, acc_t, mp, K, s_c, w_i, col);
        } else {
          bhhh_add_dense<false>(acc, acc_t, mp, K, s_c, w_i, col);
        }
      }
      nl_situation_map_clear(mp);
    } // end of i loop
  } // end parallel region

  // Return PSD information matrix (same sign convention as negated Hessian).
  acc.finish(false);
  return acc.result;
}

// Cluster-robust meat for the nested logit model (internal).
//
// Returns the n_params x n_params matrix sum_g s_g s_g' with
// s_g = sum_{i in g} w_i s_i the weighted score of cluster g, s_i the
// weight-free score of choice situation i (nl_scores_parallel()'s row) and
// w_i its weight: what crossprod(rowsum(w * S, cluster)) forms from the
// N x n_params score matrix, without that matrix. `cluster` holds one code
// per choice situation, 1..G (R passes match(as.character(labels),
// unique(...))); codes need not all occur, but the sort's arrays have G + 1
// entries, so G should not be far above the number of situations.
//
// As mnl_cluster_meat_parallel(), with the nested logit's compact score:
// the situations sorted by cluster (a stable counting sort), cut into chunks
// of similar numbers of rows processed in parallel, each cluster's weighted
// scores summed over the parameters its situations touch, clusters cut by a
// chunk boundary merged after the parallel region. A cluster whose weight,
// score or sum is not finite is updated over every parameter, as
// crossprod() spread it. acc_bytes chooses between the threads' private
// triangles and a shared result, as for nl_bhhh_parallel().
// [[Rcpp::export]]
Rcpp::NumericMatrix nl_cluster_meat_parallel(
    const arma::vec& theta,
    const arma::mat& X,
    const Rcpp::IntegerVector& alt_idx,
    const Rcpp::IntegerVector& choice_idx,
    const Rcpp::IntegerVector& nest_idx,
    const Rcpp::IntegerVector& M,
    const arma::vec& weights,
    const Rcpp::IntegerVector& cluster,
    const bool use_asc = true,
    const bool include_outside_option = false,
    const double acc_bytes = 2147483648.0
) {
  const int K = X.n_cols;
  const int n_params = theta.n_elem;
  const int n_nests = nl_n_nests(nest_idx);

  // Inputs, layout and choices validated on the primary thread as in the
  // other weighted kernels; then the cluster codes.
  arma::vec beta, lambda, delta;
  int delta_start_idx, delta_length;
  arma::ivec nest_k_to_theta_idx;
  nl_parse_theta(theta, nest_idx, K, use_asc, include_outside_option,
                 beta, lambda, delta, delta_start_idx, delta_length,
                 &nest_k_to_theta_idx);
  const ChoiceLayout lay = nl_layout_build(X, alt_idx, nest_idx, M, use_asc,
                                           delta, &weights, &choice_idx);
  const choicer_off N = lay.N;
  const int* nest = nest_idx.begin(); // alternative a is in nest nest[a] - 1
  validate_choices(lay, include_outside_option, "nl_cluster_meat_parallel");
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
  const int K_l = delta_start_idx - K;  // dissimilarities
  const bool has_delta = use_asc && delta_length > 0;
  const int J_d = has_delta ? delta_length : 0;  // free constants
  const double m_max = static_cast<double>(lay.max_m);
  const double blk_max = K + std::min<double>(K_l, m_max) + std::min<double>(J_d, m_max);
  const int n_threads = kernel_team_threads(N);
  kernel_check_memory("cluster meat", n_threads,
                      nl_map_bytes(J_d, K_l, lay.max_m) + 8.0 * (blk_max + 64.0) +
                          nl_thread_bytes(lay, K, n_nests, n_params, 0.0) +
                          9.0 * ClusterSum::bytes(n_params),
                      [&](const int t) {
                        return 8.0 * static_cast<double>(lay.n_rows) + 4.0 * K_l +
                               sym_acc_need(n_params, K + K_l, t, acc_bytes) +
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
  arma::vec base_util = compute_base_util(X, beta, lay, use_asc, delta);

  // Each dissimilarity's nest
  std::vector<int> lam_nest(static_cast<std::size_t>(K_l) + 1);
  for (int k = 0; k < n_nests; ++k) {
    const arma::sword t = nest_k_to_theta_idx[k];
    if (t >= 0) lam_nest[static_cast<std::size_t>(t - K)] = k;
  }

  // The result and the threads' accumulators, slot maps, buffers and cluster
  // sums, and the chunks' partial sums of their cut clusters (two per chunk:
  // the cluster it starts in and the one it ends in), all allocated here on
  // the primary thread.
  SymAcc acc(n_params, n_threads, K + K_l, "cluster meat", acc_bytes);
  NlSituationMaps maps = nl_situation_maps(J_d, K_l, lay.max_m, n_threads);
  std::vector<NlThreadBuffers> bufs = nl_thread_buffers(lay.max_m, n_nests, n_threads);
  std::vector<std::vector<double>> scores =
      nl_thread_vectors(static_cast<std::size_t>(blk_max), n_threads);
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
    NlSituationMap& mp = *maps[tid];
    NlThreadBuffers& b = bufs[tid];
    int* nest0 = b.nest0.data();
    ClusterSum& cs = sums[tid];
    double* acc_t = acc.base(tid);
    double* s_buf = scores[tid].data();
    arma::vec grad_vec; // resized per situation

#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
    for (choicer_off c = 0; c < n_chunks; ++c) {
      const choicer_off q0 = chunk[c], q1 = chunk[c + 1];
      for (choicer_off q = q0; q < q1; ++q) {
        const choicer_off i = order[q];
        const int g = cl[i];
        if (q == q0 || cl[order[q - 1]] != g) cs.start(g, K);

        const int m_i = lay.m(i);
        const arma::uword start_idx = static_cast<arma::uword>(lay.row_off[i]);
        nl_gather_nests(nest0, lay.alt + lay.row_off[i], nest, m_i);
        const double* V_inside = base_util.memptr() + start_idx;
        nl_individual_probs(V_inside, nest0, m_i, lambda, n_nests,
                            include_outside_option, b.pr);
        const int chosen_alt_idx = lay.choice[i];
        int chosen_nest_k = -1;
        int chosen_inside_idx = -1;
        if (!(include_outside_option && chosen_alt_idx == 0)) {
          chosen_inside_idx = chosen_alt_idx - 1;
          chosen_nest_k = nest0[chosen_inside_idx];
        }
        nl_situation_map(lay, i, nest0, b.pr.log_I_k.data(), nest_k_to_theta_idx,
                         K, has_delta, include_outside_option, mp);
        arma::vec s_c(s_buf, K + mp.n_lam + mp.a, false, true);
        nl_situation_score(X, start_idx, m_i, nest0, V_inside, b.pr, lambda,
                           lam_nest.data(), chosen_nest_k, chosen_inside_idx, mp,
                           K, b.sum_P_V.data(), grad_vec, s_c);
        cs.add(mp, K, s_c, weights[i]);
        nl_situation_map_clear(mp);

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

// Per-situation score matrix for the nested logit model (internal).
//
// Returns the N x n_params matrix whose row i is the weight-free score
// s_i = d log P_i / d theta over the beta, lambda and delta/ASC blocks.
// Same loop body as nl_bhhh_parallel with the outer-product accumulator
// replaced by a row write; weights are applied on the R side (see
// .assemble_score_vcov in R/classes.R). Singleton-nest lambdas are fixed to
// 1 and contribute no score (mirroring the gradient kernel).
// [[Rcpp::export]]
arma::mat nl_scores_parallel(
    const arma::vec& theta,
    const arma::mat& X,
    const Rcpp::IntegerVector& alt_idx,
    const Rcpp::IntegerVector& choice_idx,
    const Rcpp::IntegerVector& nest_idx,
    const Rcpp::IntegerVector& M,
    const bool use_asc = true,
    const bool include_outside_option = false
) {
  const int K = X.n_cols;
  const int n_params = theta.n_elem;
  const int n_nests = nl_n_nests(nest_idx);

  // Parameter parsing (shared helper; also returns nest-k -> theta index map)
  arma::vec beta, lambda, delta;
  int delta_start_idx, delta_length;
  arma::ivec nest_k_to_theta_idx;
  nl_parse_theta(theta, nest_idx, K, use_asc, include_outside_option,
                 beta, lambda, delta, delta_start_idx, delta_length,
                 &nest_k_to_theta_idx);
  // Layout of the stacked design: 64-bit row offsets, and the alternative
  // codes, choices and nest codes read in place
  const ChoiceLayout lay = nl_layout_build(X, alt_idx, nest_idx, M, use_asc,
                                           delta, nullptr, &choice_idx);
  const choicer_off N = lay.N;
  const int* nest = nest_idx.begin(); // alternative a is in nest nest[a] - 1

  // --- Serial pre-loop validation of chosen-alternative indices ---
  validate_choices(lay, include_outside_option, "nl_scores_parallel");

  // Memory: each thread's score and situation buffers; X beta, the N x n
  // matrix and its copy into R
  const int n_threads = kernel_team_threads(N);
  kernel_check_memory("score matrix", n_threads,
                      nl_thread_bytes(lay, K, n_nests, n_params, 1.0),
                      [&](const int) {
                        return 8.0 * static_cast<double>(lay.n_rows) +
                               16.0 * static_cast<double>(N) * n_params;
                      });

  // Pre-compute base utility for all individuals (single BLAS call)
  arma::vec base_util = compute_base_util(X, beta, lay, use_asc, delta);

  // Output: one row per choice situation (each written by exactly one
  // iteration, so no accumulator or critical section is needed).
  arma::mat scores(static_cast<arma::uword>(N), n_params);
  const choicer_off max_m = lay.max_m;

  // Each thread allocates its probability buffers before its first
  // situation: the team checked above
#ifdef _OPENMP
#pragma omp parallel num_threads(n_threads)
#endif
  {
    arma::vec grad_vec; // resized per individual
    arma::vec s_i;      // per-individual score
    NlProbs pr(max_m, n_nests);         // nl_individual_probs() buffers
    std::vector<int> nest_buf(max_m);   // each inside alternative's nest

    // Make the index map thread-private
    const arma::ivec thread_nest_k_to_theta_idx = nest_k_to_theta_idx;

#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
    for (choicer_off i = 0; i < N; ++i) {
      const int m_i         = lay.m(i);
      const arma::uword start_idx = static_cast<arma::uword>(lay.row_off[i]);
      const arma::uword end_idx   = start_idx + m_i - 1;

      const auto X_i        = X.rows(start_idx, end_idx); // m_i x K
      const AltCodes0 alt_idx0_i{lay.alt + lay.row_off[i]}; // m_i
      int* nest_idx0_i      = nest_buf.data(); // m_i
      nl_gather_nests(nest_idx0_i, lay.alt + lay.row_off[i], nest, m_i);

      // V_ij = X_ij * beta + delta_j (pre-computed)
      const double* V_inside = base_util.memptr() + start_idx;

      // Per-individual probability block (shared helper)
      nl_individual_probs(V_inside, nest_idx0_i, m_i, lambda, n_nests,
                          include_outside_option, pr);
      const arma::vec P_i(pr.P_i.data(), m_i, false, true);
      const double* P_j_given_k = pr.P_j_given_k.data();
      const double* P_k = pr.P_k.data();
      const double* log_I_k = pr.log_I_k.data();

      const int chosen_alt_idx = lay.choice[i];
      int chosen_nest_k = -1;
      int chosen_inside_idx = -1;
      if (!(include_outside_option && chosen_alt_idx == 0)) {
        chosen_inside_idx = chosen_alt_idx - 1;
        chosen_nest_k = nest_idx0_i[chosen_inside_idx];
      }

      // Assemble the (weight-free) per-individual score s_i.
      s_i.zeros(n_params);

      // --- beta / delta blocks: build grad_vec exactly as in the gradient ---
      // grad_vec[j] = -P_i[j]
      //             + P_j_given_k[j] * (1 - 1/lambda_k)  if j in chosen nest
      //             + 1/lambda_k                          if j == chosen alt
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

      // Beta block (scale = 1)
      s_i.subvec(0, K - 1) = X_i.t() * grad_vec;

      // Delta block: scatter
      if (use_asc && delta_length > 0) {
        scatter_delta_grad(s_i, delta_start_idx, grad_vec, alt_idx0_i,
                           m_i, include_outside_option, 1.0);
      }

      // --- lambda blocks: same expression as the gradient kernel ---
      // sum_{j in B_k} P(j|B_k) * V_ij
      arma::vec sum_P_V_k = arma::zeros(n_nests);
      for (int j = 0; j < m_i; ++j) {
        sum_P_V_k[nest_idx0_i[j]] += P_j_given_k[j] * V_inside[j];
      }
      for (int k = 0; k < n_nests; ++k) {
        const int theta_idx_k = thread_nest_k_to_theta_idx[k];
        if (theta_idx_k == -1) continue;
        if (!std::isfinite(log_I_k[k])) continue;

        const double lambda_k = lambda[k];
        double grad_lambda_k = 0.0;
        double term_in_brackets = log_I_k[k] - (1.0 / lambda_k) * sum_P_V_k[k];
        double P_k_k = P_k[k];

        if (k == chosen_nest_k) {
          grad_lambda_k += (1.0 - P_k_k) * term_in_brackets;
          double V_chosen = V_inside[chosen_inside_idx];
          grad_lambda_k += (1.0 / (lambda_k * lambda_k)) * (sum_P_V_k[k] - V_chosen);
        } else {
          grad_lambda_k += -P_k_k * term_in_brackets;
        }
        s_i[theta_idx_k] += grad_lambda_k;
      }

      scores.row(i) = s_i.t();
    } // end of i loop
  } // end parallel region

  return scores;
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
//' H <- choicer:::nl_loglik_numeric_hessian(theta, d$X, d$alt_idx, d$choice_idx,
//'   d$nest_idx, d$M, d$weights)
//' dim(H)
//' }
//' @keywords internal
// [[Rcpp::export]]
arma::mat nl_loglik_numeric_hessian(
  const arma::vec& theta,
  const arma::mat& X,
  const Rcpp::IntegerVector& alt_idx,
  const Rcpp::IntegerVector& choice_idx,
  const Rcpp::IntegerVector& nest_idx,
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

//' Analytical Hessian of the negated log-likelihood for the Nested Logit model
//'
//' Computes the exact (analytical) Hessian of the negated log-likelihood for the
//' Nested Logit model using OpenMP parallelisation with thread-local accumulators.
//' Covers all parameter blocks: beta-beta, beta-lambda, beta-delta,
//' lambda-lambda, lambda-delta, and delta-delta. Singleton nests (lambda fixed
//' to 1, not estimated) contribute no rows or columns to the lambda blocks.
//' Each situation's terms are formed only over the parameters its likelihood
//' involves (the coefficients, the dissimilarities of the nests it offers and
//' the constants of its alternatives), so the work per situation does not
//' grow with the number of alternatives or nests; the n_params x n_params
//' result and the threads' accumulators do (see acc_bytes).
//'
//' @param theta (K + n_non_singleton_nests + n_delta) parameter vector.
//'   Order: \code{[beta (K), lambda (n_non_singleton_nests), delta (n_delta)]}.
//'   Same layout as \code{nl_loglik_gradient_parallel}.
//' @param X sum(M) x K design matrix of covariates.
//' @param alt_idx sum(M)-length integer vector of 1-based alternative indices.
//' @param choice_idx N-length integer vector of 1-based chosen alternative
//'   indices; 0 indicates the outside option was chosen.
//' @param nest_idx J-length integer vector of 1-based nest indices for each
//'   inside alternative.
//' @param M N-length integer vector of alternative-set sizes.
//' @param weights N-length numeric vector of individual weights.
//' @param use_asc Logical; whether alternative-specific constants are included.
//' @param include_outside_option Logical; whether an outside option (V=0) is
//'   present.
//' @param acc_bytes Bytes for the threads' private upper triangles (2 GiB by
//'   default), as in \code{nl_bhhh_parallel()}; \code{0} forces the shared
//'   result (for tests).
//' @returns A symmetric (P x P) matrix: the Hessian of the negated
//'   log-likelihood evaluated at \code{theta}. Structurally identical to the
//'   output of \code{nl_loglik_numeric_hessian}; suitable for
//'   \code{invert_hessian()}.
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
//' K_x <- ncol(d$X)
//' K_l <- sum(table(d$nest_idx) > 1)   # number of non-singleton nests (= 2)
//' theta <- c(rep(0, K_x), rep(0.8, K_l), rep(0, J - 1))
//' H <- choicer:::nl_loglik_hessian_parallel(theta, d$X, d$alt_idx, d$choice_idx,
//'   d$nest_idx, d$M, d$weights)
//' dim(H)
//' }
//' @keywords internal
// [[Rcpp::export]]
Rcpp::NumericMatrix nl_loglik_hessian_parallel(
    const arma::vec& theta,
    const arma::mat& X,
    const Rcpp::IntegerVector& alt_idx,
    const Rcpp::IntegerVector& choice_idx,
    const Rcpp::IntegerVector& nest_idx,
    const Rcpp::IntegerVector& M,
    const arma::vec& weights,
    const bool use_asc = true,
    const bool include_outside_option = false,
    const double acc_bytes = 2147483648.0
) {
  // Extract dimensions
  const int K = X.n_cols;
  const int n_params = theta.n_elem;
  const int n_nests = nl_n_nests(nest_idx);

  // --- 1. Parameter Parsing ---
  arma::vec beta, lambda, delta;
  int delta_start_idx, delta_length;
  arma::ivec nest_k_to_theta_idx;
  nl_parse_theta(theta, nest_idx, K, use_asc, include_outside_option,
                 beta, lambda, delta, delta_start_idx, delta_length,
                 &nest_k_to_theta_idx);
  // Layout of the stacked design: 64-bit row offsets, and the alternative
  // codes, choices and nest codes read in place
  const ChoiceLayout lay = nl_layout_build(X, alt_idx, nest_idx, M, use_asc,
                                           delta, &weights, &choice_idx);
  const choicer_off N = lay.N;
  const int* nest = nest_idx.begin(); // alternative a is in nest nest[a] - 1

  // --- H2: Serial pre-loop validation of chosen-alternative indices ---
  validate_choices(lay, include_outside_option, nullptr);

  // Memory: each thread's slot map, X_i, the sums of the nests a situation
  // offers (four numbers, two K-vectors and two K x K matrices each), the
  // nests' slots, a column of Q_i, its weighted copy, d_V and situation
  // buffers; X beta, the result and the threads' accumulators (the
  // coefficients' and dissimilarities' rows private in a shared result)
  const int K_l = delta_start_idx - K;  // dissimilarities
  const bool has_delta = use_asc && delta_length > 0;
  const int J_d = has_delta ? delta_length : 0;  // free constants
  const double m_max = static_cast<double>(lay.max_m);
  const double off_max = std::min<double>(n_nests, m_max);  // nests offered
  const std::size_t stride = static_cast<std::size_t>(4 + 2 * K + 2 * K * K);
  const int n_threads = kernel_team_threads(N);
  kernel_check_memory(
      "Hessian", n_threads,
      nl_map_bytes(J_d, K_l, lay.max_m) +
          8.0 * (m_max * K + off_max * stride + 3.0 * m_max + 3.0 * K + 96.0) +
          sizeof(int) * (n_nests + off_max + 64.0) +
          nl_thread_bytes(lay, K, n_nests, n_params, 0.0),
      [&](const int t) {
        return 8.0 * static_cast<double>(lay.n_rows) + 4.0 * K_l +
               sym_acc_need(n_params, K + K_l, t, acc_bytes);
      });

  // Pre-compute base utility for all individuals (single BLAS call)
  arma::vec base_util = compute_base_util(X, beta, lay, use_asc, delta);

  // Each dissimilarity's nest
  std::vector<int> lam_nest(static_cast<std::size_t>(K_l) + 1);
  for (int k = 0; k < n_nests; ++k) {
    const arma::sword t = nest_k_to_theta_idx[k];
    if (t >= 0) lam_nest[static_cast<std::size_t>(t - K)] = k;
  }

  // The result and the threads' accumulators, slot maps and buffers, all
  // allocated here on the primary thread
  SymAcc acc(n_params, n_threads, K + K_l, "Hessian", acc_bytes);
  NlSituationMaps maps = nl_situation_maps(J_d, K_l, lay.max_m, n_threads);
  std::vector<NlThreadBuffers> bufs = nl_thread_buffers(lay.max_m, n_nests, n_threads);
  // per thread: X_i (max_m x K), the offered nests' sums (off_max x stride),
  // a column of Q_i, its weighted copy and d_V (max_m each), Xbar_i, X_i' v
  // and a coefficients' column (K each)
  const std::size_t mK = static_cast<std::size_t>(lay.max_m) * static_cast<std::size_t>(K);
  const std::size_t n_off = static_cast<std::size_t>(off_max);
  const std::size_t n_buf = mK + n_off * stride + 3 * static_cast<std::size_t>(lay.max_m) +
                            3 * static_cast<std::size_t>(K);
  std::vector<std::vector<double>> work = nl_thread_vectors(n_buf, n_threads);
  std::vector<std::vector<int>> nest_slots(
      static_cast<std::size_t>(n_threads),
      std::vector<int>(static_cast<std::size_t>(n_nests) + 32, -1));
  std::vector<std::vector<int>> offers(
      static_cast<std::size_t>(n_threads),
      std::vector<int>(n_off + 32));

#ifdef _OPENMP
#pragma omp parallel num_threads(n_threads)
#endif
  {
    const int tid = kernel_thread_num();
    NlSituationMap& mp = *maps[tid];
    NlThreadBuffers& b_t = bufs[tid];
    int* nest0 = b_t.nest0.data();
    int* nest_slot = nest_slots[tid].data();  // nest -> offered slot, or -1
    int* nests_off = offers[tid].data();      // offered nests, increasing
    double* acc_t = acc.base(tid);
    double* xbuf = work[tid].data();
    double* sums = xbuf + mK;
    double* qcol = sums + n_off * stride;
    double* hcol = qcol + lay.max_m;
    double* dV = hcol + lay.max_m;
    double* xbar_m = dV + lay.max_m;
    double* gvec_m = xbar_m + K;
    const std::size_t Ks = static_cast<std::size_t>(K);

    // Add a term t(r, c) of the situation's Hessian, with t(c, r) its
    // transposed term, into the upper triangle, r and c global, as their
    // sum: the result is halved once, as the dense code halved H + H'
    auto add2 = [&](const int r, const int c, const double t_rc, const double t_cr) {
      if (acc.shared) {
        acc.add<true>(acc_t, std::min(r, c), std::max(r, c), t_rc + t_cr);
      } else {
        acc.add<false>(acc_t, std::min(r, c), std::max(r, c), t_rc + t_cr);
      }
    };

#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
    for (choicer_off i = 0; i < N; ++i) {
      const int m_i       = lay.m(i);
      const arma::uword start_idx = static_cast<arma::uword>(lay.row_off[i]);
      const arma::uword end_idx   = start_idx + m_i - 1;
      const double w_i    = weights[i];
      nl_gather_nests(nest0, lay.alt + lay.row_off[i], nest, m_i);
      const double* V_inside     = base_util.memptr() + start_idx;

      // Step 2b: Compute probabilities
      NlProbs& pr = b_t.pr;
      nl_individual_probs(V_inside, nest0, m_i, lambda, n_nests,
                          include_outside_option, pr);
      const arma::vec P_i_vec(pr.P_i.data(), m_i, false, true);
      const double* P_j_given_k = pr.P_j_given_k.data();
      const double* P_k = pr.P_k.data();
      const double* log_I_k = pr.log_I_k.data();
      const double* log_P_i = pr.log_P_i.data();
      const double log_P_outside = pr.log_P_outside;

      // Resolve chosen alternative
      const int chosen_alt_idx_i  = lay.choice[i];
      int chosen_nest_k       = -1;
      int chosen_inside_idx   = -1;

      double log_P_choice;
      if (include_outside_option && chosen_alt_idx_i == 0) {
        log_P_choice = log_P_outside;
      } else {
        chosen_inside_idx = chosen_alt_idx_i - 1;
        log_P_choice      = log_P_i[chosen_inside_idx];
        chosen_nest_k     = static_cast<int>(nest0[chosen_inside_idx]);
      }

      // H1 clamp: skip individual if log_P_choice is non-finite
      if (!std::isfinite(log_P_choice)) continue;

      // X_i, copied once (the dense code's copy, the operand of its BLAS
      // products)
      arma::mat X_i(xbuf, m_i, K, false, true);
      X_i = X.rows(start_idx, end_idx);

      // The nests the situation offers, increasing, and their sums
      // (Step 2c): per offered slot o, at sums + o * stride, Vbar, S2Vbar,
      // T, VarV, m (K), c (K), S2 and Q (K x K, upper triangles). Nests it
      // does not offer added exact zeros in the dense code.
      int n_o = 0;
      for (int j = 0; j < m_i; ++j) {
        const int k = nest0[j];
        if (nest_slot[k] < 0) {
          nest_slot[k] = 0;  // seen
          nests_off[n_o++] = k;
        }
      }
      std::sort(nests_off, nests_off + n_o);
      for (int o = 0; o < n_o; ++o) nest_slot[nests_off[o]] = o;
      std::fill(sums, sums + static_cast<std::size_t>(n_o) * stride, 0.0);
      auto nest_sums = [&](const int k) { return sums + static_cast<std::size_t>(nest_slot[k]) * stride; };
      for (int j = 0; j < m_i; ++j) {
        double* sk = nest_sums(nest0[j]);
        double* m_k = sk + 4;
        double* c_k = m_k + K;
        double* S2_k = c_k + K;
        double* Q_k = S2_k + Ks * Ks;
        const double pj_given_k = P_j_given_k[j];
        const double p_ij       = P_i_vec[j];
        const double V_j        = V_inside[j];

        sk[0] += pj_given_k * V_j;          // Vbar_k
        sk[1] += pj_given_k * V_j * V_j;    // S2Vbar_k
        // m_k += p_ij x_j and c_k += P(j|k) x_j, each product rounded
        // before it is added, as the dense code's vector updates did
        for (int e = 0; e < K; ++e) {
          const double x_e = X_i(j, e);
          const double tm = x_e * p_ij;
          m_k[e] += tm;
          const double tc = x_e * pj_given_k;
          c_k[e] += tc;
        }
        // S2_k += p_ij x_j x_j' and Q_k += P(j|k) x_j x_j' (upper triangles)
        for (int c = 0; c < K; ++c) {
          for (int r = 0; r <= c; ++r) {
            const double xx = X_i(j, r) * X_i(j, c);
            const double ts = xx * p_ij;
            S2_k[r + c * Ks] += ts;
            const double tq = xx * pj_given_k;
            Q_k[r + c * Ks] += tq;
          }
        }
      }
      // Per-nest derived quantities
      for (int o = 0; o < n_o; ++o) {
        double* sk = sums + static_cast<std::size_t>(o) * stride;
        const int k = nests_off[o];
        sk[2] = log_I_k[k] - sk[0] / lambda[k];  // T_k
        sk[3] = sk[1] - sk[0] * sk[0];           // VarV_k
      }

      // P_i-weighted global mean: Xbar_i = X_i^T * P_i  (K-vector)
      arma::vec Xbar_i(xbar_m, K, false, true);
      Xbar_i = X_i.t() * P_i_vec;

      // The situation's slots: its free constants and touched dissimilarities
      nl_situation_map(lay, i, nest0, log_I_k, nest_k_to_theta_idx, K,
                       has_delta, include_outside_option, mp);
      const int d0 = K + K_l;  // global index of the first constant

      // Q_i[a, b], the V-space Hessian of the negated log-likelihood
      // (the dense code's expression)
      auto Q_ab = [&](const int a, const int b) {
        const int ka = static_cast<int>(nest0[a]);
        const double Pia = P_i_vec[a];
        const double Pa_given_ka = P_j_given_k[a];
        const double lam_ka = lambda[ka];
        const double inv_lam_ka = 1.0 / lam_ka;
        const int kb = static_cast<int>(nest0[b]);
        const double Pib = P_i_vec[b];
        const double Pb_given_kb = P_j_given_k[b];

        // Base term: -P_ia * P_ib (all pairs)
        double Qab = -Pia * Pib;

        // Within-nest correction: P_ia * [1_{a=b}/lam_k + P(b|k) * (1 - 1/lam_k)]
        if (ka == kb) {
          const double indicator_ab = (a == b) ? 1.0 : 0.0;
          Qab += Pia * (indicator_ab * inv_lam_ka + Pb_given_kb * (1.0 - inv_lam_ka));
        }

        // Chosen-nest correction C(a,b), nonzero when k(a)=k(b)=ki
        if (chosen_nest_k >= 0 && ka == chosen_nest_k && kb == chosen_nest_k) {
          const double lam_ki = lambda[chosen_nest_k];
          const double scale_C = (1.0 - 1.0 / lam_ki) / lam_ki;
          const double ind_ab = (a == b) ? 1.0 : 0.0;
          Qab -= scale_C * Pa_given_ka * (ind_ab - Pb_given_kb);
        }
        return Qab;
      };

      // Steps 2d and 2e, delta-delta and beta-delta: for each row b with a
      // free constant, the column Q_i[, b]; w_i Q_i[a, b] into the constants'
      // block for each row a with one, and X_i' (w_i Q_i[, b]) into the
      // coefficients' column of b's constant
      if (has_delta) {
        arma::vec h_col_b(hcol, m_i, false, true);
        arma::vec beta_delta_contrib(gvec_m, K, false, true);
        for (int b = 0; b < m_i; ++b) {
          const int sb = mp.row_slot[b];
          if (sb < 0) continue;
          const int gb = d0 + mp.delta[sb];
          for (int a = 0; a < m_i; ++a) qcol[a] = Q_ab(a, b);
          for (int a = 0; a < m_i; ++a) {
            const int sa = mp.row_slot[a];
            if (sa < 0) continue;
            const double t = w_i * qcol[a];
            if (sa == sb) {
              add2(gb, gb, t, t);
            } else {
              add2(d0 + mp.delta[sa], gb, t, 0.0);
            }
          }
          for (int a = 0; a < m_i; ++a) h_col_b[a] = qcol[a] * w_i;
          beta_delta_contrib = X_i.t() * h_col_b;  // K-vector
          for (int r = 0; r < K; ++r) {
            add2(r, gb, beta_delta_contrib[r], beta_delta_contrib[r]);
          }
        }
      }

      // Step 2e: Beta-beta block = w_i * X_i^T * Q_i * X_i, term by term in
      // the dense code's order, each entry as it formed it:
      // Term A (subtract): -w_i * Xbar_i * Xbar_i^T
      for (int c = 0; c < K; ++c) {
        for (int r = 0; r <= c; ++r) {
          const double xx = Xbar_i[r] * Xbar_i[c];
          const double t = xx * w_i;
          add2(r, c, -t, -t);
        }
      }
      // Term B (add): w_i * ((1 - 1/lam_k) m_k c_k^T + (1/lam_k) S2_k) for
      // each nest offered, in nest order
      for (int o = 0; o < n_o; ++o) {
        const double* sk = sums + static_cast<std::size_t>(o) * stride;
        const double* m_k = sk + 4;
        const double* c_k = m_k + K;
        const double* S2_k = c_k + K;
        const double lam_k   = lambda[nests_off[o]];
        const double inv_lam = 1.0 / lam_k;
        const double scale_B = 1.0 - inv_lam;  // = 0 for singletons (lam=1)
        for (int c = 0; c < K; ++c) {
          for (int r = 0; r <= c; ++r) {
            const double s2 = S2_k[r + c * Ks];
            const double o_rc = m_k[r] * c_k[c];
            const double b_rc = o_rc * scale_B;
            const double s_rc = s2 * inv_lam;
            double t_rc = b_rc + s_rc;
            t_rc *= w_i;
            const double o_cr = m_k[c] * c_k[r];
            const double b_cr = o_cr * scale_B;
            double t_cr = b_cr + s_rc;
            t_cr *= w_i;
            add2(r, c, t_rc, t_cr);
          }
        }
      }
      // Term C (subtract): -w_i * scale_Ci * (Q_ki - c_ki c_ki^T)
      if (chosen_nest_k >= 0) {
        const double* sk = nest_sums(chosen_nest_k);
        const double* c_k = sk + 4 + K;
        const double* Q_k = c_k + K + Ks * Ks;
        const double lam_ki = lambda[chosen_nest_k];
        const double scale_C = (1.0 - 1.0 / lam_ki) / lam_ki;
        const double w_scale_C = w_i * scale_C;
        for (int c = 0; c < K; ++c) {
          for (int r = 0; r <= c; ++r) {
            const double cc = c_k[r] * c_k[c];
            double t = Q_k[r + c * Ks] - cc;
            t *= w_scale_C;
            add2(r, c, -t, -t);
          }
        }
      }

      // Step 2f: Lambda-beta and lambda-delta cross blocks, and step 2g, the
      // lambda-lambda block, for each dissimilarity the situation touches
      arma::vec d_V_k(dV, m_i, false, true);
      arma::vec lambda_beta_contrib(gvec_m, K, false, true);
      for (int l = 0; l < mp.n_lam; ++l) {
        const int k = lam_nest[mp.lam[l]];
        const int gk = K + mp.lam[l];
        const double* sk = nest_sums(k);
        const double lam_k = lambda[k];
        const double P_k_k = P_k[k];
        const double T_kk = sk[2];
        const double Vbar_kk = sk[0];

        // d_V^(k)[b] = dg_{lambda_k}/dV_b
        if (k != chosen_nest_k) {
          // Non-chosen nest formula
          for (int b = 0; b < m_i; ++b) {
            const int kb = static_cast<int>(nest0[b]);
            d_V_k[b] = P_k_k * T_kk * P_i_vec[b];
            if (kb == k) {
              d_V_k[b] += P_k_k * P_j_given_k[b]
                          * ((V_inside[b] - Vbar_kk) / (lam_k * lam_k) - T_kk);
            }
          }
        } else {
          // Chosen nest formula (k == chosen_nest_k, ki >= 0)
          for (int b = 0; b < m_i; ++b) {
            const int kb = static_cast<int>(nest0[b]);
            d_V_k[b] = P_k_k * T_kk * P_i_vec[b];
            if (b == chosen_inside_idx) {
              d_V_k[b] -= 1.0 / (lam_k * lam_k);
            }
            if (kb == k) {
              double bracket = -P_k_k * T_kk
                + (V_inside[b] - Vbar_kk) / (lam_k * lam_k)
                  * (P_k_k - 1.0 + 1.0 / lam_k)
                + 1.0 / (lam_k * lam_k);
              d_V_k[b] += P_j_given_k[b] * bracket;
            }
          }
        }

        // lambda-beta cross: K-vector contribution
        lambda_beta_contrib = -w_i * (X_i.t() * d_V_k);
        for (int r = 0; r < K; ++r) {
          add2(r, gk, lambda_beta_contrib[r], lambda_beta_contrib[r]);
        }

        // lambda-delta cross: scatter over free-ASC alternatives
        if (has_delta) {
          for (int b = 0; b < m_i; ++b) {
            const int sb = mp.row_slot[b];
            if (sb < 0) continue;
            const double val = -w_i * d_V_k[b];
            add2(gk, d0 + mp.delta[sb], val, val);
          }
        }

        // lambda-lambda: off-diagonal with the later touched nests
        for (int l2 = l + 1; l2 < mp.n_lam; ++l2) {
          const int k2 = lam_nest[mp.lam[l2]];
          const double off_diag =
            -w_i * P_k[k] * P_k[k2] * T_kk * nest_sums(k2)[2];
          add2(gk, K + mp.lam[l2], off_diag, off_diag);
        }

        // Diagonal
        const double VarV_kk = sk[3];
        double diag_val;
        if (k == chosen_nest_k) {
          // Chosen-nest diagonal formula
          const double V_ca = V_inside[chosen_inside_idx];
          diag_val = w_i * (
              P_k[k] * (1.0 - P_k[k]) * T_kk * T_kk
            - (1.0 - P_k[k]) * VarV_kk / (lam_k * lam_k * lam_k)
            + VarV_kk / (lam_k * lam_k * lam_k * lam_k)
            + 2.0 * (Vbar_kk - V_ca) / (lam_k * lam_k * lam_k)
          );
        } else {
          // Non-chosen nest diagonal formula
          diag_val = w_i * (
              P_k[k] * (1.0 - P_k[k]) * T_kk * T_kk
            + P_k[k] * VarV_kk / (lam_k * lam_k * lam_k)
          );
        }
        add2(gk, gk, diag_val, diag_val);
      }

      for (int o = 0; o < n_o; ++o) nest_slot[nests_off[o]] = -1;
      nl_situation_map_clear(mp);
    }  // end individual loop
  }  // end parallel region

  // Step 3: Symmetrize and sanitize: the sums of the two triangles' terms,
  // halved, every non-finite entry set to 0
  acc.finish(false);
  double* out = acc.result.begin();
  const std::size_t nn = static_cast<std::size_t>(n_params) * static_cast<std::size_t>(n_params);
  for (std::size_t e = 0; e < nn; ++e) {
    const double h = out[e] * 0.5;
    out[e] = std::isfinite(h) ? h : 0.0;
  }
  return acc.result;
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
    const Rcpp::IntegerVector& nest_idx,
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
  if (nest_idx.size() == 0) {
    Rcpp::stop("nest_idx must be non-empty.");
  }
  const int n_nests = nl_n_nests(nest_idx); // 1-based nest_idx

  // Every code is at most n_nests, so a code is valid when it is at least 1
  // (NA_INTEGER, the most negative int, is not). Checked before n_nests sizes
  // anything.
  for (const int k : nest_idx) {
    if (k < 1) Rcpp::stop("Invalid nest index found in nest_idx.");
  }

  beta = theta.subvec(0, K - 1);

  // Identify singleton nests
  arma::uvec nest_counts = arma::zeros<arma::uvec>(n_nests);
  for (const int k : nest_idx) nest_counts[k - 1]++;
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
//' pred <- choicer:::nl_predict(coef(fit), fit$data$X, fit$data$alt_idx, fit$data$M,
//'   fit$data$nest_idx, use_asc = TRUE)
//' head(pred$choice_prob)
//' }
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List nl_predict(
    const arma::vec& theta,
    const arma::mat& X,
    const Rcpp::IntegerVector& alt_idx,
    const Rcpp::IntegerVector& M,
    const Rcpp::IntegerVector& nest_idx,
    const bool use_asc = true,
    const bool include_outside_option = false
) {
  const int K = X.n_cols;
  const int n_nests = nl_n_nests(nest_idx);

  arma::vec beta, lambda, delta;
  int delta_start_idx, delta_length;
  nl_parse_theta(theta, nest_idx, K, use_asc, include_outside_option,
                 beta, lambda, delta, delta_start_idx, delta_length);
  const ChoiceLayout lay = nl_layout_build(X, alt_idx, nest_idx, M, use_asc,
                                           delta);
  const choicer_off N = lay.N;
  const int* nest = nest_idx.begin(); // alternative a is in nest nest[a] - 1

  arma::vec base_util = compute_base_util(X, beta, lay, use_asc, delta);

  arma::vec V_all = arma::zeros(X.n_rows);
  arma::vec P_all = arma::zeros(X.n_rows);
  const choicer_off max_m = lay.max_m;

#ifdef _OPENMP
#pragma omp parallel
#endif
  {
    NlProbs pr(max_m, n_nests);         // nl_individual_probs() buffers
    std::vector<int> nest_buf(max_m);   // each inside alternative's nest

#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
    for (choicer_off i = 0; i < N; ++i) {
      const int m_i       = lay.m(i);
      const arma::uword start_idx = static_cast<arma::uword>(lay.row_off[i]);
      const arma::uword end_idx   = start_idx + m_i - 1;

      int* nest_idx0_i = nest_buf.data();
      nl_gather_nests(nest_idx0_i, lay.alt + lay.row_off[i], nest, m_i);
      const double* V_inside = base_util.memptr() + start_idx;

      nl_individual_probs(V_inside, nest_idx0_i, m_i, lambda, n_nests,
                          include_outside_option, pr);
      const arma::vec P_i(pr.P_i.data(), m_i, false, true);

      V_all.subvec(start_idx, end_idx) = base_util.subvec(start_idx, end_idx);
      P_all.subvec(start_idx, end_idx) = P_i;
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("choice_prob") = P_all,
    Rcpp::Named("utility")     = V_all
  );
}

// Predicted NL market shares (file-local internal), from the base utilities
// of all rows (X * beta plus the rows' ASCs). nest holds the 1-based nest
// codes of the alternatives, and denominator is shares_denominator(weights),
// which the callers check before they form the utilities. Returns weighted
// shares of length num_alts (index 0 = outside option when present, then
// inside alts in order).
static arma::vec nl_predict_shares_internal(
    const arma::vec& base_util,
    const ChoiceLayout& lay,
    const int* nest,
    const arma::vec& lambda,
    const arma::vec& weights,
    const double denominator,
    const int n_nests,
    const int num_alts,
    const bool include_outside_option
) {
  const choicer_off N = lay.N;
  arma::vec global_shares = arma::zeros(num_alts);
  const choicer_off max_m = lay.max_m;

#ifdef _OPENMP
#pragma omp parallel
#endif
  {
    arma::vec local_shares = arma::zeros(num_alts);
    NlProbs pr(max_m, n_nests);         // nl_individual_probs() buffers
    std::vector<int> nest_buf(max_m);   // each inside alternative's nest

#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
    for (choicer_off i = 0; i < N; ++i) {
      const int m_i       = lay.m(i);
      const arma::uword start_idx = static_cast<arma::uword>(lay.row_off[i]);
      const double w_i    = weights[i];

      const AltCodes0 alt_idx0_i{lay.alt + lay.row_off[i]};
      int* nest_idx0_i = nest_buf.data();
      nl_gather_nests(nest_idx0_i, lay.alt + lay.row_off[i], nest, m_i);
      const double* V_inside = base_util.memptr() + start_idx;

      nl_individual_probs(V_inside, nest_idx0_i, m_i, lambda, n_nests,
                          include_outside_option, pr);
      const arma::vec P_i(pr.P_i.data(), m_i, false, true);
      const double log_P_outside = pr.log_P_outside;

      if (include_outside_option) {
        local_shares(0) += w_i * std::exp(log_P_outside);
      }
      for (int a = 0; a < m_i; ++a) {
        if (include_outside_option) {
          local_shares(alt_idx0_i[a] + 1) += w_i * P_i(a);
        } else {
          local_shares(alt_idx0_i[a]) += w_i * P_i(a);
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
//' shares <- choicer:::nl_predict_shares(coef(fit), fit$data$X, fit$data$alt_idx,
//'   fit$data$M, fit$data$weights, fit$data$nest_idx, use_asc = TRUE)
//' shares
//' }
//' @keywords internal
// [[Rcpp::export]]
arma::vec nl_predict_shares(
    const arma::vec& theta,
    const arma::mat& X,
    const Rcpp::IntegerVector& alt_idx,
    const Rcpp::IntegerVector& M,
    const arma::vec& weights,
    const Rcpp::IntegerVector& nest_idx,
    const bool use_asc = true,
    const bool include_outside_option = false
) {
  const int K = X.n_cols;
  const int n_nests = nl_n_nests(nest_idx);

  arma::vec beta, lambda, delta;
  int delta_start_idx, delta_length;
  nl_parse_theta(theta, nest_idx, K, use_asc, include_outside_option,
                 beta, lambda, delta, delta_start_idx, delta_length);
  const ChoiceLayout lay = nl_layout_build(X, alt_idx, nest_idx, M, use_asc,
                                           delta, &weights);

  // total number of distinct alternatives (alt_idx.max() on the arma::uvec
  // argument this replaces threw on an empty alt_idx; kept)
  if (lay.n_rows == 0) {
    throw std::logic_error("Mat::max(): object has no elements");
  }
  const int num_alts = compute_J_total(lay.J, include_outside_option);

  const double denominator = shares_denominator(weights);
  const arma::vec base_util = compute_base_util(X, beta, lay, use_asc, delta);
  return nl_predict_shares_internal(
    base_util, lay, nest_idx.begin(), lambda, weights, denominator, n_nests,
    num_alts, include_outside_option
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
//' elas <- choicer:::nl_elasticities_parallel(coef(fit), fit$data$X, fit$data$alt_idx,
//'   fit$data$choice_idx, fit$data$nest_idx, fit$data$M, fit$data$weights,
//'   elast_var_idx = 1L)
//' elas
//' }
//' @keywords internal
// [[Rcpp::export]]
arma::mat nl_elasticities_parallel(
    const arma::vec& theta,
    const arma::mat& X,
    const Rcpp::IntegerVector& alt_idx,
    const Rcpp::IntegerVector& choice_idx, // kept for consistency, not used
    const Rcpp::IntegerVector& nest_idx,
    const Rcpp::IntegerVector& M,
    const arma::vec& weights,
    const int elast_var_idx,
    const bool use_asc = true,
    const bool include_outside_option = false
) {
  (void)choice_idx;
  const int K = X.n_cols;
  const int n_nests = nl_n_nests(nest_idx);

  const int var_idx = elast_var_idx - 1;
  if (var_idx < 0 || var_idx >= K) {
    Rcpp::stop("elast_var_idx is out of bounds for design matrix X.");
  }

  arma::vec beta, lambda, delta;
  int delta_start_idx, delta_length;
  nl_parse_theta(theta, nest_idx, K, use_asc, include_outside_option,
                 beta, lambda, delta, delta_start_idx, delta_length);
  const ChoiceLayout lay = nl_layout_build(X, alt_idx, nest_idx, M, use_asc,
                                           delta, &weights);
  const choicer_off N = lay.N;
  const int* nest = nest_idx.begin(); // alternative a is in nest nest[a] - 1
  const double beta_k = beta(var_idx);

  const int J_inside = compute_J_inside(use_asc, delta, lay);
  const int J_total = compute_J_total(J_inside, include_outside_option);

  arma::vec base_util = compute_base_util(X, beta, lay, use_asc, delta);

  arma::mat global_elas_matrix = arma::zeros(J_total, J_total);
  double global_total_weight = 0.0;
  const choicer_off max_m = lay.max_m;

#ifdef _OPENMP
#pragma omp parallel
#endif
  {
    arma::mat local_elas_matrix = arma::zeros(J_total, J_total);
    double local_total_weight = 0.0;
    NlProbs pr(max_m, n_nests);         // nl_individual_probs() buffers
    std::vector<int> nest_buf(max_m);   // each inside alternative's nest
    std::vector<int> global_map_buf(max_m); // each one's global index

#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
    for (choicer_off i = 0; i < N; ++i) {
      const int m_i       = lay.m(i);
      const arma::uword start_idx = static_cast<arma::uword>(lay.row_off[i]);
      const arma::uword end_idx   = start_idx + m_i - 1;
      const double w_i    = weights[i];
      const auto X_i      = X.rows(start_idx, end_idx);

      int* nest_idx0_i = nest_buf.data();
      nl_gather_nests(nest_idx0_i, lay.alt + lay.row_off[i], nest, m_i);
      const double* V_inside = base_util.memptr() + start_idx;

      nl_individual_probs(V_inside, nest_idx0_i, m_i, lambda, n_nests,
                          include_outside_option, pr);
      const arma::vec P_i(pr.P_i.data(), m_i, false, true);
      const double* P_j_given_k = pr.P_j_given_k.data();

      // x_k for each inside alternative
      arma::vec x_k_i = X_i.col(var_idx);

      // Global alternative index for each inside alt (outside option = 0)
      int* global_map = global_map_buf.data();
      fill_global_alt_map_inside(global_map, lay.alt + lay.row_off[i], m_i,
                                 include_outside_option);

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
//' dr <- choicer:::nl_diversion_ratios_parallel(coef(fit), fit$data$X, fit$data$alt_idx,
//'   fit$data$nest_idx, fit$data$M, fit$data$weights)
//' dr
//' }
//' @keywords internal
// [[Rcpp::export]]
arma::mat nl_diversion_ratios_parallel(
    const arma::vec& theta,
    const arma::mat& X,
    const Rcpp::IntegerVector& alt_idx,
    const Rcpp::IntegerVector& nest_idx,
    const Rcpp::IntegerVector& M,
    const arma::vec& weights,
    const bool use_asc = true,
    const bool include_outside_option = false
) {
  const int K = X.n_cols;
  const int n_nests = nl_n_nests(nest_idx);

  arma::vec beta, lambda, delta;
  int delta_start_idx, delta_length;
  nl_parse_theta(theta, nest_idx, K, use_asc, include_outside_option,
                 beta, lambda, delta, delta_start_idx, delta_length);
  const ChoiceLayout lay = nl_layout_build(X, alt_idx, nest_idx, M, use_asc,
                                           delta, &weights);
  const choicer_off N = lay.N;
  const int* nest = nest_idx.begin(); // alternative a is in nest nest[a] - 1

  const int J_inside = compute_J_inside(use_asc, delta, lay);
  const int J_total = compute_J_total(J_inside, include_outside_option);

  arma::vec base_util = compute_base_util(X, beta, lay, use_asc, delta);

  // numerator(k, j) = sum_i w_i * (-dP_ik/dV_ij), k != j
  // denominator(j)  = sum_i w_i * (dP_ij/dV_ij)
  arma::mat global_numerator = arma::zeros(J_total, J_total);
  arma::vec global_denominator = arma::zeros(J_total);
  const choicer_off max_m = lay.max_m;

#ifdef _OPENMP
#pragma omp parallel
#endif
  {
    arma::mat local_numerator = arma::zeros(J_total, J_total);
    arma::vec local_denominator = arma::zeros(J_total);
    NlProbs pr(max_m, n_nests);         // nl_individual_probs() buffers
    std::vector<int> nest_buf(max_m);   // each inside alternative's nest
    std::vector<int> global_map_buf(max_m); // each one's global index

#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
    for (choicer_off i = 0; i < N; ++i) {
      const int m_i       = lay.m(i);
      const arma::uword start_idx = static_cast<arma::uword>(lay.row_off[i]);
      const double w_i    = weights[i];

      int* nest_idx0_i = nest_buf.data();
      nl_gather_nests(nest_idx0_i, lay.alt + lay.row_off[i], nest, m_i);
      const double* V_inside = base_util.memptr() + start_idx;

      nl_individual_probs(V_inside, nest_idx0_i, m_i, lambda, n_nests,
                          include_outside_option, pr);
      const arma::vec P_i(pr.P_i.data(), m_i, false, true);
      const double* P_j_given_k = pr.P_j_given_k.data();
      const double log_P_outside = pr.log_P_outside;

      const double P_out = include_outside_option ? std::exp(log_P_outside) : 0.0;

      // Global alternative index for each inside alt (outside option = 0)
      int* global_map = global_map_buf.data();
      fill_global_alt_map_inside(global_map, lay.alt + lay.row_off[i], m_i,
                                 include_outside_option);

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
    const Rcpp::IntegerVector& alt_idx,
    const Rcpp::IntegerVector& nest_idx,
    const Rcpp::IntegerVector& M,
    const arma::vec& weights,
    const bool include_outside_option = false,
    const double damping  = 1.0,
    const double tol      = 1e-8,
    const int    max_iter = 1000
) {
  const bool use_asc = true;
  const int n_nests = nl_n_nests(nest_idx);

  const int num_alts =
      include_outside_option ? (delta.n_elem + 1) : delta.n_elem;
  if ((int)target_shares.n_elem != num_alts) {
    Rcpp::stop("Error: target_shares must have the same length as the total number of alternatives.");
  }
  if (arma::any(target_shares <= 0)) {
    Rcpp::stop("Error: all target_shares must be strictly positive (log(share) is undefined otherwise).");
  }
  const ChoiceLayout lay = nl_layout_build(X, alt_idx, nest_idx, M, use_asc,
                                           delta, &weights);

  // nl_parse_theta() is not run here, so check the two invariants that keep
  // nl_individual_probs() in bounds against the lambda given: nest codes from
  // 1, and one lambda per nest. A code below 1 or a lambda of another length
  // used to fail an Armadillo bounds or size check inside the parallel loop
  // (which terminates R under OpenMP); nl_individual_probs() now indexes
  // lambda directly.
  const int nest_min = *std::min_element(nest_idx.begin(), nest_idx.end());
  if (nest_min < 1) {
    Rcpp::stop("nest_idx must use 1-based nest indices (found %s).",
               nest_min == NA_INTEGER ? "NA" : std::to_string(nest_min));
  }
  if (static_cast<int>(lambda.n_elem) != n_nests) {
    Rcpp::stop("lambda must have one entry per nest of nest_idx (%d); it has "
               "%d.", n_nests, lambda.n_elem);
  }

  // X * beta, once: only delta changes across the iterations. Each one adds
  // its ASCs to a copy of it, in one buffer kept across them.
  const double denominator = shares_denominator(weights);
  const arma::vec Xb = X * beta;
  arma::vec base_util;
  auto predict_shares = [&](const arma::vec& inside_delta) {
    base_util = Xb;
    add_row_asc(base_util, lay, inside_delta);
    return nl_predict_shares_internal(base_util, lay, nest_idx.begin(), lambda,
                                      weights, denominator, n_nests, num_alts,
                                      include_outside_option);
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
    delta_old.subvec(1, num_alts - 1) = delta;
    delta_old[0] = 0.0;
  } else {
    delta_old = delta;
    delta_old -= delta_old[0];
  }

  arma::vec inside_delta_old = include_outside_option
    ? arma::vec(delta_old.subvec(1, num_alts - 1))
    : delta_old;

  arma::vec log_shares_old = predict_shares(inside_delta_old);
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
    log_shares_old = predict_shares(inside_delta_old);
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
