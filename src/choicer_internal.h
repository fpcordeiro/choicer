#ifndef CHOICER_INTERNAL_HPP
#define CHOICER_INTERNAL_HPP

#include "choicer.h"

// ============================================================================
// Internal helpers shared by mnlogit.cpp, mxlogit.cpp and nestlogit.cpp.
//
// Every numeric helper body is a verbatim transplant of the code it replaced,
// with identical floating-point operation order, so single-threaded results
// are bit-identical to the pre-refactor implementation. Helpers are
// header-only (inline / templates) so no Makevars or linkage changes are
// needed.
//
// Validation lives in two places, and is always on:
//   * theta-block validation (lengths, K > 0, lambda > 0) inside the theta
//     parsers, so no entry point can parse an inconsistent theta;
//   * data-shape validation (X/W/alt_idx/M/eta/weights consistency) in the
//     validate_*_inputs helpers below, called by every exported entry point.
// Every check is O(1) or a single O(rows) integer scan — negligible next to
// one likelihood evaluation — and turns what would otherwise be an obscure
// Armadillo bounds error (or silently wrong output) into an actionable
// message.
// ============================================================================

// Defined in mxlogit.cpp (exported via Rcpp attributes).
arma::mat build_L_mat(const arma::vec& L_params, const int K_w,
                      const bool rc_correlation);

// ----------------------------------------------------------------------------
// Data-shape validation shared by all exported entry points.
// `delta` is the *full padded* ASC vector returned by the theta parsers; its
// coverage of alt_idx is only checked when use_asc. `weights` / `choice_idx`
// are optional: pass nullptr when the entry point does not take them (or,
// for choice_idx, does not use them).
// ----------------------------------------------------------------------------
inline void validate_choice_data(const arma::mat& X, const arma::uvec& alt_idx,
                                 const Rcpp::IntegerVector& M,
                                 const bool use_asc, const arma::vec& delta,
                                 const arma::vec* weights = nullptr,
                                 const arma::uvec* choice_idx = nullptr) {
  const int N = M.size();
  long long total_rows = 0;
  for (int i = 0; i < N; ++i) {
    if (M[i] <= 0) {
      Rcpp::stop("M must be positive for every individual (M[%d] = %d).",
                 i + 1, M[i]);
    }
    total_rows += M[i];
  }
  if (total_rows != static_cast<long long>(X.n_rows)) {
    Rcpp::stop("X has %d rows but sum(M) is %d.",
               static_cast<int>(X.n_rows), static_cast<int>(total_rows));
  }
  if (alt_idx.n_elem != X.n_rows) {
    Rcpp::stop("alt_idx length (%d) does not match the number of rows of X (%d).",
               static_cast<int>(alt_idx.n_elem), static_cast<int>(X.n_rows));
  }
  if (weights && static_cast<int>(weights->n_elem) != N) {
    Rcpp::stop("weights length (%d) does not match N (%d)",
               weights->n_elem, N);
  }
  if (choice_idx && static_cast<int>(choice_idx->n_elem) != N) {
    Rcpp::stop("choice_idx length (%d) does not match N (%d)",
               static_cast<int>(choice_idx->n_elem), N);
  }
  if (alt_idx.n_elem > 0) {
    if (alt_idx.min() < 1) {
      Rcpp::stop("alt_idx must use 1-based alternative indices (found %d).",
                 static_cast<int>(alt_idx.min()));
    }
    if (use_asc && delta.n_elem < alt_idx.max()) {
      Rcpp::stop("Theta's delta (ASC) block implies %d alternatives but "
                 "alt_idx references alternative %d.",
                 static_cast<int>(delta.n_elem),
                 static_cast<int>(alt_idx.max()));
    }
  }
}

inline void validate_nl_inputs(const arma::mat& X, const arma::uvec& alt_idx,
                               const arma::uvec& nest_idx,
                               const Rcpp::IntegerVector& M,
                               const bool use_asc, const arma::vec& delta,
                               const arma::vec* weights = nullptr,
                               const arma::uvec* choice_idx = nullptr) {
  validate_choice_data(X, alt_idx, M, use_asc, delta, weights, choice_idx);
  if (alt_idx.n_elem > 0 && nest_idx.n_elem < alt_idx.max()) {
    Rcpp::stop("nest_idx has %d entries but alt_idx references alternative %d "
               "(one nest index per global alternative is required).",
               static_cast<int>(nest_idx.n_elem),
               static_cast<int>(alt_idx.max()));
  }
}

inline void validate_mxl_inputs(const arma::mat& X, const arma::mat& W,
                                const arma::uvec& alt_idx,
                                const Rcpp::IntegerVector& M,
                                const arma::cube& eta_draws,
                                const bool use_asc, const arma::vec& delta,
                                const arma::vec* weights = nullptr,
                                const arma::uvec* choice_idx = nullptr) {
  validate_choice_data(X, alt_idx, M, use_asc, delta, weights, choice_idx);
  const int N = M.size();
  const int K_w = W.n_cols;
  if (static_cast<int>(eta_draws.n_slices) != N) {
    Rcpp::stop("eta_draws 3rd dimension (%d) does not match N (%d)",
               eta_draws.n_slices, N);
  }
  if (static_cast<int>(eta_draws.n_rows) != K_w) {
    Rcpp::stop("eta_draws 1st dimension (%d) does not match K_w (%d)",
               eta_draws.n_rows, K_w);
  }
  if (W.n_rows != X.n_rows && alt_idx.n_elem > 0 &&
      W.n_rows < alt_idx.max()) {
    Rcpp::stop("W must be row-aligned with X (%d rows) or contain one row per "
               "global alternative (at least %d rows); got %d rows.",
               static_cast<int>(X.n_rows), static_cast<int>(alt_idx.max()),
               static_cast<int>(W.n_rows));
  }
}

inline void check_rc_dist_length(const arma::uvec& rc_dist, const int K_w) {
  if (static_cast<int>(rc_dist.n_elem) != K_w) {
    Rcpp::stop("rc_dist must be a vector of length K_w (%d)", K_w);
  }
}

// ----------------------------------------------------------------------------
// MNL theta parsing: theta = [beta (K), delta (J or J-1)].
// delta is returned as the *full* padded vector: when there is no outside
// option the first inside alternative's ASC is fixed to 0 and prepended.
// ----------------------------------------------------------------------------
struct MnlParams {
  arma::vec beta;
  arma::vec delta;
};

inline MnlParams parse_mnl_theta(const arma::vec& theta, const int K,
                                 const bool use_asc,
                                 const bool include_outside_option) {
  const int n_params = theta.n_elem;
  if (K <= 0) {
    Rcpp::stop("K must be positive, got %d", K);
  }
  if (n_params < K) {
    Rcpp::stop("Theta vector too short: missing beta parameters "
               "(expected at least %d, got %d).", K, n_params);
  }
  if (!use_asc && n_params != K) {
    Rcpp::stop("Theta vector too long: %d parameters given but the model "
               "expects %d. Did you mean use_asc = TRUE?", n_params, K);
  }
  MnlParams P;
  P.beta = theta.subvec(0, K - 1);

  if (use_asc) {
    const int delta_length = n_params - K;
    if (delta_length <= 0) {
      Rcpp::stop("Error: ASC parameters expected but not provided.");
    }
    if (include_outside_option) {
      // delta covers all J inside alternatives
      P.delta = theta.subvec(K, n_params - 1);
    } else {
      // delta_1 = 0 fixed
      P.delta = arma::zeros(delta_length + 1);
      P.delta.subvec(1, delta_length) = theta.subvec(K, n_params - 1);
    }
  } else {
    P.delta = arma::zeros(0);
  }
  return P;
}

// ----------------------------------------------------------------------------
// MXL theta parsing: theta = [beta (K_x), mu (K_w, if rc_mean), L (L_size),
// delta]. Returns block start indices, beta, the transformed mu (mu_final =
// exp(mu_k) for log-normal coefficients) with its first and second
// derivatives, the rebuilt Cholesky factor L, and the padded delta.
//
// Derivative-vector defaults (dmu_final_dmu = ones, dmu2_final_dmu2 = zeros)
// are standardized across callers: every read site in mxlogit.cpp (gradient,
// hessian, bhhh) sits inside an `if (rc_mean)` block, where the values below
// match what each function computed before; when !rc_mean they are never read.
//
// All theta-block validation (K_x > 0, rc_dist length, per-block theta
// lengths, no trailing unparsed parameters) happens here unconditionally, so
// every entry point — including the predict family — fails with the same
// actionable message on a malformed theta.
// ----------------------------------------------------------------------------
struct MxlParams {
  int idx_beta_start, idx_mu_start, idx_L_start, idx_delta_start, L_size;
  arma::vec beta;
  arma::vec mu_final;
  arma::vec dmu_final_dmu;
  arma::vec dmu2_final_dmu2;
  arma::mat L;
  arma::vec delta;
};

inline MxlParams parse_mxl_theta(const arma::vec& theta,
                                 const int K_x, const int K_w,
                                 const arma::uvec& rc_dist,
                                 const bool rc_correlation, const bool rc_mean,
                                 const bool use_asc,
                                 const bool include_outside_option) {
  const int n_params = theta.n_elem;
  if (K_x <= 0) {
    Rcpp::stop("K_x must be positive, got %d", K_x);
  }
  check_rc_dist_length(rc_dist, K_w);

  MxlParams P;
  P.L_size = rc_correlation ? (K_w * (K_w + 1)) / 2 : K_w;
  P.idx_beta_start = 0;
  P.idx_mu_start = K_x;
  P.idx_L_start = rc_mean ? K_x + K_w : K_x;
  P.idx_delta_start = P.idx_L_start + P.L_size;

  // beta: coefficients for design matrix X
  if (n_params < P.idx_mu_start)
    Rcpp::stop("Theta vector too short: missing beta parameters "
               "(expected at least %d, got %d).", K_x, n_params);
  P.beta = theta.subvec(P.idx_beta_start, P.idx_mu_start - 1);

  // mu: means of random coefficients (zeros when not estimated)
  arma::vec mu;
  if (rc_mean) {
    if (P.idx_L_start > n_params)
      Rcpp::stop("Theta vector too short: missing mu parameters.");
    mu = theta.subvec(P.idx_mu_start, P.idx_L_start - 1);
  } else {
    mu = arma::zeros(K_w);
  }

  // L: choleski decomposition of random coefficients matrix
  if (P.idx_delta_start > n_params)
    Rcpp::stop("Theta vector too short: missing L parameters.");
  if (!use_asc && n_params != P.idx_delta_start)
    Rcpp::stop("Theta vector too long: %d parameters given but the model "
               "expects %d. Did you mean use_asc = TRUE?",
               n_params, P.idx_delta_start);
  arma::vec L_params = theta.subvec(P.idx_L_start, P.idx_delta_start - 1);
  P.L = build_L_mat(L_params, K_w, rc_correlation);

  // mu transformations for distributions
  P.mu_final = arma::zeros(K_w);
  P.dmu_final_dmu = arma::ones(K_w);
  P.dmu2_final_dmu2 = arma::zeros(K_w);
  if (rc_mean) {
    P.mu_final = mu;
    for (int k = 0; k < K_w; ++k) {
      if (rc_dist(k) == 1) { // 1 == log-normal
        P.mu_final(k) = std::exp(mu(k));
        P.dmu_final_dmu(k) = P.mu_final(k);   // d(exp(mu))/dmu = exp(mu)
        P.dmu2_final_dmu2(k) = P.mu_final(k); // d^2(exp(mu))/dmu^2 = exp(mu)
      }
    }
  }

  // delta (ASC)
  if (use_asc) {
    const int delta_free_len = n_params - P.idx_delta_start;
    if (delta_free_len <= 0) {
      Rcpp::stop("Theta vector too short: missing delta parameters.");
    }
    if (include_outside_option) {
      // all inside alternatives are free
      P.delta = theta.subvec(P.idx_delta_start, n_params - 1);
    } else {
      // first delta is fixed to 0 -> it's not in theta
      P.delta = arma::zeros(delta_free_len + 1);
      P.delta.subvec(1, delta_free_len) =
          theta.subvec(P.idx_delta_start, n_params - 1);
    }
  } else {
    P.delta.set_size(0); // empty
  }
  return P;
}

// ----------------------------------------------------------------------------
// Base utility, pre-computed for all stacked rows with single BLAS calls:
// base_util = X*beta (+ W*mu_final) (+ delta scattered by alternative).
// The MXL overload handles both W layouts: row-aligned with X
// (sum(M) x K_w) or one row per global alternative (J x K_w).
// ----------------------------------------------------------------------------
inline arma::vec compute_base_util(const arma::mat& X, const arma::vec& beta,
                                   const arma::uvec& alt_idx0,
                                   const bool use_asc, const arma::vec& delta) {
  arma::vec base_util = X * beta;
  if (use_asc) base_util += delta.elem(alt_idx0);
  return base_util;
}

inline arma::vec compute_base_util_mxl(const arma::mat& X, const arma::mat& W,
                                       const arma::vec& beta,
                                       const arma::vec& mu_final,
                                       const arma::uvec& alt_idx0,
                                       const bool use_asc,
                                       const arma::vec& delta) {
  arma::vec base_util = X * beta;
  if (static_cast<int>(W.n_rows) == static_cast<int>(X.n_rows)) {
    base_util += W * mu_final;
  } else {
    arma::vec W_mu = W * mu_final;
    base_util += W_mu.elem(alt_idx0);
  }
  if (use_asc) base_util += delta.elem(alt_idx0);
  return base_util;
}

// ----------------------------------------------------------------------------
// Per-individual slice of the random-coefficient design matrix W:
// row-aligned with X -> contiguous row block; alt-level W -> gather rows by
// this individual's alternative indices. Templated on the index type so both
// zero-copy subviews and materialized uvec indices forward without copies.
// ----------------------------------------------------------------------------
template <typename IdxT>
inline arma::mat make_W_i(const arma::mat& W, const arma::uword x_n_rows,
                          const arma::uword start_idx,
                          const arma::uword end_idx,
                          const IdxT& alt_idx0_i) {
  if (W.n_rows == x_n_rows)            // row-aligned with X
    return W.rows(start_idx, end_idx); // m_i x K_w
  return W.rows(alt_idx0_i);           // global alt-level W
}

// ----------------------------------------------------------------------------
// Batched Cholesky draws: Gamma_final = L * eta_i in a single dgemm, then the
// log-normal transform applied row-wise where rc_dist == 1. Optional outputs
// Dgamma1/Dgamma2 receive the first/second derivative of the transform
// (ones/zeros for normal coefficients, exp(L*eta) rows for log-normal).
// ----------------------------------------------------------------------------
inline arma::mat batch_gamma_draws(const arma::mat& L, const arma::mat& eta_i,
                                   const arma::uvec& rc_dist,
                                   arma::mat* Dgamma1 = nullptr,
                                   arma::mat* Dgamma2 = nullptr) {
  arma::mat Gamma_final = L * eta_i; // single dgemm
  if (Dgamma1) Dgamma1->ones(L.n_rows, eta_i.n_cols);
  if (Dgamma2) Dgamma2->zeros(L.n_rows, eta_i.n_cols);
  for (arma::uword k = 0; k < L.n_rows; ++k) {
    if (rc_dist(k) == 1) { // log-normal
      Gamma_final.row(k) = arma::exp(Gamma_final.row(k));
      if (Dgamma1) Dgamma1->row(k) = Gamma_final.row(k);
      if (Dgamma2) Dgamma2->row(k) = Gamma_final.row(k);
    }
  }
  return Gamma_final;
}

// ----------------------------------------------------------------------------
// Utility-vector fill and numerically stable softmax.
// fill_choice_utilities places the inside utilities into the caller-owned V
// (the outside option, when present, occupies slot 0 with V = 0).
// stable_softmax stabilizes V in place (max subtraction), writes the choice
// probabilities into the caller-owned P, and returns log_denom.
// Note: arma::sum on a vector dispatches to arma::accu, so unifying the
// historical sum/accu mix on accu is bit-identical.
// ----------------------------------------------------------------------------
inline void fill_choice_utilities(arma::vec& V, const arma::vec& inside_utils,
                                  const int num_choices,
                                  const bool include_outside_option) {
  V.zeros();
  if (include_outside_option)
    V.subvec(1, num_choices - 1) = inside_utils;
  else
    V = inside_utils;
}

inline double stable_softmax(arma::vec& V, arma::vec& P) {
  V -= V.max(); // for numerical stability
  const double log_denom = std::log(arma::accu(arma::exp(V)));
  P = arma::exp(V - log_denom);
  return log_denom;
}

// ----------------------------------------------------------------------------
// Delta-block (ASC) gradient scatter, in inside-alternative space:
// diff_inside has one entry per inside alternative (callers with an outside
// option pass diff_vec.subvec(1, m_i), a zero-copy subview). When there is no
// outside option, the first inside alternative's ASC is the fixed
// normalization and receives no gradient.
// ----------------------------------------------------------------------------
template <typename DiffT, typename IdxT>
inline void scatter_delta_grad(arma::vec& g, const int delta_start,
                               const DiffT& diff_inside,
                               const IdxT& alt_idx0_i, const int m_i,
                               const bool include_outside_option,
                               const double scale) {
  for (int j = 0; j < m_i; ++j) {
    const int id = static_cast<int>(alt_idx0_i[j]);
    if (include_outside_option) {
      g[delta_start + id] += scale * diff_inside[j];
    } else if (id > 0) { // delta of first inside alt is normalised to 0
      g[delta_start + (id - 1)] += scale * diff_inside[j];
    }
  }
}

// ----------------------------------------------------------------------------
// Map local choice-set indices to global alternative indices for the
// J_total x J_total output matrices (elasticities, diversion ratios).
// Full variant: index 0 = outside option, inside alts shifted by +1.
// Inside variant (NL): m_i-length map over inside alternatives only.
// ----------------------------------------------------------------------------
inline arma::uvec build_global_alt_map(const arma::uvec& alt_idx0_i,
                                       const int m_i,
                                       const bool include_outside_option) {
  arma::uvec global_j_map(include_outside_option ? m_i + 1 : m_i);
  if (include_outside_option) {
    global_j_map[0] = 0;                       // outside option = global index 0
    global_j_map.subvec(1, m_i) = alt_idx0_i + 1; // inside alts are 1...J
  } else {
    global_j_map = alt_idx0_i;                 // no outside option: 0...J-1
  }
  return global_j_map;
}

inline arma::uvec build_global_alt_map_inside(const arma::uvec& alt_idx0_i,
                                              const bool include_outside_option) {
  arma::uvec global_map(alt_idx0_i.n_elem);
  if (include_outside_option) {
    global_map = alt_idx0_i + 1;
  } else {
    global_map = alt_idx0_i;
  }
  return global_map;
}

// ----------------------------------------------------------------------------
// Output-matrix dimensions: number of inside alternatives and total
// alternatives (including the outside option when present). arma::max is only
// evaluated when !use_asc, as in every historical copy.
// ----------------------------------------------------------------------------
inline int compute_J_inside(const bool use_asc, const arma::vec& delta,
                            const arma::uvec& alt_idx0) {
  return use_asc ? static_cast<int>(delta.n_elem)
                 : (static_cast<int>(arma::max(alt_idx0)) + 1);
}

inline int compute_J_total(const int J_inside,
                           const bool include_outside_option) {
  return include_outside_option ? J_inside + 1 : J_inside;
}

#endif // CHOICER_INTERNAL_HPP
