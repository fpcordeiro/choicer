#ifndef CHOICER_INTERNAL_HPP
#define CHOICER_INTERNAL_HPP

#include "choicer.h"
#include <algorithm>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

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
//   * data-shape validation (X/W/alt_idx/M/eta/weights consistency), called
//     by every exported entry point: in the layout builders (MNL:
//     choice_layout_build below; MXL estimation: mxl_layout_build) and, for
//     the kernels that still take arma::uvec indices, in the
//     validate_*_inputs helpers.
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

// eta_draws holds one K_w x S draw block per likelihood unit: per choice
// situation in the cross-section (n_units < 0, the default), per decision
// maker in a panel (n_units = number of decision makers).
inline void validate_mxl_inputs(const arma::mat& X, const arma::mat& W,
                                const arma::uvec& alt_idx,
                                const Rcpp::IntegerVector& M,
                                const arma::cube& eta_draws,
                                const bool use_asc, const arma::vec& delta,
                                const arma::vec* weights = nullptr,
                                const arma::uvec* choice_idx = nullptr,
                                const int n_units = -1) {
  validate_choice_data(X, alt_idx, M, use_asc, delta, weights, choice_idx);
  const int N = M.size();
  const int K_w = W.n_cols;
  if (n_units < 0 && static_cast<int>(eta_draws.n_slices) != N) {
    Rcpp::stop("eta_draws 3rd dimension (%d) does not match N (%d)",
               eta_draws.n_slices, N);
  }
  if (n_units >= 0 && static_cast<int>(eta_draws.n_slices) != n_units) {
    Rcpp::stop("eta_draws 3rd dimension (%d) does not match the number of "
               "decision makers (%d)", eta_draws.n_slices, n_units);
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
// Layout of the stacked design
//
// The theta-independent layout of a kernel call: situation offsets, and the
// per-row alternative codes and per-situation choices, read in place from the
// kernel's integer arguments (no per-call copy). Offsets are 64-bit
// (choicer_off). R caps a matrix at 2^31 - 1 rows, so the row offsets fit in
// an int; the element offsets formed from them (row + column * n_rows) are
// Armadillo's (arma::uword). ChoiceLayout is the core shared by the MNL
// kernels (choice_layout_build) and the MXL estimation kernels (MxlLayout,
// mxl_layout_build). It is valid for one kernel call only: when an argument
// was coerced from doubles, the integer copy belongs to that call's Rcpp
// wrapper. Build it on the primary thread; it is read-only afterwards.
// ----------------------------------------------------------------------------
using choicer_off = std::ptrdiff_t;

// A pass over the stacked rows runs in parallel above this many rows; below
// it a thread team costs more than the pass. Each such pass is elementwise or
// an exact min/max, so the thread count never changes its result.
constexpr choicer_off CHOICER_PARALLEL_ROWS = 1000000;

struct ChoiceLayout {
  choicer_off n_rows = 0;            // stacked alternative rows, sum(M)
  choicer_off N = 0;                 // choice situations
  choicer_off max_m = 0;             // rows of the largest situation
  int J = 0;                         // largest alternative code
  std::vector<choicer_off> row_off;  // situation t: rows [row_off[t], row_off[t+1])
  const int* alt = nullptr;          // row r: 1-based alternative code
  const int* choice = nullptr;       // situation t: the kernel's choice_idx

  // Rows of situation t (M[t]), and the 0-based alternative of row r
  int m(const choicer_off t) const {
    return static_cast<int>(row_off[t + 1] - row_off[t]);
  }
  int alt0(const choicer_off r) const { return alt[r] - 1; }
};

// 0-based alternative codes of consecutive rows, read from the 1-based codes
// in place: what scatter_delta_grad() and the per-situation loops index.
struct AltCodes0 {
  const int* alt; // 1-based codes from the first row
  int operator[](const int j) const { return alt[j] - 1; }
};

// The situations: row offsets, summed in 64 bits, and the largest choice set;
// then the design's height and the lengths of alt_idx, weights and
// choice_idx, with validate_choice_data()'s messages in its order. Pass
// nullptr for weights or choice_idx where the kernel takes none, or (MNL)
// where validate_choice_data() was not given them.
inline void choice_layout_situations(ChoiceLayout& lay, const arma::mat& X,
                                     const Rcpp::IntegerVector& alt_idx,
                                     const Rcpp::IntegerVector& M,
                                     const arma::vec* weights,
                                     const Rcpp::IntegerVector* choice_idx) {
  const choicer_off N = M.size();
  const int* m = M.begin();
  lay.N = N;
  lay.row_off.assign(static_cast<std::size_t>(N) + 1, 0);
  for (choicer_off t = 0; t < N; ++t) {
    if (m[t] <= 0) { // NA_INTEGER is negative
      Rcpp::stop("M must be positive for every individual (M[%d] = %d).",
                 t + 1, m[t]);
    }
    lay.row_off[t + 1] = lay.row_off[t] + m[t];
    lay.max_m = std::max<choicer_off>(lay.max_m, m[t]);
  }
  lay.n_rows = lay.row_off[N];
  if (lay.n_rows != static_cast<choicer_off>(X.n_rows)) {
    Rcpp::stop("X has %d rows but sum(M) is %d.", X.n_rows, lay.n_rows);
  }
  if (static_cast<choicer_off>(alt_idx.size()) != lay.n_rows) {
    Rcpp::stop("alt_idx length (%d) does not match the number of rows of X "
               "(%d).", alt_idx.size(), X.n_rows);
  }
  if (weights && static_cast<choicer_off>(weights->n_elem) != N) {
    Rcpp::stop("weights length (%d) does not match N (%d)", weights->n_elem, N);
  }
  if (choice_idx && static_cast<choicer_off>(choice_idx->size()) != N) {
    Rcpp::stop("choice_idx length (%d) does not match N (%d)",
               choice_idx->size(), N);
  }
  lay.alt = alt_idx.begin();
  lay.choice = choice_idx ? choice_idx->begin() : nullptr;
}

// The alternative codes, in one min/max pass over the rows: they are 1-based,
// and NA_INTEGER, the most negative int, fails the same check and is reported
// as NA. (The arma::uvec arguments this replaces read NA and negative codes
// through a double-to-unsigned cast, which is undefined.) J is the largest
// code, 0 without rows.
inline void choice_layout_codes(ChoiceLayout& lay) {
  int a_min = std::numeric_limits<int>::max(), a_max = 0;
  const int* alt = lay.alt;
  const choicer_off n_rows = lay.n_rows;
#ifdef _OPENMP
#pragma omp parallel for reduction(min : a_min) reduction(max : a_max) \
    if (n_rows > CHOICER_PARALLEL_ROWS)
#endif
  for (choicer_off r = 0; r < n_rows; ++r) {
    a_min = std::min(a_min, alt[r]);
    a_max = std::max(a_max, alt[r]);
  }
  if (n_rows > 0 && a_min < 1) {
    Rcpp::stop("alt_idx must use 1-based alternative indices (found %s).",
               a_min == NA_INTEGER ? "NA" : std::to_string(a_min));
  }
  lay.J = a_max;
}

// Every situation's choice is a slot of its choice set: 1..M[t], or 0 for the
// outside option when there is one. The slot is formed in 64 bits, so no
// choice, NA_INTEGER included, overflows it; NA is rejected by name. `kernel`
// (or nullptr) names the kernel in the message, as the BHHH, score and MXL
// kernels have always done.
inline void validate_choices(const ChoiceLayout& lay,
                             const bool include_outside_option,
                             const char* kernel) {
  for (choicer_off t = 0; t < lay.N; ++t) {
    const int c = lay.choice[t];
    const choicer_off slot =
        include_outside_option ? choicer_off(c) : choicer_off(c) - 1;
    const choicer_off n_choices =
        include_outside_option ? choicer_off(lay.m(t)) + 1 : lay.m(t);
    if (c == NA_INTEGER || slot < 0 || slot >= n_choices) {
      if (kernel) {
        Rcpp::stop("Invalid chosen alternative index for individual %d (%s)",
                   t, kernel);
      }
      Rcpp::stop("Invalid chosen alternative index for individual %d", t);
    }
  }
}

// The layout of the MNL kernels, with validate_choice_data()'s checks in its
// order and with its messages: the situations, the alternative codes, and
// the delta block's coverage of them. choice_idx is only length-checked here
// and kept for validate_choices(): pass nullptr where the kernel takes no
// choices or does not check them.
inline ChoiceLayout choice_layout_build(
    const arma::mat& X, const Rcpp::IntegerVector& alt_idx,
    const Rcpp::IntegerVector& M, const bool use_asc, const arma::vec& delta,
    const arma::vec* weights = nullptr,
    const Rcpp::IntegerVector* choice_idx = nullptr) {
  ChoiceLayout lay;
  choice_layout_situations(lay, X, alt_idx, M, weights, choice_idx);
  choice_layout_codes(lay);
  if (use_asc && lay.n_rows > 0 &&
      static_cast<choicer_off>(delta.n_elem) < lay.J) {
    Rcpp::stop("Theta's delta (ASC) block implies %d alternatives but "
               "alt_idx references alternative %d.", delta.n_elem, lay.J);
  }
  return lay;
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

// The same with the alternative codes read in place: X * beta in the same
// single BLAS call, then each row's ASC added to its base utility, the one
// addition per element that += delta.elem(alt_idx0) makes (in parallel above
// 10^6 rows; every element is independent, so the sums are unchanged).
// add_row_asc() serves callers that keep X * beta across calls (BLP).
inline void add_row_asc(arma::vec& base_util, const ChoiceLayout& lay,
                        const arma::vec& delta) {
  double* bu = base_util.memptr();
  const double* d = delta.memptr();
  const int* alt = lay.alt;
  const choicer_off n_rows = lay.n_rows;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) \
    if (n_rows > CHOICER_PARALLEL_ROWS)
#endif
  for (choicer_off r = 0; r < n_rows; ++r) {
    bu[r] += d[alt[r] - 1];
  }
}

inline arma::vec compute_base_util(const arma::mat& X, const arma::vec& beta,
                                   const ChoiceLayout& lay,
                                   const bool use_asc, const arma::vec& delta) {
  arma::vec base_util = X * beta;
  if (use_asc) add_row_asc(base_util, lay, delta);
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
// The _into form writes into a caller-owned (reused) Gamma_final.
// ----------------------------------------------------------------------------
inline void batch_gamma_draws_into(arma::mat& Gamma_final, const arma::mat& L,
                                   const arma::mat& eta_i,
                                   const arma::uvec& rc_dist,
                                   arma::mat* Dgamma1 = nullptr,
                                   arma::mat* Dgamma2 = nullptr) {
  Gamma_final = L * eta_i; // single dgemm
  if (Dgamma1) Dgamma1->ones(L.n_rows, eta_i.n_cols);
  if (Dgamma2) Dgamma2->zeros(L.n_rows, eta_i.n_cols);
  for (arma::uword k = 0; k < L.n_rows; ++k) {
    if (rc_dist(k) == 1) { // log-normal, in place (no aliasing temporary)
      for (arma::uword s = 0; s < Gamma_final.n_cols; ++s) {
        Gamma_final(k, s) = std::exp(Gamma_final(k, s));
      }
      if (Dgamma1) Dgamma1->row(k) = Gamma_final.row(k);
      if (Dgamma2) Dgamma2->row(k) = Gamma_final.row(k);
    }
  }
}

inline arma::mat batch_gamma_draws(const arma::mat& L, const arma::mat& eta_i,
                                   const arma::uvec& rc_dist,
                                   arma::mat* Dgamma1 = nullptr,
                                   arma::mat* Dgamma2 = nullptr) {
  arma::mat Gamma_final;
  batch_gamma_draws_into(Gamma_final, L, eta_i, rc_dist, Dgamma1, Dgamma2);
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
  V -= V.max(); // for numerical stability; after this, max(V) == 0, so exp(V) <= 1
  P = arma::exp(V);       // reuse caller-owned P as the exp buffer — no per-call alloc
  const double s = arma::accu(P);
  P /= s;                 // in-place; bitwise identical to previous `P = e / s`
  return std::log(s);
}

// ----------------------------------------------------------------------------
// stable_softmax() and logSumExp() on raw arrays, for loops that keep one
// buffer at the largest size they need. The operations are Armadillo's, in
// Armadillo's order: the paired max scan of op_max::direct_max, the shift and
// exp, the two-accumulator sum of arrayops::accumulate, the division. The
// results are therefore bitwise those of stable_softmax() and logSumExp()
// (unless Armadillo is built with -ffast-math, where it sums with a single
// accumulator), without a vector object per call and without the temporary
// that accu() makes of an exp() expression when Armadillo uses OpenMP.
// ----------------------------------------------------------------------------
inline double direct_max_n(const double* x, const int n) {
  double max_i = -arma::datum::inf, max_j = -arma::datum::inf;
  int i, j;
  for (i = 0, j = 1; j < n; i += 2, j += 2) {
    if (x[i] > max_i) max_i = x[i];
    if (x[j] > max_j) max_j = x[j];
  }
  if (i < n && x[i] > max_i) max_i = x[i];
  return (max_i > max_j) ? max_i : max_j;
}

// stable_softmax() of the n >= 1 entries of v (shifted in place), with the
// probabilities in p; returns the log of the denominator.
inline double stable_softmax_n(double* v, double* p, const int n) {
  const double v_max = direct_max_n(v, n);
  for (int i = 0; i < n; ++i) v[i] -= v_max;
  for (int i = 0; i < n; ++i) p[i] = std::exp(v[i]);
  double acc1 = 0.0, acc2 = 0.0;
  int j;
  for (j = 1; j < n; j += 2) {
    acc1 += p[j - 1];
    acc2 += p[j];
  }
  if (j - 1 < n) acc1 += p[j - 1];
  const double sum = acc1 + acc2;
  for (int i = 0; i < n; ++i) p[i] /= sum;
  return std::log(sum);
}

// logSumExp() of the n entries of x, with its handling of empty input, NaN
// and infinities.
inline double log_sum_exp_n(const double* x, const int n) {
  if (n == 0) return -arma::datum::inf;
  const double a = direct_max_n(x, n);
  if (std::isnan(a)) return arma::datum::nan;
  if (a == arma::datum::inf) return arma::datum::inf;
  if (a == -arma::datum::inf) return -arma::datum::inf;
  double acc1 = 0.0, acc2 = 0.0;
  int j;
  for (j = 1; j < n; j += 2) {
    acc1 += std::exp(x[j - 1] - a);
    acc2 += std::exp(x[j] - a);
  }
  if (j - 1 < n) acc1 += std::exp(x[j - 1] - a);
  return a + std::log(acc1 + acc2);
}

// ----------------------------------------------------------------------------
// Nested logit: one individual's probabilities
//
// Given one individual's inside utilities V (m entries), the 0-based nest of
// each inside alternative (nest, m entries), the *full* lambda vector
// (n_nests entries, singletons fixed to 1) and the outside-option flag,
// nl_individual_probs() fills in pr:
//   P_i           (m)        joint choice probability  P_ij = P(j|k) * P_k
//   P_j_given_k   (m)        within-nest conditional probability P(j|k)
//   P_k           (n_nests)  marginal nest probability P_k
//   log_I_k       (n_nests)  log inclusive value of each nest (-inf if empty)
//   log_P_i       (m)        log joint choice probability (stabilised)
//   log_P_outside            log probability of the outside option (-inf if
//                            there is none)
// with the same two-level log-sum-exp stabilisation as the likelihood kernel
// (the outside option has V = 0 and lambda = 1, so its nest term is 0).
//
// NlProbs holds the outputs and the working arrays: one per thread, sized
// once for the largest choice set and the number of nests, so the helper
// allocates nothing per individual. The arithmetic is that of the Armadillo
// expressions this replaced, operation for operation and in their order:
// elementwise quotients, products, differences and exp() of the same
// operands; the maximum nest term by op_max::direct_max's paired scan
// (direct_max_n); and the sum of their exponentials by
// arrayops::accumulate's two accumulators, which is how accu() sums an exp()
// expression whether or not it first materializes it (it does when Armadillo
// uses OpenMP; under -ffast-math it sums in one accumulator, but no order is
// fixed there). lambda_k log I_k is formed once, into nest_terms, and every
// later use reads that rounded product, never a fused multiply-add of it. The
// results are therefore bitwise those of the Armadillo version, which
// test_nl_individual_probs() keeps for comparison (kernel_test_exports.cpp).
// ----------------------------------------------------------------------------
struct NlProbs {
  std::vector<double> P_i, P_j_given_k, log_P_i; // m (sized for the largest)
  std::vector<double> P_k, log_I_k;              // n_nests
  double log_P_outside = 0.0;
  // working arrays
  std::vector<double> V_over_lambda, log_P_j_given_k;             // m
  std::vector<double> max_V_k, I_k_unscaled, nest_terms, log_P_k; // n_nests

  NlProbs(const choicer_off max_m, const int n_nests)
      : P_i(max_m), P_j_given_k(max_m), log_P_i(max_m), P_k(n_nests),
        log_I_k(n_nests), V_over_lambda(max_m), log_P_j_given_k(max_m),
        max_V_k(n_nests), I_k_unscaled(n_nests), nest_terms(n_nests),
        log_P_k(n_nests) {}
};

inline void nl_individual_probs(const double* V, const int* nest, const int m,
                                const arma::vec& lambda, const int n_nests,
                                const bool include_outside_option,
                                NlProbs& pr) {
  const double* lam = lambda.memptr();
  const double inf = arma::datum::inf;

  // V_ij / lambda_k  (lambda_k = 1 for singletons)
  double* V_over_lambda = pr.V_over_lambda.data();
  for (int j = 0; j < m; ++j) V_over_lambda[j] = V[j] / lam[nest[j]];

  // --- log_I_k (inclusive value) via log-sum-exp within each nest ---
  double* max_V_k = pr.max_V_k.data();
  for (int k = 0; k < n_nests; ++k) max_V_k[k] = -inf;
  for (int j = 0; j < m; ++j) {
    const int k = nest[j];
    if (V_over_lambda[j] > max_V_k[k]) {
      max_V_k[k] = V_over_lambda[j];
    }
  }

  double* I_k_unscaled = pr.I_k_unscaled.data();
  for (int k = 0; k < n_nests; ++k) I_k_unscaled[k] = 0.0;
  for (int j = 0; j < m; ++j) {
    const int k = nest[j];
    if (std::isfinite(max_V_k[k])) {
      I_k_unscaled[k] += std::exp(V_over_lambda[j] - max_V_k[k]);
    }
  }

  double* log_I_k = pr.log_I_k.data();
  for (int k = 0; k < n_nests; ++k) {
    log_I_k[k] = -inf;
    if (I_k_unscaled[k] > 0) {
      log_I_k[k] = max_V_k[k] + std::log(I_k_unscaled[k]);
    }
  }

  // --- log(P_k) (nest probability) ---
  double* nest_terms = pr.nest_terms.data();
  for (int k = 0; k < n_nests; ++k) nest_terms[k] = lam[k] * log_I_k[k];

  double max_nest_term = direct_max_n(nest_terms, n_nests);
  if (!std::isfinite(max_nest_term)) {
    max_nest_term = 0;
  }

  double acc1 = 0.0, acc2 = 0.0;
  int k2;
  for (k2 = 1; k2 < n_nests; k2 += 2) {
    acc1 += std::exp(nest_terms[k2 - 1] - max_nest_term);
    acc2 += std::exp(nest_terms[k2] - max_nest_term);
  }
  if (k2 - 1 < n_nests) acc1 += std::exp(nest_terms[k2 - 1] - max_nest_term);
  double sum_exp_nest_terms = acc1 + acc2;
  if (include_outside_option) {
    // Outside option: V=0, lambda=1 -> term = 0
    sum_exp_nest_terms += std::exp(0.0 - max_nest_term);
  }

  const double log_denom_P_nest = max_nest_term + std::log(sum_exp_nest_terms);

  double* log_P_k = pr.log_P_k.data();
  for (int k = 0; k < n_nests; ++k) {
    log_P_k[k] = nest_terms[k] - log_denom_P_nest;
  }
  pr.log_P_outside = include_outside_option ? (0.0 - log_denom_P_nest) : -inf;

  // --- log(P_j|k) and log(P_ij) ---
  double* log_P_j_given_k = pr.log_P_j_given_k.data();
  double* log_P_i = pr.log_P_i.data();
  for (int j = 0; j < m; ++j) {
    log_P_j_given_k[j] = V_over_lambda[j] - log_I_k[nest[j]];
  }
  for (int j = 0; j < m; ++j) {
    log_P_i[j] = log_P_j_given_k[j] + log_P_k[nest[j]];
  }

  double* P_i = pr.P_i.data();
  double* P_j_given_k = pr.P_j_given_k.data();
  double* P_k = pr.P_k.data();
  for (int j = 0; j < m; ++j) P_i[j] = std::exp(log_P_i[j]);
  for (int j = 0; j < m; ++j) P_j_given_k[j] = std::exp(log_P_j_given_k[j]);
  for (int k = 0; k < n_nests; ++k) P_k[k] = std::exp(log_P_k[k]);
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

// build_global_alt_map() from the 1-based codes of a situation's m rows, into
// a caller-owned buffer of at least m + 1 entries (one per thread, sized once
// from the layout's max_m).
inline void fill_global_alt_map(int* map, const int* alt, const int m,
                                const bool include_outside_option) {
  if (include_outside_option) {
    map[0] = 0;                              // outside option = global index 0
    for (int j = 0; j < m; ++j) map[j + 1] = alt[j]; // inside alts are 1...J
  } else {
    for (int j = 0; j < m; ++j) map[j] = alt[j] - 1; // no outside: 0...J-1
  }
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

// The same from the layout. arma::max() of an empty alt_idx0 threw
// std::logic_error "max(): object has no elements", which is kept.
inline int compute_J_inside(const bool use_asc, const arma::vec& delta,
                            const ChoiceLayout& lay) {
  if (use_asc) return static_cast<int>(delta.n_elem);
  if (lay.n_rows == 0) throw std::logic_error("max(): object has no elements");
  return lay.J;
}

// Formed in 64 bits and narrowed: the largest code, 2^31 - 1, with an
// outside option wraps to a negative count instead of overflowing int.
inline int compute_J_total(const int J_inside,
                           const bool include_outside_option) {
  return static_cast<int>(static_cast<choicer_off>(J_inside) +
                          (include_outside_option ? 1 : 0));
}

#endif // CHOICER_INTERNAL_HPP
