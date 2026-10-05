// [[Rcpp::depends(RcppArmadillo)]]
#include "choicer.h"
#include "choicer_internal.h"
#include "halton.h"
#include <algorithm>
#include <cstdint>
#include <cstring>
#include <limits>
#include <memory>
#include <vector>

// Reconstruct lower-triangular choleski factor L from L_params
// [[Rcpp::export]]
arma::mat build_L_mat(const arma::vec &L_params, const int K_w,
                      const bool rc_correlation) {
  arma::mat L = arma::zeros(K_w, K_w);
  int idx = 0;
  if (rc_correlation) { // full (lower-triangular) factor
    for (int i = 0; i < K_w; ++i) {
      for (int j = 0; j <= i; ++j, ++idx) {
        double val = L_params(idx);
        if (i == j) { // diagonal - exp()
          L(i, j) = std::exp(val);
        } else {
          L(i, j) = val; // off-diagonal stays unconstrained
        }
      }
    }
  } else { // diagonal-only (Sigma is diagonal)
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
//' @returns matrix equal to LL', where L is the choleski decomposition of random coefficient matrix
//' @examples
//' L_params <- c(log(1.0), 0.3, log(0.5))
//' Sigma <- choicer:::build_var_mat(L_params, K_w = 2, rc_correlation = TRUE)
//' Sigma  # 2x2 covariance matrix
//' @keywords internal
// [[Rcpp::export]]
arma::mat build_var_mat(const arma::vec &L_params, const int K_w,
                        const bool rc_correlation) {
  arma::mat L = build_L_mat(L_params, K_w, rc_correlation);
  // Return the variance matrix
  return L * L.t();
}

// ============================================================================
// Likelihood units: the panel structure shared by the estimation kernels
//
// A likelihood unit is one decision maker. Unit u owns the contiguous choice
// situations t = off[u], ..., off[u+1]-1 (situations are sorted by decision
// maker) and ONE K_w x S block of draws eta_u, shared by all of them. With
// Ti = NULL every situation is its own unit (off = 0, 1, ..., N): the
// cross-sectional likelihood, with the same draw block and weight per
// situation as before. With beta_us = mu_final + Gamma_us the tastes at draw
// s and j_t the alternative chosen in situation t:
//
//   lambda_us = sum_t log P_ts(j_t)               log-probability of u's choices
//   ell_u     = log (1/S) sum_s exp(lambda_us)    simulated log-likelihood
//   omega_us  = exp(lambda_us - LSE_s lambda_u.)  posterior draw weights (sum 1)
//   s_u       = sum_s omega_us sum_t g_uts        score, where
//   g_uts     = sum_a (1{a = j_t} - P_ts(a)) z_tsa,  z_tsa = dV_tsa / dtheta
//
// Everything stays in log space: each log P_ts(j_t) is the chosen stabilized
// utility minus the log-sum-exp denominator (never the log of a probability
// that may have underflowed), and the draw weights come from a max-shifted
// LSE of lambda, so neither one choice probability nor a product of many
// underflows. lambda_us is finite whenever the utilities are.
// ============================================================================

// ============================================================================
// Layout of the stacked design
//
// Offsets into the stacked design are 64-bit (mxl_off). R caps a matrix at
// 2^31 - 1 rows, so they keep offset arithmetic clear of int overflow; element
// counts and offsets are Armadillo's (arma::uword), 64-bit because choicer
// defines ARMA_64BIT_WORD (src/Makevars), so a view of a design past 2^32 - 1
// elements does not wrap. Armadillo indices are taken from the offsets where
// a matrix is sliced.
// ============================================================================
using mxl_off = choicer_off;

// The theta-independent layout of the stacked design, built and validated on
// the primary thread by every kernel call: the situation core shared with the
// MNL and NL kernels (ChoiceLayout, choicer_internal.h: situation offsets,
// and pointers to the per-row alternative codes and per-situation choices,
// read in place from the kernel's integer arguments), plus the likelihood
// units of a panel. It is valid for one kernel call only: when an argument
// was coerced from doubles, the integer copy belongs to the Rcpp wrapper of
// that call. Building it is one pass over the rows (in parallel above 10^6
// rows) and one over the situations, about 20 ms at 10^8 rows: under 1% of
// an evaluation, too little to repay caching.
struct MxlLayout : ChoiceLayout {
  mxl_off U = 0;                 // likelihood units
  mxl_off max_unit_rows = 0;     // rows of the largest unit
  bool include_outside_option = false;
  std::vector<mxl_off> unit_off; // panel: unit u's situations [unit_off[u],
                                 // unit_off[u+1]); empty in the cross-section

  // First situation of unit u (u = U gives N): unit u is situation u in the
  // cross-section.
  mxl_off unit_first(const mxl_off u) const {
    return unit_off.empty() ? u : unit_off[u];
  }
  // The slot of situation t's choice in P
  int chosen(const mxl_off t) const {
    return include_outside_option ? choice[t] : choice[t] - 1;
  }
};

// Build and validate the layout: the units (Ti), the situations and the
// alternative codes (ChoiceLayout's checks), then the choices; the delta,
// draw and W checks follow in MxlUnitData. Situations are sorted by decision
// maker; Ti = NULL makes every situation its own unit (the cross-section).
// The Ti checks mirror hmnl_gibbs (src/hmnlogit.cpp): at least one
// respondent, every Ti positive, and sum(Ti) equal to the number of choice
// situations.
inline MxlLayout mxl_layout_build(const arma::mat& X,
                                  const Rcpp::IntegerVector& alt_idx,
                                  const Rcpp::IntegerVector& choice_idx,
                                  const Rcpp::IntegerVector& M,
                                  const arma::vec* weights,
                                  const Rcpp::Nullable<Rcpp::IntegerVector>& Ti,
                                  const bool include_outside_option,
                                  const char* kernel) {
  MxlLayout lay;
  const mxl_off N = M.size();
  lay.include_outside_option = include_outside_option;

  // Likelihood units: decision makers (Ti) or choice situations (Ti = NULL).
  if (Ti.isNull()) {
    lay.U = N;
  } else {
    const Rcpp::IntegerVector T_u(Ti.get());
    if (T_u.size() == 0) {
      Rcpp::stop("Ti must contain at least one respondent.");
    }
    lay.U = T_u.size();
    lay.unit_off.assign(lay.U + 1, 0);
    for (mxl_off u = 0; u < lay.U; ++u) {
      if (T_u[u] == NA_INTEGER || T_u[u] < 1) {
        Rcpp::stop("Ti must be positive for every respondent (Ti[%d] = %s).",
                   u + 1,
                   T_u[u] == NA_INTEGER ? "NA" : std::to_string(T_u[u]));
      }
      lay.unit_off[u + 1] = lay.unit_off[u] + T_u[u];
    }
    if (lay.unit_off[lay.U] != N) {
      Rcpp::stop("sum(Ti) (%d) does not match the number of choice situations "
                 "(%d).", lay.unit_off[lay.U], N);
    }
  }

  // Situations and alternative codes (ChoiceLayout's checks, shared with the
  // MNL and NL kernels), then the choices.
  choice_layout_situations(lay, X, alt_idx, M, weights, &choice_idx);
  choice_layout_codes(lay);
  validate_choices(lay, include_outside_option, kernel);

  lay.max_unit_rows = 0;
  for (mxl_off u = 0; u < lay.U; ++u) {
    lay.max_unit_rows =
        std::max(lay.max_unit_rows, lay.row_off[lay.unit_first(u + 1)] -
                                        lay.row_off[lay.unit_first(u)]);
  }
  return lay;
}

// Draw batches. A unit's R x S matrices (WGamma, and DiffW for the score)
// are formed B draws at a time, B = S unless the kernel's n_mats R x B
// matrices would exceed MXL_BATCH_BYTES: a unit is split only when R S
// exceeds 2^18 (2^19 when n_mats = 1), i.e. a long panel, a choice set of
// hundreds of alternatives or a very large S, and a thread's scratch then
// grows with R rather than R S. The S-length and K_w x S pieces (lambda, the
// draw weights, Gamma, BW) stay whole. The budget is small because the draw
// loop reads each batch a few times in quick succession, and a batch that
// stays in cache is fastest. A split unit's draws are partitioned into
// batches whose sizes differ by at most one. draw_batch > 0 caps the batch
// size instead (tests).
constexpr double MXL_BATCH_BYTES = 4.0 * 1024.0 * 1024.0;

// Number of draw batches of a unit of R rows, for a kernel that keeps n_mats
// R x B matrices (1: the unit is not split).
inline int mxl_batch_count(const mxl_off R, const int S, const int n_mats,
                           const int draw_batch) {
  const double b_max =
      draw_batch > 0
          ? draw_batch
          : MXL_BATCH_BYTES / (sizeof(double) * n_mats *
                               static_cast<double>(std::max<mxl_off>(R, 1)));
  if (b_max >= S) return 1;
  const int b = std::max(1, static_cast<int>(b_max));
  return (S + b - 1) / b;
}

// Everything the per-unit routines read, shared by all threads: the stacked
// design and its layout, the parameters at theta, the draw source and the
// model flags. Construct it once per kernel call, on the primary thread: it
// owns what it derives (the layout, the parsed parameters and the Halton
// generator) and validates the inputs, each check O(1) or one pass over the
// rows or the situations. Holds no SEXP.
struct MxlUnitData {
  const arma::mat& X;
  const arma::mat& W;              // row-aligned with X, or J x K_w
  const arma::cube& eta_draws;     // store mode: slice u
  const arma::uvec& rc_dist;
  const MxlParams par;
  const MxlLayout lay;
  HaltonGen gen;                   // generate mode: Halton block u + 1
  const int n_params;
  const int S;                     // draws per unit
  const int draw_batch;            // > 0: draws per batch (tests); 0: automatic
  const bool use_generate, rc_correlation, rc_mean, use_asc,
      include_outside_option, alt_level_W;

  MxlUnitData(const arma::vec& theta, const arma::mat& X_, const arma::mat& W_,
              const Rcpp::IntegerVector& alt_idx,
              const Rcpp::IntegerVector& choice_idx,
              const Rcpp::IntegerVector& M, const arma::vec* weights,
              const arma::cube& eta_draws_, const arma::uvec& rc_dist_,
              const bool rc_correlation_, const bool rc_mean_,
              const bool use_asc_, const bool include_outside_option_,
              const int gen_seed, const int gen_scramble, const int gen_S,
              const Rcpp::Nullable<Rcpp::IntegerVector>& Ti,
              const int draw_batch_, const char* kernel)
      : X(X_), W(W_), eta_draws(eta_draws_), rc_dist(rc_dist_),
        // Parse theta into parameter blocks (shared helper; validates theta
        // and the length of rc_dist)
        par(parse_mxl_theta(theta, X_.n_cols, W_.n_cols, rc_dist_,
                            rc_correlation_, rc_mean_, use_asc_,
                            include_outside_option_)),
        lay(mxl_layout_build(X_, alt_idx, choice_idx, M, weights, Ti,
                             include_outside_option_, kernel)),
        n_params(theta.n_elem),
        S(gen_seed >= 0 ? gen_S : static_cast<int>(eta_draws_.n_cols)),
        draw_batch(draw_batch_),
        use_generate(gen_seed >= 0), rc_correlation(rc_correlation_),
        rc_mean(rc_mean_), use_asc(use_asc_),
        include_outside_option(include_outside_option_),
        alt_level_W(W_.n_rows != X_.n_rows) {
    const int K_w = W.n_cols;
    if (use_asc && lay.n_rows > 0 &&
        static_cast<int>(par.delta.n_elem) < lay.J) {
      Rcpp::stop("Theta's delta (ASC) block implies %d alternatives but "
                 "alt_idx references alternative %d.", par.delta.n_elem,
                 lay.J);
    }
    if (use_generate) {
      if (gen_S <= 0) Rcpp::stop("gen_S must be positive when gen_seed >= 0");
      if (K_w > HALTON_N_PRIMES) {
        Rcpp::stop("K_w exceeds the primes table size (128); reduce K_w or "
                   "extend the primes table.");
      }
      // The store-mode points (2) are the prediction kernels' alone: here they
      // would be the identity points through HaltonGen's own inverse CDF.
      if (gen_scramble != 0 && gen_scramble != 1) {
        Rcpp::stop("gen_scramble must be 0 or 1 when gen_seed >= 0; got %d.",
                   gen_scramble);
      }
    } else {
      // eta_draws holds one K_w x S draw block per likelihood unit.
      if (Ti.isNull() && static_cast<mxl_off>(eta_draws.n_slices) != lay.N) {
        Rcpp::stop("eta_draws 3rd dimension (%d) does not match N (%d)",
                   eta_draws.n_slices, lay.N);
      }
      if (Ti.isNotNull() && static_cast<mxl_off>(eta_draws.n_slices) != lay.U) {
        Rcpp::stop("eta_draws 3rd dimension (%d) does not match the number "
                   "of decision makers (%d)", eta_draws.n_slices, lay.U);
      }
      if (static_cast<int>(eta_draws.n_rows) != K_w) {
        Rcpp::stop("eta_draws 1st dimension (%d) does not match K_w (%d)",
                   eta_draws.n_rows, K_w);
      }
    }
    if (alt_level_W && lay.n_rows > 0 &&
        static_cast<mxl_off>(W.n_rows) < lay.J) {
      Rcpp::stop("W must be row-aligned with X (%d rows) or contain one row "
                 "per global alternative (at least %d rows); got %d rows.",
                 X.n_rows, lay.J, W.n_rows);
    }
    if (weights && !lay.unit_off.empty()) {
      // Weights are decision-maker weights (objective sum_u w_u ell_u):
      // constant within each unit, whose weight is weights[unit_first(u)].
      // NaN-aware: two NaN weights count as equal (finiteness is checked in
      // R).
      for (mxl_off u = 0; u < lay.U; ++u) {
        const double w_u = (*weights)[lay.unit_off[u]];
        for (mxl_off t = lay.unit_off[u] + 1; t < lay.unit_off[u + 1]; ++t) {
          const double w_t = (*weights)[t];
          if (w_t != w_u && !(std::isnan(w_t) && std::isnan(w_u))) {
            Rcpp::stop("weights must be constant within each decision maker "
                       "(unit %d).", u + 1);
          }
        }
      }
    }
    // Beyond its R x B draw batches, a unit's per-thread scratch is O(R):
    // X_u, W_u, the base utilities, d_bar and a batch of at least one draw,
    // about 8 R (K_x + K_w + 4) bytes. An allocation failure inside the
    // parallel region would terminate R, so refuse here when that would
    // exceed 2 GiB, i.e. one decision maker with tens of millions of rows.
    const double bytes = sizeof(double) *
                         static_cast<double>(lay.max_unit_rows) *
                         (X.n_cols + K_w + 4.0);
    if (bytes > 2.0 * 1024 * 1024 * 1024) {
      Rcpp::stop("The largest decision maker stacks %d alternative rows, "
                 "which needs %.1f GB of scratch memory per thread. Check "
                 "person_col: it should identify decision makers, not "
                 "markets.", lay.max_unit_rows, bytes / 1e9);
    }
    if (use_generate) {
      // Constructed here, outside the parallel region; read-only after.
      gen = HaltonGen(static_cast<uint64_t>(gen_seed), S, K_w, gen_scramble);
    }
  }
  MxlUnitData(const MxlUnitData&) = delete;
  MxlUnitData& operator=(const MxlUnitData&) = delete;
};

// Thread-private state of the unit in hand; declare it inside the parallel
// region. Members are resized per unit and keep their memory across units,
// so a unit allocates nothing once the thread has seen a larger one (except
// that Armadillo keeps a matrix of at most 16 elements in local storage, so a
// member that shrinks that far re-acquires its heap block when it grows
// again). Not copyable: `eta` may point into this object's own eta_buf.
struct MxlUnitScratch {
  MxlUnitScratch() = default;
  MxlUnitScratch(const MxlUnitScratch&) = delete;
  MxlUnitScratch& operator=(const MxlUnitScratch&) = delete;
  mxl_off t0 = 0, t1 = 0; // situations [t0, t1)
  mxl_off r0 = 0, R = 0;  // stacked rows [r0, r0 + R)
  int nb = 1, s0 = 0;   // draw batches of the unit; first draw in WGamma
  mxl_off wg_r0 = 0;    // unit row of WGamma's first row (0, or a situation's)
  const double* eta = nullptr; // K_w x S draws: eta_buf, or cube slice u
  arma::mat eta_buf;    // generate mode: the unit's Halton block
  arma::mat X_u;        // R x K_x
  arma::mat W_u;        // R x K_w
  arma::mat W_t;        // m_t x K_w: one situation's rows (Hessian, long units)
  arma::vec bu;         // R: base utilities X_u beta + W_u mu_final + delta
  arma::mat Gamma;      // K_w x S random-coefficient draws (Gamma_final)
  arma::mat Dgamma1;    // K_w x S first derivative of the RC transform
  arma::mat Dgamma2;    // K_w x S second derivative (Hessian only)
  arma::mat WGamma;     // R x B: W_u * Gamma, draws [s0, s0 + B); or m_t x S
  arma::vec V, P;       // one situation at one draw: leading m + o entries
  arma::vec lambda;     // S: sum_t log P_ts(j_t)
  arma::vec omega;      // S: posterior draw weights
  arma::vec w_b;        // B: a batch's draw weights, relative to ref
  arma::mat DiffW;      // R x B unweighted residuals 1{r = j_t} - P_ts(r)
  arma::vec d_bar;      // R: DiffW * omega, over the batches
  arma::mat BW, A;      // K_w x S and K_w x K_w: Cholesky-block collapse
  arma::rowvec omega_t; // S: omega as a row, for the Cholesky collapse
};

// Non-owning K_w x S view of the loaded unit's draws (no copy of the cube
// slice in store mode). Bind it only to a const matrix and never assign to
// it: in store mode it aliases the caller's R array.
inline const arma::mat mxl_eta_view(const MxlUnitData& ud,
                                    const MxlUnitScratch& sc) {
  return arma::mat(const_cast<double*>(sc.eta), ud.W.n_cols, ud.S, false,
                   true);
}

// Load unit u into the thread's buffers: its rows of X and W (for an
// alternative-level W, the rows of its alternatives), the base utilities
// X_u beta + W_u mu_final + delta of those rows (per unit, inside the
// parallel region: nothing of the stacked length is formed per evaluation),
// the draws eta_u (cube slice u in place, or on-the-fly Halton block u + 1),
// Gamma_u = L eta_u with the derivatives of its transform (Dgamma1 unless
// only the draws themselves are needed), and the unit's draw batch size for
// a kernel that keeps n_mats R x B matrices.
inline void mxl_unit_load(const MxlUnitData& ud, const mxl_off u,
                          MxlUnitScratch& sc, const int n_mats,
                          const bool with_Dgamma1 = true,
                          const bool with_Dgamma2 = false) {
  const MxlLayout& lay = ud.lay;
  sc.t0 = lay.unit_first(u);
  sc.t1 = lay.unit_first(u + 1);
  sc.r0 = lay.row_off[sc.t0];
  sc.R = lay.row_off[sc.t1] - sc.r0;
  const arma::uword r0 = static_cast<arma::uword>(sc.r0);
  const arma::uword r1 = static_cast<arma::uword>(sc.r0 + sc.R - 1);
  const int K_w = ud.W.n_cols;
  sc.X_u = ud.X.rows(r0, r1);
  if (!ud.alt_level_W) {
    sc.W_u = ud.W.rows(r0, r1);
  } else {
    sc.W_u.set_size(sc.R, K_w);
    for (int k = 0; k < K_w; ++k) {
      for (mxl_off i = 0; i < sc.R; ++i) {
        sc.W_u(i, k) = ud.W(lay.alt0(sc.r0 + i), k);
      }
    }
  }
  sc.bu = sc.X_u * ud.par.beta;
  if (ud.rc_mean) sc.bu += sc.W_u * ud.par.mu_final;
  if (ud.use_asc) {
    for (mxl_off i = 0; i < sc.R; ++i) {
      sc.bu[i] += ud.par.delta[lay.alt0(sc.r0 + i)];
    }
  }
  if (ud.use_generate) {
    ud.gen.fill_eta_i(sc.eta_buf, static_cast<uint64_t>(u) + 1);
    sc.eta = sc.eta_buf.memptr();
  } else {
    sc.eta = ud.eta_draws.slice_memptr(static_cast<arma::uword>(u));
  }
  batch_gamma_draws_into(sc.Gamma, ud.par.L, mxl_eta_view(ud, sc), ud.rc_dist,
                         with_Dgamma1 ? &sc.Dgamma1 : nullptr,
                         with_Dgamma2 ? &sc.Dgamma2 : nullptr);
  sc.nb = mxl_batch_count(sc.R, ud.S, n_mats, ud.draw_batch);
}

// WGamma for the draws [s0, s1) of the loaded unit: W_u Gamma(:, s0:s1-1) in
// a single dgemm.
inline void mxl_batch_wgamma(MxlUnitScratch& sc, const int s0, const int s1) {
  sc.s0 = s0;
  sc.wg_r0 = 0;
  if (s0 == 0 && s1 == static_cast<int>(sc.Gamma.n_cols)) {
    sc.WGamma = sc.W_u * sc.Gamma;
  } else {
    sc.WGamma = sc.W_u * sc.Gamma.cols(s0, s1 - 1);
  }
}

// WGamma for one situation of the loaded unit, its m rows from unit row r,
// at every draw: W_u(r:r+m-1, :) Gamma in a single dgemm.
inline void mxl_situation_wgamma(MxlUnitScratch& sc, const mxl_off r,
                                 const int m) {
  sc.s0 = 0;
  sc.wg_r0 = r;
  sc.W_t = sc.W_u.rows(static_cast<arma::uword>(r),
                       static_cast<arma::uword>(r + m - 1));
  sc.WGamma = sc.W_t * sc.Gamma;
}

// Logit probabilities of situation t at draw s, in the leading m + o entries
// of sc.P: inside utilities bu + WGamma(rows of t, s), with the outside
// option (V = 0) in slot 0 when present. Leaves the max-shifted utilities in
// sc.V and returns their log-sum-exp, so that log P_ts(a) = sc.V(a) -
// log_denom: exact wherever the utilities are finite, while P(a) itself
// turns subnormal beyond a utility gap of about 708 and zero beyond about
// 745. V and P only grow: sized per situation, they would bounce between the
// heap and Armadillo's 16-element local storage whenever consecutive
// situations straddle that size, which the draw-major loop of
// mxl_unit_simulate() makes happen at every draw. WGamma must hold
// situation t at draw s: the batch containing s (mxl_batch_wgamma()), or
// t's own rows at every draw (mxl_situation_wgamma()).
inline double mxl_situation_probs(const MxlUnitData& ud, MxlUnitScratch& sc,
                                  const mxl_off t, const int s) {
  const mxl_off first = ud.lay.row_off[t] - sc.r0;      // first row in unit
  const int m = static_cast<int>(ud.lay.row_off[t + 1] - ud.lay.row_off[t]);
  const int o = ud.include_outside_option ? 1 : 0;  // outside option: slot 0
  const int n = m + o;
  if (static_cast<int>(sc.V.n_elem) < n) {
    sc.V.set_size(n);
    sc.P.set_size(n);
  }
  double* v = sc.V.memptr();
  if (o) v[0] = 0.0;
  const double* bu = sc.bu.memptr() + first;
  const double* wg = sc.WGamma.colptr(s - sc.s0) + (first - sc.wg_r0);
  for (int a = 0; a < m; ++a) v[o + a] = bu[a] + wg[a];
  return stable_softmax_n(v, sc.P.memptr(), n);
}

// Draw loop of the loaded unit, in sc.nb batches of draws: lambda_s =
// sum_t log P_ts(j_t), each term the chosen shifted utility minus the
// log-sum-exp from mxl_situation_probs() (the multinomial logit kernels'
// V_choice - log_denom), never the log of a probability. With `score`, each
// batch's UNWEIGHTED residuals DiffW (R x B; over the rows of situation t,
// column s holds 1{r = j_t} - P_ts(r) for the inside alternatives, the
// outside option carrying no parameter) are folded into the score's pieces
// before the next batch overwrites them:
//   BW(:, batch) = W_u' DiffW                       (unweighted, K_w x B)
//   d_bar = sum_s omega_s DiffW(:, s): for a unit in one batch, from omega
//     after the loop, exactly as without batches; across batches, by a
//     streaming max shift: the batch's weights exp(lambda_s - ref) are
//     relative to the largest lambda seen so far. d_bar is rescaled by
//     exp(ref_old - ref) as ref rises, then divided by the same sum of
//     shifted weights used to normalize omega. Keeping the shift separate
//     from that sum avoids subtracting a rounded, large-magnitude lse.
// Fills normalized posterior weights omega and returns lse = log sum_s
// exp(lambda_s).
inline double mxl_unit_simulate(const MxlUnitData& ud, MxlUnitScratch& sc,
                                const bool score) {
  const MxlLayout& lay = ud.lay;
  const int S = ud.S;
  const int K_w = ud.W.n_cols;
  const int o = ud.include_outside_option ? 1 : 0; // slot of 1st inside alt
  const bool stream = score && sc.nb > 1;          // fold d_bar batch by batch
  sc.lambda.zeros(S);
  if (score) sc.BW.set_size(K_w, S);
  double ref = -arma::datum::inf; // running maximum of lambda (stream)
  for (int k = 0; k < sc.nb; ++k) {
    // Batch k: draws [s0, s1), sizes differing by at most one.
    const int s0 = static_cast<int>(static_cast<long long>(k) * S / sc.nb);
    const int s1 = static_cast<int>(static_cast<long long>(k + 1) * S / sc.nb);
    mxl_batch_wgamma(sc, s0, s1);
    if (score) sc.DiffW.set_size(sc.R, s1 - s0);
    // Draws outer, situations inner: each column of the batch is swept in
    // row order (a situation-major sweep jumps R rows between draws of a
    // long panel), and lambda_s still sums the situations in order.
    for (int s = s0; s < s1; ++s) {
      for (mxl_off t = sc.t0; t < sc.t1; ++t) {
        const arma::uword r =
            static_cast<arma::uword>(lay.row_off[t] - sc.r0); // first row of t
        const int m = static_cast<int>(lay.row_off[t + 1] - lay.row_off[t]);
        const int chosen_alt = lay.chosen(t);  // slot in P, validated
        const double log_denom = mxl_situation_probs(ud, sc, t, s);
        sc.lambda(s) += sc.V(chosen_alt) - log_denom;
        if (score) {
          sc.DiffW.col(s - s0).subvec(r, r + m - 1) =
              -sc.P.subvec(o, o + m - 1);
          if (chosen_alt >= o) sc.DiffW(r + chosen_alt - o, s - s0) += 1.0;
        }
      }
    }
    if (!score) continue;
    if (K_w > 0) { // BW(:, s0:s1-1) = W_u' DiffW: one dgemm, written in place
      arma::mat BW_b(sc.BW.colptr(s0), K_w, s1 - s0, false, true);
      BW_b = sc.W_u.t() * sc.DiffW;
    }
    if (!stream) continue;
    // Fold the batch into d_bar, relative to the running maximum.
    const double* lambda_b = sc.lambda.memptr() + s0;
    const double ref_new = std::max(ref, direct_max_n(lambda_b, s1 - s0));
    if (ref_new == -arma::datum::inf) {  // no draw with positive weight yet
      sc.d_bar.zeros(sc.R);
    } else {
      sc.w_b.set_size(s1 - s0);
      for (int j = 0; j < s1 - s0; ++j) sc.w_b[j] = std::exp(lambda_b[j] - ref_new);
      if (ref == -arma::datum::inf) {    // first batch with positive weight
        sc.d_bar = sc.DiffW * sc.w_b;
      } else {
        sc.d_bar *= std::exp(ref - ref_new);
        sc.d_bar += sc.DiffW * sc.w_b;
      }
    }
    ref = ref_new;
  }
  const double lse = log_sum_exp_n(sc.lambda.memptr(), S); // = logSumExp()
  const double lambda_max = direct_max_n(sc.lambda.memptr(), S);
  sc.omega = arma::exp(sc.lambda - lambda_max);
  const double weight_sum = arma::accu(sc.omega);
  sc.omega /= weight_sum;
  if (score) {
    if (!stream) {                       // one batch: the draw weights directly
      sc.d_bar = sc.DiffW * sc.omega;    // R x 1, one dgemv
    } else if (!std::isfinite(lse)) {    // no score, as omega has none
      sc.d_bar.fill(arma::datum::nan);
    } else {                            // ref is now the global lambda_max
      sc.d_bar /= weight_sum;
    }
  }
  return lse;
}

// Score of the loaded unit, s_u = sum_s omega_s sum_t g_uts, by the BLAS-3
// collapse of the residuals, from the pieces mxl_unit_simulate(score = true)
// folded batch by batch. Utilities are linear in beta, mu and the ASC
// dummies, so those blocks need only d_bar = DiffW omega; the Cholesky block
// couples the residuals of draw s with eta_s:
//   beta : X_u' d_bar
//   mu   : (W_u' d_bar) % dmu_final_dmu                          (rc_mean)
//   L    : A = ((W_u' DiffW % Dgamma1) diag(omega)) eta_u',  s[L_pq] = dL_pq A(p, q)
//   delta: d_bar scattered by alternative over the stacked rows
inline void mxl_unit_score(const MxlUnitData& ud, MxlUnitScratch& sc,
                           arma::vec& score) {
  const MxlParams& par = ud.par;
  const int K_w = ud.W.n_cols;
  score.zeros(ud.n_params);

  // Beta block (X_u loaded by mxl_unit_load)
  score.subvec(par.idx_beta_start, par.idx_mu_start - 1) = sc.X_u.t() * sc.d_bar;

  if (K_w > 0) {
    // Mu block (only if rc_mean)
    if (ud.rc_mean) {
      score.subvec(par.idx_mu_start, par.idx_L_start - 1) =
          (sc.W_u.t() * sc.d_bar) % par.dmu_final_dmu;
    }

    // L block: the only block with per-draw eta coupling; BW = W_u' DiffW
    sc.BW %= sc.Dgamma1;               // Dgamma1 = 1 for normal rows
    sc.omega_t = sc.omega.t();
    sc.BW.each_row() %= sc.omega_t;    // draw weights
    sc.A = sc.BW * mxl_eta_view(ud, sc).t(); // K_w x K_w, one dgemm
    if (ud.rc_correlation) {
      int lp = 0;
      for (int p = 0; p < K_w; ++p) {
        for (int q = 0; q <= p; ++q, ++lp) {
          const double dLpq = (p == q) ? par.L(p, p) : 1.0;
          score[par.idx_L_start + lp] = dLpq * sc.A(p, q);
        }
      }
    } else {
      // Diagonal L: s[L_p] = L(p,p) * A(p,p)
      for (int p = 0; p < K_w; ++p) {
        score[par.idx_L_start + p] = par.L(p, p) * sc.A(p, p);
      }
    }
  }

  // Delta block (scatter -- irregular alt-index mapping)
  if (ud.use_asc) {
    scatter_delta_grad(score, par.idx_delta_start, sc.d_bar,
                       AltCodes0{ud.lay.alt + sc.r0}, static_cast<int>(sc.R),
                       ud.include_outside_option, 1.0);
  }
}

//' Log-likelihood and gradient for Mixed Logit
//'
//' Computes the log-likelihood and its gradient for the Mixed Logit model using
//' OpenMP for parallelization. Allows for inclusion of alternative-specific
//' constants, outside option, observation weights, correlated random
//' coefficients, and panel data (one draw block per decision maker, see Ti).
//'
//' @param theta vector collecting model parameters (beta, mu, L, delta (ASCs))
//' @param X design matrix for covariates with fixed coefficients; sum(M_i) x K_x
//' @param W design matrix for covariates with random coefficients; sum(M_i) x K_w or J x K_w
//' @param alt_idx sum(M) x 1 vector with indices of alternatives within each choice set; 1-based indexing
//' @param choice_idx N x 1 vector with indices of chosen alternatives; 1-based indexing relative to X; 0 is used if include_outside_option=True
//' @param M N x 1 vector with number of alternatives for each individual
//' @param weights N x 1 vector with weights for each observation; when Ti is
//'   supplied they must be constant within each decision maker
//' @param eta_draws Array of standard-normal draws, K_w x S x U, where U is the
//'   number of decision makers when Ti is supplied and the number of choice
//'   situations otherwise
//' @param rc_dist K_w x 1 integer vector indicating distribution of random coefficients: 0 = normal, 1 = log-normal
//' @param rc_correlation whether random coefficients should be correlated
//' @param rc_mean whether to estimate means for random coefficients. If so, mean parameters (mu) should be included in theta after beta parameters.
//' @param use_asc whether to use alternative-specific constants. If so, parameters should be included in theta after beta and L (and mu, if applicable).
//' @param include_outside_option whether to include outside option normalized to 0 (if so, the outside option is not included in the data)
//' @param gen_seed Integer master seed for the on-the-fly Halton generator. \code{< 0}
//'   (default) uses the materialized \code{eta_draws} cube; \code{>= 0} generates draws
//'   on the fly from this seed.
//' @param gen_scramble Integer scramble mode for on-the-fly generation: \code{0} =
//'   identity permutations (plain Halton, compat), \code{1} = seeded position-wise
//'   digit permutations; other values are an error.
//' @param gen_S Integer number of draws per individual, used only when \code{gen_seed >= 0}.
//' @param Ti Optional integer vector with the number of choice situations of
//'   each decision maker (panel likelihood); situations must be sorted by
//'   decision maker. NULL (default): every choice situation is its own unit
//'   (cross-sectional likelihood).
//' @param draw_batch Integer; \code{0} (default) forms each decision maker's
//'   draws in batches sized to a per-thread memory budget, a positive value
//'   caps the number of draws per batch (for tests).
//' @returns List with the negated log-likelihood (\code{objective}), its
//'   \code{gradient}, and an \code{overflow} flag indicating that a
//'   non-finite objective was replaced by the finite optimizer sentinel.
//' @note For log-normal random coefficients (rc_dist=1) with rc_mean=TRUE,
//'   the distribution is a shifted log-normal: beta_k = exp(mu_k) + exp(L_k * eta),
//'   where exp(mu_k) shifts the location and exp(L_k * eta) ~ LogNormal(0, sigma_k^2).
//'   This differs from the textbook parameterization exp(mu_k + L_k * eta).
//' @examples
//' \donttest{
//' library(data.table)
//' set.seed(42)
//' N <- 50; J <- 3
//' dt <- data.table(id = rep(1:N, each = J), alt = rep(1:J, N))
//' dt[, `:=`(x1 = rnorm(.N), w1 = rnorm(.N))]
//' dt[, choice := 0L]
//' dt[, choice := sample(c(1L, rep(0L, J - 1))), by = id]
//' d <- prepare_mxl_data(dt, "id", "alt", "choice", "x1", "w1")
//' eta <- get_halton_normals(50, d$N, ncol(d$W))
//' K_x <- ncol(d$X); K_w <- ncol(d$W); J <- nrow(d$alt_mapping)
//' theta <- rep(0, K_x + K_w + J - 1)
//' result <- choicer:::mxl_loglik_gradient_parallel(theta, d$X, d$W, d$alt_idx,
//'   d$choice_idx, d$M, d$weights, eta, rc_dist = rep(0L, K_w),
//'   rc_correlation = FALSE, rc_mean = FALSE)
//' result$objective
//' }
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List mxl_loglik_gradient_parallel(
    const arma::vec &theta, const arma::mat &X, const arma::mat &W,
    const Rcpp::IntegerVector &alt_idx,
    const Rcpp::IntegerVector &choice_idx,
    const Rcpp::IntegerVector &M, const arma::vec &weights,
    const arma::cube &eta_draws, const arma::uvec &rc_dist,
    const bool rc_correlation = true, const bool rc_mean = false,
    const bool use_asc = true, const bool include_outside_option = false,
    const int gen_seed = -1, const int gen_scramble = 1, const int gen_S = 0,
    const Rcpp::Nullable<Rcpp::IntegerVector> Ti = R_NilValue,
    const int draw_batch = 0) {
  // Inputs, layout and parameters at theta, validated on the primary thread
  const MxlUnitData ud(theta, X, W, alt_idx, choice_idx, M,
                       &weights, eta_draws, rc_dist,
                       rc_correlation, rc_mean, use_asc, include_outside_option,
                       gen_seed, gen_scramble, gen_S, Ti, draw_batch,
                       "mxl_loglik_gradient_parallel");
  const MxlLayout& lay = ud.lay;
  const int n_params = ud.n_params;
  const double log_S = std::log(static_cast<double>(ud.S));

  // Prepare global accumulators
  double global_loglik = 0.0;
  arma::vec global_grad = arma::zeros(n_params);

#ifdef _OPENMP
#pragma omp parallel
#endif
  {
    // Thread-local accumulators and unit scratch
    double local_loglik = 0.0;
    arma::vec local_grad = arma::zeros(n_params);
    MxlUnitScratch sc;
    arma::vec s_u; // score of unit u

// Loop over likelihood units in parallel
#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
    for (mxl_off u = 0; u < lay.U; ++u) {
      mxl_unit_load(ud, u, sc, 2);
      const double lse = mxl_unit_simulate(ud, sc, true); // log sum_s exp(lambda_s)
      mxl_unit_score(ud, sc, s_u);

      // ell_u = lse - log(S); unit weight = weight of its first situation
      const double w_u = weights[lay.unit_first(u)];
      local_loglik += w_u * (lse - log_S);
      local_grad += w_u * s_u;
    } // end unit loop

// Combine thread-local results
#ifdef _OPENMP
#pragma omp critical
#endif
    {
      global_loglik += local_loglik;
      global_grad += local_grad;
    }
  } // end parallel region

  // Sanitize NaN/Inf so the optimizer's line-search can backtrack instead
  // of stalling at an undefined objective. The minimizer expects a finite
  // reference value; a "very bad" sentinel lets it shrink the step. Since
  // every log probability is formed in log space, the objective is finite
  // whenever the utilities are, however improbable the observed choices:
  // the sentinel marks utilities that overflowed (e.g. an exploding
  // Cholesky diagonal during a line search), not a poor fit. Finite
  // objectives are then unbounded, so run_mxlogit() lifts the sentinel above
  // every objective seen along the optimizer's path.
  double obj = -global_loglik;
  arma::vec grad = -global_grad;
  const bool overflow = !std::isfinite(obj);
  if (overflow) {
    obj = 1e10;
    grad.zeros();
  } else {
    grad.elem(arma::find_nonfinite(grad)).zeros();
  }
  return Rcpp::List::create(Rcpp::Named("objective") = obj,
                            Rcpp::Named("gradient") = grad,
                            Rcpp::Named("overflow") = overflow);
}

// vech(): lower-triangular vectorisation (including the diagonal), row-major
// Matches the row-major packing of L_params used by build_L_mat() above and
// by the outer parameter loop in jacobian_vech_Sigma() below.
inline arma::vec vech(const arma::mat &M) {
  arma::uword K = M.n_rows;
  arma::vec out(K * (K + 1) / 2);
  arma::uword idx = 0;
  for (arma::uword i = 0; i < K; ++i)
    for (arma::uword j = 0; j <= i; ++j)
      out(idx++) = M(i, j);
  return out;
}

//' Utility to compute analytical Jacobian of random coefficient matrix transformed by vech (dVech(Sigma) / dTheta)
//'
//' @param L_params flattened choleski decomposition version of the random coefficient parameters matrix
//' @param K_w dimension of the random coefficient parameter (symmetric) matrix
//' @param rc_correlation whether random coefficients are correlated
//' @returns Jacobian (dVech(Sigma) / dTheta)
//' @examples
//' L_params <- c(log(0.8), 0.2, log(0.6))
//' J_mat <- choicer:::jacobian_vech_Sigma(L_params, K_w = 2, rc_correlation = TRUE)
//' dim(J_mat)  # 3 x 3 for K_w=2 correlated
//' @keywords internal
// [[Rcpp::export]]
arma::mat jacobian_vech_Sigma(const arma::vec &L_params, const int K_w,
                              const bool rc_correlation = true) {
  // dimensions
  const int L_size = rc_correlation ? K_w * (K_w + 1) / 2 : K_w;

  arma::mat L = build_L_mat(L_params, K_w, rc_correlation);
  arma::mat J(L_size, L_size, arma::fill::zeros);

  // loop over parameters
  arma::mat E(K_w, K_w, arma::fill::zeros); // holds dL / dtheta_m
  std::size_t idx_param = 0;

  if (rc_correlation) {
    for (int i = 0; i < K_w; ++i) {
      for (int j = 0; j <= i; ++j, ++idx_param) {
        // reset E
        E.zeros();
        if (i == j) {        // diagonal: L_ii = exp(z_i)
          E(i, j) = L(i, i); // dL_ii / dz_i = exp(z_i)
        } else {             // off-diagonal parameter
          E(i, j) = 1.0;
        }
        arma::mat dSigma = E * L.t() + L * E.t(); // product rule
        J.col(idx_param) = vech(dSigma);
      }
    }

  } else {
    // diagonal Sigma only (no correlations)
    // Jacobian of Sigma wrt L_params  (diagonal-only case)
    for (int k = 0; k < K_w; ++k) {
      // Sigma_kk = L_kk^2  ,  L_kk = exp(z_k)
      // dSigma_kk/dz_k = 2 * exp(2 z_k) = 2 * L_kk^2
      double deriv = 2.0 * L(k, k) * L(k, k);
      J(k, k) = deriv;
    }
  }
  return J;
}

// ============================================================================
// Threads and per-thread buffers (estimation and prediction kernels)
// ============================================================================

// The team of a kernel's parallel regions, and so the number of scratch sets
// or accumulators it allocates: OpenMP's next team, within the thread limit,
// and at most one thread per work item (likelihood unit or choice situation),
// but two for a single item, so that the region stays active and, as before,
// an OpenMP-built BLAS runs the item's products single-threaded.
inline int mxl_team_threads(const mxl_off n_items) {
#ifdef _OPENMP
  mxl_off n = std::min(omp_get_max_threads(), omp_get_thread_limit());
#else
  mxl_off n = 1;
#endif
  n = std::min(n, std::max<mxl_off>(n_items, 2));
  return static_cast<int>(std::max<mxl_off>(n, 1));
}

inline int mxl_thread_num() {
#ifdef _OPENMP
  return omp_get_thread_num();
#else
  return 0;
#endif
}

// An uninitialized array of n doubles, the thread that uses it touching it
// first, with 128 bytes of padding (a cache line on Apple silicon, two on
// x86-64) so that no two threads' buffers share one.
inline std::unique_ptr<double[]> mxl_buffer(const std::size_t n) {
  return std::unique_ptr<double[]>(n > 0 ? new double[n + 16] : nullptr);
}

// ============================================================================
// Compact derivative blocks
//
// A unit's score and Hessian are nonzero only on the continuous parameters
// c = [beta | mu | L] (K_c = idx_delta_start of them) and on the free ASCs of
// the alternatives in its rows (a_u of them): about 75 of 2,620 parameters for
// a student who chose among 10-100 of 2,602 schools. The derivative kernels
// work in that block, m_u = K_c + a_u, and add it into the n x n result
// (n = n_params) once per unit. A unit's free ASCs, numbered by increasing
// global index, are its slots; within a situation the Hessian works on the
// situation's own slots, so that a situation costs O((K_c + a_t)^2 S) with
// a_t <= m_t its distinct free ASCs, not O(n^2 S) (nor O((K_c + a_u)^2 S): a
// patient's a_u, over all their visits, can be far above any one visit's).
// The free ASC of a row is its alternative's: none without ASCs, none for
// the first inside alternative without an outside option (the normalized
// reference), and the outside option has none. A situation that lists an
// alternative twice adds both rows into one slot, in row order, as the dense
// code added them into one delta.
// ============================================================================

// Free ASC (0-based index into the delta block) of stacked row r, or -1.
inline int mxl_free_delta(const MxlUnitData& ud, const mxl_off r) {
  if (!ud.use_asc) return -1;
  const int id = ud.lay.alt0(r);
  return ud.include_outside_option ? id : id - 1;
}

// Slot maps of the unit and the situation in hand. One per thread, sized on
// the primary thread for the largest unit and situation, so that building
// them allocates nothing; the markers (unit_of, sit_of) are -1 except while a
// unit or situation is mapped.
struct MxlUnitMap {
  std::vector<int> unit_of;   // free ASC j -> unit slot, or -1 (J_d entries)
  std::vector<int> delta;     // unit slot -> free ASC, increasing
  std::vector<int> row_slot;  // row of the unit -> unit slot, or -1
  std::vector<int> sit_of;    // unit slot -> situation slot, or -1
  std::vector<int> sit_slot;  // situation slot -> unit slot, increasing
  std::vector<int> row_sit;   // row of the situation -> situation slot, or -1
  int a = 0;                  // free ASCs of the unit
  int a_t = 0;                // free ASCs of the situation
  char pad[128];              // keeps the next thread's map off these lines

  // Each array is reserved 32 entries (128 bytes) beyond its largest size,
  // so that no two threads' arrays share a cache line.
  MxlUnitMap(const int J_d, const mxl_off max_unit_rows, const mxl_off max_m)
      : unit_of(static_cast<std::size_t>(J_d) + 32, -1) {
    const std::size_t a_max = static_cast<std::size_t>(
        std::min<mxl_off>(J_d, max_unit_rows));
    delta.reserve(a_max + 32);
    sit_of.reserve(a_max + 32);
    sit_of.assign(a_max, -1);
    sit_slot.reserve(static_cast<std::size_t>(std::min<mxl_off>(J_d, max_m)) + 32);
    row_slot.reserve(static_cast<std::size_t>(max_unit_rows) + 32);
    row_sit.reserve(static_cast<std::size_t>(max_m) + 32);
  }
  // Global parameter index of unit-local index i (K_c continuous first).
  int global(const int i, const int K_c) const {
    return i < K_c ? i : K_c + delta[i - K_c];
  }
  // Unit-local index of situation-local index i.
  int unit(const int i, const int K_c) const {
    return i < K_c ? i : K_c + sit_slot[i - K_c];
  }
};

// Every thread's maps, each in its own allocation, on the primary thread.
using MxlUnitMaps = std::vector<std::unique_ptr<MxlUnitMap>>;
inline MxlUnitMaps mxl_unit_maps(const MxlUnitData& ud, const int n_threads) {
  const int J_d = ud.n_params - ud.par.idx_delta_start;
  MxlUnitMaps maps;
  maps.reserve(n_threads);
  for (int i = 0; i < n_threads; ++i) {
    maps.emplace_back(new MxlUnitMap(J_d, ud.lay.max_unit_rows, ud.lay.max_m));
  }
  return maps;
}

// The slots of the loaded unit: its distinct free ASCs, in increasing order,
// and the slot of each of its rows. mxl_unit_map_clear() resets the markers.
inline void mxl_unit_map(const MxlUnitData& ud, const MxlUnitScratch& sc,
                         MxlUnitMap& mp) {
  mp.delta.clear();
  for (mxl_off i = 0; i < sc.R; ++i) {
    const int j = mxl_free_delta(ud, sc.r0 + i);
    if (j >= 0 && mp.unit_of[j] < 0) {
      mp.unit_of[j] = 0;  // seen
      mp.delta.push_back(j);
    }
  }
  std::sort(mp.delta.begin(), mp.delta.end());
  mp.a = static_cast<int>(mp.delta.size());
  for (int k = 0; k < mp.a; ++k) mp.unit_of[mp.delta[k]] = k;
  mp.row_slot.resize(static_cast<std::size_t>(sc.R));
  for (mxl_off i = 0; i < sc.R; ++i) {
    const int j = mxl_free_delta(ud, sc.r0 + i);
    mp.row_slot[i] = j >= 0 ? mp.unit_of[j] : -1;
  }
}

inline void mxl_unit_map_clear(MxlUnitMap& mp) {
  for (int k = 0; k < mp.a; ++k) mp.unit_of[mp.delta[k]] = -1;
}

// The slots of situation t of the loaded unit: its distinct unit slots, in
// increasing order, and the situation slot of each of its rows.
// mxl_situation_map_clear() resets the markers.
inline void mxl_situation_map(const MxlUnitData& ud, const MxlUnitScratch& sc,
                              const mxl_off t, MxlUnitMap& mp) {
  const mxl_off first = ud.lay.row_off[t] - sc.r0;  // first row in the unit
  const int m = ud.lay.m(t);
  mp.sit_slot.clear();
  for (int a = 0; a < m; ++a) {
    const int k = mp.row_slot[first + a];
    if (k >= 0 && mp.sit_of[k] < 0) {
      mp.sit_of[k] = 0;  // seen
      mp.sit_slot.push_back(k);
    }
  }
  std::sort(mp.sit_slot.begin(), mp.sit_slot.end());
  mp.a_t = static_cast<int>(mp.sit_slot.size());
  for (int q = 0; q < mp.a_t; ++q) mp.sit_of[mp.sit_slot[q]] = q;
  mp.row_sit.resize(static_cast<std::size_t>(m));
  for (int a = 0; a < m; ++a) {
    const int k = mp.row_slot[first + a];
    mp.row_sit[a] = k >= 0 ? mp.sit_of[k] : -1;
  }
}

inline void mxl_situation_map_clear(MxlUnitMap& mp) {
  for (int q = 0; q < mp.a_t; ++q) mp.sit_of[mp.sit_slot[q]] = -1;
}

// The n x n result of a derivative kernel (Hessian, BHHH, cluster meat): the
// units' compact blocks are added into its upper triangle, which is mirrored
// at the end (every block is exactly symmetric, so the mirror is the sum the
// lower triangle would have accumulated). With one thread the blocks go into
// the result itself; with more, each thread adds into its own packed upper
// triangle (column j's rows 0..j from j (j + 1) / 2), and the triangles are
// added up in thread order, starting from +0 as the critical section's
// addition into a zero matrix did. Everything is allocated here, on the
// primary thread, so that running out of memory is an R error saying how much
// each thread needs, not a failure inside a parallel region.
struct MxlSymAcc {
  const mxl_off n;
  const int T;
  Rcpp::NumericMatrix result;  // n x n, zero
  double* out = nullptr;       // its elements
  std::vector<std::unique_ptr<double[]>> part;  // T > 1: packed triangles

  static std::size_t packed_len(const mxl_off n) {
    return static_cast<std::size_t>(n) * static_cast<std::size_t>(n + 1) / 2;
  }

  MxlSymAcc(const int n_, const int T_, const char* what)
      : n(n_), T(T_), result(n_, n_) {
    out = result.begin();
    if (T > 1) {
      const std::size_t len = packed_len(n);
      try {
        part.reserve(T);
        for (int i = 0; i < T; ++i) part.push_back(mxl_buffer(len));
      } catch (const std::bad_alloc&) {
        std::vector<std::unique_ptr<double[]>>().swap(part);  // release
        Rcpp::stop("Not enough memory for the %s accumulators: %.2f GB per "
                   "thread for %d threads. Run fewer threads with "
                   "set_num_threads().", what,
                   sizeof(double) * static_cast<double>(len) / 1e9, T);
      }
      // Thread i zeros triangle i (schedule(static)), touching first the
      // memory it will add to.
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(T)
#endif
      for (int i = 0; i < T; ++i) {
        std::fill(part[i].get(), part[i].get() + len, 0.0);
      }
    }
  }
  MxlSymAcc(const MxlSymAcc&) = delete;
  MxlSymAcc& operator=(const MxlSymAcc&) = delete;

  // Thread tid's accumulator and the offset of entry (i, j), i <= j.
  double* base(const int tid) const { return T > 1 ? part[tid].get() : out; }
  std::size_t at(const mxl_off i, const mxl_off j) const {
    return T > 1 ? static_cast<std::size_t>(j) * static_cast<std::size_t>(j + 1) / 2 +
                       static_cast<std::size_t>(i)
                 : static_cast<std::size_t>(i) +
                       static_cast<std::size_t>(j) * static_cast<std::size_t>(n);
  }

  // The upper triangle of the result: the threads' triangles added up in
  // thread order (a no-op with one thread); then the lower triangle as its
  // mirror, every element negated when `negate` (the Hessian's sign).
  void finish(const bool negate) {
    const mxl_off nn = n;
    double* o = out;
    const std::vector<std::unique_ptr<double[]>>& pt = part;
    const int n_part = T > 1 ? T : 0;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 16) num_threads(T)
#endif
    for (mxl_off j = 0; j < nn; ++j) {
      const std::size_t col = static_cast<std::size_t>(j) * static_cast<std::size_t>(nn);
      if (n_part > 0) {
        // Column j's entries from +0, then each thread's in turn: each entry
        // adds the threads in thread order, reading one triangle at a time.
        const std::size_t pk = static_cast<std::size_t>(j) * static_cast<std::size_t>(j + 1) / 2;
        double* oc = o + col;
        for (mxl_off i = 0; i <= j; ++i) oc[i] = 0.0;
        for (int k = 0; k < n_part; ++k) {
          const double* pc = pt[k].get() + pk;
          for (mxl_off i = 0; i <= j; ++i) oc[i] += pc[i];
        }
      }
      for (mxl_off i = 0; i <= j; ++i) {
        const double v = negate ? -o[col + static_cast<std::size_t>(i)]
                                : o[col + static_cast<std::size_t>(i)];
        o[col + static_cast<std::size_t>(i)] = v;
        o[static_cast<std::size_t>(j) + static_cast<std::size_t>(i) * static_cast<std::size_t>(nn)] = v;
      }
    }
  }
};

// Non-finite blocks. In the dense code a unit's block sat in an n x n matrix
// of zeros, and two operations spread a non-finite value beyond the block. A
// row of the score stash F_t or of the centered G holding a NaN or an
// infinity gave NaN in that row's and that column's entries outside the
// block, through 0 * NaN and 0 * Inf in the BLAS product (as OpenBLAS forms
// it; a BLAS that skips zero operands spread less). And a non-finite unit
// weight multiplied the zeros outside the unit's block: w * 0 is NaN, or NA
// for an NA weight. The compact kernels record both in per-thread arrays
// allocated on the primary thread, a flag per parameter and a count per
// parameter of the non-finite-weight units whose block has it, so that
// nothing grows inside the parallel region however many units are affected;
// mxl_apply_nonfinite() reproduces them on the result in one pass: whole rows
// and columns of NaN for the first, and every entry outside the intersection
// of those units' blocks for the second.
struct MxlNonFinite {
  std::vector<char> row;  // parameter with a non-finite row of F_t or G
  std::vector<int> hits;  // non-finite-weight units whose block has it
  int n_blk = 0;          // non-finite-weight units
  bool na = false;        // one of their weights is NA
  explicit MxlNonFinite(const int n)
      : row(static_cast<std::size_t>(n), 0), hits(static_cast<std::size_t>(n), 0) {}
};

inline std::vector<MxlNonFinite> mxl_nonfinite(const int n, const int n_threads) {
  std::vector<MxlNonFinite> nf;
  nf.reserve(n_threads);
  for (int i = 0; i < n_threads; ++i) nf.emplace_back(n);
  return nf;
}

// Whether x is R's NA_real_: a NaN whose low-order 32 bits are 1954, the test
// of R_IsNA(), here without the R API so that threads may call it.
inline bool mxl_is_na(const double x) {
  if (!std::isnan(x)) return false;
  uint64_t bits;
  std::memcpy(&bits, &x, sizeof bits);
  return (bits & 0xFFFFFFFFu) == 1954u;
}

// Flag the parameters of the rows of a unit-local or situation-local matrix B
// (rows x S) that hold a non-finite value; `to_global` maps a row of B to its
// global index.
template <typename Map>
inline void mxl_nonfinite_rows(const arma::mat& B, Map to_global,
                               MxlNonFinite& nf) {
  if (B.is_finite()) return;
  for (arma::uword i = 0; i < B.n_rows; ++i) {
    if (!B.row(i).is_finite()) nf.row[to_global(static_cast<int>(i))] = 1;
  }
}

// Record a unit with a non-finite weight w over its block of m parameters.
template <typename Map>
inline void mxl_nonfinite_weight(const double w, const int m, Map to_global,
                                 MxlNonFinite& nf) {
  for (int i = 0; i < m; ++i) ++nf.hits[to_global(i)];
  ++nf.n_blk;
  if (mxl_is_na(w)) nf.na = true;
}

// Apply the threads' records to the finished result.
inline void mxl_apply_nonfinite(MxlSymAcc& acc,
                                const std::vector<MxlNonFinite>& nf) {
  const mxl_off n = acc.n;
  std::vector<char> row(static_cast<std::size_t>(n), 0);
  std::vector<int> hits(static_cast<std::size_t>(n), 0);
  int n_blk = 0;
  bool na = false, any_row = false;
  for (const MxlNonFinite& t : nf) {
    for (mxl_off g = 0; g < n; ++g) {
      row[g] = row[g] | t.row[g];
      hits[g] += t.hits[g];
      any_row = any_row || t.row[g];
    }
    n_blk += t.n_blk;
    na = na || t.na;
  }
  if (!any_row && n_blk == 0) return;
  const double nan = std::numeric_limits<double>::quiet_NaN();
  const double spread = na ? NA_REAL : nan;
  double* o = acc.out;
  for (mxl_off j = 0; j < n; ++j) {
    const std::size_t col = static_cast<std::size_t>(j) * static_cast<std::size_t>(n);
    for (mxl_off i = 0; i < n; ++i) {
      if (row[i] || row[j]) {
        o[col + static_cast<std::size_t>(i)] = nan;
      } else if (n_blk > 0 && !(hits[i] == n_blk && hits[j] == n_blk)) {
        o[col + static_cast<std::size_t>(i)] = spread;
      }
    }
  }
}

//' Analytical Hessian of the log-likelihood v2
//'
//' Computes the Hessian of the log-likelihood for the Mixed Logit model using
//' OpenMP for parallelization. Mirrors the parameters of mxl_loglik_gradient_parallel.
//'
//' @param theta vector collecting model parameters (beta, mu, L, delta (ASCs))
//' @param X design matrix for covariates with fixed coefficients; sum(M_i) x K_x
//' @param W design matrix for covariates with random coefficients; sum(M_i) x K_w or J x K_w
//' @param alt_idx sum(M) x 1 vector with indices of alternatives within each choice set; 1-based indexing
//' @param choice_idx N x 1 vector with indices of chosen alternatives; 1-based indexing relative to X; 0 is used if include_outside_option=True
//' @param M N x 1 vector with number of alternatives for each individual
//' @param weights N x 1 vector with weights for each observation; when Ti is
//'   supplied they must be constant within each decision maker
//' @param eta_draws Array of standard-normal draws, K_w x S x U, where U is the
//'   number of decision makers when Ti is supplied and the number of choice
//'   situations otherwise
//' @param rc_dist K_w x 1 integer vector indicating distribution of random coefficients: 0 = normal, 1 = log-normal
//' @param rc_correlation whether random coefficients should be correlated
//' @param rc_mean whether to estimate means for random coefficients.
//' @param use_asc whether to use alternative-specific constants.
//' @param include_outside_option whether to include outside option normalized to 0 (if so, the outside option is not included in the data)
//' @param gen_seed Integer master seed for the on-the-fly Halton generator. \code{< 0}
//'   (default) uses the materialized \code{eta_draws} cube; \code{>= 0} generates draws
//'   on the fly from this seed.
//' @param gen_scramble Integer scramble mode for on-the-fly generation: \code{0} =
//'   identity permutations (plain Halton, compat), \code{1} = seeded position-wise
//'   digit permutations; other values are an error.
//' @param gen_S Integer number of draws per individual, used only when \code{gen_seed >= 0}.
//' @param Ti Optional integer vector with the number of choice situations of
//'   each decision maker (panel likelihood); situations must be sorted by
//'   decision maker. NULL (default): every choice situation is its own unit
//'   (cross-sectional likelihood).
//' @param draw_batch Integer; \code{0} (default) forms each decision maker's
//'   draws in batches sized to a per-thread memory budget, a positive value
//'   caps the number of draws per batch (for tests).
//' @returns Hessian evaluated at input arguments
//' @note For log-normal random coefficients (rc_dist=1) with rc_mean=TRUE,
//'   the distribution is a shifted log-normal: beta_k = exp(mu_k) + exp(L_k * eta),
//'   where exp(mu_k) shifts the location and exp(L_k * eta) ~ LogNormal(0, sigma_k^2).
//'   This differs from the textbook parameterization exp(mu_k + L_k * eta).
//' @examples
//' \donttest{
//' library(data.table)
//' set.seed(42)
//' N <- 50; J <- 3
//' dt <- data.table(id = rep(1:N, each = J), alt = rep(1:J, N))
//' dt[, `:=`(x1 = rnorm(.N), w1 = rnorm(.N))]
//' dt[, choice := 0L]
//' dt[, choice := sample(c(1L, rep(0L, J - 1))), by = id]
//' d <- prepare_mxl_data(dt, "id", "alt", "choice", "x1", "w1")
//' eta <- get_halton_normals(50, d$N, ncol(d$W))
//' theta <- rep(0, ncol(d$X) + ncol(d$W) + nrow(d$alt_mapping) - 1)
//' H <- choicer:::mxl_hessian_parallel(theta, d$X, d$W, d$alt_idx, d$choice_idx,
//'   d$M, d$weights, eta, rc_dist = rep(0L, ncol(d$W)),
//'   rc_correlation = FALSE, rc_mean = FALSE)
//' dim(H)
//' }
//' @keywords internal
// [[Rcpp::export]]
Rcpp::NumericMatrix mxl_hessian_parallel(
    const arma::vec &theta, const arma::mat &X, const arma::mat &W,
    const Rcpp::IntegerVector &alt_idx,
    const Rcpp::IntegerVector &choice_idx,
    const Rcpp::IntegerVector &M, const arma::vec &weights,
    const arma::cube &eta_draws, const arma::uvec &rc_dist,
    const bool rc_correlation = true, const bool rc_mean = false,
    const bool use_asc = true, const bool include_outside_option = false,
    const int gen_seed = -1, const int gen_scramble = 1, const int gen_S = 0,
    const Rcpp::Nullable<Rcpp::IntegerVector> Ti = R_NilValue,
    const int draw_batch = 0) {
  // Inputs, layout and parameters at theta, validated on the primary thread
  const MxlUnitData ud(theta, X, W, alt_idx, choice_idx, M,
                       &weights, eta_draws, rc_dist,
                       rc_correlation, rc_mean, use_asc, include_outside_option,
                       gen_seed, gen_scramble, gen_S, Ti, draw_batch,
                       "mxl_hessian_parallel");
  const MxlLayout& lay = ud.lay;
  const int n_params = ud.n_params;
  const int K_w = W.n_cols;
  const int Sdraw = ud.S;
  // Reuse the parser's layout when assembling derivatives.
  const MxlParams& par = ud.par;
  const int idx_beta_start = par.idx_beta_start;
  const int idx_mu_start = par.idx_mu_start;
  const int idx_L_start = par.idx_L_start;
  const int idx_delta_start = par.idx_delta_start;
  const arma::mat& L = par.L;
  const arma::vec& dmu_final_dmu = par.dmu_final_dmu;
  const arma::vec& dmu2_final_dmu2 = par.dmu2_final_dmu2;

  // Block layout: continuous block c = [beta | mu | L], size Kc = idx_delta_start,
  // then the unit's free ASCs (its slots, see "Compact derivative blocks"):
  // a unit's Hessian lives on its Kc + a_u parameters, a situation's
  // per-draw pieces on its Kc + a_t. H_V (second derivative of utility
  // w.r.t. theta) is nonzero only in the (mu, L) x (mu, L) sub-block, which
  // lives entirely within the continuous block.
  const int Kc = idx_delta_start;  // size of continuous block

  // Per unit u (Louis 1982), with omega_s the posterior draw weights and
  // g_s = sum_t g_ts, H_s = sum_t H_ts the per-draw score and Hessian of
  // lambda_s = sum_t log P_ts(j_t):
  //   H_u = sum_s omega_s (H_s + (g_s - g_bar)(g_s - g_bar)'),  g_bar = sum_s omega_s g_s,
  //   H_ts = -sum_Pzz_ts + sum_Pz_ts sum_Pz_ts' + sum_diff_H_V_ts.
  // Pass 1 (the shared draw loop) yields omega; pass 2 accumulates the O3
  // buffers with omega_s where the cross-section used P_choice_s / P_i_hat.

  // The result and the threads' accumulators, slot maps and non-finite
  // records, all allocated here on the primary thread.
  const int n_threads = mxl_team_threads(lay.U);
  MxlSymAcc acc(n_params, n_threads, "Hessian");
  MxlUnitMaps maps = mxl_unit_maps(ud, n_threads);
  std::vector<MxlNonFinite> nonfinite = mxl_nonfinite(n_params, n_threads);

#ifdef _OPENMP
#pragma omp parallel num_threads(n_threads)
#endif
  {
    const int tid = mxl_thread_num();
    MxlUnitMap& mp = *maps[tid];
    MxlNonFinite& nf = nonfinite[tid];
    double* acc_t = acc.base(tid);

    // Thread-private scratch — sized per unit or situation, reset inside
    // loops as needed; Armadillo keeps a member's memory as it shrinks.
    MxlUnitScratch sc; // unit slices, draws, probabilities P_ts, weights omega

    // Per-(t, s) block accumulators (reset at the start of each situation-draw).
    // cc: Kc x Kc dense outer-product sum
    arma::mat sum_Pzz_cc(Kc, Kc);
    // cd: Kc x a_t; column q accumulates P_a * zc_a for the alternatives in
    // situation slot q
    arma::mat sum_Pzz_cd;
    // dd: diagonal only — stored as a length-a_t vector
    arma::vec sum_Pzz_dd;
    // P*z sums
    arma::vec sum_Pz_c(Kc);
    arma::vec sum_Pz_d;
    // gradient components
    arma::vec g_c(Kc);
    arma::vec g_d;
    // H_V in the continuous block (mu,L sub-block; beta and delta rows/cols are zero)
    arma::mat sum_diff_H_V_cc(Kc, Kc);

    // Per-alt scratch for the continuous-block z vector and its outer product
    arma::vec zc_a(Kc);
    arma::mat zz(Kc, Kc);

    // O3: Per-unit block buffers — accumulate the pieces linear in omega_s
    // across the unit's situations and draws, in unit slots.
    arma::mat buf_Pzz_cc(Kc, Kc);  // Identity C: Σ_ts ω_s sum_Pzz_cc_ts
    arma::mat buf_Pzz_cd;          // Identity C: Σ_ts ω_s sum_Pzz_cd_ts (Kc x a_u)
    arma::vec buf_Pzz_dd;          // Identity C: Σ_ts ω_s sum_Pzz_dd_ts (a_u)
    arma::mat buf_diff_HV_cc(Kc, Kc);  // Identity C: Σ_ts ω_s sum_diff_H_V_cc_ts
    // O3: Column stashes for BLAS-3 batching.
    // G_stash(:,s) = sqrt(ω_s) * (g_s - g_bar), g_s = Σ_t g_ts (unit slots)
    //   → G Gᵀ = Σ_s ω_s (g_s - g_bar)(g_s - g_bar)ᵀ              (Identity A)
    // F_stash(:,s) = sqrt(ω_s) * [sum_Pz_c; sum_Pz_d]_ts (situation slots)
    //   → F Fᵀ = Σ_s ω_s sum_Pz sum_PzT for situation t           (Identity B)
    // F column s is fully written in draw s of every situation; G is zeroed
    // per unit and accumulates over the unit's situations.
    arma::mat G_stash;
    arma::mat F_stash;
    arma::mat opg_pz;         // G Gᵀ + Σ_t F_t F_tᵀ, then the unit's Hessian
    arma::vec g_bar;          // reference score, then mean deviation
    // Per-unit assembly buffers, reused so that no unit or situation
    // allocates once the thread has seen a larger one: the products F Fᵀ and
    // G Gᵀ, and sqrt(ω). The unit's Hessian is assembled in place in opg_pz.
    arma::mat prod_buf;
    arma::rowvec sqrt_omega(Sdraw);

#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
    for (mxl_off u = 0; u < lay.U; ++u) {
      const double w_u = weights[lay.unit_first(u)]; // unit weight

      // --- Pass 1: posterior draw weights omega_s of the unit's choices ---
      mxl_unit_load(ud, u, sc, 1, true, true);
      const double lse = mxl_unit_simulate(ud, sc, false);
      // lse is finite whenever the utilities are (log-space probabilities),
      // so this skips only a unit whose utilities overflowed, where its
      // log-likelihood is undefined.
      if (!std::isfinite(lse)) continue;
      const arma::mat eta_u = mxl_eta_view(ud, sc);

      // The unit's slots: m_u = Kc + a_u local parameters.
      mxl_unit_map(ud, sc, mp);
      const int a_u = mp.a;
      const int m_u = Kc + a_u;

      // O3: Initialize per-unit block buffers.
      buf_Pzz_cc.zeros();
      buf_diff_HV_cc.zeros();
      buf_Pzz_cd.zeros(Kc, a_u);
      buf_Pzz_dd.zeros(a_u);
      G_stash.zeros(m_u, Sdraw);
      opg_pz.zeros(m_u, m_u);

      // --- Pass 2: situations t (outer) x draws s (inner). A unit in one
      // draw batch reuses pass 1's WGamma; for a longer one each situation's
      // rows of W_u Gamma_u (m_t x S) are formed in turn, so the thread's
      // memory grows with the largest choice set rather than the unit.
      const bool one_batch = sc.nb == 1;
      for (mxl_off t = sc.t0; t < sc.t1; ++t) {
        const int m_t = static_cast<int>(lay.row_off[t + 1] - lay.row_off[t]);
        const int num_choices = include_outside_option ? m_t + 1 : m_t;
        if (!one_batch) mxl_situation_wgamma(sc, lay.row_off[t] - sc.r0, m_t);

        // The situation's slots: its per-draw pieces live on Kc + a_t.
        mxl_situation_map(ud, sc, t, mp);
        const int a_t = mp.a_t;
        sum_Pzz_cd.set_size(Kc, a_t);
        sum_Pzz_dd.set_size(a_t);
        sum_Pz_d.set_size(a_t);
        g_d.set_size(a_t);
        F_stash.set_size(Kc + a_t, Sdraw);

        // chosen alternative index (validated with the layout)
        const int chosen_alt = lay.chosen(t);

        // Loop over simulations (draws) s
        for (int s = 0; s < Sdraw; ++s) {
          // Column views into the batched matrices (zero-copy)
          const auto eta_i_s                = eta_u.col(s);
          const auto dgamma_final_dgamma    = sc.Dgamma1.col(s);
          const auto d2gamma_final_dgamma2  = sc.Dgamma2.col(s);
          const double omega_s = sc.omega(s);

          // === 3. Utility and probabilities P_ts (as in pass 1) ===
          mxl_situation_probs(ud, sc, t, s);
          const arma::vec& P_s = sc.P;

          // === 4. Per-(t, s) Gradient (g_ts) and Hessian (H_ts) — block form ===
          // Reset per-draw block accumulators.
          sum_Pzz_cc.zeros();
          sum_Pz_c.zeros();
          g_c.zeros();
          sum_diff_H_V_cc.zeros();
          if (a_t > 0) {
            sum_Pzz_cd.zeros();
            sum_Pzz_dd.zeros();
            sum_Pz_d.zeros();
            g_d.zeros();
          }

          for (int a = 0; a < num_choices; ++a) {
            // Determine if this alt contributes a nonzero z (outside option at a==0 skips)
            if (include_outside_option && a == 0) {
              // Outside option: z_as == 0 and H_V_as == 0, so this alternative
              // contributes nothing to the Hessian accumulators — skip entirely.
              continue;
            }

            const int current_a_idx = include_outside_option ? a - 1 : a;
            const mxl_off row = lay.row_off[t] + current_a_idx; // stacked row of alt a
            const arma::rowvec w_ap_row = sc.W_u.row(row - sc.r0);

            // --- Build continuous-block vector zc_a ---
            zc_a.zeros();

            // beta sub-block
            zc_a.subvec(idx_beta_start, idx_mu_start - 1) =
                sc.X_u.row(row - sc.r0).t();

            // mu sub-block  (nonzero only when rc_mean)
            if (rc_mean) {
              for (int p = 0; p < K_w; ++p) {
                zc_a[idx_mu_start + p] = w_ap_row(p) * dmu_final_dmu(p);
              }
            }

            // L sub-block
            if (rc_correlation) {
              int lp_idx_r = 0;
              for (int p = 0; p < K_w; ++p) {
                const double f_p_prime = dgamma_final_dgamma(p);
                for (int q = 0; q <= p; ++q, ++lp_idx_r) {
                  const double dLpq_dparam_r = (p == q) ? L(p, p) : 1.0;
                  zc_a[idx_L_start + lp_idx_r] =
                      w_ap_row(p) * f_p_prime * dLpq_dparam_r * eta_i_s(q);
                }
              }
            } else { // Diagonal L
              for (int p = 0; p < K_w; ++p) {
                zc_a[idx_L_start + p] =
                    w_ap_row(p) * dgamma_final_dgamma(p) * L(p, p) * eta_i_s(p);
              }
            }

            // --- Situation slot of this alt's free ASC (-1: none) ---
            const int k_sit = mp.row_sit[current_a_idx];

            // --- Accumulate block sums ---
            const double P_a = P_s(a);
            const double diff = (a == chosen_alt ? 1.0 : 0.0) - P_a;

            // cc block: P_a * zc_a * zc_a^T  (rank-1 update)
            zz = zc_a * zc_a.t();
            sum_Pzz_cc += P_a * zz;
            sum_Pz_c   += P_a * zc_a;
            g_c        += diff * zc_a;

            // cd and dd blocks (delta scatter)
            if (k_sit >= 0) {
              sum_Pzz_cd.col(k_sit) += P_a * zc_a;
              sum_Pzz_dd(k_sit)     += P_a;
              sum_Pz_d(k_sit)       += P_a;
              g_d(k_sit)            += diff;
            }

            // H_V accumulation — nonzero only in (mu,L)x(mu,L) within cc block
            // Reuse same logic as before but write directly into sum_diff_H_V_cc
            if (rc_mean) {
              for (int p = 0; p < K_w; ++p) {
                sum_diff_H_V_cc(idx_mu_start + p, idx_mu_start + p) +=
                    diff * w_ap_row(p) * dmu2_final_dmu2(p);
              }
            }
            if (rc_correlation) {
              int lp_idx_r = 0;
              for (int p = 0; p < K_w; ++p) {
                const double f_p_prime        = dgamma_final_dgamma(p);
                const double f_p_double_prime = d2gamma_final_dgamma2(p);
                for (int q = 0; q <= p; ++q, ++lp_idx_r) {
                  const int r = idx_L_start + lp_idx_r;
                  const double dLpq_dparam_r = (p == q) ? L(p, p) : 1.0;
                  const double d2Lpq_dparam2_r = (p == q) ? L(p, p) : 0.0;
                  const double dgamma_p_dparam_r = dLpq_dparam_r * eta_i_s(q);
                  const double d2gamma_p_dparam2_r = d2Lpq_dparam2_r * eta_i_s(q);

                  // Diagonal H_V entry
                  sum_diff_H_V_cc(r, r) +=
                      diff * w_ap_row(p) *
                      (f_p_double_prime * std::pow(dgamma_p_dparam_r, 2) +
                       f_p_prime * d2gamma_p_dparam2_r);

                  // Off-diagonal L-L entries (same row p, column < q)
                  for (int q_s = 0; q_s < q; ++q_s) {
                    const int param_s = idx_L_start + lp_idx_r - (q - q_s);
                    const double dLpq_s_dparam = (p == q_s) ? L(p, p) : 1.0;
                    const double dgamma_p_dparam_s = dLpq_s_dparam * eta_i_s(q_s);
                    const double d2V = diff * w_ap_row(p) * f_p_double_prime *
                                       dgamma_p_dparam_r * dgamma_p_dparam_s;
                    sum_diff_H_V_cc(r, param_s) += d2V;
                    sum_diff_H_V_cc(param_s, r) += d2V;
                  }
                }
              }
            } else { // Diagonal L
              for (int p = 0; p < K_w; ++p) {
                const int r = idx_L_start + p;
                const double f_p_prime        = dgamma_final_dgamma(p);
                const double f_p_double_prime = d2gamma_final_dgamma2(p);
                const double dgamma_p_dparam  = L(p, p) * eta_i_s(p);
                const double d2gamma_p_dparam = L(p, p) * eta_i_s(p);
                sum_diff_H_V_cc(r, r) +=
                    diff * w_ap_row(p) *
                    (f_p_double_prime * std::pow(dgamma_p_dparam, 2) +
                     f_p_prime * d2gamma_p_dparam);
              }
            }
          } // end alt loop

          // === 5. O3: Accumulate per-unit buffers (Identity C) and fill
          //          column stashes G_stash/F_stash (Identities A + B). ===
          //
          // Instead of assembling H_ts and accumulating
          //   hess_term1 += omega_s * ((g_s - g_bar)(g_s - g_bar)T + sum_t H_ts),
          // we split the linear-in-omega_s and the outer-product pieces:
          //
          //   Linear (Identity C): buf_Pzz_cc     += omega_s * sum_Pzz_cc
          //                        buf_diff_HV_cc += omega_s * sum_diff_H_V_cc
          //                        (and cd/dd variants, situation slots to
          //                        unit slots)
          //   Outer-product (Identity A): G_stash(:,s) += g_ts, centered at
          //     g_bar and scaled by sqrt(omega_s) after the unit, so that
          //     G GT = Σ_s omega_s (g_s - g_bar)(g_s - g_bar)T.
          //   Outer-product (Identity B): F_stash(:,s) = sqrt(omega_s) * [sum_Pz_c; sum_Pz_d]
          //     so that F FT = Σ_s omega_s sum_Pz sum_PzT after situation t.
          // Each addition into a unit buffer is the one the dense buffers made
          // (their other columns gained omega_s * 0), its product formed first
          // as Armadillo's element-wise update forms it.

          const double sqrt_ws = std::sqrt(omega_s);

          // Identity C — scalar-times-matrix accumulation into per-unit buffers.
          buf_Pzz_cc    += omega_s * sum_Pzz_cc;
          buf_diff_HV_cc += omega_s * sum_diff_H_V_cc;
          for (int q = 0; q < a_t; ++q) {
            const int k_u = mp.sit_slot[q];
            buf_Pzz_cd.col(k_u) += omega_s * sum_Pzz_cd.col(q);
            const double v = sum_Pzz_dd(q) * omega_s;
            buf_Pzz_dd(k_u) += v;
          }

          // Identity A — accumulate g_s = Σ_t g_ts in column s (g_ts = [g_c; g_d]).
          G_stash.col(s).head(Kc) += g_c;
          for (int q = 0; q < a_t; ++q) G_stash(Kc + mp.sit_slot[q], s) += g_d(q);

          // Identity B — fill column s of F_stash with sqrt(omega_s) * [sum_Pz_c; sum_Pz_d].
          F_stash.col(s).head(Kc) = sqrt_ws * sum_Pz_c;
          if (a_t > 0) {
            F_stash.col(s).tail(a_t) = sqrt_ws * sum_Pz_d;
          }
        } // end S loop

        // Identity B for situation t (one BLAS-3 product per situation: the
        // outer products do not combine across situations before the product),
        // added into the unit block at the situation's slots.
        mxl_nonfinite_rows(F_stash, [&](const int i) {
          return mp.global(mp.unit(i, Kc), Kc);
        }, nf);
        prod_buf = F_stash * F_stash.t();
        const int m_ts = Kc + a_t;
        for (int j = 0; j < m_ts; ++j) {
          const int uj = mp.unit(j, Kc);
          for (int i = 0; i < m_ts; ++i) {
            opg_pz(mp.unit(i, Kc), uj) += prod_buf(i, j);
          }
        }
        mxl_situation_map_clear(mp);
      } // end situation loop

      // Identity A, centered: G Gᵀ = Σ_s omega_s (g_s - g_bar)(g_s - g_bar)ᵀ.
      // Remove the highest-weight draw's score before computing the mean,
      // so neither a large common score nor a zero-weight outlier can
      // amplify rounding in the centered covariance.
      g_bar = G_stash.col(sc.omega.index_max());
      G_stash.each_col() -= g_bar;
      g_bar = G_stash * sc.omega;
      G_stash.each_col() -= g_bar;
      for (int s = 0; s < Sdraw; ++s) sqrt_omega[s] = std::sqrt(sc.omega[s]);
      G_stash.each_row() %= sqrt_omega;
      mxl_nonfinite_rows(G_stash, [&](const int i) { return mp.global(i, Kc); },
                         nf);
      prod_buf = G_stash * G_stash.t();
      opg_pz += prod_buf;

      // === 6. O3: Per-unit finalization — assemble Hessian once from buffers ===
      // 6a. Assemble hess_term1 block by block, in place in opg_pz:
      //     hess_term1 = (G Gᵀ + Σ_t F_t F_tᵀ) + (-buf_Pzz) + buf_diff_HV_cc (cc block only)
      //     This is the batched equivalent of Σ_s omega_s (H_s +
      //     (g_s - g_bar)(g_s - g_bar)T), the unit's whole Hessian H_u.

      // cc block: contributions from OPG, sum-Pz outer product, -Pzz, and H_V.
      opg_pz.submat(0, 0, Kc - 1, Kc - 1) -= buf_Pzz_cc;
      opg_pz.submat(0, 0, Kc - 1, Kc - 1) += buf_diff_HV_cc;

      if (a_u > 0) {
        // cd block: -buf_Pzz_cd + opg_pz cd block (NOT symmetric — full rectangular).
        opg_pz.submat(0, Kc, Kc - 1, m_u - 1) -= buf_Pzz_cd;
        for (int j = 0; j < a_u; ++j) {        // dc = (cd)ᵀ
          for (int i = 0; i < Kc; ++i) opg_pz(Kc + j, i) = opg_pz(i, Kc + j);
        }

        // dd block: -diag(buf_Pzz_dd) + opg_pz dd block (sum_Pz_d outer products).
        opg_pz.submat(Kc, Kc, m_u - 1, m_u - 1).diag() -= buf_Pzz_dd;
      }

      // 6b. Louis identity: H_u = hess_term1, whose centered Identity A
      //     already subtracts g_bar g_barᵀ. The draw weights are normalized,
      //     so there is no division by the simulated P_u. Add w_u H_u into
      //     the thread's upper triangle (opg_pz is exactly symmetric), the
      //     product formed first as Armadillo's element-wise update forms it.
      for (int j = 0; j < m_u; ++j) {
        const int gj = mp.global(j, Kc);
        for (int i = 0; i <= j; ++i) {
          double v = opg_pz(i, j);
          v *= w_u;
          acc_t[acc.at(mp.global(i, Kc), gj)] += v;
        }
      }
      if (!std::isfinite(w_u)) {
        mxl_nonfinite_weight(w_u, m_u, [&](const int i) { return mp.global(i, Kc); },
                             nf);
      }
      mxl_unit_map_clear(mp);
    } // end unit loop
  } // end parallel region

  acc.finish(true);  // the negated Hessian
  mxl_apply_nonfinite(acc, nonfinite);
  return acc.result;
}

//' BHHH (outer product of gradients) information matrix for Mixed Logit
//'
//' Computes the BHHH approximation to the observed information matrix for the
//' Mixed Logit model: \eqn{H_{BHHH} = \sum_i w_i \cdot s_i s_i^\top}, where
//' \eqn{s_i} is the score of likelihood unit i (gradient of its simulated
//' log-likelihood \eqn{\log \bar{P}_i}; a decision maker when Ti is supplied,
//' a choice situation otherwise).
//' This outer product of gradients (OPG) estimator provides an alternative to
//' the analytical Hessian for standard error computation that scales to large
//' problems where the analytical Hessian is infeasible (e.g., many alternatives
//' or simulation draws).
//'
//' @param theta vector collecting model parameters (beta, mu, L, delta (ASCs))
//' @param X design matrix for covariates with fixed coefficients; sum(M_i) x K_x
//' @param W design matrix for covariates with random coefficients; sum(M_i) x K_w or J x K_w
//' @param alt_idx sum(M) x 1 vector with indices of alternatives within each choice set; 1-based indexing
//' @param choice_idx N x 1 vector with indices of chosen alternatives; 1-based indexing relative to X; 0 is used if include_outside_option=True
//' @param M N x 1 vector with number of alternatives for each individual
//' @param weights N x 1 vector with weights for each observation; when Ti is
//'   supplied they must be constant within each decision maker
//' @param eta_draws Array of standard-normal draws, K_w x S x U, where U is the
//'   number of decision makers when Ti is supplied and the number of choice
//'   situations otherwise
//' @param rc_dist K_w x 1 integer vector indicating distribution of random coefficients: 0 = normal, 1 = log-normal
//' @param rc_correlation whether random coefficients should be correlated
//' @param rc_mean whether to estimate means for random coefficients.
//' @param use_asc whether to use alternative-specific constants.
//' @param include_outside_option whether to include outside option normalized to 0 (if so, the outside option is not included in the data)
//' @param gen_seed Integer master seed for the on-the-fly Halton generator. \code{< 0}
//'   (default) uses the materialized \code{eta_draws} cube; \code{>= 0} generates draws
//'   on the fly from this seed.
//' @param gen_scramble Integer scramble mode for on-the-fly generation: \code{0} =
//'   identity permutations (plain Halton, compat), \code{1} = seeded position-wise
//'   digit permutations; other values are an error.
//' @param gen_S Integer number of draws per individual, used only when \code{gen_seed >= 0}.
//' @param Ti Optional integer vector with the number of choice situations of
//'   each decision maker (panel likelihood); situations must be sorted by
//'   decision maker. NULL (default): every choice situation is its own unit
//'   (cross-sectional likelihood).
//' @param draw_batch Integer; \code{0} (default) forms each decision maker's
//'   draws in batches sized to a per-thread memory budget, a positive value
//'   caps the number of draws per batch (for tests).
//' @returns n_params x n_params PSD matrix representing the observed information
//'   matrix estimated by the outer product of gradients (same sign convention
//'   as the negated Hessian returned by \code{mxl_hessian_parallel}, so it can
//'   be inverted directly to obtain vcov).
//' @note The BHHH/OPG estimator is only asymptotically equivalent to the
//'   Hessian-based information matrix at the true MLE. In finite samples it can
//'   underestimate standard errors, particularly when the model is mis-specified
//'   or away from the optimum.
//' @examples
//' \donttest{
//' library(data.table)
//' set.seed(42)
//' N <- 50; J <- 3
//' dt <- data.table(id = rep(1:N, each = J), alt = rep(1:J, N))
//' dt[, `:=`(x1 = rnorm(.N), w1 = rnorm(.N))]
//' dt[, choice := 0L]
//' dt[, choice := sample(c(1L, rep(0L, J - 1))), by = id]
//' d <- prepare_mxl_data(dt, "id", "alt", "choice", "x1", "w1")
//' eta <- get_halton_normals(50, d$N, ncol(d$W))
//' theta <- rep(0, ncol(d$X) + ncol(d$W) + nrow(d$alt_mapping) - 1)
//' H <- choicer:::mxl_bhhh_parallel(theta, d$X, d$W, d$alt_idx, d$choice_idx,
//'   d$M, d$weights, eta, rc_dist = rep(0L, ncol(d$W)),
//'   rc_correlation = FALSE, rc_mean = FALSE)
//' dim(H)
//' }
//' @keywords internal
// [[Rcpp::export]]
arma::mat mxl_bhhh_parallel(
    const arma::vec &theta, const arma::mat &X, const arma::mat &W,
    const Rcpp::IntegerVector &alt_idx,
    const Rcpp::IntegerVector &choice_idx,
    const Rcpp::IntegerVector &M, const arma::vec &weights,
    const arma::cube &eta_draws, const arma::uvec &rc_dist,
    const bool rc_correlation = true, const bool rc_mean = false,
    const bool use_asc = true, const bool include_outside_option = false,
    const int gen_seed = -1, const int gen_scramble = 1, const int gen_S = 0,
    const Rcpp::Nullable<Rcpp::IntegerVector> Ti = R_NilValue,
    const int draw_batch = 0) {
  // Inputs, layout and parameters at theta, validated on the primary thread
  const MxlUnitData ud(theta, X, W, alt_idx, choice_idx, M,
                       &weights, eta_draws, rc_dist,
                       rc_correlation, rc_mean, use_asc, include_outside_option,
                       gen_seed, gen_scramble, gen_S, Ti, draw_batch,
                       "mxl_bhhh_parallel");
  const MxlLayout& lay = ud.lay;
  const int n_params = ud.n_params;

  // Global BHHH accumulator
  arma::mat global_bhhh = arma::zeros(n_params, n_params);

#ifdef _OPENMP
#pragma omp parallel
#endif
  {
    // Thread-local accumulator and unit scratch
    arma::mat local_bhhh = arma::zeros(n_params, n_params);
    MxlUnitScratch sc;
    arma::vec s_u; // score of unit u

// Loop over likelihood units in parallel
#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
    for (mxl_off u = 0; u < lay.U; ++u) {
      mxl_unit_load(ud, u, sc, 2);
      mxl_unit_simulate(ud, sc, true);
      mxl_unit_score(ud, sc, s_u);
      local_bhhh += weights[lay.unit_first(u)] * s_u * s_u.t();
    } // end unit loop

#ifdef _OPENMP
#pragma omp critical
#endif
    {
      global_bhhh += local_bhhh;
    }
  } // end parallel region

  // Return PSD information matrix (same sign convention as negated Hessian).
  return global_bhhh;
}

// Per-unit score matrix for the mixed logit model (internal).
//
// Returns the U x n_params matrix whose row u is the weight-free score
// s_u = d ell_u / d theta of likelihood unit u (a decision maker when Ti is
// supplied, a choice situation otherwise) over the beta, mu, L and delta/ASC
// blocks. Same per-unit routine as mxl_bhhh_parallel with the outer-product
// accumulator replaced by a row write; weights are applied on the R side (see
// .assemble_score_vcov in R/classes.R).
// [[Rcpp::export]]
arma::mat mxl_scores_parallel(
    const arma::vec &theta, const arma::mat &X, const arma::mat &W,
    const Rcpp::IntegerVector &alt_idx,
    const Rcpp::IntegerVector &choice_idx,
    const Rcpp::IntegerVector &M,
    const arma::cube &eta_draws, const arma::uvec &rc_dist,
    const bool rc_correlation = true, const bool rc_mean = false,
    const bool use_asc = true, const bool include_outside_option = false,
    const int gen_seed = -1, const int gen_scramble = 1, const int gen_S = 0,
    const Rcpp::Nullable<Rcpp::IntegerVector> Ti = R_NilValue,
    const int draw_batch = 0) {
  // Inputs, layout and parameters at theta, validated on the primary thread
  const MxlUnitData ud(theta, X, W, alt_idx, choice_idx, M,
                       nullptr, eta_draws, rc_dist,
                       rc_correlation, rc_mean, use_asc, include_outside_option,
                       gen_seed, gen_scramble, gen_S, Ti, draw_batch,
                       "mxl_scores_parallel");
  const MxlLayout& lay = ud.lay;
  const int n_params = ud.n_params;

  // Output: one row per likelihood unit (each written by exactly one
  // iteration, so no accumulator or critical section is needed).
  arma::mat scores(lay.U, n_params);

#ifdef _OPENMP
#pragma omp parallel
#endif
  {
    MxlUnitScratch sc;
    arma::vec s_u; // score of unit u

// Loop over likelihood units in parallel
#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
    for (mxl_off u = 0; u < lay.U; ++u) {
      mxl_unit_load(ud, u, sc, 2);
      mxl_unit_simulate(ud, sc, true);
      mxl_unit_score(ud, sc, s_u);
      scores.row(u) = s_u.t();
    } // end unit loop
  } // end parallel region

  return scores;
}

// Conditional tastes of each likelihood unit (internal).
//
// Simulated mean and standard deviation of the random coefficients given the
// unit's observed choices (Revelt & Train 2000; Train 2009, ch. 11), on the
// scale the coefficients enter utility, beta_us = mu_final + Gamma_us (the
// log-normal rows of Gamma exponentiated). With the posterior draw weights
// omega_us of the shared draw loop:
//   mean_u = mu_final + Gamma_u omega_u
//   sd_u   = sqrt(sum_s omega_us (Gamma_us - Gamma_u omega_u)^2)   (elementwise)
// Takes the inputs of mxl_scores_parallel and returns list(mean, sd), each
// K_w x U (U = decision makers when Ti is supplied, choice situations
// otherwise). A unit whose simulated log-likelihood is not finite, which
// takes utilities that overflowed, gets NA.
// [[Rcpp::export]]
Rcpp::List mxl_conditional_tastes_parallel(
    const arma::vec &theta, const arma::mat &X, const arma::mat &W,
    const Rcpp::IntegerVector &alt_idx,
    const Rcpp::IntegerVector &choice_idx,
    const Rcpp::IntegerVector &M,
    const arma::cube &eta_draws, const arma::uvec &rc_dist,
    const bool rc_correlation = true, const bool rc_mean = false,
    const bool use_asc = true, const bool include_outside_option = false,
    const int gen_seed = -1, const int gen_scramble = 1, const int gen_S = 0,
    const Rcpp::Nullable<Rcpp::IntegerVector> Ti = R_NilValue,
    const int draw_batch = 0) {
  // Inputs, layout and parameters at theta, validated on the primary thread
  const MxlUnitData ud(theta, X, W, alt_idx, choice_idx, M,
                       nullptr, eta_draws, rc_dist,
                       rc_correlation, rc_mean, use_asc, include_outside_option,
                       gen_seed, gen_scramble, gen_S, Ti, draw_batch,
                       "mxl_conditional_tastes_parallel");
  const MxlLayout& lay = ud.lay;
  const arma::vec& mu_final = ud.par.mu_final;
  const double na = NA_REAL; // read on the master thread

  // Output: one column per likelihood unit (disjoint writes, no reduction).
  arma::mat taste_mean(W.n_cols, lay.U);
  arma::mat taste_sd(W.n_cols, lay.U);

#ifdef _OPENMP
#pragma omp parallel
#endif
  {
    MxlUnitScratch sc;
    arma::vec gamma_bar; // Gamma_u omega_u: conditional mean of Gamma
    arma::mat dev;       // (Gamma_u - gamma_bar)^2, K_w x S

// Loop over likelihood units in parallel
#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
    for (mxl_off u = 0; u < lay.U; ++u) {
      mxl_unit_load(ud, u, sc, 1, false);
      const double lse = mxl_unit_simulate(ud, sc, false);
      if (!std::isfinite(lse)) {
        taste_mean.col(u).fill(na);
        taste_sd.col(u).fill(na);
        continue;
      }
      gamma_bar = sc.Gamma * sc.omega;
      taste_mean.col(u) = mu_final + gamma_bar;
      dev = sc.Gamma;                 // copied into the reused buffer
      dev.each_col() -= gamma_bar;
      dev %= dev;
      taste_sd.col(u) = arma::sqrt(dev * sc.omega);
    } // end unit loop
  } // end parallel region

  return Rcpp::List::create(Rcpp::Named("mean") = taste_mean,
                            Rcpp::Named("sd") = taste_sd);
}

// ============================================================================
// Prediction: one choice situation at a time
//
// The prediction kernels integrate each choice situation over the population
// taste distribution with its own K_w x S block of draws: situation t reads
// slice t of a given cube or forms Halton block t + 1, in a panel as in the
// cross-section (predictions are unconditional). A kernel validates
// its inputs and lays out the stacked design on the primary thread (the
// alternative codes are read in place, row offsets are 64-bit), forms the base
// utilities once and allocates every thread's buffers there, then loops over
// the situations in parallel without allocating.
//
// Store-mode fits pass no cube either: with gen_scramble = 2
// (MXL_STORE_POINTS) a kernel forms each situation's block on the fly from
// the generator's identity-permutation uniforms, randtoolbox::halton()'s
// points, mapped to normals by R's qnorm(), as get_halton_normals() maps
// them. The draws are the cube's, bit for bit, where randtoolbox and choicer
// are compiled with the same floating-point contraction (see
// HaltonGen::fill_uniforms()); otherwise some uniforms differ in their last
// bit, which moves those draws by at most a few parts in 10^9 (in the tails,
// through qnorm()'s slope). blp() can also keep every situation's
// Gamma = L eta across the contraction's passes, within a byte budget
// (mxl_pred_keep_gamma()), so that only its first pass forms the draws.
//
// Every result is bitwise that of the per-situation Armadillo code this
// replaced, with BLAS libraries whose results do not depend on operand
// addresses, such as the reference BLAS and OpenBLAS: per-situation outputs at
// any thread count (a multithreaded BLAS may split the base product, as
// before), sums over situations at one thread (at more, the situations reach
// the threads in a varying order, as before):
//   * the base utilities are the same full products, X beta, then
//     += W mu_final for a row-aligned W (one dgemv with beta = 1 into
//     X beta), with an alternative-level W's W mu_final and then the ASCs
//     added per row, so no BLAS call changes its shape;
//   * Gamma = L eta and W_t Gamma are the same Armadillo products on
//     matrices of the same shapes (Armadillo picks its BLAS call, or its own
//     code for tiny matrices, from the shapes alone), and the draws reach them
//     through a copy into a thread buffer, as they did;
//   * the draw loops use stable_softmax_n() and max_shifted_lse_n()
//     (Armadillo's operations in Armadillo's order) and element-wise sums;
//   * every statement that adds a product to an accumulator is kept as it
//     was, its Armadillo operator() calls (and their bounds checks)
//     included, so each compiler fuses the same statements into multiply-adds
//     as before: clang contracts each such statement as it parses it, and
//     GCC decides after optimizing, from the code around the statement
//     (g++-16 fuses the same statements here as in the code this replaced,
//     at -O2 and -O3, with and without OpenMP). Written with [], .at() or a
//     raw pointer, such a statement can fuse where it did not, which changes
//     results under GCC with FMA.
// ============================================================================

// An n x 1 R matrix, the form in which RcppArmadillo returns an arma::vec
// (choicer.h includes <RcppArmadillo.h>, which leaves
// RCPP_ARMADILLO_RETURN_COLVEC_AS_VECTOR undefined), not initialized: the
// kernel writes every element.
inline Rcpp::NumericVector mxl_col_result(const mxl_off n) {
  Rcpp::NumericVector x(Rcpp::no_init(static_cast<R_xlen_t>(n)));
  x.attr("dim") = Rcpp::Dimension(static_cast<std::size_t>(n), 1);
  return x;
}

// The gen_scramble value of store-mode draws formed on the fly.
constexpr int MXL_STORE_POINTS = 2;

// Where situation t's K_w x S draws come from: slice t of a cube in memory,
// or Halton block t + 1, formed on the fly. On the fly, the generator's
// uniforms go through its own inverse normal CDF (generate mode, gen_scramble
// 0 or 1), or, for store-mode points (MXL_STORE_POINTS), its
// identity-permutation uniforms go through R's qnorm(). The threads may call
// R::qnorm, as bayes_samplers.h does: Rf_qnorm5 keeps no state, allocates
// nothing and never reaches R's warning or error machinery (its one
// diagnostic, for p outside [0, 1], is an ME_DOMAIN warning that Rmath
// suppresses), and the uniforms lie in (0, 1). Rmath functions that can warn
// (ME_RANGE, ME_PRECISION) would not be safe there.
struct MxlPredDraws {
  int K_w = 0, S = 0;
  bool generate = false;   // formed on the fly
  bool r_qnorm = false;    // store-mode points: identity uniforms, R's qnorm()
  HaltonGen gen;
  const double* cube = nullptr;  // a given cube: every situation's draws

  // Write situation t's draws (column-major K_w x S) to eta.
  void fill(double* eta, const mxl_off t) const {
    const std::size_t n = static_cast<std::size_t>(K_w) * static_cast<std::size_t>(S);
    if (!generate) {
      if (n > 0) {
        std::memcpy(eta, cube + static_cast<std::size_t>(t) * n,
                    n * sizeof(double));
      }
      return;
    }
    const uint64_t n0 = static_cast<uint64_t>(t) * static_cast<uint64_t>(S) + 1;
    if (!r_qnorm) {
      gen.fill_block(eta, n0);
      return;
    }
    gen.fill_uniforms(eta, n0);
    for (std::size_t j = 0; j < n; ++j) eta[j] = R::qnorm(eta[j], 0.0, 1.0, 1, 0);
  }
};

// What the per-situation routine reads, shared by all threads: the stacked
// design and its layout, the parameters, the base utilities and the draws.
// A kernel fills it on the primary thread in the order of its own checks;
// it holds no SEXP. The routine writes to it in one place: in blp()'s first
// pass, each situation's slot of gamma_cache (through the const reference;
// the slots are disjoint).
struct MxlPredData {
  const arma::mat& X;
  const arma::mat& W;               // row-aligned with X, or J x K_w
  const arma::uvec& rc_dist;
  ChoiceLayout lay;
  arma::vec mu_final, delta;        // delta: padded, all J; empty without ASCs
  arma::mat L;
  const double* base = nullptr;     // per row: X beta (+ W mu_final, row-aligned W)
  arma::vec W_mu;                   // alternative-level W: W mu_final per alternative
  MxlPredDraws draws;
  // blp(): every situation's Gamma (K_w x S, slice t at t K_w S), kept across
  // the contraction's passes, in which L and the draws do not change: the
  // first pass writes it, gamma_filled is then set and the later passes read
  // it.
  std::unique_ptr<double[]> gamma_cache;
  bool gamma_filled = false;
  int K_w = 0, S = 0;
  bool use_asc = false, include_outside_option = false, alt_level_W = false;

  MxlPredData(const arma::mat& X_, const arma::mat& W_, const arma::uvec& rc_dist_)
      : X(X_), W(W_), rc_dist(rc_dist_) {}
  MxlPredData(const MxlPredData&) = delete;
  MxlPredData& operator=(const MxlPredData&) = delete;
};

// Checks of a prediction kernel's draws formed on the fly, before its layout.
inline void mxl_pred_check_generate(const int gen_S, const int K_w,
                                    const int gen_scramble) {
  if (gen_S <= 0) Rcpp::stop("gen_S must be positive when gen_seed >= 0");
  if (K_w > HALTON_N_PRIMES) {
    Rcpp::stop("K_w exceeds the primes table size (128); reduce K_w or "
               "extend the primes table.");
  }
  if (gen_scramble < 0 || gen_scramble > MXL_STORE_POINTS) {
    Rcpp::stop("gen_scramble must be 0, 1 or 2 when gen_seed >= 0; got %d.",
               gen_scramble);
  }
}

// The layout of a prediction kernel and its checks, in the order of the code
// this replaced: the situations, the alternative codes and the ASCs' coverage
// of them (choice_layout_build(), which names NA and negative codes, where
// the old cast to an unsigned index reported an arbitrary value), then the
// dimensions of a given cube (cube: no draws formed on the fly), then W's
// rows. Draws formed on the fly check W's rows too (generate mode used to
// stop at an Armadillo bounds error).
inline void mxl_pred_layout(MxlPredData& pd, const Rcpp::IntegerVector& alt_idx,
                            const Rcpp::IntegerVector& M,
                            const arma::cube& eta_draws, const bool cube,
                            const bool use_asc, const arma::vec& delta,
                            const arma::vec* weights) {
  pd.lay = choice_layout_build(pd.X, alt_idx, M, use_asc, delta, weights);
  const ChoiceLayout& lay = pd.lay;
  if (cube) {
    if (static_cast<mxl_off>(eta_draws.n_slices) != lay.N) {
      Rcpp::stop("eta_draws 3rd dimension (%d) does not match N (%d)",
                 eta_draws.n_slices, lay.N);
    }
    if (eta_draws.n_rows != pd.W.n_cols) {
      Rcpp::stop("eta_draws 1st dimension (%d) does not match K_w (%d)",
                 eta_draws.n_rows, pd.W.n_cols);
    }
  }
  if (pd.W.n_rows != pd.X.n_rows && lay.n_rows > 0 &&
      static_cast<mxl_off>(pd.W.n_rows) < lay.J) {
    Rcpp::stop("W must be row-aligned with X (%d rows) or contain one row per "
               "global alternative (at least %d rows); got %d rows.",
               pd.X.n_rows, lay.J, pd.W.n_rows);
  }
}

// The parameters, flags and draw source of a prediction kernel, after its
// checks.
inline void mxl_pred_setup(MxlPredData& pd, const arma::vec& mu_final,
                           const arma::mat& L, const arma::vec& delta,
                           const bool use_asc, const bool include_outside_option,
                           const arma::cube& eta_draws, const int gen_seed,
                           const int gen_scramble, const int gen_S) {
  pd.mu_final = mu_final;
  pd.L = L;
  pd.gamma_filled = false;  // a cache holds the Gamma of this L only
  pd.delta = delta;
  pd.use_asc = use_asc;
  pd.include_outside_option = include_outside_option;
  pd.alt_level_W = pd.W.n_rows != pd.X.n_rows;
  pd.K_w = static_cast<int>(pd.W.n_cols);
  pd.S = gen_seed >= 0 ? gen_S : static_cast<int>(eta_draws.n_cols);
  pd.draws.K_w = pd.K_w;
  pd.draws.S = pd.S;
  pd.draws.generate = gen_seed >= 0;
  if (pd.draws.generate) {
    pd.draws.r_qnorm = gen_scramble == MXL_STORE_POINTS;
    pd.draws.gen = HaltonGen(static_cast<uint64_t>(gen_seed), pd.S, pd.K_w,
                             pd.draws.r_qnorm ? 0 : gen_scramble);
  } else {
    pd.draws.cube = eta_draws.memptr();
  }
}

// The base utilities of every row into `base` (n rows; it may view a
// kernel's output): X beta, then += W mu_final for a row-aligned W, the
// expressions the prediction kernels have always formed before the ASCs. An
// alternative-level W's W mu_final goes to W_mu instead; mxl_pred_load() adds
// it and then the ASCs row by row, in the order they always were.
inline void mxl_pred_base(MxlPredData& pd, arma::vec& base,
                          const arma::vec& beta) {
  base = pd.X * beta;
  if (!pd.alt_level_W) {
    base += pd.W * pd.mu_final;
  } else {
    pd.W_mu = pd.W * pd.mu_final;
  }
  pd.base = base.memptr();
}

// A thread's buffers for the situation in hand, sized once on the primary
// thread for the largest choice set and viewed per situation by Armadillo
// matrices of the exact shape; never resized. aux holds a kernel's own
// per-situation arrays; tid is the thread that uses them.
struct MxlPredScratch {
  std::unique_ptr<double[]> eta, gamma, W_t, WGamma, bu, v, p, aux;
  std::unique_ptr<int[]> map;  // local-to-global alternatives (m + 1)
  int tid = 0;
  std::size_t gamma_reads = 0;  // situations whose Gamma came from blp()'s cache

  MxlPredScratch(const MxlPredData& pd, const std::size_t n_aux, const int tid_)
      : tid(tid_) {
    const std::size_t K_w = pd.K_w, S = pd.S;
    const std::size_t m = static_cast<std::size_t>(pd.lay.max_m);
    eta = mxl_buffer(K_w * S);
    gamma = mxl_buffer(K_w * S);
    W_t = mxl_buffer(m * K_w);
    WGamma = mxl_buffer(m * S);
    bu = mxl_buffer(m);
    v = mxl_buffer(m + 1);
    p = mxl_buffer(m + 1);
    aux = mxl_buffer(n_aux);
    map.reset(new int[m + 1 + 32]);  // padded as mxl_buffer()
  }
};

// Every thread's scratch, allocated on the primary thread, so that running
// out of memory is an R error saying how much each thread needs, not a
// failure inside the parallel region.
inline std::vector<MxlPredScratch> mxl_pred_scratch(const MxlPredData& pd,
                                                    const int n_threads,
                                                    const std::size_t n_aux) {
  std::vector<MxlPredScratch> sc;
  try {
    sc.reserve(n_threads);
    for (int i = 0; i < n_threads; ++i) sc.emplace_back(pd, n_aux, i);
  } catch (const std::bad_alloc&) {
    std::vector<MxlPredScratch>().swap(sc);  // release before reporting
    const double K_w = pd.K_w, S = pd.S, m = static_cast<double>(pd.lay.max_m);
    const double bytes = sizeof(double) * (2.0 * K_w * S + m * (K_w + S + 3.0) +
                                           static_cast<double>(n_aux)) +
                         sizeof(int) * m;
    if (n_threads > 1) {
      Rcpp::stop("Not enough memory for the prediction's working arrays: "
                 "%.2f GB per thread for %d threads (the largest choice "
                 "situation stacks %d alternative rows; S = %d draws). Run "
                 "fewer threads with set_num_threads().", bytes / 1e9,
                 n_threads, pd.lay.max_m, pd.S);
    }
    Rcpp::stop("Not enough memory for the prediction's working arrays: %.2f GB "
               "(the largest choice situation stacks %d alternative rows; "
               "S = %d draws).", bytes / 1e9, pd.lay.max_m, pd.S);
  }
  return sc;
}

// Run body(t, sc) for every situation t, in parallel, each thread with its
// own scratch.
template <typename Body>
inline void mxl_pred_run(const MxlPredData& pd,
                         std::vector<MxlPredScratch>& scratch, Body body) {
  const int n_threads = static_cast<int>(scratch.size());
  const mxl_off N = pd.lay.N;
  if (N == 0) return;
#ifdef _OPENMP
#pragma omp parallel num_threads(n_threads)
#endif
  {
    MxlPredScratch& sc = scratch[mxl_thread_num()];
#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
    for (mxl_off t = 0; t < N; ++t) body(t, sc);
  }
}

// Load situation t into the thread's buffers: its base utilities (with an
// alternative-level W's W mu_final, then the ASCs, added row by row), its
// draws, Gamma = L eta with the log-normal transform (or Gamma copied from
// the cache once blp()'s first pass has filled it), its rows of W and
// W_t Gamma (m x S). Returns its number of rows m.
inline int mxl_pred_load(const MxlPredData& pd, const mxl_off t,
                         MxlPredScratch& sc) {
  const ChoiceLayout& lay = pd.lay;
  const mxl_off r0 = lay.row_off[t];
  const int m = lay.m(t);
  const AltCodes0 alt0{lay.alt + r0};
  double* bu = sc.bu.get();
  for (int a = 0; a < m; ++a) bu[a] = pd.base[r0 + a];
  if (pd.alt_level_W) {
    for (int a = 0; a < m; ++a) bu[a] += pd.W_mu[alt0[a]];
  }
  if (pd.use_asc) {
    for (int a = 0; a < m; ++a) bu[a] += pd.delta[alt0[a]];
  }

  const int K_w = pd.K_w, S = pd.S;
  const std::size_t n_eta = static_cast<std::size_t>(K_w) * static_cast<std::size_t>(S);
  double* slot = pd.gamma_cache
                 ? pd.gamma_cache.get() + static_cast<std::size_t>(t) * n_eta
                 : nullptr;
  arma::mat Gamma(sc.gamma.get(), K_w, S, false, true);
  if (slot && pd.gamma_filled) {
    std::memcpy(sc.gamma.get(), slot, n_eta * sizeof(double));
    ++sc.gamma_reads;
  } else {
    pd.draws.fill(sc.eta.get(), t);
    const arma::mat eta(sc.eta.get(), K_w, S, false, true);
    batch_gamma_draws_into(Gamma, pd.L, eta, pd.rc_dist);
    if (slot) std::memcpy(slot, sc.gamma.get(), n_eta * sizeof(double));
  }

  arma::mat W_t(sc.W_t.get(), m, K_w, false, true);
  for (int k = 0; k < K_w; ++k) {
    for (int a = 0; a < m; ++a) {
      W_t.at(a, k) = pd.alt_level_W ? pd.W.at(alt0[a], k) : pd.W.at(r0 + a, k);
    }
  }
  arma::mat WGamma(sc.WGamma.get(), m, S, false, true);
  WGamma = W_t * Gamma;
  return m;
}

// Per-thread accumulators of a kernel that sums over situations: n doubles
// per scratch set, in padded buffers allocated on the primary thread and not
// initialized (mxl_pred_zero() zeros them); running out of memory is an R
// error saying how much each needs. Each thread adds into its own through an
// Armadillo view, whose operator() keeps the bounds check of the code this
// replaced, and the kernel adds them up in thread order afterwards, which at
// one thread is the single addition the critical section used to make.
using MxlPredAcc = std::vector<std::unique_ptr<double[]>>;

inline MxlPredAcc mxl_pred_accumulators(const std::vector<MxlPredScratch>& scratch,
                                        const std::size_t n, const char* what) {
  const int n_threads = static_cast<int>(scratch.size());
  MxlPredAcc acc;
  try {
    acc.reserve(n_threads);
    for (int i = 0; i < n_threads; ++i) acc.push_back(mxl_buffer(n));
  } catch (const std::bad_alloc&) {
    MxlPredAcc().swap(acc);  // release before reporting
    const double gb = sizeof(double) * static_cast<double>(n) / 1e9;
    if (n_threads > 1) {
      Rcpp::stop("Not enough memory for the %s accumulators: %.2f GB per "
                 "thread for %d threads. Run fewer threads with "
                 "set_num_threads().", what, gb, n_threads);
    }
    Rcpp::stop("Not enough memory for the %s accumulators: %.2f GB.", what, gb);
  }
  return acc;
}

// The count of alternatives the accumulators index: compute_J_total() wraps
// the largest code, 2^31 - 1, plus the outside option to a negative int (it
// used to reach Armadillo's allocation as an enormous size).
inline void mxl_pred_check_alternatives(const int J_total) {
  if (J_total <= 0) {
    Rcpp::stop("alt_idx references alternative %d, which with the outside "
               "option is more alternatives than an int can count.",
               std::numeric_limits<int>::max());
  }
}

// Zero the accumulators, thread i zeroing accumulator i (schedule(static)),
// so that it touches first the memory it will add to.
inline void mxl_pred_zero(MxlPredAcc& acc, const std::size_t n) {
  const int n_threads = static_cast<int>(acc.size());
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(n_threads)
#endif
  for (int i = 0; i < n_threads; ++i) {
    std::fill(acc[i].get(), acc[i].get() + n, 0.0);
  }
}

// Simulated market shares over num_alts alternatives: each situation's
// probabilities, averaged over draws, weighted and added by alternative into
// its thread's accumulator (the outside option in slot 0 when present), then
// the accumulators added up, over the sum of the weights. The additions into
// local_shares are the statements of the code this replaced (see the top of
// this section).
inline arma::vec mxl_pred_shares(const MxlPredData& pd, const arma::vec& weights,
                                 const double weight_sum, const int num_alts,
                                 std::vector<MxlPredScratch>& scratch,
                                 MxlPredAcc& acc) {
  if (acc.size() != scratch.size()) {
    Rcpp::stop("Internal error: one accumulator per scratch set expected.");
  }
  mxl_pred_zero(acc, num_alts);
  const int S = pd.S;
  const bool include_outside_option = pd.include_outside_option;
  const int o = include_outside_option ? 1 : 0;  // outside option: slot 0

  mxl_pred_run(pd, scratch, [&](const mxl_off t, MxlPredScratch& sc) {
    arma::vec local_shares(acc[sc.tid].get(), num_alts, false, true);
    const int m_i = mxl_pred_load(pd, t, sc);
    const int num_choices = m_i + o;
    const double w_i = weights[t];
    const AltCodes0 alt0{pd.lay.alt + pd.lay.row_off[t]};
    const double* bu = sc.bu.get();
    double* v = sc.v.get();
    double* p = sc.p.get();

    // Accumulate probabilities over draws
    arma::vec P_bar_i(sc.aux.get(), num_choices, false, true);
    P_bar_i.zeros();
    for (int s = 0; s < S; ++s) {
      const double* wg = sc.WGamma.get() + static_cast<std::size_t>(s) * m_i;
      if (o) v[0] = 0.0;
      for (int a = 0; a < m_i; ++a) v[o + a] = bu[a] + wg[a];
      stable_softmax_n(v, p, num_choices);
      for (int j = 0; j < num_choices; ++j) P_bar_i[j] += p[j];
    }
    P_bar_i /= static_cast<double>(S);

    // Accumulate shares by alternative
    if (include_outside_option) {
      local_shares(0) += w_i * P_bar_i(0);
    }
    for (int a = 0; a < m_i; ++a) {
      if (include_outside_option) {
        local_shares(alt0[a] + 1) += w_i * P_bar_i(a + 1);
      } else {
        local_shares(alt0[a]) += w_i * P_bar_i(a);
      }
    }
  });

  arma::vec global_shares = arma::zeros(num_alts);
  for (const std::unique_ptr<double[]>& a : acc) {
    global_shares += arma::vec(a.get(), num_alts, false, true);
  }
  return global_shares / weight_sum;
}

// The values x_k of the perturbed variable over situation t's choice set, in
// the slots of its probabilities (the outside option's, 0, first): column
// var_idx of its rows of W (random coefficient) or of X (fixed).
inline void mxl_pred_x_k(const MxlPredData& pd, const mxl_off t,
                         const MxlPredScratch& sc, const int m_i,
                         const int var_idx, const bool is_random_coef,
                         arma::vec& x_k_i) {
  const int o = pd.include_outside_option ? 1 : 0;
  const mxl_off r0 = pd.lay.row_off[t];
  const arma::mat W_t(sc.W_t.get(), m_i, pd.K_w, false, true);
  x_k_i.zeros();
  for (int a = 0; a < m_i; ++a) {
    x_k_i.at(o + a) = is_random_coef ? W_t.at(a, var_idx)
                                     : pd.X.at(r0 + a, var_idx);
  }
}

// ============================================================================
// Mixed Logit: Share Prediction and BLP Contraction
// ============================================================================

//' Per-observation simulated choice probabilities for Mixed Logit
//'
//' Returns the simulated choice probability for each (individual, alternative)
//' row of `X`, averaged over the supplied Halton draws. Mirrors `mnl_predict`.
//'
//' @param theta parameter vector (beta, \[mu\], L, delta)
//' @param X design matrix for fixed coefficients; sum(M_i) x K_x
//' @param W design matrix for random coefficients; sum(M_i) x K_w or J x K_w
//' @param alt_idx sum(M) x 1 vector with indices of alternatives; 1-based indexing
//' @param M N x 1 vector with number of alternatives for each individual
//' @param eta_draws Array with draws; K_w x S x N
//' @param rc_dist K_w vector indicating distribution (0=normal, 1=log-normal)
//' @param rc_correlation whether random coefficients are correlated
//' @param rc_mean whether mu parameters are estimated
//' @param use_asc whether ASCs are included
//' @param include_outside_option whether the outside option is present
//' @param gen_seed Integer. \code{< 0} (default) uses the materialized
//'   \code{eta_draws} cube; \code{>= 0} forms the draws on the fly (see
//'   \code{gen_scramble}), with this master seed for the digit permutations of
//'   \code{gen_scramble = 1}.
//' @param gen_scramble Integer mode of the draws formed on the fly: \code{0} =
//'   identity permutations (plain Halton) through the generator's inverse normal
//'   CDF, \code{1} = seeded position-wise digit permutations, \code{2} = identity
//'   permutations through R's \code{qnorm()}, the store-mode draws of
//'   \code{\link{get_halton_normals}} formed without its cube (bit for bit where
//'   randtoolbox and choicer are compiled with the same floating-point
//'   contraction).
//' @param gen_S Integer number of draws per individual, used only when \code{gen_seed >= 0}.
//' @returns List with `choice_prob` (length sum(M)), `utility` (length sum(M),
//'   simulated mean of the deterministic + W*gamma component), and, when
//'   `include_outside_option = TRUE`, `choice_prob_outside` (length N).
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List mxl_predict(
    const arma::vec& theta,
    const arma::mat& X,
    const arma::mat& W,
    const Rcpp::IntegerVector& alt_idx,
    const Rcpp::IntegerVector& M,
    const arma::cube& eta_draws,
    const arma::uvec& rc_dist,
    const bool rc_correlation = true,
    const bool rc_mean = false,
    const bool use_asc = true,
    const bool include_outside_option = false,
    const int gen_seed = -1, const int gen_scramble = 1, const int gen_S = 0
) {
  // Parse theta into parameter blocks (shared helper; validates theta), then
  // the inputs and their layout, on the primary thread
  const MxlParams par = parse_mxl_theta(theta, X.n_cols, W.n_cols, rc_dist,
                                        rc_correlation, rc_mean, use_asc,
                                        include_outside_option);
  const bool generate = gen_seed >= 0;
  if (generate) mxl_pred_check_generate(gen_S, W.n_cols, gen_scramble);
  MxlPredData pd(X, W, rc_dist);
  mxl_pred_layout(pd, alt_idx, M, eta_draws, !generate, use_asc, par.delta,
                  nullptr);
  mxl_pred_setup(pd, par.mu_final, par.L, par.delta, use_asc,
                 include_outside_option, eta_draws, gen_seed, gen_scramble,
                 gen_S);
  const ChoiceLayout& lay = pd.lay;
  const std::size_t max_m = static_cast<std::size_t>(lay.max_m);
  std::vector<MxlPredScratch> scratch =
      mxl_pred_scratch(pd, mxl_team_threads(lay.N), 2 * max_m);

  // Outputs, each written once by the situation that owns its slots. The
  // base utilities are formed in choice_prob, which a situation reads into
  // its buffers before writing its probabilities there.
  Rcpp::NumericVector choice_prob = mxl_col_result(lay.n_rows);
  Rcpp::NumericVector utility = mxl_col_result(lay.n_rows);
  Rcpp::NumericVector choice_prob_outside;
  if (include_outside_option) choice_prob_outside = mxl_col_result(lay.N);
  arma::vec base(choice_prob.begin(), lay.n_rows, false, true);
  mxl_pred_base(pd, base, par.beta);
  double* prob = choice_prob.begin();
  double* util = utility.begin();
  double* prob_outside = include_outside_option ? choice_prob_outside.begin()
                                                : nullptr;
  const int S = pd.S;
  const double S_d = static_cast<double>(S);
  const int o = include_outside_option ? 1 : 0;  // outside option: slot 0

  mxl_pred_run(pd, scratch, [&](const mxl_off t, MxlPredScratch& sc) {
    const int m = mxl_pred_load(pd, t, sc);
    const double* bu = sc.bu.get();
    double* v = sc.v.get();
    double* p = sc.p.get();
    double* P_inside_avg = sc.aux.get();          // m: averaged over draws
    double* util_inside_avg = P_inside_avg + max_m;
    for (int a = 0; a < m; ++a) {
      P_inside_avg[a] = 0.0;
      util_inside_avg[a] = 0.0;
    }
    double P_outside_avg = 0.0;

    for (int s = 0; s < S; ++s) {
      const double* wg = sc.WGamma.get() + static_cast<std::size_t>(s) * m;
      if (o) v[0] = 0.0;
      for (int a = 0; a < m; ++a) {
        v[o + a] = bu[a] + wg[a];          // inside utility at draw s
        util_inside_avg[a] += v[o + a];
      }
      stable_softmax_n(v, p, m + o);       // shifts v in place
      if (o) P_outside_avg += p[0];
      for (int a = 0; a < m; ++a) P_inside_avg[a] += p[o + a];
    }

    // Average over draws; disjoint writes by situation — no race
    const mxl_off r0 = lay.row_off[t];
    for (int a = 0; a < m; ++a) {
      prob[r0 + a] = P_inside_avg[a] / S_d;
      util[r0 + a] = util_inside_avg[a] / S_d;
    }
    if (o) prob_outside[t] = P_outside_avg / S_d;
  });

  Rcpp::List out;
  out["choice_prob"] = choice_prob;
  out["utility"] = utility;
  if (include_outside_option) {
    out["choice_prob_outside"] = choice_prob_outside;
  }
  return out;
}

//' Simulated expected logsum (inclusive value) for Mixed Logit
//'
//' Computes the simulated expected logsum (expected maximum utility, up to an
//' additive constant) for each choice situation:
//' \deqn{logsum_i = (1/S) \sum_s \log \sum_j \exp(V_{ij}^s),}
//' where the inner sum runs over individual i's alternatives and includes the
//' outside option's \eqn{\exp(0)} term when `include_outside_option = TRUE`.
//' The log-sum-exp must be averaged *across draws*: applying log-sum-exp to
//' the draw-averaged utilities returned by `mxl_predict` understates the
//' expectation because log-sum-exp is convex (Jensen's inequality).
//'
//' @param theta parameter vector (beta, \[mu\], L, delta)
//' @param X design matrix for fixed coefficients; sum(M_i) x K_x
//' @param W design matrix for random coefficients; sum(M_i) x K_w or J x K_w
//' @param alt_idx sum(M) x 1 vector with indices of alternatives; 1-based indexing
//' @param M N x 1 vector with number of alternatives for each individual
//' @param eta_draws Array with draws; K_w x S x N
//' @param rc_dist K_w vector indicating distribution (0=normal, 1=log-normal)
//' @param rc_correlation whether random coefficients are correlated
//' @param rc_mean whether mu parameters are estimated
//' @param use_asc whether ASCs are included
//' @param include_outside_option whether the outside option is present
//' @param gen_seed Integer. \code{< 0} (default) uses the materialized
//'   \code{eta_draws} cube; \code{>= 0} forms the draws on the fly (see
//'   \code{gen_scramble}), with this master seed for the digit permutations of
//'   \code{gen_scramble = 1}.
//' @param gen_scramble Integer mode of the draws formed on the fly: \code{0} =
//'   identity permutations (plain Halton) through the generator's inverse normal
//'   CDF, \code{1} = seeded position-wise digit permutations, \code{2} = identity
//'   permutations through R's \code{qnorm()}, the store-mode draws of
//'   \code{\link{get_halton_normals}} formed without its cube (bit for bit where
//'   randtoolbox and choicer are compiled with the same floating-point
//'   contraction).
//' @param gen_S Integer number of draws per individual, used only when \code{gen_seed >= 0}.
//' @returns Vector of length N with the simulated expected logsum per choice
//'   situation.
//' @note For log-normal random coefficients (rc_dist=1) with rc_mean=TRUE,
//'   the distribution is a shifted log-normal: beta_k = exp(mu_k) + exp(L_k * eta),
//'   where exp(mu_k) shifts the location and exp(L_k * eta) ~ LogNormal(0, sigma_k^2).
//'   This differs from the textbook parameterization exp(mu_k + L_k * eta).
//' @examples
//' \donttest{
//' library(data.table)
//' set.seed(42)
//' N <- 50; J <- 3
//' dt <- data.table(id = rep(1:N, each = J), alt = rep(1:J, N))
//' dt[, `:=`(x1 = rnorm(.N), w1 = rnorm(.N))]
//' dt[, choice := 0L]
//' dt[, choice := sample(c(1L, rep(0L, J - 1))), by = id]
//' d <- prepare_mxl_data(dt, "id", "alt", "choice", "x1", "w1")
//' eta <- get_halton_normals(50, d$N, ncol(d$W))
//' fit <- run_mxlogit(input_data = d, eta_draws = eta)
//' ls <- choicer:::mxl_logsum(coef(fit), d$X, d$W, d$alt_idx, d$M, eta,
//'   rc_dist = rep(0L, ncol(d$W)), rc_correlation = FALSE, rc_mean = FALSE)
//' head(ls)
//' }
//' @keywords internal
// [[Rcpp::export]]
Rcpp::NumericVector mxl_logsum(const arma::vec &theta, const arma::mat &X, const arma::mat &W,
                     const Rcpp::IntegerVector &alt_idx, const Rcpp::IntegerVector &M,
                     const arma::cube &eta_draws, const arma::uvec &rc_dist,
                     const bool rc_correlation = true, const bool rc_mean = false,
                     const bool use_asc = true, const bool include_outside_option = false,
                     const int gen_seed = -1, const int gen_scramble = 1, const int gen_S = 0) {
  // Parse theta into parameter blocks (shared helper; validates theta), then
  // the inputs and their layout, on the primary thread
  const MxlParams par = parse_mxl_theta(theta, X.n_cols, W.n_cols, rc_dist,
                                        rc_correlation, rc_mean, use_asc,
                                        include_outside_option);
  const bool generate = gen_seed >= 0;
  if (generate) mxl_pred_check_generate(gen_S, W.n_cols, gen_scramble);
  MxlPredData pd(X, W, rc_dist);
  mxl_pred_layout(pd, alt_idx, M, eta_draws, !generate, use_asc, par.delta,
                  nullptr);
  mxl_pred_setup(pd, par.mu_final, par.L, par.delta, use_asc,
                 include_outside_option, eta_draws, gen_seed, gen_scramble,
                 gen_S);
  const ChoiceLayout& lay = pd.lay;

  std::vector<MxlPredScratch> scratch =
      mxl_pred_scratch(pd, mxl_team_threads(lay.N), 0);
  arma::vec base(lay.n_rows, arma::fill::none);
  mxl_pred_base(pd, base, par.beta);
  const int S = pd.S;
  const int o = include_outside_option ? 1 : 0;  // outside option: slot 0

  // Output (each situation writes only its own slot)
  Rcpp::NumericVector logsum = mxl_col_result(lay.N);
  double* ls = logsum.begin();

  mxl_pred_run(pd, scratch, [&](const mxl_off t, MxlPredScratch& sc) {
    const int m = mxl_pred_load(pd, t, sc);
    const double* bu = sc.bu.get();
    double* v = sc.v.get();

    // Accumulate the per-draw log-sum-exp (NOT the logsum of averaged
    // utilities; see the Jensen note in the docs above).
    double logsum_acc = 0.0;
    for (int s = 0; s < S; ++s) {
      const double* wg = sc.WGamma.get() + static_cast<std::size_t>(s) * m;
      if (o) v[0] = 0.0;  // outside option fixed at 0
      for (int a = 0; a < m; ++a) v[o + a] = bu[a] + wg[a];
      logsum_acc += max_shifted_lse_n(v, m + o);
    }

    // Average over draws; disjoint write by situation — no race
    ls[t] = logsum_acc / static_cast<double>(S);
  });

  return logsum;
}

//' Predicted aggregate market shares for Mixed Logit
//'
//' Parses `theta` using the standard parameter ordering and returns the
//' simulated weighted-average market shares.
//'
//' @param theta parameter vector (beta, \[mu\], L, delta)
//' @param X design matrix for fixed coefficients; sum(M_i) x K_x
//' @param W design matrix for random coefficients; sum(M_i) x K_w or J x K_w
//' @param alt_idx sum(M) x 1 vector with indices of alternatives; 1-based indexing
//' @param M N x 1 vector with number of alternatives for each individual
//' @param weights N x 1 vector with weights for each observation
//' @param eta_draws Array with draws; K_w x S x N
//' @param rc_dist K_w vector indicating distribution (0=normal, 1=log-normal)
//' @param rc_correlation whether random coefficients are correlated
//' @param rc_mean whether mu parameters are estimated
//' @param use_asc whether ASCs are included
//' @param include_outside_option whether outside option is included
//' @param gen_seed Integer. \code{< 0} (default) uses the materialized
//'   \code{eta_draws} cube; \code{>= 0} forms the draws on the fly (see
//'   \code{gen_scramble}), with this master seed for the digit permutations of
//'   \code{gen_scramble = 1}.
//' @param gen_scramble Integer mode of the draws formed on the fly: \code{0} =
//'   identity permutations (plain Halton) through the generator's inverse normal
//'   CDF, \code{1} = seeded position-wise digit permutations, \code{2} = identity
//'   permutations through R's \code{qnorm()}, the store-mode draws of
//'   \code{\link{get_halton_normals}} formed without its cube (bit for bit where
//'   randtoolbox and choicer are compiled with the same floating-point
//'   contraction).
//' @param gen_S Integer number of draws per individual, used only when \code{gen_seed >= 0}.
//' @returns Vector of length J (or J+1 with outside option) of predicted shares.
//' @keywords internal
// [[Rcpp::export]]
arma::vec mxl_predict_shares(
    const arma::vec& theta,
    const arma::mat& X,
    const arma::mat& W,
    const Rcpp::IntegerVector& alt_idx,
    const Rcpp::IntegerVector& M,
    const arma::vec& weights,
    const arma::cube& eta_draws,
    const arma::uvec& rc_dist,
    const bool rc_correlation = true,
    const bool rc_mean = false,
    const bool use_asc = true,
    const bool include_outside_option = false,
    const int gen_seed = -1, const int gen_scramble = 1, const int gen_S = 0
) {
  // Parse theta into parameter blocks (shared helper; validates theta), then
  // the inputs and their layout, on the primary thread
  const MxlParams par = parse_mxl_theta(theta, X.n_cols, W.n_cols, rc_dist,
                                        rc_correlation, rc_mean, use_asc,
                                        include_outside_option);
  const bool generate = gen_seed >= 0;
  if (generate) mxl_pred_check_generate(gen_S, W.n_cols, gen_scramble);
  MxlPredData pd(X, W, rc_dist);
  mxl_pred_layout(pd, alt_idx, M, eta_draws, !generate, use_asc, par.delta,
                  &weights);
  mxl_pred_setup(pd, par.mu_final, par.L, par.delta, use_asc,
                 include_outside_option, eta_draws, gen_seed, gen_scramble,
                 gen_S);
  const ChoiceLayout& lay = pd.lay;

  // Number of alternatives for the output
  const int J_inside = compute_J_inside(use_asc, par.delta, lay);
  const int num_alts = compute_J_total(J_inside, include_outside_option);
  const double weight_sum = shares_denominator(weights);
  mxl_pred_check_alternatives(num_alts);

  const int n_threads = mxl_team_threads(lay.N);
  std::vector<MxlPredScratch> scratch = mxl_pred_scratch(
      pd, n_threads, static_cast<std::size_t>(lay.max_m) + 1);
  MxlPredAcc acc = mxl_pred_accumulators(scratch, num_alts, "shares");
  arma::vec base(lay.n_rows, arma::fill::none);
  mxl_pred_base(pd, base, par.beta);
  return mxl_pred_shares(pd, weights, weight_sum, num_alts, scratch, acc);
}

//' Diversion ratios for Mixed Logit (simulated, derivative-based)
//'
//' Computes the matrix of attribute-based diversion ratios for a fitted
//' Mixed Logit model. DR(k, j) is the fraction of demand lost by alternative
//' `j` that is captured by alternative `k` when a marginal change in
//' alternative j's `elast_var` attribute reduces s_j.
//'
//' In MNL the per-draw realized coefficient is a constant, so it cancels in
//' the ratio and the result is independent of the variable chosen. In MXL,
//' the realized coefficient \eqn{\beta_{ik}^s} varies across individuals
//' and draws, so the diversion ratio depends on which attribute is perturbed.
//' For a variable with a fixed coefficient the dependence again vanishes
//' (the constant cancels); for a random-coefficient variable it does not.
//'
//' @param theta parameter vector (beta, \[mu\], L, delta)
//' @param X design matrix for fixed coefficients; sum(M_i) x K_x
//' @param W design matrix for random coefficients; sum(M_i) x K_w or J x K_w
//' @param alt_idx sum(M) x 1 vector with indices of alternatives; 1-based indexing
//' @param M N x 1 vector with number of alternatives for each individual
//' @param weights N x 1 vector with weights for each observation
//' @param eta_draws Array with draws; K_w x S x N
//' @param rc_dist K_w vector indicating distribution (0=normal, 1=log-normal)
//' @param elast_var_idx 1-based index of the perturbed variable
//' @param is_random_coef TRUE if the variable is in W (random coef), FALSE if in X (fixed)
//' @param rc_correlation whether random coefficients are correlated
//' @param rc_mean whether mu parameters are estimated
//' @param use_asc whether ASCs are included
//' @param include_outside_option whether outside option is included
//' @param gen_seed Integer. \code{< 0} (default) uses the materialized
//'   \code{eta_draws} cube; \code{>= 0} forms the draws on the fly (see
//'   \code{gen_scramble}), with this master seed for the digit permutations of
//'   \code{gen_scramble = 1}.
//' @param gen_scramble Integer mode of the draws formed on the fly: \code{0} =
//'   identity permutations (plain Halton) through the generator's inverse normal
//'   CDF, \code{1} = seeded position-wise digit permutations, \code{2} = identity
//'   permutations through R's \code{qnorm()}, the store-mode draws of
//'   \code{\link{get_halton_normals}} formed without its cube (bit for bit where
//'   randtoolbox and choicer are compiled with the same floating-point
//'   contraction).
//' @param gen_S Integer number of draws per individual, used only when \code{gen_seed >= 0}.
//' @returns J x J (or (J+1) x (J+1)) matrix of diversion ratios with zero diagonal.
//' @keywords internal
// [[Rcpp::export]]
arma::mat mxl_diversion_ratios_parallel(
    const arma::vec& theta,
    const arma::mat& X,
    const arma::mat& W,
    const Rcpp::IntegerVector& alt_idx,
    const Rcpp::IntegerVector& M,
    const arma::vec& weights,
    const arma::cube& eta_draws,
    const arma::uvec& rc_dist,
    const int elast_var_idx,
    const bool is_random_coef,
    const bool rc_correlation = true,
    const bool rc_mean = false,
    const bool use_asc = true,
    const bool include_outside_option = false,
    const int gen_seed = -1, const int gen_scramble = 1, const int gen_S = 0
) {
  // Attribute-based diversion ratio (simulated):
  //   DR(k, j) = E_i[ w_i * (1/S) sum_s beta_{ik}^s * P_ij(s) * P_ik(s) ]
  //           /  E_i[ w_i * (1/S) sum_s beta_{ik}^s * P_ij(s) * (1 - P_ij(s)) ]
  // where beta_{ik}^s is the realized coefficient on the perturbed variable
  // for individual i and draw s. For a fixed coefficient beta_{ik}^s is a
  // constant scalar and cancels (MNL property). For a random coefficient it
  // varies and does not cancel. Cross-products P_ij(s)*P_ik(s) must be
  // accumulated INSIDE the draw loop; averaging across draws first and
  // multiplying later is biased. See docs/mixed_logit_math.md.

  // Basic dimensions
  const int K_x = X.n_cols;
  const int K_w = W.n_cols;

  // Validate the perturbed variable index. Catch the empty-block cases
  // first (K_x=0 with is_random_coef=FALSE, or K_w=0 with is_random_coef=TRUE)
  // with an actionable message before the index-range check.
  if (!is_random_coef && K_x == 0) {
    Rcpp::stop("Cannot compute diversion ratios w.r.t. a fixed-coefficient "
               "variable: the model has no fixed coefficients (K_x = 0). "
               "Did you mean is_random_coef = TRUE?");
  }
  if (is_random_coef && K_w == 0) {
    Rcpp::stop("Cannot compute diversion ratios w.r.t. a random-coefficient "
               "variable: the model has no random coefficients (K_w = 0).");
  }
  const int var_idx = elast_var_idx - 1;
  if (is_random_coef) {
    if (var_idx < 0 || var_idx >= K_w) {
      Rcpp::stop("elast_var_idx (%d) is out of bounds for W matrix (K_w=%d).",
                 elast_var_idx, K_w);
    }
  } else {
    if (var_idx < 0 || var_idx >= K_x) {
      Rcpp::stop("elast_var_idx (%d) is out of bounds for X matrix (K_x=%d).",
                 elast_var_idx, K_x);
    }
  }

  // Parse theta into parameter blocks (shared helper; validates theta), then
  // the inputs and their layout, on the primary thread
  const MxlParams par = parse_mxl_theta(theta, K_x, K_w, rc_dist,
                                        rc_correlation, rc_mean, use_asc,
                                        include_outside_option);
  const bool generate = gen_seed >= 0;
  if (generate) mxl_pred_check_generate(gen_S, K_w, gen_scramble);
  MxlPredData pd(X, W, rc_dist);
  mxl_pred_layout(pd, alt_idx, M, eta_draws, !generate, use_asc, par.delta,
                  &weights);
  mxl_pred_setup(pd, par.mu_final, par.L, par.delta, use_asc,
                 include_outside_option, eta_draws, gen_seed, gen_scramble,
                 gen_S);
  const ChoiceLayout& lay = pd.lay;
  const double beta_k = is_random_coef ? 0.0 : par.beta(var_idx);
  const arma::vec& mu_final = par.mu_final;

  // Total alternatives for output matrix
  const int J_inside = compute_J_inside(use_asc, par.delta, lay);
  const int J_total = compute_J_total(J_inside, include_outside_option);
  mxl_pred_check_alternatives(J_total);

  // Every thread's buffers (aux: a situation's cross-products and
  // denominators, sized for the largest choice set) and accumulators, then
  // the base utilities
  const int n_threads = mxl_team_threads(lay.N);
  const std::size_t n_max = static_cast<std::size_t>(lay.max_m) + 1;
  std::vector<MxlPredScratch> scratch =
      mxl_pred_scratch(pd, n_threads, n_max * n_max + n_max);
  const std::size_t J_n = static_cast<std::size_t>(J_total);
  MxlPredAcc numerators =
      mxl_pred_accumulators(scratch, J_n * J_n, "diversion-ratio");
  MxlPredAcc denominators =
      mxl_pred_accumulators(scratch, J_n, "diversion-ratio");
  arma::vec base(lay.n_rows, arma::fill::none);
  mxl_pred_base(pd, base, par.beta);

  mxl_pred_zero(numerators, J_n * J_n);
  mxl_pred_zero(denominators, J_n);
  const int S = pd.S;
  const int o = include_outside_option ? 1 : 0;  // outside option: slot 0

  // The accumulations below are the statements of the code this replaced, on
  // Armadillo matrices over the threads' buffers (see mxl_pred_shares()).
  mxl_pred_run(pd, scratch, [&](const mxl_off t, MxlPredScratch& sc) {
    arma::mat local_numerator(numerators[sc.tid].get(), J_total, J_total,
                              false, true);
    arma::vec local_denominator(denominators[sc.tid].get(), J_total, false,
                                true);
    const int m_i = mxl_pred_load(pd, t, sc);
    const int num_choices = m_i + o;
    const double w_i = weights[t];
    const double* bu = sc.bu.get();
    double* v = sc.v.get();
    const arma::mat Gamma_final(sc.gamma.get(), pd.K_w, S, false, true);

    // Map local indices to global alternative indices
    int* global_j_map = sc.map.get();
    fill_global_alt_map(global_j_map, lay.alt + lay.row_off[t], m_i,
                        include_outside_option);

    // Per-individual accumulators (sum across draws, divided by S below)
    arma::mat ind_num(sc.aux.get(), num_choices, num_choices, false, true);
    arma::vec ind_den(sc.aux.get() + n_max * n_max, num_choices, false, true);
    ind_num.zeros();
    ind_den.zeros();
    arma::vec P_s(sc.p.get(), num_choices, false, true);

    // Loop over draws — cross-products MUST be accumulated INSIDE this loop
    for (int s = 0; s < S; ++s) {
      // Realized coefficient on the perturbed variable for this (i, s).
      // For a fixed coef this is the constant beta_k; for a random coef
      // it is mu_final(var_idx) + Gamma_final(var_idx, s), already
      // transformed (exp(.) applied upstream when rc_dist == 1).
      const double beta_k_eff = is_random_coef
          ? (mu_final(var_idx) + Gamma_final(var_idx, s))
          : beta_k;

      const double* wg = sc.WGamma.get() + static_cast<std::size_t>(s) * m_i;
      if (o) v[0] = 0.0;
      for (int a = 0; a < m_i; ++a) v[o + a] = bu[a] + wg[a];
      stable_softmax_n(v, P_s.memptr(), num_choices);

      // Accumulate cross-products weighted by beta_k_eff inside the draw loop
      for (int j_local = 0; j_local < num_choices; ++j_local) {
        const double P_j = P_s(j_local);
        ind_den(j_local) += beta_k_eff * P_j * (1.0 - P_j);
        for (int k_local = 0; k_local < num_choices; ++k_local) {
          if (k_local == j_local) continue;
          ind_num(k_local, j_local) += beta_k_eff * P_j * P_s(k_local);
        }
      }
    }  // end s loop

    // Average over draws
    const double S_d = static_cast<double>(S);
    ind_num /= S_d;
    ind_den /= S_d;

    // Scatter individual contribution into thread-local globals
    for (int j_local = 0; j_local < num_choices; ++j_local) {
      const int global_j = global_j_map[j_local];
      local_denominator(global_j) += w_i * ind_den(j_local);
      for (int k_local = 0; k_local < num_choices; ++k_local) {
        if (k_local == j_local) continue;
        const int global_k = global_j_map[k_local];
        local_numerator(global_k, global_j) += w_i * ind_num(k_local, j_local);
      }
    }
  });

  // Global accumulators: the threads' added up in thread order
  arma::mat global_numerator = arma::zeros(J_total, J_total);
  arma::vec global_denominator = arma::zeros(J_total);
  for (int i = 0; i < n_threads; ++i) {
    global_numerator += arma::mat(numerators[i].get(), J_total, J_total, false,
                                  true);
    global_denominator += arma::vec(denominators[i].get(), J_total, false, true);
  }
  // Free the threads' accumulators and buffers before the ratios' matrix
  MxlPredAcc().swap(numerators);
  MxlPredAcc().swap(denominators);
  std::vector<MxlPredScratch>().swap(scratch);

  // Final ratios with numerical guard (denominator can be negative when
  // beta_k_eff is negative, e.g. price; check magnitude, not sign)
  arma::mat DR = arma::zeros(J_total, J_total);
  for (int j = 0; j < J_total; ++j) {
    if (std::abs(global_denominator(j)) > 1e-15) {
      for (int k = 0; k < J_total; ++k) {
        if (k != j) {
          DR(k, j) = global_numerator(k, j) / global_denominator(j);
        }
      }
    }
  }

  return DR;
}

// Allocate blp()'s Gamma cache (MxlPredData::gamma_cache) when the draws are
// formed on the fly and all situations' Gamma together take at most
// cache_bytes; without the room, or if the allocation fails, the passes form
// the draws again, as without a budget. The cache holds the values each pass
// would compute, so the shares are the same either way. (Where memory is
// committed only when touched, a failure shows in the first pass instead.)
inline void mxl_pred_keep_gamma(MxlPredData& pd, const double cache_bytes) {
  if (!pd.draws.generate || pd.lay.N <= 0) return;
  const std::size_t per = static_cast<std::size_t>(pd.K_w) * static_cast<std::size_t>(pd.S);
  const std::size_t N = static_cast<std::size_t>(pd.lay.N);
  if (per == 0 || N > std::numeric_limits<std::size_t>::max() / sizeof(double) / per) {
    return;
  }
  const std::size_t n = N * per;
  if (!(static_cast<double>(sizeof(double)) * static_cast<double>(n) <= cache_bytes)) {
    return;
  }
  try {
    pd.gamma_cache.reset(new double[n]);
  } catch (const std::bad_alloc&) {
    pd.gamma_cache.reset();
  }
}

// What blp()'s contraction did with its draws: whether it allocated the
// cache, how many evaluations of the shares it made, and how many situations
// those evaluations read from the cache (all N in each after the first).
struct MxlBlpKeep {
  bool kept = false;
  double passes = 0, reads = 0;
};

// The BLP contraction behind mxl_blp_contraction() and
// mxl_blp_contraction_cached(), which keeps the draws across passes within
// cache_bytes (0: never; nor with max_iter <= 0, which leaves no later pass);
// *keep, when given, reports what it did.
static arma::vec mxl_blp_run(
    const arma::vec& delta,
    const arma::vec& target_shares,
    const arma::mat& X,
    const arma::mat& W,
    const arma::vec& beta,
    const arma::vec& mu,
    const arma::vec& L_params,
    const Rcpp::IntegerVector& alt_idx,
    const Rcpp::IntegerVector& M,
    const arma::vec& weights,
    const arma::cube& eta_draws,
    const arma::uvec& rc_dist,
    const bool rc_correlation,
    const bool rc_mean,
    const bool include_outside_option,
    const double tol,
    const int max_iter,
    const int gen_seed,
    const int gen_scramble,
    const int gen_S,
    const double cache_bytes,
    MxlBlpKeep* keep
) {
  const int K_w = W.n_cols;

  // delta is harmonized to cover every referenced alternative below, so the
  // ASC-coverage check is skipped here (use_asc = false, empty delta).
  check_rc_dist_length(rc_dist, K_w);
  const bool generate = gen_seed >= 0;
  if (generate) mxl_pred_check_generate(gen_S, K_w, gen_scramble);
  MxlPredData pd(X, W, rc_dist);
  mxl_pred_layout(pd, alt_idx, M, eta_draws, !generate, /*use_asc=*/false,
                  arma::vec(), &weights);

  // Build L matrix
  arma::mat L = build_L_mat(L_params, K_w, rc_correlation);

  // Transform mu for log-normal coefficients
  arma::vec mu_final = mu;
  if (rc_mean) {
    for (int k = 0; k < K_w; ++k) {
      if (rc_dist(k) == 1) {  // log-normal
        mu_final(k) = std::exp(mu(k));
      }
    }
  }

  // Deduce number of alternatives from data: the largest code
  const int J_inside = compute_J_inside(false, arma::vec(), pd.lay);      // inside options
  const int num_alts = compute_J_total(J_inside, include_outside_option); // total options incl. outside

  // Validate target shares length
  if (static_cast<int>(target_shares.n_elem) != num_alts) {
    Rcpp::stop("Error: target_shares must have length %d (total alternatives, incl. outside if present).", num_alts);
  }

  // ---------------------------------------------------------------------------
  // Harmonize delta input:
  // - include_outside_option = TRUE: delta must have length J_inside (inside
  //   alternatives only). The outside option has no ASC (utility normalized to 0).
  // - include_outside_option = FALSE: allow J-1 (free) or J (full) length. Pad
  //   a leading zero if only free deltas are provided; keep baseline anchored.
  // ---------------------------------------------------------------------------
  arma::vec delta_current;  // length J_inside
  if (include_outside_option) {
    if (static_cast<int>(delta.n_elem) == J_inside) {
      delta_current = delta;
    } else {
      Rcpp::stop("Error: delta must have length %d (inside alternatives only).", J_inside);
    }
  } else {
    if (static_cast<int>(delta.n_elem) == J_inside) {
      delta_current = delta;
    } else if (static_cast<int>(delta.n_elem) == J_inside - 1) {
      delta_current = arma::zeros(J_inside);
      delta_current.subvec(1, J_inside - 1) = delta;  // pad baseline with 0
    } else {
      Rcpp::stop("Error: delta must have length %d (full) or %d (free, with baseline omitted).",
                 J_inside, J_inside - 1);
    }
    // Anchor baseline at zero for identification
    delta_current -= delta_current(0);
  }

  // The shares at a delta. X beta (+ W mu_final) is formed once, after the
  // check of the weights, where the shares function formed it at every call
  // (so a beta or mu of the wrong length stops at the same point); each pass
  // adds the current ASCs row by row.
  mxl_pred_setup(pd, mu_final, L, delta_current, /*use_asc=*/true,
                 include_outside_option, eta_draws, gen_seed, gen_scramble,
                 gen_S);
  const double weight_sum = shares_denominator(weights);
  const int n_threads = mxl_team_threads(pd.lay.N);
  std::vector<MxlPredScratch> scratch = mxl_pred_scratch(
      pd, n_threads, static_cast<std::size_t>(pd.lay.max_m) + 1);
  MxlPredAcc acc = mxl_pred_accumulators(scratch, num_alts, "shares");
  arma::vec base(pd.lay.n_rows, arma::fill::none);
  mxl_pred_base(pd, base, beta);
  if (max_iter > 0) mxl_pred_keep_gamma(pd, cache_bytes);
  double passes = 0;
  auto shares_at = [&](const arma::vec& d) {
    pd.delta = d;
    arma::vec shares = mxl_pred_shares(pd, weights, weight_sum, num_alts,
                                       scratch, acc);
    pd.gamma_filled = pd.gamma_cache != nullptr;  // filled by the first pass
    ++passes;
    return shares;
  };

  // Compute initial predicted shares
  arma::vec shares_pred = shares_at(delta_current);

  // Work with inside shares only for the contraction step
  arma::vec shares_pred_inside = include_outside_option
                                 ? shares_pred.subvec(1, num_alts - 1)
                                 : shares_pred;
  arma::vec target_shares_inside = include_outside_option
                                   ? target_shares.subvec(1, num_alts - 1)
                                   : target_shares;

  // Guard against zeros before taking logs
  auto validate_positive = [](const arma::vec& v, const char* name) {
    if (v.min() <= 0.0) {
      Rcpp::stop("%s contains non-positive entries; cannot take logarithm.", name);
    }
  };
  validate_positive(shares_pred_inside, "Predicted shares");
  validate_positive(target_shares_inside, "Target shares");

  arma::vec log_shares_old = arma::log(shares_pred_inside);
  arma::vec log_shares_target = arma::log(target_shares_inside);

  // Iteration
  int iter = 0;
  double residual = 10.0;

  while (iter < max_iter) {
    Rcpp::checkUserInterrupt(); // H4: allow user to interrupt long-running contraction
    arma::vec delta_new = delta_current + (log_shares_target - log_shares_old);

    // Re-anchor baseline each iteration when there is no outside option
    if (!include_outside_option) {
      delta_new -= delta_new(0);
    }

    residual = arma::max(arma::abs(delta_new - delta_current));

    if (residual < tol) {
      break;
    }

    delta_current = delta_new;
    shares_pred = shares_at(delta_current);

    shares_pred_inside = include_outside_option
                         ? shares_pred.subvec(1, num_alts - 1)
                         : shares_pred;
    validate_positive(shares_pred_inside, "Predicted shares");
    log_shares_old = arma::log(shares_pred_inside);
    ++iter;
  }

  if (iter >= max_iter) {
    Rcpp::Rcout << "Warning: Maximum iterations reached without convergence. Residual: "
                << residual << std::endl;
  }

  if (keep) {
    keep->kept = pd.gamma_cache != nullptr;
    keep->passes = passes;
    keep->reads = 0;
    for (const MxlPredScratch& sc : scratch) {
      keep->reads += static_cast<double>(sc.gamma_reads);
    }
  }
  return delta_current;
}

//' BLP contraction mapping for mixed logit
//'
//' Finds the ASC (delta) parameters such that predicted market shares
//' match target shares, using the contraction mapping of Berry, Levinsohn,
//' and Pakes (1995).
//'
//' @param delta J-1 or J vector with initial guess for deltas (ASCs)
//' @param target_shares J vector with target market shares
//' @param X design matrix for fixed coefficients; sum(M_i) x K_x
//' @param W design matrix for random coefficients; sum(M_i) x K_w or J x K_w
//' @param beta K_x vector with fixed coefficients
//' @param mu K_w vector with mean parameters (raw, will be transformed if log-normal)
//' @param L_params Cholesky parameters vector
//' @param alt_idx sum(M) x 1 vector with indices of alternatives; 1-based indexing
//' @param M N x 1 vector with number of alternatives for each individual
//' @param weights N x 1 vector with weights for each observation
//' @param eta_draws Array with draws; K_w x S x N
//' @param rc_dist K_w vector indicating distribution (0=normal, 1=log-normal)
//' @param rc_correlation whether random coefficients are correlated
//' @param rc_mean whether mu parameters represent means (TRUE) or are zero (FALSE)
//' @param include_outside_option whether outside option is included
//' @param tol convergence tolerance (default 1e-8)
//' @param max_iter maximum iterations (default 1000)
//' @param gen_seed Integer. \code{< 0} (default) uses the materialized
//'   \code{eta_draws} cube; \code{>= 0} forms the draws on the fly (see
//'   \code{gen_scramble}), with this master seed for the digit permutations of
//'   \code{gen_scramble = 1}.
//' @param gen_scramble Integer mode of the draws formed on the fly: \code{0} =
//'   identity permutations (plain Halton) through the generator's inverse normal
//'   CDF, \code{1} = seeded position-wise digit permutations, \code{2} = identity
//'   permutations through R's \code{qnorm()}, the store-mode draws of
//'   \code{\link{get_halton_normals}} formed without its cube (bit for bit where
//'   randtoolbox and choicer are compiled with the same floating-point
//'   contraction). Other values are an error.
//' @param gen_S Integer number of draws per individual, used only when \code{gen_seed >= 0}.
//' @details With \code{gen_seed >= 0} the draws are formed again at every
//'   evaluation of the shares. \code{\link{blp}()} on a fitted model runs the
//'   same contraction but keeps them across iterations within its
//'   \code{keep_draws_bytes} (see \code{\link{blp.choicer_mxl}}).
//' @returns vector with converged delta (ASC) values
//' @examples
//' \donttest{
//' library(data.table)
//' set.seed(42)
//' N <- 50; J <- 3
//' dt <- data.table(id = rep(1:N, each = J), alt = rep(1:J, N))
//' dt[, `:=`(x1 = rnorm(.N), w1 = rnorm(.N))]
//' dt[, choice := 0L]
//' dt[, choice := sample(c(1L, rep(0L, J - 1))), by = id]
//' d <- prepare_mxl_data(dt, "id", "alt", "choice", "x1", "w1")
//' eta <- get_halton_normals(50, d$N, ncol(d$W))
//' fit <- run_mxlogit(input_data = d, eta_draws = eta)
//' pm <- fit$param_map
//' delta <- mxl_blp_contraction(rep(0, J), rep(1/J, J), d$X, d$W,
//'   coef(fit)[pm$beta], rep(0, ncol(d$W)), coef(fit)[pm$sigma],
//'   d$alt_idx, d$M, d$weights, eta, rc_dist = rep(0L, ncol(d$W)),
//'   rc_correlation = FALSE, rc_mean = FALSE)
//' delta
//' }
//' @export
// [[Rcpp::export]]
arma::vec mxl_blp_contraction(
    const arma::vec& delta,
    const arma::vec& target_shares,
    const arma::mat& X,
    const arma::mat& W,
    const arma::vec& beta,
    const arma::vec& mu,
    const arma::vec& L_params,
    const Rcpp::IntegerVector& alt_idx,
    const Rcpp::IntegerVector& M,
    const arma::vec& weights,
    const arma::cube& eta_draws,
    const arma::uvec& rc_dist,
    const bool rc_correlation = true,
    const bool rc_mean = false,
    const bool include_outside_option = false,
    const double tol = 1e-8,
    const int max_iter = 1000,
    const int gen_seed = -1,
    const int gen_scramble = 1,
    const int gen_S = 0
) {
  return mxl_blp_run(delta, target_shares, X, W, beta, mu, L_params, alt_idx,
                     M, weights, eta_draws, rc_dist, rc_correlation, rc_mean,
                     include_outside_option, tol, max_iter, gen_seed,
                     gen_scramble, gen_S, 0.0, nullptr);
}

//' BLP contraction for blp(): mxl_blp_contraction(), keeping the draws
//'
//' The contraction of mxl_blp_contraction(), which, when the draws are formed
//' on the fly (gen_seed >= 0) and all situations' Gamma = L eta together take
//' at most cache_bytes bytes (8 K_w S N), keeps them from the first
//' evaluation of the shares on, instead of forming them at every evaluation.
//' The result is the same.
//'
//' @inheritParams mxl_blp_contraction
//' @param cache_bytes Largest cache, in bytes; 0 keeps nothing.
//' @return List: delta (as mxl_blp_contraction() returns it); kept, whether
//'   the cache was allocated; passes, the evaluations of the shares; and
//'   reads, the situations those evaluations read from the cache.
//' @noRd
// [[Rcpp::export]]
Rcpp::List mxl_blp_contraction_cached(
    const arma::vec& delta,
    const arma::vec& target_shares,
    const arma::mat& X,
    const arma::mat& W,
    const arma::vec& beta,
    const arma::vec& mu,
    const arma::vec& L_params,
    const Rcpp::IntegerVector& alt_idx,
    const Rcpp::IntegerVector& M,
    const arma::vec& weights,
    const arma::cube& eta_draws,
    const arma::uvec& rc_dist,
    const bool rc_correlation,
    const bool rc_mean,
    const bool include_outside_option,
    const double tol,
    const int max_iter,
    const int gen_seed,
    const int gen_scramble,
    const int gen_S,
    const double cache_bytes
) {
  MxlBlpKeep keep;
  const arma::vec d = mxl_blp_run(delta, target_shares, X, W, beta, mu,
                                  L_params, alt_idx, M, weights, eta_draws,
                                  rc_dist, rc_correlation, rc_mean,
                                  include_outside_option, tol, max_iter,
                                  gen_seed, gen_scramble, gen_S, cache_bytes,
                                  &keep);
  return Rcpp::List::create(Rcpp::Named("delta") = d,
                            Rcpp::Named("kept") = keep.kept,
                            Rcpp::Named("passes") = keep.passes,
                            Rcpp::Named("reads") = keep.reads);
}

//' Compute aggregate elasticities for mixed logit model
//'
//' Computes the aggregate elasticity matrix (weighted average of individual
//' elasticities) for the Mixed Logit model. The elasticity E(i,j) represents
//' the percentage change in the probability of choosing alternative i when
//' the attribute of alternative j changes by 1%.
//'
//' @param theta parameter vector (beta, \[mu\], L, delta)
//' @param X design matrix for fixed coefficients; sum(M_i) x K_x
//' @param W design matrix for random coefficients; sum(M_i) x K_w or J x K_w
//' @param alt_idx sum(M) x 1 vector with indices of alternatives; 1-based indexing
//' @param choice_idx N x 1 vector (kept for API consistency, not used)
//' @param M N x 1 vector with number of alternatives for each individual
//' @param weights N x 1 vector with weights for each observation
//' @param eta_draws Array with draws; K_w x S x N
//' @param rc_dist K_w vector indicating distribution (0=normal, 1=log-normal)
//' @param elast_var_idx 1-based index of the variable for elasticity computation
//' @param is_random_coef TRUE if variable is in W (random coef), FALSE if in X (fixed coef)
//' @param rc_correlation whether random coefficients are correlated
//' @param rc_mean whether mu parameters are estimated
//' @param use_asc whether ASCs are included
//' @param include_outside_option whether outside option is included
//' @param gen_seed Integer. \code{< 0} (default) uses the materialized
//'   \code{eta_draws} cube; \code{>= 0} forms the draws on the fly (see
//'   \code{gen_scramble}), with this master seed for the digit permutations of
//'   \code{gen_scramble = 1}.
//' @param gen_scramble Integer mode of the draws formed on the fly: \code{0} =
//'   identity permutations (plain Halton) through the generator's inverse normal
//'   CDF, \code{1} = seeded position-wise digit permutations, \code{2} = identity
//'   permutations through R's \code{qnorm()}, the store-mode draws of
//'   \code{\link{get_halton_normals}} formed without its cube (bit for bit where
//'   randtoolbox and choicer are compiled with the same floating-point
//'   contraction).
//' @param gen_S Integer number of draws per individual, used only when \code{gen_seed >= 0}.
//' @returns J x J matrix of aggregate elasticities
//' @examples
//' \donttest{
//' library(data.table)
//' set.seed(42)
//' N <- 50; J <- 3
//' dt <- data.table(id = rep(1:N, each = J), alt = rep(1:J, N))
//' dt[, `:=`(x1 = rnorm(.N), w1 = rnorm(.N))]
//' dt[, choice := 0L]
//' dt[, choice := sample(c(1L, rep(0L, J - 1))), by = id]
//' d <- prepare_mxl_data(dt, "id", "alt", "choice", "x1", "w1")
//' eta <- get_halton_normals(50, d$N, ncol(d$W))
//' fit <- run_mxlogit(input_data = d, eta_draws = eta)
//' elas <- choicer:::mxl_elasticities_parallel(coef(fit), d$X, d$W, d$alt_idx,
//'   d$choice_idx, d$M, d$weights, eta, rc_dist = rep(0L, ncol(d$W)),
//'   elast_var_idx = 1L, is_random_coef = FALSE,
//'   rc_correlation = FALSE, rc_mean = FALSE)
//' elas
//' }
//' @keywords internal
// [[Rcpp::export]]
arma::mat mxl_elasticities_parallel(
    const arma::vec& theta,
    const arma::mat& X,
    const arma::mat& W,
    const Rcpp::IntegerVector& alt_idx,
    SEXP choice_idx,  // unused, kept for API consistency; never converted
    const Rcpp::IntegerVector& M,
    const arma::vec& weights,
    const arma::cube& eta_draws,
    const arma::uvec& rc_dist,
    const int elast_var_idx,
    const bool is_random_coef,
    const bool rc_correlation = true,
    const bool rc_mean = false,
    const bool use_asc = true,
    const bool include_outside_option = false,
    const int gen_seed = -1, const int gen_scramble = 1, const int gen_S = 0
) {
  (void)choice_idx;  // unused, kept for API consistency

  // Basic dimensions
  const int K_x = X.n_cols;
  const int K_w = W.n_cols;

  // Convert 1-based R index to 0-based C++ index
  const int var_idx = elast_var_idx - 1;

  // Validate variable index
  if (is_random_coef) {
    if (var_idx < 0 || var_idx >= K_w) {
      Rcpp::stop("elast_var_idx (%d) is out of bounds for W matrix (K_w=%d).",
                 elast_var_idx, K_w);
    }
  } else {
    if (var_idx < 0 || var_idx >= K_x) {
      Rcpp::stop("elast_var_idx (%d) is out of bounds for X matrix (K_x=%d).",
                 elast_var_idx, K_x);
    }
  }

  // Parse theta into parameter blocks (shared helper; validates theta), then
  // the inputs and their layout, on the primary thread
  const MxlParams par = parse_mxl_theta(theta, K_x, K_w, rc_dist,
                                        rc_correlation, rc_mean, use_asc,
                                        include_outside_option);
  const bool generate = gen_seed >= 0;
  if (generate) mxl_pred_check_generate(gen_S, K_w, gen_scramble);
  MxlPredData pd(X, W, rc_dist);
  mxl_pred_layout(pd, alt_idx, M, eta_draws, !generate, use_asc, par.delta,
                  &weights);
  mxl_pred_setup(pd, par.mu_final, par.L, par.delta, use_asc,
                 include_outside_option, eta_draws, gen_seed, gen_scramble,
                 gen_S);
  const ChoiceLayout& lay = pd.lay;
  const double beta_k = is_random_coef ? 0.0 : par.beta(var_idx);
  const arma::vec& mu_final = par.mu_final;

  // Determine total number of alternatives
  const int J_inside = compute_J_inside(use_asc, par.delta, lay);
  const int J_total = compute_J_total(J_inside, include_outside_option);
  mxl_pred_check_alternatives(J_total);

  // Every thread's buffers (aux: a situation's elasticity terms, mean
  // probabilities and x_k, sized for the largest choice set) and
  // accumulators, then the base utilities
  const int n_threads = mxl_team_threads(lay.N);
  const std::size_t n_max = static_cast<std::size_t>(lay.max_m) + 1;
  std::vector<MxlPredScratch> scratch =
      mxl_pred_scratch(pd, n_threads, n_max * n_max + 2 * n_max);
  const std::size_t J_n = static_cast<std::size_t>(J_total);
  MxlPredAcc elas = mxl_pred_accumulators(scratch, J_n * J_n, "elasticity");
  MxlPredAcc total_weight = mxl_pred_accumulators(scratch, 1, "elasticity");
  arma::vec base(lay.n_rows, arma::fill::none);
  mxl_pred_base(pd, base, par.beta);

  mxl_pred_zero(elas, J_n * J_n);
  mxl_pred_zero(total_weight, 1);
  const int S = pd.S;
  const int o = include_outside_option ? 1 : 0;  // outside option: slot 0

  // The accumulations below are the statements of the code this replaced, on
  // Armadillo matrices over the threads' buffers (see mxl_pred_shares()).
  mxl_pred_run(pd, scratch, [&](const mxl_off t, MxlPredScratch& sc) {
    arma::mat local_elas_matrix(elas[sc.tid].get(), J_total, J_total, false,
                                true);
    const int m_i = mxl_pred_load(pd, t, sc);
    const int num_choices = m_i + o;
    const double w_i = weights[t];
    const double* bu = sc.bu.get();
    double* v = sc.v.get();
    const arma::mat Gamma_final(sc.gamma.get(), pd.K_w, S, false, true);

    // Map local indices to global alternative indices
    int* global_j_map = sc.map.get();
    fill_global_alt_map(global_j_map, lay.alt + lay.row_off[t], m_i,
                        include_outside_option);

    // Get attribute values for the elasticity variable
    arma::vec x_k_i(sc.aux.get() + n_max * n_max + n_max, num_choices, false,
                    true);
    mxl_pred_x_k(pd, t, sc, m_i, var_idx, is_random_coef, x_k_i);

    // Compute P_bar (average probabilities) and accumulate elasticity terms
    arma::vec P_bar_i(sc.aux.get() + n_max * n_max, num_choices, false, true);
    arma::mat elas_accum(sc.aux.get(), num_choices, num_choices, false, true);
    P_bar_i.zeros();
    elas_accum.zeros();
    arma::vec P_s(sc.p.get(), num_choices, false, true);

    for (int s = 0; s < S; ++s) {
      // Get effective coefficient for this draw
      double beta_k_eff;
      if (is_random_coef) {
        beta_k_eff = mu_final(var_idx) + Gamma_final(var_idx, s);
      } else {
        beta_k_eff = beta_k;
      }

      const double* wg = sc.WGamma.get() + static_cast<std::size_t>(s) * m_i;
      if (o) v[0] = 0.0;
      for (int a = 0; a < m_i; ++a) v[o + a] = bu[a] + wg[a];
      stable_softmax_n(v, P_s.memptr(), num_choices);

      P_bar_i += P_s;

      // Accumulate elasticity terms for this draw
      for (int j_local = 0; j_local < num_choices; ++j_local) {
        const double P_j = P_s(j_local);

        for (int m_local = 0; m_local < num_choices; ++m_local) {
          const double P_m = P_s(m_local);
          const double x_km = x_k_i(m_local);

          double elas_term;
          if (j_local == m_local) {
            // Own-elasticity: beta_k * x_k * P_j * (1 - P_j)
            elas_term = beta_k_eff * x_km * P_j * (1.0 - P_j);
          } else {
            // Cross-elasticity: -beta_k * x_km * P_j * P_m
            elas_term = -beta_k_eff * x_km * P_j * P_m;
          }

          elas_accum(j_local, m_local) += elas_term;
        }
      }
    }  // end s loop

    P_bar_i /= static_cast<double>(S);
    elas_accum /= static_cast<double>(S);

    // Compute final elasticities: E = elas_accum / P_bar
    // and map to global indices
    for (int j_local = 0; j_local < num_choices; ++j_local) {
      const int global_j = global_j_map[j_local];
      const double P_bar_j = P_bar_i(j_local);

      if (P_bar_j > 1e-12) {  // Avoid division by zero
        for (int m_local = 0; m_local < num_choices; ++m_local) {
          const int global_m = global_j_map[m_local];
          double elasticity = elas_accum(j_local, m_local) / P_bar_j;
          local_elas_matrix(global_j, global_m) += w_i * elasticity;
        }
      }
    }

    total_weight[sc.tid][0] += w_i;
  });

  // Global accumulators: the threads' added up in thread order
  arma::mat global_elas_matrix = arma::zeros(J_total, J_total);
  double global_total_weight = 0.0;
  for (int i = 0; i < n_threads; ++i) {
    global_elas_matrix += arma::mat(elas[i].get(), J_total, J_total, false, true);
    global_total_weight += total_weight[i][0];
  }

  // Compute weighted average
  if (global_total_weight > 1e-10) {
    global_elas_matrix /= global_total_weight;
  }

  return global_elas_matrix;
}
