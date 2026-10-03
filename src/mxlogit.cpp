// [[Rcpp::depends(RcppArmadillo)]]
#include "choicer.h"
#include "choicer_internal.h"
#include "halton.h"
#include <cstring>
#include <memory>

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

// Build and validate the layout, with the messages of mxl_unit_offsets() and
// validate_choice_data(), in that order, followed by the choices' check; the
// delta, draw and W checks follow in MxlUnitData, in validate_mxl_inputs()'s
// order. Situations are sorted by decision maker; Ti = NULL makes every
// situation its own unit (the cross-section).
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
//'   identity permutations (plain Halton, compat), \code{1} = seeded position-wise digit permutations.
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
//'   identity permutations (plain Halton, compat), \code{1} = seeded position-wise digit permutations.
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
arma::mat mxl_hessian_parallel(
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

  // Block layout: continuous block c = [beta | mu | L], size Kc = idx_delta_start
  //               delta block d = [delta ASCs],         size Jd = n_params - idx_delta_start
  // H_V (second derivative of utility w.r.t. theta) is nonzero only in the
  // (mu, L) x (mu, L) sub-block, which lives entirely within the continuous block.
  const int Kc = idx_delta_start;           // size of continuous block
  const int Jd = n_params - idx_delta_start; // size of delta block (0 when !use_asc)

  // Per unit u (Louis 1982), with omega_s the posterior draw weights and
  // g_s = sum_t g_ts, H_s = sum_t H_ts the per-draw score and Hessian of
  // lambda_s = sum_t log P_ts(j_t):
  //   H_u = sum_s omega_s (H_s + (g_s - g_bar)(g_s - g_bar)'),  g_bar = sum_s omega_s g_s,
  //   H_ts = -sum_Pzz_ts + sum_Pz_ts sum_Pz_ts' + sum_diff_H_V_ts.
  // Pass 1 (the shared draw loop) yields omega; pass 2 accumulates the O3
  // buffers with omega_s where the cross-section used P_choice_s / P_i_hat.

  // Global accumulator
  arma::mat global_hess = arma::zeros(n_params, n_params);

#ifdef _OPENMP
#pragma omp parallel
#endif
  {
    arma::mat local_hess = arma::zeros(n_params, n_params);

    // Thread-private scratch — sized once, reset inside loops as needed.
    MxlUnitScratch sc; // unit slices, draws, probabilities P_ts, weights omega

    // Per-(t, s) block accumulators (reset at the start of each situation-draw).
    // cc: Kc x Kc dense outer-product sum
    arma::mat sum_Pzz_cc(Kc, Kc);
    // cd: Kc x Jd; each col j accumulates P_a * zc_a for the alt whose delta is j
    // Note: the "Jd > 0 ? Jd : 1" sentinel below (and for all Jd-sized buffers) avoids
    // zero-size allocation when there are no ASC parameters; these buffers are never read
    // when Jd == 0 because every access is guarded by "if (Jd > 0)" / "if (j_delta >= 0)".
    arma::mat sum_Pzz_cd(Kc, Jd > 0 ? Jd : 1);
    // dd: diagonal only — stored as length-Jd vector
    arma::vec sum_Pzz_dd(Jd > 0 ? Jd : 1);
    // P*z sums
    arma::vec sum_Pz_c(Kc);
    arma::vec sum_Pz_d(Jd > 0 ? Jd : 1);
    // gradient components
    arma::vec g_c(Kc);
    arma::vec g_d(Jd > 0 ? Jd : 1);
    // H_V in the continuous block (mu,L sub-block; beta and delta rows/cols are zero)
    arma::mat sum_diff_H_V_cc(Kc, Kc);

    // Per-alt scratch for the continuous-block z vector and its outer product
    arma::vec zc_a(Kc);
    arma::mat zz(Kc, Kc);

    // O3: Per-unit block buffers — accumulate the pieces linear in omega_s
    // across the unit's situations and draws. These replace the per-draw H_is
    // allocation and a running sum_s omega_s (g_s g_sT + H_s) accumulation.
    arma::mat buf_Pzz_cc(Kc, Kc);               // Identity C: Σ_ts ω_s sum_Pzz_cc_ts
    arma::mat buf_Pzz_cd(Kc, Jd > 0 ? Jd : 1); // Identity C: Σ_ts ω_s sum_Pzz_cd_ts
    arma::vec buf_Pzz_dd(Jd > 0 ? Jd : 1);      // Identity C: Σ_ts ω_s sum_Pzz_dd_ts
    arma::mat buf_diff_HV_cc(Kc, Kc);           // Identity C: Σ_ts ω_s sum_diff_H_V_cc_ts
    // O3: Column stashes for BLAS-3 batching.
    // G_stash(:,s) = sqrt(ω_s) * (g_s - g_bar), g_s = Σ_t g_ts
    //   → G Gᵀ = Σ_s ω_s (g_s - g_bar)(g_s - g_bar)ᵀ              (Identity A)
    // F_stash(:,s) = sqrt(ω_s) * [sum_Pz_c; sum_Pz_d]_ts → F Fᵀ = Σ_s ω_s sum_Pz sum_PzT
    //                                                     for situation t (Identity B)
    // These are allocated once per thread at max size (n_params x Sdraw).
    // F column s is fully written in draw s of every situation; G is zeroed
    // per unit and accumulates over the unit's situations.
    arma::mat G_stash(n_params, Sdraw);
    arma::mat F_stash(n_params, Sdraw);
    arma::mat opg_pz(n_params, n_params); // G Gᵀ + Σ_t F_t F_tᵀ
    arma::vec g_bar(n_params);            // reference score, then mean deviation
    // Per-unit assembly buffers, reused so that no unit or situation
    // allocates: the products F Fᵀ and G Gᵀ, and sqrt(ω). The unit's
    // Hessian is then assembled in place in opg_pz.
    arma::mat prod_buf(n_params, n_params);
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

      // O3: Initialize per-unit block buffers.
      buf_Pzz_cc.zeros();
      buf_diff_HV_cc.zeros();
      if (Jd > 0) {
        buf_Pzz_cd.zeros();
        buf_Pzz_dd.zeros();
      }
      G_stash.zeros();
      opg_pz.zeros();

      // --- Pass 2: situations t (outer) x draws s (inner). A unit in one
      // draw batch reuses pass 1's WGamma; for a longer one each situation's
      // rows of W_u Gamma_u (m_t x S) are formed in turn, so the thread's
      // memory grows with the largest choice set rather than the unit.
      const bool one_batch = sc.nb == 1;
      for (mxl_off t = sc.t0; t < sc.t1; ++t) {
        const int m_t = static_cast<int>(lay.row_off[t + 1] - lay.row_off[t]);
        const int num_choices = include_outside_option ? m_t + 1 : m_t;
        if (!one_batch) mxl_situation_wgamma(sc, lay.row_off[t] - sc.r0, m_t);

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
          if (Jd > 0) {
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

            // --- Delta index for this alt (may be -1 meaning no delta entry) ---
            int j_delta = -1; // index into delta block; -1 = no entry
            if (use_asc && Jd > 0) {
              const int id = lay.alt0(row);
              if (include_outside_option) {
                j_delta = id;           // always valid (id >= 0)
              } else if (id > 0) {
                j_delta = id - 1;       // first inside alt (id==0) has no ASC
              }
              // id==0 and !include_outside_option => j_delta stays -1
            }

            // --- Accumulate block sums ---
            const double P_a = P_s(a);
            const double diff = (a == chosen_alt ? 1.0 : 0.0) - P_a;

            // cc block: P_a * zc_a * zc_a^T  (rank-1 update)
            zz = zc_a * zc_a.t();
            sum_Pzz_cc += P_a * zz;
            sum_Pz_c   += P_a * zc_a;
            g_c        += diff * zc_a;

            // cd and dd blocks (delta scatter)
            if (j_delta >= 0) {
              sum_Pzz_cd.col(j_delta) += P_a * zc_a;
              sum_Pzz_dd(j_delta)     += P_a;
              sum_Pz_d(j_delta)       += P_a;
              g_d(j_delta)            += diff;
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
          // Instead of assembling H_ts (n_params x n_params) and accumulating
          //   hess_term1 += omega_s * ((g_s - g_bar)(g_s - g_bar)T + sum_t H_ts),
          // we split the linear-in-omega_s and the outer-product pieces:
          //
          //   Linear (Identity C): buf_Pzz_cc     += omega_s * sum_Pzz_cc
          //                        buf_diff_HV_cc += omega_s * sum_diff_H_V_cc
          //                        (and cd/dd variants)
          //   Outer-product (Identity A): G_stash(:,s) += g_ts, centered at
          //     g_bar and scaled by sqrt(omega_s) after the unit, so that
          //     G GT = Σ_s omega_s (g_s - g_bar)(g_s - g_bar)T.
          //   Outer-product (Identity B): F_stash(:,s) = sqrt(omega_s) * [sum_Pz_c; sum_Pz_d]
          //     so that F FT = Σ_s omega_s sum_Pz sum_PzT after situation t.

          const double sqrt_ws = std::sqrt(omega_s);

          // Identity C — scalar-times-matrix accumulation into per-unit buffers.
          buf_Pzz_cc    += omega_s * sum_Pzz_cc;
          buf_diff_HV_cc += omega_s * sum_diff_H_V_cc;
          if (Jd > 0) {
            buf_Pzz_cd += omega_s * sum_Pzz_cd;
            buf_Pzz_dd += omega_s * sum_Pzz_dd;
          }

          // Identity A — accumulate g_s = Σ_t g_ts in column s (g_ts = [g_c; g_d]).
          G_stash.col(s).head(Kc) += g_c;
          if (Jd > 0) {
            G_stash.col(s).tail(Jd) += g_d;
          }

          // Identity B — fill column s of F_stash with sqrt(omega_s) * [sum_Pz_c; sum_Pz_d].
          F_stash.col(s).head(Kc) = sqrt_ws * sum_Pz_c;
          if (Jd > 0) {
            F_stash.col(s).tail(Jd) = sqrt_ws * sum_Pz_d;
          }
        } // end S loop

        // Identity B for situation t (one BLAS-3 product per situation: the
        // outer products do not combine across situations before the product).
        prod_buf = F_stash * F_stash.t();
        opg_pz += prod_buf;
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
      prod_buf = G_stash * G_stash.t();
      opg_pz += prod_buf;

      // === 6. O3: Per-unit finalization — assemble Hessian once from buffers ===
      // 6a. Assemble hess_term1 block by block, in place in opg_pz (zeroed
      //     again for the next unit):
      //     hess_term1 = (G Gᵀ + Σ_t F_t F_tᵀ) + (-buf_Pzz) + buf_diff_HV_cc (cc block only)
      //     This is the batched equivalent of Σ_s omega_s (H_s +
      //     (g_s - g_bar)(g_s - g_bar)T), the unit's whole Hessian H_u.

      // cc block: contributions from OPG, sum-Pz outer product, -Pzz, and H_V.
      opg_pz.submat(0, 0, Kc - 1, Kc - 1) -= buf_Pzz_cc;
      opg_pz.submat(0, 0, Kc - 1, Kc - 1) += buf_diff_HV_cc;

      if (Jd > 0) {
        // cd block: -buf_Pzz_cd + opg_pz cd block (NOT symmetric — full rectangular).
        opg_pz.submat(0, Kc, Kc - 1, n_params - 1) -= buf_Pzz_cd;
        for (int j = 0; j < Jd; ++j) {        // dc = (cd)ᵀ
          for (int i = 0; i < Kc; ++i) opg_pz(Kc + j, i) = opg_pz(i, Kc + j);
        }

        // dd block: -diag(buf_Pzz_dd) + opg_pz dd block (sum_Pz_d outer products).
        opg_pz.submat(Kc, Kc, n_params - 1, n_params - 1).diag() -= buf_Pzz_dd;
      }

      // 6b. Louis identity: H_u = hess_term1, whose centered Identity A
      //     already subtracts g_bar g_barᵀ. The draw weights are normalized,
      //     so there is no division by the simulated P_u.
      local_hess += w_u * opg_pz;
    } // end unit loop

#ifdef _OPENMP
#pragma omp critical
#endif
    {
      global_hess += local_hess;
    }
  } // end parallel region

  return -global_hess;
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
//'   identity permutations (plain Halton, compat), \code{1} = seeded position-wise digit permutations.
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
// cube slice t (store mode) or Halton block t + 1 (generate mode), in a panel
// as in the cross-section (predictions are unconditional). A kernel validates
// its inputs and lays out the stacked design on the primary thread (the
// alternative codes are read in place, row offsets are 64-bit), forms the base
// utilities once and allocates every thread's buffers there, then loops over
// the situations in parallel without allocating.
//
// Every result is bitwise that of the per-situation Armadillo code this
// replaced at the same thread count (choicer's own per-situation work does not
// depend on it; a multithreaded BLAS may split the base product, as before),
// with BLAS libraries whose results do not depend on operand addresses, such
// as the reference BLAS and OpenBLAS:
//   * the base utilities are the same full products, X beta, then
//     += W mu_final for a row-aligned W (one dgemv with beta = 1 into
//     X beta), with an alternative-level W's W mu_final and then the ASCs
//     added per row, so no BLAS call changes its shape;
//   * Gamma = L eta and W_t Gamma are the same Armadillo products on
//     matrices of the same shapes (Armadillo picks its BLAS call, or its own
//     code for tiny matrices, from the shapes alone), and the draws reach them
//     through a copy into a thread buffer, as they did;
//   * the draw loops use stable_softmax_n() and max_shifted_lse_n()
//     (Armadillo's operations in Armadillo's order) and element-wise sums.
// ============================================================================

// The team of a prediction kernel's parallel regions, and so the number of
// scratch sets it allocates: OpenMP's next team, within the thread limit, and
// at most one thread per situation, but two for a single situation, so that
// its region stays active and, as before, an OpenMP-built BLAS runs the
// situation's products single-threaded.
inline int mxl_pred_threads(const mxl_off n_situations) {
#ifdef _OPENMP
  mxl_off n = std::min(omp_get_max_threads(), omp_get_thread_limit());
#else
  mxl_off n = 1;
#endif
  n = std::min(n, std::max<mxl_off>(n_situations, 2));
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

// An n x 1 R matrix, the form in which RcppArmadillo returns an arma::vec
// (choicer.h includes <RcppArmadillo.h>, which leaves
// RCPP_ARMADILLO_RETURN_COLVEC_AS_VECTOR undefined), not initialized: the
// kernel writes every element.
inline Rcpp::NumericVector mxl_col_result(const mxl_off n) {
  Rcpp::NumericVector x(Rcpp::no_init(static_cast<R_xlen_t>(n)));
  x.attr("dim") = Rcpp::Dimension(static_cast<std::size_t>(n), 1);
  return x;
}

// Where situation t's K_w x S draws come from: Halton block t + 1, formed on
// the fly (generate mode), or slice t - c0 of the store-mode draws in memory,
// which hold situations [c0, c1). The primary thread loads them between
// parallel regions (load()); they are read-only inside.
struct MxlPredDraws {
  int K_w = 0, S = 0;
  bool generate = false;
  HaltonGen gen;                 // generate mode
  const double* cube = nullptr;  // store mode: the draws of situations [c0, c1)
  mxl_off c0 = 0;

  // Make the draws of situations [first, c1) available and return c1 > first:
  // all of them, from the generator or the whole cube.
  mxl_off load(const mxl_off first, const mxl_off N) {
    c0 = 0;
    (void)first;
    return N;
  }

  // Write situation t's draws (column-major K_w x S) to eta.
  void fill(double* eta, const mxl_off t) const {
    if (generate) {
      gen.fill_block(eta, static_cast<uint64_t>(t) * static_cast<uint64_t>(S) + 1);
      return;
    }
    const std::size_t n = static_cast<std::size_t>(K_w) * static_cast<std::size_t>(S);
    if (n > 0) {
      std::memcpy(eta, cube + static_cast<std::size_t>(t - c0) * n,
                  n * sizeof(double));
    }
  }
};

// What the per-situation routine reads, shared by all threads: the stacked
// design and its layout, the parameters, the base utilities and the draws.
// A kernel fills it on the primary thread in the order of its own checks;
// it holds no SEXP.
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
  int K_w = 0, S = 0;
  bool use_asc = false, include_outside_option = false, alt_level_W = false;

  MxlPredData(const arma::mat& X_, const arma::mat& W_, const arma::uvec& rc_dist_)
      : X(X_), W(W_), rc_dist(rc_dist_) {}
  MxlPredData(const MxlPredData&) = delete;
  MxlPredData& operator=(const MxlPredData&) = delete;
};

// Generate-mode checks of a prediction kernel, before its layout.
inline void mxl_pred_check_generate(const int gen_S, const int K_w) {
  if (gen_S <= 0) Rcpp::stop("gen_S must be positive when gen_seed >= 0");
  if (K_w > HALTON_N_PRIMES) {
    Rcpp::stop("K_w exceeds the primes table size (128); reduce K_w or "
               "extend the primes table.");
  }
}

// The layout of a prediction kernel, with the messages of
// validate_mxl_inputs() and in its order: the situations, the alternative
// codes and the ASCs' coverage of them (choice_layout_build()), then in store
// mode the cube's dimensions, then W's rows. Generate mode checks W's rows
// too (it used to stop at an Armadillo bounds error).
inline void mxl_pred_layout(MxlPredData& pd, const Rcpp::IntegerVector& alt_idx,
                            const Rcpp::IntegerVector& M,
                            const arma::cube& eta_draws, const bool store,
                            const bool use_asc, const arma::vec& delta,
                            const arma::vec* weights) {
  pd.lay = choice_layout_build(pd.X, alt_idx, M, use_asc, delta, weights);
  const ChoiceLayout& lay = pd.lay;
  if (store) {
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
    pd.draws.gen = HaltonGen(static_cast<uint64_t>(gen_seed), pd.S, pd.K_w,
                             gen_scramble);
  } else {
    pd.draws.cube = eta_draws.memptr();
  }
}

// The base utilities of every row into `base` (n rows; it may view a
// kernel's output): X beta, then += W mu_final for a row-aligned W, the
// expressions compute_base_util_mxl() formed before its ASCs. An
// alternative-level W's W mu_final goes to W_mu instead; mxl_pred_load() adds
// it and then the ASCs row by row, as that function did.
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
  int tid = 0;

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
                                           static_cast<double>(n_aux));
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
// own scratch: chunk by chunk of the draw source (a single chunk unless the
// store-mode draws come in chunks), each chunk one parallel region, so the
// primary thread loads a chunk's draws between regions.
template <typename Body>
inline void mxl_pred_run(MxlPredData& pd, std::vector<MxlPredScratch>& scratch,
                         Body body) {
  const int n_threads = static_cast<int>(scratch.size());
  for (mxl_off c0 = 0; c0 < pd.lay.N;) {
    const mxl_off c1 = pd.draws.load(c0, pd.lay.N);
    if (c1 <= c0) Rcpp::stop("Internal error: the draws made no progress.");
#ifdef _OPENMP
#pragma omp parallel num_threads(n_threads)
#endif
    {
      MxlPredScratch& sc = scratch[mxl_thread_num()];
#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
      for (mxl_off t = c0; t < c1; ++t) body(t, sc);
    }
    c0 = c1;
  }
}

// Load situation t into the thread's buffers: its base utilities (with an
// alternative-level W's W mu_final, then the ASCs, added row by row), its
// draws, Gamma = L eta with the log-normal transform, its rows of W and
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
  pd.draws.fill(sc.eta.get(), t);
  const arma::mat eta(sc.eta.get(), K_w, S, false, true);
  arma::mat Gamma(sc.gamma.get(), K_w, S, false, true);
  batch_gamma_draws_into(Gamma, pd.L, eta, pd.rc_dist);

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

// ============================================================================
// Mixed Logit: Share Prediction and BLP Contraction
// ============================================================================

// Internal function for computing simulated market shares
arma::vec mxl_predict_shares_internal(
    const arma::mat& X,
    const arma::mat& W,
    const arma::vec& beta,
    const arma::vec& mu_final,         // Transformed mu (exp(mu) for log-normal)
    const arma::mat& L,                // Cholesky factor
    const arma::uvec& alt_idx0,        // 0-based indexing
    const Rcpp::IntegerVector& M,
    const Rcpp::IntegerVector& S_prefix,
    const arma::vec& weights,
    const arma::vec& delta,            // Full J-element delta
    const arma::cube& eta_draws,
    const arma::uvec& rc_dist,
    const int num_alts,
    const bool use_asc,
    const bool include_outside_option,
    const int gen_seed = -1, const int gen_scramble = 1, const int gen_S = 0
) {
  const int N = M.size();
  const int K_w = W.n_cols;
  const int Sdraw = (gen_seed >= 0) ? gen_S : static_cast<int>(eta_draws.n_cols);
  // C++ runtime guards for generate mode
  if (gen_seed >= 0 && gen_S <= 0) Rcpp::stop("gen_S must be positive when gen_seed >= 0");
  if (gen_seed >= 0 && K_w > HALTON_N_PRIMES) Rcpp::stop("K_w exceeds the primes table size (128); reduce K_w or extend the primes table.");

  const bool use_generate_s = (gen_seed >= 0);
  HaltonGen halton_gen_s;
  if (use_generate_s) {
    halton_gen_s = HaltonGen(static_cast<uint64_t>(gen_seed), Sdraw, K_w, gen_scramble);
  }

  const double weight_sum = arma::accu(weights);

  if (weight_sum <= 0) {
    Rcpp::stop("Error: Sum of weights must be positive.");
  }

  // Pre-compute base utility for all individuals (single BLAS call)
  arma::vec base_util_s = compute_base_util_mxl(X, W, beta, mu_final,
                                          alt_idx0, use_asc, delta);

  // Initialize global accumulator for predicted shares
  arma::vec global_shares = arma::zeros(num_alts);

#ifdef _OPENMP
#pragma omp parallel
#endif
  {
    // Thread-local accumulator
    arma::vec local_shares = arma::zeros(num_alts);
    arma::mat eta_i_buf_s;    // generate mode scratch
    arma::mat eta_i_store_s;  // cube mode materialised slice

#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
    for (int i = 0; i < N; ++i) {
      const int m_i = M[i];
      const int num_choices = include_outside_option ? m_i + 1 : m_i;
      const int start_idx = S_prefix[i];
      const int end_idx = start_idx + m_i - 1;
      const double w_i = weights[i];
      const arma::uvec alt_idx0_i = alt_idx0.subvec(start_idx, end_idx);

      arma::mat W_i =
          make_W_i(W, X.n_rows, start_idx, end_idx, alt_idx0_i);

      // Pre-computed base utility for this individual
      const arma::vec base_util_i = base_util_s.subvec(start_idx, end_idx);

      // Accumulate probabilities over draws
      arma::vec P_bar_i = arma::zeros(num_choices);

      // --- Batch Cholesky: compute L * eta for all draws in one dgemm ---
      const arma::mat* eta_i_ptr_s;
      if (use_generate_s) {
        halton_gen_s.fill_eta_i(eta_i_buf_s, i + 1);
        eta_i_ptr_s = &eta_i_buf_s;
      } else {
        eta_i_store_s = eta_draws.slice(i);
        eta_i_ptr_s = &eta_i_store_s;
      }
      const arma::mat& eta_i_s_ref = *eta_i_ptr_s;
      arma::mat Gamma_final = batch_gamma_draws(L, eta_i_s_ref, rc_dist);

      // Batch W_i * Gamma_final into a single dgemm (m_i x Sdraw)
      const arma::mat WGamma = W_i * Gamma_final;

      arma::vec V_s(num_choices);
      arma::vec inside_utils(m_i);
      arma::vec P_s;

      for (int s = 0; s < Sdraw; ++s) {
        inside_utils = base_util_i + WGamma.col(s);
        fill_choice_utilities(V_s, inside_utils, num_choices,
                              include_outside_option);

        // Compute probabilities with numerical stability
        stable_softmax(V_s, P_s);

        P_bar_i += P_s;
      }  // end s loop

      P_bar_i /= static_cast<double>(Sdraw);

      // Accumulate shares by alternative
      if (include_outside_option) {
        local_shares(0) += w_i * P_bar_i(0);
      }
      for (int a = 0; a < m_i; ++a) {
        if (include_outside_option) {
          local_shares(alt_idx0_i(a) + 1) += w_i * P_bar_i(a + 1);
        } else {
          local_shares(alt_idx0_i(a)) += w_i * P_bar_i(a);
        }
      }
    }  // end i loop

#ifdef _OPENMP
#pragma omp critical
#endif
    {
      global_shares += local_shares;
    }
  }  // end parallel region

  return global_shares / weight_sum;
}

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
//' @param gen_seed Integer master seed for the on-the-fly Halton generator. \code{< 0}
//'   (default) uses the materialized \code{eta_draws} cube; \code{>= 0} generates draws
//'   on the fly from this seed.
//' @param gen_scramble Integer scramble mode for on-the-fly generation: \code{0} =
//'   identity permutations (plain Halton, compat), \code{1} = seeded position-wise digit permutations.
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
  if (generate) mxl_pred_check_generate(gen_S, W.n_cols);
  MxlPredData pd(X, W, rc_dist);
  mxl_pred_layout(pd, alt_idx, M, eta_draws, !generate, use_asc, par.delta,
                  nullptr);
  mxl_pred_setup(pd, par.mu_final, par.L, par.delta, use_asc,
                 include_outside_option, eta_draws, gen_seed, gen_scramble,
                 gen_S);
  const ChoiceLayout& lay = pd.lay;
  const std::size_t max_m = static_cast<std::size_t>(lay.max_m);
  std::vector<MxlPredScratch> scratch =
      mxl_pred_scratch(pd, mxl_pred_threads(lay.N), 2 * max_m);

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
//' @param gen_seed Integer master seed for the on-the-fly Halton generator. \code{< 0}
//'   (default) uses the materialized \code{eta_draws} cube; \code{>= 0} generates draws
//'   on the fly from this seed.
//' @param gen_scramble Integer scramble mode for on-the-fly generation: \code{0} =
//'   identity permutations (plain Halton, compat), \code{1} = seeded position-wise digit permutations.
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
  if (generate) mxl_pred_check_generate(gen_S, W.n_cols);
  MxlPredData pd(X, W, rc_dist);
  mxl_pred_layout(pd, alt_idx, M, eta_draws, !generate, use_asc, par.delta,
                  nullptr);
  mxl_pred_setup(pd, par.mu_final, par.L, par.delta, use_asc,
                 include_outside_option, eta_draws, gen_seed, gen_scramble,
                 gen_S);
  const ChoiceLayout& lay = pd.lay;

  std::vector<MxlPredScratch> scratch =
      mxl_pred_scratch(pd, mxl_pred_threads(lay.N), 0);
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
//' Exported wrapper around the internal `mxl_predict_shares_internal`. Parses
//' `theta` using the standard parameter ordering and returns the simulated
//' weighted-average market shares.
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
//' @param gen_seed Integer master seed for the on-the-fly Halton generator. \code{< 0}
//'   (default) uses the materialized \code{eta_draws} cube; \code{>= 0} generates draws
//'   on the fly from this seed.
//' @param gen_scramble Integer scramble mode for on-the-fly generation: \code{0} =
//'   identity permutations (plain Halton, compat), \code{1} = seeded position-wise digit permutations.
//' @param gen_S Integer number of draws per individual, used only when \code{gen_seed >= 0}.
//' @returns Vector of length J (or J+1 with outside option) of predicted shares.
//' @keywords internal
// [[Rcpp::export]]
arma::vec mxl_predict_shares(
    const arma::vec& theta,
    const arma::mat& X,
    const arma::mat& W,
    const arma::uvec& alt_idx,
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
  // Basic dimensions
  const int K_x = X.n_cols;
  const int K_w = W.n_cols;

  // Parse theta into parameter blocks (shared helper; validates theta)
  const MxlParams par = parse_mxl_theta(theta, K_x, K_w, rc_dist,
                                        rc_correlation, rc_mean, use_asc,
                                        include_outside_option);
  if (gen_seed < 0) {
    validate_mxl_inputs(X, W, alt_idx, M, eta_draws, use_asc, par.delta,
                        &weights);
  } else {
    if (gen_S <= 0) Rcpp::stop("gen_S must be positive when gen_seed >= 0");
    if (K_w > HALTON_N_PRIMES) Rcpp::stop("K_w exceeds the primes table size (128); reduce K_w or extend the primes table.");
    validate_choice_data(X, alt_idx, M, use_asc, par.delta, &weights);
    check_rc_dist_length(rc_dist, K_w);
  }
  const arma::vec& beta = par.beta;
  const arma::vec& mu_final = par.mu_final;
  const arma::mat& L = par.L;
  const arma::vec& delta = par.delta;

  // 0-based indexing and prefix sums
  arma::uvec alt_idx0 = alt_idx - 1;
  Rcpp::IntegerVector S_prefix = compute_prefix_sum(M);

  // Number of alternatives for the output
  const int J_inside = compute_J_inside(use_asc, delta, alt_idx0);
  const int num_alts = compute_J_total(J_inside, include_outside_option);

  return mxl_predict_shares_internal(
    X, W, beta, mu_final, L, alt_idx0, M, S_prefix, weights,
    delta, eta_draws, rc_dist, num_alts, use_asc, include_outside_option,
    gen_seed, gen_scramble, gen_S
  );
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
//' @param gen_seed Integer master seed for the on-the-fly Halton generator. \code{< 0}
//'   (default) uses the materialized \code{eta_draws} cube; \code{>= 0} generates draws
//'   on the fly from this seed.
//' @param gen_scramble Integer scramble mode for on-the-fly generation: \code{0} =
//'   identity permutations (plain Halton, compat), \code{1} = seeded position-wise digit permutations.
//' @param gen_S Integer number of draws per individual, used only when \code{gen_seed >= 0}.
//' @returns J x J (or (J+1) x (J+1)) matrix of diversion ratios with zero diagonal.
//' @keywords internal
// [[Rcpp::export]]
arma::mat mxl_diversion_ratios_parallel(
    const arma::vec& theta,
    const arma::mat& X,
    const arma::mat& W,
    const arma::uvec& alt_idx,
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
  const int N = M.size();
  const int K_x = X.n_cols;
  const int K_w = W.n_cols;
  const int Sdraw = (gen_seed >= 0) ? gen_S : static_cast<int>(eta_draws.n_cols);

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

  // Parse theta into parameter blocks (shared helper; validates theta)
  const MxlParams par = parse_mxl_theta(theta, K_x, K_w, rc_dist,
                                        rc_correlation, rc_mean, use_asc,
                                        include_outside_option);
  if (gen_seed < 0) {
    validate_mxl_inputs(X, W, alt_idx, M, eta_draws, use_asc, par.delta,
                        &weights);
  } else {
    if (gen_S <= 0) Rcpp::stop("gen_S must be positive when gen_seed >= 0");
    if (K_w > HALTON_N_PRIMES) Rcpp::stop("K_w exceeds the primes table size (128); reduce K_w or extend the primes table.");
    validate_choice_data(X, alt_idx, M, use_asc, par.delta, &weights);
    check_rc_dist_length(rc_dist, K_w);
  }
  const arma::vec& beta = par.beta;
  const double beta_k = is_random_coef ? 0.0 : beta(var_idx);
  const arma::vec& mu_final = par.mu_final;
  const arma::mat& L = par.L;
  const arma::vec& delta = par.delta;

  // 0-based alt indices
  arma::uvec alt_idx0 = alt_idx - 1;

  // Total alternatives for output matrix
  const int J_inside = compute_J_inside(use_asc, delta, alt_idx0);
  const int J_total = compute_J_total(J_inside, include_outside_option);

  // Prefix sums
  const Rcpp::IntegerVector S_prefix = compute_prefix_sum(M);

  // Pre-compute base utility for all individuals (single BLAS call)
  arma::vec base_util = compute_base_util_mxl(X, W, beta, mu_final,
                                          alt_idx0, use_asc, delta);

  // Construct on-the-fly generator outside parallel region.
  const bool use_generate_d = (gen_seed >= 0);
  HaltonGen halton_gen_d;
  if (use_generate_d) {
    halton_gen_d = HaltonGen(static_cast<uint64_t>(gen_seed), Sdraw, K_w, gen_scramble);
  }

  // Global accumulators
  arma::mat global_numerator = arma::zeros(J_total, J_total);
  arma::vec global_denominator = arma::zeros(J_total);

#ifdef _OPENMP
#pragma omp parallel
#endif
  {
    // Thread-local accumulators
    arma::mat local_numerator = arma::zeros(J_total, J_total);
    arma::vec local_denominator = arma::zeros(J_total);
    arma::mat eta_i_buf_d;
    arma::mat eta_i_store_d;

#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
    for (int i = 0; i < N; ++i) {
      const int m_i = M[i];
      const int num_choices = include_outside_option ? m_i + 1 : m_i;
      const int start_idx = S_prefix[i];
      const int end_idx = start_idx + m_i - 1;
      const double w_i = weights[i];
      const arma::uvec alt_idx0_i = alt_idx0.subvec(start_idx, end_idx);

      arma::mat W_i =
          make_W_i(W, X.n_rows, start_idx, end_idx, alt_idx0_i);

      // Pre-computed base utility for this individual
      const arma::vec base_util_i = base_util.subvec(start_idx, end_idx);

      // Map local indices to global alternative indices
      arma::uvec global_j_map =
          build_global_alt_map(alt_idx0_i, m_i, include_outside_option);

      // Per-individual accumulators (sum across draws, divided by S below)
      arma::mat ind_num = arma::zeros(num_choices, num_choices);
      arma::vec ind_den = arma::zeros(num_choices);

      // --- Batch Cholesky: compute L * eta for all draws in one dgemm ---
      const arma::mat* eta_i_ptr_d;
      if (use_generate_d) {
        halton_gen_d.fill_eta_i(eta_i_buf_d, i + 1);
        eta_i_ptr_d = &eta_i_buf_d;
      } else {
        eta_i_store_d = eta_draws.slice(i);
        eta_i_ptr_d = &eta_i_store_d;
      }
      const arma::mat& eta_i_d_ref = *eta_i_ptr_d;
      arma::mat Gamma_final = batch_gamma_draws(L, eta_i_d_ref, rc_dist);

      // Batch W_i * Gamma_final into a single dgemm (m_i x Sdraw)
      const arma::mat WGamma = W_i * Gamma_final;

      arma::vec V_s(num_choices);
      arma::vec inside_utils(m_i);
      arma::vec P_s;

      // Loop over draws — cross-products MUST be accumulated INSIDE this loop
      for (int s = 0; s < Sdraw; ++s) {
        const auto gamma_i_s_final = Gamma_final.col(s); // still needed for random coef value

        // Realized coefficient on the perturbed variable for this (i, s).
        // For a fixed coef this is the constant beta_k; for a random coef
        // it is mu_final(var_idx) + gamma_i_s_final(var_idx), already
        // transformed (exp(.) applied upstream when rc_dist == 1).
        const double beta_k_eff = is_random_coef
            ? (mu_final(var_idx) + gamma_i_s_final(var_idx))
            : beta_k;

        // CHANGE #2: use pre-computed WGamma column instead of W_i * gamma_i_s_final
        inside_utils = base_util_i + WGamma.col(s);
        fill_choice_utilities(V_s, inside_utils, num_choices,
                              include_outside_option);

        // Stable softmax
        stable_softmax(V_s, P_s);

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
      const double S_d = static_cast<double>(Sdraw);
      ind_num /= S_d;
      ind_den /= S_d;

      // Scatter individual contribution into thread-local globals
      for (int j_local = 0; j_local < num_choices; ++j_local) {
        const int global_j = global_j_map(j_local);
        local_denominator(global_j) += w_i * ind_den(j_local);
        for (int k_local = 0; k_local < num_choices; ++k_local) {
          if (k_local == j_local) continue;
          const int global_k = global_j_map(k_local);
          local_numerator(global_k, global_j) += w_i * ind_num(k_local, j_local);
        }
      }
    }  // end i loop

#ifdef _OPENMP
#pragma omp critical
#endif
    {
      global_numerator += local_numerator;
      global_denominator += local_denominator;
    }
  }  // end parallel region

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
//' @param gen_seed Integer master seed for the on-the-fly Halton generator. \code{< 0}
//'   (default) uses the materialized \code{eta_draws} cube; \code{>= 0} generates draws
//'   on the fly from this seed.
//' @param gen_scramble Integer scramble mode for on-the-fly generation: \code{0} =
//'   identity permutations, \code{1} = seeded position-wise digit permutations.
//' @param gen_S Integer number of draws per individual, used only when \code{gen_seed >= 0}.
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
    const arma::uvec& alt_idx,
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
  const int K_w = W.n_cols;
  const bool use_asc = true;

  // delta is harmonized to cover every referenced alternative below, so the
  // ASC-coverage check is skipped here (use_asc = false, empty delta).
  check_rc_dist_length(rc_dist, K_w);
  if (gen_seed < 0) {
    validate_mxl_inputs(X, W, alt_idx, M, eta_draws,
                        /*use_asc=*/false, arma::vec(), &weights);
  } else {
    if (gen_S <= 0) Rcpp::stop("gen_S must be positive when gen_seed >= 0");
    if (K_w > HALTON_N_PRIMES) Rcpp::stop("K_w exceeds the primes table size (128); reduce K_w or extend the primes table.");
    validate_choice_data(X, alt_idx, M, /*use_asc=*/false, arma::vec(),
                         &weights);
    check_rc_dist_length(rc_dist, K_w);
  }

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

  // Convert to 0-based indexing
  arma::uvec alt_idx0 = alt_idx - 1;

  // Deduce number of alternatives from data
  const int J_inside = static_cast<int>(arma::max(alt_idx0)) + 1;  // inside options
  const int num_alts = include_outside_option ? (J_inside + 1) : J_inside;  // total options incl. outside

  // Validate target shares length
  if (static_cast<int>(target_shares.n_elem) != num_alts) {
    Rcpp::stop("Error: target_shares must have length %d (total alternatives, incl. outside if present).", num_alts);
  }

  // Compute prefix sums
  Rcpp::IntegerVector S_prefix = compute_prefix_sum(M);

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

  // Compute initial predicted shares
  arma::vec shares_pred = mxl_predict_shares_internal(
    X, W, beta, mu_final, L, alt_idx0, M, S_prefix, weights,
    delta_current, eta_draws, rc_dist, num_alts, use_asc, include_outside_option,
    gen_seed, gen_scramble, gen_S
  );

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
    shares_pred = mxl_predict_shares_internal(
      X, W, beta, mu_final, L, alt_idx0, M, S_prefix, weights,
      delta_current, eta_draws, rc_dist, num_alts, use_asc, include_outside_option,
      gen_seed, gen_scramble, gen_S
    );

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

  return delta_current;
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
//' @param gen_seed Integer master seed for the on-the-fly Halton generator. \code{< 0}
//'   (default) uses the materialized \code{eta_draws} cube; \code{>= 0} generates draws
//'   on the fly from this seed.
//' @param gen_scramble Integer scramble mode for on-the-fly generation: \code{0} =
//'   identity permutations (plain Halton, compat), \code{1} = seeded position-wise digit permutations.
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
    const arma::uvec& alt_idx,
    const arma::uvec& choice_idx,
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
  const int N = M.size();
  const int K_x = X.n_cols;
  const int K_w = W.n_cols;
  const int Sdraw = (gen_seed >= 0) ? gen_S : static_cast<int>(eta_draws.n_cols);

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

  // Parse theta into parameter blocks (shared helper; validates theta)
  const MxlParams par = parse_mxl_theta(theta, K_x, K_w, rc_dist,
                                        rc_correlation, rc_mean, use_asc,
                                        include_outside_option);
  if (gen_seed < 0) {
    validate_mxl_inputs(X, W, alt_idx, M, eta_draws, use_asc, par.delta,
                        &weights);
  } else {
    if (gen_S <= 0) Rcpp::stop("gen_S must be positive when gen_seed >= 0");
    if (K_w > HALTON_N_PRIMES) Rcpp::stop("K_w exceeds the primes table size (128); reduce K_w or extend the primes table.");
    validate_choice_data(X, alt_idx, M, use_asc, par.delta, &weights);
    check_rc_dist_length(rc_dist, K_w);
  }
  const arma::vec& beta = par.beta;
  const double beta_k = is_random_coef ? 0.0 : beta(var_idx);
  const arma::vec& mu_final = par.mu_final;
  const arma::mat& L = par.L;
  const arma::vec& delta = par.delta;

  // Convert to 0-based indexing
  arma::uvec alt_idx0 = alt_idx - 1;

  // Determine total number of alternatives
  const int J_inside = compute_J_inside(use_asc, delta, alt_idx0);
  const int J_total = compute_J_total(J_inside, include_outside_option);

  // Compute prefix sums
  const Rcpp::IntegerVector S_prefix = compute_prefix_sum(M);

  // Pre-compute base utility for all individuals (single BLAS call)
  arma::vec base_util_e = compute_base_util_mxl(X, W, beta, mu_final,
                                          alt_idx0, use_asc, delta);

  // Construct on-the-fly generator outside parallel region.
  const bool use_generate_e = (gen_seed >= 0);
  HaltonGen halton_gen_e;
  if (use_generate_e) {
    halton_gen_e = HaltonGen(static_cast<uint64_t>(gen_seed), Sdraw, K_w, gen_scramble);
  }

  // Global accumulators
  arma::mat global_elas_matrix = arma::zeros(J_total, J_total);
  double global_total_weight = 0.0;

#ifdef _OPENMP
#pragma omp parallel
#endif
  {
    // Thread-local accumulators
    arma::mat local_elas_matrix = arma::zeros(J_total, J_total);
    double local_total_weight = 0.0;
    arma::mat eta_i_buf_e;
    arma::mat eta_i_store_e;

#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
    for (int i = 0; i < N; ++i) {
      const int m_i = M[i];
      const int num_choices = include_outside_option ? m_i + 1 : m_i;
      const int start_idx = S_prefix[i];
      const int end_idx = start_idx + m_i - 1;
      const double w_i = weights[i];
      const auto X_i = X.rows(start_idx, end_idx);
      const arma::uvec alt_idx0_i = alt_idx0.subvec(start_idx, end_idx);

      arma::mat W_i =
          make_W_i(W, X.n_rows, start_idx, end_idx, alt_idx0_i);

      // Pre-computed base utility for this individual
      const arma::vec base_util_i = base_util_e.subvec(start_idx, end_idx);

      // Map local indices to global alternative indices
      arma::uvec global_j_map =
          build_global_alt_map(alt_idx0_i, m_i, include_outside_option);

      // Get attribute values for the elasticity variable
      arma::vec x_k_i = arma::zeros(num_choices);
      if (is_random_coef) {
        if (include_outside_option) {
          x_k_i.subvec(1, num_choices - 1) = W_i.col(var_idx);
        } else {
          x_k_i = W_i.col(var_idx);
        }
      } else {
        if (include_outside_option) {
          x_k_i.subvec(1, num_choices - 1) = X_i.col(var_idx);
        } else {
          x_k_i = X_i.col(var_idx);
        }
      }

      // Compute P_bar (average probabilities) and accumulate elasticity terms
      arma::vec P_bar_i = arma::zeros(num_choices);
      arma::mat elas_accum = arma::zeros(num_choices, num_choices);

      // --- Batch Cholesky: compute L * eta for all draws in one dgemm ---
      const arma::mat* eta_i_ptr_e;
      if (use_generate_e) {
        halton_gen_e.fill_eta_i(eta_i_buf_e, i + 1);
        eta_i_ptr_e = &eta_i_buf_e;
      } else {
        eta_i_store_e = eta_draws.slice(i);
        eta_i_ptr_e = &eta_i_store_e;
      }
      const arma::mat& eta_i_e_ref = *eta_i_ptr_e;
      arma::mat Gamma_final = batch_gamma_draws(L, eta_i_e_ref, rc_dist);

      // Batch W_i * Gamma_final into a single dgemm (m_i x Sdraw)
      const arma::mat WGamma = W_i * Gamma_final;

      arma::vec V_s(num_choices);
      arma::vec inside_utils(m_i);
      arma::vec P_s;

      for (int s = 0; s < Sdraw; ++s) {
        const auto gamma_i_s_final = Gamma_final.col(s); // still needed for random coef value

        // Get effective coefficient for this draw
        double beta_k_eff;
        if (is_random_coef) {
          beta_k_eff = mu_final(var_idx) + gamma_i_s_final(var_idx);
        } else {
          beta_k_eff = beta_k;
        }

        // CHANGE #2: use pre-computed WGamma column instead of W_i * gamma_i_s_final
        inside_utils = base_util_i + WGamma.col(s);
        fill_choice_utilities(V_s, inside_utils, num_choices,
                              include_outside_option);

        // Compute probabilities
        stable_softmax(V_s, P_s);

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

      P_bar_i /= static_cast<double>(Sdraw);
      elas_accum /= static_cast<double>(Sdraw);

      // Compute final elasticities: E = elas_accum / P_bar
      // and map to global indices
      for (int j_local = 0; j_local < num_choices; ++j_local) {
        const int global_j = global_j_map(j_local);
        const double P_bar_j = P_bar_i(j_local);

        if (P_bar_j > 1e-12) {  // Avoid division by zero
          for (int m_local = 0; m_local < num_choices; ++m_local) {
            const int global_m = global_j_map(m_local);
            double elasticity = elas_accum(j_local, m_local) / P_bar_j;
            local_elas_matrix(global_j, global_m) += w_i * elasticity;
          }
        }
      }

      local_total_weight += w_i;
    }  // end i loop

#ifdef _OPENMP
#pragma omp critical
#endif
    {
      global_elas_matrix += local_elas_matrix;
      global_total_weight += local_total_weight;
    }
  }  // end parallel region

  // Compute weighted average
  if (global_total_weight > 1e-10) {
    global_elas_matrix /= global_total_weight;
  }

  return global_elas_matrix;
}
