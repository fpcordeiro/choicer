#ifndef CHOICER_HB_INTERNAL_HPP
#define CHOICER_HB_INTERNAL_HPP

#include "choicer.h"
#include "rng.h"
#include "bayes_samplers.h"
#include <R_ext/Utils.h>   // R_CheckUserInterrupt / R_ToplevelExec
#include <vector>

// ============================================================================
// Shared infrastructure for the hierarchical Bayes engines (hmnl_gibbs /
// hmnp_gibbs). Header-only; every helper below is either (a) callable from
// OpenMP worker threads because it uses hand-rolled fixed-order loops and no
// LAPACK/BLAS/R API, or (b) explicitly master-only (see each block). The
// engine safety contract is the one documented in src/mnprobit.cpp:15-36:
// no BLAS/LAPACK and no R RNG off the master thread, all cross-unit
// accumulations in a fixed order so draws are bitwise independent of the
// thread count.
//
// RNG tag partition (per iteration r, master seed `seed`, via
// make_stream(seed, r, tag) — src/rng.h). With N respondents and J inside
// alternatives:
//
//   tag i          (i = 0..N-1)  respondent i's latent sweep / RW-MH work
//   tag N + i      (i = 0..N-1)  respondent i's beta_i regression draw (HMNP)
//   tag 2N + j     (j = 0..J-1)  delta_j draw
//   tag 2N + J + 0               b draw            (master)
//   tag 2N + J + 1               W draw            (master)
//   tag 2N + J + 2               theta draw        (master)
//   tag 2N + J + 3               sigma_d2 (and, HMNP, sigma2) draws (master)
//
// The master block is based past the delta block (2N + J, not 3N), so the
// partition is collision-free for ANY J — including the J > N regime this
// design targets.
//
// delta-update asymmetry (documented here because both kernels share this
// header): the HMNL delta_j conditionals are coupled through the shared
// per-task softmax denominators, so the delta phase MUST be a SERIAL sweep
// over j — each delta_j update has to see the latest delta and the updated
// denominator caches, or the update does not leave the posterior invariant.
// The HMNP delta_j conditionals are Gaussian and conditionally INDEPENDENT
// across j given the augmented utilities U, so they may validly be drawn in
// one work-shared parallel pass. Same target distribution, opposite
// parallelization contract.
// ============================================================================

// ----------------------------------------------------------------------------
// Panel indexing: two-level (person, task) structure over the row-major data
// produced by .prepare_hb_panel() (R/hb_data.R). All offsets are CSR-style
// half-open ranges built with compute_prefix_sum (src/utils.cpp:98).
//
// build() allocates R vectors (through compute_prefix_sum), so it must be
// called on the master thread BEFORE the parallel region opens; afterwards
// the struct is read-only and safe to share across threads.
// ----------------------------------------------------------------------------
struct HbPanel {
  int N_persons = 0;      // respondents
  int n_tasks = 0;        // choice situations (sum of Ti)
  int J = 0;              // inside alternatives
  int total_rows = 0;     // inside rows (sum of M)
  bool include_outside = true;

  Rcpp::IntegerVector row_offsets;    // n_tasks + 1: task t owns rows
                                      //   [row_offsets[t], row_offsets[t+1])
  Rcpp::IntegerVector task_offsets;   // N_persons + 1: person i owns tasks
                                      //   [task_offsets[i], task_offsets[i+1])
  std::vector<int> alt_row_offsets;   // J + 1: alternative a owns the slots
                                      //   [alt_row_offsets[a], alt_row_offsets[a+1])
                                      //   of alt_rows
  std::vector<int> alt_rows;          // total_rows: 0-based row ids grouped by
                                      //   alternative, ascending row order
                                      //   within each alternative (fixed-order
                                      //   contract for delta sufficient stats)
  std::vector<int> alt_of_row;        // total_rows: 0-based alternative per row
  std::vector<int> task_of_row;       // total_rows: 0-based task per row
  std::vector<int> chosen_row;        // n_tasks: 0-based global row of the
                                      //   chosen inside alternative; -1 when
                                      //   the outside option was chosen
  std::vector<int> choice_pos;        // n_tasks: 1-based within-task position
                                      //   of the choice; 0 = outside chosen

  // Build from the prep outputs (M, Ti, alt_of_row, choice_pos are the
  // 1-based R-side vectors). Returns false on any inconsistency (sizes that
  // do not add up, out-of-range codes, a 0 choice without an outside option)
  // so kernels can abort cleanly instead of indexing out of bounds.
  bool build(const Rcpp::IntegerVector& M, const Rcpp::IntegerVector& Ti,
             const Rcpp::IntegerVector& alt_of_row_1b,
             const Rcpp::IntegerVector& choice_pos_1b,
             const int J_in, const bool include_outside_option) {
    n_tasks = M.size();
    N_persons = Ti.size();
    J = J_in;
    include_outside = include_outside_option;
    if (n_tasks < 1 || N_persons < 1 || J < 1) return false;

    row_offsets = compute_prefix_sum(M);
    task_offsets = compute_prefix_sum(Ti);
    total_rows = row_offsets[n_tasks];
    if (task_offsets[N_persons] != n_tasks) return false;
    if (alt_of_row_1b.size() != total_rows) return false;
    if (choice_pos_1b.size() != n_tasks) return false;

    // Row -> task map and per-alternative counts (fixed row order).
    task_of_row.assign(total_rows, 0);
    for (int t = 0; t < n_tasks; ++t) {
      if (M[t] < 1) return false;
      for (int r = row_offsets[t]; r < row_offsets[t + 1]; ++r) {
        task_of_row[r] = t;
      }
    }
    alt_of_row.assign(total_rows, 0);
    std::vector<int> counts(J, 0);
    for (int r = 0; r < total_rows; ++r) {
      const int a = alt_of_row_1b[r] - 1;
      if (a < 0 || a >= J) return false;
      alt_of_row[r] = a;
      ++counts[a];
    }

    // CSR row lists over alternatives; the fill loop runs in ascending row
    // order, so rows are ascending within each alternative — the fixed
    // summation order for the per-delta_j sufficient statistics.
    alt_row_offsets.assign(J + 1, 0);
    for (int a = 0; a < J; ++a) {
      alt_row_offsets[a + 1] = alt_row_offsets[a] + counts[a];
    }
    alt_rows.assign(total_rows, 0);
    std::vector<int> cursor(alt_row_offsets.begin(), alt_row_offsets.end() - 1);
    for (int r = 0; r < total_rows; ++r) {
      alt_rows[cursor[alt_of_row[r]]++] = r;
    }

    // Chosen-row lookup (0 = outside chosen, valid only with an outside good).
    chosen_row.assign(n_tasks, -1);
    choice_pos.assign(n_tasks, 0);
    for (int t = 0; t < n_tasks; ++t) {
      const int cp = choice_pos_1b[t];
      const int m_t = row_offsets[t + 1] - row_offsets[t];
      if (cp < 0 || cp > m_t) return false;
      if (cp == 0 && !include_outside) return false;
      choice_pos[t] = cp;
      chosen_row[t] = (cp == 0) ? -1 : row_offsets[t] + cp - 1;
    }
    return true;
  }
};

// ----------------------------------------------------------------------------
// Hand-rolled dense Cholesky and triangular solves.
//
// Plain fixed-order loops over arma element storage: no arma::chol,
// arma::solve, LAPACK, or BLAS, so these are callable from OpenMP worker
// threads (per-respondent proposal/regression precisions). Failure (a
// non-positive or non-finite pivot) is reported through the return value —
// never an exception, which must not cross an OpenMP region boundary.
// ----------------------------------------------------------------------------

// A = L L' with L lower triangular. Returns false when A is not (numerically)
// SPD. `!(s > 0)` also catches NaN pivots.
inline bool hb_chol_lower(const arma::mat& A, arma::mat& L) {
  const arma::uword K = A.n_rows;
  if (A.n_cols != K || K == 0) return false;
  L.zeros(K, K);
  for (arma::uword j = 0; j < K; ++j) {
    double s = A(j, j);
    for (arma::uword k = 0; k < j; ++k) s -= L(j, k) * L(j, k);
    if (!(s > 0.0) || !std::isfinite(s)) return false;
    const double ljj = std::sqrt(s);
    L(j, j) = ljj;
    for (arma::uword i = j + 1; i < K; ++i) {
      double t = A(i, j);
      for (arma::uword k = 0; k < j; ++k) t -= L(i, k) * L(j, k);
      L(i, j) = t / ljj;
    }
  }
  return true;
}

// Forward substitution: x = L^{-1} b (L lower triangular).
inline void hb_forward_solve(const arma::mat& L, const arma::vec& b,
                             arma::vec& x) {
  const arma::uword K = L.n_rows;
  x.set_size(K);
  for (arma::uword i = 0; i < K; ++i) {
    double s = b(i);
    for (arma::uword k = 0; k < i; ++k) s -= L(i, k) * x(k);
    x(i) = s / L(i, i);
  }
}

// Back substitution on the transpose: x = L^{-T} b (L lower triangular).
inline void hb_back_solve(const arma::mat& L, const arma::vec& b,
                          arma::vec& x) {
  const int K = static_cast<int>(L.n_rows);
  x.set_size(K);
  for (int i = K - 1; i >= 0; --i) {
    double s = b(i);
    for (int k = i + 1; k < K; ++k) s -= L(k, i) * x(k);
    x(i) = s / L(i, i);
  }
}

// x = (L L')^{-1} b given the Cholesky factor: two trisolves, no inverse and
// no downdate ever formed.
inline void hb_chol_solve(const arma::mat& L, const arma::vec& b,
                          arma::vec& x) {
  arma::vec y;
  hb_forward_solve(L, b, y);
  hb_back_solve(L, y, x);
}

// Solve A x = b for SPD A via hb_chol_lower + two trisolves. Returns false
// (x untouched) when A is not SPD.
inline bool hb_spd_solve(const arma::mat& A, const arma::vec& b,
                         arma::vec& x) {
  arma::mat L;
  if (!hb_chol_lower(A, L)) return false;
  hb_chol_solve(L, b, x);
  return true;
}

// Draw x ~ N(mean, P^{-1}) given the lower Cholesky factor Lprec of the
// PRECISION P = Lprec Lprec': x = mean + Lprec^{-T} z, z ~ N(0, I). Draws
// exclusively from the Xoshiro stream, so it is worker-thread safe with a
// per-task stream.
inline void hb_mvn_precision_draw(Xoshiro256pp& rng, const arma::vec& mean,
                                  const arma::mat& Lprec, arma::vec& out) {
  const arma::uword K = mean.n_elem;
  arma::vec z(K);
  for (arma::uword k = 0; k < K; ++k) z(k) = rng.rnorm();
  arma::vec t;
  hb_back_solve(Lprec, z, t);
  out = mean + t;
}

// ----------------------------------------------------------------------------
// Fixed-order scalar helpers (worker-thread safe).
// ----------------------------------------------------------------------------

// Log-sum-exp over the contiguous range v[0..n) with an optional implicit
// "+ exp(0)" outside-option term (the implicit-outside convention of
// .prepare_hb_panel: the outside good has systematic utility 0 and no
// physical row). Semantics match logSumExp in choicer.h: NaN propagates,
// any +inf gives +inf, an empty range without the outside term gives -inf.
inline double hb_logsumexp(const double* v, const int n,
                           const bool include_outside) {
  double a = include_outside ? 0.0 : -arma::datum::inf;
  for (int i = 0; i < n; ++i) {
    if (std::isnan(v[i])) return arma::datum::nan;
    if (v[i] > a) a = v[i];
  }
  if (a == arma::datum::inf) return arma::datum::inf;
  if (a == -arma::datum::inf) return -arma::datum::inf;  // empty, no outside
  double s = include_outside ? std::exp(0.0 - a) : 0.0;
  for (int i = 0; i < n; ++i) s += std::exp(v[i] - a);
  return a + std::log(s);
}

// Quadratic form x' A x, accumulated column-by-column in fixed order.
inline double hb_quad_form(const arma::mat& A, const arma::vec& x) {
  const arma::uword K = x.n_elem;
  double out = 0.0;
  for (arma::uword j = 0; j < K; ++j) {
    const double* acol = A.colptr(j);
    double s = 0.0;
    for (arma::uword i = 0; i < K; ++i) s += acol[i] * x(i);
    out += s * x(j);
  }
  return out;
}

// Symmetric rank-1 accumulate A += w * x x', column-major fixed order.
inline void hb_rank1_accum(arma::mat& A, const arma::vec& x,
                           const double w = 1.0) {
  const arma::uword K = x.n_elem;
  for (arma::uword j = 0; j < K; ++j) {
    double* acol = A.colptr(j);
    const double xj = w * x(j);
    for (arma::uword i = 0; i < K; ++i) acol[i] += x(i) * xj;
  }
}

// Inverse-gamma(shape, scale): X ~ IG iff 1/X ~ Gamma(shape, rate = scale),
// drawn as scale / rgamma(shape) from the Xoshiro stream — NEVER R's RNG.
inline double hb_rinvgamma(Xoshiro256pp& rng, const double shape,
                           const double scale) {
  return scale / rng.rgamma(shape);
}

// ----------------------------------------------------------------------------
// Master-only conditional draws for the shared hierarchy (b, W, theta,
// sigma_d2). These run on the master thread between barriers, so they may
// use armadillo / bayes_samplers freely; the sufficient statistics are still
// accumulated in fixed i = 0..N-1 / j = 0..J-1 order for bitwise
// thread-count invariance of everything feeding them.
// ----------------------------------------------------------------------------

// b | beta_1..N, W  ~  MVN with precision (A + N W^{-1}) and mean
// (A + N W^{-1})^{-1} (A b_bar + W^{-1} sum_i beta_i).
// beta_i is K x N, one column per respondent. Returns false when the
// posterior precision is not SPD.
inline bool draw_b_conditional(Xoshiro256pp& rng, const arma::mat& beta_i,
                               const arma::vec& b_bar, const arma::mat& A,
                               const arma::mat& Winv, arma::vec& b_out) {
  const arma::uword K = beta_i.n_rows;
  const arma::uword N = beta_i.n_cols;
  arma::vec beta_sum(K, arma::fill::zeros);
  for (arma::uword i = 0; i < N; ++i) {           // fixed respondent order
    const double* bc = beta_i.colptr(i);
    for (arma::uword k = 0; k < K; ++k) beta_sum(k) += bc[k];
  }
  const arma::mat P = A + static_cast<double>(N) * Winv;
  arma::vec rhs(K);
  for (arma::uword k = 0; k < K; ++k) {
    double s = 0.0;
    for (arma::uword l = 0; l < K; ++l) {
      s += A(k, l) * b_bar(l) + Winv(k, l) * beta_sum(l);
    }
    rhs(k) = s;
  }
  arma::mat L;
  if (!hb_chol_lower(P, L)) return false;
  arma::vec mean;
  hb_chol_solve(L, rhs, mean);
  hb_mvn_precision_draw(rng, mean, L, b_out);
  return true;
}

// W | beta_1..N, b  ~  IW(nu + N, V + S), S = sum_i (beta_i - b)(beta_i - b)'
// accumulated in fixed respondent order; drawn via riwishart_nothrow
// (src/bayes_samplers.h:146). Returns false on a non-PD scale matrix.
inline bool draw_W_conditional(Xoshiro256pp& rng, const arma::mat& beta_i,
                               const arma::vec& b, const double nu,
                               const arma::mat& V, arma::mat& W_out) {
  const arma::uword K = beta_i.n_rows;
  const arma::uword N = beta_i.n_cols;
  arma::mat S = V;
  arma::vec d(K);
  for (arma::uword i = 0; i < N; ++i) {           // fixed respondent order
    const double* bc = beta_i.colptr(i);
    for (arma::uword k = 0; k < K; ++k) d(k) = bc[k] - b(k);
    hb_rank1_accum(S, d, 1.0);
  }
  return riwishart_nothrow(rng, nu + static_cast<double>(N), S, W_out);
}

// theta | delta, sigma_d2  ~  MVN — the Bayesian regression of delta on Z
// with known error variance sigma_d2 and prior N(theta_bar, A_theta^{-1}):
// precision A_theta + Z'Z / sigma_d2, mean solve of
// (A_theta theta_bar + Z' delta / sigma_d2). Z is J x P. Returns false when
// the posterior precision is not SPD or sigma_d2 <= 0.
inline bool draw_theta_conditional(Xoshiro256pp& rng, const arma::mat& Z,
                                   const arma::vec& delta,
                                   const double sigma_d2,
                                   const arma::vec& theta_bar,
                                   const arma::mat& A_theta,
                                   arma::vec& theta_out) {
  const arma::uword J = Z.n_rows;
  const arma::uword P = Z.n_cols;
  if (!(sigma_d2 > 0.0) || !std::isfinite(sigma_d2)) return false;
  const double inv_s = 1.0 / sigma_d2;
  arma::mat Prec = A_theta;
  arma::vec rhs(P);
  for (arma::uword p = 0; p < P; ++p) {
    double s = 0.0;
    for (arma::uword l = 0; l < P; ++l) s += A_theta(p, l) * theta_bar(l);
    rhs(p) = s;
  }
  arma::vec zrow(P);
  for (arma::uword j = 0; j < J; ++j) {           // fixed alternative order
    for (arma::uword p = 0; p < P; ++p) zrow(p) = Z(j, p);
    hb_rank1_accum(Prec, zrow, inv_s);
    const double dj = delta(j) * inv_s;
    for (arma::uword p = 0; p < P; ++p) rhs(p) += zrow(p) * dj;
  }
  arma::mat L;
  if (!hb_chol_lower(Prec, L)) return false;
  arma::vec mean;
  hb_chol_solve(L, rhs, mean);
  hb_mvn_precision_draw(rng, mean, L, theta_out);
  return true;
}

// ----------------------------------------------------------------------------
// sigma_d prior: half-Cauchy(0, s_d) by default, implemented as the
// Makalic-Schmidt (2016) conjugate inverse-gamma scale mixture
//
//   sigma_d2 | a_d ~ IG(1/2, 1/a_d),   a_d ~ IG(1/2, 1/s_d^2)
//
// which marginally gives sigma_d ~ half-Cauchy(0, s_d) and avoids the
// IG(eps, eps) near-zero pathology (Gelman 2006). With xi_j ~ N(0, sigma_d2),
// j = 1..J, the two Gibbs blocks are (derivation in the comments below):
//
//   sigma_d2 | a_d, xi ~ IG(1/2 + J/2, 1/a_d + sum_j xi_j^2 / 2)
//   a_d | sigma_d2     ~ IG(1,         1/s_d^2 + 1/sigma_d2)
//
// [p(s2|a) ∝ s2^{-3/2} e^{-1/(a s2)} times the N(0, s2) likelihood
//  ∝ s2^{-J/2} e^{-ss/(2 s2)} gives the first line; collecting the a-terms
//  a^{-3/2} e^{-1/(s_d^2 a)} · a^{-1/2} e^{-1/(a s2)} = a^{-2} e^{-(1/s_d^2 +
//  1/s2)/a} gives the second.]
//
// half_cauchy = false switches to a plain conjugate IG(c0, d0) prior:
// sigma_d2 | xi ~ IG(c0 + J/2, d0 + sum_j xi_j^2 / 2) (a_d untouched).
// ----------------------------------------------------------------------------
struct HbSigmaDPrior {
  bool half_cauchy = true;  // half-Cauchy(0, s_d) scale mixture vs IG(c0, d0)
  double s_d = 1.0;         // half-Cauchy scale
  double c0 = 3.0;          // IG fallback shape
  double d0 = 3.0;          // IG fallback scale
};

// One Gibbs update of (sigma_d2, a_d) given xi. Updates both in place;
// deterministic block order (sigma_d2 first, then a_d). Master-only by
// convention (single hierarchy-level draw), though nothing in it is
// thread-unsafe given its own stream.
inline void draw_sigma_d2_conditional(Xoshiro256pp& rng, const arma::vec& xi,
                                      const HbSigmaDPrior& prior,
                                      double& sigma_d2, double& a_d) {
  const arma::uword J = xi.n_elem;
  double ss = 0.0;
  for (arma::uword j = 0; j < J; ++j) ss += xi(j) * xi(j);  // fixed order
  if (prior.half_cauchy) {
    sigma_d2 = hb_rinvgamma(rng, 0.5 + 0.5 * static_cast<double>(J),
                            1.0 / a_d + 0.5 * ss);
    a_d = hb_rinvgamma(rng, 1.0,
                       1.0 / (prior.s_d * prior.s_d) + 1.0 / sigma_d2);
  } else {
    sigma_d2 = hb_rinvgamma(rng, prior.c0 + 0.5 * static_cast<double>(J),
                            prior.d0 + 0.5 * ss);
  }
}

// ----------------------------------------------------------------------------
// Welford online mean/SD accumulator for the keep_beta_i = "means" and the
// always-on delta_j / xi_j summaries: one update per kept draw, O(1) memory
// in the number of draws. Element-wise over an arbitrary matrix shape
// (K x N for beta_i, J x 1 for delta). Master-only (updated at recording
// time between barriers).
// ----------------------------------------------------------------------------
struct HbWelford {
  long long n = 0;
  arma::mat mean;   // running mean
  arma::mat m2;     // running sum of squared deviations

  void reset(const arma::uword n_rows, const arma::uword n_cols) {
    n = 0;
    mean.zeros(n_rows, n_cols);
    m2.zeros(n_rows, n_cols);
  }

  void update(const arma::mat& x) {
    ++n;
    const arma::uword nel = x.n_elem;
    const double* xv = x.memptr();
    double* mv = mean.memptr();
    double* sv = m2.memptr();
    const double inv_n = 1.0 / static_cast<double>(n);
    for (arma::uword t = 0; t < nel; ++t) {
      const double d = xv[t] - mv[t];
      mv[t] += d * inv_n;
      sv[t] += d * (xv[t] - mv[t]);
    }
  }

  // Sample SD (n - 1 denominator); all-zero until two updates arrive.
  arma::mat sd() const {
    if (n < 2) return arma::mat(arma::size(mean), arma::fill::zeros);
    return arma::sqrt(m2 / static_cast<double>(n - 1));
  }
};

// ----------------------------------------------------------------------------
// Interrupt poller (clone of src/mnprobit.cpp:40-47). True if a user
// interrupt is pending. Safe inside the parallel region as long as only the
// master (R's main) thread calls it.
// ----------------------------------------------------------------------------
inline void hb_chk_int_fn(void*) { R_CheckUserInterrupt(); }

inline bool hb_pending_interrupt() {
  return R_ToplevelExec(hb_chk_int_fn, nullptr) == FALSE;
}

#endif
