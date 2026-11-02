// ============================================================================
// The geographic generator's C++ port (plan §4: D3 hospital and D4 school at
// population scale, 5 x 10^7 rows). lib/dgp_geo.R is the reference and its
// header documents the designs; this file holds only the parts whose cost
// grows with the long table: the choice sets (D4: each unit's K nearest
// alternatives of its market; D3: the alternatives within a radius, at
// least K_min and at most K_max of the nearest, by a grid search), each
// set's rows in choicer's order, and the logit probabilities and choices.
// Randomness stays in R (geo_draw()): this file draws nothing, so at small
// sizes it reproduces the R reference exactly (identical sets, distances and
// choices; probabilities to rounding: X beta is BLAS's dot product in R).
//
// Arithmetic as R's: a distance is sqrt(dx * dx + dy * dy), each product
// and sum rounded on its own (R evaluates them as separate vector
// operations), so floating-point contraction into fused multiply-adds must
// be off. The one guard is the compile flag: geo_gen_compile() passes
// -ffp-contract=off (a pragma does not override -ffp-contract=fast in Apple
// clang 17), and since a user's Makevars can override it, geo_gen_compile()
// checks geo_knn()'s distances against R's bit for bit after every load and
// stops if they differ. Compiled by geo_gen_compile() (lib/dgp_geo.R) with
// Rcpp::sourceCpp(). Single-threaded: generation is a small share of a
// replication's time.
// ============================================================================

#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <utility>
#include <vector>

using Rcpp::IntegerVector;
using Rcpp::List;
using Rcpp::NumericMatrix;
using Rcpp::NumericVector;

typedef std::pair<double, int> DistIdx;   // (distance, 0-based alternative)

// The Euclidean distance, rounded as R rounds it.
static inline double geo_dist(double x1, double y1, double x2, double y2) {
  const double dx = x1 - x2;
  const double dy = y1 - y2;
  const double dxx = dx * dx;
  const double dyy = dy * dy;
  const double s = dxx + dyy;
  return std::sqrt(s);
}

static bool by_index(const DistIdx& a, const DistIdx& b) { return a.second < b.second; }

// D4: each unit's K nearest alternatives of its own market (ties by index),
// rows in (unit, alternative index) order. Units and alternatives carry
// their market (1..M); units are taken in their order, which is id order.
// Returns list(id, j, dist): the unit (1-based), the alternative (1-based)
// and the distance, one row per unit and set member.
// [[Rcpp::export]]
List geo_knn(NumericVector ux, NumericVector uy, IntegerVector um, NumericVector ax,
             NumericVector ay, IntegerVector am, int K) {
  const R_xlen_t n = ux.size();
  const int J = ax.size();
  if (uy.size() != n || um.size() != n || ay.size() != J || am.size() != J) {
    Rcpp::stop("geo_knn(): coordinates and markets differ in length");
  }
  if (K < 1) Rcpp::stop("geo_knn(): K must be positive");
  int M = 0;
  for (int j = 0; j < J; ++j) {
    if (am[j] < 1) Rcpp::stop("geo_knn(): markets are coded 1..M");
    M = std::max(M, am[j]);
  }
  std::vector<std::vector<int>> members(M + 1);
  for (int j = 0; j < J; ++j) members[am[j]].push_back(j);
  IntegerVector id(n * K), jj(n * K);
  NumericVector dd(n * K);
  std::vector<DistIdx> cand;
  R_xlen_t row = 0;
  for (R_xlen_t i = 0; i < n; ++i) {
    const int m = um[i];
    if (m < 1 || m > M) Rcpp::stop("geo_knn(): a unit's market holds no alternatives");
    const std::vector<int>& sc = members[m];
    if ((int)sc.size() < K) {
      Rcpp::stop("market %d has %d schools, fewer than K = %d", m, (int)sc.size(), K);
    }
    cand.clear();
    for (int j : sc) cand.push_back(DistIdx(geo_dist(ux[i], uy[i], ax[j], ay[j]), j));
    std::partial_sort(cand.begin(), cand.begin() + K, cand.end());
    std::sort(cand.begin(), cand.begin() + K, by_index);
    for (int k = 0; k < K; ++k, ++row) {
      id[row] = (int)(i + 1);
      jj[row] = cand[k].second + 1;
      dd[row] = cand[k].first;
    }
  }
  return List::create(Rcpp::_["id"] = id, Rcpp::_["j"] = jj, Rcpp::_["dist"] = dd);
}

// D3: each patient's set, the alternatives within distance r, at least the
// K_min and at most the K_max nearest (ties by index), repeated over the
// patient's T admissions; rows in (patient, admission, alternative index)
// order, admissions numbered 1..sum(T) patient by patient. The alternatives
// lie in [0, side]^2, bucketed in square cells of side r / 2; a patient's
// search scans rings of cells around its own until every alternative within
// r is found and, if fewer than K_min are, until the K_min nearest are.
// After rings 0..k every alternative within k * cell of the patient has
// been scanned in exact arithmetic, but floor(v / cell) can round a point
// on a cell edge into the next cell (lattice points tied at 1.4 - 1 ulp with
// cell 0.7), so the search counts only (k - 1) * cell as covered: one ring
// of slack. Returns list(id, patient, j, dist).
// [[Rcpp::export]]
List geo_radius(NumericVector ux, NumericVector uy, IntegerVector T, NumericVector ax,
                NumericVector ay, double side, double r, int K_min, int K_max) {
  const R_xlen_t n = ux.size();
  const int J = ax.size();
  if (uy.size() != n || T.size() != n || ay.size() != J) {
    Rcpp::stop("geo_radius(): lengths differ");
  }
  if (!(r > 0) || !(side > 0) || K_min < 1 || K_max < K_min) {
    Rcpp::stop("geo_radius(): needs r > 0, side > 0 and 1 <= K_min <= K_max");
  }
  const double cell = r / 2.0;
  const int G = std::max(1, (int)std::ceil(side / cell));
  auto cell_of = [&](double v) {
    int c = (int)std::floor(v / cell);
    return std::min(G - 1, std::max(0, c));
  };
  std::vector<std::vector<int>> grid((size_t)G * G);
  for (int j = 0; j < J; ++j) grid[(size_t)cell_of(ay[j]) * G + cell_of(ax[j])].push_back(j);
  const int k_min = std::min(K_min, J);
  // Pass 1: each patient's set, kept compactly (start offsets into flat arrays).
  std::vector<R_xlen_t> start(n + 1, 0);
  std::vector<int> set_j;
  std::vector<double> set_d;
  std::vector<DistIdx> cand;
  R_xlen_t A = 0, rows = 0;
  for (R_xlen_t i = 0; i < n; ++i) {
    if (T[i] < 1) Rcpp::stop("geo_radius(): every patient has at least one admission");
    const double x = ux[i], y = uy[i];
    const int cx = cell_of(x), cy = cell_of(y);
    cand.clear();
    for (int k = 0;; ++k) {
      // Ring k: the cells at Chebyshev distance k from the patient's.
      for (int iy = cy - k; iy <= cy + k; ++iy) {
        if (iy < 0 || iy >= G) continue;
        const int step = (iy == cy - k || iy == cy + k) ? 1 : 2 * k;   // the ring's rows
        for (int ix = cx - k; ix <= cx + k; ix += step) {
          if (ix < 0 || ix >= G) continue;
          for (int j : grid[(size_t)iy * G + ix]) {
            cand.push_back(DistIdx(geo_dist(x, y, ax[j], ay[j]), j));
          }
        }
      }
      if (k >= G - 1) break;                               // every cell scanned
      const double certain = (k - 1) * cell;               // one ring of slack
      if (certain < r) continue;
      int n_r = 0, n_c = 0;
      for (const DistIdx& c : cand) {
        n_r += (c.first <= r);
        n_c += (c.first <= certain);
      }
      if (n_r >= k_min || n_c >= k_min) break;
    }
    std::sort(cand.begin(), cand.end());
    int n_r = 0;
    for (const DistIdx& c : cand) n_r += (c.first <= r);
    const int take = std::min({K_max, std::max(k_min, n_r), (int)cand.size()});
    std::sort(cand.begin(), cand.begin() + take, by_index);
    for (int k = 0; k < take; ++k) {
      set_j.push_back(cand[k].second);
      set_d.push_back(cand[k].first);
    }
    start[i + 1] = start[i] + take;
    A += T[i];
    rows += (R_xlen_t)take * T[i];
  }
  // Pass 2: the long rows, each set repeated over the patient's admissions.
  IntegerVector id(rows), patient(rows), jj(rows);
  NumericVector dd(rows);
  R_xlen_t row = 0;
  int adm = 0;
  for (R_xlen_t i = 0; i < n; ++i) {
    for (int t = 0; t < T[i]; ++t) {
      ++adm;
      for (R_xlen_t s = start[i]; s < start[i + 1]; ++s, ++row) {
        id[row] = adm;
        patient[row] = (int)(i + 1);
        jj[row] = set_j[s] + 1;
        dd[row] = set_d[s];
      }
    }
  }
  return List::create(Rcpp::_["id"] = id, Rcpp::_["patient"] = patient, Rcpp::_["j"] = jj,
                      Rcpp::_["dist"] = dd);
}

// The logit probabilities the choices are drawn from and the choices, as
// lib/simulate.R's mnl_probs() / mxl_cond_probs() and draw_choice(u =).
// Rows of a situation are contiguous (`id` changes between situations);
// `code` is each row's alternative code (1..J, indexing `delta`); `X` and
// `W` are lists of columns (`beta` the coefficients of X; `B` a U x K_w
// matrix of each unit's realized random coefficients, `unit` each
// situation's unit, 1..U; W empty for the MNL); `u` one uniform per
// situation. Utility X beta + delta[code] (+ W b_unit); the outside option,
// if any, has utility 0. The X and W columns are read in place, so they
// must be double (a converted copy would not outlive its scope). Returns
// list(p (per row, only with keep_p; else empty), p0 (per situation),
// choice (0/1 per row), expected (per code: the sum of p)).
// [[Rcpp::export]]
List geo_choose(IntegerVector id, IntegerVector code, List X, NumericVector beta,
                NumericVector delta, List W, NumericMatrix B, IntegerVector unit, bool outside,
                NumericVector u, bool keep_p) {
  const R_xlen_t n = id.size();
  const int J = delta.size();
  const int K_x = X.size(), K_w = W.size();
  if (code.size() != n) Rcpp::stop("geo_choose(): `code` needs one value per row");
  if (beta.size() != K_x) Rcpp::stop("geo_choose(): `beta` needs one value per X column");
  // A column's data, read in place: double and one value per row.
  auto column = [&](List cols, int k, const char* what) -> const double* {
    SEXP s = cols[k];
    if (TYPEOF(s) != REALSXP) Rcpp::stop("geo_choose(): %s column %d is not double", what, k + 1);
    if (XLENGTH(s) != n) Rcpp::stop("geo_choose(): %s column %d has the wrong length", what, k + 1);
    return REAL(s);
  };
  std::vector<const double*> xc(K_x), wc(K_w);
  for (int k = 0; k < K_x; ++k) xc[k] = column(X, k, "an X");
  for (int k = 0; k < K_w; ++k) wc[k] = column(W, k, "a W");
  if (K_w > 0 && B.ncol() != K_w) Rcpp::stop("geo_choose(): `B` needs one column per W column");
  // Situations: runs of equal id.
  std::vector<R_xlen_t> first;
  for (R_xlen_t r = 0; r < n; ++r) {
    if (r == 0 || id[r] != id[r - 1]) first.push_back(r);
  }
  const R_xlen_t N = (R_xlen_t)first.size();
  first.push_back(n);
  if (u.size() != N) Rcpp::stop("geo_choose(): `u` needs one uniform per situation");
  if (K_w > 0 && unit.size() != N) Rcpp::stop("geo_choose(): `unit` needs one value per situation");
  NumericVector p(keep_p ? n : 0), p0(N), expected(J);
  IntegerVector choice(n);
  std::vector<double> V, cum;
  for (R_xlen_t s = 0; s < N; ++s) {
    const R_xlen_t a = first[s], b = first[s + 1];
    const R_xlen_t m_s = b - a;
    V.assign(m_s, 0.0);
    const int un = K_w > 0 ? unit[s] - 1 : 0;
    if (K_w > 0 && (un < 0 || un >= B.nrow())) Rcpp::stop("geo_choose(): `unit` out of range");
    double vmax = R_NegInf;
    for (R_xlen_t r = a; r < b; ++r) {
      const int c = code[r];
      if (c < 1 || c > J) Rcpp::stop("geo_choose(): `code` out of range");
      double xb = 0.0;
      for (int k = 0; k < K_x; ++k) {
        const double t = xc[k][r] * beta[k];
        xb = (k == 0) ? t : xb + t;
      }
      double v = xb + delta[c - 1];
      if (K_w > 0) {
        double rp = 0.0;
        for (int k = 0; k < K_w; ++k) {
          const double t = wc[k][r] * B(un, k);
          rp = (k == 0) ? t : rp + t;
        }
        v = v + rp;
      }
      V[r - a] = v;
      if (v > vmax) vmax = v;
    }
    // Log-sum-exp as lib/simulate.R's .group_lse(): the outside option a
    // term exp(0), added last.
    double m = vmax;
    if (outside && m < 0) m = 0;
    if (!std::isfinite(m)) m = 0;
    double sum = 0.0;
    for (R_xlen_t k = 0; k < m_s; ++k) sum += std::exp(V[k] - m);
    if (outside) sum += std::exp(-m);
    const double D = m + std::log(sum);
    // Probabilities and the choice by inverse CDF (draw_choice()).
    cum.assign(m_s, 0.0);
    R_xlen_t last = 0;
    for (R_xlen_t k = 0; k < m_s; ++k) {
      const double pk = std::exp(V[k] - D);
      if (keep_p) p[a + k] = pk;
      expected[code[a + k] - 1] += pk;
      cum[k] = (k == 0) ? pk : cum[k - 1] + pk;
      if (pk > 0) last = k + 1;
    }
    const double q0 = outside ? std::exp(-D) : 0.0;
    p0[s] = q0;
    const double total = q0 + cum[m_s - 1];
    if (!(total > 0)) Rcpp::stop("A situation has no positive probability.");
    const double thr = u[s] * total - q0;
    if (thr >= 0) {
      R_xlen_t k = 0;
      for (R_xlen_t q = 0; q < m_s; ++q) k += (cum[q] <= thr);
      k = std::min(k + 1, last);
      choice[a + k - 1] = 1;
    }
  }
  return List::create(Rcpp::_["p"] = p, Rcpp::_["p0"] = p0, Rcpp::_["choice"] = choice,
                      Rcpp::_["expected"] = expected);
}
