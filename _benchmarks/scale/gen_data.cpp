// Synthetic population-scale mixed logit data in choicer's prepared (stacked)
// layout, generated directly (no prepare_mxl_data), deterministic per unit.
//
// Unit u (decision maker) belongs to region r(u) with hospital/school set
// A_r (size log-uniform in [m_min, m_max], drawn from J_total alternatives).
// Each of its T_u choice situations faces A_r, each alternative kept with
// probability keep_p (at least 2 kept). T_u:
//   "pareto": T = min(T_max, floor(t_min * U^(-1/alpha)))   heavy tail
//   "short" : T = 1 with prob p_single, else 2 or 3
// Utility of row r in situation t of unit u:
//   V = X_r beta + W_r (mu + L eta_u) + delta_alt(r) + Gumbel
// X columns are N(0,1), except column 1 which is scaled Student-t(nu) when
// x1_df > 0 (fat tails). W columns are N(0,1). The choice is the argmax.
// Alternatives are relabelled 1..J in order of first appearance sorted, so
// alt_idx is dense and every present alternative has an index.
// [[Rcpp::plugins(openmp)]]
#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <numeric>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

struct Rng {
  uint64_t s;
  explicit Rng(uint64_t seed) : s(seed) {}
  uint64_t next() {  // splitmix64
    uint64_t z = (s += 0x9e3779b97f4a7c15ULL);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
  }
  double unif() { return ((next() >> 11) + 0.5) * 0x1.0p-53; }  // (0, 1)
  double norm() {
    const double u1 = unif(), u2 = unif();
    return std::sqrt(-2.0 * std::log(u1)) * std::cos(6.283185307179586 * u2);
  }
  double chisq(int df) {
    double s = 0.0;
    for (int i = 0; i < df; ++i) { const double z = norm(); s += z * z; }
    return s;
  }
};

static uint64_t mix(uint64_t a, uint64_t b) {
  Rng r(a * 0x100000001b3ULL ^ (b + 0x632be59bd9b4e019ULL));
  r.next();
  return r.next();
}

// [[Rcpp::export]]
Rcpp::List gen_mxl_data(int U, int K_x, int K_w, int J_total, int n_regions,
                        int m_min, int m_max, double keep_p,
                        std::string t_dist, double alpha, int t_min,
                        int T_max, double p_single,
                        Rcpp::NumericVector beta, Rcpp::NumericVector mu,
                        Rcpp::NumericMatrix L, double delta_sd,
                        double x1_df, double x1_scale, double seed_d,
                        int threads) {
  const uint64_t seed = static_cast<uint64_t>(seed_d);
  // --- regions: alternative sets --------------------------------------------
  std::vector<std::vector<int>> region(n_regions);
  {
    Rng rng(mix(seed, 1));
    std::vector<int> pool(J_total);
    for (int r = 0; r < n_regions; ++r) {
      const double lm = std::log((double)m_min) +
                        rng.unif() * (std::log((double)m_max) - std::log((double)m_min));
      int m = std::min(J_total, std::max(2, (int)std::lround(std::exp(lm))));
      std::iota(pool.begin(), pool.end(), 0);
      for (int k = 0; k < m; ++k) {  // partial Fisher-Yates
        const int j = k + (int)(rng.unif() * (J_total - k));
        std::swap(pool[k], pool[std::min(j, J_total - 1)]);
      }
      region[r].assign(pool.begin(), pool.begin() + m);
      std::sort(region[r].begin(), region[r].end());
    }
  }
  // --- units: T_u and region --------------------------------------------------
  std::vector<int> Tu(U), ru(U);
  {
    Rng rng(mix(seed, 2));
    for (int u = 0; u < U; ++u) {
      int T;
      if (t_dist == "pareto") {
        T = (int)std::floor(t_min * std::pow(rng.unif(), -1.0 / alpha));
        T = std::max(1, std::min(T_max, T));
      } else {
        const double v = rng.unif();
        T = v < p_single ? 1 : (v < p_single + (1 - p_single) / 2 ? 2 : 3);
      }
      Tu[u] = T;
      ru[u] = (int)(rng.unif() * n_regions) % n_regions;
    }
  }
  std::vector<int64_t> sit_off(U + 1, 0);
  for (int u = 0; u < U; ++u) sit_off[u + 1] = sit_off[u] + Tu[u];
  const int64_t N = sit_off[U];
  if (N > 2147483647LL) Rcpp::stop("too many situations");
  // --- choice sets: kept alternatives per situation (serial, per-unit RNG) ---
  std::vector<int> M(N);
  std::vector<uint64_t> keep_mask_seed(N);
  for (int u = 0; u < U; ++u) {
    const int m_r = (int)region[ru[u]].size();
    for (int64_t t = sit_off[u]; t < sit_off[u + 1]; ++t) {
      const uint64_t s = mix(seed ^ 0xabcdefULL, (uint64_t)t);
      keep_mask_seed[t] = s;
      if (keep_p >= 1.0) { M[t] = m_r; continue; }
      Rng rng(s);
      int kept = 0;
      for (int k = 0; k < m_r; ++k) kept += rng.unif() < keep_p;
      M[t] = std::max(2, kept);
    }
  }
  std::vector<int64_t> row_off(N + 1, 0);
  for (int64_t t = 0; t < N; ++t) row_off[t + 1] = row_off[t] + M[t];
  const int64_t R = row_off[N];
  if (R > 2147483647LL) Rcpp::stop("too many rows for an R matrix index here");
  // --- alternative relabelling: present alternatives -> 1..J ------------------
  std::vector<int> present(J_total, 0);
  std::vector<char> region_used(n_regions, 0);
  for (int u = 0; u < U; ++u) region_used[ru[u]] = 1;
  for (int r = 0; r < n_regions; ++r)
    if (region_used[r]) for (int j : region[r]) present[j] = 1;
  std::vector<int> label(J_total, 0);
  int J = 0;
  for (int j = 0; j < J_total; ++j) if (present[j]) label[j] = ++J;
  // ASCs of the present alternatives (delta_1 = 0)
  std::vector<double> delta(J + 1, 0.0);
  {
    Rng rng(mix(seed, 3));
    for (int j = 2; j <= J; ++j) delta[j] = delta_sd * rng.norm();
  }
  // --- fill ---------------------------------------------------------------
  Rcpp::NumericMatrix X((int)R, K_x), W((int)R, K_w);
  Rcpp::IntegerVector alt_idx((int)R), choice_idx((int)N), Mv((int)N),
      Tv(U);
  Rcpp::NumericMatrix gamma_true(K_w, U);
  double* Xp = X.begin();
  double* Wp = W.begin();
  int* ap = alt_idx.begin();
  int* cp = choice_idx.begin();
  double* gp = gamma_true.begin();
  for (int64_t t = 0; t < N; ++t) Mv[t] = M[t];
  for (int u = 0; u < U; ++u) Tv[u] = Tu[u];
  const double* bp = beta.begin();
  const double* mp = mu.begin();
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 64) num_threads(threads)
#endif
  for (int u = 0; u < U; ++u) {
    Rng rng(mix(seed ^ 0x51ed2701ULL, (uint64_t)u));
    std::vector<double> g(K_w), eta(K_w);
    for (int k = 0; k < K_w; ++k) eta[k] = rng.norm();
    for (int k = 0; k < K_w; ++k) {
      double s = mp[k];
      for (int q = 0; q <= k; ++q) s += L(k, q) * eta[q];
      g[k] = s;
      gp[(int64_t)u * K_w + k] = s;
    }
    const std::vector<int>& A = region[ru[u]];
    const int m_r = (int)A.size();
    for (int64_t t = sit_off[u]; t < sit_off[u + 1]; ++t) {
      // kept alternatives
      std::vector<int> kept;
      kept.reserve(m_r);
      if (keep_p >= 1.0) {
        kept = A;
      } else {
        Rng kr(keep_mask_seed[t]);
        for (int k = 0; k < m_r; ++k) if (kr.unif() < keep_p) kept.push_back(A[k]);
        for (int k = 0; (int)kept.size() < 2 && k < m_r; ++k)
          if (std::find(kept.begin(), kept.end(), A[k]) == kept.end()) kept.push_back(A[k]);
        std::sort(kept.begin(), kept.end());
      }
      const int64_t r0 = row_off[t];
      int best = 0;
      double best_v = -INFINITY;
      for (int a = 0; a < (int)kept.size(); ++a) {
        const int64_t r = r0 + a;
        double v = delta[label[kept[a]]];
        for (int k = 0; k < K_x; ++k) {
          double x = rng.norm();
          if (k == 0 && x1_df > 0) x = x1_scale * x / std::sqrt(rng.chisq((int)x1_df) / x1_df);
          Xp[r + (int64_t)k * R] = x;
          v += bp[k] * x;
        }
        for (int k = 0; k < K_w; ++k) {
          const double w = rng.norm();
          Wp[r + (int64_t)k * R] = w;
          v += g[k] * w;
        }
        v += -std::log(-std::log(rng.unif()));  // Gumbel
        ap[r] = label[kept[a]];
        if (v > best_v) { best_v = v; best = a; }
      }
      cp[t] = best + 1;
    }
  }
  Rcpp::CharacterVector xn(K_x), wn(K_w);
  for (int k = 0; k < K_x; ++k) xn[k] = "x" + std::to_string(k + 1);
  for (int k = 0; k < K_w; ++k) wn[k] = "w" + std::to_string(k + 1);
  Rcpp::colnames(X) = xn;
  Rcpp::colnames(W) = wn;
  Rcpp::NumericVector delta_out(J);
  for (int j = 0; j < J; ++j) delta_out[j] = delta[j + 1];
  return Rcpp::List::create(
      Rcpp::Named("X") = X, Rcpp::Named("W") = W,
      Rcpp::Named("alt_idx") = alt_idx, Rcpp::Named("choice_idx") = choice_idx,
      Rcpp::Named("M") = Mv, Rcpp::Named("Ti") = Tv,
      Rcpp::Named("J") = J, Rcpp::Named("N") = (double)N,
      Rcpp::Named("rows") = (double)R, Rcpp::Named("delta") = delta_out,
      Rcpp::Named("gamma") = gamma_true);
}
