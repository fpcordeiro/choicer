#ifndef CHOICER_RNG_HPP
#define CHOICER_RNG_HPP

#include <cstdint>
#include <cmath>

// ============================================================================
// Self-contained RNG core for the Bayesian (MCMC) engines.
//
// R's RNG is single-threaded only, so the Gibbs samplers use their own
// generator: xoshiro256++ (Blackman & Vigna 2021) seeded via splitmix64, the
// initialisation the xoshiro authors recommend. All scalar samplers are exact
// (no approximate densities): uniform (53-bit), normal (Marsaglia polar),
// exponential (inversion), gamma (Marsaglia & Tsang 2000), chi-squared.
//
// Reproducibility / thread-safety contract:
//   * Each consumer derives an independent stream with make_stream(seed,
//     iter, tag) — for the MNP latent sweep, one stream per (iteration,
//     observation). Streams are deterministic functions of their key, so
//     results are bitwise reproducible regardless of the number of OpenMP
//     threads or the loop schedule.
//   * An Xoshiro256pp instance must never be shared across threads; create
//     one per task with make_stream() instead.
//   * The master seed is supplied from R (drawn from R's RNG when the user
//     does not set one), so set.seed() governs reproducibility end to end.
// ============================================================================

// Seeding / stream-derivation primitive (Vigna's splitmix64).
inline uint64_t splitmix64_next(uint64_t& s) {
  uint64_t z = (s += 0x9E3779B97F4A7C15ULL);
  z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
  z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
  return z ^ (z >> 31);
}

struct Xoshiro256pp {
  uint64_t s[4];
  bool has_spare;     // cached second draw from the polar method
  double spare;

  explicit Xoshiro256pp(uint64_t seed) : has_spare(false), spare(0.0) {
    for (int i = 0; i < 4; ++i) s[i] = splitmix64_next(seed);
  }

  static uint64_t rotl(const uint64_t x, const int k) {
    return (x << k) | (x >> (64 - k));
  }

  uint64_t next() {
    const uint64_t result = rotl(s[0] + s[3], 23) + s[0];
    const uint64_t t = s[1] << 17;
    s[2] ^= s[0];
    s[3] ^= s[1];
    s[1] ^= s[2];
    s[0] ^= s[3];
    s[2] ^= t;
    s[3] = rotl(s[3], 45);
    return result;
  }

  // Uniform on the open interval (0, 1): 53-bit resolution, never 0 or 1,
  // so log(runif()) and qnorm(runif()) are always finite.
  double runif() {
    return (static_cast<double>(next() >> 11) + 0.5) * (1.0 / 9007199254740992.0);
  }

  // Standard normal via the Marsaglia polar method (exact); generates pairs
  // and caches the spare.
  double rnorm() {
    if (has_spare) {
      has_spare = false;
      return spare;
    }
    double u, v, q;
    do {
      u = 2.0 * runif() - 1.0;
      v = 2.0 * runif() - 1.0;
      q = u * u + v * v;
    } while (q >= 1.0 || q == 0.0);
    const double f = std::sqrt(-2.0 * std::log(q) / q);
    spare = v * f;
    has_spare = true;
    return u * f;
  }

  // Exponential(rate 1) via inversion.
  double rexp() {
    return -std::log(runif());
  }

  // Gamma(shape a, scale 1) via Marsaglia & Tsang (2000); exact for any
  // a > 0 (the a < 1 case uses the boosting identity X = Gamma(a+1) U^{1/a}).
  double rgamma(double a) {
    double boost = 1.0;
    if (a < 1.0) {
      boost = std::pow(runif(), 1.0 / a);
      a += 1.0;
    }
    const double d = a - 1.0 / 3.0;
    const double c = 1.0 / std::sqrt(9.0 * d);
    for (;;) {
      double x, v;
      do {
        x = rnorm();
        v = 1.0 + c * x;
      } while (v <= 0.0);
      v = v * v * v;
      const double u = runif();
      if (u < 1.0 - 0.0331 * x * x * x * x) return boost * d * v;
      if (std::log(u) < 0.5 * x * x + d * (1.0 - v + std::log(v))) return boost * d * v;
    }
  }

  // Chi-squared(df) = 2 * Gamma(df / 2, 1); df > 0, non-integer allowed.
  double rchisq(double df) {
    return 2.0 * rgamma(0.5 * df);
  }
};

// Independent stream for task (iter, tag) under the master seed. The three
// keys are folded through splitmix64 so distinct keys give decorrelated
// streams; the Xoshiro256pp constructor then expands the mixed value into
// the full 256-bit state.
inline uint64_t mix_seed(uint64_t a, const uint64_t b) {
  a += 0x9E3779B97F4A7C15ULL * (b + 1);
  return splitmix64_next(a);
}

inline Xoshiro256pp make_stream(const uint64_t seed, const uint64_t iter,
                                const uint64_t tag) {
  return Xoshiro256pp(mix_seed(mix_seed(seed, iter), tag));
}

#endif
