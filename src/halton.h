// halton.h — Header-only on-the-fly Halton draw generator
//
// Provides:
//   HALTON_PRIMES[]       — 128-prime table (primes[k] is the base for dimension k)
//   radical_inverse()     — van der Corput radical-inverse function
//   inv_normal_cdf()      — Wichura (1988) AS241 PPND16 inverse normal CDF
//   struct HaltonGen      — Owen digit-scrambled Halton generator with fill_eta_i()
//
// DESIGN CONTRACT:
//   - Header-only: all functions are inline or templated; no .cpp counterpart.
//   - NO R API: does not include <Rcpp.h>, <R.h>, <Rmath.h>, or call Rf_* functions.
//     Safe to use inside OpenMP parallel regions.
//   - Include order: this header must be included from a .cpp translation unit that
//     has already included choicer.h (which brings in RcppArmadillo.h), so arma::
//     types are available here even though halton.h does not include them directly.
//   - Thread safety: HaltonGen is const after construction; fill_eta_i() takes a
//     mutable eta_i buffer owned by the calling thread. Multiple threads may call
//     fill_eta_i() simultaneously on the same const HaltonGen without data races.
//   - Bitwise reproducibility: n = (i-1)*S + s + 1 is a deterministic function of
//     (i, s, S) and the permutation table is a deterministic function of the seed,
//     so results are identical regardless of OpenMP thread count or schedule.

#ifndef CHOICER_HALTON_HPP
#define CHOICER_HALTON_HPP

#include <cstdint>
#include <cmath>
#include <vector>
#include "rng.h"   // for splitmix64_next(uint64_t&) and mix_seed(uint64_t, uint64_t)

// ============================================================================
// §2.2 Primes table — 128 primes; dimension k uses HALTON_PRIMES[k]
// ============================================================================

static constexpr uint32_t HALTON_PRIMES[] = {
  2,   3,   5,   7,  11,  13,  17,  19,  23,  29,
 31,  37,  41,  43,  47,  53,  59,  61,  67,  71,
 73,  79,  83,  89,  97, 101, 103, 107, 109, 113,
127, 131, 137, 139, 149, 151, 157, 163, 167, 173,
179, 181, 191, 193, 197, 199, 211, 223, 227, 229,
233, 239, 241, 251, 257, 263, 269, 271, 277, 281,
283, 293, 307, 311, 313, 317, 331, 337, 347, 349,
353, 359, 367, 373, 379, 383, 389, 397, 401, 409,
419, 421, 431, 433, 439, 443, 449, 457, 461, 463,
467, 479, 487, 491, 499, 503, 509, 521, 523, 541,
547, 557, 563, 569, 571, 577, 587, 593, 599, 601,
607, 613, 617, 619, 631, 641, 643, 647, 653, 659,
661, 673, 677, 683, 691, 701, 709, 719
};

static constexpr int HALTON_N_PRIMES =
    static_cast<int>(sizeof(HALTON_PRIMES) / sizeof(HALTON_PRIMES[0])); // 128

// ============================================================================
// §2.3 radical_inverse — van der Corput function
//
// Returns the base-b radical inverse of n in [0, 1).
// radical_inverse(0, b) = 0.0
// radical_inverse(1, 2) = 0.5  (matches randtoolbox row 1 with start=1)
// ============================================================================

inline double radical_inverse(uint64_t n, uint32_t base) {
    double result = 0.0;
    double f = 1.0 / static_cast<double>(base);
    while (n > 0) {
        result += static_cast<double>(n % base) * f;
        n /= base;
        f /= static_cast<double>(base);
    }
    return result;
}

// ============================================================================
// §2.4 inv_normal_cdf — Wichura (1988) Algorithm AS241 PPND16
//
// Standalone C++ inverse normal CDF; no R API calls.
// Achieves max absolute error < 1e-12 vs stats::qnorm on p in [1e-10, 1-1e-10].
// (Measured max error on 1600-point dense grid: 8e-15.)
//
// Coefficients verified against R 4.6.0 nmath/qnorm.o binary (ARM64 LE) and
// against a dense grid of stats::qnorm() reference values.
//
// Edge cases: p <= 0 returns -8.29; p >= 1 returns +8.29 (sentinel, never NaN).
// ============================================================================

inline double inv_normal_cdf(double p) {
    // Algorithm split constants
    static const double SPLIT1 = 0.425;    // |q| <= SPLIT1 -> central region
    static const double SPLIT2 = 5.0;      // r  <= SPLIT2  -> intermediate tail
    static const double CONST1 = 0.180625; // = SPLIT1^2; central poly shift
    static const double CONST2 = 1.6;      // r -= CONST2 before C/D evaluation

    // Central region: |p - 0.5| <= SPLIT1
    // Rational poly in (p-0.5)^2; A[0] is constant term, A[7] is degree-7 term.
    static const double A[8] = {
        3.38713287279637e+00,  1.33141667891784e+02,
        1.97159095030655e+03,  1.37316937655095e+04,
        4.59219539315499e+04,  6.72657709270087e+04,
        3.34305755835881e+04,  2.50908092873012e+03
    };
    static const double B[8] = {
        1.0,                   4.23133307016009e+01,
        6.87187007492058e+02,  5.39419602142475e+03,
        2.12137943015866e+04,  3.93078958000927e+04,
        2.87290857357219e+04,  5.22649527885285e+03
    };

    // Intermediate tail: sqrt(-log(min(p,1-p))) - CONST2 in [0, SPLIT2-CONST2]
    static const double C[8] = {
        1.42343711074968e+00,  4.63033784615655e+00,
        5.76949722146069e+00,  3.64784832476320e+00,
        1.27045825245237e+00,  2.41780725177451e-01,
        2.27238449892692e-02,  7.74545014278341e-04
    };
    static const double D[8] = {
        1.0,                   2.05319162663776e+00,
        1.67638483018380e+00,  6.89767334985100e-01,
        1.48103976427480e-01,  1.51986665636165e-02,
        5.47593808499535e-04,  1.05075007164442e-09
    };

    // Far tail: sqrt(-log(min(p,1-p))) - SPLIT2 > 0
    static const double E[8] = {
        6.65790464350110e+00,  5.46378491116411e+00,
        1.78482653991729e+00,  2.96560571828505e-01,
        2.65321895265761e-02,  1.24266094738808e-03,
        2.71155556874349e-05,  2.01033439929229e-07
    };
    static const double F[8] = {
        1.0,                   5.99832206555888e-01,
        1.36929880922736e-01,  1.48753612908506e-02,
        7.86869131145613e-04,  1.84631831751005e-05,
        1.42151175831645e-07,  2.04426310338994e-15
    };

    if (p <= 0.0) return -8.29;
    if (p >= 1.0) return  8.29;

    const double q = p - 0.5;
    double r, result;

    if (std::abs(q) <= SPLIT1) {
        r = CONST1 - q * q;
        result = q * (((((((A[7]*r+A[6])*r+A[5])*r+A[4])*r+A[3])*r+A[2])*r+A[1])*r+A[0]) /
                     (((((((B[7]*r+B[6])*r+B[5])*r+B[4])*r+B[3])*r+B[2])*r+B[1])*r+B[0]);
    } else {
        r = (q < 0.0) ? p : (1.0 - p);  // min(p, 1-p)
        r = std::sqrt(-std::log(r));
        if (r <= SPLIT2) {
            r -= CONST2;
            result = (((((((C[7]*r+C[6])*r+C[5])*r+C[4])*r+C[3])*r+C[2])*r+C[1])*r+C[0]) /
                     (((((((D[7]*r+D[6])*r+D[5])*r+D[4])*r+D[3])*r+D[2])*r+D[1])*r+D[0]);
        } else {
            r -= SPLIT2;
            result = (((((((E[7]*r+E[6])*r+E[5])*r+E[4])*r+E[3])*r+E[2])*r+E[1])*r+E[0]) /
                     (((((((F[7]*r+F[6])*r+F[5])*r+F[4])*r+F[3])*r+F[2])*r+F[1])*r+F[0]);
        }
        if (q < 0.0) result = -result;
    }
    return result;
}

// ============================================================================
// §2.5–2.6 Owen digit scrambling and HaltonGen
//
// scramble_mode = 0: identity permutations (compat mode, matches randtoolbox exactly)
// scramble_mode = 1: Owen (2017) digit scrambling via Fisher–Yates + splitmix64
//
// HALTON_MAX_DIGITS = 64 is a conservative upper bound on the number of base-b
// digits needed for any practical sequence index (covers indices up to b^64).
// ============================================================================

static const int HALTON_MAX_DIGITS = 64;

struct HaltonGen {
    int S;
    int K_w;
    int scramble_mode;  // 0 = identity (none/compat), 1 = Owen digit scrambling

    // perm[k][d][j]: permuted digit j at digit-position d for dimension k
    // perm[k][d] is a std::vector<uint32_t> of size HALTON_PRIMES[k]
    std::vector<std::vector<std::vector<uint32_t>>> perm;

    // Default constructor: placeholder (non-functional); used by two-step init:
    //   HaltonGen gen;
    //   if (use_generate) gen = HaltonGen(seed, S, K_w, scramble);
    HaltonGen() : S(0), K_w(0), scramble_mode(0) {}

    // Main constructor: build permutation tables from master seed.
    //   seed          — master seed (uint64_t)
    //   S_            — number of draws per individual
    //   K_w_          — number of random-coefficient dimensions (< HALTON_N_PRIMES)
    //   scramble_mode_— 0 = identity, 1 = Owen digit scrambling
    HaltonGen(uint64_t seed, int S_, int K_w_, int scramble_mode_)
        : S(S_), K_w(K_w_), scramble_mode(scramble_mode_),
          perm(K_w_)
    {
        for (int k = 0; k < K_w_; ++k) {
            uint32_t b = HALTON_PRIMES[k];
            perm[k].resize(HALTON_MAX_DIGITS, std::vector<uint32_t>(b));
            for (int d = 0; d < HALTON_MAX_DIGITS; ++d) {
                // Initialize to identity
                for (uint32_t j = 0; j < b; ++j) perm[k][d][j] = j;
                if (scramble_mode_ == 1) {
                    // Fisher-Yates shuffle with a splitmix64 stream derived
                    // from the triple (seed, k, d) via two mix_seed folds.
                    // mix_seed() and splitmix64_next() are from rng.h.
                    // splitmix64_next takes uint64_t& (modifies in place).
                    uint64_t local_seed = mix_seed(
                        mix_seed(seed, static_cast<uint64_t>(k)),
                        static_cast<uint64_t>(d));
                    for (uint32_t j = b - 1; j > 0; --j) {
                        uint64_t rv = splitmix64_next(local_seed);
                        uint32_t swap_idx = static_cast<uint32_t>(rv % (j + 1));
                        uint32_t tmp = perm[k][d][j];
                        perm[k][d][j] = perm[k][d][swap_idx];
                        perm[k][d][swap_idx] = tmp;
                    }
                }
            }
        }
    }

    // Compute one scrambled Halton draw in [0, 1) for dimension k, index n (1-based).
    // In identity mode (scramble_mode = 0), reduces to radical_inverse(n, HALTON_PRIMES[k]).
    inline double scrambled_halton_uniform(uint64_t n, int k) const {
        uint32_t b = HALTON_PRIMES[k];
        double result = 0.0;
        double f = 1.0 / static_cast<double>(b);
        int d = 0;
        while (n > 0) {
            uint32_t digit = static_cast<uint32_t>(n % b);
            uint32_t sd = perm[k][d][digit];  // identity permutation when scramble_mode=0
            result += static_cast<double>(sd) * f;
            n /= b;
            f /= static_cast<double>(b);
            ++d;
        }
        return result;
    }

    // Fill eta_i: write K_w × S standard-normal draws into eta_i for individual i (1-based).
    //
    // Global Halton index: n = (i-1)*S + s + 1   (1-based)
    //   i=1, s=0 → n=1 → phi_2(1) = 0.5, matching randtoolbox start=1 in compat mode.
    //
    // eta_i is resized to K_w × S on entry; the caller owns this buffer (thread-private).
    void fill_eta_i(arma::mat& eta_i, int i) const {
        eta_i.set_size(K_w, S);
        for (int s = 0; s < S; ++s) {
            const uint64_t n = static_cast<uint64_t>(i - 1) *
                               static_cast<uint64_t>(S) +
                               static_cast<uint64_t>(s) + 1ULL;
            for (int k = 0; k < K_w; ++k) {
                double u = scrambled_halton_uniform(n, k);
                eta_i(k, s) = inv_normal_cdf(u);
            }
        }
    }

};

#endif // CHOICER_HALTON_HPP
