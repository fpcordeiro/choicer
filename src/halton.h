// halton.h — Header-only on-the-fly Halton draw generator
//
// Provides:
//   HALTON_PRIMES[]       — 128-prime table (primes[k] is the base for dimension k)
//   radical_inverse()     — van der Corput radical-inverse function
//   inv_normal_cdf()      — Wichura (1988) AS241 PPND16 inverse normal CDF
//   struct HaltonGen      — position-wise digit-permuted Halton generator
//
// DESIGN CONTRACT:
//   - Header-only: all functions are inline or templated; no .cpp counterpart.
//   - NO R API: does not include <Rcpp.h>, <R.h>, <Rmath.h>, or call Rf_* functions.
//     Safe to use inside OpenMP parallel regions.
//   - Include order: this header must be included from a .cpp translation unit that
//     has already included choicer.h (which brings in RcppArmadillo.h), so arma::
//     types are available here even though halton.h does not include them directly.
//   - Thread safety: HaltonGen is const after construction; fill_block() and
//     fill_eta_i() write only a buffer owned by the calling thread and keep their
//     working state on the stack. Multiple threads may call them simultaneously on
//     the same const HaltonGen without data races.
//   - Bitwise reproducibility: n = (i-1)*S + s + 1 is a deterministic function of
//     (i, s, S) and the permutation table is a deterministic function of the seed,
//     so results are identical regardless of OpenMP thread count or schedule.

#ifndef CHOICER_HALTON_HPP
#define CHOICER_HALTON_HPP

#include <cstdint>
#include <cmath>
#include <limits>
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
// §2.5–2.6 Position-wise digit permutations and HaltonGen
//
// scramble_mode = 0: identity permutations (compat mode, matches randtoolbox exactly)
// scramble_mode = 1: seeded Fisher–Yates digit permutation for every
// (dimension k, digit position d), shared across all sequence indices.
//
// This is NOT Owen's nested-uniform scramble: Owen's permutation at position d
// depends on the preceding d-digit prefix, whereas the permutation of
// (k, d) does not. Also, scrambled_halton_uniform() stops when n == 0, so
// implicit trailing zero digits are left unchanged instead of being permuted.
// Consequently this construction must not be advertised with standard RQMC
// marginal-uniformity, unbiasedness, variance-rate, or replicate-error
// guarantees. The public R value "owen" is retained only as a deprecated
// compatibility alias for mode 1.
//
// HALTON_MAX_DIGITS = 64 is the number of base-2 digits of a uint64_t index,
// hence an upper bound on the digits of any index in any base; the stride of
// the place-value table. The permutation table keeps, per base b, only the
// positions an index can reach: halton_index_digits(b).
// ============================================================================

static const int HALTON_MAX_DIGITS = 64;

// Largest index of a block whose base-2 uniforms come from the integer
// odometer (HaltonGen::uniforms_odometer2): 2^53 - 1, so that an index has at
// most 53 binary digits.
static const uint64_t HALTON_ODOMETER_MAX = (static_cast<uint64_t>(1) << 53) - 1;

// Number of base-b digits of the largest index, 2^64 - 1: the digit positions
// any uint64_t index can have (64 in base 2, 41 in base 3, 7 in base 719).
inline int halton_index_digits(uint32_t b) {
    int digits = 0;
    for (uint64_t m = std::numeric_limits<uint64_t>::max(); m > 0; m /= b) ++digits;
    return digits;
}

// ============================================================================
// Block generation: fill_block() and why its draws are bit-identical to the
// per-index loop
//
//   for s = 0..S-1, k = 0..K_w-1:
//     e[s * K_w + k] = inv_normal_cdf(scrambled_halton_uniform(n0 + s, k))
//
// that fill_eta_i() ran before (tests: halton_fill_block() against
// halton_reference_block()).
//   1. Two passes. Pass 1 writes every uniform of the block, dimension by
//      dimension; pass 2 maps each through inv_normal_cdf(). Each value goes
//      through the same operations as in the per-index loop; only the order in
//      which different values are computed changes, and no value depends on
//      another.
//   2. Digits. uniforms_from_digits<B>() takes digit d as n - (n / b) * b with
//      the base a compile-time constant for b <= 31: exact integer arithmetic,
//      so the digits are those of n % b.
//   3. Terms and order. Each permuted digit is the same table entry, and the
//      place value place[k * HALTON_MAX_DIGITS + d] is the f of
//      scrambled_halton_uniform() after d divisions, computed by the same
//      divisions. The terms are added from digit 0 upward in the same
//      expression form, u += (double)digit * f, so a compiler that contracts it
//      into a fused multiply-add does so in both. This assumes doubles are
//      evaluated in double precision (FLT_EVAL_METHOD == 0, true on every
//      platform R supports but 32-bit x87 builds, where f may be held in an
//      80-bit register while place[d] is rounded).
//   4. Base 2 by an integer odometer (uniforms_odometer2(), blocks whose last
//      index is below 2^53). A base-2 term is 0 or 2^-(d+1). An index below
//      2^53 has at most 53 binary digits, so every partial sum of the
//      floating-point loop is exact, and its result is U 2^-63 with the
//      integer U = sum_d P_d(digit_d) 2^(62-d). U is a multiple of 2^10 below
//      2^63, so it has at most 53 significant bits and converts to double
//      exactly, and scaling by 2^-63 is exact. The odometer maintains U from
//      one index to the next: see uniforms_odometer2(). Blocks that reach 2^53
//      use the digit loop.
// ============================================================================

struct HaltonGen {
    int S;
    int K_w;
    int scramble_mode;  // 0 = identity; 1 = position-wise digit permutation

    // Digit permutations, flat: at digit position d of dimension k (base
    // b = HALTON_PRIMES[k]), digit j becomes perm[perm_off[k] + d * b + j],
    // for d < halton_index_digits(b). At K_w = 128 this is 1.3 MB.
    std::vector<uint32_t> perm;
    std::vector<size_t> perm_off;
    // Place values: place[k * HALTON_MAX_DIGITS + d] = b^-(d+1), each the f
    // of scrambled_halton_uniform() after d divisions.
    std::vector<double> place;

    // Default constructor: placeholder (non-functional); used by two-step init:
    //   HaltonGen gen;
    //   if (use_generate) gen = HaltonGen(seed, S, K_w, scramble);
    HaltonGen() : S(0), K_w(0), scramble_mode(0) {}

    // Main constructor: build permutation tables from master seed.
    //   seed          — master seed (uint64_t)
    //   S_            — number of draws per individual
    //   K_w_          — number of random-coefficient dimensions (<= HALTON_N_PRIMES)
    //   scramble_mode_— 0 = identity, 1 = position-wise digit permutation
    HaltonGen(uint64_t seed, int S_, int K_w_, int scramble_mode_)
        : S(S_), K_w(K_w_), scramble_mode(scramble_mode_),
          perm_off(K_w_), place(static_cast<size_t>(K_w_) * HALTON_MAX_DIGITS)
    {
        size_t n_perm = 0;
        for (int k = 0; k < K_w_; ++k) {
            perm_off[k] = n_perm;
            n_perm += static_cast<size_t>(halton_index_digits(HALTON_PRIMES[k])) *
                      HALTON_PRIMES[k];
        }
        perm.resize(n_perm);
        for (int k = 0; k < K_w_; ++k) {
            const uint32_t b = HALTON_PRIMES[k];
            const int n_digits = halton_index_digits(b);
            for (int d = 0; d < n_digits; ++d) {
                uint32_t* p = &perm[perm_off[k] + static_cast<size_t>(d) * b];
                // Initialize to identity
                for (uint32_t j = 0; j < b; ++j) p[j] = j;
                if (scramble_mode_ == 1) {
                    // Fisher-Yates shuffle with a splitmix64 stream derived
                    // from the triple (seed, k, d) via two mix_seed folds.
                    // mix_seed() and splitmix64_next() are from rng.h.
                    // splitmix64_next takes uint64_t& (modifies in place).
                    // Each position has its own stream, so the positions no
                    // index reaches can be left out without changing the
                    // others.
                    uint64_t local_seed = mix_seed(
                        mix_seed(seed, static_cast<uint64_t>(k)),
                        static_cast<uint64_t>(d));
                    for (uint32_t j = b - 1; j > 0; --j) {
                        uint64_t rv = splitmix64_next(local_seed);
                        uint32_t swap_idx = static_cast<uint32_t>(rv % (j + 1));
                        uint32_t tmp = p[j];
                        p[j] = p[swap_idx];
                        p[swap_idx] = tmp;
                    }
                }
            }
            double f = 1.0 / static_cast<double>(b);
            for (int d = 0; d < HALTON_MAX_DIGITS; ++d) {
                place[static_cast<size_t>(k) * HALTON_MAX_DIGITS + d] = f;
                f /= static_cast<double>(b);
            }
        }
    }

    // Compute one digit-permuted Halton draw in [0, 1) for dimension k, index n.
    // Only the explicit base-b digits of n are processed; trailing zeros are not.
    // In identity mode (scramble_mode = 0), reduces to radical_inverse(n, HALTON_PRIMES[k]).
    // The per-index reference for fill_block().
    inline double scrambled_halton_uniform(uint64_t n, int k) const {
        uint32_t b = HALTON_PRIMES[k];
        const uint32_t* P = perm.data() + perm_off[k];
        double result = 0.0;
        double f = 1.0 / static_cast<double>(b);
        int d = 0;
        while (n > 0) {
            uint32_t digit = static_cast<uint32_t>(n % b);
            uint32_t sd = P[static_cast<size_t>(d) * b + digit];  // identity when scramble_mode=0
            result += static_cast<double>(sd) * f;
            n /= b;
            f /= static_cast<double>(b);
            ++d;
        }
        return result;
    }

    // Fill e[s * K_w + k] (column-major K_w x S) with the standard-normal draws
    // of indices n0, ..., n0 + S - 1: all uniforms first, then the inverse
    // normal CDF (see the proof above). e must hold K_w * S doubles, and the
    // last index n0 + S - 1 must not pass 2^64 - 1.
    void fill_block(double* e, const uint64_t n0) const {
        if (S <= 0 || K_w <= 0) return;
        const uint64_t span = static_cast<uint64_t>(S) - 1;  // last index n0 + span
        const bool odometer = n0 <= HALTON_ODOMETER_MAX &&
                              span <= HALTON_ODOMETER_MAX - n0;
        for (int k = 0; k < K_w; ++k) {
            switch (HALTON_PRIMES[k]) {
            case 2:
                if (odometer) uniforms_odometer2(e, k, n0);
                else uniforms_from_digits<2>(e, k, n0);
                break;
            case 3:  uniforms_from_digits<3>(e, k, n0);  break;
            case 5:  uniforms_from_digits<5>(e, k, n0);  break;
            case 7:  uniforms_from_digits<7>(e, k, n0);  break;
            case 11: uniforms_from_digits<11>(e, k, n0); break;
            case 13: uniforms_from_digits<13>(e, k, n0); break;
            case 17: uniforms_from_digits<17>(e, k, n0); break;
            case 19: uniforms_from_digits<19>(e, k, n0); break;
            case 23: uniforms_from_digits<23>(e, k, n0); break;
            case 29: uniforms_from_digits<29>(e, k, n0); break;
            case 31: uniforms_from_digits<31>(e, k, n0); break;
            default: uniforms_from_digits<0>(e, k, n0);
            }
        }
        const size_t n_eta = static_cast<size_t>(K_w) * static_cast<size_t>(S);
        for (size_t j = 0; j < n_eta; ++j) e[j] = inv_normal_cdf(e[j]);
    }

    // Fill eta_i: write K_w × S standard-normal draws into eta_i for individual i (1-based).
    //
    // Global Halton index: n = (i-1)*S + s + 1   (1-based)
    //   i=1, s=0 → n=1 → phi_2(1) = 0.5, matching randtoolbox start=1 in compat mode.
    //
    // eta_i is resized to K_w × S on entry; the caller owns this buffer (thread-private).
    // i is 64-bit, so the start of a block past the 2^31 - 1st does not overflow.
    void fill_eta_i(arma::mat& eta_i, const uint64_t i) const {
        eta_i.set_size(K_w, S);
        fill_block(eta_i.memptr(), (i - 1) * static_cast<uint64_t>(S) + 1ULL);
    }

private:
    // Pass-1 helpers of fill_block(), which establishes their preconditions.

    // Uniforms of dimension k for indices n0, ..., n0 + S - 1, written to
    // e[s * K_w + k], each from its own digits as in scrambled_halton_uniform().
    // B > 0 fixes the base at compile time, so the digit division is by a
    // constant; B = 0 reads it from the table.
    template <uint32_t B>
    void uniforms_from_digits(double* e, const int k, const uint64_t n0) const {
        const uint32_t b = B ? B : HALTON_PRIMES[k];
        const uint32_t* P = perm.data() + perm_off[k];
        const double* F = place.data() + static_cast<size_t>(k) * HALTON_MAX_DIGITS;
        for (int s = 0; s < S; ++s) {
            uint64_t n = n0 + static_cast<uint64_t>(s);
            double u = 0.0;
            int d = 0;
            while (n > 0) {
                const uint64_t q = n / b;
                const uint32_t digit = static_cast<uint32_t>(n - q * b);
                u += static_cast<double>(P[static_cast<size_t>(d) * b + digit]) * F[d];
                n = q;
                ++d;
            }
            e[static_cast<size_t>(s) * K_w + k] = u;
        }
    }

    // Uniforms of base-2 dimension k for indices n0, ..., n0 + S - 1 by an
    // exact integer odometer; requires n0 + S - 1 <= HALTON_ODOMETER_MAX.
    // U = sum_d P_d(digit_d) 2^(62-d) over the explicit binary digits of the
    // current index, and u = U 2^-63 (see the proof above). The next index is
    // n + 1: its explicit 1-digits from the bottom turn into explicit 0-digits,
    // which are still permuted (P_d(0), not dropped), and the first 0-digit
    // becomes 1, or, past the top digit, a new top digit 1 appears. U changes
    // by the permuted values of the changed digits only.
    void uniforms_odometer2(double* e, const int k, const uint64_t n0) const {
        const uint32_t* P = perm.data() + perm_off[k];  // P[2 * d + digit]
        const double scale = 1.0 / 9223372036854775808.0;  // 2^-63, exact
        uint32_t digit[HALTON_MAX_DIGITS];
        int n_digits = 0;  // explicit digits of the current index
        uint64_t U = 0;
        for (uint64_t m = n0; m > 0; m >>= 1, ++n_digits) {
            digit[n_digits] = static_cast<uint32_t>(m & 1u);
            U += static_cast<uint64_t>(P[2 * n_digits + digit[n_digits]])
                 << (62 - n_digits);
        }
        e[k] = static_cast<double>(U) * scale;
        for (int s = 1; s < S; ++s) {
            int d = 0;
            while (d < n_digits && digit[d] == 1u) {  // 1 -> explicit 0, carry
                U = U - (static_cast<uint64_t>(P[2 * d + 1]) << (62 - d)) +
                    (static_cast<uint64_t>(P[2 * d]) << (62 - d));
                digit[d] = 0u;
                ++d;
            }
            if (d == n_digits) {  // new top digit 1
                ++n_digits;
                U += static_cast<uint64_t>(P[2 * d + 1]) << (62 - d);
            } else {              // explicit 0 -> 1
                U = U - (static_cast<uint64_t>(P[2 * d]) << (62 - d)) +
                    (static_cast<uint64_t>(P[2 * d + 1]) << (62 - d));
            }
            digit[d] = 1u;
            e[static_cast<size_t>(s) * K_w + k] = static_cast<double>(U) * scale;
        }
    }
};

#endif // CHOICER_HALTON_HPP
