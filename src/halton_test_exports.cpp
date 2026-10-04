// halton_test_exports.cpp — Thin Rcpp wrappers exposing halton.h internals for
// unit testing.  These functions are NOT user-facing API; they are @noRd and
// only exported to make the Phase-A correctness gates callable from R tests.
//
// DO NOT add any of these to the public documentation or NAMESPACE.
// They are kept internal through @noRd roxygen tags.

// [[Rcpp::depends(RcppArmadillo)]]
#include "choicer.h"         // brings in RcppArmadillo.h; arma:: available
#include "halton.h"          // halton.h included AFTER choicer.h (per contract)

//' Radical inverse (van der Corput) for testing halton.h
//'
//' @param n Sequence index (coerced to uint64_t).
//' @param base Prime base (coerced to uint32_t).
//' @return Radical inverse value in [0, 1).
//' @noRd
// [[Rcpp::export]]
double halton_radical_inverse(double n, double base) {
    return radical_inverse(
        static_cast<uint64_t>(n),
        static_cast<uint32_t>(base)
    );
}

//' Wichura AS241 inverse normal CDF for testing halton.h
//'
//' @param p Probability in (0, 1).
//' @return Quantile value.
//' @noRd
// [[Rcpp::export]]
double halton_inv_normal_cdf(double p) {
    return inv_normal_cdf(p);
}

//' Generate an n x dim matrix of uniform digit-permuted Halton draws for testing
//'
//' Returns an n x dim matrix using global indices 1..n (one row per index,
//' one column per dimension). For scramble=0 (compat mode) the result
//' reproduces randtoolbox::halton(n, dim, normal=FALSE) (bit for bit where both
//' are compiled with the same floating-point contraction).
//'
//' @param n Number of Halton points (rows).
//' @param dim Number of dimensions (columns).
//' @param seed Master seed for position-wise digit permutations (coerced to uint64_t).
//'   Ignored when scramble=0.
//' @param scramble 0 = identity (compat), 1 = position-wise digit permutation.
//' @return n x dim arma::mat of uniform [0,1) values.
//' @noRd
// [[Rcpp::export]]
arma::mat halton_generate_uniform(int n, int dim, double seed, int scramble) {
    if (dim < 1 || dim > HALTON_N_PRIMES) {
        Rcpp::stop("dim must be between 1 and %d.", HALTON_N_PRIMES);
    }
    HaltonGen gen(static_cast<uint64_t>(seed), 1, dim, scramble);
    arma::mat out(n, dim);
    for (int row = 0; row < n; ++row) {
        uint64_t idx = static_cast<uint64_t>(row) + 1ULL;  // 1-based global index
        for (int k = 0; k < dim; ++k) {
            out(row, k) = gen.scrambled_halton_uniform(idx, k);
        }
    }
    return out;
}

//' Generate a K_w x (S*N) matrix of normal draws for testing halton.h
//'
//' Built via HaltonGen::fill_eta_i for i=1..N.
//'
//' Layout: columns `[(i-1)*S, i*S)` hold eta_i (K_w x S) for individual i.
//' Within individual i, column s holds the K_w variates for draw s (0-based),
//' so `out(k, (i-1)*S + s) = inv_normal_cdf(phi_{PRIMES[k]}((i-1)*S + s + 1))`.
//'
//' @param S   Number of draws per individual.
//' @param N   Number of individuals.
//' @param K_w Number of random-coefficient dimensions.
//' @param seed Master seed for position-wise digit permutations (coerced to uint64_t).
//' @param scramble 0 = identity (compat), 1 = position-wise digit permutation.
//' @return K_w x (S*N) arma::mat of standard-normal draws.
//' @noRd
// [[Rcpp::export]]
arma::mat halton_generate_normal(int S, int N, int K_w, double seed, int scramble) {
    if (K_w < 1 || K_w > HALTON_N_PRIMES) {
        Rcpp::stop("K_w must be between 1 and %d.", HALTON_N_PRIMES);
    }
    HaltonGen gen(static_cast<uint64_t>(seed), S, K_w, scramble);
    arma::mat out(K_w, static_cast<arma::uword>(S) * static_cast<arma::uword>(N));
    arma::mat eta_i;
    for (int i = 1; i <= N; ++i) {
        gen.fill_eta_i(eta_i, i);
        arma::uword col_start = static_cast<arma::uword>(i - 1) *
                                static_cast<arma::uword>(S);
        out.cols(col_start, col_start + static_cast<arma::uword>(S) - 1) = eta_i;
    }
    return out;
}

// Validate the arguments of the block exports below and return the block's
// first index. n0 and seed are integer-valued doubles: every such double below
// 2^64 converts to uint64_t exactly, including starts past 2^53.
static uint64_t halton_block_start(double n0, int S, int K_w, double seed) {
    const double two64 = 18446744073709551616.0;  // 2^64, exact
    if (K_w < 1 || K_w > HALTON_N_PRIMES) {
        Rcpp::stop("K_w must be between 1 and %d.", HALTON_N_PRIMES);
    }
    if (S < 1) Rcpp::stop("S must be positive.");
    if (!(seed >= 0.0 && seed < two64 && seed == std::floor(seed))) {
        Rcpp::stop("seed must be an integer in [0, 2^64).");
    }
    if (!(n0 >= 1.0 && n0 < two64 && n0 == std::floor(n0))) {
        Rcpp::stop("n0 must be an integer in [1, 2^64).");
    }
    const uint64_t start = static_cast<uint64_t>(n0);
    if (static_cast<uint64_t>(S) - 1 >
        std::numeric_limits<uint64_t>::max() - start) {
        Rcpp::stop("The block n0, ..., n0 + S - 1 passes 2^64 - 1.");
    }
    return start;
}

//' Block of normal draws from HaltonGen::fill_block for testing halton.h
//'
//' @param n0 First global Halton index of the block: an integer-valued double
//'   in [1, 2^64), so blocks past 2^53 are addressable.
//' @param S Number of draws (columns).
//' @param K_w Number of random-coefficient dimensions (rows).
//' @param seed Master seed for position-wise digit permutations: an
//'   integer-valued double in [0, 2^64).
//' @param scramble 0 = identity (compat), 1 = position-wise digit permutation.
//' @return K_w x S arma::mat; column s holds the draw of index n0 + s.
//' @noRd
// [[Rcpp::export]]
arma::mat halton_fill_block(double n0, int S, int K_w, double seed, int scramble) {
    const uint64_t start = halton_block_start(n0, S, K_w, seed);
    HaltonGen gen(static_cast<uint64_t>(seed), S, K_w, scramble);
    arma::mat out(K_w, S);
    gen.fill_block(out.memptr(), start);
    return out;
}

//' Block of uniforms from HaltonGen::fill_uniforms (pass 1 of fill_block)
//'
//' The uniforms the MXL prediction kernels map to normals with R's qnorm()
//' for store-mode draws on the fly (gen_scramble = 2, identity permutations).
//'
//' @param n0 First global Halton index of the block, as for halton_fill_block.
//' @param S Number of draws (columns).
//' @param K_w Number of random-coefficient dimensions (rows).
//' @param seed Master seed, as for halton_fill_block.
//' @param scramble 0 = identity (compat), 1 = position-wise digit permutation.
//' @return K_w x S arma::mat; column s holds the uniforms of index n0 + s.
//' @noRd
// [[Rcpp::export]]
arma::mat halton_fill_uniforms(double n0, int S, int K_w, double seed, int scramble) {
    const uint64_t start = halton_block_start(n0, S, K_w, seed);
    HaltonGen gen(static_cast<uint64_t>(seed), S, K_w, scramble);
    arma::mat out(K_w, S);
    gen.fill_uniforms(out.memptr(), start);
    return out;
}

//' Per-index reference for halton_fill_block
//'
//' The per-index computation that HaltonGen::fill_eta_i ran before the block
//' generator: `out(k, s) = inv_normal_cdf(scrambled_halton_uniform(n0 + s, k))`,
//' draw by draw, with the indices formed in uint64_t.
//'
//' @param n0 First global Halton index of the block, as for halton_fill_block.
//' @param S Number of draws (columns).
//' @param K_w Number of random-coefficient dimensions (rows).
//' @param seed Master seed for position-wise digit permutations, as for
//'   halton_fill_block.
//' @param scramble 0 = identity (compat), 1 = position-wise digit permutation.
//' @return K_w x S arma::mat of standard-normal draws, laid out as in
//'   halton_fill_block.
//' @noRd
// [[Rcpp::export]]
arma::mat halton_reference_block(double n0, int S, int K_w, double seed, int scramble) {
    const uint64_t start = halton_block_start(n0, S, K_w, seed);
    HaltonGen gen(static_cast<uint64_t>(seed), S, K_w, scramble);
    arma::mat out(K_w, S);
    for (int s = 0; s < S; ++s) {
        const uint64_t n = start + static_cast<uint64_t>(s);
        for (int k = 0; k < K_w; ++k) {
            out(k, s) = inv_normal_cdf(gen.scrambled_halton_uniform(n, k));
        }
    }
    return out;
}

//' Layout of HaltonGen's digit-permutation table for testing halton.h
//'
//' The table keeps, per dimension k with base b, the digit positions a 64-bit
//' index can have; the block tests read it on both sides, so its extent is
//' checked here against the base-b digit count of 2^64 - 1.
//'
//' @param K_w Number of random-coefficient dimensions.
//' @return List with `digits` (positions kept per dimension), `offsets`
//'   (start of each dimension's permutations) and `size` (entries in all).
//' @noRd
// [[Rcpp::export]]
Rcpp::List halton_table_layout(int K_w) {
    if (K_w < 1 || K_w > HALTON_N_PRIMES) {
        Rcpp::stop("K_w must be between 1 and %d.", HALTON_N_PRIMES);
    }
    HaltonGen gen(0, 1, K_w, 0);
    Rcpp::IntegerVector digits(K_w);
    Rcpp::NumericVector offsets(K_w);
    for (int k = 0; k < K_w; ++k) {
        digits[k] = halton_index_digits(HALTON_PRIMES[k]);
        offsets[k] = static_cast<double>(gen.perm_off[k]);
    }
    return Rcpp::List::create(
        Rcpp::Named("digits") = digits, Rcpp::Named("offsets") = offsets,
        Rcpp::Named("size") = static_cast<double>(gen.perm.size()));
}
