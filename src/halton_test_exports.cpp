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

//' Generate an n x dim matrix of uniform scrambled-Halton draws for testing
//'
//' Returns an n x dim matrix using global indices 1..n (one row per index,
//' one column per dimension). For scramble=0 (compat mode) the result
//' reproduces randtoolbox::halton(n, dim, normal=FALSE) exactly.
//'
//' @param n Number of Halton points (rows).
//' @param dim Number of dimensions (columns).
//' @param seed Master seed for Owen scrambling (coerced to uint64_t). Ignored when scramble=0.
//' @param scramble 0 = identity (compat), 1 = Owen digit scrambling.
//' @return n x dim arma::mat of uniform [0,1) values.
//' @noRd
// [[Rcpp::export]]
arma::mat halton_generate_uniform(int n, int dim, double seed, int scramble) {
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
//' @param seed Master seed for Owen scrambling (coerced to uint64_t).
//' @param scramble 0 = identity (compat), 1 = Owen digit scrambling.
//' @return K_w x (S*N) arma::mat of standard-normal draws.
//' @noRd
// [[Rcpp::export]]
arma::mat halton_generate_normal(int S, int N, int K_w, double seed, int scramble) {
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
