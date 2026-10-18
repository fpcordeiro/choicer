// The inversions behind a fit's variance (invert_hessian() and
// .sandwich_combine()): an information matrix A, or the sandwich
// A^{-1} B A^{-1}, computed from the equilibrated D A D (and D B D), D a
// diagonal of powers of two (.jacobi_scale()). The equilibration is done in
// the working copy that LAPACK factors in place, so it costs no memory, and
// the products with D are exact.
//
// Memory: every R allocation comes first, LAPACK's small arrays from
// R_alloc() and then the result, uninitialized (an Rcpp object, whose
// protection its destructor releases), because an R allocation that fails
// jumps past C++ destructors. The n x n working copy then comes from new[],
// whose failure is a C++ exception, and is freed as soon as the call
// returns; from R_alloc() it would stay in R's heap until the next garbage
// collection, one p x p more at the peak of a fit. Nothing after it calls
// R: LAPACK's xerbla, which R makes an R error, is unreachable with these
// arguments. A singular matrix leaves the result's allocation untouched.
//
// The passes over whole matrices run over columns, or over columns of
// tiles, in parallel from order 256 (each pragma behind #ifdef _OPENMP, as
// in the kernels; the threshold literal, so that a build without OpenMP
// holds no unused constant); none sums across threads, so these passes do
// not depend on the thread count.
#include <RcppArmadillo.h>
#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstddef>
#include <memory>
#include <new>
#include <utility>

namespace {

using arma::blas_int;

// The order of the square double matrix A, checked against its scale vector
// D. An empty or non-square matrix gets solve()'s message.
blas_int vcov_order(SEXP A, SEXP D) {
  if (TYPEOF(A) != REALSXP || !Rf_isMatrix(A)) {
    Rcpp::stop("`A` must be a double matrix");
  }
  const int n = Rf_nrows(A);
  if (n == 0) Rcpp::stop("'a' is 0-diml");
  if (Rf_ncols(A) != n) {
    Rcpp::stop("'a' (%d x %d) must be square", n, Rf_ncols(A));
  }
  if (TYPEOF(D) != REALSXP || Rf_xlength(D) != n) {
    Rcpp::stop("`D` must be a double vector of length %d", n);
  }
  return n;
}

// LAPACK's small arrays, from R_alloc(): the pivots, and gecon's workspace.
struct LuArrays {
  blas_int* ipiv;
  double* work;
  blas_int* iwork;
  explicit LuArrays(blas_int n)
    : ipiv(reinterpret_cast<blas_int*>(R_alloc(n, sizeof(blas_int)))),
      work(reinterpret_cast<double*>(
        R_alloc(4 * static_cast<std::size_t>(n), sizeof(double)))),
      iwork(reinterpret_cast<blas_int*>(R_alloc(n, sizeof(blas_int)))) {}
};

// The n x n working copy, from new[].
std::unique_ptr<double[]> working_copy(std::size_t n) {
  try {
    return std::unique_ptr<double[]>(new double[n * n]);
  } catch (const std::bad_alloc&) {
    Rcpp::stop("cannot allocate a working copy of %.1f Gb for the inversion",
               8.0 * static_cast<double>(n) * static_cast<double>(n) / 1073741824.0);
  }
}

// out[i, j] = x[i, j] * (d[i] * d[j]) over column-major n x n matrices, the
// product R forms as x * tcrossprod(d): d[i] * d[j] is exact for powers of
// two in [2^-500, 2^500], and the product with x rounds only on underflow or
// overflow. `out` may be `x`.
void scale_by(double* out, const double* x, const double* d, std::ptrdiff_t n) {
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (n >= 256)
#endif
  for (std::ptrdiff_t j = 0; j < n; ++j) {
    const double dj = d[j];
    const double* xj = x + j * n;
    double* oj = out + j * n;
    for (std::ptrdiff_t i = 0; i < n; ++i) oj[i] = xj[i] * (d[i] * dj);
  }
}

// Whether every entry of a column-major n x n matrix is finite.
bool all_finite(const double* a, std::ptrdiff_t n) {
  int bad = 0;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) reduction(|:bad) if (n >= 256)
#endif
  for (std::ptrdiff_t j = 0; j < n; ++j) {
    const double* aj = a + j * n;
    for (std::ptrdiff_t i = 0; i < n; ++i) bad |= !std::isfinite(aj[i]);
  }
  return bad == 0;
}

// The LU factors of the n x n matrix in `a`, in place, with solve()'s
// checks: a pivot that is exactly zero is an error, and so, when every
// factor is finite, is a reciprocal condition number (LAPACK's estimate, in
// the 1-norm) below the machine epsilon. Like solve(), it skips the estimate
// when a factor is not finite (LAPACK 3.11.0 fails on a NaN there).
void lu_checked(double* a, blas_int n, const LuArrays& w) {
  char norm = '1';
  // The norm of the matrix itself, before its factors overwrite it.
  double anorm = arma::lapack::lange<double>(&norm, &n, &n, a, &n, nullptr);
  blas_int info = 0;
  arma::lapack::getrf<double>(&n, &n, a, &n, w.ipiv, &info);
  if (info < 0) {
    Rcpp::stop("argument %d of Lapack routine dgetrf had invalid value", -info);
  }
  if (info > 0) {
    Rcpp::stop("Lapack routine dgetrf: system is exactly singular: U[%d,%d] = 0",
               info, info);
  }
  if (!all_finite(a, n)) return;
  double rcond = 0.0;
  arma::lapack::gecon<double>(&norm, &n, a, &n, &anorm, &rcond, w.work, w.iwork,
                              &info);
  if (rcond < DBL_EPSILON) {
    Rcpp::stop("system is computationally singular: reciprocal condition number = %g",
               rcond);
  }
}

// f(a[i, j], a[j, i]) for every i < j of a column-major n x n matrix, by
// tiles of 64 x 64, so that both elements of every pair stay in cache. A
// column of tiles touches the upper triangle in its own columns and the
// lower triangle in its own rows, so columns of tiles run in parallel, the
// longest (the last, with the most tiles) handed out first.
template <typename F>
void for_each_pair(double* a, std::ptrdiff_t n, F f) {
  constexpr std::ptrdiff_t b = 64;
  const std::ptrdiff_t n_tiles = (n + b - 1) / b;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) if (n >= 256)
#endif
  for (std::ptrdiff_t k = 0; k < n_tiles; ++k) {
    const std::ptrdiff_t t = n_tiles - 1 - k;
    const std::ptrdiff_t j0 = t * b, j1 = std::min(j0 + b, n);
    for (std::ptrdiff_t i0 = 0; i0 <= j0; i0 += b) {
      const std::ptrdiff_t i1 = std::min(i0 + b, n);
      for (std::ptrdiff_t j = j0; j < j1; ++j) {
        const std::ptrdiff_t i_end = i0 == j0 ? j : i1;  // a diagonal tile's upper part
        for (std::ptrdiff_t i = i0; i < i_end; ++i) f(a[i + j * n], a[j + i * n]);
      }
    }
  }
}

}  // namespace

//' The inverse of an information matrix, equilibrated
//'
//' \eqn{V = D (D A D)^{-1} D}: `A` scaled by `D` on both sides into a
//' working copy, LU-factored there in place (`dgetrf`), the identity in the
//' result solved against the factors (`dgetrs`), and the result scaled back.
//' That factorization and solve are what reference LAPACK's `dgesv` performs
//' for `solve(A * tcrossprod(D)) * tcrossprod(D)`, on the same values, so
//' with R's reference LAPACK (or OpenBLAS's) the two agree bit for bit, but
//' without R's temporaries: `solve(A)` holds three n x n matrices besides `A`
//' (the identity, its copy that becomes the result, and the copy LAPACK
//' factors); this holds two, the factors and the result, and frees the
//' factors as it returns. The checks are `solve()`'s: an error for a pivot
//' that is exactly zero, and, when every factor is finite, for a reciprocal
//' condition number below the machine epsilon. (LAPACK's
//' in-place inverse, `dgetri`, would need one, but its many thin BLAS calls
//' ran several times slower than the two solves with threads.)
//'
//' With `D` from `.jacobi_scale()` (powers of two), every product with `D` is
//' exact unless it underflows or overflows, so only the factorization and the
//' solves round, on a matrix whose diagonal is within a factor of two of one
//' in magnitude.
//'
//' @param A The information matrix, a square double matrix (read only).
//' @param D The scales, a double vector of length `nrow(A)`.
//' @returns The double matrix \eqn{V}, without dimnames.
//' @noRd
// [[Rcpp::export(rng = false)]]
SEXP equilibrated_inverse(SEXP A, SEXP D) {
  blas_int n = vcov_order(A, D);
  const std::ptrdiff_t ns = n;
  const double* a = REAL(A);
  const double* d = REAL(D);
  const LuArrays w(n);
  Rcpp::NumericMatrix V = Rcpp::no_init(n, n);
  double* v = V.begin();
  const std::unique_ptr<double[]> lu = working_copy(ns);
  scale_by(lu.get(), a, d, ns);
  lu_checked(lu.get(), n, w);

#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (ns >= 256)
#endif
  for (std::ptrdiff_t j = 0; j < ns; ++j) {
    double* vj = v + j * ns;
    std::fill(vj, vj + ns, 0.0);
    vj[j] = 1.0;
  }
  char no_trans = 'N';
  blas_int info = 0;
  arma::lapack::getrs<double>(&no_trans, &n, &n, lu.get(), &n, w.ipiv, v, &n,
                              &info);
  if (info != 0) {
    Rcpp::stop("Lapack routine dgetrs gave error code %d", info);
  }
  scale_by(v, v, d, ns);
  return V;
}

//' The sandwich of an information matrix and a meat, equilibrated
//'
//' \eqn{V = D \, \mathrm{sym}[(DAD)^{-1} (DBD) (DAD)^{-1}] \, D}, with
//' \eqn{\mathrm{sym}(M) = (M + M^\top) / 2}, without an explicit inverse:
//' `A` is scaled into a working copy and LU-factored there (`dgetrf`, with
//' `equilibrated_inverse()`'s checks); `B` is scaled into the result, which
//' becomes \eqn{X = (DAD)^{-1} (DBD)} by one solve (`dgetrs`), is transposed
//' in place, and becomes \eqn{(DAD)^{-\top} X^\top = [X (DAD)^{-1}]^\top} by
//' a second (the transposed system). Symmetrized and scaled back, that is the
//' sandwich. Two n x n matrices besides `A` and `B`, the working copy and the
//' result. Forming \eqn{(DAD)^{-1}} as `solve()` does and multiplying twice
//' takes 43% more arithmetic, and three such matrices in place (four in R).
//'
//' @param A The bread (the information), a square double matrix.
//' @param B The meat, a double matrix of the same order.
//' @param D The scales, a double vector of length `nrow(A)` (from `A`).
//' @returns The symmetric double matrix \eqn{V}, without dimnames.
//' @noRd
// [[Rcpp::export(rng = false)]]
SEXP equilibrated_sandwich(SEXP A, SEXP B, SEXP D) {
  blas_int n = vcov_order(A, D);
  if (TYPEOF(B) != REALSXP || !Rf_isMatrix(B) || Rf_nrows(B) != n ||
      Rf_ncols(B) != n) {
    Rcpp::stop("`B` must be a double matrix of the order of `A`");
  }
  const std::ptrdiff_t ns = n;
  const double* a = REAL(A);
  const double* bm = REAL(B);
  const double* d = REAL(D);
  const LuArrays w(n);
  Rcpp::NumericMatrix V = Rcpp::no_init(n, n);
  double* v = V.begin();
  const std::unique_ptr<double[]> lu = working_copy(ns);
  scale_by(lu.get(), a, d, ns);
  lu_checked(lu.get(), n, w);

  scale_by(v, bm, d, ns);
  char no_trans = 'N', trans = 'T';
  blas_int info = 0;
  arma::lapack::getrs<double>(&no_trans, &n, &n, lu.get(), &n, w.ipiv, v, &n,
                              &info);
  if (info == 0) {
    for_each_pair(v, ns, [](double& x, double& y) { std::swap(x, y); });
    arma::lapack::getrs<double>(&trans, &n, &n, lu.get(), &n, w.ipiv, v, &n,
                                &info);
  }
  if (info != 0) {
    Rcpp::stop("Lapack routine dgetrs gave error code %d", info);
  }
  // sym(M) from its transpose: (M[i, j] + M[j, i]) / 2 is the same sum; the
  // diagonal is its own.
  for_each_pair(v, ns, [](double& x, double& y) { x = y = (x + y) / 2; });
  scale_by(v, v, d, ns);
  return V;
}
