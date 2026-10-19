// The damped Newton step behind run_polish(): for the free parameters F,
// solve (D H_FF D + lambda I) y = -D g_F by Cholesky and return the step
// D_F y on F, 0 elsewhere. H is an information matrix (the Hessian of the
// negated log-likelihood, or BHHH), g the gradient of the negated
// log-likelihood, D a diagonal of powers of two (.jacobi_scale()), so that
// the scaled block is formed exactly and lambda I on it is Marquardt's
// lambda |diag H| to within a factor of two. Only the upper triangle of H is
// read (the kernels return it exactly symmetric). A block that is positive
// definite but singular to working precision (a reciprocal condition number,
// estimated in the 1-norm as solve() does, below the machine epsilon) gives
// no step either: its step would be dominated by rounding.
//
// Memory, as in vcov.cpp: every R allocation comes first (the inputs'
// pointers, the result list with its step and info, and LAPACK's small
// arrays from R_alloc()), because an R allocation that fails jumps past C++
// destructors. The n_F x n_F working copy, with the right-hand side after it,
// then comes from new[], whose failure is a C++ exception, and is freed as
// the call returns; nothing after it calls R. LAPACK's xerbla is unreachable
// with these arguments.
#include <RcppArmadillo.h>
#include <cfloat>
#include <cmath>
#include <cstddef>
#include <memory>
#include <new>

namespace {

using arma::blas_int;

// The working copy, count doubles, from new[].
std::unique_ptr<double[]> newton_working_copy(std::size_t count) {
  try {
    return std::unique_ptr<double[]>(new double[count]);
  } catch (const std::bad_alloc&) {
    Rcpp::stop("cannot allocate a working copy of %.1f Gb for the Newton step",
               8.0 * static_cast<double>(count) / 1073741824.0);
  }
}

}  // namespace

//' One damped Newton step on the free parameters (internal)
//'
//' Solves \eqn{(D H_{FF} D + \lambda I) y = -D g_F} by Cholesky (`dpotrf`,
//' `dpotrs`) in one working copy, and returns the step \eqn{D_F y} on the
//' free set \eqn{F}, 0 elsewhere.
//'
//' @param H Square double matrix: the information (Hessian of the negated
//'   log-likelihood, or BHHH). Only its upper triangle is read.
//' @param g Double vector: the gradient of the negated log-likelihood.
//' @param D Double vector of scales, powers of two (`.jacobi_scale(H)`).
//' @param free_idx Integer vector: the free parameters' positions, 1-based
//'   and strictly increasing.
//' @param lambda Damping added to the scaled block's diagonal, finite and
//'   nonnegative.
//' @returns List with `step` (double, one entry per parameter, 0 off the free
//'   set) and `info`: 0 when solved; k > 0 when the scaled block is not
//'   positive definite (its leading minor of order k), or `length(free_idx)
//'   + 1` when it is but its reciprocal condition number is below the
//'   machine epsilon; -1 when the scaled block or the scaled gradient holds a
//'   value that is not finite (an overflow of the scaling included). Without
//'   free parameters the step is 0 and `info` 0.
//' @noRd
// [[Rcpp::export(rng = false)]]
Rcpp::List newton_step(SEXP H, SEXP g, SEXP D, SEXP free_idx, double lambda) {
  if (TYPEOF(H) != REALSXP || !Rf_isMatrix(H)) {
    Rcpp::stop("`H` must be a double matrix");
  }
  const int n = Rf_nrows(H);
  if (Rf_ncols(H) != n) Rcpp::stop("`H` must be square");
  if (TYPEOF(g) != REALSXP || Rf_xlength(g) != n) {
    Rcpp::stop("`g` must be a double vector of length %d", n);
  }
  if (TYPEOF(D) != REALSXP || Rf_xlength(D) != n) {
    Rcpp::stop("`D` must be a double vector of length %d", n);
  }
  if (TYPEOF(free_idx) != INTSXP) {
    Rcpp::stop("`free_idx` must be an integer vector");
  }
  const R_xlen_t nf = Rf_xlength(free_idx);
  const int* fr = INTEGER(free_idx);
  for (R_xlen_t k = 0; k < nf; ++k) {
    if (fr[k] == NA_INTEGER || fr[k] < 1 || fr[k] > n ||
        (k > 0 && fr[k] <= fr[k - 1])) {
      Rcpp::stop("`free_idx` must hold strictly increasing positions in 1..%d",
                 n);
    }
  }
  if (!std::isfinite(lambda) || lambda < 0) {
    Rcpp::stop("`lambda` must be finite and nonnegative");
  }

  // Every R allocation first: the inputs' data (which an ALTREP vector
  // would materialize), the result, LAPACK's small arrays
  const double* h = REAL(H);
  const double* gv = REAL(g);
  const double* d = REAL(D);
  Rcpp::NumericVector step(n);
  Rcpp::IntegerVector info(1);
  Rcpp::List out = Rcpp::List::create(Rcpp::Named("step") = step,
                                      Rcpp::Named("info") = info);
  if (nf == 0) return out;
  const std::size_t m = static_cast<std::size_t>(nf);
  double* lwork = reinterpret_cast<double*>(R_alloc(3 * m, sizeof(double)));
  blas_int* iwork = reinterpret_cast<blas_int*>(R_alloc(m, sizeof(blas_int)));

  std::unique_ptr<double[]> work = newton_working_copy(m * m + m);
  double* a = work.get();
  double* y = a + m * m;
  const std::ptrdiff_t ld = n;

  // The scaled upper triangle and right-hand side. d[i] * d[j], a power of
  // two within [2^-1000, 2^1000], is exact, and so is its product with an
  // entry of H unless that overflows or underflows.
  bool finite = true;
  for (std::size_t jj = 0; jj < m; ++jj) {
    const std::ptrdiff_t j = fr[jj] - 1;
    const double* hj = h + j * ld;
    double* aj = a + jj * m;
    for (std::size_t ii = 0; ii <= jj; ++ii) {
      const std::ptrdiff_t i = fr[ii] - 1;
      double v = hj[i] * (d[i] * d[j]);
      if (ii == jj) v += lambda;
      if (!std::isfinite(v)) finite = false;
      aj[ii] = v;
    }
    y[jj] = -d[j] * gv[j];
    if (!std::isfinite(y[jj])) finite = false;
  }
  if (!finite) {
    info[0] = -1;
    return out;
  }

  char uplo = 'U';
  char norm = '1';
  blas_int bn = static_cast<blas_int>(m);
  blas_int nrhs = 1;
  blas_int status = 0;
  double anorm = arma::lapack::lansy<double>(&norm, &uplo, &bn, a, &bn, lwork);
  arma::lapack::potrf<double>(&uplo, &bn, a, &bn, &status);
  if (status != 0) {
    // status < 0 (an invalid argument) cannot arise with these arguments
    info[0] = status > 0 ? static_cast<int>(status) : -1;
    return out;
  }
  double rcond = 0.0;
  arma::lapack::pocon<double>(&uplo, &bn, a, &bn, &anorm, &rcond, lwork, iwork,
                              &status);
  if (status != 0 || !(rcond >= DBL_EPSILON)) {
    info[0] = static_cast<int>(m) + 1;
    return out;
  }
  arma::lapack::potrs<double>(&uplo, &bn, &nrhs, a, &bn, y, &bn, &status);
  if (status != 0) {
    info[0] = -1;
    return out;
  }
  for (std::size_t jj = 0; jj < m; ++jj) {
    const std::ptrdiff_t j = fr[jj] - 1;
    step[j] = d[j] * y[jj];
  }
  return out;
}
