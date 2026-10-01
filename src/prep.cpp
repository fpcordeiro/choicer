// Data-preparation helpers behind prepare_mnl_data(), prepare_mxl_data() and
// prepare_mnp_data() (prepare_nl_data() calls the first) and the hierarchical
// preparations.
#include <Rcpp.h>
#include <climits>

//' Gather rows of numeric columns into a double design matrix
//'
//' `X[i, k] = cols[[k]][rows[i]]`, written column by column straight from the
//' caller's columns (read only), so the preps never copy the covariates into
//' their working table. Integer columns are converted exactly (NA to
//' NA_real_), as as.matrix() does after coercion to double; classed numerics
//' are read from their raw storage, as as.matrix()'s unlist() reads them.
//'
//' With `base`, the same pass writes the differenced design of the
//' multinomial probit, `X[i, k] = cols[[k]][rows[i]] - cols[[k]][base[i]]`:
//' the difference of the two values in double (NA when either is NA), as R's
//' arithmetic forms it, without a second design-sized matrix.
//'
//' @param cols List of integer or double vectors.
//' @param rows Integer vector of 1-based indices into them.
//' @param base `NULL`, or an integer vector of 1-based indices, one per
//'   element of `rows`, of the values to subtract.
//' @returns A `length(rows)` x `length(cols)` double matrix, no dimnames.
//' @noRd
// [[Rcpp::export(rng = false)]]
SEXP prep_gather_design(SEXP cols, SEXP rows, SEXP base = R_NilValue) {
  if (TYPEOF(cols) != VECSXP) Rcpp::stop("`cols` must be a list");
  if (TYPEOF(rows) != INTSXP) Rcpp::stop("`rows` must be an integer vector");
  const bool diff = !Rf_isNull(base);
  if (diff && TYPEOF(base) != INTSXP) {
    Rcpp::stop("`base` must be NULL or an integer vector");
  }
  const R_xlen_t n = XLENGTH(rows);
  const R_xlen_t K = XLENGTH(cols);
  if (n > INT_MAX || K > INT_MAX) Rcpp::stop("too many rows or columns");
  if (diff && XLENGTH(base) != n) {
    Rcpp::stop("`base` must have one index per element of `rows`");
  }
  const int* r = INTEGER_RO(rows);
  const int* b = diff ? INTEGER_RO(base) : nullptr;
  // An index below 1 (NA_INTEGER included) or past a column's end is an
  // error; r_min only needs to fall below 1 to catch the first kind.
  int r_min = 1, r_max = 0;
  for (R_xlen_t i = 0; i < n; ++i) {
    if (r[i] < r_min) r_min = r[i];
    if (r[i] > r_max) r_max = r[i];
    if (diff) {
      if (b[i] < r_min) r_min = b[i];
      if (b[i] > r_max) r_max = b[i];
    }
  }
  for (R_xlen_t k = 0; k < K; ++k) {
    SEXP c = VECTOR_ELT(cols, k);
    if (TYPEOF(c) != REALSXP && TYPEOF(c) != INTSXP) {
      Rcpp::stop("column %d is neither integer nor double",
                 static_cast<int>(k + 1));
    }
    if (n > 0 && (r_min < 1 || static_cast<R_xlen_t>(r_max) > XLENGTH(c))) {
      Rcpp::stop("row indices out of range for column %d",
                 static_cast<int>(k + 1));
    }
  }
  // All checks are done: nothing below throws.
  SEXP X = PROTECT(Rf_allocMatrix(REALSXP, static_cast<int>(n),
                                  static_cast<int>(K)));
  double* out = REAL(X);
  for (R_xlen_t k = 0; k < K; ++k) {
    SEXP c = VECTOR_ELT(cols, k);
    double* o = out + k * n;
    if (TYPEOF(c) == REALSXP) {
      const double* in = REAL_RO(c);
      if (diff) {
        #pragma omp parallel for schedule(static) if (n > 100000)
        for (R_xlen_t i = 0; i < n; ++i) o[i] = in[r[i] - 1] - in[b[i] - 1];
      } else {
        #pragma omp parallel for schedule(static) if (n > 100000)
        for (R_xlen_t i = 0; i < n; ++i) o[i] = in[r[i] - 1];
      }
    } else {
      const int* in = INTEGER_RO(c);
      if (diff) {
        #pragma omp parallel for schedule(static) if (n > 100000)
        for (R_xlen_t i = 0; i < n; ++i) {
          const int v = in[r[i] - 1], w = in[b[i] - 1];
          o[i] = (v == NA_INTEGER || w == NA_INTEGER)
            ? NA_REAL
            : static_cast<double>(v) - static_cast<double>(w);
        }
      } else {
        #pragma omp parallel for schedule(static) if (n > 100000)
        for (R_xlen_t i = 0; i < n; ++i) {
          const int v = in[r[i] - 1];
          o[i] = (v == NA_INTEGER) ? NA_REAL : static_cast<double>(v);
        }
      }
    }
  }
  UNPROTECT(1);
  return X;
}
