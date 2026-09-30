// Data-preparation helpers behind prepare_mnl_data() and prepare_mxl_data()
// (and prepare_nl_data(), which calls the former).
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
//' @param cols List of integer or double vectors.
//' @param rows Integer vector of 1-based indices into them.
//' @returns A `length(rows)` x `length(cols)` double matrix, no dimnames.
//' @noRd
// [[Rcpp::export(rng = false)]]
SEXP prep_gather_design(SEXP cols, SEXP rows) {
  if (TYPEOF(cols) != VECSXP) Rcpp::stop("`cols` must be a list");
  if (TYPEOF(rows) != INTSXP) Rcpp::stop("`rows` must be an integer vector");
  const R_xlen_t n = XLENGTH(rows);
  const R_xlen_t K = XLENGTH(cols);
  if (n > INT_MAX || K > INT_MAX) Rcpp::stop("too many rows or columns");
  const int* r = INTEGER_RO(rows);
  int r_min = 1, r_max = 0;
  for (R_xlen_t i = 0; i < n; ++i) {
    if (i == 0 || r[i] < r_min) r_min = r[i];
    if (r[i] > r_max) r_max = r[i];
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
      #pragma omp parallel for schedule(static) if (n > 100000)
      for (R_xlen_t i = 0; i < n; ++i) o[i] = in[r[i] - 1];
    } else {
      const int* in = INTEGER_RO(c);
      #pragma omp parallel for schedule(static) if (n > 100000)
      for (R_xlen_t i = 0; i < n; ++i) {
        const int v = in[r[i] - 1];
        o[i] = (v == NA_INTEGER) ? NA_REAL : static_cast<double>(v);
      }
    }
  }
  UNPROTECT(1);
  return X;
}
