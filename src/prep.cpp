// Data-preparation helpers behind prepare_mnl_data(), prepare_mxl_data() and
// prepare_mnp_data() (prepare_nl_data() calls the first) and the hierarchical
// preparations, and the column steps behind run_mxlogit()'s default start.
#include <Rcpp.h>
#include <algorithm>
#include <climits>
#include <cmath>
#include <cstdint>
#include <cstring>

namespace {

// bit64's integer64 keeps 64-bit integers in the storage of a double vector,
// NA_integer64_ as INT64_MIN. Element i as bit64's as.double() reads it:
// NA_real_ for NA, otherwise the nearest double (exact below 2^53).
inline double int64_value(const double* in, R_xlen_t i) {
  std::int64_t v;
  std::memcpy(&v, in + i, sizeof v);
  return v == INT64_MIN ? NA_REAL : static_cast<double>(v);
}

// 2^53: from this magnitude up, doubles no longer hold every integer.
constexpr double kExactInt = 9007199254740992.0;

}  // namespace

//' Gather rows of numeric columns into a double design matrix
//'
//' `X[i, k] = cols[[k]][rows[i]]`, written column by column straight from the
//' caller's columns (read only), so the preps never copy the covariates into
//' their working table. Integer columns are converted exactly (NA to
//' NA_real_), as as.matrix() does after coercion to double. integer64
//' (bit64) columns are read as their values, as bit64's as.double() reads
//' them, where as.matrix() would read the raw bits; other classed numerics
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
//'   If integer64 columns hold values of magnitude 2^53 or more in the rows
//'   read, a logical attribute `int64_big`, one element per column, marks
//'   them.
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
  SEXP big = PROTECT(Rf_allocVector(LGLSXP, K));  // see @returns
  int* bg = LOGICAL(big);
  bool any_big = false;
  SEXP X = PROTECT(Rf_allocMatrix(REALSXP, static_cast<int>(n),
                                  static_cast<int>(K)));
  double* out = REAL(X);
  for (R_xlen_t k = 0; k < K; ++k) {
    SEXP c = VECTOR_ELT(cols, k);
    double* o = out + k * n;
    bg[k] = FALSE;
    if (TYPEOF(c) == REALSXP && Rf_inherits(c, "integer64")) {
      const double* in = REAL_RO(c);
      int large = 0;
      if (diff) {
        #pragma omp parallel for schedule(static) reduction(|:large) if (n > 100000)
        for (R_xlen_t i = 0; i < n; ++i) {
          const double v = int64_value(in, r[i] - 1);
          const double w = int64_value(in, b[i] - 1);
          o[i] = (std::isnan(v) || std::isnan(w)) ? NA_REAL : v - w;
          large |= std::fabs(v) >= kExactInt || std::fabs(w) >= kExactInt;
        }
      } else {
        #pragma omp parallel for schedule(static) reduction(|:large) if (n > 100000)
        for (R_xlen_t i = 0; i < n; ++i) {
          o[i] = int64_value(in, r[i] - 1);
          large |= std::fabs(o[i]) >= kExactInt;
        }
      }
      bg[k] = large != 0;
      any_big = any_big || large != 0;
    } else if (TYPEOF(c) == REALSXP) {
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
  if (any_big) Rf_setAttrib(X, Rf_install("int64_big"), big);
  UNPROTECT(2);
  return X;
}

namespace {

// The step of one column of n values (see design_column_step()).
double column_step(const double* x, const R_xlen_t n) {
  if (n < 1) return 0.0;
  // One pass: Boyer and Moore's majority vote (a value that fills more than
  // half of the column is the candidate left standing), and whether the
  // column holds two values only (a = x[0], b the first value apart from it).
  double v = x[0], a = x[0], b = x[0];
  bool second = false, third = false;
  R_xlen_t lead = 0;
  for (R_xlen_t i = 0; i < n; ++i) {
    const double xi = x[i];
    if (lead == 0) {
      v = xi;
      lead = 1;
    } else if (xi == v) {
      ++lead;
    } else {
      --lead;
    }
    if (!(xi == a)) {
      if (!second) {
        b = xi;
        second = true;
      } else if (!(xi == b)) {
        third = true;
      }
    }
  }
  double s, level;
  if (second && !third) {  // two values: the distance between them
    s = std::fabs(a - b);
    level = std::max(std::fabs(a), std::fabs(b));
  } else {
    R_xlen_t count = 0;
    double gap = 0.0, sum = 0.0;
    for (R_xlen_t i = 0; i < n; ++i) {
      sum += x[i];
      if (x[i] == v) {
        ++count;
      } else {
        gap += std::fabs(x[i] - v);
      }
    }
    if (count > n - count) {  // v fills more than half of the column
      if (count == n) return 0.0;
      s = gap / static_cast<double>(n - count);
      level = std::fabs(v);
    } else {
      // The corrected two-pass sum of squares: the sum of the deviations
      // removes the rounding error of the mean.
      const double mean = sum / static_cast<double>(n);
      double d1 = 0.0, d2 = 0.0;
      for (R_xlen_t i = 0; i < n; ++i) {
        const double d = x[i] - mean;
        d1 += d;
        d2 += d * d;
      }
      s = std::sqrt((d2 - d1 * d1 / static_cast<double>(n)) /
                    static_cast<double>(n - 1));
      level = std::fabs(mean);
    }
  }
  // A step below 1e-12 of the level is rounding around a constant.
  return std::isfinite(s) && s > 1e-12 * level ? s : 0.0;
}

}  // namespace

//' The step of each column of a double matrix, read in place
//'
//' The typical step a random coefficient multiplies, behind `run_mxlogit()`'s
//' default start (`.mxl_default_start()`): for a column with two values (a
//' dummy, a two-level attribute), the distance between them; for a column in
//' which one value fills more than half of the rows (a mostly-zero column),
//' the mean distance from that value over the other rows; otherwise the
//' sample standard deviation. A 0/1 dummy's step is 1 whatever its share q of
//' ones, where its standard deviation, sqrt(q (1 - q)), would understate it.
//' Zero for a column whose step is below 1e-12 of its level (its larger value,
//' the dominant value or the mean): a constant column, or rounding around a
//' constant; also for a step that is not finite (a column holding a value that
//' is not finite) and a matrix with no rows. Squared deviations underflow for
//' a column whose values are below about 1e-154 in magnitude, and overflow
//' above about 1e154, which gives 0 too.
//'
//' One pass runs Boyer and Moore's majority vote and counts the distinct
//' values up to three; a second checks the vote and sums; a standard
//' deviation takes a third, the corrected two-pass sum of squares. The matrix
//' is read in place, a column at a time, with double sums in a fixed order:
//' `stats::sd()` on a column would copy it and allocate an index of its rows
//' (1.2 GB more at 9.8e7 rows).
//'
//' @param M A double matrix.
//' @returns A double vector, one step per column.
//' @noRd
// [[Rcpp::export(rng = false)]]
SEXP design_column_step(SEXP M) {
  if (TYPEOF(M) != REALSXP || !Rf_isMatrix(M)) {
    Rcpp::stop("`M` must be a double matrix");
  }
  const R_xlen_t n = Rf_nrows(M);
  const R_xlen_t K = Rf_ncols(M);
  SEXP out = PROTECT(Rf_allocVector(REALSXP, K));
  double* s = REAL(out);
  const double* x = REAL_RO(M);
  for (R_xlen_t k = 0; k < K; ++k) s[k] = column_step(x + k * n, n);
  UNPROTECT(1);
  return out;
}
