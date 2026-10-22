#ifndef CHOICER_DERIVS_HPP
#define CHOICER_DERIVS_HPP

#include "choicer.h"
#include "choicer_internal.h"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <limits>
#include <memory>
#include <new>
#include <string>
#include <vector>

// ============================================================================
// Machinery shared by the derivative kernels (Hessian, BHHH, cluster meat)
// and the memory check of every estimation kernel: the team of threads and
// their padded buffers, the memory guard, the symmetric n x n accumulator
// that units' compact blocks are added into, the cluster score sums, and the
// records that reproduce the dense products' spread of non-finite values.
// The mixed logit's kernels (mxlogit.cpp) use it with likelihood units
// (decision makers, or choice situations) and their slot maps (MxlUnitMap).
// A slot map describes the unit in hand: global(i, K_c) gives the global
// parameter of local index i, strictly increasing in i over the unit's
// block, and local(g, K_c) its inverse, the local index of global parameter
// g or -1 outside the block. The adds below rely on the increase: they pass
// (global(i), global(j)) for local i <= j to SymAcc::add(), which requires
// row <= column.
// ============================================================================

// ============================================================================
// Threads and per-thread buffers (estimation and prediction kernels)
// ============================================================================

// The team of the parallel regions that size per-thread buffers by it, and so
// the number of scratch sets or accumulators they allocate, and the thread
// count every estimation kernel's memory check assumes: OpenMP's next team,
// within the thread limit, and at most one thread per work item (likelihood
// unit or choice situation), but two for a single item, so that the region
// stays active and, as before, an OpenMP-built BLAS runs the item's products
// single-threaded. Kernels without per-thread accumulators (the gradients,
// the score matrices, the mixed logit's conditional tastes) run on OpenMP's
// default team; where it has more threads than items, the extra ones hold
// only their few n-vectors.
inline int kernel_team_threads(const choicer_off n_items) {
#ifdef _OPENMP
  choicer_off n = std::min(omp_get_max_threads(), omp_get_thread_limit());
#else
  choicer_off n = 1;
#endif
  n = std::min(n, std::max<choicer_off>(n_items, 2));
  return static_cast<int>(std::max<choicer_off>(n, 1));
}

inline int kernel_thread_num() {
#ifdef _OPENMP
  return omp_get_thread_num();
#else
  return 0;
#endif
}

// An uninitialized array of n doubles, the thread that uses it touching it
// first, with 128 bytes of padding (a cache line on Apple silicon, two on
// x86-64) so that no two threads' buffers share one.
inline std::unique_ptr<double[]> padded_buffer(const std::size_t n) {
  return std::unique_ptr<double[]>(n > 0 ? new double[n + 16] : nullptr);
}

// ============================================================================
// Memory guard
//
// An estimation kernel allocates scratch for each of its threads and buffers
// they share. Before it does, it compares an upper estimate of the total with
// the memory available: options(choicer.max_memory), in bytes, when set (a
// container's or a job's limit can be below the machine's; Inf lifts the
// check), otherwise the machine's physical memory (no check when it cannot be
// read). Running out inside a parallel region would terminate R; the guard
// stops with an error naming how many threads fit instead. Threads are never
// capped silently. Read on the primary thread.
// ============================================================================

// The limit in bytes (0: none known) and where it comes from.
inline double kernel_memory_limit(const char*& source) {
  static SEXP name = Rf_install("choicer.max_memory");
  SEXP opt = Rf_GetOption1(name);
  if (!Rf_isNull(opt)) {
    double v = NA_REAL;
    if (Rf_length(opt) == 1 && TYPEOF(opt) == REALSXP &&
        Rf_inherits(opt, "integer64")) {
      std::int64_t b;  // bit64's storage; its NA is the smallest int64
      std::memcpy(&b, REAL(opt), sizeof(b));
      if (b != std::numeric_limits<std::int64_t>::min()) v = static_cast<double>(b);
    } else if ((Rf_isReal(opt) || Rf_isInteger(opt)) && Rf_length(opt) == 1) {
      v = Rf_asReal(opt);
    }
    if (!(v > 0)) {
      Rcpp::stop("options(choicer.max_memory =) must be a positive number of "
                 "bytes (Inf to lift the memory check).");
    }
    source = "options(choicer.max_memory)";
    return v;
  }
  static const double physical = choicer_physical_memory();
  source = "physical memory";
  return physical;
}

// "12.34 GB", "567.8 MB" or "89 kB", the unit chosen after rounding.
inline std::string bytes_str(const double bytes) {
  char buf[64];
  if (bytes >= 999.95e6) {
    std::snprintf(buf, sizeof(buf), "%.2f GB", bytes / 1e9);
  } else if (bytes >= 999.5e3) {
    std::snprintf(buf, sizeof(buf), "%.1f MB", bytes / 1e6);
  } else {
    std::snprintf(buf, sizeof(buf), "%.0f kB", bytes / 1e3);
  }
  return buf;
}

// Stop when a call on T threads would take more than the limit: per_thread
// bytes for each thread plus shared(t) bytes with t threads. The total need
// not fall with t (below the accumulators' switch the per-thread triangles
// return), so the message names the largest t that fits, not a range.
template <typename Shared>
inline void kernel_check_memory(const char* what, const int T,
                                const double per_thread, const Shared& shared) {
  const char* source = "";
  const double limit = kernel_memory_limit(source);
  if (!(limit > 0) || std::isinf(limit)) return;
  const auto total = [&](const int t) { return t * per_thread + shared(t); };
  if (total(T) <= limit) return;
  int fit = T - 1;
  while (fit >= 1 && total(fit) > limit) --fit;
  if (fit >= 1) {
    Rcpp::stop("The %s needs about %s with %d threads (%s per thread and %s "
               "shared), more than the %s of %s. Use set_num_threads(%d), or "
               "raise options(choicer.max_memory =).",
               what, bytes_str(total(T)), T, bytes_str(per_thread),
               bytes_str(shared(T)), bytes_str(limit), source, fit);
  }
  Rcpp::stop("The %s needs about %s even with one thread (%s per thread and %s "
             "shared), more than the %s of %s. The estimate is an upper bound: "
             "options(choicer.max_memory = Inf) lifts the check.",
             what, bytes_str(total(1)), bytes_str(per_thread),
             bytes_str(shared(1)), bytes_str(limit), source);
}

// ----------------------------------------------------------------------------
// The derivative kernels' n x n accumulators (SymAcc, below):
// the rows of each column that a thread keeps privately. With T > 1 threads,
// n (its own packed upper triangle) while the T triangles, T n (n + 1) / 2
// doubles, fit in acc_bytes; above it, or for a NaN or negative acc_bytes,
// K_c (the continuous parameters' rows, the threads sharing the result).
// With every parameter continuous the two coincide. One thread adds into the
// result itself (n).
// ----------------------------------------------------------------------------
inline choicer_off sym_acc_rows(const choicer_off n, const choicer_off K_c,
                                const int T, const double acc_bytes) {
  if (T <= 1) return n;
  const double packed = static_cast<double>(n) * static_cast<double>(n + 1) / 2;
  return sizeof(double) * static_cast<double>(T) * packed <= acc_bytes
             ? n
             : std::min(K_c, n);
}

// Offset of column j in a thread's private buffer, which packs each column's
// rows 0..min(s, j + 1) - 1: s = n gives the packed upper triangle, s = K_c
// the continuous rows of a shared result. The buffer holds
// sym_acc_col_off(n, s) doubles, never more than the triangle's.
inline std::size_t sym_acc_col_off(const choicer_off j, const choicer_off s) {
  const std::size_t jj = static_cast<std::size_t>(j);
  const std::size_t ss = static_cast<std::size_t>(s);
  return j <= s ? jj * (jj + 1) / 2 : jj * ss - ss * (ss - 1) / 2;
}

// The n x n result of a derivative kernel (Hessian, BHHH, cluster meat): the
// units' compact blocks are added into its upper triangle, which is mirrored
// at the end (every block is exactly symmetric, so the mirror is the sum the
// lower triangle would have accumulated). With one thread the blocks go into
// the result itself. With more, while the threads' packed upper triangles fit
// in acc_bytes (sym_acc_rows(): 2 GiB by default, n of about 7,000 at eleven
// threads), each thread adds into its own, and the triangles are added up in
// thread order, starting from +0 as the critical section's addition into a
// zero matrix did. They are kept where they fit for speed as much as memory:
// no atomic additions (a clustered meat at census scale took 1.37 s with
// them, 1.82 s shared). Above it the threads share the result: each keeps
// private only the rows of the first K_c parameters (the continuous ones,
// which every unit touches), packed like the triangle's first K_c rows
// (sym_acc_col_off()) and added up in thread order at the end, and adds the
// rest, the constants' block, into the result atomically (with K_c = 0,
// every addition).
// Threads seldom meet there when units face few alternatives (a
// unit touches a_u (a_u + 1) / 2 of its (n - K_c) (n - K_c + 1) / 2 entries);
// a unit with a non-finite score or weight, or a cluster with a wide
// support, adds over many of them. Everything is allocated here, on the
// primary thread, so that running out of memory is an R error saying how
// much each thread needs, not a failure inside a parallel region.
struct SymAcc {
  const choicer_off n;
  const int T;
  const choicer_off K_c;
  const choicer_off s;         // rows a thread keeps per column: n, or K_c
  const bool shared;           // s < n: the threads add into the result
  Rcpp::NumericMatrix result;  // n x n, zero
  double* out = nullptr;       // its elements
  std::vector<std::unique_ptr<double[]>> part;  // T > 1: each thread's rows

  SymAcc(const int n_, const int T_, const int K_c_, const char* what,
         const double acc_bytes)
      : n(n_), T(T_), K_c(K_c_), s(sym_acc_rows(n_, K_c_, T_, acc_bytes)),
        shared(s < n), result(n_, n_) {
    out = result.begin();
    if (T > 1) {
      const std::size_t len = sym_acc_col_off(n, s);
      try {
        part.reserve(T);
        for (int i = 0; i < T; ++i) part.push_back(padded_buffer(len));
      } catch (const std::bad_alloc&) {
        std::vector<std::unique_ptr<double[]>>().swap(part);  // release
        Rcpp::stop("Not enough memory for the %s accumulators: %.2f GB per "
                   "thread for %d threads. Run fewer threads with "
                   "set_num_threads().", what,
                   sizeof(double) * static_cast<double>(len) / 1e9, T);
      }
      // Thread i zeros buffer i (schedule(static)), touching first the
      // memory it will add to.
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(T)
#endif
      for (int i = 0; i < T; ++i) {
        std::fill(part[i].get(), part[i].get() + len, 0.0);
      }
    }
  }
  SymAcc(const SymAcc&) = delete;
  SymAcc& operator=(const SymAcc&) = delete;

  // Thread tid's buffer, to pass to add() only: its layout depends on the
  // mode.
  double* base(const int tid) const { return part.empty() ? out : part[tid].get(); }
  // Add v to entry (i, j), i <= j, from the thread whose buffer is `buf`.
  // Shared must be `shared` (the other instantiation writes outside the
  // buffer): the callers' loops are templates instantiated for each and
  // chosen once per unit (sym_add_block(), bhhh_add(), bhhh_add_dense(),
  // ClusterSum::update_into()), so that the per-thread loops hold no atomic
  // operation (a branch per entry made the per-thread BHHH 12% slower).
  // Atomic = false is for a shared result's serial pass after the parallel
  // region.
  template <bool Shared, bool Atomic = Shared>
  void add(double* buf, const choicer_off i, const choicer_off j, const double v) const {
    const std::size_t ii = static_cast<std::size_t>(i);
    const std::size_t jj = static_cast<std::size_t>(j);
    if (!Shared) {
      buf[T > 1 ? jj * (jj + 1) / 2 + ii : ii + jj * static_cast<std::size_t>(n)] += v;
    } else if (i < s) {
      buf[sym_acc_col_off(j, s) + ii] += v;
    } else {
      double& t = out[ii + jj * static_cast<std::size_t>(n)];
      if (Atomic) {
#ifdef _OPENMP
#pragma omp atomic update
#endif
        t += v;
      } else {
        t += v;
      }
    }
  }

  // The upper triangle of the result: the entries the threads kept, added
  // up in thread order (a no-op with one thread); then the lower triangle as
  // its mirror, every element negated when `negate` (the Hessian's sign).
  void finish(const bool negate) {
    const choicer_off nn = n;
    double* o = out;
    const std::vector<std::unique_ptr<double[]>>& pt = part;
    const int n_part = static_cast<int>(part.size());
    const choicer_off ss = s;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 16) num_threads(T)
#endif
    for (choicer_off j = 0; j < nn; ++j) {
      const std::size_t col = static_cast<std::size_t>(j) * static_cast<std::size_t>(nn);
      if (n_part > 0) {
        // Column j's entries the threads kept (every entry, or with a shared
        // result the continuous rows) from +0, then each thread's in turn:
        // each entry adds the threads in thread order, reading one buffer at
        // a time.
        const choicer_off top = std::min<choicer_off>(ss, j + 1);
        const std::size_t pk = sym_acc_col_off(j, ss);
        double* oc = o + col;
        for (choicer_off i = 0; i < top; ++i) oc[i] = 0.0;
        for (int k = 0; k < n_part; ++k) {
          const double* pc = pt[k].get() + pk;
          for (choicer_off i = 0; i < top; ++i) oc[i] += pc[i];
        }
      }
      for (choicer_off i = 0; i <= j; ++i) {
        const double v = negate ? -o[col + static_cast<std::size_t>(i)]
                                : o[col + static_cast<std::size_t>(i)];
        o[col + static_cast<std::size_t>(i)] = v;
        o[static_cast<std::size_t>(j) + static_cast<std::size_t>(i) * static_cast<std::size_t>(nn)] = v;
      }
    }
  }
};

// Bytes of the derivative kernels' result and of the threads' private
// accumulators with t threads (SymAcc).
inline double sym_acc_need(const choicer_off n, const int K_c, const int t,
                           const double acc_bytes) {
  const double result = sizeof(double) * static_cast<double>(n) * static_cast<double>(n);
  if (t <= 1) return result;
  const choicer_off s = sym_acc_rows(n, K_c, t, acc_bytes);
  return result + sizeof(double) * t * static_cast<double>(sym_acc_col_off(n, s));
}

// The weighted score sum s_g = sum_{u in g} w_u s_u of the cluster in hand,
// as a dense n-vector, zero outside the parameters its units touch (its
// support: the first K_c parameters and the rest of its units' blocks). One
// per thread, and one for each part of a cluster a chunk boundary cuts; every
// vector is allocated on the primary thread, and handing a part over swaps
// vectors, so nothing allocates inside the parallel region. A unit with a non-finite
// weight or score makes the sum dense: its weighted score is added over
// every parameter (w * 0 included) and the update covers every entry, as
// crossprod() of rowsum(w * S) spread NaN and Inf.
struct ClusterSum {
  std::vector<double> s;     // n_params, zero outside the support
  std::vector<char> in;      // support marker
  std::vector<int> supp;     // support, in the order first touched
  bool dense = false;        // support is every parameter
  bool used = false;         // a cut part holds a sum
  int g = 0;                 // its cluster
  char pad[128];             // keeps another thread's sum off these lines

  // Each array carries 128 bytes of slack, so that no two threads' arrays
  // share a cache line.
  explicit ClusterSum(const int n)
      : s(static_cast<std::size_t>(n) + 16, 0.0),
        in(static_cast<std::size_t>(n) + 128, 0) {
    s.resize(static_cast<std::size_t>(n));
    in.resize(static_cast<std::size_t>(n));
    supp.reserve(static_cast<std::size_t>(n) + 32);
  }
  void touch(const int k) {
    if (!in[k]) {
      in[k] = 1;
      supp.push_back(k);
    }
  }
  // Begin cluster g's sum: the continuous parameters are in every support.
  void start(const int g_, const int K_c) {
    g = g_;
    for (int k = 0; k < K_c; ++k) touch(k);
  }
  void make_dense() {
    if (dense) return;
    dense = true;
    for (int k = 0; k < static_cast<int>(s.size()); ++k) touch(k);
  }
  // Add w_u s_u (the unit's compact score), each product formed before the
  // addition, as R's w * S and then rowsum() form it (GCC may contract the
  // pair into one multiply-add where clang does not; the meat is class T).
  template <typename Map>
  void add(const Map& mp, const int K_c, const arma::vec& s_u,
           const double w_u) {
    const int m_u = static_cast<int>(s_u.n_elem);
    if (!std::isfinite(w_u) || !s_u.is_finite()) make_dense();
    if (dense && !std::isfinite(w_u)) {
      // w_u times every entry of the unit's n-vector score, zeros included
      for (int k = 0; k < static_cast<int>(s.size()); ++k) {
        const int l = mp.local(k, K_c);
        const double sk = l >= 0 ? s_u[l] : 0.0;
        const double v = w_u * sk;
        s[k] += v;
      }
      return;
    }
    for (int i = 0; i < m_u; ++i) {
      const int k = mp.global(i, K_c);
      touch(k);
      const double v = w_u * s_u[i];
      s[k] += v;
    }
  }
  // Add another part of the same cluster.
  void merge(const ClusterSum& part) {
    if (part.dense) make_dense();
    for (const int k : part.supp) {
      touch(k);
      s[k] += part.s[k];
    }
  }
  // s_g s_g' over the support into the accumulator's upper triangle. A sum
  // that overflowed (finite weights and scores can add up past the largest
  // double) makes it dense first, as crossprod() spread its Inf (Inf * 0).
  // `serial`: no other thread runs (after the parallel region).
  void update(const SymAcc& acc, double* buf, const bool serial = false) {
    if (!dense && !std::all_of(supp.begin(), supp.end(), [this](const int k) {
          return std::isfinite(s[k]);
        }))
      make_dense();
    std::sort(supp.begin(), supp.end());
    if (!acc.shared) {
      update_into<false, false>(acc, buf);
    } else if (serial) {
      update_into<true, false>(acc, buf);
    } else {
      update_into<true, true>(acc, buf);
    }
  }
  template <bool Shared, bool Atomic>
  void update_into(const SymAcc& acc, double* buf) const {
    const int n = static_cast<int>(supp.size());
    for (int b = 0; b < n; ++b) {
      const int g_b = supp[b];
      const double s_b = s[g_b];
      for (int a = 0; a <= b; ++a) {
        const double v = s[supp[a]] * s_b;
        acc.add<Shared, Atomic>(buf, supp[a], g_b, v);
      }
    }
  }
  void clear() {
    for (const int k : supp) {
      s[k] = 0.0;
      in[k] = 0;
    }
    supp.clear();
    dense = false;
    used = false;
    g = 0;
  }
  // Hand this (cut) sum to an empty part: swap the vectors.
  void hand_to(ClusterSum& part) {
    s.swap(part.s);
    in.swap(part.in);
    supp.swap(part.supp);
    std::swap(dense, part.dense);
    part.g = g;
    part.used = true;
    used = false;
    dense = false;
    g = 0;
  }
};

// Non-finite blocks. In the dense code a unit's block sat in an n x n matrix
// of zeros, and two operations spread a non-finite value beyond the block. A
// non-finite value in a parameter's row of a product gave NaN in that row's
// and that column's entries outside the block, through 0 * NaN and 0 * Inf
// (in the mixed logit, a row of the score stash F_t or of the centered G, in
// the BLAS product as OpenBLAS forms it; a BLAS that skips zero operands
// spread less). And a non-finite unit weight multiplied the zeros outside the
// unit's block: w * 0 is NaN, or NA for an NA weight. The compact kernels
// record both in per-thread arrays allocated on the primary thread, a flag
// per parameter and a count per parameter of the non-finite-weight units
// whose block has it, so that nothing grows inside the parallel region
// however many units are affected; apply_nonfinite() reproduces them on the
// result in one pass: whole rows and columns of NaN for the first, and every
// entry outside the intersection of those units' blocks for the second.
struct NonFiniteRecord {
  std::vector<char> row;  // parameter whose row and column get NaN
  std::vector<int> hits;  // non-finite-weight units whose block has it
  // With row_na (the MNL and NL kernels): the parameter's non-finite value was
  // R's NA, which the dense products spread as NA (empty otherwise)
  std::vector<char> row_na;
  int n_blk = 0;          // non-finite-weight units
  bool na = false;        // one of their weights is NA
  explicit NonFiniteRecord(const int n, const bool with_na = false)
      : row(static_cast<std::size_t>(n), 0), hits(static_cast<std::size_t>(n), 0),
        row_na(with_na ? static_cast<std::size_t>(n) : 0u, 0) {}
};

inline std::vector<NonFiniteRecord> nonfinite_records(const int n, const int n_threads,
                                                      const bool with_na = false) {
  std::vector<NonFiniteRecord> nf;
  nf.reserve(n_threads);
  for (int i = 0; i < n_threads; ++i) nf.emplace_back(n, with_na);
  return nf;
}

// Whether x is R's NA_real_: a NaN whose low-order 32 bits are 1954, the test
// of R_IsNA(), here without the R API so that threads may call it.
inline bool is_na_real(const double x) {
  if (!std::isnan(x)) return false;
  std::uint64_t bits;
  std::memcpy(&bits, &x, sizeof bits);
  return (bits & 0xFFFFFFFFu) == 1954u;
}

// Flag the parameters of the rows of a unit-local or situation-local matrix B
// (rows x S) that hold a non-finite value; `to_global` maps a row of B to its
// global index.
template <typename Map>
inline void nonfinite_rows(const arma::mat& B, Map to_global,
                           NonFiniteRecord& nf) {
  if (B.is_finite()) return;
  for (arma::uword i = 0; i < B.n_rows; ++i) {
    if (!B.row(i).is_finite()) nf.row[to_global(static_cast<int>(i))] = 1;
  }
}

// Flag parameter g when v, the value that spreads over its row and column in
// the dense products, is not finite (and whether it is NA, given row_na).
inline void nonfinite_value(const int g, const double v, NonFiniteRecord& nf) {
  if (std::isfinite(v)) return;
  nf.row[g] = 1;
  if (!nf.row_na.empty() && is_na_real(v)) nf.row_na[g] = 1;
}

// Record a unit with a non-finite weight w over its block of m parameters.
template <typename Map>
inline void nonfinite_weight(const double w, const int m, Map to_global,
                             NonFiniteRecord& nf) {
  for (int i = 0; i < m; ++i) ++nf.hits[to_global(i)];
  ++nf.n_blk;
  if (is_na_real(w)) nf.na = true;
}

// Apply the threads' records to the finished result.
inline void apply_nonfinite(SymAcc& acc,
                            const std::vector<NonFiniteRecord>& nf) {
  const choicer_off n = acc.n;
  std::vector<char> row(static_cast<std::size_t>(n), 0);
  std::vector<char> row_na(static_cast<std::size_t>(n), 0);
  std::vector<int> hits(static_cast<std::size_t>(n), 0);
  int n_blk = 0;
  bool na = false, any_row = false;
  for (const NonFiniteRecord& t : nf) {
    for (choicer_off g = 0; g < n; ++g) {
      row[g] = row[g] | t.row[g];
      if (!t.row_na.empty()) row_na[g] = row_na[g] | t.row_na[g];
      hits[g] += t.hits[g];
      any_row = any_row || t.row[g];
    }
    n_blk += t.n_blk;
    na = na || t.na;
  }
  if (!any_row && n_blk == 0) return;
  const double nan = std::numeric_limits<double>::quiet_NaN();
  const double spread = na ? NA_REAL : nan;
  double* o = acc.out;
  for (choicer_off j = 0; j < n; ++j) {
    const std::size_t col = static_cast<std::size_t>(j) * static_cast<std::size_t>(n);
    for (choicer_off i = 0; i < n; ++i) {
      if (row[i] || row[j]) {
        o[col + static_cast<std::size_t>(i)] = (row_na[i] || row_na[j]) ? NA_REAL : nan;
      } else if (n_blk > 0 && !(hits[i] == n_blk && hits[j] == n_blk)) {
        o[col + static_cast<std::size_t>(i)] = spread;
      }
    }
  }
}

// Add w_u H_u, the unit's m_u x m_u Hessian block (exactly symmetric), into
// the accumulator's upper triangle; each term is the block's entry times the
// weight, formed before it is added, as Armadillo's element-wise update of
// the dense code formed it (see mxl_hessian_parallel()).
template <bool Shared, typename Map>
inline void sym_add_block(const SymAcc& acc, double* buf, const Map& mp,
                          const int Kc, const arma::mat& H, const double w_u) {
  const int m_u = static_cast<int>(H.n_rows);
  for (int j = 0; j < m_u; ++j) {
    const int gj = mp.global(j, Kc);
    for (int i = 0; i <= j; ++i) {
      double v = H(i, j);
      v *= w_u;
      acc.add<Shared>(buf, mp.global(i, Kc), gj, v);
    }
  }
}

// Add w_u s_u s_u', s_u the unit's compact score, into the accumulator's
// upper triangle, column by column; c holds a column (m_u entries). Each term
// is formed as Armadillo's syrk_vec formed the dense rank-one update
// w * s * s.t() (the product of the two scores, times the weight), stored,
// then added, so that the weight's product is never fused into the addition
// (see mxl_bhhh_parallel()).
template <bool Shared, typename Map>
inline void bhhh_add(const SymAcc& acc, double* buf, const Map& mp,
                     const int K_c, const arma::vec& s_u, const double w_u,
                     double* c) {
  const int m_u = static_cast<int>(s_u.n_elem);
  for (int b = 0; b < m_u; ++b) {
    const double s_b = s_u[b];
    for (int a = 0; a <= b; ++a) {
      const double acc1 = s_u[a] * s_b;
      c[a] = w_u * acc1;
    }
    const int g_b = mp.global(b, K_c);
    for (int a = 0; a <= b; ++a) acc.add<Shared>(buf, mp.global(a, K_c), g_b, c[a]);
  }
}

// The same update over every entry, the score read from the slot map (zero
// outside the block), for a unit whose score or weight is not finite; c holds
// a column (n_params entries).
template <bool Shared, typename Map>
inline void bhhh_add_dense(const SymAcc& acc, double* buf, const Map& mp,
                           const int K_c, const arma::vec& s_u,
                           const double w_u, double* c) {
  const auto s_of = [&](const int i) -> double {
    const int l = mp.local(i, K_c);
    return l >= 0 ? s_u[l] : 0.0;
  };
  const int n = static_cast<int>(acc.n);
  for (int i = 0; i < n; ++i) {
    const double s_i = s_of(i);
    for (int k = 0; k <= i; ++k) {
      const double acc1 = s_of(k) * s_i;
      c[k] = w_u * acc1;
    }
    for (int k = 0; k <= i; ++k) acc.add<Shared>(buf, k, i, c[k]);
  }
}

#endif // CHOICER_DERIVS_HPP
