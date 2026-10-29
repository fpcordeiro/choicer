#ifndef CHOICER_FIXED_POINT_HPP
#define CHOICER_FIXED_POINT_HPP

#include "choicer.h"
#include <algorithm>
#include <cmath>
#include <limits>

// ============================================================================
// Fixed-point iteration x <- x + damping * step(x), on a map that reports, at
// every point it evaluates, a stopping residual and a base step: the share
// inversion's (blp_engine.h), and any later map of the same form (a markup
// equation, an EM step). A map provides
//   arma::uword n() const;  // its coordinates
//   bool eval(const arma::vec& x, double& residual, arma::vec& step);
//   bool has_fallback() const;  // a second step, from the same evaluation
//   void use_fallback(bool on);  // eval() writes that step while on
// eval() evaluates the map at x (in the share inversion, one evaluation of the
// shares), writes the residual the solver compares with tol (it may be +Inf
// where the step is still usable) and the step, and returns false where the
// evaluation failed: shares or a step that are not finite.
//
// The solver runs on the primary thread (eval() may open parallel passes),
// checks for a user interrupt before every evaluation and never raises an R
// warning (under options(warn = 2) a warning would jump past the destructors
// of the working vectors): its outcome is a code, which the R callers turn
// into a classed warning. Its working vectors are allocated once per call
// (the fallback's when it starts), none per iteration.
//
// Methods (FpMethod), each an iteration that takes one accepted evaluation
// per turn (besides the trial points that fail) and races the fallback step
// in the same way (fp_race()):
//   FP_SIMPLE    the damped base step, with stalls (FpSimpleLane);
//   FP_ANDERSON  Anderson acceleration of the base step (FpAndersonLane);
//   FP_SQUAREM   SQUAREM's squared extrapolation of it (FpSquaremLane).
// ============================================================================

enum FpMethod { FP_SIMPLE = 0, FP_ANDERSON = 1, FP_SQUAREM = 2 };

enum FpCode {
  FP_CONVERGED = 0,       // an evaluated point with residual <= tol, returned
  FP_MAX_EVALS = 1,       // max_evals evaluations without one
  FP_NONFINITE = 2,       // every iteration stopped at a trial point that
                          // failed even after backtracking
  FP_NONFINITE_START = 3  // the start's evaluation failed
};

struct FpControl {
  int method = FP_SIMPLE;
  int max_evals = 1000;     // evaluations of the map, the start's included
  double tol = 1e-8;        // on the map's residual
  double damping = 1.0;     // multiplies the step
  int max_backtracks = 10;  // halvings of a step whose trial point failed
  int patience = 0;         // accepted steps without a new best residual that
                            // make a stall (0: no stalls)
  int window = 0;           // accepted base steps within which the best
                            // residual must halve, or the fallback step joins
                            // (0: it never does) or, stopped, goes on (twice
                            // as many while it is stopped); and the
                            // fallback's evaluations within which its own must
                            // halve while it is behind, or it stops
  // Anderson acceleration (FP_ANDERSON)
  int memory = 10;               // differences kept
  double max_condition = 1e8;    // of their least squares: the oldest are
                                 // dropped while its estimate exceeds this
  int restart_patience = 5;      // iterations without a residual below
  double restart_progress = 0.9; // this fraction of the best: a restart
  bool restart_best = false;     // restart from the best point always (Phase
                                 // R's rule; research), not only where the
                                 // residual has grown past twice the best
  // SQUAREM (FP_SQUAREM)
  int scheme = 3;                // the step length's formula: 1, 2 or 3
  double step_min = 1.0;         // the step length's bounds at the start
  double step_max = 1.0;
  double step_factor = 4.0;      // the upper bound's factor: it multiplies the
                                 // bound when alpha reaches it at a point
                                 // accepted, and divides it (not below its
                                 // start) when such a point fails
};

struct FpResult {
  int code = FP_MAX_EVALS;
  int iterations = 0;       // accepted steps, the fallback's included
  int evaluations = 0;      // of the map, the failed ones included
  int backtracks = 0;       // halvings
  int stalls = 0;           // of the base step
  int race = 0;             // evaluations before the fallback joined (0: it did not)
  int resumes = 0;          // times the fallback, stopped as behind, went on
  int restarts = 0;         // of an Anderson iteration; SQUAREM's extrapolated
                            // points that failed
  bool fallback = false;    // the point returned is the fallback's
  double residual = std::numeric_limits<double>::quiet_NaN();  // at the point returned
};

// An iteration's state: its point and step, the best point it has evaluated,
// and its stall rule's.
struct FpLane {
  arma::vec x, step, trial, trial_step, best, best_step;
  double res = 0.0;
  double best_res = std::numeric_limits<double>::infinity();  // none yet
  double relax = 1.0;        // halved at every stall
  int patience = 0;          // doubled at a stall that does not diverge
  int since = 0;             // accepted steps since the best residual improved
  bool returned = false;     // to the current best point, at a stall
  int window = 0;            // the current window's accepted steps (the
                             // fallback's: its evaluations)
  double window_best = 0.0;  // the best residual when it began
  void init(const arma::uword n) {
    x.set_size(n);
    step.set_size(n);
    trial.set_size(n);
    trial_step.set_size(n);
    best.set_size(n);
    best_step.set_size(n);
  }
};

enum FpStep { FP_STEP_OK, FP_STEP_FAILED, FP_STEP_BUDGET };

// One step of a lane: its trial point, retried with the step halved while it
// fails (each an evaluation), then accepted; and, with stalls, the stall rule
// (fp_race()). FP_STEP_FAILED: every trial failed; FP_STEP_BUDGET: no
// evaluation was left for one.
template <class Map>
inline int fp_lane_step(Map& map, FpLane& L, const FpControl& ctl,
                        const bool stalls, FpResult& out) {
  if (out.evaluations >= ctl.max_evals) return FP_STEP_BUDGET;
  double t = ctl.damping * L.relax;
  double trial_res = 0.0;
  for (int b = 0;; ++b) {
    L.trial = L.x + t * L.step;
    Rcpp::checkUserInterrupt();
    ++out.evaluations;
    if (map.eval(L.trial, trial_res, L.trial_step)) break;
    // a spent budget first, then the halvings (each counted when its trial
    // point is evaluated)
    if (out.evaluations >= ctl.max_evals) return FP_STEP_BUDGET;
    if (b >= ctl.max_backtracks) return FP_STEP_FAILED;
    t *= 0.5;
    ++out.backtracks;
  }
  L.x = L.trial;
  L.step = L.trial_step;
  L.res = trial_res;
  ++out.iterations;
  if (L.res < L.best_res) {
    L.best_res = L.res;
    L.best = L.x;
    L.best_step = L.step;
    L.since = 0;
    L.returned = false;
  } else if (stalls && std::isfinite(L.best_res) && ++L.since >= L.patience) {
    if (L.res > 2.0 * L.best_res && !L.returned) {
      L.x = L.best;  // diverging: back to the best point, once
      L.step = L.best_step;
      L.res = L.best_res;
      L.returned = true;
    } else if (L.patience <= std::numeric_limits<int>::max() / 2) {
      L.patience *= 2;
    }
    L.relax *= 0.5;
    L.since = 0;
    ++out.stalls;
  }
  return FP_STEP_OK;
}

// The iterations that race the fallback step (fp_race()), each the method's
// lane a: setup() sizes its vectors, and turn() takes one step, an
// evaluation accepted (more where trial points fail), returning an FpStep.

// FP_SIMPLE: the damped base step, with the stall rule (fp_lane_step()).
struct FpSimpleLane : FpLane {
  void setup(const arma::uword n, const FpControl& ctl) {
    init(n);
    patience = ctl.patience;
  }
  template <class Map>
  int turn(Map& map, const FpControl& ctl, FpResult& out) {
    return fp_lane_step(map, *this, ctl, ctl.patience > 0, out);
  }
};

// FP_ANDERSON: Anderson acceleration, type II (Walker and Ni 2011), of the
// base step f. From the differences of the last memory iterates and of their
// steps, dX = [x_k - x_{k-1}, ...] and dF = [f_k - f_{k-1}, ...] (newest
// first), the least-squares coefficients gamma = argmin ||f_k - dF gamma||
// give
//   x_{k+1} = x_k - dX gamma + damping (f_k - dF gamma),
// the damped base step while there are none. The least squares go through a
// Householder QR of dF, newest column first, so that dropping the oldest
// columns while the triangular factor's condition estimate, max |R_ii| /
// min |R_ii|, exceeds max_condition (or a column depends on the newer ones)
// truncates the factorization: one per iteration, O(n memory^2), on vectors
// allocated once. Safeguards on the residual, never on a merit function
// (after Phase R's solve_anderson_safe() in the BLP engine plan's research):
// a restart when restart_patience iterations bring no residual below
// restart_progress times the best, or when an extrapolated point fails. It
// clears the differences and takes the damped base step from where the
// iteration is, or, where its residual has grown past twice the best (it
// diverges), from the best point, once per best point. A base step (the
// restart's, or the first one, while there are no differences) is halved
// while its trial point fails, and the turn ends in FP_STEP_FAILED when it
// still fails after max_backtracks halvings. Restarting from the best point
// every time (Phase R's rule, restart_best) holds the iteration where a path
// to the solution climbs through larger residuals first, as the plain
// iteration's stalls once did (fp_race()): on Phase A's sweeps of small hard
// designs it failed on four to five times as many.
struct FpAndersonLane : FpLane {
  arma::mat dX, dF, qr;  // n x memory: the differences, slot `head` the
                         // newest, and the factorization's work
  arma::vec beta, gam;   // the reflectors' scalars and the coefficients
  arma::vec qtf;         // Q' f
  int m = 0;             // columns: memory, at most n
  int k = 0;             // differences kept
  int head = 0;          // slot of the newest
  int quiet = 0;         // iterations without the restart's progress
  bool due = false;      // a restart, at the next turn

  void setup(const arma::uword n, const FpControl& ctl) {
    init(n);
    m = static_cast<int>(std::max<arma::uword>(
        1, std::min(n, static_cast<arma::uword>(std::max(ctl.memory, 1)))));
    dX.set_size(n, m);
    dF.set_size(n, m);
    qr.set_size(n, m);
    beta.set_size(m);
    gam.set_size(m);
    qtf.set_size(n);
  }

  int slot(const int c) const { return (head - c + m) % m; }  // c-th newest

  // gamma for the k newest differences, the older ones dropped where their
  // least squares are ill-conditioned (k may fall to 0)
  void solve(const FpControl& ctl) {
    const arma::uword n = x.n_elem;
    for (int c = 0; c < k; ++c) {
      const double* src = dF.colptr(slot(c));
      std::copy(src, src + n, qr.colptr(c));
    }
    int kept = k;
    double rmax = 0.0, rmin = 0.0;
    for (int c = 0; c < k; ++c) {
      // the reflector that zeroes column c below its diagonal (LAPACK's
      // dlarfg): R_cc = r, the rest of v stored in place, v_c = 1
      double* a = qr.colptr(c);
      double sq = 0.0;
      for (arma::uword i = c + 1; i < n; ++i) sq += a[i] * a[i];
      const double alpha = a[c];
      double r = alpha, bc = 0.0;
      if (sq > 0) {
        r = -std::copysign(std::hypot(alpha, std::sqrt(sq)), alpha);
        bc = (r - alpha) / r;
        const double scal = 1.0 / (alpha - r);
        for (arma::uword i = c + 1; i < n; ++i) a[i] *= scal;
        a[c] = r;
      }
      beta[c] = bc;
      const double ar = std::fabs(r);
      rmax = c ? std::max(rmax, ar) : ar;
      rmin = c ? std::min(rmin, ar) : ar;
      if (!(ar > 0) || !std::isfinite(ar) || rmax > ctl.max_condition * rmin) {
        kept = c;
        break;
      }
      for (int d = c + 1; d < k; ++d) {
        double* b = qr.colptr(d);
        double w = b[c];
        for (arma::uword i = c + 1; i < n; ++i) w += a[i] * b[i];
        w *= bc;
        b[c] -= w;
        for (arma::uword i = c + 1; i < n; ++i) b[i] -= w * a[i];
      }
    }
    k = kept;  // the older differences are dropped for good
    if (k == 0) return;
    // gamma = R^-1 (Q' f), the first k entries
    double* q = qtf.memptr();
    std::copy(step.memptr(), step.memptr() + n, q);
    for (int c = 0; c < k; ++c) {
      const double* a = qr.colptr(c);
      double w = q[c];
      for (arma::uword i = c + 1; i < n; ++i) w += a[i] * q[i];
      w *= beta[c];
      q[c] -= w;
      for (arma::uword i = c + 1; i < n; ++i) q[i] -= w * a[i];
    }
    for (int c = k - 1; c >= 0; --c) {
      double v = q[c];
      for (int d = c + 1; d < k; ++d) v -= qr(c, d) * gam[d];
      gam[c] = v / qr(c, c);
    }
  }

  // The trial point accepted: the differences from the current point (none
  // after a restart), the point, and the restart's rule
  void accept(const double trial_res, const bool push, const FpControl& ctl,
              FpResult& out) {
    const arma::uword n = x.n_elem;
    if (push) {
      head = (head + 1) % m;
      double* dx = dX.colptr(head);
      double* df = dF.colptr(head);
      for (arma::uword i = 0; i < n; ++i) {
        dx[i] = trial[i] - x[i];
        df[i] = trial_step[i] - step[i];
      }
      if (k < m) ++k;
    }
    x = trial;
    step = trial_step;
    res = trial_res;
    ++out.iterations;
    // (no count while every residual so far is infinite: no best point yet)
    if (std::isfinite(best_res)) {
      if (res < ctl.restart_progress * best_res) {
        quiet = 0;
      } else {
        ++quiet;
      }
    }
    if (res < best_res) {
      best_res = res;
      best = x;
      best_step = step;
      returned = false;
    }
    if (quiet >= ctl.restart_patience) due = true;
  }

  template <class Map>
  int turn(Map& map, const FpControl& ctl, FpResult& out) {
    if (out.evaluations >= ctl.max_evals) return FP_STEP_BUDGET;
    const arma::uword n = x.n_elem;
    const double t0 = ctl.damping;
    double trial_res = 0.0;
    bool restart = due;
    if (!restart) {
      solve(ctl);
      if (k > 0) {  // the extrapolated point
        double* tr = trial.memptr();
        for (arma::uword i = 0; i < n; ++i) tr[i] = x[i] + t0 * step[i];
        for (int c = 0; c < k; ++c) {
          const double* dx = dX.colptr(slot(c));
          const double* df = dF.colptr(slot(c));
          const double g = gam[c];
          for (arma::uword i = 0; i < n; ++i) tr[i] -= g * (dx[i] + t0 * df[i]);
        }
        Rcpp::checkUserInterrupt();
        ++out.evaluations;
        if (map.eval(trial, trial_res, trial_step)) {
          accept(trial_res, true, ctl, out);
          return FP_STEP_OK;
        }
        if (out.evaluations >= ctl.max_evals) return FP_STEP_BUDGET;
        restart = true;  // it failed
      }
    }
    // The damped base step, halved while its trial point fails: while there
    // are no differences, from where the iteration is; at a restart, which
    // clears them, from there too, or from the best point where the
    // iteration diverges
    bool back = false;
    if (restart) {
      ++out.restarts;
      due = false;
      quiet = 0;
      k = 0;
      back = ctl.restart_best || (res > 2.0 * best_res && !returned);
      if (back) returned = true;
    }
    const arma::vec& from = back ? best : x;
    const arma::vec& from_step = back ? best_step : step;
    double t = t0;
    for (int b = 0;; ++b) {
      trial = from + t * from_step;
      Rcpp::checkUserInterrupt();
      ++out.evaluations;
      if (map.eval(trial, trial_res, trial_step)) break;
      if (out.evaluations >= ctl.max_evals) return FP_STEP_BUDGET;
      if (b >= ctl.max_backtracks) return FP_STEP_FAILED;
      t *= 0.5;
      ++out.backtracks;
    }
    accept(trial_res, !restart, ctl, out);
    return FP_STEP_OK;
  }
};

// FP_SQUAREM: SQUAREM (Varadhan and Roland 2008), the squared extrapolation
// of the damped base step, G(x) = x + damping f(x). From a cycle's start x0,
// x1 = G(x0) and x2 = G(x1) give r = x1 - x0 and v = x2 - 2 x1 + x0, and the
// cycle moves to the extrapolated point
//   x_e = x0 + 2 alpha r + alpha^2 v,
// alpha = ||r|| / ||v|| (scheme 3; scheme 1, -r'v / v'v; scheme 2,
// -r'r / r'v) bounded to [step_min, step_max]: with PyBLP's defaults, 1 and
// 1, never shorter than x2 itself, the point of alpha = 1. The next cycle
// starts at G(x_e), a stabilizing step. When alpha reaches step_max and x_e
// is accepted with a finite residual (x2 only if its step was not halved),
// step_max is multiplied by step_factor: with the defaults, after the first
// cycle's x2, and the bound grows while the extrapolations succeed. An x_e
// whose shares are not finite gives way to x2, and where alpha had reached
// step_max it is divided by step_factor, not below its start; step_min stays
// (Phase R's solve_squarem() in the BLP engine plan's research; residual
// alone, no merit function). One evaluation per turn: x1, then x_e, then the
// stabilizing step; a base step whose trial point fails is halved, and a
// halved step from x0 starts a new cycle where it lands.
struct FpSquaremLane : FpLane {
  arma::vec x0, r, v;       // the cycle's start, x1 - x0 and x2 - 2 x1 + x0
  int phase = 0;            // where the lane is: a cycle's start (0), x1 (1),
                            // the extrapolated point (2)
  double step_lo = 1.0, step_hi = 1.0;  // alpha's bounds now

  void setup(const arma::uword n, const FpControl& ctl) {
    init(n);
    x0.set_size(n);
    r.set_size(n);
    v.set_size(n);
    step_lo = ctl.step_min;
    step_hi = ctl.step_max;
  }

  void accept(const double trial_res, FpResult& out) {
    x = trial;
    step = trial_step;
    res = trial_res;
    ++out.iterations;
    if (res < best_res) {
      best_res = res;
      best = x;
      best_step = step;
    }
  }

  // The damped base step from where the lane is, halved while its trial
  // point fails; `full` whether it was not halved.
  template <class Map>
  int base_step(Map& map, const FpControl& ctl, FpResult& out, bool& full) {
    double t = ctl.damping;
    double trial_res = 0.0;
    full = true;
    for (int b = 0;; ++b) {
      trial = x + t * step;
      Rcpp::checkUserInterrupt();
      ++out.evaluations;
      if (map.eval(trial, trial_res, trial_step)) break;
      if (out.evaluations >= ctl.max_evals) return FP_STEP_BUDGET;
      if (b >= ctl.max_backtracks) return FP_STEP_FAILED;
      t *= 0.5;
      ++out.backtracks;
      full = false;
    }
    accept(trial_res, out);
    return FP_STEP_OK;
  }

  template <class Map>
  int turn(Map& map, const FpControl& ctl, FpResult& out) {
    if (out.evaluations >= ctl.max_evals) return FP_STEP_BUDGET;
    const arma::uword n = x.n_elem;
    bool full = true;
    if (phase == 0) {  // x1 = G(x0)
      x0 = x;
      const int st = base_step(map, ctl, out, full);
      if (st == FP_STEP_OK && full) {
        for (arma::uword i = 0; i < n; ++i) r[i] = x[i] - x0[i];
        phase = 1;
      }
      return st;
    }
    if (phase == 2) {  // the stabilizing step, the next cycle's start
      const int st = base_step(map, ctl, out, full);
      if (st == FP_STEP_OK) phase = 0;
      return st;
    }
    // At x1: the step length from r and v = G(x1) - x1 - r
    double rr = 0.0, vv = 0.0, rv = 0.0;
    for (arma::uword i = 0; i < n; ++i) {
      v[i] = ctl.damping * step[i] - r[i];
      rr += r[i] * r[i];
      vv += v[i] * v[i];
      rv += r[i] * v[i];
    }
    double alpha = ctl.scheme == 1 ? -rv / vv
                   : ctl.scheme == 2 ? -rr / rv
                                     : std::sqrt(rr / vv);
    if (!std::isfinite(alpha)) alpha = 1.0;
    alpha = std::max(step_lo, std::min(step_hi, alpha));
    if (alpha == 1.0) {  // x2 itself
      const int st = base_step(map, ctl, out, full);
      if (st == FP_STEP_OK) {
        if (full && alpha == step_hi && std::isfinite(res)) step_hi *= ctl.step_factor;
        phase = 2;
      }
      return st;
    }
    double trial_res = 0.0;
    for (arma::uword i = 0; i < n; ++i) {
      trial[i] = x0[i] + 2.0 * alpha * r[i] + alpha * alpha * v[i];
    }
    Rcpp::checkUserInterrupt();
    ++out.evaluations;
    if (map.eval(trial, trial_res, trial_step)) {
      accept(trial_res, out);
      // an overshoot to where every outside share underflows (an infinite
      // residual) is taken, but does not lengthen the next one
      if (alpha == step_hi && std::isfinite(trial_res)) step_hi *= ctl.step_factor;
      phase = 2;
      return FP_STEP_OK;
    }
    // The extrapolated point failed: x2 instead, and a shorter bound
    ++out.restarts;
    if (alpha == step_hi) step_hi = std::max(ctl.step_max, step_hi / ctl.step_factor);
    if (out.evaluations >= ctl.max_evals) return FP_STEP_BUDGET;
    const int st = base_step(map, ctl, out, full);
    if (st == FP_STEP_OK) phase = 2;
    return st;
  }
};

// The method's iteration (Lane: FpSimpleLane, the plain iteration,
// FpAndersonLane or FpSquaremLane), safeguarded:
//   * a trial point that fails is retried with the step halved, up to
//     max_backtracks times (each an evaluation; Anderson's extrapolated point
//     restarts instead, and its restart's step is halved; SQUAREM's gives way
//     to the double step x2);
//   * the plain iteration, with patience > 0: a stall (patience accepted
//     steps without a new best residual) halves the relaxation of the later
//     steps; if the residual has grown past twice the best (the steps
//     diverge), the iteration also returns to its best point, once per best
//     point, and otherwise it continues from where it is and doubles its
//     patience: a path whose residual rises on its way (constants climbing
//     back from where their shares barely respond) is shortened but never
//     discarded, so a stall cannot hold the iteration at a point. (While
//     every residual so far is infinite, as when an outside share underflows
//     at the start, there is no best point to return to, and no step
//     counts.) Anderson's iteration restarts instead (FpAndersonLane);
//   * with window > 0 and a map with a fallback step, once window accepted
//     steps of the method fail to halve its best residual (or a trial point
//     fails after every halving), the fallback step joins from the start:
//     the two iterations take steps in turn (each with its halvings), the
//     fallback's without stalls, and the first to meet tol returns. The
//     fallback stops where a window of its evaluations fails to halve its
//     best residual while the method's is lower (it is behind, and slow);
//     the method then has windows of its own, twice as long, and the
//     fallback goes on from where it stopped whenever the method fails to
//     halve its best residual within one (or its trial point fails).
//     Wherever the fallback alone converges from the start in E evaluations,
//     this then converges within T + 2E if the fallback halves its best
//     residual within every window of its evaluations, and within T + 4E
//     otherwise (each stop follows a window of its evaluations, matched by
//     at most as many of the method's steps, and is followed by at most one
//     failing window of the method's, twice as long), T being the
//     evaluations before it joined, however the method fares; not counted
//     are the method's failed trial points after it joined (halvings, and
//     the extrapolated points that failed), and the windows in which
//     the method alone halved its best residual (with tol > 0, at most log2
//     of its residual over tol). Every method takes one evaluation per
//     accepted step besides those failed trial points, so the bound holds
//     for each.
// x holds the start and receives the result: the first evaluated point whose
// residual is at most tol, or, when none is, the evaluated point of smallest
// residual of either iteration (the start itself when max_evals is 0,
// unevaluated).
template <class Lane, class Map>
inline FpResult fp_race(Map& map, arma::vec& x, const FpControl& ctl) {
  FpResult out;
  if (ctl.max_evals <= 0) return out;
  const arma::uword n = map.n();
  Lane a;    // the method's iteration
  FpLane b;  // the fallback's
  a.setup(n, ctl);
  a.x = x;
  Rcpp::checkUserInterrupt();
  ++out.evaluations;
  if (!map.eval(a.x, a.res, a.step)) {
    out.code = FP_NONFINITE_START;
    out.residual = a.res;  // not finite
    return out;
  }
  a.best = a.x;
  a.best_step = a.step;
  a.best_res = a.res;
  a.window_best = a.res;
  const bool can_race = ctl.window > 0 && map.has_fallback();
  // the base step's windows while the fallback is stopped
  const int solo_window = ctl.window <= std::numeric_limits<int>::max() / 2
                              ? 2 * ctl.window
                              : std::numeric_limits<int>::max();
  bool a_on = true, b_on = false, b_failed = false, raced = false;
  for (;;) {
    if (a_on && a.res <= ctl.tol) {
      x = a.x;
      out.code = FP_CONVERGED;
      out.residual = a.res;
      return out;
    }
    if (b_on && b.res <= ctl.tol) {
      x = b.x;
      out.code = FP_CONVERGED;
      out.residual = b.res;
      out.fallback = true;
      return out;
    }
    if (out.evaluations >= ctl.max_evals) break;
    bool join = false;
    if (a_on) {
      const int a_window =
          raced && !b_on && !b_failed ? solo_window : ctl.window;
      const int st = a.turn(map, ctl, out);
      if (st == FP_STEP_BUDGET) break;
      if (st == FP_STEP_FAILED) {
        a_on = false;
        join = true;
      } else if (a.res <= ctl.tol) {
        continue;
      } else if (ctl.window > 0 && ++a.window >= a_window) {
        if (std::isfinite(a.best_res) && std::isfinite(a.window_best) &&
            !(a.best_res <= 0.5 * a.window_best)) {
          join = true;
        }
        a.window = 0;
        a.window_best = a.best_res;
      }
    }
    if (join && can_race && !raced) {
      // The fallback's iteration from the start (which x still holds)
      if (out.evaluations >= ctl.max_evals) break;
      raced = true;
      out.race = out.evaluations;
      b.init(n);
      b.x = x;
      map.use_fallback(true);
      Rcpp::checkUserInterrupt();
      ++out.evaluations;
      b_on = map.eval(b.x, b.res, b.step);
      b_failed = !b_on;
      map.use_fallback(false);
      b.best = b.x;
      b.best_step = b.step;
      b.best_res = b_on ? b.res : arma::datum::inf;
      b.window_best = b.best_res;
      continue;
    }
    if (join && raced && !b_on && !b_failed) {
      if (out.evaluations >= ctl.max_evals) break;  // none left for it
      b_on = true;  // stopped as behind: it goes on from where it stopped
      ++out.resumes;
    }
    if (b_on) {
      map.use_fallback(true);
      const int before = out.evaluations;
      const int st = fp_lane_step(map, b, ctl, false, out);
      map.use_fallback(false);
      if (st == FP_STEP_BUDGET) break;
      if (st == FP_STEP_FAILED) {
        b_on = false;
        b_failed = true;
      } else if (b.res > ctl.tol &&
                 (b.window += out.evaluations - before) >= ctl.window) {
        // the fallback stops where it falls behind: a window of its
        // evaluations that failed to halve its best residual, the base
        // step's being lower; the base step then has a full window, twice
        // as long, to halve its own
        if (a_on && !(b.best_res <= 0.5 * b.window_best) &&
            a.best_res < b.best_res) {
          b_on = false;
          a.window = 0;
          a.window_best = a.best_res;
        }
        b.window = 0;
        b.window_best = b.best_res;
      }
    }
    if (!a_on && !b_on) {
      out.code = FP_NONFINITE;
      break;
    }
  }
  out.fallback = raced && b.best_res < a.best_res;
  x = out.fallback ? b.best : a.best;
  out.residual = out.fallback ? b.best_res : a.best_res;
  return out;
}

template <class Map>
inline FpResult fp_solve(Map& map, arma::vec& x, const FpControl& ctl) {
  switch (ctl.method) {
  case FP_ANDERSON:
    return fp_race<FpAndersonLane>(map, x, ctl);
  case FP_SQUAREM:
    return fp_race<FpSquaremLane>(map, x, ctl);
  case FP_SIMPLE:
  default:
    return fp_race<FpSimpleLane>(map, x, ctl);
  }
}

#endif // CHOICER_FIXED_POINT_HPP
