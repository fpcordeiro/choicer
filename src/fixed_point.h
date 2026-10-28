#ifndef CHOICER_FIXED_POINT_HPP
#define CHOICER_FIXED_POINT_HPP

#include "choicer.h"
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
// (the fallback's when it starts), none per iteration. Phase B of the BLP
// engine plan adds Anderson acceleration and SQUAREM on the same maps
// (FpMethod).
// ============================================================================

enum FpMethod { FP_SIMPLE = 0 };

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
};

struct FpResult {
  int code = FP_MAX_EVALS;
  int iterations = 0;       // accepted steps, the fallback's included
  int evaluations = 0;      // of the map, the failed ones included
  int backtracks = 0;       // halvings
  int stalls = 0;           // of the base step
  int race = 0;             // evaluations before the fallback joined (0: it did not)
  int resumes = 0;          // times the fallback, stopped as behind, went on
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
// (fp_simple). FP_STEP_FAILED: every trial failed; FP_STEP_BUDGET: no
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

// The plain iteration, safeguarded:
//   * a trial point that fails is retried with the step halved, up to
//     max_backtracks times (each an evaluation);
//   * with patience > 0, a stall (patience accepted steps without a new best
//     residual) halves the relaxation of the later steps; if the residual has
//     grown past twice the best (the steps diverge), the iteration also
//     returns to its best point, once per best point, and otherwise it
//     continues from where it is and doubles its patience: a path whose
//     residual rises on its way (constants climbing back from where their
//     shares barely respond) is shortened but never discarded, so a stall
//     cannot hold the iteration at a point. (While every residual so far is
//     infinite, as when an outside share underflows at the start, there is
//     no best point to return to, and no step counts.);
//   * with window > 0 and a map with a fallback step, once window accepted
//     steps fail to halve the best residual (or a trial point fails after
//     every halving), the fallback step joins from the start: the two
//     iterations take steps in turn (each with its halvings), the
//     fallback's without stalls, and the first to meet tol returns. The
//     fallback stops where a window of its evaluations fails to halve its
//     best residual while the base step's is lower (it is behind, and slow);
//     the base step then has windows of its own, twice as long, and the
//     fallback goes on from where it stopped whenever the base step fails to
//     halve its best residual within one (or its trial point fails).
//     Wherever the fallback alone converges from the start in E evaluations,
//     this then converges within T + 2E if the fallback halves its best
//     residual within every window of its evaluations, and within T + 4E
//     otherwise (each stop follows a window of its evaluations, matched by
//     at most as many of the base step's steps, and is followed by at most
//     one failing window of the base step's, twice as long), T being the
//     evaluations before it joined, however the base step fares; not counted
//     are halvings of the base step's trial points after it joined, and the
//     windows in which the base step alone halved its best residual (with
//     tol > 0, at most log2 of its residual over tol).
// x holds the start and receives the result: the first evaluated point whose
// residual is at most tol, or, when none is, the evaluated point of smallest
// residual of either iteration (the start itself when max_evals is 0,
// unevaluated).
template <class Map>
inline FpResult fp_simple(Map& map, arma::vec& x, const FpControl& ctl) {
  FpResult out;
  if (ctl.max_evals <= 0) return out;
  const arma::uword n = map.n();
  FpLane a, b;  // the base step's iteration and the fallback's
  a.init(n);
  a.patience = ctl.patience;
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
      const int st = fp_lane_step(map, a, ctl, ctl.patience > 0, out);
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
  case FP_SIMPLE:
  default:
    return fp_simple(map, x, ctl);
  }
}

#endif // CHOICER_FIXED_POINT_HPP
