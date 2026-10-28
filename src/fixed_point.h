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
// eval() evaluates the map at x (in the share inversion, one evaluation of the
// shares), writes the residual the solver compares with tol (it may be +Inf
// where the step is still usable) and the base step, and returns false where
// the evaluation failed: shares or a step that are not finite.
//
// The solver runs on the primary thread (eval() may open parallel passes),
// checks for a user interrupt before every evaluation and never raises an R
// warning (under options(warn = 2) a warning would jump past the destructors
// of the working vectors): its outcome is a code, which the R callers turn
// into a classed warning. Its working vectors are allocated once per call,
// none per iteration. Phase B of the BLP engine plan adds Anderson
// acceleration and SQUAREM on the same maps (FpMethod).
// ============================================================================

enum FpMethod { FP_SIMPLE = 0 };

enum FpCode {
  FP_CONVERGED = 0,       // an evaluated point with residual <= tol, returned
  FP_MAX_EVALS = 1,       // max_evals evaluations without one
  FP_NONFINITE = 2,       // a trial point failed even after backtracking
  FP_NONFINITE_START = 3  // the start's evaluation failed
};

struct FpControl {
  int method = FP_SIMPLE;
  int max_evals = 1000;     // evaluations of the map, the start's included
  double tol = 1e-8;        // on the map's residual
  double damping = 1.0;     // multiplies the base step
  int max_backtracks = 10;  // halvings of a step whose trial point failed
  int patience = 0;         // accepted steps without a new best residual
                            // before a restart from the best point with the
                            // step's relaxation halved; 0: no restarts
};

struct FpResult {
  int code = FP_MAX_EVALS;
  int iterations = 0;       // accepted steps
  int evaluations = 0;      // of the map, the failed ones included
  int backtracks = 0;       // halvings
  int restarts = 0;         // from the best point, after patience ran out
  double residual = std::numeric_limits<double>::quiet_NaN();  // at the point returned
};

// The plain iteration, safeguarded: a trial point that fails is retried with
// the step halved, up to max_backtracks times (each an evaluation), before
// the solver stops; and, with patience > 0, an iteration whose best residual
// has not improved in patience accepted steps (it diverges, or cycles)
// restarts from its best point with the relaxation of its steps halved for
// the rest of the call (again at each further restart): a base step that
// overshoots some direction contracts once it is short enough. (While every residual so far is infinite, as when
// an outside share underflows at the start, there is no best point to return
// to, and no step counts.) x holds the start and receives the result: the
// first evaluated point whose residual is at most tol, or, when none is, the
// evaluated point of smallest residual (the start itself when max_evals is
// 0, unevaluated).
template <class Map>
inline FpResult fp_simple(Map& map, arma::vec& x, const FpControl& ctl) {
  FpResult out;
  if (ctl.max_evals <= 0) return out;
  const arma::uword n = map.n();
  arma::vec step(n), trial(n), trial_step(n), best(n), best_step(n);
  double res = 0.0;
  Rcpp::checkUserInterrupt();
  ++out.evaluations;
  if (!map.eval(x, res, step)) {
    out.code = FP_NONFINITE_START;
    out.residual = res;  // not finite
    return out;
  }
  best = x;
  best_step = step;
  double best_res = res;
  double relax = 1.0;  // halved at every restart
  int since = 0;       // accepted steps since the best residual improved
  for (;;) {
    if (res <= ctl.tol) {
      out.code = FP_CONVERGED;
      out.residual = res;
      return out;
    }
    if (out.evaluations >= ctl.max_evals) {
      out.code = FP_MAX_EVALS;
      break;
    }
    double t = ctl.damping * relax;
    double trial_res = 0.0;
    bool ok = false;
    for (int b = 0;; ++b) {
      trial = x + t * step;
      Rcpp::checkUserInterrupt();
      ++out.evaluations;
      ok = map.eval(trial, trial_res, trial_step);
      if (ok || b >= ctl.max_backtracks || out.evaluations >= ctl.max_evals) break;
      t *= 0.5;
      ++out.backtracks;
    }
    if (!ok) {
      out.code = out.evaluations >= ctl.max_evals ? FP_MAX_EVALS : FP_NONFINITE;
      break;
    }
    x = trial;
    step = trial_step;
    res = trial_res;
    ++out.iterations;
    if (res < best_res) {
      best_res = res;
      best = x;
      best_step = step;
      since = 0;
    } else if (ctl.patience > 0 && std::isfinite(best_res) &&
               ++since >= ctl.patience) {
      x = best;
      step = best_step;
      res = best_res;
      relax *= 0.5;
      since = 0;
      ++out.restarts;
    }
  }
  x = best;
  out.residual = best_res;
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
