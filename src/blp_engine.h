#ifndef CHOICER_BLP_ENGINE_HPP
#define CHOICER_BLP_ENGINE_HPP

#include "choicer.h"
#include "choicer_internal.h"
#include "fixed_point.h"
#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

// ============================================================================
// The share inversion (Berry 1994): the constants delta at which a model's
// market shares equal target shares, s(delta) = S, for the multinomial,
// nested and mixed logits (.claude/plans/BLP_ENGINE_PLAN.md). This header
// holds what does not depend on the model: the components of the market and
// their references, the checks of the targets, the residual and stopping
// rule, and the base step; each model file holds an adapter whose pass
// evaluates, at a delta, the shares s_j, the diagonal A_jj of their Jacobian
// ds/ddelta and the co-choice C_jb = sum_u w_u P_uj P_ub with the step's base
// b (u a choice situation, or a situation and draw in the mixed logit).
//
// Identification. With an outside option every inside constant is free.
// Without one, the shares of the alternatives of a connected component of
// the situation-alternative graph (the situations of positive weight) do not
// change when their constants move together, so each component's reference,
// its smallest code, is held at its start (after the start is normalized,
// delta <- delta - delta_1, as before), and the targets of a component must
// add up to the weight of its situations.
//
// The base map (method "simple"): with r = log S - log s, b the outside
// option (or, without one, the component's reference) and r_b its residual,
// num_j = r_j - r_b is the residual of the equations log(s_j / s_b) =
// log(S_j / S_b) (with an outside option, Fukasawa's (2024) outside-share
// correction). Their Jacobian on the free constants is
//   G = diag(1 / s) J + 1 c' / s_b,
// J the share Jacobian and c_j = C_jb = -ds_b / ddelta_j the co-choice of j
// with the base: a part with a dominant diagonal, plus a rank-one coupling
// through the base that moves every free constant's residual together. The
// step treats the two in turn (a two-level step):
//   * the common shift of the component's free constants by tau, from the
//     share-weighted equation sum_j s_j (num_j - h_j tau) = 0, where
//     h = G 1, h_j = C_jb / s_j + sum_k C_kb / s_b, is the exact rate at
//     which num responds to that shift (the share Jacobian's rows sum to
//     zero over every alternative, the base's included, and it is
//     symmetric, so a row's sum over the free constants is the co-choice
//     with the base): it fits the base's share to first order;
//   * the rest by the exact diagonal of G, D_j = A_jj / s_j + C_jb / s_b:
//   step_j = tau + (num_j - h_j tau) / D_j.
// The diagonal alone (the internal map BLP_MAP_JACOBI, step_j = num_j / D_j)
// misjudges the common shift: where the base's share is small and each
// alternative dominates choice situations of its own, its multiplier on that
// shift, 1 - h_j / D_j, falls below -1, and the map diverges where Berry's
// converges.
// Berry's contraction (BLP_MAP_PLAIN, the method "original": Berry, Levinsohn
// and Pakes 1995, eq. 6.8) takes r_j itself (num_j without an outside option,
// the reference held).
// Safeguards of the diagonal (the first and third for the nested logit with
// lambda > 1, whose nest-mates are complements where the nest's choice
// probability is below 1 - 1 / lambda):
//   * an alternative whose co-choice with the base is negative (the base in
//     its nest, without an outside option) has a denominator of at least 1:
//     the diagonal then nets out complementarity and can be near 0 or
//     negative;
//   * a denominator that is not positive and finite takes 1;
//   * in a nest with lambda > 1 the denominator is at least
//     (1 + lambda) / (2 lambda) times the one the same probabilities give
//     with lambda = 1 (A1_jj = sum_u w_u P_uj (1 - P_uj), C1_jb = sum_u w_u
//     P_uj P_ub), the model's floor_scale: a common shift of the nest's
//     constants moves its shares as in a logit, while the exact diagonal is
//     smaller by the within-nest 1 / lambda, so its step would overshoot that
//     shift by up to lambda (the map diverges for lambda above 2); the factor
//     relaxes the two directions equally;
//   * a step is capped at max(|r_j - r_b|, BLP_STEP_CAP) in magnitude: the
//     scaling may lengthen the unscaled step r_j - r_b to that many log
//     units, not beyond (a near-captive alternative far from its target has
//     a denominator near 0);
//   * a component whose outside share underflows to 0 takes the step without
//     the outside share's correction, and without the common shift, that
//     iteration (gamma = 0: r_j over A_jj / s_j; r_b is then infinite).
// And the iteration's own (fixed_point.h): a trial point whose shares are not
// finite halves the step; with method "simple", a stall (no new best residual
// in patience steps) halves the later steps, returning to the best point only
// when the residual has grown past twice the best, once per best point; and
// once a window of steps fails to halve the best residual, Berry's
// contraction joins from the start, the two taking steps in turn, and the
// first to meet tol returns (Berry's stops when a window of its evaluations
// leaves it behind and slow, and goes on whenever the step alone fails to
// halve its best residual within a window twice as long). A diagonal step
// can amplify a small residual of an alternative whose share barely
// responds to its own constant (one that nearly owns its choice situations)
// into a step that wrecks the shares of its neighbours there, and then
// oscillate, and it can carry the constants of alternatives that make up
// the whole choice set of some situations far down, where the residual
// barely responds to their common level; Berry's step is never amplified,
// so wherever Berry's contraction converges the race bounds the method's
// evaluations by four times its own (twice where it halves its best
// residual within every window), plus those before it joined, leaving out
// the step's halvings after it joined and the windows in which the step
// alone halved its best residual (fixed_point.h).
//
// Stopping rule (every method): the largest absolute residual |r_j - r_ref|
// over the inside alternatives without an outside option, |r_j| with one,
// and then also each component's outside share's |r_0c|, at most tol. The
// returned delta is the evaluated point that met it.
// ============================================================================

// The amplification cap of a scaled step, in log units (above).
constexpr double BLP_STEP_CAP = 5.0;

// Targets within this relative distance of their sums (1, or a component's
// weight) are rescaled to them; farther, they are an error.
constexpr double BLP_SUM_TOL = 1.4901161193847656e-08;  // sqrt(DBL_EPSILON)

enum BlpMap { BLP_MAP_REDUCED = 0, BLP_MAP_PLAIN = 1, BLP_MAP_JACOBI = 2 };

struct BlpControl {
  int method = FP_SIMPLE;
  int map = BLP_MAP_REDUCED;
  int max_evals = 1000;
  double tol = 1e-8;
  double damping = 1.0;
  int patience = 0;         // steps without a new best residual: a stall
  int window = 0;           // steps that must halve the best residual, or
                            // Berry's contraction joins, or goes on after a
                            // stop (twice as many then); and its evaluations
                            // that must halve its own (FpControl)
  int memory = 10;          // Anderson acceleration's settings (FpControl)
  double max_condition = 1e8;
  int restart_patience = 5;
  double restart_progress = 0.9;
  bool restart_best = false;
  bool keep_pass = false;   // tests: return the last pass
};

// Components of the situation-alternative graph through the situations of
// positive weight, and what the checks need of them; every share is over the
// total weight W.
struct BlpComponents {
  int n = 0;                   // components
  std::vector<int> of;         // alternative (0-based) -> component, -1 if offered nowhere
  std::vector<int> ref;        // component -> its smallest code (0-based)
  std::vector<int> size;       // component -> alternatives
  std::vector<double> mass;    // component -> weight of its situations
  std::vector<double> offer;   // alternative -> weight of the situations offering it
  std::vector<double> captive; // alternative -> of those offering only it (no outside option)
  bool repeated = false;       // a situation lists a code twice
};

// One serial pass over the rows: union-find over the alternatives (each set's
// root its smallest code), each alternative's offer and captive weight (a code
// listed twice counts its situation once), then each component's weight.
inline BlpComponents blp_components(const ChoiceLayout& lay, const double* w,
                                    const double W, const int J,
                                    const bool outside) {
  BlpComponents c;
  std::vector<int> parent(static_cast<std::size_t>(J));
  for (int j = 0; j < J; ++j) parent[j] = j;
  auto find = [&parent](int a) {
    while (parent[a] != a) {
      parent[a] = parent[parent[a]];
      a = parent[a];
    }
    return a;
  };
  std::vector<double> offer(static_cast<std::size_t>(J), 0.0),
      captive(static_cast<std::size_t>(J), 0.0);
  std::vector<choicer_off> seen(static_cast<std::size_t>(J), -1);
  for (choicer_off t = 0; t < lay.N; ++t) {
    const double wt = w[t];
    if (!(wt > 0)) continue;  // a situation without weight neither offers nor links
    const int* alt = lay.alt + lay.row_off[t];
    const int m = lay.m(t);
    const int a0 = alt[0] - 1;
    bool alone = true;
    for (int r = 0; r < m; ++r) {
      const int a = alt[r] - 1;
      if (seen[a] == t) {
        c.repeated = true;
        continue;
      }
      seen[a] = t;
      offer[a] += wt;
      if (a != a0) {
        alone = false;
        const int ra = find(a), r0 = find(a0);
        if (ra < r0) {
          parent[r0] = ra;
        } else if (r0 < ra) {
          parent[ra] = r0;
        }
      }
    }
    if (!outside && alone) captive[a0] += wt;
  }
  c.of.assign(static_cast<std::size_t>(J), -1);
  for (int j = 0; j < J; ++j) {
    if (!(offer[j] > 0)) continue;
    const int root = find(j);
    if (root == j) {
      c.ref.push_back(j);
      c.size.push_back(0);
      c.of[j] = c.n++;
    } else {
      c.of[j] = c.of[root];  // roots are smallest, so met first
    }
    ++c.size[c.of[j]];
  }
  c.mass.assign(static_cast<std::size_t>(c.n), 0.0);
  for (choicer_off t = 0; t < lay.N; ++t) {
    if (!(w[t] > 0)) continue;
    c.mass[c.of[lay.alt[lay.row_off[t]] - 1]] += w[t];
  }
  for (int k = 0; k < c.n; ++k) c.mass[k] /= W;
  if (c.n == 1) c.mass[0] = 1.0;  // every situation of positive weight
  c.offer.resize(static_cast<std::size_t>(J));
  c.captive.resize(static_cast<std::size_t>(J));
  for (int j = 0; j < J; ++j) {
    c.offer[j] = offer[j] / W;
    c.captive[j] = captive[j] / W;
  }
  return c;
}

// Up to ten 1-based codes, for messages.
inline std::string blp_codes(const std::vector<int>& codes0) {
  std::string s;
  const std::size_t n = std::min<std::size_t>(codes0.size(), 10);
  for (std::size_t k = 0; k < n; ++k) {
    if (k) s += ", ";
    s += std::to_string(codes0[k] + 1);
  }
  if (codes0.size() > n) s += ", ...";
  return s;
}

// The problem the solver works on: the free constants and the targets'
// logarithms, after the checks.
struct BlpProblem {
  int J = 0;
  bool outside = false;
  std::vector<int> free;  // free coordinate k -> alternative (0-based)
  arma::vec log_S;        // inside alternatives
  arma::vec log_S0;       // with an outside option, each component's outside share
};

// The targets' checks, before any pass, and their rescaling: target holds
// the J inside shares, after the outside option's when there is one (as the
// callers take them), already checked for length and for being positive.
inline BlpProblem blp_problem(const BlpComponents& c, const arma::vec& target,
                              const bool outside, const int J) {
  std::vector<int> nowhere;
  for (int j = 0; j < J; ++j) {
    if (c.of[j] < 0) nowhere.push_back(j);
  }
  if (!nowhere.empty()) {
    const bool many = nowhere.size() > 1;
    Rcpp::stop("Alternative%s %s %s offered in no choice situation of positive "
               "weight: %s constant%s not identified, and %s target share%s "
               "cannot be met.", many ? "s" : "", blp_codes(nowhere),
               many ? "are" : "is", many ? "their" : "its",
               many ? "s are" : " is", many ? "their" : "its", many ? "s" : "");
  }
  const int o = outside ? 1 : 0;
  arma::vec S = target.subvec(o, o + J - 1);
  double S0 = outside ? target[0] : 0.0;
  if (outside) {
    const double total = arma::accu(target);
    if (!(std::fabs(total - 1.0) <= BLP_SUM_TOL)) {
      Rcpp::stop("target_shares must sum to 1, the outside option's share "
                 "included; they sum to %.15g.", total);
    }
    S /= total;
    S0 /= total;
  } else {
    std::vector<double> sum(static_cast<std::size_t>(c.n), 0.0);
    for (int j = 0; j < J; ++j) sum[c.of[j]] += S[j];
    for (int k = 0; k < c.n; ++k) {
      if (!(std::fabs(sum[k] - c.mass[k]) <= BLP_SUM_TOL * c.mass[k])) {
        if (c.n == 1) {
          Rcpp::stop("target_shares must sum to 1; they sum to %.15g.", sum[k]);
        }
        Rcpp::stop("The alternatives that share choice situations with "
                   "alternative %d have target shares summing to %.15g, but "
                   "their choice situations weigh %.15g of the total: without "
                   "an outside option the shares of such a group of "
                   "alternatives add up to its situations' weight.",
                   c.ref[k] + 1, sum[k], c.mass[k]);
      }
    }
    for (int j = 0; j < J; ++j) S[j] *= c.mass[c.of[j]] / sum[c.of[j]];
  }
  BlpProblem pb;
  pb.J = J;
  pb.outside = outside;
  if (outside) {
    pb.log_S0.set_size(static_cast<arma::uword>(c.n));
    if (c.n == 1) {
      if (!(S0 > 0)) {
        Rcpp::stop("The target shares leave no share to the outside option.");
      }
      pb.log_S0[0] = std::log(S0);
    } else {
      std::vector<double> inside(static_cast<std::size_t>(c.n), 0.0);
      for (int j = 0; j < J; ++j) inside[c.of[j]] += S[j];
      for (int k = 0; k < c.n; ++k) {
        const double s0 = c.mass[k] - inside[k];
        if (!(s0 > 0)) {
          Rcpp::stop("The target shares leave no share to the outside option in "
                     "the choice situations of alternative %d: their inside "
                     "targets sum to %.15g of the %.15g those situations weigh.",
                     c.ref[k] + 1, inside[k], c.mass[k]);
        }
        pb.log_S0[k] = std::log(s0);
      }
    }
  }
  for (int j = 0; j < J; ++j) {
    const int k = c.of[j];
    if (!outside && c.size[k] == 1) continue;  // alone in its component: S_j is its mass
    if (!(S[j] < c.offer[j])) {
      Rcpp::stop("Alternative %d's target share (%.15g) is not below the weight "
                 "of the choice situations that offer it (%.15g): no constants "
                 "attain it.", j + 1, S[j], c.offer[j]);
    }
    if (!outside && !(S[j] > c.captive[j])) {
      Rcpp::stop("Alternative %d's target share (%.15g) is not above the weight "
                 "of the choice situations that offer it alone (%.15g): no "
                 "constants attain it.", j + 1, S[j], c.captive[j]);
    }
  }
  pb.log_S = arma::log(S);
  for (int j = 0; j < J; ++j) {
    if (outside || c.ref[c.of[j]] != j) pb.free.push_back(j);
  }
  return pb;
}

// A pass's output, over the total weight: inside shares, the diagonal of
// ds/ddelta and the co-choice with the base (J each), each component's
// outside share (with an outside option), and, from a model with a
// floor_scale, the diagonal and co-choice its probabilities give with lambda
// = 1 (A1, C1; empty otherwise).
struct BlpPass {
  arma::vec s, A, C, s0, A1, C1;
  BlpPass(const int J, const int n_comp)
      : s(J, arma::fill::zeros), A(J, arma::fill::zeros), C(J, arma::fill::zeros),
        s0(static_cast<arma::uword>(n_comp), arma::fill::zeros) {}
};

struct BlpResult {
  FpResult fp;
  double fallbacks = 0;  // safeguards applied: denominators raised or replaced,
                         // steps capped, gamma = 0 steps taken (a count that
                         // can pass INT_MAX: coordinates times evaluations)
};

// The share inversion as a map of the free constants (fixed_point.h). Every
// evaluation is one pass of the model at the full delta.
template <class Model>
struct BlpMapEval {
  Model& model;
  const BlpComponents& comp;
  const BlpProblem& pb;
  const int map;
  arma::vec& delta;  // J; the free constants written in, the references held
  const arma::vec* floor_scale;  // J, 0 where no floor; nullptr: none
  BlpPass pass;
  arma::vec r, r0;
  arma::vec num, den, h;                  // per free constant
  arma::vec sb, csum, wnum, wh, tau;      // per component (the common shift)
  double fallbacks = 0;
  bool berry = false;  // the fallback step, Berry's, while the solver asks

  BlpMapEval(Model& model_, const BlpComponents& comp_, const BlpProblem& pb_,
             const int map_, arma::vec& delta_, const arma::vec* floor_scale_)
      : model(model_), comp(comp_), pb(pb_), map(map_), delta(delta_),
        floor_scale(floor_scale_), pass(pb_.J, comp_.n), r(pb_.J),
        r0(pb_.outside ? comp_.n : 0), num(pb_.free.size()),
        den(pb_.free.size()), h(pb_.free.size()), sb(comp_.n), csum(comp_.n),
        wnum(comp_.n), wh(comp_.n), tau(comp_.n) {}

  arma::uword n() const { return pb.free.size(); }
  bool has_fallback() const { return map != BLP_MAP_PLAIN; }
  void use_fallback(const bool on) { berry = on; }

  bool eval(const arma::vec& x, double& res, arma::vec& step) {
    const int J = pb.J;
    const std::size_t nf = pb.free.size();
    for (std::size_t k = 0; k < nf; ++k) delta[pb.free[k]] = x[k];
    model.pass(delta, pass);
    bool ok = true;
    res = 0.0;
    for (int j = 0; j < J; ++j) r[j] = pb.log_S[j] - std::log(pass.s[j]);
    for (int j = 0; j < J; ++j) {
      const double rn = pb.outside ? r[j] : r[j] - r[comp.ref[comp.of[j]]];
      if (!std::isfinite(rn)) ok = false;
      res = std::max(res, std::fabs(rn));
    }
    if (pb.outside) {
      for (int k = 0; k < comp.n; ++k) {
        r0[k] = pb.log_S0[k] - std::log(pass.s0[k]);
        if (std::isnan(r0[k])) ok = false;
        res = std::max(res, std::fabs(r0[k]));  // +Inf where s0 underflowed
      }
    }
    if (!ok) {
      res = arma::datum::inf;
      return false;
    }
    if (map == BLP_MAP_PLAIN || berry) {
      for (std::size_t k = 0; k < nf; ++k) {
        const int j = pb.free[k];
        step[k] = pb.outside ? r[j] : r[j] - r[comp.ref[comp.of[j]]];
      }
      return true;
    }
    // The residuals against the base, the diagonal and its safeguards; and
    // each component's base share (0: none this iteration) and co-choices
    sb.zeros();
    csum.zeros();
    for (int c = 0; c < comp.n; ++c) {
      if (pb.outside) {
        if (pass.s0[c] > 0 && std::isfinite(r0[c])) sb[c] = pass.s0[c];
      } else {
        sb[c] = pass.s[comp.ref[c]];
      }
    }
    for (std::size_t k = 0; k < nf; ++k) {
      const int j = pb.free[k];
      const int c = comp.of[j];
      const double a_s = pass.A[j] / pass.s[j];
      double d;
      if (sb[c] > 0) {
        num[k] = r[j] - (pb.outside ? r0[c] : r[comp.ref[c]]);
        d = a_s + pass.C[j] / sb[c];
        csum[c] += pass.C[j];
      } else {
        num[k] = r[j];  // gamma = 0: the outside share underflowed
        d = a_s;
        ++fallbacks;
      }
      if (sb[c] > 0 && pass.C[j] < 0 && !(d >= 1.0)) {
        d = 1.0;
        ++fallbacks;
      }
      if (!(d > 0) || !std::isfinite(d)) {
        d = 1.0;
        ++fallbacks;
      }
      if (sb[c] > 0 && floor_scale && (*floor_scale)[j] > 0) {
        const double d1 =
            (*floor_scale)[j] * (pass.A1[j] / pass.s[j] + pass.C1[j] / sb[c]);
        if (d1 > d) {
          d = d1;
          ++fallbacks;
        }
      }
      den[k] = d;
    }
    // The common shift: the rate h_j at which residual j responds to it, and
    // its size tau from sum_j s_j (num_j - h_j tau) = 0
    tau.zeros();
    if (map == BLP_MAP_REDUCED) {
      wnum.zeros();
      wh.zeros();
      for (std::size_t k = 0; k < nf; ++k) {
        const int j = pb.free[k];
        const int c = comp.of[j];
        if (!(sb[c] > 0)) continue;
        h[k] = pass.C[j] / pass.s[j] + csum[c] / sb[c];
        wnum[c] += pass.s[j] * num[k];
        wh[c] += pass.s[j] * h[k];
      }
      for (int c = 0; c < comp.n; ++c) {
        // (none without a base, or without free constants: an alternative
        // alone in its component, with no outside option)
        if (!(sb[c] > 0) || (!pb.outside && comp.size[c] < 2)) continue;
        // (wh = csum (1 + sum_j s_j / s_b) > 0, csum being the base's own
        // derivative: lost only to underflow, and the shift with it)
        const double t = wnum[c] / wh[c];
        if (wh[c] > 0 && std::isfinite(t)) {
          tau[c] = t;
        } else {
          ++fallbacks;
        }
      }
    }
    for (std::size_t k = 0; k < nf; ++k) {
      const int c = comp.of[pb.free[k]];
      double st = tau[c] != 0.0 ? tau[c] + (num[k] - h[k] * tau[c]) / den[k]
                                : num[k] / den[k];
      const double cap = std::max(std::fabs(num[k]), BLP_STEP_CAP);
      if (std::fabs(st) > cap) {
        st = std::copysign(cap, st);
        ++fallbacks;
      }
      if (!std::isfinite(st)) {
        res = arma::datum::inf;
        return false;
      }
      step[k] = st;
    }
    return true;
  }
};

// The solver's settings (fixed_point.h) of a control list's.
inline FpControl blp_fp_control(const BlpControl& ctl) {
  FpControl fc;
  fc.method = ctl.method;
  fc.max_evals = ctl.max_evals;
  fc.tol = ctl.tol;
  fc.damping = ctl.damping;
  fc.patience = ctl.patience;
  fc.window = ctl.window;
  fc.memory = ctl.memory;
  fc.max_condition = ctl.max_condition;
  fc.restart_patience = ctl.restart_patience;
  fc.restart_progress = ctl.restart_progress;
  fc.restart_best = ctl.restart_best;
  return fc;
}

// Solve for delta (J entries, in: the start, already normalized; out: the
// result) with the model's passes.
template <class Model>
inline BlpResult blp_solve(Model& model, const BlpComponents& comp,
                           const BlpProblem& pb, arma::vec& delta,
                           const BlpControl& ctl, BlpPass* last_pass = nullptr,
                           const arma::vec* floor_scale = nullptr) {
  BlpMapEval<Model> map(model, comp, pb, ctl.map, delta, floor_scale);
  arma::vec x(pb.free.size());
  for (std::size_t k = 0; k < pb.free.size(); ++k) x[k] = delta[pb.free[k]];
  const FpControl fc = blp_fp_control(ctl);
  BlpResult out;
  if (pb.free.empty()) {
    out.fp.code = FP_CONVERGED;  // every constant held: nothing to solve
    out.fp.residual = 0.0;
    return out;
  }
  out.fp = fp_solve(map, x, fc);
  for (std::size_t k = 0; k < pb.free.size(); ++k) delta[pb.free[k]] = x[k];
  out.fallbacks = map.fallbacks;
  if (last_pass) *last_pass = map.pass;
  return out;
}

// The control list of the internal entry points (.blp_control() in R).
inline BlpControl blp_control(const Rcpp::List& control) {
  BlpControl c;
  const std::string iteration = Rcpp::as<std::string>(control["iteration"]);
  if (iteration == "simple") {
    c.method = FP_SIMPLE;
  } else if (iteration == "anderson") {
    c.method = FP_ANDERSON;
  } else {
    Rcpp::stop("Internal error: unknown BLP iteration \"%s\".", iteration);
  }
  const std::string map = Rcpp::as<std::string>(control["map"]);
  if (map == "reduced") {
    c.map = BLP_MAP_REDUCED;
  } else if (map == "plain") {
    c.map = BLP_MAP_PLAIN;
  } else if (map == "jacobi") {
    c.map = BLP_MAP_JACOBI;
  } else {
    Rcpp::stop("Internal error: unknown BLP map \"%s\".", map);
  }
  c.max_evals = Rcpp::as<int>(control["max_evals"]);
  c.tol = Rcpp::as<double>(control["tol"]);
  c.damping = Rcpp::as<double>(control["damping"]);
  c.patience = Rcpp::as<int>(control["patience"]);
  c.window = Rcpp::as<int>(control["window"]);
  c.memory = Rcpp::as<int>(control["memory"]);
  c.max_condition = Rcpp::as<double>(control["max_condition"]);
  c.restart_patience = Rcpp::as<int>(control["restart_patience"]);
  c.restart_progress = Rcpp::as<double>(control["restart_progress"]);
  c.restart_best = Rcpp::as<bool>(control["restart_best"]);
  c.keep_pass = Rcpp::as<bool>(control["keep_pass"]);
  return c;
}

// The internal entry points' result: delta as the n x 1 matrix RcppArmadillo
// returns, and the solver's record.
inline Rcpp::List blp_result(const arma::vec& delta, const BlpResult& res,
                             const BlpControl& ctl, const BlpPass& last) {
  Rcpp::List out = Rcpp::List::create(
      Rcpp::Named("delta") = delta,
      Rcpp::Named("code") = res.fp.code,
      Rcpp::Named("iterations") = res.fp.iterations,
      Rcpp::Named("evaluations") = res.fp.evaluations,
      Rcpp::Named("backtracks") = res.fp.backtracks,
      Rcpp::Named("stalls") = res.fp.stalls,
      Rcpp::Named("race") = res.fp.race,
      Rcpp::Named("resumes") = res.fp.resumes,
      Rcpp::Named("restarts") = res.fp.restarts,
      Rcpp::Named("berry") = res.fp.fallback,
      Rcpp::Named("residual") = res.fp.residual,
      Rcpp::Named("fallbacks") = res.fallbacks);
  if (ctl.keep_pass) {
    Rcpp::List pass = Rcpp::List::create(Rcpp::Named("s") = last.s,
                                         Rcpp::Named("A") = last.A,
                                         Rcpp::Named("C") = last.C,
                                         Rcpp::Named("s0") = last.s0);
    if (last.A1.n_elem) {
      pass["A1"] = last.A1;
      pass["C1"] = last.C1;
    }
    out["pass"] = pass;
  }
  return out;
}

// delta finite (the start); weights finite and non-negative.
inline void blp_check_start(const arma::vec& delta, const arma::vec& weights) {
  if (!delta.is_finite()) {
    Rcpp::stop("delta must be finite (the start of the share inversion).");
  }
  for (arma::uword i = 0; i < weights.n_elem; ++i) {
    if (!(weights[i] >= 0) || !std::isfinite(weights[i])) {
      Rcpp::stop("weights must be finite and non-negative (weight %d is %g).",
                 static_cast<int>(i) + 1, weights[i]);
    }
  }
}

#endif // CHOICER_BLP_ENGINE_HPP
