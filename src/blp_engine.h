#ifndef CHOICER_BLP_ENGINE_HPP
#define CHOICER_BLP_ENGINE_HPP

#include "choicer.h"
#include "choicer_internal.h"
#include "choicer_derivs.h"
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

// The exp-space pass (BlpExpCache, below): the largest distance from the
// cache's origin, in log units, at which an evaluation uses the cache; and
// the share of the memory limit (kernel_memory_limit()) that caches may take,
// less what the call already holds.
constexpr double BLP_EXP_SPAN = 50.0;
constexpr double BLP_CACHE_SHARE = 0.5;

enum BlpMap { BLP_MAP_REDUCED = 0, BLP_MAP_PLAIN = 1, BLP_MAP_JACOBI = 2 };
enum BlpPassKind { BLP_PASS_LOG = 0, BLP_PASS_EXP = 1 };
enum BlpCompensate { BLP_COMP_OFF = 0, BLP_COMP_ON = 1, BLP_COMP_AUTO = 2 };

// Compensated sums (BlpControl::compensate, "auto") above this many choice
// situations per thread: a share summed over 2.5e7 situations on one thread
// had a rounding floor of 1.7e-13 (BLP engine plan, Phase B), growing about
// as N^0.7, so about 5e-14 at this count.
constexpr double BLP_COMP_SITUATIONS = 4194304.0;  // 2^22

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
  int scheme = 3;           // SQUAREM's settings (FpControl)
  double step_min = 1.0;
  double step_max = 1.0;
  double step_factor = 4.0;
  bool keep_pass = false;   // tests: return the last pass
  int pass = BLP_PASS_LOG;  // the passes' arithmetic (BlpExpCache)
  double cache_bytes = arma::datum::inf;  // the multinomial and nested logits'
                                          // exp-space cache: most bytes
  arma::vec origin;         // tests: where the exp-space cache is formed first
                            // (empty: at the first evaluation)
  int compensate = BLP_COMP_OFF;  // the share sums: plain, Kahan's, or Kahan's
                                  // above BLP_COMP_SITUATIONS per thread
};

// ----------------------------------------------------------------------------
// Compensated sums
//
// A pass adds each situation's weighted probabilities into its thread's share
// accumulators, so a share sums N / T terms, and its rounding grows with
// them: on one thread, near the default tol at 10^9 rows (BLP engine plan,
// Phase B). With compensate, the passes add the shares (the outside options'
// included) by Kahan's compensated summation, the threads' totals too, and
// the inversion's own sums over the situations (the weights' total, the
// components' weights) likewise, which leaves a share's error at a few
// rounding units whatever N. The diagonal and the co-choices, which only scale
// the steps, keep plain sums.
// ----------------------------------------------------------------------------

// x added to s, whose running compensation is c (Kahan 1965): the sum is s - c.
// (A build with -ffast-math may fold the compensation away: the sum is then
// the plain one.)
inline void kahan_add(double& s, double& c, const double x) {
  const double y = x - c;
  const double t = s + y;
  c = (t - s) - y;
  s = t;
}

// Whether a pass over N situations on T threads compensates its share sums
// (each thread's sum runs over its N / T situations).
inline bool blp_compensate(const BlpControl& ctl, const double N, const int T) {
  if (ctl.compensate == BLP_COMP_ON) return true;
  return ctl.compensate == BLP_COMP_AUTO &&
         N > BLP_COMP_SITUATIONS * static_cast<double>(std::max(T, 1));
}

// Whether the inversion's serial sums over the N situations (the weights'
// total, the components' weights, and the targets' bounds: the weights of
// the situations offering each alternative, and offering it alone) are
// compensated. They run on one thread whatever the passes', and with an
// outside option an error in the total moves every share alike, a floor
// under the residual.
inline bool blp_compensate_serial(const BlpControl& ctl, const double N) {
  return blp_compensate(ctl, N, 1);
}

// The weights' total by compensated summation in two interleaved chains
// (the dependent adds of one chain would cost four times the plain sum's),
// with shares_denominator()'s check.
inline double blp_weight_sum(const arma::vec& weights) {
  double s1 = 0.0, c1 = 0.0, s2 = 0.0, c2 = 0.0;
  const arma::uword n = weights.n_elem;
  arma::uword i = 0;
  for (; i + 1 < n; i += 2) {
    kahan_add(s1, c1, weights[i]);
    kahan_add(s2, c2, weights[i + 1]);
  }
  if (i < n) kahan_add(s1, c1, weights[i]);
  double s = 0.0, c = 0.0;
  kahan_add(s, c, s1 - c1);
  kahan_add(s, c, s2 - c2);
  s -= c;
  if (s <= 0) {
    Rcpp::stop("Error: Sum of weights must be positive.");
  }
  return s;
}

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
// Compensated (kahan), the offer and captive weights are too, and each is
// taken as a share of its component's weight summed in the same order (the
// total of positive weights for a single component): an alternative that
// every situation of its component offers then bounds its target by the
// component's weight exactly, as the plain sums, added in one order, do. (A
// plain offer weight over a compensated total fell below a target that its
// constants attain.)
inline BlpComponents blp_components(const ChoiceLayout& lay, const double* w,
                                    const double W, const int J,
                                    const bool outside, const bool kahan = false) {
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
  const std::size_t J_k = kahan ? static_cast<std::size_t>(J) : 0;
  std::vector<double> offer_c(J_k, 0.0), captive_c(J_k, 0.0);
  double total = 0.0, total_c = 0.0;  // (compensated: the positive weights)
  std::vector<choicer_off> seen(static_cast<std::size_t>(J), -1);
  for (choicer_off t = 0; t < lay.N; ++t) {
    const double wt = w[t];
    if (!(wt > 0)) continue;  // a situation without weight neither offers nor links
    if (kahan) kahan_add(total, total_c, wt);
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
      if (kahan) {
        kahan_add(offer[a], offer_c[a], wt);
      } else {
        offer[a] += wt;
      }
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
    if (!outside && alone) {
      if (kahan) {
        kahan_add(captive[a0], captive_c[a0], wt);
      } else {
        captive[a0] += wt;
      }
    }
  }
  if (kahan) {
    for (int j = 0; j < J; ++j) {
      offer[j] -= offer_c[j];
      captive[j] -= captive_c[j];
    }
    total -= total_c;
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
  std::vector<double> mass_c(kahan ? static_cast<std::size_t>(c.n) : 0, 0.0);
  for (choicer_off t = 0; t < lay.N; ++t) {
    if (kahan && c.n == 1) break;  // (one component's weight is 1, set below)
    if (!(w[t] > 0)) continue;
    if (kahan) {
      const int k = c.of[lay.alt[lay.row_off[t]] - 1];
      kahan_add(c.mass[k], mass_c[k], w[t]);
    } else {
      c.mass[c.of[lay.alt[lay.row_off[t]] - 1]] += w[t];
    }
  }
  if (kahan) {
    for (int k = 0; k < c.n; ++k) c.mass[k] -= mass_c[k];
  }
  // (compensated: each component's weight as its offer weights were summed)
  std::vector<double> raw;
  if (kahan) raw = c.n == 1 ? std::vector<double>{total} : c.mass;
  for (int k = 0; k < c.n; ++k) c.mass[k] /= W;
  if (c.n == 1) c.mass[0] = 1.0;  // every situation of positive weight
  c.offer.resize(static_cast<std::size_t>(J));
  c.captive.resize(static_cast<std::size_t>(J));
  for (int j = 0; j < J; ++j) {
    const int k = c.of[j];
    if (kahan && k >= 0) {
      c.offer[j] = c.mass[k] * (offer[j] / raw[k]);
      c.captive[j] = c.mass[k] * (captive[j] / raw[k]);
    } else {
      c.offer[j] = offer[j] / W;
      c.captive[j] = captive[j] / W;
    }
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

// ----------------------------------------------------------------------------
// The exp-space pass
//
// A model's pass forms every choice probability as a quotient of
// exponentials, one per row (and draw), each exp(utility - shift) for a shift
// per choice situation (per offered nest in the nested logit, per draw in the
// mixed logit). Only delta changes between evaluations, so a model may keep
// those exponentials from its pass at an origin delta0, its cache, and
// evaluate at delta by multiplying the row's by a_j = exp((delta_j - delta0_j)
// / scale_j) (scale_j = 1, or the dissimilarity of j's nest): J exponentials
// per evaluation instead of one or more per row, and the same quotients, so
// the same probabilities up to rounding.
//
// The cache is used only within BLP_EXP_SPAN of its origin, d = max_j
// |delta_j - delta0_j| / scale_j <= 50, where the factors and the sums stay
// far from overflow and a row whose exponential was subnormal at the origin
// (known to a few bits) carries at most DBL_MIN e^(2 d) into its probability.
// An evaluation from the cache is kept only if every share it gives, the
// outside options' included, is finite and at least F(d) = 2^53 DBL_MIN
// e^(2 d) (about 5.4e-249 at the span), where that error, summed over a
// situation's rows, is below rounding. Otherwise, and beyond the span, the
// model runs its usual pass at delta and stores its exponentials on the way:
// delta becomes the origin, and that evaluation is the usual pass bit for
// bit. A cache whose own pass gave a share below F(span) is not used until a
// later usual pass gives none (conservative: near such an origin, shares are
// tiny, and most evaluations from it would end in the usual pass). A delta
// that is not finite takes the usual pass and leaves the cache as it was.
// ----------------------------------------------------------------------------

// F(d): the smallest share an evaluation at distance d from the cache's
// origin may give (above).
inline double blp_exp_floor(const double d) {
  return std::ldexp(std::exp(2.0 * d), -969);  // 2^53 DBL_MIN e^(2 d)
}

struct BlpExpCache {
  arma::vec d0;        // the origin
  arma::vec factor;    // exp((delta - delta0) / scale) of the evaluation in hand
  arma::vec scale;     // empty: 1 for every alternative
  bool outside = false;
  bool valid = false;  // a cache formed at d0
  bool tiny = false;   // its pass gave a share below F(BLP_EXP_SPAN)
  double evals = 0;    // evaluations from the cache
  double forms = 0;    // usual passes that formed it
  double floors = 0;   // evaluations from it discarded by the floor

  BlpExpCache(const int J, const bool outside_)
      : d0(static_cast<arma::uword>(J)), factor(static_cast<arma::uword>(J)),
        outside(outside_) {}

  // Whether the evaluation at delta uses the cache; then factor holds its
  // factors and d its distance from the origin
  bool use(const arma::vec& delta, double& d) {
    if (!valid || tiny) return false;
    const arma::uword J = delta.n_elem;
    const bool scaled = scale.n_elem == J;
    d = 0.0;
    for (arma::uword j = 0; j < J; ++j) {
      double x = delta[j] - d0[j];
      if (scaled) x /= scale[j];
      const double ax = std::fabs(x);
      if (!(ax <= BLP_EXP_SPAN)) return false;  // (NaN included)
      d = std::max(d, ax);
      factor[j] = x;
    }
    for (arma::uword j = 0; j < J; ++j) factor[j] = std::exp(factor[j]);
    return true;
  }

  // Whether a pass's shares are all finite and at least lo
  bool above(const BlpPass& out, const double lo) const {
    for (arma::uword j = 0; j < out.s.n_elem; ++j) {
      if (!(out.s[j] >= lo) || !std::isfinite(out.s[j])) return false;
    }
    if (outside) {
      for (arma::uword k = 0; k < out.s0.n_elem; ++k) {
        if (!(out.s0[k] >= lo) || !std::isfinite(out.s0[k])) return false;
      }
    }
    return true;
  }

  // The usual pass at delta formed the cache
  void formed(const arma::vec& delta, const BlpPass& out) {
    d0 = delta;
    valid = true;
    tiny = !above(out, blp_exp_floor(BLP_EXP_SPAN));
    ++forms;
  }
};

// One evaluation of a model that keeps an exp-space cache (above): its
// sweep<true, false>() evaluates from the cache, sweep<false, true>() is its
// usual pass storing the cache, and sweep<false, false>() its usual pass.
template <class Model>
inline void blp_exp_pass(Model& model, BlpExpCache& xc, const arma::vec& delta,
                         BlpPass& out) {
  double d = 0.0;
  if (xc.use(delta, d)) {
    model.template sweep<true, false>(delta, out);
    if (xc.above(out, blp_exp_floor(d))) {
      ++xc.evals;
      return;
    }
    ++xc.floors;
  }
  if (!delta.is_finite()) {
    model.template sweep<false, false>(delta, out);
    return;
  }
  model.template sweep<false, true>(delta, out);
  xc.formed(delta, out);
}

// The bytes the share inversion may keep in a cache: budget, capped by
// BLP_CACHE_SHARE of the memory limit (kernel_memory_limit(): the option
// choicer.max_memory, else the physical memory) less held, the bytes the call
// already holds; no cap where the limit is unknown (0) or infinite. A budget
// that is not positive (NaN included) keeps nothing. It reads an option (and
// stops on a malformed one): call it on the primary thread, and only where a
// cache may be kept.
inline double blp_cache_budget(const double budget, const double held) {
  if (!(budget > 0)) return 0.0;
  const char* source = "";
  const double limit = kernel_memory_limit(source);
  if (!(limit > 0) || std::isinf(limit)) return budget;
  const double cap = BLP_CACHE_SHARE * limit - held;
  return cap > 0 ? std::min(budget, cap) : 0.0;
}

// The exp-space cache's origin given for tests (control$origin), normalized
// as the solvers normalize a start: J finite constants, less the first's
// without an outside option.
inline arma::vec blp_origin(const BlpControl& ctl, const int J, const bool outside) {
  if (ctl.origin.n_elem != static_cast<arma::uword>(J) || !ctl.origin.is_finite()) {
    Rcpp::stop("Internal error: control$origin must hold %d finite constants.", J);
  }
  arma::vec o = ctl.origin;
  if (!outside && J > 0) o -= ctl.origin[0];
  return o;
}

// Whether the passes compensated their sums, on an internal entry point's
// result (with a compensate setting other than "off").
inline void blp_comp_record(Rcpp::List& out, const BlpControl& ctl, const bool kahan) {
  if (ctl.compensate != BLP_COMP_OFF) out["compensated"] = kahan;
}

// The exp-space cache's record on an internal entry point's result.
inline void blp_exp_record(Rcpp::List& out, const BlpControl& ctl, const bool kept,
                           const BlpExpCache& xc) {
  if (ctl.pass != BLP_PASS_EXP) return;
  out["exp_kept"] = kept;
  out["exp_evals"] = xc.evals;
  out["exp_forms"] = xc.forms;
  out["exp_floors"] = xc.floors;
}

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
  fc.scheme = ctl.scheme;
  fc.step_min = ctl.step_min;
  fc.step_max = ctl.step_max;
  fc.step_factor = ctl.step_factor;
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
  } else if (iteration == "squarem") {
    c.method = FP_SQUAREM;
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
  c.scheme = Rcpp::as<int>(control["scheme"]);
  c.step_min = Rcpp::as<double>(control["step_min"]);
  c.step_max = Rcpp::as<double>(control["step_max"]);
  c.step_factor = Rcpp::as<double>(control["step_factor"]);
  c.keep_pass = Rcpp::as<bool>(control["keep_pass"]);
  // (optional: lists made before these settings existed lack them)
  if (control.containsElementNamed("pass")) {
    const std::string pass = Rcpp::as<std::string>(control["pass"]);
    if (pass == "log") {
      c.pass = BLP_PASS_LOG;
    } else if (pass == "exp") {
      c.pass = BLP_PASS_EXP;
    } else {
      Rcpp::stop("Internal error: unknown BLP pass \"%s\".", pass);
    }
  }
  if (control.containsElementNamed("cache_bytes")) {
    c.cache_bytes = Rcpp::as<double>(control["cache_bytes"]);
  }
  if (control.containsElementNamed("origin") && !Rf_isNull(control["origin"])) {
    c.origin = Rcpp::as<arma::vec>(control["origin"]);
  }
  if (control.containsElementNamed("compensate")) {
    const std::string comp = Rcpp::as<std::string>(control["compensate"]);
    if (comp == "off") {
      c.compensate = BLP_COMP_OFF;
    } else if (comp == "on") {
      c.compensate = BLP_COMP_ON;
    } else if (comp == "auto") {
      c.compensate = BLP_COMP_AUTO;
    } else {
      Rcpp::stop("Internal error: unknown BLP compensate \"%s\".", comp);
    }
  }
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
