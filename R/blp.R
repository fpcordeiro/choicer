# The share inversion (Berry 1994; Berry, Levinsohn and Pakes 1995) behind
# blp() and the expert functions blp_contraction(), nl_blp_contraction() and
# mxl_blp_contraction(): one engine with one convention for every model
# (src/blp_engine.h, src/fixed_point.h; the models' passes in their files).
# The internal entry points mnl_blp_solve(), nl_blp_solve() and
# mxl_blp_solve() return the solver's record; the functions here validate the
# scalar settings, call them, and turn an unconverged record into the classed
# warning choicer_blp_unconverged, raised from R.

#' Invert market shares for the constants of a multinomial logit
#'
#' Finds the alternative-specific constants \eqn{\delta} at which the model's
#' market shares \eqn{s_j(\delta) = \sum_i w_i P_{ij}(\delta) / \sum_i w_i}
#' equal target shares \eqn{S_j}, holding the other parameters fixed: the
#' inversion of Berry (1994) and Berry, Levinsohn and Pakes (1995)
#' \doi{10.2307/2171802} behind calibration to aggregate shares and
#' counterfactual share matching. \code{\link{blp}()} runs it on a fitted
#' model; these functions take the model's inputs directly.
#'
#' @param delta Start: one constant per inside alternative (\eqn{J}). Without
#'   an outside option the first alternative's constant is the normalization:
#'   the start is shifted so that it is 0.
#' @param target_shares Target shares: \eqn{J} of them, or \eqn{J + 1} with the
#'   outside option's first. They must be positive and sum to 1 (to within
#'   about \eqn{1.5 \times 10^{-8}}, the remainder rescaled away).
#' @param X sum(M) x K design matrix with covariates. M\[i] x K matrix for
#'   individual i.
#' @param beta K x 1 vector with model parameters.
#' @param alt_idx sum(M) x 1 vector with indices of alternatives within each
#'   choice set; 1-based indexing.
#' @param M N x 1 vector with number of alternatives for each individual.
#' @param weights N x 1 vector with weights for each observation (finite and
#'   non-negative).
#' @param include_outside_option Whether there is an outside option, with
#'   utility normalized to 0 (its rows are not in the data).
#' @param tol Convergence tolerance on the largest absolute log-share
#'   residual \eqn{|\log S_j - \log s_j(\delta)|} over the alternatives (the
#'   outside option's included; see Details). Default 1e-12. Rounding in the
#'   sums of the shares over choice situations sets a floor under the
#'   residual that grows with their number: on one thread about
#'   \eqn{2 \times 10^{-13}} at \eqn{10^8} rows and near \eqn{10^{-12}} at
#'   \eqn{10^9} (lower with more threads). Where it reaches \code{tol}, raise
#'   \code{tol}.
#' @param max_iter The most evaluations of the shares (default 1000).
#' @param method The iteration: \code{"anderson"} (the default), Anderson
#'   acceleration of the two-level step of Details; \code{"simple"}, that
#'   step alone; \code{"original"}, Berry, Levinsohn and Pakes's contraction;
#'   or \code{"squarem"}, SQUAREM on the step.
#' @param control A list of the method's settings. \code{"anderson"} takes
#'   \code{memory}, how many differences between successive iterates it keeps
#'   (a whole number from 1 to 100, cut to the number of free constants;
#'   default 10), and \code{max_condition}, the largest condition estimate it
#'   accepts for their least squares (at least 1; default \eqn{10^8}).
#'   \code{"squarem"} takes \code{scheme}, its step length's formula (1, 2 or
#'   3; default 3), \code{step_min} and \code{step_max}, the step length's
#'   bounds at the start (positive; default 1 and 1), and
#'   \code{step_factor}, what the upper bound is multiplied by when the step
#'   length reaches it and divided by when such a step fails (at least 1;
#'   default 4): PyBLP's defaults. The other methods take none.
#' @details With \eqn{r = \log S - \log s(\delta)}, the iteration solves
#'   \eqn{\rho_j = r_j - r_b = 0} for every free constant, the log-share
#'   residual against a base: the outside option (Fukasawa's 2024
#'   outside-share correction) or, without one, the reference of the
#'   alternative's group. (The alternatives that share choice situations,
#'   directly or through others, form a connected component, whose constants
#'   the shares identify only up to a common level, so its smallest code's
#'   constant is held at its start.) The default, \code{method = "anderson"},
#'   extrapolates the step below, which \code{method = "simple"} takes alone.
#'   Each iteration evaluates the shares once and takes a step in two parts.
#'   First a common shift \eqn{\tau} of the group's free constants, to which
#'   each \eqn{\rho_j} responds at the rate
#'   \eqn{-h_j}, \eqn{h_j = C_{jb}/s_j + \sum_k C_{kb}/s_b} exactly, where
#'   \eqn{C_{jb} = \sum_i w_i P_{ij} P_{ib} / \sum_i w_i} is the share-weighted
#'   co-choice of \eqn{j} and the base; \eqn{\tau} solves
#'   \eqn{\sum_j s_j (\rho_j - h_j \tau) = 0}, to first order the shift that
#'   matches the base's share. Then the rest of each residual, divided by the
#'   exact diagonal of its derivative, \eqn{A_{jj}/s_j + C_{jb}/s_b}, where
#'   \eqn{A_{jj}} is the diagonal of the share Jacobian
#'   \eqn{\partial s/\partial \delta}. The diagonal alone would misjudge the
#'   common shift, and diverge where the base's share is small and each
#'   alternative dominates choice situations of its own. Berry, Levinsohn and
#'   Pakes's (1995) contraction (\code{blp(method = "original")}) moves each
#'   constant by its own residual, \eqn{r_j} (without an outside option
#'   \eqn{r_j - r_b}, the reference held). A step is at most 5 log units long,
#'   or \eqn{|\rho_j|} where that is longer; a step whose shares are not finite
#'   is halved, up to ten times. Under \code{"simple"}, an iteration whose
#'   best residual has not improved for 30 iterations (a stall) halves its
#'   later steps; it returns
#'   to its best point only when its residual at the stall exceeds twice the
#'   best (it diverges), once per best point, and otherwise goes on from where
#'   it is, each such stall doubling the iterations the next one waits for,
#'   so that constants climbing back from where their shares barely respond
#'   are not pulled back. And once 100 iterations fail to halve the best
#'   residual, or a trial point's shares stay non-finite after every halving
#'   (the step then stops, and the contraction goes on alone), Berry,
#'   Levinsohn and Pakes's contraction joins from the start: the two take
#'   steps in turn, each with its halvings, and the first to meet \code{tol}
#'   returns. The contraction stops if 100 of its own evaluations fail to
#'   halve its best residual while the step's best is lower; the step then
#'   has 200 iterations at a time to halve its own best residual, and
#'   whenever it fails to (or its trial point's shares stay non-finite), the
#'   contraction goes on from where it stopped. A divided step can amplify
#'   the small residual of an alternative that nearly owns its choice
#'   situations into a step that wrecks its neighbours' shares there, and
#'   then oscillate, or carry the constants of alternatives that make up the
#'   whole choice set of some choice situations so far down that the
#'   residual barely responds to their common level; the contraction's step
#'   is never amplified. So wherever the contraction converges in \eqn{E}
#'   evaluations, the step converges within \eqn{T + 4E}, and within
#'   \eqn{T + 2E} where the contraction halves its best residual within every
#'   100 of its evaluations, \eqn{T} being the evaluations before the
#'   contraction joined; the bound leaves out the step's halvings after the
#'   contraction joined and any 200 iterations in which the step alone
#'   halved its best residual.
#'
#'   \code{method = "anderson"} accelerates the same step (Anderson
#'   acceleration of type II; Walker and Ni 2011). Write \eqn{f(\delta)} for
#'   the step at \eqn{\delta}, and \eqn{\Delta X} and \eqn{\Delta F} for the
#'   last \code{memory} differences between successive iterates and between
#'   their steps. Each iteration finds the coefficients \eqn{\gamma} whose
#'   combination of past changes in the step best cancels the current step,
#'   \eqn{\min_\gamma \| f - \Delta F \gamma \|_2}, and moves to
#'   \eqn{\delta - \Delta X \gamma + f - \Delta F \gamma}, with the step's part
#'   multiplied by \code{damping} for the nested logit: the capped step
#'   above, extrapolated. On a linear map, undamped and with unlimited
#'   memory, this is essentially GMRES. The oldest differences are dropped
#'   while the least squares' condition estimate exceeds
#'   \code{max_condition}. When five iterations bring no residual below 0.9
#'   times the best so far, or the shares at an extrapolated point are not
#'   finite, the iteration clears its differences and takes the step from
#'   where it is, or from its best point where its residual has grown past
#'   twice the best (once per best point), halving that step while its
#'   shares are not finite. Where a component's outside share underflows,
#'   it takes the step above unextrapolated until the share returns. It
#'   races Berry, Levinsohn and Pakes's
#'   contraction as the step does, with the same windows; each of its
#'   iterations takes one evaluation of the shares besides the trial points
#'   whose shares are not finite, which the bound leaves out as it leaves out
#'   the step's halvings, so the same bound holds.
#'
#'   \code{method = "squarem"} extrapolates the same step by SQUAREM
#'   (Varadhan and Roland 2008), PyBLP's default (Conlon and Gortmaker 2020).
#'   With \eqn{G(\delta) = \delta + f(\delta)} (\eqn{f} multiplied by
#'   \code{damping} for the nested logit), each cycle takes
#'   \eqn{\delta_1 = G(\delta_0)} and \eqn{\delta_2 = G(\delta_1)}, sets
#'   \eqn{r = \delta_1 - \delta_0} and \eqn{v = \delta_2 - 2\delta_1 + \delta_0},
#'   moves to \eqn{\delta_0 + 2\alpha r + \alpha^2 v} with
#'   \eqn{\alpha = \|r\|/\|v\|} (\code{scheme} 3) bounded below by
#'   \code{step_min} and above by \code{step_max}, and takes one more step
#'   \eqn{G} from there: three evaluations of the shares per cycle. With the
#'   default bounds the first cycle's point is \eqn{\delta_2} itself. When
#'   \eqn{\alpha} reaches \code{step_max} at a point whose residual is
#'   finite (for \eqn{\delta_2}, a step not halved), the bound is
#'   multiplied by \code{step_factor}, so the extrapolation lengthens while
#'   it succeeds; where the shares at an extrapolated point are not finite,
#'   the cycle takes \eqn{\delta_2} instead, and a bound that \eqn{\alpha}
#'   had reached is divided by \code{step_factor}, not below its start.
#'   \code{step_min} stays. It races the contraction too, one evaluation per
#'   turn besides the trial points that fail, so the same bound holds.
#'
#'   The iteration stops at the first evaluation whose largest absolute
#'   log-share residual is at most \code{tol} (residuals relative to the
#'   component's reference without an outside option; with one, the outside
#'   option's own residual included), and returns that point. The error in
#'   \eqn{\delta} is then of order \code{tol} times the conditioning of the
#'   share equations, which segmented markets (choice sets with few
#'   alternatives in common) can raise to hundreds or more. When it does not
#'   converge within \code{max_iter} evaluations, or stops at shares that are
#'   not finite, it warns (class \code{choicer_blp_unconverged}; the
#'   condition's \code{convergence} element holds the status, iterations,
#'   evaluations and residual) and returns the evaluated point of smallest
#'   residual (the start when \code{max_iter = 0}).
#'
#'   Before any evaluation it checks that the targets can be met: every
#'   alternative is offered in a choice situation of positive weight; the
#'   targets of each component add up to the weight of its situations (to 1
#'   with an outside option, and then leave a positive share to the outside
#'   option in each component); and each target lies strictly between the
#'   weight of the situations that offer the alternative alone and the weight
#'   of those that offer it. These conditions are necessary, not sufficient: a
#'   target that passes them may still not be attainable, and then the
#'   iteration does not converge.
#' @returns The constants, an n x 1 matrix (\eqn{J} rows, inside alternatives
#'   only); without an outside option the first is 0.
#' @seealso \code{\link{blp}()}, \code{\link{nl_blp_contraction}()},
#'   \code{\link{mxl_blp_contraction}()}
#' @references Berry, S. T. (1994). Estimating discrete-choice models of
#'   product differentiation. \emph{RAND Journal of Economics}, 25(2),
#'   242--262.
#'
#'   Berry, S., Levinsohn, J., and Pakes, A. (1995). Automobile prices in
#'   market equilibrium. \emph{Econometrica}, 63(4), 841--890.
#'
#'   Fukasawa, T. (2024). Fast and simple inner-loop algorithms of static /
#'   dynamic BLP estimations. arXiv:2404.04494.
#'
#'   Conlon, C., and Gortmaker, J. (2020). Best practices for differentiated
#'   products demand estimation with PyBLP. \emph{RAND Journal of Economics},
#'   51(4), 1108--1161.
#'
#'   Varadhan, R., and Roland, C. (2008). Simple and globally convergent
#'   methods for accelerating the convergence of any EM algorithm.
#'   \emph{Scandinavian Journal of Statistics}, 35(2), 335--353.
#'
#'   Walker, H. F., and Ni, P. (2011). Anderson acceleration for fixed-point
#'   iterations. \emph{SIAM Journal on Numerical Analysis}, 49(4), 1715--1735.
#' @examples
#' \donttest{
#' library(data.table)
#' set.seed(42)
#' N <- 50; J <- 3
#' dt <- data.table(id = rep(1:N, each = J), alt = rep(1:J, N))
#' dt[, `:=`(x1 = rnorm(.N), x2 = rnorm(.N))]
#' dt[, choice := 0L]
#' dt[, choice := sample(c(1L, rep(0L, J - 1))), by = id]
#' fit <- run_mnlogit(dt, "id", "alt", "choice", c("x1", "x2"))
#' beta <- coef(fit)[fit$param_map$beta]
#' delta <- blp_contraction(rep(0, J), rep(1/J, J), fit$data$X,
#'   beta, fit$data$alt_idx, fit$data$M, fit$data$weights)
#' delta
#' }
#' @export
blp_contraction <- function(delta, target_shares, X, beta, alt_idx, M, weights,
                            include_outside_option = FALSE, tol = 1e-12,
                            max_iter = 1000L, method = "anderson",
                            control = list()) {
  m <- .blp_method(method, control)
  ctl <- .blp_control(tol, max_iter, method = m$method, settings = m$settings)
  res <- mnl_blp_solve(delta, target_shares, X, beta, alt_idx, M, weights,
                       include_outside_option, ctl)
  .blp_unconverged(res, ctl, "blp_contraction()")
  res$delta
}

#' Invert market shares for the constants of a nested logit
#'
#' \code{\link{blp_contraction}()} for the nested logit: the constants at which
#' the nested logit's shares equal the targets, with \eqn{\beta} and the
#' dissimilarities held fixed.
#'
#' @inheritParams blp_contraction
#' @param beta K x 1 vector with fixed coefficients.
#' @param lambda The dissimilarity of every nest of \code{nest_idx}
#'   (\code{max(nest_idx)} entries, a singleton nest's 1), each positive and
#'   finite.
#' @param alt_idx sum(M) x 1 vector with indices of alternatives; 1-based
#'   indexing.
#' @param nest_idx J x 1 vector with nest indices for each alternative;
#'   1-based indexing.
#' @param include_outside_option Whether there is an outside option, with
#'   utility normalized to 0 and dissimilarity 1.
#' @param damping Multiplies each step, in (0, 1] (default 1). The steps are
#'   scaled by the share Jacobian's diagonal, whose within-nest term
#'   \eqn{(1 - P(j|g))/\lambda_g} keeps them stable under strong nesting, so
#'   damping is seldom needed.
#' @details The step's diagonal is built from the nested logit's derivative,
#'   \eqn{A_{jj} = \sum_i w_i P_{ij}((1 - P_i(j|g))/\lambda_g + P_i(j|g) -
#'   P_{ij}) / \sum_i w_i} for alternative \eqn{j} in nest \eqn{g}, and,
#'   without an outside option, its exact co-choice with the reference
#'   \eqn{b}, \eqn{\sum_i w_i P_{ib}((1 - \lambda_g)/\lambda_g\, P_i(j|g)\,
#'   1\{b \in g\} + P_{ij}) / \sum_i w_i}, which also give the exact rate at
#'   which the residuals respond to the common shift (the nested logit's share
#'   Jacobian is symmetric, with rows summing to zero). In a nest with
#'   \eqn{\lambda_g > 1} nest-mates are complements where the nest's choice
#'   probability is below \eqn{1 - 1/\lambda_g}, and two safeguards apply. An alternative whose co-choice with the reference is
#'   negative (the reference in its nest) has a diagonal of at least 1: its
#'   exact one nets out the complementarity and can be near 0 or negative.
#'   And the diagonal is at least \eqn{(1 + \lambda_g)/(2\lambda_g)} times the
#'   one the same probabilities give with \eqn{\lambda_g = 1}: a common shift
#'   of a nest's constants moves its shares as in a logit, while the exact
#'   diagonal is smaller by the within-nest \eqn{1/\lambda_g}, so its step
#'   would overshoot that shift by up to \eqn{\lambda_g} and, for
#'   \eqn{\lambda_g} above 2, diverge. See \code{\link{blp_contraction}()} for
#'   the iteration, its stopping rule and its checks.
#' @returns The constants, an n x 1 matrix (\eqn{J} rows); without an outside
#'   option the first is 0.
#' @seealso \code{\link{blp}()}, \code{\link{blp_contraction}()}
#' @examples
#' \donttest{
#' library(data.table)
#' set.seed(42)
#' N <- 50; J <- 4
#' dt <- data.table(id = rep(1:N, each = J), alt = rep(1:J, N))
#' dt[, `:=`(x1 = rnorm(.N), x2 = rnorm(.N))]
#' dt[, nest := ifelse(alt <= 2, "A", "B")]
#' dt[, choice := 0L]
#' dt[, choice := sample(c(1L, rep(0L, J - 1))), by = id]
#' fit <- run_nestlogit(dt, "id", "alt", "choice", c("x1", "x2"), "nest")
#' beta <- coef(fit)[fit$param_map$beta]
#' lambda <- rep(1, length(unique(fit$data$nest_idx)))
#' lambda[as.integer(names(which(table(fit$data$nest_idx) > 1)))] <-
#'   coef(fit)[fit$param_map$lambda]
#' target <- as.numeric(predict(fit, type = "shares"))
#' delta <- nl_blp_contraction(rep(0, J), rev(target), fit$data$X, beta,
#'   lambda, fit$data$alt_idx, fit$data$nest_idx, fit$data$M, fit$data$weights)
#' delta
#' }
#' @export
nl_blp_contraction <- function(delta, target_shares, X, beta, lambda, alt_idx,
                               nest_idx, M, weights,
                               include_outside_option = FALSE, damping = 1.0,
                               tol = 1e-12, max_iter = 1000L,
                               method = "anderson", control = list()) {
  m <- .blp_method(method, control)
  ctl <- .blp_control(tol, max_iter, damping, method = m$method,
                      settings = m$settings)
  res <- nl_blp_solve(delta, target_shares, X, beta, lambda, alt_idx, nest_idx,
                      M, weights, include_outside_option, ctl)
  .blp_unconverged(res, ctl, "nl_blp_contraction()")
  res$delta
}

#' Invert market shares for the constants of a mixed logit
#'
#' \code{\link{blp_contraction}()} for the mixed logit: the constants at which
#' the simulated shares \eqn{s_j(\delta) = \sum_i w_i S^{-1} \sum_r
#' P_{ijr}(\delta) / \sum_i w_i} equal the targets, over the draws given, with
#' the other parameters fixed. The step's \eqn{A_{jj}} and \eqn{C_{jb}}
#' average each draw's terms, \eqn{A_{jj} = \sum_i w_i S^{-1} \sum_r
#' P_{ijr}(1 - P_{ijr}) / \sum_i w_i} and \eqn{C_{jb} = \sum_i w_i S^{-1}
#' \sum_r P_{ijr} P_{ibr} / \sum_i w_i}.
#'
#' @inheritParams blp_contraction
#' @param delta Start: one constant per inside alternative (\eqn{J}, the
#'   largest code in \code{alt_idx}); without an outside option also the
#'   \eqn{J - 1} free ones, after the first's normalized 0.
#' @param target_shares Target shares: \eqn{J}, or \eqn{J + 1} with the outside
#'   option's first; positive, summing to 1.
#' @param W design matrix for random coefficients; sum(M_i) x K_w or J x K_w
#' @param beta K_x vector with fixed coefficients
#' @param mu K_w vector with mean parameters (raw, transformed if log-normal);
#'   used as given when \code{rc_mean = FALSE}
#' @param L_params Cholesky parameters vector
#' @param alt_idx sum(M) x 1 vector with indices of alternatives; 1-based
#'   indexing
#' @param M N x 1 vector with number of alternatives for each individual
#' @param weights N x 1 vector with weights for each observation
#' @param eta_draws Array with draws, K_w x S x N (an integer or logical array
#'   is converted to double first)
#' @param rc_dist K_w vector indicating distribution (0=normal, 1=log-normal)
#' @param rc_correlation whether random coefficients are correlated
#' @param rc_mean whether mu parameters represent means (TRUE) or are zero
#'   (FALSE)
#' @param include_outside_option whether outside option is included
#' @param gen_seed Integer. \code{< 0} (default) uses the materialized
#'   \code{eta_draws} cube; \code{>= 0} forms the draws on the fly (see
#'   \code{gen_scramble}), with this master seed for the digit permutations of
#'   \code{gen_scramble = 1}.
#' @param gen_scramble Integer mode of the draws formed on the fly: \code{0} =
#'   identity permutations (plain Halton) through the generator's inverse normal
#'   CDF, \code{1} = seeded position-wise digit permutations, \code{2} = identity
#'   permutations through R's \code{qnorm()}, the store-mode draws of
#'   \code{\link{get_halton_normals}} formed without its cube (bit for bit where
#'   randtoolbox and choicer are compiled with the same floating-point
#'   contraction). Other values are an error.
#' @param gen_S Integer number of draws per individual, used only when
#'   \code{gen_seed >= 0}.
#' @details With \code{gen_seed >= 0} the draws are formed again at every
#'   evaluation of the shares. \code{\link{blp}()} on a fitted model runs the
#'   same inversion but keeps them across evaluations within its
#'   \code{keep_draws_bytes} (see \code{\link{blp.choicer_mxl}}). See
#'   \code{\link{blp_contraction}()} for the iteration, its stopping rule and
#'   its checks.
#' @returns The constants, an n x 1 matrix (\eqn{J} rows); without an outside
#'   option the first is 0.
#' @seealso \code{\link{blp}()}, \code{\link{blp_contraction}()}
#' @examples
#' \donttest{
#' library(data.table)
#' set.seed(42)
#' N <- 50; J <- 3
#' dt <- data.table(id = rep(1:N, each = J), alt = rep(1:J, N))
#' dt[, `:=`(x1 = rnorm(.N), w1 = rnorm(.N))]
#' dt[, choice := 0L]
#' dt[, choice := sample(c(1L, rep(0L, J - 1))), by = id]
#' d <- prepare_mxl_data(dt, "id", "alt", "choice", "x1", "w1")
#' eta <- get_halton_normals(50, d$N, ncol(d$W))
#' fit <- run_mxlogit(input_data = d, eta_draws = eta)
#' pm <- fit$param_map
#' delta <- mxl_blp_contraction(rep(0, J), rep(1/J, J), d$X, d$W,
#'   coef(fit)[pm$beta], rep(0, ncol(d$W)), coef(fit)[pm$sigma],
#'   d$alt_idx, d$M, d$weights, eta, rc_dist = rep(0L, ncol(d$W)),
#'   rc_correlation = FALSE, rc_mean = FALSE)
#' delta
#' }
#' @export
mxl_blp_contraction <- function(delta, target_shares, X, W, beta, mu, L_params,
                                alt_idx, M, weights, eta_draws, rc_dist,
                                rc_correlation = TRUE, rc_mean = FALSE,
                                include_outside_option = FALSE, tol = 1e-12,
                                max_iter = 1000L, gen_seed = -1L,
                                gen_scramble = 1L, gen_S = 0L,
                                method = "anderson", control = list()) {
  m <- .blp_method(method, control)
  ctl <- .blp_control(tol, max_iter, method = m$method, settings = m$settings)
  # The kernels read a cube of doubles in place; a cube of another type would
  # be read from a temporary that R may free during the call
  if (is.integer(eta_draws) || is.logical(eta_draws)) {
    storage.mode(eta_draws) <- "double"
  }
  res <- mxl_blp_solve(delta, target_shares, X, W, beta, mu, L_params, alt_idx,
                       M, weights, eta_draws, rc_dist, rc_correlation, rc_mean,
                       include_outside_option, ctl, gen_seed, gen_scramble,
                       gen_S, 0)
  .blp_unconverged(res, ctl, "mxl_blp_contraction()")
  res$delta
}

#' The share inversion's settings, checked, as the internal solvers take them
#'
#' @param tol,max_iter,damping As in \code{blp_contraction()} and
#'   \code{nl_blp_contraction()}, but for \code{tol}'s default: 1e-8 here,
#'   the internal callers' and the tests', where the exports' is 1e-12.
#' @param method The method, as blp() names it: \code{"simple"} (the
#'   two-level step, its stalls and its race with Berry's contraction),
#'   \code{"original"} (Berry, Levinsohn and Pakes's contraction),
#'   \code{"anderson"} (Anderson acceleration of the two-level step) or
#'   \code{"squarem"} (SQUAREM on it), the last two racing the contraction as
#'   \code{"simple"} does. (This default, \code{"simple"}, is the internal
#'   callers' and the tests'; the exports and blp(), whose default is
#'   \code{"anderson"}, pass theirs.) The first two run the engine's simple
#'   iteration (\code{iteration}) and differ in the map, the stalls and the
#'   race.
#' @param settings The method's settings, checked (\code{.blp_method()});
#'   the defaults of \code{.BLP_SETTINGS} for those not given.
#' @param map The base map, by default the method's: \code{"reduced"} (the
#'   two-level step), \code{"plain"} (Berry's) or \code{"jacobi"} (the
#'   diagonal alone; tests).
#' @param patience Accepted steps without a new best residual that make a
#'   stall, which halves the later steps (0: none); by default
#'   \code{.BLP_PATIENCE} for the simple iteration on the two-level map, else
#'   0.
#' @param window Accepted steps within which the best residual must halve, or
#'   Berry's contraction joins from the start (0: never), or, stopped as
#'   behind, goes on (twice as many while it is stopped); and the
#'   contraction's evaluations within which its own must halve while the
#'   method's is lower, or it stops. By default \code{.BLP_WINDOW} on the
#'   two-level map, else 0.
#' @param keep_pass Whether the solver returns its last pass (tests).
#' @param restart_patience,restart_progress Anderson's restart: this many
#'   iterations without a residual below \code{restart_progress} times the
#'   best so far.
#' @param restart_best Whether Anderson restarts from its best point every
#'   time (Phase R's rule; research), rather than where it is unless its
#'   residual has grown past twice the best.
#' @returns A list for \code{*_blp_solve()}'s \code{control}.
#' @noRd
.blp_control <- function(tol = 1e-8, max_iter = 1000L, damping = 1,
                         method = "simple", settings = list(), map = NULL,
                         patience = NULL, window = NULL, keep_pass = FALSE,
                         restart_patience = 5L, restart_progress = 0.9,
                         restart_best = FALSE) {
  tol <- .int64_to_double(tol, "`tol`")
  max_iter <- .int64_to_double(max_iter, "`max_iter`")
  if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol < 0) {
    stop("`tol` must be a single non-negative number: the largest absolute ",
         "log-share residual accepted.", call. = FALSE)
  }
  if (!is.numeric(max_iter) || length(max_iter) != 1L || !is.finite(max_iter) ||
      max_iter < 0 || max_iter != round(max_iter) ||
      max_iter > .Machine$integer.max) {
    stop("`max_iter` must be a whole number between 0 and ",
         ".Machine$integer.max: the most evaluations of the shares.", call. = FALSE)
  }
  if (!is.numeric(damping) || length(damping) != 1L || !is.finite(damping) ||
      damping <= 0 || damping > 1) {
    stop("`damping` must be a single number in (0, 1].", call. = FALSE)
  }
  iteration <- if (method %in% c("anderson", "squarem")) method else "simple"
  if (is.null(map)) map <- if (identical(method, "original")) "plain" else "reduced"
  if (is.null(patience)) {
    patience <- if (map == "reduced" && iteration == "simple") .BLP_PATIENCE else 0L
  }
  if (is.null(window)) window <- if (map == "reduced") .BLP_WINDOW else 0L
  st <- c(.BLP_SETTINGS$anderson, .BLP_SETTINGS$squarem)
  st[names(settings)] <- settings
  list(method = method, iteration = iteration, map = map,
       max_evals = as.integer(max_iter), tol = as.numeric(tol),
       damping = as.numeric(damping), patience = as.integer(patience),
       window = as.integer(window), memory = as.integer(st$memory),
       max_condition = as.numeric(st$max_condition),
       restart_patience = as.integer(restart_patience),
       restart_progress = as.numeric(restart_progress),
       restart_best = isTRUE(restart_best), scheme = as.integer(st$scheme),
       step_min = as.numeric(st$step_min), step_max = as.numeric(st$step_max),
       step_factor = as.numeric(st$step_factor), keep_pass = isTRUE(keep_pass))
}

# Method "simple"'s stalls: this many accepted steps without a new best
# residual (it diverges, cycles or climbs; a converging two-level iteration
# improves far sooner)
.BLP_PATIENCE <- 30L

# And its race, every method's on the two-level step: Berry's contraction
# joins from the start once this many accepted steps fail to halve the best
# residual, stops when this many of its own evaluations fail to halve its best
# residual while the method's is lower, and goes on whenever twice this many
# of the method's, alone, fail to halve the method's
.BLP_WINDOW <- 100L

# The methods, with their settings' defaults ("simple" and "original" take
# none): Anderson acceleration keeps the differences of its last `memory`
# iterates, the oldest dropped while their least squares' condition estimate
# exceeds `max_condition`; SQUAREM's step length follows `scheme` (1, 2 or
# 3), bounded to [step_min, step_max], the upper bound multiplied by
# `step_factor` when an accepted step reaches it and divided by it, not below
# its start, when one that reached it fails (PyBLP's defaults)
.BLP_SETTINGS <- list(
  simple = list(), original = list(),
  anderson = list(memory = 10L, max_condition = 1e8),
  squarem = list(scheme = 3L, step_min = 1, step_max = 1, step_factor = 4))

# The solver's statuses, by code (src/fixed_point.h)
.BLP_STATUS <- c("converged", "max_iter", "non_finite", "non_finite_start")

#' The share inversion's convergence record
#'
#' @param res An internal solver's result.
#' @param ctl Its control list.
#' @returns List: converged, status, iterations, evaluations, residual (the
#'   largest absolute log-share residual at the point returned; NA without an
#'   evaluation), tol and method.
#' @noRd
.blp_convergence <- function(res, ctl) {
  list(converged = identical(res$code, 0L), status = .BLP_STATUS[res$code + 1L],
       iterations = res$iterations, evaluations = res$evaluations,
       residual = res$residual, tol = ctl$tol, method = ctl$method)
}

#' Warn when the share inversion did not converge
#'
#' The warning has class \code{choicer_blp_unconverged}; its
#' \code{convergence} element is \code{.blp_convergence()}'s record.
#'
#' @param res An internal solver's result.
#' @param ctl Its control list.
#' @param what The function named in the message.
#' @returns NULL, invisibly.
#' @noRd
.blp_unconverged <- function(res, ctl, what) {
  if (identical(res$code, 0L)) return(invisible(NULL))
  conv <- .blp_convergence(res, ctl)
  fmt <- function(x) format(x, digits = 3)
  msg <- switch(conv$status,
    max_iter = if (conv$evaluations == 0L) {
      paste0(what, " made no evaluation of the shares (max_iter = 0) and ",
             "returns the start.")
    } else {
      paste0(what, " did not converge: after ", conv$evaluations,
             " evaluations of the shares (max_iter) the largest log-share ",
             "residual is ", fmt(conv$residual), ", above tol = ", fmt(conv$tol),
             ". It returns the point of smallest residual. Raise max_iter, or ",
             "check that the target shares are attainable.")
    },
    non_finite = paste0(
      what, " stopped after ", conv$evaluations, " evaluations of the ",
      "shares: at a trial point they were not finite even with its step ",
      "halved ten times. It returns the point of smallest residual (",
      fmt(conv$residual), "; tol = ", fmt(conv$tol), ")."),
    non_finite_start = paste0(
      what, " could not start: at the start a share is zero or not finite ",
      "(a constant far below the others', or utilities that are not finite). ",
      "It returns the start; start from other constants."))
  cond <- warningCondition(msg, class = "choicer_blp_unconverged", call = NULL)
  cond$convergence <- conv
  warning(cond)
  invisible(NULL)
}

#' The target shares of blp(), checked
#'
#' Numeric, a vector or a one-column matrix, one share per alternative (the
#' outside option's first), finite, positive and summing to 1 within the
#' engine's tolerance: checked in R, before any draw or evaluation of the
#' shares.
#'
#' @param target_shares The user's targets.
#' @param J Inside alternatives.
#' @param outside Whether the fit has an outside option.
#' @returns The targets as a plain numeric vector.
#' @noRd
.blp_targets <- function(target_shares, J, outside) {
  ts <- .int64_to_double(target_shares, "`target_shares`")
  n <- J + outside
  what <- if (outside) {
    sprintf("%d: the outside option's share, then one per inside alternative", n)
  } else {
    sprintf("%d, one per alternative", n)
  }
  if (!is.numeric(ts) || (!is.null(dim(ts)) && (length(dim(ts)) != 2L || ncol(ts) != 1L))) {
    stop("`target_shares` must be a numeric vector (or a one-column matrix) of ",
         what, ".", call. = FALSE)
  }
  ts <- as.numeric(ts)
  if (length(ts) != n) {
    stop(sprintf("`target_shares` has %d entries; it must have %s.", length(ts), what),
         call. = FALSE)
  }
  bad <- which(!is.finite(ts) | !(ts > 0))
  if (length(bad)) {
    stop(sprintf(paste0("`target_shares` must be finite and positive (entry %d is %s): ",
                        "the inversion matches log shares."), bad[1L], format(ts[bad[1L]])),
         call. = FALSE)
  }
  if (!(abs(sum(ts) - 1) <= sqrt(.Machine$double.eps))) {
    stop(sprintf("`target_shares` must sum to 1%s; they sum to %.15g.",
                 if (outside) ", the outside option's share included" else "", sum(ts)),
         call. = FALSE)
  }
  ts
}

#' The start of blp(): the fit's constants, or delta_init, checked
#'
#' @param object The fit.
#' @param delta_init NULL (the fitted constants; zeros for a fit without
#'   them) or the user's start.
#' @param J Inside alternatives.
#' @param outside Whether the fit has an outside option.
#' @param free_ok Whether, without an outside option, a start of J - 1 free
#'   constants is accepted (the mixed logit's solver prepends the zero).
#' @returns The start, a numeric vector.
#' @noRd
.blp_delta_init <- function(object, delta_init, J, outside, free_ok = FALSE) {
  if (is.null(delta_init)) {
    asc <- object$param_map$asc
    if (is.null(asc)) return(rep(0, J))
    co <- as.numeric(object$coefficients[asc])
    return(if (outside) co else c(0, co))
  }
  di <- .int64_to_double(delta_init, "`delta_init`")
  lens <- if (free_ok && !outside) c(J, J - 1L) else J
  if (!is.numeric(di) || !is.null(dim(di)) && NCOL(di) != 1L ||
      !(length(di) %in% lens)) {
    stop(sprintf(paste0("`delta_init` must be a numeric vector of %d constants, one ",
                        "per inside alternative%s."), J,
                 if (length(lens) > 1L) sprintf(" (or %d, without the first's zero)", J - 1L) else ""),
         call. = FALSE)
  }
  di <- as.numeric(di)
  if (!all(is.finite(di))) {
    stop("`delta_init` must be finite.", call. = FALSE)
  }
  di
}

#' The method blp() and the exported functions run, and its settings
#'
#' "simple" (the two-level step) and "original" (Berry, Levinsohn and
#' Pakes's contraction) take no settings in `control`; "anderson" takes
#' `memory` and `max_condition`, "squarem" `scheme`, `step_min`, `step_max`
#' and `step_factor` (.BLP_SETTINGS).
#'
#' @param method The method's name.
#' @param control A list of the method's settings.
#' @returns List: method, and settings (every one of the method's, the
#'   defaults for those not given).
#' @noRd
.blp_method <- function(method, control) {
  if (!is.character(method) || length(method) != 1L || is.na(method) ||
      !method %in% names(.BLP_SETTINGS)) {
    stop("`method` must be \"anderson\" (the default: Anderson acceleration of ",
         "the two-level step of ?blp_contraction), \"simple\" (that step alone), ",
         "\"original\" (Berry, Levinsohn and Pakes's contraction) or \"squarem\" ",
         "(SQUAREM on the step).", call. = FALSE)
  }
  if (!is.list(control) || is.object(control)) {
    stop("`control` must be a list of the method's settings.", call. = FALSE)
  }
  st <- .BLP_SETTINGS[[method]]
  if (length(control)) {
    nm <- names(control)
    if (!length(st)) {
      shown <- if (is.null(nm)) "unnamed entries" else paste0("`", nm, "`", collapse = ", ")
      stop(sprintf("`control` takes no settings with method \"%s\" (it has %s).",
                   method, shown), call. = FALSE)
    }
    if (is.null(nm) || anyNA(nm) || !all(nzchar(nm)) || anyDuplicated(nm)) {
      stop("`control` must name each of its settings once.", call. = FALSE)
    }
    bad <- setdiff(nm, names(st))
    if (length(bad)) {
      shown <- paste0("`", names(st), "`")
      stop(sprintf("`control` takes %s with method \"%s\"; it has %s.",
                   if (length(shown) > 2L) {
                     paste0(paste(shown[-length(shown)], collapse = ", "), " and ",
                            shown[length(shown)])
                   } else {
                     paste(shown, collapse = " and ")
                   }, method, paste0("`", bad, "`", collapse = ", ")), call. = FALSE)
    }
    st[nm] <- control
  }
  if (identical(method, "anderson")) {
    m <- .int64_to_double(st$memory, "`control$memory`")
    if (!is.numeric(m) || length(m) != 1L || !is.finite(m) || m != round(m) ||
        m < 1 || m > 100) {
      stop("`control$memory` must be a whole number from 1 to 100: the past ",
           "iterates whose differences Anderson acceleration keeps.", call. = FALSE)
    }
    mc <- .int64_to_double(st$max_condition, "`control$max_condition`")
    if (!is.numeric(mc) || length(mc) != 1L || is.na(mc) || mc < 1) {
      stop("`control$max_condition` must be a single number of at least 1: the ",
           "largest condition estimate of the least squares on the differences ",
           "kept.", call. = FALSE)
    }
    st$memory <- as.integer(m)
    st$max_condition <- as.numeric(mc)
  }
  if (identical(method, "squarem")) {
    num1 <- function(x, what) {
      x <- .int64_to_double(x, what)
      if (!is.numeric(x) || length(x) != 1L || !is.finite(x)) NA_real_ else as.numeric(x)
    }
    sc <- num1(st$scheme, "`control$scheme`")
    if (!sc %in% 1:3) {
      stop("`control$scheme` must be 1, 2 or 3: SQUAREM's step length (Varadhan ",
           "and Roland 2008).", call. = FALSE)
    }
    lo <- num1(st$step_min, "`control$step_min`")
    hi <- num1(st$step_max, "`control$step_max`")
    if (is.na(lo) || is.na(hi) || !(lo > 0) || hi < lo) {
      stop("`control$step_min` and `control$step_max` must be finite numbers with ",
           "0 < step_min <= step_max: the first bounds of SQUAREM's step length.",
           call. = FALSE)
    }
    fa <- num1(st$step_factor, "`control$step_factor`")
    if (is.na(fa) || fa < 1) {
      stop("`control$step_factor` must be a finite number of at least 1: what ",
           "SQUAREM multiplies the upper bound of its step length by when a step ",
           "reaches it.", call. = FALSE)
    }
    st$scheme <- as.integer(sc)
    st$step_min <- lo
    st$step_max <- hi
    st$step_factor <- fa
  }
  list(method = method, settings = st)
}

#' The share inversion of blp(), common to its methods
#'
#' Checks, in order, the stored data, the targets, the start and the
#' settings, all before the model's solver (and so before a mixed logit's
#' draws or any evaluation of the shares); then solves, warns when the
#' iteration did not converge, and labels the result.
#'
#' @param object The fit.
#' @param target_shares,delta_init,tol,max_iter,damping,method,control As in
#'   blp().
#' @param free_ok See .blp_delta_init().
#' @param solve function(start, target, ctl) running the model's internal
#'   solver.
#' @returns The constants: an n x 1 matrix with the inside alternatives as
#'   row names and the attribute "convergence" (.blp_convergence()).
#' @noRd
.blp_common <- function(object, target_shares, delta_init, tol, max_iter,
                        damping = 1, method = "anderson", control = list(),
                        free_ok = FALSE, solve) {
  if (is.null(object[["data"]])) {
    stop("blp() requires stored data. Refit with keep_data = TRUE.")
  }
  outside <- isTRUE(object$include_outside_option)
  J <- nrow(object$alt_mapping) - outside
  target <- .blp_targets(target_shares, J, outside)
  start <- .blp_delta_init(object, delta_init, J, outside, free_ok)
  m <- .blp_method(method, control)
  ctl <- .blp_control(tol, max_iter, damping, method = m$method,
                      settings = m$settings)
  res <- solve(start, target, ctl)
  .blp_unconverged(res, ctl, "blp()")
  out <- res$delta
  am <- object$alt_mapping
  inside <- am[["alt_int"]] > 0L
  lab <- as.character(am[[2L]][inside])[order(am[["alt_int"]][inside])]
  if (length(lab) == nrow(out)) rownames(out) <- lab
  attr(out, "convergence") <- .blp_convergence(res, ctl)
  out
}
