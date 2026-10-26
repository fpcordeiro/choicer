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
#'   outside option's included; see Details). Default 1e-8.
#' @param max_iter The most evaluations of the shares (default 1000).
#' @details Each iteration evaluates the shares once and moves every constant
#'   by its log-share residual against a base, \eqn{r_j - r_b} with
#'   \eqn{r = \log S - \log s(\delta)}, scaled by the exact diagonal of that
#'   residual's derivative, \eqn{A_{jj}/s_j + C_{jb}/s_b}, where \eqn{A_{jj}}
#'   is the diagonal of the share Jacobian \eqn{\partial s/\partial \delta} and
#'   \eqn{C_{jb} = \sum_i w_i P_{ij} P_{ib}} the share-weighted co-choice of
#'   \eqn{j} and the base. The base is the outside option (Fukasawa's 2024
#'   outside-share correction) or, without one, the reference of the
#'   alternative's group: the alternatives that share choice situations
#'   (directly or through others) form a connected component, whose constants
#'   the shares identify only up to a common level, so its smallest code's
#'   constant is held at its start. Berry's contraction takes the same steps
#'   without the scaling. A scaled step is at most 5 log units long, or as long
#'   as the unscaled one where that is longer, and a step whose shares are not
#'   finite is halved, up to ten times.
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
#'   Fukasawa, T. (2024). Fast and simple inner-loop algorithms of static /
#'   dynamic BLP estimations. arXiv:2404.04494.
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
                            include_outside_option = FALSE, tol = 1e-8,
                            max_iter = 1000L) {
  ctl <- .blp_control(tol, max_iter)
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
#' @details The step's scale is the nested logit's derivative,
#'   \eqn{A_{jj} = \sum_i w_i P_{ij}((1 - P_i(j|g))/\lambda_g + P_i(j|g) -
#'   P_{ij})} for alternative \eqn{j} in nest \eqn{g}, and, without an outside
#'   option, its exact co-choice with the reference \eqn{b}, \eqn{\sum_i w_i
#'   P_{ib}((1 - \lambda_g)/\lambda_g\, P_i(j|g)\, 1\{b \in g\} + P_{ij})}.
#'   In a nest with \eqn{\lambda_g > 1} nest-mates are complements, and two
#'   safeguards apply. An alternative whose co-choice with the reference is
#'   negative (the reference in its nest) takes at most the unscaled step:
#'   its scale nets out the complementarity and can be near 0 or negative.
#'   And the scale is at least \eqn{(1 + \lambda_g)/(2\lambda_g)} times the
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
                               tol = 1e-8, max_iter = 1000L) {
  ctl <- .blp_control(tol, max_iter, damping)
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
#' the other parameters fixed. The step's scale averages each draw's terms,
#' \eqn{A_{jj} = \sum_i w_i S^{-1} \sum_r P_{ijr}(1 - P_{ijr})} and
#' \eqn{C_{jb} = \sum_i w_i S^{-1} \sum_r P_{ijr} P_{ibr}}.
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
                                include_outside_option = FALSE, tol = 1e-8,
                                max_iter = 1000L, gen_seed = -1L,
                                gen_scramble = 1L, gen_S = 0L) {
  ctl <- .blp_control(tol, max_iter)
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
#'   \code{nl_blp_contraction()}.
#' @param method The solver: \code{"simple"} (Phase B adds others).
#' @param map The base map: \code{"reduced"} (the scaled step) or
#'   \code{"plain"} (Berry's contraction; tests).
#' @param keep_pass Whether the solver returns its last pass (tests).
#' @returns A list for \code{*_blp_solve()}'s \code{control}.
#' @noRd
.blp_control <- function(tol = 1e-8, max_iter = 1000L, damping = 1,
                         method = "simple", map = "reduced", keep_pass = FALSE) {
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
  list(method = method, map = map, max_evals = as.integer(max_iter),
       tol = as.numeric(tol), damping = as.numeric(damping),
       keep_pass = isTRUE(keep_pass))
}

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
