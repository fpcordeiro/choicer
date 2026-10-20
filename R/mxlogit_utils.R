
# Normalize a bounds vector to nloptr's full-length representation. Accepts
# NULL (returns the default), a full-length unnamed numeric (returned as-is
# after a length check), or a named partial numeric whose names must be a
# subset of `param_names`. Used by run_mxlogit() for both `lower` and `upper`.
.normalize_bound <- function(b, param_names, default, which) {
  if (is.null(b) || length(b) == 0L) {
    return(rep(default, length(param_names)))
  }
  if (!is.numeric(b)) stop("`", which, "` must be numeric or NULL.")
  if (anyNA(b)) stop("`", which, "` must not contain NA/NaN.")
  if (is.null(names(b))) {
    if (length(b) != length(param_names)) {
      stop("Unnamed `", which, "` must have length n_params (",
           length(param_names), "); got ", length(b), ".")
    }
    return(b)
  }
  dups <- unique(names(b)[duplicated(names(b))])
  if (length(dups)) {
    stop("Duplicate name(s) in `", which, "`: ",
         paste(dups, collapse = ", "), ".")
  }
  bad <- setdiff(names(b), param_names)
  if (length(bad)) {
    stop("Unknown parameter name(s) in `", which, "`: ",
         paste(bad, collapse = ", "),
         ". Valid names: ", paste(param_names, collapse = ", "), ".")
  }
  out <- rep(default, length(param_names))
  names(out) <- param_names
  out[names(b)] <- b
  unname(out)
}

# Labels of the random-coefficient covariance block, in the order of its
# Cholesky parameters: <prefix>_<i><j> over the lower triangle by rows, or
# the diagonal <prefix>_<k><k> without correlation. "L" names the Cholesky
# factor in the parameter vector; summary() shows "Sigma", the covariance
# L L' it implies.
.mxl_cov_names <- function(prefix, K_w, rc_correlation) {
  if (!rc_correlation) {
    return(paste0(prefix, "_", seq_len(K_w), seq_len(K_w)))
  }
  i <- rep(seq_len(K_w), seq_len(K_w))
  sprintf("%s_%d%d", prefix, i, sequence(seq_len(K_w)))
}

# Labels summary() prints for a mixed logit's parameters: the covariance
# block as Sigma_<i><j>, and the mean of a log-normal coefficient as
# exp(Mu_<variable>). run_mxlogit() checks them for repeats before the
# kernels have validated rc_dist, so only its first K_w entries are read.
.mxl_summary_labels <- function(param_names, param_map, K_w, rc_dist,
                                rc_correlation, rc_mean) {
  labels <- param_names
  if (!is.null(param_map$sigma)) {
    labels[param_map$sigma] <- .mxl_cov_names("Sigma", K_w, rc_correlation)
  }
  if (rc_mean && !is.null(param_map$mu)) {
    idx <- param_map$mu[rc_dist[seq_len(K_w)] %in% 1L]
    labels[idx] <- paste0("exp(", labels[idx], ")")
  }
  labels
}

# run_mxlogit()'s default start, in the columns' units: zero for beta, mu, the
# Cholesky factor's off-diagonal entries and the constants; on its diagonal
# (which holds log L_pp), L_pp = 1 / h_p for a normal random coefficient, where
# h_p is the typical step of its column (design_column_step(), reading W in
# place): the distance between its values for a column with two (a dummy, a
# two-level attribute), the mean distance from a value that fills more than
# half of the column over the other rows (a mostly-zero column), the standard
# deviation otherwise. A step of h_p in the column then moves the utility by a
# random amount with standard deviation 1 at the start. For a dummy with a
# share q of ones the standard deviation, sqrt(q (1 - q)), would understate the
# step: with one of J alternatives marked per situation, L_pp would start about
# sqrt(J) times higher, which took about 1.7 times the evaluations on a
# census-style design. Rescaling a column, W_p -> a W_p, gives the same model
# with row p of L (and mu_p) divided by a, and this start moves with it, so the
# optimizer's problem does not depend on the units: in exact arithmetic under
# scale_vars = "sd" (whose map undoes the rescaling), up to the power-of-two
# rounding of its scales under "bhhh". From L_pp = 0.5 in each column's own
# units, fits of rescaled designs stopped hundreds to thousands of
# log-likelihood units short under every scale_vars. The level errs high:
# "bhhh" fixes its scales at the start, and too little heterogeneity there
# makes them too generous for the Cholesky diagonal, whose curvature grows with
# L_pp (from 0.25 / sd(W_p) a fit reported convergence 554 log-likelihood
# units short). A log-normal coefficient keeps L_pp = 0.5, the standard
# deviation of the logarithm of its random part, which has no units; that
# random part, exp((L eta)_p), has median 1 in the column's units, so no start
# makes a log-normal column's units irrelevant. A column with no variation
# (relative to its level), or a step that is not finite, starts at L_pp = 1;
# steps are clamped to [2^-500, 2^500]. A W that is not double (a hand-built
# input_data) is read through a double copy.
.mxl_default_start <- function(n_params, param_map, W, rc_dist,
                               rc_correlation) {
  theta <- rep(0, n_params)
  K_w <- ncol(W)
  if (!K_w) return(theta)
  diag_idx <- if (rc_correlation) {
    param_map$sigma[cumsum(seq_len(K_w))]
  } else {
    param_map$sigma
  }
  if (!is.double(W)) storage.mode(W) <- "double"
  s <- design_column_step(W)
  ok <- is.finite(s) & s > 0
  log_L <- numeric(K_w)
  log_L[ok] <- -log(pmin(pmax(s[ok], 2^-500), 2^500))
  log_L[rc_dist[seq_len(K_w)] %in% 1L] <- log(0.5)
  theta[diag_idx] <- log_L
  theta
}

#' The profile likelihood's requirements on the constants
#'
#' `concentrate_out_asc = "profile_likelihood"` solves the constants'
#' first-order condition at every value of the other parameters, without
#' bounds, factoring their block of the information densely: it needs
#' constants, none of them bounded (or pinned), each with a finite
#' maximum-likelihood estimate (every case `.unbounded_constants()` flags is
#' an error, with the warning's remedies), and at most `max_dense` of them.
#' An alternative listed twice in a situation counts twice in
#' `alt_mapping$N_OBS`, which can hide one chosen wherever it is offered, so
#' the inside alternatives' counts are taken from the data instead
#' (`.unchosen_offers()`). Checked before the draws.
#'
#' @param n_asc Number of free constants.
#' @param asc_lower,asc_upper Their bounds.
#' @param asc_names Their names.
#' @param alt_mapping The preparation's `alt_mapping`.
#' @param ctl From `.concentrate_control()`.
#' @param d The prepared `input_data`.
#' @returns NULL, invisibly.
#' @noRd
.check_profile_inputs <- function(n_asc, asc_lower, asc_upper, asc_names,
                                  alt_mapping, ctl, d) {
  what <- "concentrate_out_asc = \"profile_likelihood\""
  if (n_asc == 0L) {
    stop(what, " concentrates the alternative-specific constants out of the ",
         "likelihood, so it needs them: use_asc = TRUE and more than one ",
         "alternative.", call. = FALSE)
  }
  b <- is.finite(asc_lower) | is.finite(asc_upper)
  if (any(b)) {
    stop(what, " solves the constants' first-order condition without bounds, ",
         "but `lower`/`upper` bound ",
         paste(utils::head(asc_names[b], 5L), collapse = ", "),
         if (sum(b) > 5L) paste0(" and ", sum(b) - 5L, " more"),
         ". Remove those bounds, or fit without concentrate_out_asc.",
         call. = FALSE)
  }
  am <- alt_mapping
  if (all(c("alt_int", "N_OBS", "N_CHOICES") %in% names(am))) {
    am <- as.list(am)
    code <- am$alt_int
    inside <- code >= 1L & !is.na(code)
    un <- .unchosen_offers(d, max(c(0L, code[inside])))
    # chosen wherever offered: no row in a situation that chose another
    am$N_OBS[inside] <- am$N_CHOICES[inside] + un[code[inside]]
  }
  u <- .unbounded_constants(am)
  if (!is.null(u)) {
    msg <- .unbounded_constants_text(u, alt_mapping)
    stop(errorCondition(paste(c(paste0(
      what, " needs every alternative-specific constant to have a finite ",
      "maximum-likelihood estimate, so that the constants' first-order ",
      "condition has a solution at every value of the other parameters."),
      msg[-1L]), collapse = "\n"),
      class = "choicer_unbounded_constants", call = NULL))
  }
  if (n_asc > ctl$max_dense) {
    stop(what, " factors the constants' block of the information densely, up ",
         "to concentrate_control$max_dense = ", ctl$max_dense, " constants; ",
         "this model has ", n_asc, ". A solve holds two J x J matrices (8 J^2 ",
         "bytes each) and each Newton step factors one (about J^3 / 3 ",
         "operations): raise max_dense, or fit without concentrate_out_asc.",
         call. = FALSE)
  }
  invisible(NULL)
}

#' Rows of each alternative in situations that chose another
#'
#' For each inside alternative code, the rows that carry it in the choice
#' situations whose chosen alternative is another (or the outside option): 0
#' exactly when the alternative is chosen in every situation that offers it,
#' however often a situation lists it. In blocks of about 2^22 rows.
#'
#' @param d A prepared `input_data`.
#' @param J The largest inside code.
#' @returns Numeric vector of length J.
#' @noRd
.unchosen_offers <- function(d, J) {
  M <- d$M
  n <- length(M)
  out <- numeric(J)
  if (!n || !J) return(out)
  off <- c(0, cumsum(as.numeric(M)))
  inside <- d$choice_idx > 0L
  chosen <- integer(n)
  chosen[inside] <- d$alt_idx[off[which(inside)] + d$choice_idx[inside]]
  step <- max(1, floor(2^22 / max(1, off[n + 1L] / n)))
  lo <- 1
  while (lo <= n) {
    hi <- min(n, lo + step - 1)
    code <- d$alt_idx[seq.int(off[lo] + 1, off[hi + 1L])]
    other <- code != rep.int(chosen[lo:hi], M[lo:hi])
    out <- out + tabulate(code[other], J)
    lo <- hi + 1
  }
  out
}

#' Weighted chosen counts of the free constants
#'
#' The data side of the constants' first-order condition: for each free
#' constant (`param_map$asc`'s order), the summed weights of the choice
#' situations whose chosen alternative carries it, and the same sum for the
#' reference (the outside option, or the base alternative, coded 1), whose
#' constant is normalized.
#'
#' @param d A prepared `input_data`.
#' @param n_asc Number of free constants.
#' @returns List with `obs` (n_asc) and `ref`.
#' @noRd
.chosen_counts <- function(d, n_asc) {
  M <- d$M
  first <- c(0, cumsum(as.numeric(M)))[seq_along(M)]
  inside <- d$choice_idx > 0L
  code <- integer(length(M))
  code[inside] <- d$alt_idx[first[inside] + d$choice_idx[inside]]
  # the free constant of the chosen alternative; 0 for the reference
  k <- if (isTRUE(d$include_outside_option)) code else code - 1L
  w <- as.numeric(d$weights)
  free <- k >= 1L
  obs <- numeric(n_asc)
  if (any(free)) {
    s <- rowsum(w[free], k[free])
    obs[as.integer(rownames(s))] <- s[, 1L]
  }
  list(obs = obs, ref = sum(w[!free]))
}

#' The inner solve of the mixed logit's profile likelihood
#'
#' The maximum over the free constants at fixed continuous parameters
#' theta_c, by `run_polish()` on delta: Newton steps on the constants' block
#' of the information, Levenberg-Marquardt damping where it is not positive
#' definite, and the line search that, near the maximum of a likelihood
#' summed over millions of units, also accepts a step whose gain is below
#' the objective's rounding by reading the gradient. The kernel's one pass
#' gives the objective, the full gradient and the block, so a trial point
#' the line search accepts is the next step's curvature; the pass is kept
#' (one point, the previous result released before the next call) and never
#' carried from one solve to the next. A solve stops converged when the
#' Newton decrement over the mean weight falls below `tol_decrement` (so that
#' the rule does not depend on the weights' scale), or when every constant's
#' score, and the reference's (minus their sum: the scores of all the
#' alternatives' constants sum to zero), is below `tol_residual` relative to
#' its weighted chosen count plus the mean weight; within `max_iter` steps
#' and three kernel passes per step.
#'
#' @param kernel_f Function of the full parameter vector returning the
#'   kernel's `list(objective, gradient, overflow, nonfinite, hessian)`.
#' @param c_idx,d_idx Positions of theta_c and of the free constants.
#' @param obs Weighted chosen counts of the free constants.
#' @param ref The reference's weighted chosen count.
#' @param w_bar Mean weight.
#' @param ctl From `.concentrate_control()`.
#' @returns Function(theta_c, delta0, max_iter) returning `list(converged,
#'   value, delta, gradient, iterations, passes, residual, decrement,
#'   message, error)`: the negated log-likelihood, the constants and the full
#'   gradient at the point returned (value and gradient only for a converged
#'   solve), and the error that ended the solve, if one did.
#' @noRd
.mxl_profile_solver <- function(kernel_f, c_idx, d_idx, obs, ref, w_bar, ctl) {
  n <- length(c_idx) + length(d_idx)
  force(kernel_f)
  residual <- function(g) {
    max(abs(g) / (obs + w_bar), abs(sum(g)) / (ref + w_bar))
  }
  function(theta_c, delta0, max_iter) {
    full <- numeric(n)
    full[c_idx] <- theta_c
    stash <- NULL
    passes <- 0L
    budget <- 3 * max_iter + 1
    spent <- FALSE
    at <- function(dl, within_budget = TRUE) {
      if (!is.null(stash) && identical(dl, stash$delta, num.eq = FALSE)) {
        return(stash$r)
      }
      stash <<- NULL
      if (within_budget && passes >= budget) {
        spent <<- TRUE
        stop("the budget is spent")
      }
      passes <<- passes + 1L
      th <- full
      th[d_idx] <- dl
      r <- kernel_f(th)
      stash <<- list(delta = dl, r = r)
      r
    }
    eval_f <- function(dl) {
      r <- at(dl)
      list(objective = r$objective, gradient = r$gradient[d_idx],
           overflow = isTRUE(r$overflow) || isTRUE(r$nonfinite))
    }
    curv_f <- function(dl) at(dl)$hessian
    stop_f <- function(dl, g) {
      if (residual(g) < ctl$tol_residual) {
        "every constant's score is below tol_residual"
      }
    }
    pol <- run_polish(delta0, eval_f, curv_f,
                      control = .polish_control(list(
                        max_iter = max_iter,
                        tol_decrement = ctl$tol_decrement * w_bar)),
                      method = "newton", stop_f = stop_f)
    # the budget's stop is a failure to converge, not an error
    out <- list(converged = FALSE, value = NA_real_, delta = pol$par,
                gradient = NULL, iterations = pol$iterations, passes = passes,
                residual = NA_real_, decrement = pol$decrement / w_bar,
                message = if (spent) {
                  paste0("the solve used its ", budget, " passes")
                } else pol$message,
                error = if (!spent) pol$error)
    if (!isTRUE(pol$converged) || !is.null(pol$error) || !is.finite(pol$value)) {
      return(out)
    }
    r <- at(pol$par, within_budget = FALSE)
    out$converged <- !isTRUE(r$overflow) && !isTRUE(r$nonfinite)
    out$value <- r$objective
    out$gradient <- r$gradient
    out$passes <- passes
    out$residual <- residual(r$gradient[d_idx])
    out
  }
}

#' Runs mixed logit estimation
#'
#' Estimates a mixed logit model via simulated maximum likelihood.
#'
#' Two workflows are supported:
#' \describe{
#'   \item{Convenience}{Supply \code{data} and column names. Data preparation
#'     (\code{\link{prepare_mxl_data}}) and Halton draw generation
#'     (\code{\link{get_halton_normals}}) are handled automatically.}
#'   \item{Advanced}{Call \code{\link{prepare_mxl_data}} and
#'     \code{\link{get_halton_normals}} yourself, then pass the results via
#'     \code{input_data} and \code{eta_draws}.}
#' }
#'
#' \strong{Cross-section vs panel.} With \code{person_col}, each decision
#' maker \eqn{n} draws one taste vector \eqn{\beta_n} from the mixing
#' distribution \eqn{f(\beta \mid \theta)} and keeps it in all of their choice
#' situations \eqn{t = 1, \dots, T_n}. The likelihood integrates the product
#' of their conditional logit probabilities (Revelt and Train 1998),
#' \deqn{L_n(\theta) = \int \prod_{t=1}^{T_n} P_{nt}(j_{nt} \mid \beta)
#'   f(\beta \mid \theta) \, d\beta,}
#' simulated with \code{S} draws per decision maker, so repeated choices
#' speak to each decision maker's tastes, and
#' \code{\link{conditional_tastes}} recovers what they reveal. Without
#' \code{person_col} (the default), each choice situation is an independent
#' draw from the mixing distribution: \eqn{L_i(\theta) = \int P_i(j_i \mid
#' \beta) f(\beta \mid \theta) \, d\beta}. That is the right model when tastes
#' are redrawn in every situation (or each decision maker chooses once). On
#' panel data with stable tastes it is a composite (pseudo-) likelihood: it
#' still targets the population mixing distribution, but it ignores the
#' within-person correlation of choices, so it is less efficient, its
#' standard errors must be clustered by decision maker (\code{cluster_col}),
#' and its log-likelihood, AIC and BIC are not comparable with those of the
#' panel fit.
#'
#' @param data Data frame containing choice data (convenience workflow).
#'   Mutually exclusive with \code{input_data}.
#' @param id_col Name of the column identifying choice situations.
#' @param alt_col Name of the column identifying alternatives.
#' @param choice_col Name of the column indicating chosen alternative (1/0).
#' @param covariate_cols Vector of column names for fixed covariates. None
#'   may take a name the model gives its parameters, or prints for them in
#'   \code{summary()}: \code{L_<i><j>} and \code{Sigma_<i><j>}, with
#'   \code{use_asc = TRUE} \code{ASC_<label>}, and with
#'   \code{rc_mean = TRUE} \code{Mu_<variable>} and, for a log-normal
#'   coefficient, \code{exp(Mu_<variable>)}.
#' @param random_var_cols Vector of column names for random coefficients.
#' @param input_data List output from \code{\link{prepare_mxl_data}} (advanced
#'   workflow). Mutually exclusive with \code{data}.
#' @param eta_draws Array of shape K_w x S x U with standard normal draws, one
#'   K_w x S block per likelihood unit: U is the number of decision makers
#'   when \code{input_data} was prepared with \code{person_col}
#'   (\code{length(input_data$Ti)}), and the number of choice situations
#'   otherwise. Required for the advanced workflow; auto-generated from
#'   \code{S} in the convenience workflow. Post-hoc methods
#'   (\code{vcov(fit, type = )}, \code{\link{wesml_vcov}},
#'   \code{\link{conditional_tastes}}, prediction)
#'   regenerate the Halton draws \code{\link{get_halton_normals}} builds from
#'   the recorded draw count, so build \code{eta_draws} with it for them to
#'   reproduce the estimation draws.
#' @param S Integer number of Halton draws per decision maker, or per choice
#'   situation without \code{person_col} (convenience workflow only). Default
#'   100.
#' @param rc_dist Integer vector indicating distribution of random coefficients
#'   (0 = normal, 1 = log-normal). Default: all normal.
#' @param rc_mean Logical indicating whether to estimate means for random
#'   coefficients.
#' @param rc_correlation Logical indicating whether random coefficients are
#'   correlated (convenience workflow). Ignored when \code{input_data} is used
#'   (taken from the prepared data).
#' @param use_asc Logical indicating whether to include alternative-specific
#'   constants: one per inside alternative with the outside option, otherwise
#'   one per alternative but the first (\code{outside_opt_label} when it is
#'   given, else the first label in \code{sort()} order, which for a factor
#'   is the first of its levels present), whose constant is normalized to 0.
#'   With \code{FALSE} the utilities have no constants, so alternatives
#'   differ only through their covariates. An alternative's constant has no
#'   finite maximum-likelihood estimate when the alternative is never chosen
#'   in the data, or is chosen in every choice situation that offers it (the
#'   likelihood keeps rising as the constant falls, or grows); when the
#'   outside option or the base alternative is never chosen, or chosen
#'   wherever it is offered, the other constants have none together. The fit
#'   then warns (class \code{choicer_unbounded_constants}) before the
#'   optimizer runs, naming them; a finite bound in the direction a constant
#'   diverges (\code{lower} for a never-chosen alternative's, \code{upper}
#'   for an always-chosen one's) holds it, and that alternative is not named
#'   (the outside option's and the base's cases, which move every constant
#'   together, are named whatever the bounds).
#' @param theta_init Initial parameter vector in natural-scale units. If
#'   \code{NULL} (default), the start is in the columns' units: zero for the
#'   \eqn{\beta}, \eqn{\mu} and ASC blocks and the Cholesky factor's
#'   off-diagonal entries; on its diagonal, \eqn{L_{pp} = 1 / h_p} for a
#'   normal random coefficient, where \eqn{h_p} is the typical step of its
#'   column of \eqn{W}: the distance between its values for a column with two
#'   (a dummy: 1), the mean distance from a value that fills more than half of
#'   the column over the other rows (a mostly-zero column), and the standard
#'   deviation otherwise; a column with no variation starts at
#'   \eqn{L_{pp} = 1}. At the start, a step of \eqn{h_p} in the column then
#'   moves the utility by a random amount with standard deviation 1, whatever
#'   the column's units of measurement. A log-normal coefficient starts at
#'   \eqn{L_{pp} = 0.5}, the standard deviation of the logarithm of its random
#'   part, which has no units (that random part has median 1 in the column's
#'   units, so a log-normal covariate belongs in units where its coefficient is
#'   of order one). Where the default start falls outside
#'   \code{lower}/\code{upper}, it is moved to the nearest bound, with a
#'   message. A supplied \code{theta_init} is used as given.
#' @param lower,upper Optional parameter bounds for the optimizer, in
#'   natural-scale units (mapped to the optimizer's coordinates when
#'   \code{scale_vars != "none"}); \code{lower} must not exceed \code{upper},
#'   nor be \code{Inf} (nor \code{upper} \code{-Inf}). Each accepts three
#'   forms:
#'   \describe{
#'     \item{\code{NULL}}{(default) Unbounded (\code{-Inf}/\code{Inf}).}
#'     \item{Unnamed numeric vector of length \code{n_params}}{Full-length
#'       vector ordered exactly like \code{theta_init} (the nloptr-native form).}
#'     \item{Named numeric vector}{Names must be a subset of the parameter
#'       names (\eqn{\beta} block: column names of \code{X};
#'       \eqn{\mu} block: \code{Mu_<col>} (if \code{rc_mean = TRUE});
#'       Cholesky block: \code{L_<i><j>} for \eqn{i \ge j}; ASC block:
#'       \code{ASC_<level>}). Unlisted parameters default to \eqn{\pm\infty}.
#'       This is the recommended form for typical use, e.g.
#'       \code{lower = c(L_11 = -5, L_22 = -5)} to clip Cholesky diagonals.}
#'   }
#' @param optimizer Optimizer to use: \code{"nloptr"} (default), \code{"optim"},
#'   or a custom function. See \code{\link{run_mnlogit}} for details.
#' @param control List of optimizer-specific control parameters. With
#'   \code{scale_vars} other than \code{"none"}, settings in parameter units
#'   apply to the scaled coordinates.
#' @param se_method Method for computing standard errors. One of
#'   \code{"hessian"} (default) for the analytical Hessian of the simulated
#'   log-likelihood, \code{"bhhh"} for the BHHH/outer-product-of-gradients
#'   (OPG) estimator, \code{"sandwich"} for the robust (Huber-White)
#'   variance \eqn{V = A^{-1} B A^{-1}} (bread \eqn{A} = weighted negated
#'   Hessian, meat \eqn{B} = weight-squared OPG), or \code{"cluster"} for the
#'   cluster-robust sandwich (requires \code{cluster_col} or a prepared
#'   \code{input_data} with a \code{cluster} field). Use \code{"sandwich"} for
#'   valid inference under choice-based / WESML weighting, where the
#'   inverse-Hessian and ordinary BHHH are invalid; it reduces to the usual
#'   robust variance under uniform weights. BHHH needs no second derivatives
#'   and so costs less than the Hessian. For both, the work per decision maker
#'   involves only the parameters that enter that decision maker's likelihood
#'   (the fixed coefficients, the random coefficients' means and Cholesky
#'   factor, and the constants of the alternatives they faced), so it does not
#'   grow with the number of alternatives in the data; the
#'   parameters-by-parameters matrix and its inversion do. BHHH may
#'   underestimate standard errors in finite samples or away from the
#'   optimum. Any of these variances can also be recomputed post hoc via
#'   \code{vcov(fit, type = )}. Without \code{person_col}, clustering repairs
#'   the inference, not the likelihood: the cross-sectional likelihood still
#'   treats each choice situation as an
#'   independent draw from the mixing distribution, and on panel data is a
#'   less efficient composite likelihood for the same taste distribution.
#'   With \code{person_col} the likelihood unit is the decision maker, so
#'   scores are per decision maker and
#'   \code{"sandwich"} is already robust to dependence across a decision
#'   maker's choice situations; \code{"cluster"} is then needed only for
#'   coarser groups (e.g. households or markets) that nest decision makers.
#' @param cluster_col Optional name of a column in \code{data} holding cluster
#'   labels for cluster-robust standard errors (e.g. a person id when the same
#'   decision maker contributes several choice situations to a
#'   cross-sectional fit). Must be constant within each \code{id_col}, and
#'   within each decision maker when \code{person_col} is used (clusters must
#'   nest decision makers). Supplying \code{cluster_col} without an explicit
#'   \code{se_method} selects \code{se_method = "cluster"}.
#' @param scale_vars How the optimizer's coordinates are scaled. The choice
#'   does not change the estimator; it changes how fast, and whether, the
#'   optimizer reaches the maximum. \code{"bhhh"} (default) scales each
#'   parameter by the inverse square root of the diagonal of the BHHH (outer
#'   product of gradients) matrix at the start values, rounded to a power of
#'   two, from one extra gradient pass, so that the BHHH diagonal, an estimate
#'   of the curvature, is near one at the start values in every coordinate: it
#'   reaches the constants and the Cholesky factor, whose curvatures can differ
#'   by orders of magnitude that no column scale reaches (the constants of
#'   large and small alternatives, say). A start near the estimates
#'   (\code{theta_init}) gives it a better diagonal. Where the diagonal is
#'   unusable (no finite, positive entry, as when the likelihood overflows at
#'   the start values), the parameters are left unscaled, with a message.
#'   \code{"none"} hands the optimizer the parameters as they are (the default
#'   in choicer 0.2.x). With \code{"sd"} (sample standard deviation),
#'   \code{"mad"} (\code{stats::mad}, i.e. 1.4826 \eqn{\times} median absolute
#'   deviation; SD-equivalent under normality) or \code{"iqr"}
#'   (\code{stats::IQR(x) / 1.349}; also SD-equivalent under normality), the
#'   optimizer works on the parameters the design would have with every column
#'   of \code{X} and \code{W} divided by the chosen scale; the constants keep
#'   their scale, the data are not divided, and the fit holds no scaled copy of
#'   them. Robust scales (\code{"mad"}/\code{"iqr"}) better capture the bulk
#'   for heavy-tailed columns where SD is dominated by outliers, but
#'   \code{stats::mad} can return zero when more than half of a column's
#'   entries are identical (e.g., a sparse 0/1 dummy) and will then trigger the
#'   same near-constant-column error as \code{"sd"}. The parameters of
#'   log-normal random coefficients (\code{rc_dist == 1}) keep their scale
#'   under these three, since the shifted log-normal parameterization does not
#'   admit a closed-form rescaling under multiplicative column scaling. With
#'   any choice but \code{"none"} (so by default), the optimizer, a custom one
#'   included, works in the scaled coordinates (\code{theta_init}, \code{lower}
#'   and \code{upper} are given in natural units and mapped), so settings in
#'   parameter units in \code{control} (nloptr's \code{xtol_abs}; optim's
#'   \code{parscale}, and \code{pgtol}, which then applies to the scaled
#'   gradient) apply to the scaled coordinates; \code{"none"} keeps the
#'   parameters' units. Coefficients are reported in natural units and the
#'   standard errors computed in them, from the information matrix
#'   equilibrated before it is inverted, whatever the choice, so, where the
#'   likelihood is well identified and every choice reaches its maximum,
#'   reported quantities do not depend on this choice beyond the optimizer's
#'   tolerance.
#' @param weights Optional weight vector (convenience workflow), one weight per
#'   choice situation in ascending-id order (see
#'   \code{\link{prepare_mxl_data}}). With \code{person_col}, situations are
#'   reordered by decision maker, so a positional vector is accepted only when
#'   that leaves the id order unchanged; otherwise supply
#'   \code{weights_col}. If \code{NULL},
#'   equal weights are used. All weights must be finite and strictly positive.
#'   With \code{person_col} they are decision-maker weights, constant within
#'   each decision maker: the objective is \eqn{\sum_n w_n \log L_n}.
#' @param weights_col Optional name of a column in \code{data} holding a per-row
#'   weight (constant within each choice situation, finite and strictly positive).
#'   Mutually exclusive with
#'   \code{weights}; the recommended way to pass WESML weights from
#'   \code{\link{sample_by_choice}} / \code{\link{wesml_weights}}, since
#'   alignment is by id rather than by position. Convenience workflow only. If
#'   \code{data} carries choice-based-sampling provenance (a
#'   \code{"choice_sampling"} attribute, as attached by
#'   \code{\link{sample_by_choice}} / \code{\link{wesml_weights}}) and neither
#'   \code{weights} nor \code{weights_col} is supplied, the recorded weight
#'   column is auto-detected and applied (with a message); if that column is
#'   absent the call errors rather than silently fitting an unweighted model
#'   under a WESML label. With \code{person_col} the weight must also be
#'   constant within each decision maker, and WESML is not supported:
#'   choice-based weights vary with the chosen alternative across a decision
#'   maker's situations, so they cannot weight the panel likelihood.
#' @param outside_opt_label (Convenience workflow.) The label of the outside
#'   option or of the base alternative, given in the class of the alternative
#'   column. With \code{include_outside_option = TRUE} it names the outside
#'   option's row of \code{alt_mapping}. Without the outside option it names the
#'   base alternative, coded first (with constants, its constant is normalized
#'   to 0); it must be one of the alternatives, and \code{NULL} (the default)
#'   makes the first alternative in sort order the base.
#' @param include_outside_option Logical whether to include an outside option
#'   (convenience workflow).
#' @param draws Draw storage mode. One of \code{"store"} (default) or \code{"generate"}.
#'   \code{"store"} pre-materializes the full \eqn{K_w \times S \times U} Halton cube
#'   (\eqn{8 K_w S U} bytes), one block per likelihood unit (\code{U} decision
#'   makers with \code{person_col}, choice situations otherwise; existing
#'   behavior, exact reproducibility). It supports at most \eqn{2^{31} - 1}
#'   points (\eqn{S \times U}). \code{predict()}, \code{logsum()},
#'   \code{consumer_surplus()}, \code{elasticities()},
#'   \code{diversion_ratios()}, \code{blp()} and \code{gof()} (hence
#'   \code{summary()}'s fit statistics) integrate over one block per choice
#'   situation instead (see \code{\link{get_halton_normals}}). With up to 128
#'   random coefficients they form each block on the fly, in parallel, from
#'   the points \code{get_halton_normals(S, N, K_w)} would hold for the
#'   \eqn{N} choice situations, mapped by R's \code{qnorm()}, so they never
#'   hold that cube and have no limit on the number of points; with more,
#'   they build the cube. Their draws, and so their results, are the cube's,
#'   bit for bit where randtoolbox and choicer are compiled with the same
#'   floating-point contraction, as on CRAN's platforms (otherwise some of the
#'   underlying uniforms differ in their last bit, which moves those draws by
#'   at most a few parts in \eqn{10^9}, in the far tails); sums over choice
#'   situations are bit for bit at one thread, and at more threads vary in
#'   their last bits with the threads' order, as they always did.
#'   \code{blp()} may keep the draws it forms across its iterations, within
#'   its \code{keep_draws_bytes} (4 GiB by default), with either storage
#'   mode. In the convenience workflow, \code{run_mxlogit()} warns before it
#'   builds an estimation cube above 1 GiB: the cube is held while the model is
#'   fitted, and \code{vcov(type = )}, \code{\link{wesml_vcov}} and
#'   \code{\link{conditional_tastes}} rebuild it whole. \code{"generate"}
#'   computes each unit's draws on-the-fly in C++ from a stored seed, eliminating the O(U)
#'   cube; recommended for memory-constrained or large-N settings. With
#'   \code{scramble = "permuted"}, each base-\eqn{b} digit position in each
#'   dimension receives a seeded permutation shared across sequence indices.
#'   This is not Owen's nested-uniform scramble and does not carry standard
#'   randomized-QMC unbiasedness or error-estimation guarantees. Only supported
#'   in the convenience workflow.
#' @param seed Integer master seed for the on-the-fly generator. Used only when
#'   \code{draws = "generate"}. If \code{NULL} (default), a seed is drawn from R's
#'   RNG at call time (so \code{set.seed()} governs reproducibility). Ignored when
#'   \code{draws = "store"}.
#' @param scramble Scrambling mode for on-the-fly Halton draws. One of
#'   \code{"permuted"} (default) for seeded position-wise digit permutations or
#'   \code{"none"} for plain Halton (identity permutations). The historical value
#'   \code{"owen"} is accepted with a deprecation warning as an alias for
#'   \code{"permuted"}; the implementation is not Owen's nested-uniform scramble.
#'   \code{"none"} uses the points of the randtoolbox sequence, mapped to
#'   normals by choicer's own inverse normal CDF, which is within about
#'   \eqn{10^{-14}} of \code{qnorm()}: its draws agree with store mode's to
#'   that precision, not bit for bit. Simulation-draw
#'   sensitivity should be assessed by increasing \code{S} and, for
#'   \code{"permuted"}, varying \code{seed}. Used only when
#'   \code{draws = "generate"}.
#' @param keep_data Logical. If \code{TRUE} (default), stores prepared data in
#'   the returned object for post-estimation functions.
#' @param nloptr_opts Deprecated. Use \code{optimizer} and \code{control}
#'   instead.
#' @param person_col Optional name of the column in \code{data} identifying
#'   decision makers (respondents). When supplied, all choice situations of a
#'   decision maker share one draw of the random coefficients (panel
#'   likelihood); \code{id_col} must still identify choice situations uniquely
#'   across the data set. \code{NULL} (default) makes each choice situation its
#'   own decision maker (cross-sectional likelihood). Convenience workflow
#'   only; in the advanced workflow pass \code{person_col} to
#'   \code{\link{prepare_mxl_data}}.
#' @param polish A Newton phase after the optimizer (see the section
#'   below): \code{"none"} (default) has none; \code{"newton"} steps on the
#'   analytic Hessian of the simulated log-likelihood, \code{"bhhh"} on the
#'   BHHH (outer product of the gradients) matrix, both with the fit's draws
#'   and weights.
#' @param polish_control List of settings of the Newton phase, any of
#' \code{max_iter} (20; 0 evaluates the decrement at the optimizer's solution
#' without stepping), \code{tol_decrement} (\eqn{10^{-8}}{1e-8}, in
#' log-likelihood units) and \code{tol_step} (\eqn{10^{-12}}{1e-12}: a step
#' shorter than this, relative to the parameters' magnitudes plus their
#' scales in the optimizer's coordinates, stops the polish unconverged).
#' @param concentrate_out_asc How the alternative-specific constants are
#'   estimated (see the section below). \code{"none"} (default): with the
#'   other parameters, by the optimizer. \code{"profile_likelihood"}: the
#'   optimizer moves the other parameters only, and at each of its points the
#'   constants are solved from their first-order condition; the estimator, and
#'   so the estimates and their standard errors, are the same. \code{"blp"} is
#'   reserved for matching aggregate shares instead, a different estimator,
#'   and is not implemented.
#' @param concentrate_control List of settings of the constants' inner
#'   solve, any of \code{max_iter} (25 Newton steps per solve; the first, from
#'   the start values, may take four times as many), \code{tol_residual}
#'   (\eqn{10^{-10}}{1e-10}: the solve stops when every constant's score, and
#'   the reference's, is below it relative to its alternative's weighted
#'   chosen count plus the mean weight), \code{tol_decrement}
#'   (\eqn{10^{-12}}{1e-12}: or when the Newton decrement over the mean
#'   weight falls below it, in log-likelihood units of weights with mean one,
#'   so that neither rule depends on the weights' scale) and
#'   \code{max_dense} (10000: the most constants whose block of the
#'   information is formed and factored densely; for \eqn{J}{J} constants a
#'   solve holds two \eqn{J \times J}{J x J} matrices, \eqn{8 J^2}{8 J^2} bytes
#'   each (0.8 GB at \eqn{J = 10^4}{J = 1e4}), besides the kernel's
#'   accumulators, and each Newton step factors one, about
#'   \eqn{J^3 / 3}{J^3 / 3} operations).
#' @section Newton phase (polish):
#' The polish does not change the estimator. From wherever the optimizer
#' stopped (at its tolerance or at \code{maxeval}), it finishes the climb to
#' the maximum the optimizer approached, which first-order methods such as
#' L-BFGS reach slowly when curvatures differ by orders of magnitude
#' (thousands of constants, say), and certifies it. It converges when the
#' Newton decrement
#' \eqn{\nu = \frac{1}{2} g^\top C^{-1} g}{nu = g' C^-1 g / 2}, with \eqn{g}
#' the gradient and \eqn{C} the curvature on the parameters it moves, falls
#' below \code{tol_decrement}. On a positive definite Hessian,
#' \eqn{\sqrt{2\nu}}{sqrt(2 nu)} is the distance from the estimates to the
#' maximum of the log-likelihood's quadratic model in Hessian standard
#' errors: at the default \eqn{10^{-8}}{1e-8}, every linear combination of
#' the estimates lies within \eqn{1.4 \times 10^{-4}}{1.4e-4} of its standard
#' error of that maximum. BHHH is positive semi-definite everywhere and
#' cheaper per step, but converges linearly; where the Hessian \eqn{H} is
#' positive definite, its \eqn{\sqrt{2\nu}}{sqrt(2 nu)} is that distance in
#' the standard errors of the sandwich \eqn{H^{-1} B H^{-1}}{H^-1 B H^-1}
#' (the robust standard errors of an unweighted fit), but on its own it shows
#' a stationary point, not a maximum. \code{"bhhh"} needs at least as many
#' likelihood units as the parameters it moves.
#'
#' Where the curvature is not positive definite the steps are damped
#' (Levenberg-Marquardt); they respect \code{lower}/\code{upper} and pass a
#' line search. The constants of inside alternatives without a finite
#' estimate (see \code{use_asc}) are held at the optimizer's values, and the
#' certificate concerns the other parameters given them. Where the base
#' alternative or the outside option is the one without, the other constants
#' diverge together, and the polish does not report convergence, unless the
#' reference is never chosen and the constant of an alternative that is
#' chosen is held or bounded above, which stops them. A random coefficient
#' whose variance's maximum is at zero drives its log-Cholesky diagonal
#' toward \eqn{-\infty} by Newton steps of \eqn{\frac{1}{2}}{1/2} to 1, each
#' shrinking the decrement by a factor of about \eqn{e}: \code{"newton"} may
#' certify a point where the variance is negligible but that diagonal, and
#' its standard error, depend on \code{tol_decrement}; \code{"bhhh"}, whose
#' curvature in that direction vanishes faster, comes to hold the diagonal
#' and stops unconverged. A parameter without usable curvature (no
#' information at the point, or a curvature that has all but vanished, giving
#' a step more than 100 times its scale in the optimizer's coordinates) is
#' held for the step and stops the polish unconverged. Under
#' \code{scale_vars = "none"} those scales are the parameters' natural units,
#' so a badly scaled column far from its maximum can be held as if its
#' curvature had vanished; the default \code{"bhhh"} scaling avoids this. The
#' standard errors are computed at the polished estimates, by
#' \code{se_method}. The polish evaluates its curvature at the point it
#' returns before it stops (unless it returns the start, or the curvature
#' itself failed), so that matrix serves the variance whenever it is the one
#' \code{se_method} needs, instead of being evaluated again: after
#' \code{"newton"}, the Hessian of \code{"hessian"} and the bread of
#' \code{"sandwich"} and \code{"cluster"}; after \code{"bhhh"}, the BHHH
#' matrix of \code{"bhhh"} (the sandwich's meat squares the weights, so it is
#' computed). Where it serves so, a polish that certifies the optimizer's
#' solution at once costs one evaluation of the log-likelihood and its
#' gradient and one factorization of the curvature more than the fit without
#' it.
#'
#' Each step costs a curvature evaluation and a line search. A smaller
#' \code{maxeval} (nloptr) or \code{maxit} (optim) in \code{control} hands
#' the polish a point sooner, which saves time where the optimizer's last
#' iterations creep; from farther away the steps are shortened (and damped
#' where the curvature is not positive definite), and more of them are
#' needed.
#' @section Profile likelihood of the constants (concentrate_out_asc):
#' With many alternatives the constants \eqn{\delta}{delta} are most of the
#' parameters. \code{concentrate_out_asc = "profile_likelihood"} maximizes the
#' profile log-likelihood
#' \eqn{\ell^*(\theta_c) = \max_\delta \ell(\theta_c, \delta)}{l*(theta_c) =
#' max_delta l(theta_c, delta)} over the other parameters \eqn{\theta_c}{theta_c},
#' so the optimizer works on those alone. The estimator does not change: the
#' profile's maximizer is the joint one, and the variance is the joint fit's,
#' computed at \eqn{(\hat\theta_c, \hat\delta)}{(theta_c_hat, delta_hat)}
#' from the full matrices. The full Hessian's inverse gives
#' \eqn{\theta_c}{theta_c} the inverse of the Schur complement
#' \eqn{H_{cc} - H_{c\delta} H_{\delta\delta}^{-1} H_{\delta c}}{H_cc - H_cd
#' H_dd^-1 H_dc}, the curvature of the profile; the inverse of
#' \eqn{H_{cc}}{H_cc} alone would treat the constants as known and
#' understate the standard errors. At each point the optimizer tries, the
#' constants solve their first-order condition: each alternative's weighted
#' chosen count equals its predicted count, the probabilities weighted by
#' each decision maker's posterior over the draws given their choices (not
#' the unconditional shares that BLP share matching equates). The solve takes
#' damped Newton steps on the constants' block of the information,
#' warm-started from the best point the optimizer has reached, and the
#' optimizer's gradient is the joint gradient at the solution (the envelope
#' theorem). A solve that does not converge is an overflow to the optimizer,
#' which then backtracks, and is counted; the fit warns (class
#' \code{choicer_profile_unconverged}) when such solves come after the
#' optimizer's best point, which may have stopped it, or when the constants
#' cannot be solved at its solution.
#'
#' Every constant needs a finite maximum-likelihood estimate: with an
#' alternative never chosen, or chosen wherever it is offered, or such a
#' reference (see \code{use_asc}), the fit stops before the draws (class
#' \code{choicer_unbounded_constants}), and so it does with bounds on the
#' constants (bounds on the other parameters go to the optimizer). The counts
#' do not reveal every such pattern (a group of alternatives chosen only
#' where nothing outside the group is offered, say): the constants then drift
#' as in a joint fit, and their standard errors are huge or missing. Each of
#' the optimizer's evaluations costs an inner solve, a few passes over the
#' data and a dense factorization of the constants' block, so profiling pays
#' only where the joint optimizer needs many more evaluations; with
#' \code{polish = "newton"} the joint maximum is certified after it.
#' @returns A \code{choicer_mxl} object (inherits from \code{choicer_fit}).
#'   Standard S3 methods available: \code{summary()}, \code{coef()},
#'   \code{vcov()}, \code{logLik()}, \code{AIC()}, \code{BIC()},
#'   \code{nobs()}. A panel fit (\code{person_col}) also carries
#'   \code{n_persons}, the number of decision makers; \code{nobs()} remains
#'   the number of choice situations, which is also the \eqn{n} in BIC's
#'   \eqn{\log n} penalty. The variance is computed once the optimizer has
#'   stopped. If that step fails with an error (the kernels' memory check,
#'   see \code{\link{set_num_threads}}, say), the fit returns its estimates
#'   without it (\code{vcov} \code{NULL}, \code{NA} standard errors, the
#'   error's message in \code{variance_error}) and warns, naming the
#'   \code{vcov(fit, type = )} that computes the variance later from a fit
#'   that keeps its data (\code{keep_data = TRUE}); the memory check's message
#'   names the threads that fit or how to lift its limit. With \code{polish},
#'   \code{fit$optimizer$polish} records the Newton phase: \code{method},
#'   \code{converged}, \code{iterations} (steps taken) and \code{damped},
#'   \code{decrement} (\code{NA} where the curvature on the parameters moved
#'   was not positive definite) and \code{lambda} (its damping) at the last
#'   curvature evaluated, \code{message} (why it stopped),
#'   \code{evaluations} (of the objective and the curvature), \code{active},
#'   \code{held} and \code{flat} (the parameters on a bound, held, and left
#'   out of the last steps for want of usable curvature), \code{improvement}
#'   (the log-likelihood gained), \code{error}, \code{elapsed_time} (the
#'   polish's seconds, which \code{fit$optimizer$elapsed_time} includes), and
#'   \code{optimizer}, the optimizer's own \code{convergence},
#'   \code{message}, \code{iterations} and \code{elapsed_time}. A converged
#'   polish sets \code{convergence} to the optimizer's success code (1 for
#'   nloptr, 0 for optim and a custom optimizer) and leads \code{message}
#'   with its verdict; otherwise both keep the optimizer's, the polish's
#'   reason appended to the message. A polish that stops with an error (the
#'   Hessian's memory check, say) keeps its last accepted point and warns
#'   (class \code{choicer_polish_failed}). With \code{concentrate_out_asc},
#'   \code{fit$optimizer} describes the optimizer over the other parameters,
#'   and \code{fit$optimizer$concentration} records the inner solves:
#'   \code{method}, \code{evaluations} (the optimizer's), \code{solves},
#'   \code{failed} (solves at the optimizer's points that did not converge),
#'   \code{inner_iterations} (Newton steps), \code{passes} (of the
#'   constants' kernel), \code{max_residual} (the largest relative score of
#'   the constants at a converged solve, above \code{tol_residual} where the
#'   decrement stopped the solve), \code{residual} and \code{decrement} (over
#'   the mean weight; NA where the residual rule stopped the solve) at the
#'   profile's solution, before any polish, and \code{elapsed_time}.
#' @references Revelt, D. and Train, K. (1998). Mixed logit with repeated
#'   choices: households' choices of appliance efficiency level.
#'   \emph{Review of Economics and Statistics} 80(4), 647-657.
#' @seealso \code{\link{conditional_tastes}}, \code{\link{run_hmnlogit}} (the
#'   hierarchical Bayes counterpart)
#' @inheritSection prepare_mnl_data Column names
#' @examples
#' \donttest{
#' library(data.table)
#' set.seed(42)
#' N <- 100; J <- 3
#' dt <- data.table(id = rep(1:N, each = J), alt = rep(1:J, N))
#' dt[, `:=`(x1 = rnorm(.N), w1 = rnorm(.N), w2 = rnorm(.N))]
#' dt[, choice := 0L]
#' dt[, choice := sample(c(1L, rep(0L, J - 1))), by = id]
#'
#' fit <- run_mxlogit(
#'   data = dt, id_col = "id", alt_col = "alt", choice_col = "choice",
#'   covariate_cols = "x1", random_var_cols = c("w1", "w2"), S = 50L
#' )
#' summary(fit)
#' }
#' @importFrom nloptr nloptr
#' @export
run_mxlogit <- function(
    data = NULL,
    id_col = NULL,
    alt_col = NULL,
    choice_col = NULL,
    covariate_cols = NULL,
    random_var_cols = NULL,
    input_data = NULL,
    eta_draws = NULL,
    S = 100L,
    rc_dist = NULL,
    rc_mean = FALSE,
    rc_correlation = FALSE,
    use_asc = TRUE,
    theta_init = NULL,
    lower = NULL,
    upper = NULL,
    optimizer = NULL,
    control = list(),
    se_method = c("hessian", "bhhh", "sandwich", "cluster"),
    scale_vars = c("bhhh", "none", "sd", "mad", "iqr"),
    weights = NULL,
    outside_opt_label = NULL,
    include_outside_option = FALSE,
    draws       = c("store", "generate"),
    seed        = NULL,
    scramble    = c("permuted", "none", "owen"),
    keep_data = TRUE,
    nloptr_opts = NULL,
    weights_col = NULL,
    cluster_col = NULL,
    person_col = NULL,
    polish = c("none", "newton", "bhhh"),
    polish_control = list(),
    concentrate_out_asc = c("none", "profile_likelihood", "blp"),
    concentrate_control = list()
) {
  cl <- match.call()

  se_method_default <- missing(se_method)
  se_method <- match.arg(se_method)
  if (!is.null(cluster_col) && se_method_default) se_method <- "cluster"
  scale_vars <- match.arg(scale_vars)
  polish <- match.arg(polish)
  polish_control <- .polish_control(polish_control)
  concentrate_out_asc <- match.arg(concentrate_out_asc)
  concentrate_control <- .concentrate_control(concentrate_control)
  # The concentration of the constants, by method. "profile_likelihood"
  # maximizes over them at every value of the other parameters (an inner
  # Newton solve of their first-order condition: observed counts equal to
  # posterior-weighted predicted counts), which leaves the estimator
  # unchanged, and differentiates the profile by the envelope theorem.
  # "blp" is reserved for matching aggregate shares instead (Berry,
  # Levinsohn and Pakes 1995), which changes the estimator: its inner solve
  # would be a share contraction such as mxl_blp_contraction()'s, but on the
  # estimation's draws (a panel fit's are per decision maker, the
  # contraction's per choice situation) and to a tight tolerance (Dube, Fox
  # and Su 2012); the outer gradient would need the implicit-function term
  # -(ds/d delta)^-1 ds/d theta_c, since the shares' solution does not set
  # the likelihood's score of the constants to zero; and the variance would
  # be that of the share-constrained (micro-BLP) estimator.
  if (identical(concentrate_out_asc, "blp")) {
    stop("concentrate_out_asc = \"blp\" is reserved for a share-matching ",
         "(BLP) concentration of the constants, which is not implemented.",
         call. = FALSE)
  }
  profile <- identical(concentrate_out_asc, "profile_likelihood")
  draws   <- match.arg(draws)
  scramble <- match.arg(scramble)
  if (identical(scramble, "owen")) {
    warning("scramble = \"owen\" is a deprecated alias for ",
            "scramble = \"permuted\". The implemented position-wise digit ",
            "permutation is not Owen's nested-uniform scramble.",
            call. = FALSE)
    scramble <- "permuted"
  }

  # Validate seed parameter
  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) || seed < 0) {
      stop("'seed' must be NULL or a single non-negative integer.")
    }
    seed <- as.integer(seed)
    if (draws != "generate") {
      message("'seed' is ignored when draws = 'store'.")
    }
  }

  # Backward compatibility: nloptr_opts -> optimizer + control
  if (!is.null(nloptr_opts)) {
    message("'nloptr_opts' is deprecated. Use 'optimizer' and 'control' instead.")
    optimizer <- optimizer %||% "nloptr"
    control <- nloptr_opts
  }

  # --- Resolve input pathway --------------------------------------------------
  has_data <- !is.null(data)
  has_input <- !is.null(input_data)
  cs_meta <- if (has_data) attr(data, "choice_sampling") else attr(input_data, "choice_sampling")

  if (has_data && has_input) {
    stop("Supply either 'data' (convenience) or 'input_data' (advanced), not both.")
  }
  if (!has_data && !has_input) {
    stop("Supply either 'data' (convenience) or 'input_data' (advanced).")
  }
  if (has_input && !is.null(weights_col)) {
    stop("`weights_col` is only supported in the convenience (data) workflow. ",
         "Bake weights into `input_data` via prepare_mxl_data(weights_col = ) ",
         "or supply `weights` to prepare_mxl_data().")
  }
  if (has_input && !is.null(cluster_col)) {
    stop("`cluster_col` is only supported in the convenience (data) workflow. ",
         "Bake cluster labels into `input_data` via ",
         "prepare_mxl_data(cluster_col = ).")
  }
  if (has_input && !is.null(person_col)) {
    stop("`person_col` is only supported in the convenience (data) workflow. ",
         "Bake the panel structure into `input_data` via ",
         "prepare_mxl_data(person_col = ).")
  }

  if (has_data) {
    # Convenience workflow: validate required column-name arguments
    if (is.null(id_col) || is.null(alt_col) || is.null(choice_col) ||
        is.null(covariate_cols) || is.null(random_var_cols)) {
      stop("Convenience workflow requires: id_col, alt_col, choice_col, ",
           "covariate_cols, and random_var_cols.")
    }
    # WESML provenance present but no weights supplied: auto-adopt the recorded
    # weight column, or error -- never silently fit unweighted under a WESML label.
    # (With person_col, prepare_mxl_data() rejects WESML provenance outright.)
    if (has_data && !is.null(cs_meta) && is.null(person_col) &&
        is.null(weights) && is.null(weights_col)) {
      wn <- cs_meta$weight_name
      if (!is.null(wn) && wn %in% names(data)) {
        weights_col <- wn
        message("Detected WESML choice-based-sampling provenance; applying attached ",
                "weights from column '", wn, "'.")
      } else {
        stop("Data carries WESML choice-based-sampling provenance but no weights were ",
             "supplied, and the recorded weight column (",
             if (is.null(wn)) "unknown" else paste0("'", wn, "'"),
             ") is not present in `data`. Pass `weights_col=` or `weights=` explicitly.",
             call. = FALSE)
      }
    }
    input_data <- prepare_mxl_data(
      data = data,
      id_col = id_col,
      alt_col = alt_col,
      choice_col = choice_col,
      covariate_cols = covariate_cols,
      random_var_cols = random_var_cols,
      weights = weights,
      weights_col = weights_col,
      outside_opt_label = outside_opt_label,
      include_outside_option = include_outside_option,
      rc_correlation = rc_correlation,
      cluster_col = cluster_col,
      person_col = person_col
    )
  } else {
    # Advanced workflow
    if (draws == "generate") {
      stop("draws=generate is only supported in the convenience workflow. ",
           "In the advanced workflow, supply 'eta_draws' directly.")
    }
    if (is.null(eta_draws)) {
      stop("'eta_draws' is required when using 'input_data' (advanced workflow).")
    }
    # `S` is a convenience-workflow argument; record the draws actually used so
    # post-hoc methods regenerate the same Halton blocks.
    S <- dim(eta_draws)[2L]
  }

  if (se_method == "cluster" && is.null(input_data$cluster)) {
    stop("se_method = \"cluster\" needs cluster labels: pass `cluster_col=` ",
         "(convenience workflow) or prepare `input_data` with ",
         "prepare_mxl_data(cluster_col = ).", call. = FALSE)
  }

  # Parameter dimensions. alt_mapping has one row per alternative, the outside
  # option's first when there is one, so use_asc frees J - 1 constants: every
  # inside alternative's with the outside option, all but the first's
  # (normalized to 0) without it.
  J <- nrow(input_data$alt_mapping)
  K_x <- ncol(input_data$X)
  K_w <- ncol(input_data$W)
  rc_correlation <- input_data$rc_correlation
  L_size <- if (rc_correlation) K_w * (K_w + 1) / 2 else K_w
  mu_size <- if (rc_mean) K_w else 0
  n_asc <- if (use_asc) J - 1 else 0
  n_params <- K_x + mu_size + L_size + n_asc

  if (is.null(rc_dist)) rc_dist <- rep(0L, K_w)

  # Parameter index map (built early so the scaling layer can address blocks)
  pos <- 0
  param_map <- list(beta = seq_len(K_x))
  pos <- K_x
  if (mu_size > 0) {
    param_map$mu <- pos + seq_len(mu_size)
    pos <- pos + mu_size
  }
  param_map$sigma <- pos + seq_len(L_size)
  pos <- pos + L_size
  if (n_asc > 0) param_map$asc <- pos + seq_len(n_asc)

  # Parameter names (built early; reused for theta_hat, vcov, se downstream)
  beta_names <- colnames(input_data$X)
  mu_names <- if (rc_mean) paste0("Mu_", colnames(input_data$W)) else character(0)
  sigma_names <- .mxl_cov_names("L", K_w, rc_correlation)
  alt_col <- names(input_data$alt_mapping)[2]
  asc_names <- if (use_asc) {
    paste0("ASC_", input_data$alt_mapping[[alt_col]][2:J])
  } else {
    character(0)
  }
  param_names <- c(beta_names, mu_names, sigma_names, asc_names)
  # summary() relabels some parameters: its labels may not repeat either
  .check_param_names(
    list(param_names,
         .mxl_summary_labels(param_names, param_map, K_w, rc_dist,
                             rc_correlation, rc_mean)),
    c("ASC_<label>", "Mu_<variable>", "L_<i><j>", "Sigma_<i><j>",
      "exp(Mu_<variable>)")
  )

  # Normalize lower/upper bounds (natural units); a box that holds no start
  # is an error here, before the draws, rather than at the optimizer.
  lower <- .normalize_bound(lower, param_names, -Inf, "lower")
  upper <- .normalize_bound(upper, param_names,  Inf, "upper")
  bad <- lower > upper | lower == Inf | upper == -Inf
  if (any(bad)) {
    stop("Invalid bounds for ", paste(param_names[bad], collapse = ", "),
         ": `lower` must not exceed `upper`, `lower` must not be Inf, and ",
         "`upper` must not be -Inf.", call. = FALSE)
  }

  # The variance's inputs, checked here, before the draws, the column scales,
  # the "bhhh" pass and the optimizer, as the variance step checks them: the
  # weights against se_method and the choice-based-sampling provenance, then
  # the cluster labels, one per likelihood unit (nesting decision makers with
  # person_col). A Ti that does not split the situations into decision makers
  # is left to the kernels, which report it at their first call, as before.
  choice_sampling <- .resolve_choice_sampling(
    weights = input_data$weights, se_method = se_method,
    cs_meta = cs_meta, has_input = has_input, prepare_fn = "prepare_mxl_data"
  )
  Ti <- input_data$Ti
  if (se_method == "cluster" &&
      (is.null(Ti) || (is.numeric(Ti) && !anyNA(Ti) && all(Ti >= 1) &&
                       isTRUE(sum(as.numeric(Ti)) == length(input_data$M))))) {
    .check_cluster_labels(
      .to_units(input_data$cluster, input_data, "`cluster_col`"),
      length(.unit_first(input_data))
    )
  }

  # Constants without a finite estimate (alternatives never chosen, or chosen
  # wherever offered), warned before a fit that can take hours; a finite bound
  # in the direction a constant diverges holds it. A polish holds the constants
  # of the alternatives named, the reference's case aside, which moves every
  # constant together (row r >= 2 of alt_mapping has the constant asc[r - 1]).
  # The profile likelihood needs every constant to have a finite estimate:
  # those cases are errors there, and so are bounds on the constants.
  unbounded <- NULL
  if (profile) {
    .check_profile_inputs(n_asc, lower[param_map$asc], upper[param_map$asc],
                          param_names[param_map$asc], input_data$alt_mapping,
                          concentrate_control, input_data)
    if (!is.null(theta_init)) .check_theta_init(theta_init, n_params)
  } else if (n_asc > 0) {
    unbounded <- .warn_unbounded_constants(input_data$alt_mapping,
                                           asc_lower = lower[param_map$asc],
                                           asc_upper = upper[param_map$asc])
  }
  held <- integer()
  if (polish != "none" && !is.null(unbounded)) {
    rows <- sort(c(unbounded$never, unbounded$always))
    held <- param_map$asc[rows[rows >= 2L] - 1L]
  }
  # The BHHH matrix of U likelihood units has rank at most U; the parameters
  # it moves exclude the held and the pinned (lower == upper)
  if (polish == "bhhh") {
    n_units <- length(.unit_first(input_data))
    n_moved <- sum(!seq_len(n_params) %in% held & lower < upper)
    if (n_units < n_moved) {
      stop("polish = \"bhhh\" needs at least as many likelihood units as the ",
           "parameters it moves: the BHHH matrix of ", n_units, " ",
           if (is.null(input_data$Ti)) "choice situations" else "decision makers",
           " has rank at most ", n_units, ", below the ", n_moved,
           " parameters. Use polish = \"newton\".", call. = FALSE)
    }
  }

  # Draws for the convenience workflow, built once the inputs have passed
  # the checks above: a store-mode cube can take gigabytes.
  if (has_data) {
    if (draws == "store") {
      # One K_w x S draw block per likelihood unit: per decision maker with
      # person_col, per choice situation otherwise.
      n_units <- length(.unit_first(input_data))
      .warn_store_cube(S, n_units, K_w, !is.null(input_data$Ti))
      eta_draws <- get_halton_normals(S, n_units, K_w)
    } else {
      # generate mode: no cube ever materialized; empty placeholder
      eta_draws <- array(0, dim = c(K_w, 0L, 0L))
      # Draw seed from R RNG when not supplied (like run_mnprobit)
      if (is.null(seed)) {
        seed <- sample.int(.Machine$integer.max, 1L)
      }
    }
  }

  # --- The optimizer's coordinates: column scales (scale_vars) ----------------
  # scale_vars chooses the optimizer's coordinates, theta = m * theta_t + c
  # (.coordinate_map()); the kernels evaluate the model at theta on the
  # natural design, and estimates and variances are reported in natural
  # units. "sd", "mad" and "iqr" take the map from the columns' scales
  # (computed here; log-normal W columns at 1), so that the optimizer sees
  # the likelihood of the design with its columns divided by them; "bhhh"
  # takes it from the BHHH diagonal at the start values. theta_init and the
  # bounds are in natural units and mapped below; sX and sW are returned as
  # 1s unless the columns are scaled.
  col_scaled <- scale_vars %in% c("sd", "mad", "iqr")
  sX <- rep(1, K_x); names(sX) <- colnames(input_data$X)
  sW <- rep(1, K_w); names(sW) <- colnames(input_data$W)
  if (col_scaled) {
    if (K_x > 0) {
      sX_raw <- .column_scales(input_data$X, scale_vars)
      .assert_scales_ok(sX_raw, scale_vars, "fixed-coefficient")
      sX <- sX_raw
    }
    if (K_w > 0) {
      sW_raw <- .column_scales(input_data$W, scale_vars)
      normal_cols <- which(rc_dist == 0L)
      if (length(normal_cols) > 0L) {
        .assert_scales_ok(sW_raw, scale_vars, "normal random-coefficient",
                          idx = normal_cols)
      }
      # Preserve names from sW_raw; carve out log-normal columns (their
      # parameters keep the identity map).
      sW <- sW_raw
      sW[rc_dist == 1L] <- 1
      n_lognormal <- sum(rc_dist == 1L)
      if (K_w > 0L && n_lognormal == K_w) {
        message("scale_vars='", scale_vars,
                "': all random coefficients are log-normal; their parameters ",
                "keep their scale.")
      } else if (n_lognormal > 0L) {
        message("scale_vars='", scale_vars,
                "': the parameters of log-normal random coefficient(s) keep ",
                "their scale (no closed-form rescaling).")
      }
    }
  }
  # Resolve theta_init (natural units). The default start is in the columns'
  # units (.mxl_default_start()), moved to the nearest bound where it falls
  # outside the box (the optimizer refuses a start outside it), with a
  # message. A supplied theta_init is used as given.
  start_default <- is.null(theta_init)
  if (start_default) {
    theta_init <- .mxl_default_start(n_params, param_map, input_data$W,
                                     rc_dist, rc_correlation)
    moved <- theta_init < lower | theta_init > upper
    if (any(moved)) {
      # Parameters pinned by lower == upper are moved without a word.
      named <- param_names[moved & lower < upper]
      if (length(named)) {
        message("The default start of ", paste(utils::head(named, 5L), collapse = ", "),
                if (length(named) > 5L) paste0(" and ", length(named) - 5L, " more"),
                " lies outside `lower`/`upper`; it starts at the nearest bound.")
      }
      theta_init <- pmin(pmax(theta_init, lower), upper)
    }
  }

  # Resolve generate-mode parameters for C++ kernels.
  # In store mode (draws="store"): gen_seed_cpp = -1L triggers cube path (unchanged behavior).
  gen_seed_cpp     <- if (draws == "generate") seed else -1L
  gen_scramble_cpp <- if (draws == "generate") (if (scramble == "permuted") 1L else 0L) else 1L
  gen_S_cpp        <- if (draws == "generate") S else 0L

  # The model's objective and gradient at theta; `...` takes opg_diag.
  model_f <- function(theta, ...) {
    mxl_loglik_gradient_parallel(
      theta = theta,
      X = input_data$X,
      W = input_data$W,
      alt_idx = input_data$alt_idx,
      choice_idx = input_data$choice_idx,
      M = input_data$M,
      weights = input_data$weights,
      rc_dist = rc_dist,
      rc_correlation = rc_correlation,
      rc_mean = rc_mean,
      eta_draws = eta_draws,
      use_asc = use_asc,
      include_outside_option = input_data$include_outside_option,
      gen_seed = gen_seed_cpp,
      gen_scramble = gen_scramble_cpp,
      gen_S = gen_S_cpp,
      Ti = input_data$Ti,
      ...
    )
  }

  # --- The profile likelihood of the constants (concentrate_out_asc) ---------
  # The constants delta are solved at the start's continuous parameters
  # theta_c, from the log ratio of each alternative's weighted chosen count to
  # the reference's under the default start (theta_init's constants
  # otherwise), so that the "bhhh" pass and the outer optimizer start on the
  # profile.
  if (profile) {
    c_idx <- setdiff(seq_len(n_params), param_map$asc)
    d_idx <- param_map$asc
    cc <- .chosen_counts(input_data, n_asc)
    if (start_default) theta_init[d_idx] <- log(cc$obs / cc$ref)
    delta_kernel <- function(theta) {
      mxl_delta_hessian_parallel(
        theta = theta, X = input_data$X, W = input_data$W,
        alt_idx = input_data$alt_idx, choice_idx = input_data$choice_idx,
        M = input_data$M, weights = input_data$weights, eta_draws = eta_draws,
        rc_dist = rc_dist, rc_correlation = rc_correlation, rc_mean = rc_mean,
        use_asc = use_asc,
        include_outside_option = input_data$include_outside_option,
        gen_seed = gen_seed_cpp, gen_scramble = gen_scramble_cpp,
        gen_S = gen_S_cpp, Ti = input_data$Ti
      )
    }
    profile_solve <- .mxl_profile_solver(delta_kernel, c_idx, d_idx, cc$obs,
                                         cc$ref, mean(input_data$weights),
                                         concentrate_control)
    t_first <- system.time(
      s0 <- profile_solve(theta_init[c_idx], theta_init[d_idx],
                          min(4 * concentrate_control$max_iter,
                              .Machine$integer.max)),
      gcFirst = FALSE)
    if (!is.null(s0$error)) {
      # an error, not a failure to converge (a kernel's memory check, say)
      stop("concentrate_out_asc = \"profile_likelihood\": the constants' ",
           "solve at the start values stopped with an error: ", s0$error,
           call. = FALSE)
    }
    if (!s0$converged) {
      stop(errorCondition(paste0(
        "concentrate_out_asc = \"profile_likelihood\": the constants could ",
        "not be solved at the start values (", s0$message, "). Start ",
        "elsewhere (theta_init), allow more steps ",
        "(concentrate_control$max_iter), or fit without ",
        "concentrate_out_asc."), class = "choicer_profile_unconverged",
        call = NULL))
    }
    theta_init[d_idx] <- s0$delta
  }

  # --- The optimizer's coordinates --------------------------------------------
  #   theta_natural = map$scale * theta_t + map$shift
  # "bhhh": the BHHH diagonal at the start values, from one gradient pass;
  # unusable where the likelihood overflows there (unit entries: the identity).
  if (scale_vars != "none") .check_theta_init(theta_init, n_params)
  opg0 <- NULL
  t_bpass <- system.time(if (scale_vars == "bhhh") {
    g0 <- model_f(theta_init, opg_diag = TRUE)
    if (isTRUE(g0$overflow)) {
      message("scale_vars = \"bhhh\": the likelihood overflows at the start ",
              "values (check theta_init); the parameters are left unscaled.")
      opg0 <- rep(1, n_params)
    } else {
      opg0 <- g0$opg_diag
    }
  }, gcFirst = FALSE)
  map <- .coordinate_map(scale_vars, param_map, n_params, sX = sX, sW = sW,
                         rc_correlation = rc_correlation, opg_diag = opg0)

  # Run the optimizer in the map's coordinates on the objective with the
  # kernel's overflow sentinel kept above the path
  elapsed <- system.time({
    if (!profile) {
      opt <- .optimize_in_coordinates(
        map = if (scale_vars != "none") map, optimizer = optimizer,
        theta_init = theta_init, eval_f = .lift_sentinel(model_f),
        lower = lower, upper = upper, control = control
      )
    } else {
      # The outer problem: the optimizer on theta_c alone, in the map's
      # coordinates restricted to it, on the negated profile log-likelihood
      # and its envelope gradient (.profile_objective()), warm-started from
      # the best converged point.
      prof <- .profile_objective(
        function(theta_c, delta0) {
          profile_solve(theta_c, delta0, concentrate_control$max_iter)
        },
        c_idx, list(theta_c = theta_init[c_idx], value = s0$value,
                    delta = s0$delta, gradient = s0$gradient,
                    residual = s0$residual, decrement = s0$decrement))
      opt <- .optimize_in_coordinates(
        map = if (scale_vars != "none") {
          list(scale = map$scale[c_idx], shift = map$shift[c_idx])
        },
        optimizer = optimizer, theta_init = theta_init[c_idx],
        eval_f = .lift_sentinel(prof$f, f_max = max(0, s0$value)),
        lower = lower[c_idx], upper = upper[c_idx], control = control
      )
      # The constants at the optimizer's solution: the best point's, when it
      # is that point; solved from it otherwise; the best point itself when
      # that solve fails (the optimizer then returned a point it never
      # evaluated well)
      st <- prof$state
      best <- st$best
      if (identical(unname(opt$par), unname(best$theta_c), num.eq = FALSE)) {
        s_hat <- c(best, list(converged = TRUE))
      } else {
        s_hat <- profile_solve(opt$par, best$delta, concentrate_control$max_iter)
        st$solves <- st$solves + 1L
        st$inner_iterations <- st$inner_iterations + s_hat$iterations
        st$passes <- st$passes + s_hat$passes
        if (isTRUE(s_hat$converged)) {
          st$max_residual <- max(st$max_residual, s_hat$residual)
        }
      }
      final_ok <- isTRUE(s_hat$converged)
      final_error <- if (!final_ok) s_hat$error
      theta_c_hat <- if (final_ok) opt$par else best$theta_c
      if (!final_ok) s_hat <- c(best, list(converged = TRUE))
      opt$par <- numeric(n_params)
      opt$par[c_idx] <- theta_c_hat
      opt$par[d_idx] <- s_hat$delta
      opt$value <- s_hat$value
    }
  })

  elapsed <- elapsed + t_bpass  # "bhhh"'s gradient pass counts
  conc_record <- NULL
  if (profile) {
    elapsed <- elapsed + t_first
    conc_record <- list(
      method = "profile_likelihood", evaluations = st$evaluations,
      solves = st$solves + 1L, failed = st$failed,
      inner_iterations = st$inner_iterations + s0$iterations,
      passes = st$passes + s0$passes,
      max_residual = max(st$max_residual, s0$residual),
      residual = s_hat$residual, decrement = s_hat$decrement,
      elapsed_time = elapsed[["elapsed"]])
    message("Profile likelihood: ", st$evaluations, " evaluations of the ",
            "optimizer, ", conc_record$solves, " inner solves (",
            conc_record$inner_iterations, " Newton steps, ",
            conc_record$passes, " passes)",
            if (st$failed > 0L) {
              paste0(", ", st$failed, " of them unconverged at trial points ",
                     "(overflows to the optimizer)")
            })
    # Failed solves before the best point's last improvement are the line
    # search's ordinary backtracking; later ones may have stopped the optimizer
    if (st$failed_late > 0L || !final_ok) {
      why <- paste0(
        if (st$failed_late > 0L) {
          paste0(st$failed_late, " of the inner solves after the optimizer's ",
                 "best point did not converge (it saw them as overflows, so ",
                 "its success code can mislead",
                 if (!is.null(st$error)) {
                   paste0("; the last error: ", st$error)
                 }, ")")
        },
        if (st$failed_late > 0L && !final_ok) "; ",
        if (!final_ok) {
          paste0("the constants could not be solved at the optimizer's ",
                 "solution",
                 if (!is.null(final_error)) paste0(" (", final_error, ")"),
                 ", so the best point solved is returned")
        })
      opt$message <- paste0(paste(as.character(opt$message), collapse = "; "),
                            "; profile likelihood: ", why)
      warning(warningCondition(paste0(
        "concentrate_out_asc = \"profile_likelihood\": ", why,
        ". polish = \"newton\" certifies the maximum, or not."),
        class = "choicer_profile_unconverged"))
    }
  }
  message("Optimization run time ", convertTime(elapsed))

  # --- The Newton phase (polish) ----------------------------------------------
  # From the optimizer's solution, damped Newton steps on the analytic Hessian
  # ("newton") or the BHHH matrix ("bhhh", its weights unsquared, the
  # objective's curvature) of the same draws, in natural units
  # (run_polish()), on the raw objective, whose overflow flag fails a trial
  # point. The estimator is unchanged. A converged polish reports the
  # optimizer's success code; the optimizer's own status is kept in the record.
  polish_record <- NULL
  curv <- NULL
  if (polish != "none") {
    curv_kernel <- if (polish == "newton") mxl_hessian_parallel else mxl_bhhh_parallel
    curv_f <- function(theta) {
      curv_kernel(
        theta = theta, X = input_data$X, W = input_data$W,
        alt_idx = input_data$alt_idx, choice_idx = input_data$choice_idx,
        M = input_data$M, weights = input_data$weights, eta_draws = eta_draws,
        rc_dist = rc_dist, rc_correlation = rc_correlation, rc_mean = rc_mean,
        use_asc = use_asc,
        include_outside_option = input_data$include_outside_option,
        gen_seed = gen_seed_cpp, gen_scramble = gen_scramble_cpp,
        gen_S = gen_S_cpp, Ti = input_data$Ti
      )
    }
    # The polish's last curvature, evaluated at the point it returns (unless
    # it returns the start, or its curvature fails), is the variance's own
    # matrix there when se_method needs that one: the Hessian as the bread of
    # "hessian", "sandwich" and "cluster", the BHHH matrix with the weights
    # unsquared for "bhhh" (the robust meat squares them). Only then is it
    # kept, and taken out of the result at once; otherwise it is released
    # with the polish. The variance's peak then does not grow for "hessian"
    # and "bhhh"; for "sandwich" and "cluster" the bread is held while the
    # meat is computed, which adds at most one p x p where the meat's own
    # buffers (the squared weights; the cluster kernel's 8 (U + 2G) bytes and
    # labels) outweigh the Hessian kernel's.
    reuse <- switch(polish, newton = se_method != "bhhh", bhhh = se_method == "bhhh")
    t_pol <- system.time(
      pol <- run_polish(stats::setNames(opt$par, param_names), model_f, curv_f,
                        lower = lower, upper = upper, scale = map$scale,
                        control = polish_control, method = polish,
                        fixed = held, keep_curvature = reuse),
      gcFirst = FALSE)
    curv <- pol$curvature_at_par
    pol$curvature_at_par <- NULL
    # Where the reference (the outside option, or the base alternative) is
    # never chosen, or chosen wherever offered, the other constants diverge
    # together (up, or down): there is no maximum to certify. A rise stops
    # where the constant of an alternative that is chosen (row r >= 2 of
    # alt_mapping, constant asc[r - 1]) is held or bounded above; a fall is
    # not examined, and is never certified.
    if (pol$converged && !is.null(unbounded) &&
        (unbounded$ref_never || unbounded$ref_always)) {
      chosen <- input_data$alt_mapping[["N_CHOICES"]][-1L] > 0
      stops <- param_map$asc %in% held | is.finite(upper[param_map$asc])
      if (!(unbounded$ref_never && any(chosen & stops))) {
        pol$converged <- FALSE
        pol$message <- paste0(
          pol$message, ", but the constants diverge together (the ",
          if (unbounded$outside) "outside option" else "base alternative",
          if (unbounded$ref_never) " is never chosen" else
            " is chosen wherever it is offered",
          "), so there is no maximum to certify")
      }
    }
    label <- if (polish == "newton") "Newton polish" else "BHHH polish"
    its <- paste(pol$iterations,
                 if (pol$iterations == 1L) "iteration" else "iterations")
    dec <- paste0(if (polish == "newton") "Newton" else "BHHH", " decrement ",
                  format(signif(pol$decrement, 3)))
    polish_record <- c(
      pol[c("method", "converged", "iterations", "damped", "decrement",
            "lambda", "message", "evaluations", "active", "held", "flat",
            "improvement")],
      list(error = pol$error, elapsed_time = t_pol[["elapsed"]],
           optimizer = list(convergence = opt$convergence,
                            message = opt$message,
                            iterations = opt$iterations,
                            elapsed_time = elapsed[["elapsed"]])))
    if (is.finite(pol$value)) {
      opt$par <- unname(pol$par)
      opt$value <- pol$value
    }
    opt_message <- paste(as.character(opt$message), collapse = "; ")
    if (!nzchar(opt_message)) opt_message <- NULL
    if (pol$converged) {
      opt$convergence <- if (identical(optimizer %||% "nloptr", "nloptr")) 1L else 0L
      opt$message <- paste0(label, " converged after ", its, " (", dec, ")",
                            if (!is.null(opt_message)) {
                              paste0("; the optimizer: ", opt_message)
                            })
    } else {
      opt$message <- paste0(if (!is.null(opt_message)) paste0(opt_message, "; "),
                            label, " not converged after ", its, ": ",
                            pol$message, " (", dec, ")")
    }
    message(label, ": ", if (pol$converged) "converged" else "not converged",
            " after ", its, " (", dec, ")",
            if (!pol$converged) paste0(": ", pol$message), ", ",
            convertTime(t_pol))
    if (!is.null(pol$error)) {
      warning(warningCondition(paste0(
        "The ", label, " stopped with an error, so the estimates are its ",
        "last accepted point (the optimizer's solution when it took no ",
        "step): ", pol$error), class = "choicer_polish_failed"))
    }
    elapsed <- elapsed + t_pol
  }

  # Estimate at the optimum (natural units)
  theta_hat <- opt$par
  names(theta_hat) <- param_names

  # Reconstruct Sigma for display (from the natural-scale L params)
  L_params <- theta_hat[param_map$sigma]
  sigma_mat <- build_var_mat(L_params, K_w, rc_correlation)
  w_names <- colnames(input_data$W)
  if (!is.null(w_names)) {
    rownames(sigma_mat) <- w_names
    colnames(sigma_mat) <- w_names
  }

  # Draws info (metadata only, not the full array). N is the number of choice
  # situations, the draw-block count of the per-situation prediction sites;
  # estimation-type sites use one block per likelihood unit (.unit_first()).
  draws_info <- list(
    S       = S,
    N       = input_data$N,
    K_w     = K_w,
    mode    = draws,
    seed    = if (draws == "generate") seed else NULL,
    scramble = if (draws == "generate") scramble else NULL
  )

  # The optimizer's record, with the polish's only when there is one
  opt_info <- list(
    name = if (is.function(optimizer)) "custom" else (optimizer %||% "nloptr"),
    control = control,
    elapsed_time = elapsed[["elapsed"]],
    iterations = opt$iterations
  )
  if (!is.null(polish_record)) opt_info$polish <- polish_record
  if (!is.null(conc_record)) opt_info$concentration <- conc_record

  # Build S3 object. It holds the data until its variance is computed below.
  fit <- new_choicer_mxl(
    call = cl,
    coefficients = theta_hat,
    loglik = -opt$value,
    nobs = input_data$N,
    n_params = n_params,
    convergence = opt$convergence,
    message = opt$message,
    data_spec = input_data$data_spec,
    alt_mapping = input_data$alt_mapping,
    param_map = param_map,
    use_asc = use_asc,
    include_outside_option = input_data$include_outside_option,
    optimizer = opt_info,
    vcov = NULL,
    se = NULL,
    data = list(
      X = input_data$X,
      W = input_data$W,
      alt_idx = input_data$alt_idx,
      choice_idx = input_data$choice_idx,
      M = input_data$M,
      weights = input_data$weights,
      cluster = input_data$cluster,
      situation_ids = input_data$situation_ids,
      Ti = input_data$Ti,
      person_ids = input_data$person_ids
    ),
    draws_info = draws_info,
    rc_dist = rc_dist,
    rc_correlation = rc_correlation,
    rc_mean = rc_mean,
    sigma = sigma_mat,
    se_method = se_method,
    scale_vars = scale_vars,
    sX = sX,
    sW = sW,
    choice_sampling = choice_sampling,
    n_persons = if (!is.null(input_data$Ti)) length(input_data$Ti),
    param_scale = stats::setNames(map$scale, param_names),
    param_shift = stats::setNames(map$shift, param_names)
  )

  # Compute vcov eagerly using the selected SE method, by the route
  # vcov(fit, type = ) takes post hoc, with the fit's own draws: at the
  # natural estimates, in natural units, inverted equilibrated. For
  # "sandwich" (robust / WESML) standard errors, V = A^{-1} B A^{-1} with
  # bread A = weighted negated Hessian and meat B = weight-squared OPG (the
  # BHHH kernel with weights^2, its per-unit score being weight-free). For
  # "cluster", the meat is the outer product of within-cluster sums of
  # weighted scores, summed in C++ without the U x p score matrix. Scores,
  # weights and cluster labels are per likelihood unit (decision maker with
  # person_col). After a polish, its last curvature stands in for the matrix
  # the variance would otherwise evaluate at the same point (curv, above). An
  # error here (a kernel's memory check, say) returns the estimates without a
  # variance, with a warning (.vcov_failed()).
  vcov_result <- tryCatch(
    .assemble_score_vcov(
      fit, type = if (se_method == "sandwich") "robust" else se_method,
      gp = list(eta_draws = eta_draws, gen_seed = gen_seed_cpp,
                gen_scramble = gen_scramble_cpp, gen_S = gen_S_cpp),
      what = "`cluster_col`",
      bread = if (polish == "newton") curv,
      bhhh = if (polish == "bhhh") curv
    ),
    error = function(e) {
      .vcov_failed(e, se_method, param_names, keep_data)
    }
  )
  curv <- NULL
  if (!is.null(vcov_result$vcov)) {
    rownames(vcov_result$vcov) <- param_names
    colnames(vcov_result$vcov) <- param_names
    names(vcov_result$se) <- param_names
  }
  fit["vcov"] <- list(vcov_result$vcov)
  fit["se"] <- list(vcov_result$se)
  # A variance step that failed is recorded, on such fits only
  if (!is.null(vcov_result$error)) {
    fit["variance_error"] <- list(vcov_result$error)
  }
  if (keep_data) return(fit)
  fit["data"] <- list(NULL)  # the slot is kept, NULL
  fit
}


#' Prepare inputs for mixed logit estimation
#'
#' Prepares and validates inputs for mixed logit estimation routine.
#'
#' Rows are ordered by choice-situation id and, within a situation, by
#' alternative. With \code{person_col}, rows are ordered by decision maker
#' first, then by situation id and alternative, so that each decision maker's
#' situations are contiguous. A positional \code{weights} vector is read in
#' the prepared order, which is ascending id in the cross-section; with
#' \code{person_col} it is accepted only when the prepared order is also
#' ascending id, and an error otherwise points to \code{weights_col}, which is
#' aligned by id and is the safer interface.
#'
#' @param data Data frame containing choice data
#' @param id_col Name of the column identifying choice situations
#' @param alt_col Name of the column identifying alternatives
#' @param choice_col Name of the column indicating chosen alternative (1 = chosen, 0 = not chosen)
#' @param covariate_cols Vector of names of columns to be used as covariates
#' @param random_var_cols Vector of names of columns to be used as random variables
#' @param weights Optional vector of weights, one per choice situation in the
#'   prepared order (ascending id; see Details). If NULL, equal weights are
#'   used. All weights
#'   must be finite and strictly positive. With \code{person_col} they are
#'   decision-maker weights and must be constant within each decision maker.
#' @param weights_col Optional name of a column in \code{data} holding a
#'   per-row weight (constant within each choice situation, and within each
#'   decision maker when \code{person_col} is used; finite and strictly
#'   positive). Mutually exclusive with \code{weights}.
#' @param outside_opt_label The label of the outside option or of the base
#'   alternative, given in the class of the alternative column. With
#'   \code{include_outside_option = TRUE} it names the outside option's row of
#'   \code{alt_mapping}. Without the outside option it names the base
#'   alternative, coded first (with constants, its constant is normalized to 0);
#'   it must be one of the alternatives, and \code{NULL} (the default) makes the
#'   first alternative in sort order the base.
#' @param include_outside_option Logical indicating whether to include an outside option in the model.
#' @param rc_correlation Logical indicating whether random coefficients are correlated. Default is FALSE.
#' @param cluster_col Optional name of a column in \code{data} holding cluster
#'   labels for cluster-robust standard errors. Must be constant within each
#'   \code{id_col}, and within each decision maker when \code{person_col} is
#'   used (clusters must nest decision makers); collapsed to one label per
#'   choice situation and returned as \code{cluster}.
#' @param person_col Optional name of the column identifying decision makers
#'   (respondents). When supplied, all choice situations of a decision maker
#'   share one draw of the random coefficients (panel likelihood);
#'   \code{id_col} must still identify choice situations uniquely across the
#'   data set. \code{NULL} (default) makes each choice situation its own
#'   decision maker (cross-sectional likelihood). Choice-based (WESML) samples
#'   are not supported with \code{person_col}: their weights vary with the
#'   chosen alternative, not by decision maker.
#' @returns A `choicer_data_mxl` object (list) containing:
#'   \itemize{
#'     \item `X`: Fixed-coefficient design matrix (sum(M) x K_x, double).
#'     \item `W`: Random-coefficient design matrix (sum(M) x K_w, double).
#'     \item `alt_idx`: Integer vector of alternative indices.
#'     \item `choice_idx`: Integer vector of chosen alternative indices.
#'     \item `M`: Integer vector with number of alternatives per choice situation.
#'     \item `N`: Number of choice situations.
#'     \item `weights`: Vector of weights.
#'     \item `cluster`: Vector of cluster labels (or `NULL`).
#'     \item `situation_ids`: Choice-situation ids in prepared (sorted) order.
#'     \item `Ti`: Number of choice situations of each decision maker, in
#'       prepared order (`NULL` without `person_col`).
#'     \item `person_ids`: Decision-maker ids in prepared order (`NULL`
#'       without `person_col`).
#'     \item `include_outside_option`: Logical flag.
#'     \item `rc_correlation`: Logical flag.
#'     \item `alt_mapping`: data.table mapping alternatives to summary statistics.
#'     \item `dropped_cols`: Names of columns dropped due to collinearity, if any.
#'     \item `data_spec`: List with column-name metadata (incl. `person_col`).
#'   }
#' @inheritSection prepare_mnl_data Column names
#' @examples
#' library(data.table)
#' set.seed(42)
#' N <- 50; J <- 3
#' dt <- data.table(id = rep(1:N, each = J), alt = rep(1:J, N))
#' dt[, `:=`(x1 = rnorm(.N), w1 = rnorm(.N), w2 = rnorm(.N))]
#' dt[, choice := 0L]
#' dt[, choice := sample(c(1L, rep(0L, J - 1))), by = id]
#' input <- prepare_mxl_data(dt, "id", "alt", "choice", "x1", c("w1", "w2"))
#' str(input$X)
#' str(input$W)
#' @export
prepare_mxl_data <- function(
    data,
    id_col,
    alt_col,
    choice_col,
    covariate_cols,
    random_var_cols,
    weights = NULL,
    outside_opt_label = NULL,
    include_outside_option = FALSE,
    rc_correlation = FALSE,
    weights_col = NULL,
    cluster_col = NULL,
    person_col = NULL
) {

  ## Preliminary housekeeping --------------------------------------------------
  # Capture any choice-based-sampling provenance before column drops / coercion,
  # so it can be carried onto the returned object for the advanced pathway.
  cs_provenance <- attr(data, "choice_sampling")
  if (!is.null(person_col) && !is.null(cs_provenance)) {
    stop("`person_col` cannot be used with a choice-based (WESML) sample: ",
         "WESML weights are defined per choice situation (they vary with the ",
         "chosen alternative), so they cannot weight the panel likelihood ",
         "sum_n w_n log L_n, which has one weight per decision maker. Fit the ",
         "cross-sectional likelihood with the WESML weights instead, clustering ",
         "on the decision maker via `cluster_col`.", call. = FALSE)
  }

  # Check if all relevant variables are available
  needed <- c(id_col, alt_col, choice_col, covariate_cols, random_var_cols)
  if (!is.null(weights) && !is.null(weights_col)) {
    stop("Supply only one of `weights` or `weights_col`.")
  }
  if (!is.null(weights_col)) needed <- c(needed, weights_col)
  if (!is.null(cluster_col)) needed <- c(needed, cluster_col)
  if (!is.null(person_col)) needed <- c(needed, person_col)
  # Names first: they need no data, and reading the data can stop for a
  # missing bit64.
  .check_col_names(needed, alt_col)
  # The covariates are read in place from `src` (the caller's data, never
  # modified, when possible); only the index columns are copied into `dt`,
  # and X and W are gathered by source row at the end.
  prep_src <- .prep_source(data, needed)
  src <- prep_src$src
  if (!all(needed %in% names(src)))
    stop("Missing columns: ",
         paste(setdiff(needed, names(src)), collapse = ", "))
  dt <- .copy_cols(src, c(id_col, alt_col, choice_col, weights_col,
                          cluster_col, person_col))
  dt[, .choicer_row := seq_len(.N)]
  .int64_cols_to_double(dt, weights_col)

  ## Drop ids with missing observations ----------------------------------------
  has_na <- .rows_with_na(src, prep_src$scan)
  ids_to_drop <- unique(dt[[id_col]][has_na[dt$.choicer_row]])
  if (length(ids_to_drop) > 0) {
    # computed outside dt[...], where a column named ids_to_drop would mask it
    keep <- !(dt[[id_col]] %in% ids_to_drop)
    dt <- dt[keep]
    warning("Removed ", length(ids_to_drop),
            " choice situations containing missing values.")
  }
  if (nrow(dt) == 0) {
    stop("All choice situations removed due to missing values.")
  }

  ## Sanity checks ---------------------------------------------------------

  ## covariates must be numeric
  if (!all(vapply(match(covariate_cols, names(src)),
                  function(j) is.numeric(.subset2(src, j)), NA)))
    stop("All covariates must be numeric.")
  if (!all(vapply(match(random_var_cols, names(src)),
                  function(j) is.numeric(.subset2(src, j)), NA)))
    stop("All covariates must be numeric.")

  ## choice column must be 0 or 1
  if (!.is_zero_one(dt[[choice_col]]))
    stop("`", choice_col, "` must contain only 0 and 1.")

  ## Panel: choice situations nest in decision makers. Checked before the
  ## per-id choice counts, which would misfire if ids restart within persons.
  if (!is.null(person_col)) {
    n_persons_per_id <- .n_distinct_by(dt, person_col, id_col)
    if (any(n_persons_per_id != 1L)) {
      stop("Each '", id_col, "' must belong to exactly one '", person_col,
           "': `id_col` must identify choice situations uniquely across ",
           "decision makers, but ", sum(n_persons_per_id != 1L), " id(s) ",
           "appear under several. If situation ids restart within each ",
           "decision maker, build a unique id, e.g. paste(", person_col, ", ",
           id_col, ").", call. = FALSE)
    }
  }

  ## Exactly one '1' per choice situation
  n_chosen <- .n_chosen_by(dt, choice_col, id_col)
  if (include_outside_option == FALSE && any(n_chosen != 1)) {
    stop("Each ", id_col, " must have exactly one chosen alternative (one '1' in ",
         choice_col, ").")
  }
  if (include_outside_option && any(n_chosen > 1)) {
    stop("Each ", id_col, " must have at most one chosen alternative (one '1' in ",
         choice_col, "). An id with no explicit choice is assumed to be outside option.")
  }

  ## Create integer alternative codes ------------------------------------------

  base_label <- !is.null(outside_opt_label) && include_outside_option==FALSE
  if (base_label) {
    levels <- c(outside_opt_label, sort(setdiff(unique(dt[[alt_col]]), outside_opt_label)))
  } else {
    levels <- sort(unique(dt[[alt_col]]))
  }

  dt <- .code_alternatives(dt, alt_col, levels)
  if (base_label) {
    .check_base_codes(dt$.choicer_alt_int, outside_opt_label, dt[[alt_col]],
                      alt_col, length(ids_to_drop))
  }

  ## Order rows ----------------------------------------------------------------
  ##   within each id: ascending alternative id
  ##   between ids   : ascending id
  ##   with person_col, first by decision maker, so that each decision maker's
  ##   situations are contiguous (the kernels' unit layout)
  data.table::setorderv(dt, c(person_col, id_col, ".choicer_alt_int"))

  ## index of each row within its choice set (computed outside dt[...], where
  ## a column named id_col would take the local's place)
  .choicer_idx <- data.table::rowidv(dt, cols = id_col)
  dt[, .choicer_idx_in_group := .choicer_idx]
  rm(.choicer_idx)

  ## Build objects -------------------------------------------------------------
  ## design matrix
  warn_once <- .int64_warn_once()  # a column can be in both X and W
  X <- warn_once(.gather_matrix(src, covariate_cols, dt$.choicer_row))
  X_res <- check_collinearity(X)
  X <- X_res$mat
  if (!is.null(X_res$dropped)) dropped_vars <- X_res$dropped # accumulate dropped vars if we had multiple checks


  W <- warn_once(.gather_matrix(src, random_var_cols, dt$.choicer_row))
  W_res <- check_collinearity(W)
  W <- W_res$mat
  if (!is.null(W_res$dropped)) {
     if(exists("dropped_vars")) dropped_vars <- c(dropped_vars, W_res$dropped)
     else dropped_vars <- W_res$dropped
  }

  dt[, .choicer_row := NULL]
  .check_index_covariates(c(covariate_cols, random_var_cols),
                          c(id_col, alt_col, choice_col, weights_col,
                            cluster_col, person_col))

  ## alternative ids used for delta coefficients
  alt_idx <- as.integer(dt$.choicer_alt_int)                    # length == sum(M)

  ## M[i] - # alternatives per choice situation
  M <- .n_rows_by(dt, id_col)                                   # length N

  ## N: number of individuals / choice situations
  ids <- dt[[id_col]][!duplicated(dt[[id_col]])]  # vector of ids in *current* order
  N   <- length(ids)

  ## Panel layout: Ti[u] consecutive situations for decision maker u, in
  ## sorted person order; NULL in the cross-section.
  Ti <- person_ids <- NULL
  if (!is.null(person_col)) {
    person <- .collapse_situation_col(dt, person_col, id_col, ids)
    person_ids <- unique(person)
    person_pos <- match(person, person_ids)
    Ti <- tabulate(person_pos, length(person_ids))
    stopifnot(!is.unsorted(person_pos))  # contiguous by construction (sort)
    # A positional `weights` vector is read in the prepared order, which is
    # now decision maker first; unless that is also ascending-id order, a
    # vector built in id order would be silently misaligned.
    if (!is.null(weights) && !.ids_sorted(ids)) {
      stop("With `person_col`, choice situations are reordered by decision ",
           "maker, so a positional `weights` vector is ambiguous. Supply the ",
           "weights as a column via `weights_col` (aligned by id).",
           call. = FALSE)
    }
  }

  ## Collapse a row-level weight column to one weight per choice situation.
  ## Done AFTER ordering/filtering so alignment is by id, never by position.
  if (!is.null(weights_col)) {
    if (!is.numeric(dt[[weights_col]])) {
      stop("`", weights_col, "` must be numeric.")
    }
    nuniq <- .n_distinct_by(dt, weights_col, id_col)
    if (any(nuniq != 1L)) {
      stop("`", weights_col, "` must be constant within each '", id_col,
           "' (one weight per choice situation).")
    }
    weights <- .first_by(dt, weights_col, id_col, ids)
    if (any(!is.finite(weights))) {
      stop("`", weights_col, "` produced non-finite weights.")
    }
  }

  ## Collapse a row-level cluster column to one label per choice situation
  ## (same alignment discipline as weights_col).
  cluster <- if (!is.null(cluster_col)) {
    .collapse_situation_col(dt, cluster_col, id_col, ids)
  }

  ## choice_idx[i] - 1-based index *within* the choice set data
  ## 0 == outside option (only if chosen = 0 for all inside options & include_outside_option == TRUE)
  chosen <- dt[[choice_col]] == 1
  if (include_outside_option) {
    # start with all-zero (everyone assumed to pick the outside good), then
    # match the chosen rows' ids (at most one per id) to the master index
    choice_idx <- integer(N)
    choice_idx[match(dt[[id_col]][chosen], ids)] <-
      dt$.choicer_idx_in_group[chosen]
  } else {
    # exactly one explicit choice per id
    choice_idx <- dt$.choicer_idx_in_group[chosen]
  }
  rm(chosen)

  # Weights default = 1
  if (is.null(weights)) weights <- rep(1, N)
  weights <- .int64_to_double(weights, "`weights`")

  ## Weights must be finite and strictly positive. Zero/negative weights would
  ## silently invalidate weighted and WESML sandwich inference (w in the bread,
  ## w^2 in the meat). Validated here so every resolution path (weights=,
  ## weights_col=, and the uniform default) is covered.
  if (any(!is.finite(weights))) {
    stop("Weights must be finite, but non-finite values (NA/NaN/Inf) were found.",
         call. = FALSE)
  }
  if (any(weights <= 0)) {
    stop("Weights must be strictly positive, but values <= 0 were found.",
         call. = FALSE)
  }

  ## Panel: the likelihood has one term per decision maker, so weights and
  ## cluster labels must be constant within each (identity otherwise).
  units <- list(M = M, Ti = Ti)
  .to_units(weights, units, "Weights")
  .to_units(cluster, units, "`cluster_col`")

  ## Alternative summary -------------------------------------------------------
  ## One inside-alternative aggregation; the outside branch only prepends its
  ## synthetic alt_int = 0 row.
  alt_mapping <- .alt_counts(dt, alt_col, choice_col)
  if (include_outside_option) {
    outside_alt_mapping <- data.table::data.table(alt_int=0L, N_OBS = N, N_CHOICES = sum(choice_idx == 0L))
    outside_alt_mapping[[alt_col]] <- outside_opt_label
    alt_mapping <- list(outside_alt_mapping, alt_mapping) |>
      data.table::rbindlist(use.names = TRUE, fill = TRUE)
    data.table::setcolorder(alt_mapping, c("alt_int", alt_col, "N_OBS", "N_CHOICES"))
  }

  alt_mapping[, `:=`(
    TAKE_RATE = N_CHOICES / N_OBS,
    MKT_SHARE = N_CHOICES / sum(N_CHOICES)
  )]

  ## Final validity checks -----------------------------------------------------
  stopifnot(
    length(alt_idx)    == nrow(X),
    length(choice_idx) == N,
    length(M)          == N,
    length(weights)    == N,
    is.null(Ti) || sum(Ti) == N,
    all(is.finite(X)),
    all(is.finite(W))
  )

  ## Return output -------------------------------------------------------------
  out <- structure(
    list(
      X           = X,
      W           = W,
      alt_idx     = alt_idx,
      choice_idx  = as.integer(choice_idx),
      M           = M,
      N           = N,
      weights     = weights,
      cluster     = cluster,
      situation_ids = ids,
      Ti          = Ti,
      person_ids  = person_ids,
      include_outside_option = include_outside_option,
      rc_correlation = rc_correlation,
      alt_mapping = alt_mapping[],
      dropped_cols = if(exists("dropped_vars")) dropped_vars else NULL,
      data_spec = list(
        id_col = id_col,
        alt_col = alt_col,
        choice_col = choice_col,
        covariate_cols = covariate_cols,
        random_var_cols = random_var_cols,
        outside_opt_label = outside_opt_label,
        person_col = person_col
      )
    ),
    class = "choicer_data_mxl"
  )
  if (!is.null(cs_provenance)) {
    attr(out, "choice_sampling") <- cs_provenance
  }
  out
}

#' Halton draws for mixed logit
#'
#' Create halton normal draws in appropriate format for mixed logit estimation
#'
#' Draw unit \eqn{i} receives the \eqn{S} consecutive points
#' \eqn{(i - 1)S + 1, \ldots, iS} of the \eqn{K_w}-dimensional Halton sequence
#' of \code{randtoolbox::halton()}, mapped to standard normals. The cube
#' takes \eqn{8 K_w S N} bytes and is filled a block of units at a time, so it
#' never holds a second copy of the sequence. Store mode supports at most
#' \eqn{2^{31} - 1} points (\eqn{S \times N}), the largest starting index
#' \code{halton(start = )} accepts; \code{draws = "generate"} in
#' \code{\link{run_mxlogit}} has no such limit and never materializes the cube.
#' \code{predict()}, \code{logsum()}, \code{consumer_surplus()},
#' \code{elasticities()}, \code{diversion_ratios()}, \code{blp()} and
#' \code{gof()} on a store-mode fit with up to 128 random coefficients do not
#' build their cube (one block per choice situation) with this function: they
#' form the same draws on the fly, from the same uniforms mapped by R's
#' \code{qnorm()} (bit for bit where randtoolbox and choicer are compiled
#' with the same floating-point contraction, as on CRAN's platforms).
#'
#' @param S Number of draws per draw unit
#' @param N number of draw units: choice situations, or decision makers for a
#'   panel fit (\code{person_col} in \code{\link{run_mxlogit}})
#' @param K_w dimension of random coefficients (number of columns in W matrix)
#' @returns K_w x S x N array with halton standard normal draws
#' @examples
#' draws <- get_halton_normals(S = 50, N = 10, K_w = 2)
#' dim(draws)  # 2 x 50 x 10
#' @importFrom randtoolbox halton
#' @export
get_halton_normals <- function(S, N, K_w) {
  if (!.whole_number(S)) stop("`S` must be a single positive whole number.")
  if (!.whole_number(N)) stop("`N` must be a single positive whole number.")
  if (!.whole_number(K_w)) stop("`K_w` must be a single positive whole number.")
  if (!.halton_points_ok(S, N)) stop(.halton_points_msg(as.numeric(S) * N))
  .halton_cube(S, N, K_w)
}

#' A single positive whole number, as get_halton_normals() takes S, N, K_w
#' @noRd
.whole_number <- function(x) {
  is.numeric(x) && length(x) == 1L && is.finite(x) && x >= 1 && x == round(x)
}

#' Whether S * N points of the Halton sequence are within randtoolbox's reach
#'
#' `randtoolbox::halton(start = )` takes an int start, so store-mode draws
#' reach at most 2^31 - 1 points.
#' @noRd
.halton_points_ok <- function(S, N) as.numeric(S) * N <= .Machine$integer.max

#' The error for store-mode draws past randtoolbox's largest start
#' @param n_points S * N.
#' @noRd
.halton_points_msg <- function(n_points) {
  paste0("A store-mode draw cube needs S * N = ",
         format(n_points, big.mark = ",", scientific = FALSE),
         " points of the Halton sequence, more than 2^31 - 1, the largest ",
         "starting index randtoolbox::halton(start = ) accepts. Refit with ",
         "run_mxlogit(draws = \"generate\").")
}

#' Points start, ..., start + n - 1 of the Halton sequence, as normals
#'
#' `randtoolbox::halton()`'s n x K_w matrix (a vector when K_w = 1).
#' `halton(start = )` (randtoolbox >= 1.31.0) reproduces any slice of the
#' sequence bit for bit, so every way of cutting it into calls gives the same
#' draws.
#' @noRd
.halton_points <- function(start, n, K_w) {
  randtoolbox::halton(n = n, dim = K_w, normal = TRUE, start = start)
}

#' Fill the store-mode Halton cube a block of draw units at a time
#'
#' Unit i's draws are points (i - 1) * S + 1, ..., i * S of the sequence.
#' `randtoolbox::halton(start = )` (randtoolbox >= 1.31.0) reproduces any slice
#' of the full sequence bit for bit, so the blocks give exactly the draws of a
#' single call over all S * N points. They avoid that call's copies of the
#' whole sequence and its 32-bit offsets into the output, which overflow past
#' 2^31 - 1 values.
#'
#' @param S,N,K_w As in [get_halton_normals()], already validated.
#' @param block Target number of values generated per `halton()` call; a
#'   block holds at least one unit.
#' @returns K_w x S x N array.
#' @noRd
.halton_cube <- function(S, N, K_w, block = 2^22) {
  eta <- array(0, dim = c(K_w, S, N))
  units_per_block <- max(1, floor(block / (as.numeric(S) * K_w)))
  i0 <- 1
  while (i0 <= N) {
    i1 <- min(N, i0 + units_per_block - 1)
    # (S * (i1 - i0 + 1)) x K_w, or a vector when K_w = 1; its transpose is
    # eta[, , i0:i1] in memory order either way.
    h <- .halton_points((i0 - 1) * S + 1, S * (i1 - i0 + 1), K_w)
    eta[, , i0:i1] <- t(h)
    i0 <- i1 + 1
  }
  eta
}


#' Keep the likelihood kernel's overflow sentinel above the optimizer's path
#'
#' `mxl_loglik_gradient_parallel()` returns `objective = 1e10` with a zero
#' gradient and `overflow = TRUE` where the utilities overflow, so that a
#' line search backtracks. A finite objective can also equal 1e10, so only
#' the explicit flag identifies the sentinel.
#' Its objective is otherwise finite however poor the fit, and unbounded, so
#' after a start above 1e10 (a badly scaled warm start at population scale)
#' the fixed sentinel would look like an improvement. Wraps `eval_f` to report
#' the sentinel as ten times the largest objective seen so far, and never
#' below 1e10; objectives are negated log-likelihoods, hence nonnegative.
#'
#' @param eval_f Function of theta returning `list(objective, gradient, overflow)`.
#' @param f_max The largest objective known before the first call (the
#'   profile likelihood's solved start), 0 by default.
#' @returns The wrapped function.
#' @noRd
.lift_sentinel <- function(eval_f, f_max = 0) {
  force(f_max)
  function(theta) {
    res <- eval_f(theta)
    if (isTRUE(res$overflow)) {
      res$objective <- max(1e10, min(10 * f_max, .Machine$double.xmax))
    } else if (is.finite(res$objective)) {
      f_max <<- max(f_max, res$objective)
    }
    res$overflow <- NULL
    res
  }
}


#' Resolve draw parameters for post-estimation regeneration sites
#'
#' When the fitted object used generate mode, returns an empty placeholder cube
#' plus the three gen_* integers. When in store mode, materialises the Halton
#' cube from the stored metadata. The estimation-type sites (Hessian, scores,
#' conditional tastes) use it; the prediction sites use .mxl_pred_draws(),
#' whose kernels form the store-mode draws on the fly.
#'
#' @param draws_info List from a fitted choicer_mxl object.
#' @param N Number of draw blocks: the number of likelihood units,
#'   `length(.unit_first(d))`, which is the number of decision makers for a
#'   panel fit.
#' @noRd
.mxl_gen_params <- function(draws_info, N = draws_info$N) {
  mode <- draws_info$mode %||% "store"
  if (mode == "generate") {
    list(
      eta_draws    = array(0, dim = c(draws_info$K_w, 0L, 0L)),
      gen_seed     = as.integer(draws_info$seed),
      # "owen" is the legacy serialized label for the same position-wise
      # permutation. Read it silently so old fitted objects keep working.
      gen_scramble = if (draws_info$scramble %in% c("permuted", "owen")) 1L else 0L,
      gen_S        = as.integer(draws_info$S)
    )
  } else {
    list(
      eta_draws    = get_halton_normals(draws_info$S, N, draws_info$K_w),
      gen_seed     = -1L,
      gen_scramble = 1L,
      gen_S        = 0L
    )
  }
}

#' Bytes as a size in binary units (1 GiB = 2^30 bytes), for messages
#' @noRd
.format_bytes <- function(bytes) {
  format(structure(bytes, class = "object_size"), units = "auto",
         standard = "IEC", digits = 2L)
}

#' Warn before a store-mode fit builds a draw cube above the budget
#'
#' The cube is held for the whole fit, and the estimation-type post-estimation
#' functions rebuild it; generate mode stores none. Called only where
#' get_halton_normals() will build it: above 2^31 - 1 points it stops instead.
#'
#' @param S,n_units,K_w Draws per unit, likelihood units, random coefficients.
#' @param panel Whether the units are decision makers (`person_col`).
#' @noRd
.warn_store_cube <- function(S, n_units, K_w, panel) {
  # get_halton_normals() reports invalid arguments, and stops past 2^31 - 1
  # points
  if (!.whole_number(S) || !.whole_number(n_units) || !.whole_number(K_w) ||
      !.halton_points_ok(S, n_units)) {
    return(invisible())
  }
  bytes <- 8 * as.numeric(K_w) * S * n_units
  if (bytes <= .mxl_cube_budget()) return(invisible())
  big <- function(x) format(x, big.mark = ",", scientific = FALSE)
  # immediate.: the fit that follows can take hours, or run out of memory
  # (unless the user ignores warnings, options(warn = -1), which it would
  # override)
  warning(sprintf(paste0(
    "draws = \"store\" will hold %s of Halton draws (%s %s x %s %s x %d ",
    "random %s) while the model is fitted, and vcov(type = ), wesml_vcov() ",
    "and conditional_tastes() rebuild them whole. draws = \"generate\" forms ",
    "the draws on the fly and stores none (with scramble = \"none\", the same ",
    "Halton points)."),
    .format_bytes(bytes), big(n_units),
    if (panel) "decision makers" else "choice situations", big(S),
    ngettext(S, "draw", "draws"), K_w,
    ngettext(K_w, "coefficient", "coefficients")),
    call. = FALSE, immediate. = !isTRUE(getOption("warn") < 0))
}

#' Byte budget of a store-mode estimation cube
#'
#' `run_mxlogit()` warns before it builds an estimation cube larger than it
#' (`.warn_store_cube()`). 1 GiB; a function so that tests can lower it with
#' `local_mocked_bindings()`.
#' @noRd
.mxl_cube_budget <- function() 2^30

#' The `gen_scramble` code of store-mode draws formed on the fly
#'
#' The prediction kernels map the generator's identity-permutation uniforms,
#' which are `randtoolbox::halton()`'s points (bit for bit when both are
#' compiled with the same floating-point contraction; see
#' `HaltonGen::fill_uniforms()`), to normals with R's `qnorm()`, as
#' [get_halton_normals()] does (`MXL_STORE_POINTS` in src/mxlogit.cpp).
#' @noRd
.MXL_STORE_POINTS <- 2L

#' Dimensions the on-the-fly generator supports (its primes table,
#' `HALTON_N_PRIMES` in src/halton.h)
#' @noRd
.HALTON_MAX_DIM <- 128L

#' Draw arguments of an MXL prediction kernel
#'
#' Generate mode: an empty cube and the generator's seed, scramble and S.
#' Store mode: an empty cube and the generator's arguments for the store-mode
#' points (`.MXL_STORE_POINTS`): the kernel forms each choice situation's
#' block of S draws on the fly, those of the K_w x S x N cube
#' [get_halton_normals()] builds (see `.MXL_STORE_POINTS` for the
#' bit-for-bit condition), which it never holds. With more random
#' coefficients than the generator's primes table covers, the cube itself.
#' S and K_w are checked as get_halton_normals() checks them.
#'
#' @param draws_info The fit's draw metadata.
#' @param N Number of choice situations.
#' @param what The calling function, named in the error when `draws_info` is
#'   missing.
#' @returns List with `eta_draws`, `gen_seed`, `gen_scramble` and `gen_S`.
#' @noRd
.mxl_pred_draws <- function(draws_info, N, what) {
  if (is.null(draws_info)) {
    stop(what, " requires draws_info from a fitted MXL model.", call. = FALSE)
  }
  if ((draws_info$mode %||% "store") == "generate") {
    return(.mxl_gen_params(draws_info))
  }
  K_w <- draws_info$K_w
  S <- draws_info$S
  if (!.whole_number(S)) {
    stop("`S` must be a single positive whole number.", call. = FALSE)
  }
  if (!.whole_number(K_w)) {
    stop("`K_w` must be a single positive whole number.", call. = FALSE)
  }
  if (K_w > .HALTON_MAX_DIM) {
    return(list(eta_draws = get_halton_normals(S, N, K_w), gen_seed = -1L,
                gen_scramble = 1L, gen_S = 0L))
  }
  list(eta_draws = array(0, dim = c(K_w, 0L, 0L)), gen_seed = 0L,
       gen_scramble = .MXL_STORE_POINTS, gen_S = as.integer(S))
}
