
# Expected logsum and consumer surplus (Train 2009, Ch. 3)

# --- Internal helpers ---------------------------------------------------------

#' Resolve the data inputs for logsum/consumer surplus
#'
#' Mirrors the data resolution of the predict() methods: stored fit data when
#' `newdata` is NULL (with the same keep_data error), or
#' `resolve_predict_newdata()` otherwise.
#'
#' @param object A choicer_fit object.
#' @param newdata NULL, a long-format data.frame, or a list (X, alt_idx, M, ...).
#' @param needs_draws Whether MXL draw metadata (`draws_info`) is required.
#' @param weights Optional prediction weights for the newdata path (one per
#'   choice situation, as in `predict()`); ignored when `newdata` is NULL.
#' @returns List with `X`, `W`, `alt_idx`, `M`, `N`, `weights`.
#' @noRd
.surplus_resolve_data <- function(object, newdata, needs_draws = FALSE,
                                  weights = NULL) {
  if (is.null(newdata)) {
    # [["data"]]: $data would partial-match data_spec if the element is absent
    d <- object[["data"]]
    if (is.null(d) || (needs_draws && is.null(object$draws_info))) {
      if (needs_draws) {
        stop("Prediction requires stored data and draws. ",
             "Refit with keep_data = TRUE.")
      }
      stop("Prediction requires stored data. Refit with keep_data = TRUE.")
    }
    d$N <- length(d$M)
    d
  } else {
    if (needs_draws && is.null(object$draws_info)) {
      stop("Prediction with newdata requires stored draw metadata ",
           "('draws_info'). Refit to enable newdata prediction.")
    }
    resolve_predict_newdata(object, newdata, weights = weights)
  }
}

#' Blockwise stable log-sum-exp over choice situations
#'
#' Computes `log(sum_j exp(V_ij))` per id block (plus the outside option's
#' exp(0) = 1 term when present) with max-subtraction for stability.
#'
#' @param V Stacked utility vector (one entry per inside-alternative row).
#' @param M Integer vector with the number of inside alternatives per id.
#' @param include_outside_option Whether to append the outside option's V = 0.
#' @returns Numeric vector of length `length(M)`.
#' @noRd
.blockwise_logsumexp <- function(V, M, include_outside_option) {
  ends <- cumsum(M)
  starts <- ends - M + 1L
  vapply(seq_along(M), function(i) {
    v <- V[starts[i]:ends[i]]
    if (include_outside_option) v <- c(v, 0)
    v_max <- max(v)
    v_max + log(sum(exp(v - v_max)))
  }, numeric(1))
}

#' Reconstruct the full lambda vector for a nested logit fit
#'
#' Singleton nests have lambda fixed to 1; non-singleton nests receive the
#' estimated lambdas, scattered by ascending nest index (matching the C++
#' parsing in `nl_parse_theta()`; same logic as `blp.choicer_nl`).
#'
#' @param object A choicer_nl object.
#' @param nest_idx Per-alternative (length J) nest index vector.
#' @returns Numeric vector of length `max(nest_idx)`.
#' @noRd
.nl_lambda_full <- function(object, nest_idx) {
  n_nests <- max(nest_idx)
  lambda_full <- rep(1, n_nests)
  nest_counts <- tabulate(nest_idx, n_nests)
  non_singleton <- which(nest_counts > 1)
  lambda_full[non_singleton] <- object$coefficients[object$param_map$lambda]
  lambda_full
}

# --- logsum: generic and methods ----------------------------------------------

#' Expected logsum (inclusive value) per choice situation
#'
#' Computes the expected maximum utility ("logsum" or inclusive value) for each
#' choice situation, up to the additive constant of the extreme-value error:
#' \itemize{
#'   \item MNL: \eqn{\log \sum_j \exp(V_{ij})}.
#'   \item MXL: \eqn{E_\beta[\log \sum_j \exp(V_{ij}(\beta))]}, simulated by
#'     averaging the log-sum-exp \emph{across} the deterministic Halton draws
#'     (regenerated from \code{object$draws_info}). Taking the log-sum-exp of
#'     draw-averaged utilities would understate the expectation (Jensen's
#'     inequality), so a dedicated kernel (\code{\link{mxl_logsum}}) is used.
#'   \item NL: \eqn{\log \sum_b \exp(\lambda_b I_{ib})} with the nest inclusive
#'     value \eqn{I_{ib} = \log \sum_{j \in b} \exp(V_{ij}/\lambda_b)}
#'     (singleton nests have \eqn{\lambda_b = 1}).
#' }
#' When the model includes an outside option, its normalized utility
#' \eqn{V = 0} contributes an \eqn{\exp(0)} term to the sum.
#'
#' Logsum \emph{levels} depend on the ASC normalization (and, more generally,
#' on any additive utility normalization), so only logsum \emph{differences}
#' between scenarios (e.g. via \code{newdata}) are meaningful.
#'
#' @param object A fitted model object (\code{choicer_mnl}, \code{choicer_mxl},
#'   or \code{choicer_nl}).
#' @param newdata Optional counterfactual data: a data.frame in the fit-time
#'   long format or a list with \code{X}, \code{alt_idx}, \code{M} (plus
#'   \code{W} for MXL), as in \code{predict()}. When \code{NULL} (default),
#'   the data stored at fit time is used (requires \code{keep_data = TRUE}).
#' @param ... Additional arguments passed to methods.
#' @returns Numeric vector with one logsum per choice situation. With a
#'   data.frame \code{newdata}, choice situations are ordered by id (as in
#'   \code{predict()}).
#' @seealso \code{\link{consumer_surplus}}
#' @examples
#' \donttest{
#' library(data.table)
#' sim <- simulate_mnl_data(N = 500, J = 3, beta = c(0.8, -0.6), seed = 1,
#'                          outside_option = FALSE, vary_choice_set = FALSE)
#' fit <- run_mnlogit(sim$data, "id", "alt", "choice", c("x1", "x2"))
#' head(logsum(fit))
#' }
#' @export
logsum <- function(object, newdata = NULL, ...) {
  UseMethod("logsum")
}

#' @rdname logsum
#' @export
logsum.choicer_mnl <- function(object, newdata = NULL, ...) {
  d <- .surplus_resolve_data(object, newdata)
  .logsum_mnl_core(object, d)
}

#' MNL logsum from resolved kernel inputs
#' @noRd
.logsum_mnl_core <- function(object, d) {
  pred <- mnl_predict(
    theta = object$coefficients,
    X = d$X,
    alt_idx = d$alt_idx,
    M = d$M,
    use_asc = object$use_asc,
    include_outside_option = object$include_outside_option
  )
  .blockwise_logsumexp(as.numeric(pred$utility), d$M,
                       object$include_outside_option)
}

#' @rdname logsum
#' @export
logsum.choicer_mxl <- function(object, newdata = NULL, ...) {
  d <- .surplus_resolve_data(object, newdata, needs_draws = TRUE)
  .logsum_mxl_core(object, d)
}

#' MXL logsum from resolved kernel inputs (regenerates the Halton draws)
#' @noRd
.logsum_mxl_core <- function(object, d) {
  eta_draws <- get_halton_normals(
    S   = object$draws_info$S,
    N   = d$N,
    K_w = object$draws_info$K_w
  )

  as.numeric(mxl_logsum(
    theta                  = object$coefficients,
    X                      = d$X,
    W                      = d$W,
    alt_idx                = d$alt_idx,
    M                      = d$M,
    eta_draws              = eta_draws,
    rc_dist                = object$rc_dist,
    rc_correlation         = object$rc_correlation,
    rc_mean                = object$rc_mean,
    use_asc                = object$use_asc,
    include_outside_option = object$include_outside_option
  ))
}

#' Resolve the authoritative alternative-to-nest mapping for an NL fit
#' @noRd
.nl_resolve_nest_idx <- function(object) {
  # The per-alternative (length J) nest mapping is authoritative from the fit
  # since it indexes the estimated lambda parameters (see predict.choicer_nl).
  nest_idx <- object$nest_idx %||% object[["data"]]$nest_idx
  if (is.null(nest_idx)) {
    stop("Cannot resolve the alternative-to-nest mapping: this fit stores ",
         "no 'nest_idx'. Refit to enable newdata prediction.")
  }
  nest_idx
}

#' @rdname logsum
#' @export
logsum.choicer_nl <- function(object, newdata = NULL, ...) {
  nest_idx <- .nl_resolve_nest_idx(object)
  d <- .surplus_resolve_data(object, newdata)
  .logsum_nl_core(object, d, nest_idx)
}

#' NL logsum from resolved kernel inputs
#' @noRd
.logsum_nl_core <- function(object, d, nest_idx) {
  pred <- nl_predict(
    theta = object$coefficients,
    X = d$X,
    alt_idx = d$alt_idx,
    M = d$M,
    nest_idx = nest_idx,
    use_asc = object$use_asc,
    include_outside_option = object$include_outside_option
  )
  V <- as.numeric(pred$utility)
  lambda_full <- .nl_lambda_full(object, nest_idx)
  nest_of_row <- nest_idx[d$alt_idx]

  # logsum_i = log( sum_b exp(lambda_b * I_ib) [+ exp(0)] ) with
  # I_ib = log sum_{j in b} exp(V_ij / lambda_b); mirrors nl_individual_probs.
  ends <- cumsum(d$M)
  starts <- ends - d$M + 1L
  vapply(seq_along(d$M), function(i) {
    rows <- starts[i]:ends[i]
    v_i <- V[rows]
    b_i <- nest_of_row[rows]
    terms <- vapply(unique(b_i), function(b) {
      v_b <- v_i[b_i == b] / lambda_full[b]
      v_max <- max(v_b)
      lambda_full[b] * (v_max + log(sum(exp(v_b - v_max))))
    }, numeric(1))
    if (object$include_outside_option) terms <- c(terms, 0)
    t_max <- max(terms)
    t_max + log(sum(exp(terms - t_max)))
  }, numeric(1))
}

# --- consumer surplus: internal -----------------------------------------------

#' Delta-method SE of the weighted mean consumer surplus (MNL)
#'
#' For \eqn{m(\theta) = \sum_i w_i CS_i / \sum_i w_i} with
#' \eqn{CS_i = logsum_i / (-\alpha)}, the gradient blocks are
#' \itemize{
#'   \item \eqn{\partial logsum_i / \partial \theta_r =
#'     \sum_j P_{ij} \, \partial V_{ij} / \partial \theta_r}, with
#'     \eqn{\partial V_{ij}/\partial \beta_r = x_{ijr}} and
#'     \eqn{\partial V_{ij}/\partial ASC_a = 1\{j\ has\ ASC\ a\}} (the outside
#'     option row has \eqn{V = 0} and contributes nothing);
#'   \item \eqn{\partial CS_i/\partial \theta_r =
#'     \partial logsum_i/\partial \theta_r / (-\alpha)} for non-price
#'     parameters;
#'   \item for the price coefficient \eqn{\alpha} itself,
#'     \eqn{\partial CS_i/\partial \alpha = logsum_i / \alpha^2 +
#'     (1/(-\alpha)) \sum_j P_{ij} x_{ij,price}} (it enters both the logsum
#'     and the \eqn{1/(-\alpha)} factor).
#' }
#' The SE is \eqn{\sqrt{G' V G}} with \eqn{G} the weighted average of the
#' per-i gradients and \eqn{V} the coefficient variance.
#'
#' @param object A choicer_mnl object with `vcov` populated.
#' @param d Resolved data list (X, alt_idx, M, weights).
#' @param ls Logsum vector (length N).
#' @param alpha Price coefficient.
#' @param price_idx Coefficient index of the price variable.
#' @returns Scalar SE (NA when the variance is unavailable or negative).
#' @noRd
.cs_mnl_delta_se <- function(object, d, ls, alpha, price_idx) {
  V_mat <- object$vcov
  if (is.null(V_mat)) return(NA_real_)

  pred <- mnl_predict(
    theta = object$coefficients,
    X = d$X,
    alt_idx = d$alt_idx,
    M = d$M,
    use_asc = object$use_asc,
    include_outside_option = object$include_outside_option
  )
  P <- as.numeric(pred$choice_prob)

  w <- d$weights
  w_total <- sum(w)
  w_row <- rep(w, d$M)
  wp <- w_row * P
  pm <- object$param_map

  G <- numeric(length(object$coefficients))

  # Beta block: G_r = (1/(W * -alpha)) * sum_rows w_i P_ij x_ijr
  G[pm$beta] <- colSums(d$X * wp) / (w_total * (-alpha))

  # ASC block: with an outside option every inside alt has a free ASC
  # (theta slot pm$asc[alt]); without, the first inside alt is normalized
  # to 0 (slot pm$asc[alt - 1] for alt >= 2). See src/mnlogit.cpp.
  if (object$use_asc && !is.null(pm$asc)) {
    slot <- if (object$include_outside_option) {
      pm$asc[d$alt_idx]
    } else {
      # alt 1 is normalized: index it explicitly to NA (subsetting with
      # alt_idx - 1 would drop the zero indices and misalign the rows)
      s <- rep(NA_integer_, length(d$alt_idx))
      inside <- d$alt_idx >= 2L
      s[inside] <- pm$asc[d$alt_idx[inside] - 1L]
      s
    }
    ok <- !is.na(slot)
    asc_sums <- tapply(wp[ok], slot[ok], sum)
    G[as.integer(names(asc_sums))] <- asc_sums / (w_total * (-alpha))
  }

  # Price coefficient appears in both the logsum and the 1/(-alpha) factor:
  # add the weighted mean of logsum_i / alpha^2 to its (already filled) slot.
  G[price_idx] <- G[price_idx] + sum(w * ls) / (w_total * alpha^2)

  var_m <- as.numeric(t(G) %*% V_mat %*% G)
  if (is.finite(var_m) && var_m >= 0) sqrt(var_m) else NA_real_
}

#' Assemble a choicer_cs object
#' @noRd
.cs_build <- function(cs, weights, se, price_var, level) {
  if (!is.numeric(level) || length(level) != 1L || is.na(level) ||
      level <= 0 || level >= 1) {
    stop("'level' must be a single number strictly between 0 and 1.")
  }
  mean_cs <- sum(weights * cs) / sum(weights)
  q <- stats::qnorm(1 - (1 - level) / 2)
  ci <- if (is.finite(se)) {
    c(mean_cs - q * se, mean_cs + q * se)
  } else {
    c(NA_real_, NA_real_)
  }
  structure(
    list(
      cs = cs,
      mean_cs = mean_cs,
      se_mean_cs = se,
      ci = ci,
      price_var = price_var,
      level = level,
      n = length(cs)
    ),
    class = "choicer_cs"
  )
}

# --- consumer surplus: generic and methods -------------------------------------

#' Expected consumer surplus
#'
#' Computes the expected consumer surplus per choice situation (Train 2009,
#' Ch. 3):
#' \deqn{E[CS_i] = \frac{logsum_i}{-\alpha},}
#' where \eqn{logsum_i} is the expected maximum utility (see
#' \code{\link{logsum}}) and \eqn{\alpha} is the (fixed) price coefficient,
#' so that \eqn{-\alpha} is the marginal utility of income. The formula
#' assumes \emph{no income effects}: utility is linear in price, and the
#' marginal utility of income is constant across the price changes considered.
#'
#' Consumer surplus \emph{levels} inherit the additive utility normalization
#' (in particular the ASC normalization), so the level is only defined up to a
#' constant; \emph{differences} in CS between scenarios — e.g.
#' \code{consumer_surplus(fit, "price", newdata = scenario)} minus the baseline
#' — are the economically meaningful quantity.
#'
#' For MNL fits, a delta-method standard error of the weighted mean CS is
#' reported (weights are the stored fit weights, or the resolved
#' \code{newdata} weights). For MXL and NL fits only point estimates are
#' returned (\code{se_mean_cs = NA}): the delta method for the simulated MXL
#' logsum and the nested logsum is deferred; simulation-based intervals
#' (Krinsky-Robb: resample coefficients from their asymptotic distribution
#' and recompute the mean CS) are a practical alternative.
#'
#' The price variable must have a \emph{fixed} coefficient. For mixed logit a
#' random price coefficient is rejected (as in \code{\link{wtp}}): with a
#' random denominator \eqn{1/(-\alpha)} generally has no finite moments.
#'
#' @param object A fitted model object (\code{choicer_mnl}, \code{choicer_mxl},
#'   or \code{choicer_nl}).
#' @param price_var Name of the price variable. Must be a fixed-coefficient
#'   variable (a column of the design matrix \code{X}).
#' @param newdata Optional counterfactual data (data.frame or list), as in
#'   \code{\link{logsum}} and \code{predict()}. When \code{NULL} (default),
#'   the data stored at fit time is used (requires \code{keep_data = TRUE}).
#' @param level Confidence level for the normal-approximation interval around
#'   the mean CS (MNL only). Default 0.95.
#' @param weights Optional numeric vector with one weight per choice situation,
#'   used for the mean CS (and its SE), as in `predict()`: for a data.frame
#'   `newdata`, one weight per id in order of first appearance. Defaults to
#'   equal weights. Ignored when `newdata` is `NULL` (the stored fit weights
#'   apply).
#' @param ... Additional arguments passed to methods.
#' @returns A \code{choicer_cs} object: a list with \code{cs} (per-choice-
#'   situation surplus, length N), \code{mean_cs} (weighted mean),
#'   \code{se_mean_cs} (delta-method SE; NA for MXL/NL or when the
#'   variance-covariance matrix is unavailable), \code{ci} (confidence
#'   interval for the mean), \code{price_var}, \code{level}, and \code{n}.
#' @references Train, K. (2009). \emph{Discrete Choice Methods with
#'   Simulation}, 2nd ed., Ch. 3. Cambridge University Press.
#' @seealso \code{\link{logsum}}, \code{\link{wtp}}
#' @examples
#' \donttest{
#' library(data.table)
#' sim <- simulate_mnl_data(N = 1000, J = 3, beta = c(0.8, -0.6), seed = 123,
#'                          outside_option = FALSE, vary_choice_set = FALSE)
#' fit <- run_mnlogit(sim$data, "id", "alt", "choice", c("x1", "x2"))
#'
#' # treat x2 as the price variable
#' cs0 <- consumer_surplus(fit, price_var = "x2")
#' cs0
#'
#' # Change in consumer surplus from a price increase on alternative 2:
#' # levels depend on the ASC normalization, differences do not.
#' dt_cf <- copy(sim$data)[alt == 2, x2 := x2 + 0.5]
#' cs1 <- consumer_surplus(fit, price_var = "x2", newdata = dt_cf)
#' delta_cs <- cs1$mean_cs - cs0$mean_cs
#' delta_cs  # negative: the price increase lowers expected surplus
#' }
#' @export
consumer_surplus <- function(object, price_var, newdata = NULL,
                             level = 0.95, weights = NULL, ...) {
  UseMethod("consumer_surplus")
}

#' @rdname consumer_surplus
#' @export
consumer_surplus.choicer_mnl <- function(object, price_var, newdata = NULL,
                                         level = 0.95, weights = NULL, ...) {
  price_idx <- .wtp_price_index(object, price_var)
  alpha <- unname(object$coefficients[price_idx])

  d <- .surplus_resolve_data(object, newdata, weights = weights)
  ls <- .logsum_mnl_core(object, d)
  cs <- ls / (-alpha)

  object <- ensure_vcov(object)
  se <- .cs_mnl_delta_se(object, d, ls, alpha, price_idx)

  .cs_build(cs, d$weights, se, price_var, level)
}

#' @rdname consumer_surplus
#' @export
consumer_surplus.choicer_mxl <- function(object, price_var, newdata = NULL,
                                         level = 0.95, weights = NULL, ...) {
  w_names <- names(object$sW) %||% colnames(object[["data"]]$W) %||% character(0)
  if (is.character(price_var) && length(price_var) == 1L &&
      price_var %in% w_names) {
    stop("Random price coefficients are not supported: 1/(-alpha) with a ",
         "random alpha generally has no finite moments. Use a fixed price ",
         "coefficient (a 'covariate_cols' variable) or estimate the model ",
         "in WTP space.")
  }
  price_idx <- .wtp_price_index(object, price_var)
  alpha <- unname(object$coefficients[price_idx])

  d <- .surplus_resolve_data(object, newdata, needs_draws = TRUE,
                             weights = weights)
  cs <- .logsum_mxl_core(object, d) / (-alpha)

  # Delta-method SE for the simulated logsum is deferred; use Krinsky-Robb
  # resampling of the coefficients for intervals.
  .cs_build(cs, d$weights, NA_real_, price_var, level)
}

#' @rdname consumer_surplus
#' @export
consumer_surplus.choicer_nl <- function(object, price_var, newdata = NULL,
                                        level = 0.95, weights = NULL, ...) {
  nest_idx <- .nl_resolve_nest_idx(object)
  price_idx <- .wtp_price_index(object, price_var)
  alpha <- unname(object$coefficients[price_idx])

  d <- .surplus_resolve_data(object, newdata, weights = weights)
  cs <- .logsum_nl_core(object, d, nest_idx) / (-alpha)

  # Delta-method SE for the nested logsum is deferred; use Krinsky-Robb
  # resampling of the coefficients for intervals.
  .cs_build(cs, d$weights, NA_real_, price_var, level)
}

# --- print --------------------------------------------------------------------

#' Print a consumer surplus summary
#'
#' @param x A \code{choicer_cs} object.
#' @param digits Number of significant digits to print.
#' @param ... Additional arguments (ignored).
#' @returns The object invisibly.
#' @export
print.choicer_cs <- function(x, digits = 4, ...) {
  cat(sprintf("Consumer surplus, price variable: '%s'\n", x$price_var))
  cat("  Mean CS:", format(x$mean_cs, digits = digits), "\n")
  if (is.finite(x$se_mean_cs)) {
    cat("  SE (delta method):", format(x$se_mean_cs, digits = digits), "\n")
    cat(sprintf("  %g%% CI: [%s, %s]\n", 100 * x$level,
                format(x$ci[1], digits = digits),
                format(x$ci[2], digits = digits)))
  } else {
    cat("  SE: NA (delta-method SE is available for MNL fits with a ",
        "variance-covariance matrix; use Krinsky-Robb resampling otherwise)\n",
        sep = "")
  }
  cat("  N:", x$n, "\n")
  cat("Note: CS levels depend on the utility (ASC) normalization; ",
      "differences between scenarios are the meaningful quantity.\n",
      sep = "")
  invisible(x)
}
