
# Goodness-of-fit measures for fitted choice models

# --- Internal helpers ---------------------------------------------------------

#' Closed-form market-shares null log-likelihood
#'
#' \eqn{LL_0 = \sum_j N\_CHOICES_j \log(MKT\_SHARE_j)} from
#' \code{object$alt_mapping}. This is the maximized log-likelihood of an
#' ASC-only model, which has a closed form only when every individual faces
#' the same choice set (every inside alternative present in every choice
#' situation, not merely equal set sizes) and weights are uniform.
#' When \code{include_outside_option = TRUE}, \code{alt_mapping} contains an
#' outside-option row (alt_int = 0) whose term is included. Alternatives that
#' are never chosen (\code{N_CHOICES = 0}) contribute 0 and are dropped to
#' avoid \code{0 * log(0)}.
#'
#' @param object A choicer_fit object with stored data.
#' @returns Scalar null log-likelihood.
#' @noRd
.gof_market_shares_ll0 <- function(object) {
  d <- object$data
  am <- object$alt_mapping
  # Identical composition, not just identical size: every inside alternative
  # must appear in every choice situation (equal-size heterogeneous sets do
  # not admit the closed form either).
  inside_n_obs <- am$N_OBS[am$alt_int > 0L]
  balanced <- length(unique(d$M)) == 1L && all(inside_n_obs == length(d$M))
  uniform_w <- length(unique(d$weights)) == 1L
  if (!balanced || !uniform_w) {
    stop("The 'market_shares' null log-likelihood has a closed form only for ",
         "balanced designs (every alternative available in every choice ",
         "situation) and uniform weights. Refit an ASC-only model to ",
         "obtain the market-shares null log-likelihood for this design.")
  }
  n_choices <- am$N_CHOICES
  shares <- am$MKT_SHARE
  keep <- n_choices > 0
  # Uniform weights scale every individual's contribution by the same factor.
  d$weights[1] * sum(n_choices[keep] * log(shares[keep]))
}

#' Weighted hit rate from in-sample predicted probabilities
#'
#' The hit rate is the (weighted) share of individuals whose observed choice
#' coincides with the highest predicted probability. With an outside option,
#' the outside good competes too: its probability is taken from
#' \code{choice_prob_outside} when the predictor returns it (MXL), and from
#' \eqn{1 - \sum_j p_{ij}} otherwise (MNL/NL); a win by the outside good is a
#' predicted index of 0, matching \code{choice_idx = 0}.
#'
#' @param object A choicer_fit object with stored data.
#' @returns Scalar weighted hit rate in \[0, 1\].
#' @noRd
.gof_hit_rate <- function(object) {
  d <- object$data
  preds <- predict(object, type = "probabilities")
  p <- preds$choice_prob

  M <- d$M
  ends <- cumsum(M)
  starts <- ends - M + 1L
  ioo <- isTRUE(object$include_outside_option)
  p_out <- preds$choice_prob_outside  # NULL unless returned (MXL with outside)

  hits <- vapply(seq_along(M), function(i) {
    pb <- p[starts[i]:ends[i]]
    pred <- which.max(pb)
    if (ioo) {
      po <- if (!is.null(p_out)) p_out[i] else 1 - sum(pb)
      if (po > pb[pred]) pred <- 0L
    }
    as.numeric(pred == d$choice_idx[i])
  }, numeric(1))

  sum(d$weights * hits) / sum(d$weights)
}

#' Print compact goodness-of-fit footer line
#'
#' Cats a single footer-style line, e.g.
#' \code{McFadden R2: 0.241 (adj: 0.232) | Hit rate: 0.516}.
#' Silently no-ops when the relevant fields are NA (e.g. when the model was
#' fitted with \code{keep_data = FALSE}).
#'
#' @param g A \code{choicer_gof} object.
#' @returns NULL invisibly.
#' @noRd
print_gof_lines <- function(g) {
  if (is.null(g)) return(invisible(NULL))
  parts <- character(0)
  if (!is.null(g$mcfadden_r2) && is.finite(g$mcfadden_r2)) {
    parts <- c(parts, sprintf("McFadden R2: %.3f (adj: %.3f)",
                              g$mcfadden_r2, g$mcfadden_r2_adj))
  }
  if (!is.null(g$hit_rate) && is.finite(g$hit_rate)) {
    parts <- c(parts, sprintf("Hit rate: %.3f", g$hit_rate))
  }
  if (length(parts) > 0) {
    cat(paste(parts, collapse = " | "), "\n")
  }
  invisible(NULL)
}

# --- Generic and method -------------------------------------------------------

#' Goodness of fit for a fitted choice model
#'
#' Computes McFadden's pseudo R-squared (plain and adjusted) and the in-sample
#' hit rate for a fitted model.
#'
#' Two null models are available for the pseudo R-squared
#' \eqn{R^2 = 1 - LL / LL_0} (adjusted:
#' \eqn{R^2_{adj} = 1 - (LL - K) / LL_0} with \eqn{K} the number of estimated
#' parameters):
#' \itemize{
#'   \item \code{"equal_shares"} (default): every alternative in individual
#'     \eqn{i}'s choice set is equally likely, so
#'     \eqn{LL_0 = -\sum_i w_i \log(M_i + 1_{outside})}. This is exact for
#'     unbalanced choice sets and arbitrary weights.
#'   \item \code{"market_shares"}: the maximized log-likelihood of an
#'     ASC-only model, \eqn{LL_0 = \sum_j N_j \log(s_j)} with \eqn{N_j} the
#'     choice counts and \eqn{s_j} the observed market shares (including the
#'     outside option when present). This closed form is valid only for
#'     balanced choice sets and uniform weights; otherwise an error suggests
#'     refitting an ASC-only model.
#' }
#'
#' The hit rate is the weighted share of individuals whose observed choice has
#' the highest predicted probability. When the model includes an outside
#' option, the outside good competes for the predicted maximum (its
#' probability is \eqn{1 - \sum_j p_{ij}}), and an individual predicted to
#' choose the outside good is a hit when they actually did.
#'
#' Both the null log-likelihood and the hit rate require the stored estimation
#' data; models fitted with \code{keep_data = FALSE} return NA fields with a
#' message.
#'
#' @param object A fitted model object (\code{choicer_mnl}, \code{choicer_mxl},
#'   or \code{choicer_nl}).
#' @param null Null model for the pseudo R-squared: \code{"equal_shares"}
#'   (default) or \code{"market_shares"}.
#' @param ... Additional arguments passed to methods.
#' @returns A \code{choicer_gof} object: a list with \code{loglik},
#'   \code{loglik_null}, \code{null}, \code{mcfadden_r2},
#'   \code{mcfadden_r2_adj}, \code{hit_rate}, \code{nobs}, and
#'   \code{n_params}.
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
#' gof(fit)
#' gof(fit, null = "market_shares")
#' }
#' @export
gof <- function(object, null = c("equal_shares", "market_shares"), ...) {
  UseMethod("gof")
}

#' @rdname gof
#' @export
gof.choicer_fit <- function(object, null = c("equal_shares", "market_shares"),
                            ...) {
  null <- match.arg(null)

  out <- list(
    loglik = object$loglik,
    loglik_null = NA_real_,
    null = null,
    mcfadden_r2 = NA_real_,
    mcfadden_r2_adj = NA_real_,
    hit_rate = NA_real_,
    nobs = object$nobs,
    n_params = object$n_params
  )

  if (is.null(object$data)) {
    message("gof() requires stored data for the null log-likelihood and hit ",
            "rate. Refit with keep_data = TRUE.")
    return(structure(out, class = "choicer_gof"))
  }

  d <- object$data
  ioo <- as.integer(object$include_outside_option)

  out$loglik_null <- if (null == "equal_shares") {
    -sum(d$weights * log(d$M + ioo))
  } else {
    .gof_market_shares_ll0(object)
  }

  out$mcfadden_r2 <- 1 - object$loglik / out$loglik_null
  out$mcfadden_r2_adj <- 1 - (object$loglik - object$n_params) / out$loglik_null
  out$hit_rate <- .gof_hit_rate(object)

  structure(out, class = "choicer_gof")
}

# --- print --------------------------------------------------------------------

#' Print goodness-of-fit measures
#'
#' @param x A \code{choicer_gof} object.
#' @param ... Additional arguments (ignored).
#' @returns The object invisibly.
#' @export
print.choicer_gof <- function(x, ...) {
  cat("Goodness of fit\n")
  cat("  Log-likelihood:      ", format(x$loglik, digits = 6), "\n")
  if (is.na(x$loglik_null)) {
    cat("  Null log-likelihood, McFadden R2 and hit rate unavailable ",
        "(refit with keep_data = TRUE).\n", sep = "")
  } else {
    null_label <- switch(x$null,
      equal_shares = "equal shares",
      market_shares = "market shares",
      x$null
    )
    cat("  Null log-likelihood: ", format(x$loglik_null, digits = 6),
        " (", null_label, ")\n", sep = "")
    cat(sprintf("  McFadden R2:          %.4f (adj: %.4f)\n",
                x$mcfadden_r2, x$mcfadden_r2_adj))
    cat(sprintf("  Hit rate:             %.4f\n", x$hit_rate))
  }
  cat("  N:", x$nobs, " | Parameters:", x$n_params, "\n")
  invisible(x)
}
