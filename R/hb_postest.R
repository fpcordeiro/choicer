# Post-estimation methods for the hierarchical Bayes fits (choicer_hmnl /
# choicer_hmnp), registered on the package's existing generics. Everything
# integrates over posterior draws (and, at the population level, over the
# random-coefficient distribution), returning posterior summaries rather
# than delta-method points. See _plans/hierarchical_bayes_plan.md,
# "Post-estimation".
#
# Probability engines:
#   * HMNL: closed-form softmax against the implicit outside option,
#     P(j) = exp(V_j) / (1 + sum_k exp(V_k)).
#   * HMNP: the iid-normal shocks give the 1-D integral
#     P(j) = int phi(u) prod_{k != j} Phi(V_j - V_k + u) du (outside V = 0,
#     identified scale sigma = 1), evaluated by fixed-node Gauss-Hermite
#     quadrature (.gauss_hermite) — deterministic, no simulation noise.
#     logsum / consumer_surplus stay HMNL-only in v1: the logsum formula is
#     EV1-specific and probit E[max U] has no closed form (simulated-Emax is
#     roadmapped).

# --- internal helpers ---------------------------------------------------------

#' Gauss-Hermite nodes and weights via Golub-Welsch
#'
#' Physicists' convention: int exp(-t^2) f(t) dt = sum w_m f(t_m). For
#' N(0, 1) expectations use int phi(u) f(u) du = sum (w_m / sqrt(pi))
#' f(sqrt(2) t_m).
#' @noRd
.gauss_hermite <- function(n) {
  i <- seq_len(n - 1L)
  Jm <- matrix(0, n, n)
  off <- sqrt(i / 2)
  Jm[cbind(i, i + 1L)] <- off
  Jm[cbind(i + 1L, i)] <- off
  eg <- eigen(Jm, symmetric = TRUE)
  idx <- order(eg$values)
  list(
    nodes = eg$values[idx],
    weights = (eg$vectors[1L, idx])^2 * sqrt(pi)
  )
}

#' Thinned indices into the kept posterior draws
#' @noRd
.hb_draw_index <- function(object, n_draws) {
  R_keep <- nrow(object$draws$b)
  if (n_draws >= R_keep) seq_len(R_keep)
  else unique(round(seq(1L, R_keep, length.out = n_draws)))
}

#' Random-coefficient transform: chain scale -> utility scale
#' @noRd
.hb_gamma <- function(beta, rc_dist) {
  ln <- which(rc_dist == 1L)
  if (length(ln)) beta[ln] <- exp(beta[ln])
  beta
}

#' rc_dist lookup with cross-version fallbacks (HMNP fits are all-normal)
#' @noRd
.hb_rc_dist <- function(object) {
  rc <- object$rc_dist %||% object$data_spec$rc_dist %||%
    object$data$rc_dist
  if (is.null(rc)) rc <- rep(0L, object$K_struct)
  as.integer(rc)
}

#' Resolve prediction data for a hierarchical Bayes fit
#'
#' Builds the prediction panel from `newdata` (or the stored estimation
#' data): the structural design matrix in estimation column order, the task
#' grouping, and the per-row alternative resolution. Alternatives never seen
#' in estimation are matched by their `alt_covariate_cols` values and get a
#' posterior-predictive delta, N(z_new' theta_r, sigma_d2_r) per draw — the
#' entry-counterfactual mechanism of the random-effects delta.
#' @noRd
.hb_resolve_newdata <- function(object, newdata) {
  spec <- object$data_spec
  if (is.null(newdata)) {
    if (is.null(object$data)) {
      stop("No `newdata` supplied and the fit was built with keep_data = ",
           "FALSE; refit with keep_data = TRUE or pass `newdata`.")
    }
    # Rebuild from the stored prep: rows are already sorted and encoded.
    d <- object$data
    return(list(
      X = d$X,
      task_of_row = rep(seq_along(d$M), times = d$M),
      n_tasks = d$n_tasks,
      alt_label = as.character(
        d$alt_mapping[d$alt_mapping$alt_int > 0, ][[spec$alt_col]]
      )[d$alt_of_row],
      known_idx = d$alt_of_row,
      z_new = NULL,
      new_labels = character(0),
      person = if (!is.null(spec$person_col)) {
        # one entry per ROW: person per task, expanded by the task sizes
        rep(rep(d$person_ids, times = d$Ti), times = d$M)
      }
    ))
  }

  dt <- data.table::as.data.table(newdata)
  needed <- unique(c(spec$person_col, spec$id_col, spec$alt_col,
                     spec$covariate_cols, spec$alt_covariate_cols,
                     spec$cf_residual_col))
  missing_cols <- setdiff(needed, names(dt))
  if (length(missing_cols)) {
    stop("`newdata` is missing columns used in estimation: ",
         paste(missing_cols, collapse = ", "))
  }
  x_cols <- colnames(object$draws$b)      # estimation column order
  if (!all(x_cols %in% names(dt))) {
    stop("`newdata` is missing structural covariate columns: ",
         paste(setdiff(x_cols, names(dt)), collapse = ", "))
  }

  ord_cols <- c(spec$person_col, spec$id_col, spec$alt_col)
  data.table::setorderv(dt, ord_cols)
  task_key <- if (is.null(spec$person_col)) {
    as.character(dt[[spec$id_col]])
  } else {
    paste(dt[[spec$person_col]], dt[[spec$id_col]], sep = "\r")
  }
  task_of_row <- as.integer(factor(task_key, levels = unique(task_key)))

  X <- as.matrix(dt[, x_cols, with = FALSE])
  if (!all(is.finite(X))) stop("`newdata` covariates must be finite.")

  am <- object$alt_mapping
  known_labels <- as.character(am[am$alt_int > 0, ][[spec$alt_col]])
  alt_label <- as.character(dt[[spec$alt_col]])
  known_idx <- match(alt_label, known_labels)   # NA = new alternative

  new_labels <- sort(unique(alt_label[is.na(known_idx)]))
  z_new <- NULL
  if (length(new_labels)) {
    z_cols <- setdiff(colnames(object$draws$theta), "(Intercept)")
    if (length(z_cols) && !all(z_cols %in% names(dt))) {
      stop("New alternatives in `newdata` (",
           paste(new_labels, collapse = ", "),
           ") need their alternative-level covariates: ",
           paste(setdiff(z_cols, names(dt)), collapse = ", "))
    }
    z_new <- matrix(1, nrow = length(new_labels), ncol = 1L,
                    dimnames = list(new_labels, "(Intercept)"))
    for (zc in z_cols) {
      vals <- vapply(new_labels, function(lb) {
        v <- unique(dt[[zc]][alt_label == lb])
        if (length(v) != 1L) {
          stop("Alternative-level covariate `", zc,
               "` is not constant within new alternative `", lb, "`.")
        }
        as.numeric(v)
      }, numeric(1L))
      z_new <- cbind(z_new, vals)
      colnames(z_new)[ncol(z_new)] <- zc
    }
  }

  list(
    X = X,
    task_of_row = task_of_row,
    n_tasks = max(task_of_row),
    alt_label = alt_label,
    known_idx = known_idx,
    z_new = z_new,
    new_labels = new_labels,
    person = if (!is.null(spec$person_col)) dt[[spec$person_col]]
  )
}

#' Per-row delta for one posterior draw (posterior-predictive for new alts)
#' @noRd
.hb_delta_row <- function(object, rd, r, delta_new_draws) {
  out <- numeric(length(rd$alt_label))
  known <- !is.na(rd$known_idx)
  out[known] <- object$draws$delta[r, rd$known_idx[known]]
  if (any(!known)) {
    out[!known] <- delta_new_draws[[r]][rd$alt_label[!known]]
  }
  out
}

#' Posterior-predictive delta draws for alternatives outside the estimation
#' sample: delta_new ~ N(z_new' theta_r, sigma_d2_r), one draw per kept
#' posterior draw (set.seed() governs reproducibility).
#' @noRd
.hb_delta_new_draws <- function(object, rd, idx) {
  if (is.null(rd$z_new)) return(NULL)
  theta_cols <- colnames(object$draws$theta)
  Zn <- rd$z_new[, theta_cols, drop = FALSE]
  out <- vector("list", nrow(object$draws$b))
  for (r in idx) {
    mu <- drop(Zn %*% object$draws$theta[r, ])
    out[[r]] <- stats::setNames(
      stats::rnorm(length(mu), mu, sqrt(object$draws$sigma_d2[r])),
      rownames(rd$z_new)
    )
  }
  out
}

#' Task-level choice probabilities for one draw (both engines)
#'
#' V is the inside-row utility index; returns a list with `p_inside` (per
#' row) and `p_outside` (per task).
#' @noRd
.hb_task_probs <- function(V, task_of_row, n_tasks, model, gh = NULL) {
  if (model == "hmnl") {
    eV <- exp(V)
    denom <- 1 + as.numeric(rowsum(eV, task_of_row, reorder = TRUE))
    p_inside <- eV / denom[task_of_row]
    return(list(p_inside = p_inside, p_outside = 1 / denom))
  }
  # HMNP: Gauss-Hermite over the common shock argument, per task.
  p_inside <- numeric(length(V))
  p_outside <- numeric(n_tasks)
  u <- sqrt(2) * gh$nodes
  w <- gh$weights / sqrt(pi)
  rows_by_task <- split(seq_along(V), task_of_row)
  for (t in seq_len(n_tasks)) {
    rows <- rows_by_task[[t]]
    v <- c(V[rows], 0)                    # outside last
    m <- length(v)
    p <- numeric(m)
    for (j in seq_len(m)) {
      f <- rep(1, length(u))
      for (k in seq_len(m)) {
        if (k != j) f <- f * stats::pnorm(v[j] - v[k] + u)
      }
      p[j] <- sum(w * f)
    }
    p <- p / sum(p)                       # guard quadrature roundoff
    p_inside[rows] <- p[seq_len(m - 1L)]
    p_outside[t] <- p[m]
  }
  list(p_inside = p_inside, p_outside = p_outside)
}

#' Shared prediction engine: aggregate shares (or per-row means) over
#' posterior draws
#' @noRd
.hb_predict_core <- function(object, rd, idx, level, n_gh = 20L) {
  K <- object$K_struct
  rc <- .hb_rc_dist(object)
  gh <- if (object$model == "hmnp") .gauss_hermite(n_gh)
  delta_new <- .hb_delta_new_draws(object, rd, idx)

  use_beta_draws <- identical(level, "individual") &&
    !is.null(object$beta_i$draws)
  if (identical(level, "individual")) {
    if (is.null(rd$person)) {
      stop("level = \"individual\" requires a person column in the ",
           "prediction data (person_col was NULL in estimation).")
    }
    person_pos <- match(as.character(rd$person),
                        colnames(object$beta_i$mean))
    if (anyNA(person_pos)) {
      stop("level = \"individual\": `newdata` contains respondents not in ",
           "the estimation sample.")
    }
    if (!use_beta_draws) {
      message("Using posterior-mean beta_i (keep_beta_i = \"draws\" was not ",
              "set); individual probabilities are plug-in, not fully ",
              "posterior-integrated.")
    }
  }

  labels <- c(sort(unique(rd$alt_label)), "(outside)")
  share_draws <- matrix(0, length(idx), length(labels),
                        dimnames = list(NULL, labels))
  row_prob_sum <- numeric(length(rd$alt_label))

  for (s in seq_along(idx)) {
    r <- idx[s]
    if (identical(level, "population")) {
      # one beta draw from N(b_r, W_r) per posterior draw
      W <- matrix(0, K, K)
      W[lower.tri(W, diag = TRUE)] <- object$draws$w_vech[r, ]
      W <- W + t(W) - diag(diag(W))
      L <- t(chol(W + diag(1e-10, K)))
      beta <- object$draws$b[r, ] + drop(L %*% stats::rnorm(K))
      gamma_row <- matrix(.hb_gamma(beta, rc), nrow = 1L)
      Vx <- drop(rd$X %*% t(gamma_row))
    } else {
      if (use_beta_draws) {
        bmat <- object$beta_i$draws[, person_pos, r, drop = FALSE]
        bmat <- matrix(bmat, nrow = K)
      } else {
        bmat <- object$beta_i$mean[, person_pos, drop = FALSE]
      }
      for (k in which(rc == 1L)) bmat[k, ] <- exp(bmat[k, ])
      Vx <- rowSums(rd$X * t(bmat))
    }
    V <- Vx + .hb_delta_row(object, rd, r, delta_new)
    pr <- .hb_task_probs(V, rd$task_of_row, rd$n_tasks, object$model, gh)

    agg <- rowsum(pr$p_inside, rd$alt_label, reorder = TRUE)
    share_draws[s, rownames(agg)] <- agg[, 1L] / rd$n_tasks
    share_draws[s, "(outside)"] <- mean(pr$p_outside)
    row_prob_sum <- row_prob_sum + pr$p_inside
  }

  list(share_draws = share_draws,
       row_prob_mean = row_prob_sum / length(idx))
}

# --- predict ------------------------------------------------------------------

#' Posterior choice probabilities and shares for hierarchical Bayes fits
#'
#' Computes counterfactual choice probabilities, integrating over the
#' posterior draws and (at the population level) over the random-coefficient
#' distribution: for each kept draw \eqn{(b_r, W_r, \delta_r, \ldots)} one
#' \eqn{\beta \sim N(b_r, W_r)} is drawn and the model probabilities are
#' averaged. Alternatives in `newdata` that were not in the estimation
#' sample receive a posterior-predictive
#' \eqn{\delta_{new} \sim N(z_{new}'\theta_r, \sigma_{d,r}^2)} — the entry
#' counterfactual unlocked by the random-effects \eqn{\delta}. Price or
#' subsidy counterfactuals are just modified covariate columns in `newdata`.
#'
#' HMNL probabilities are closed-form logit; HMNP probabilities use the
#' 1-D Gauss-Hermite representation of the iid-probit integral
#' \eqn{P(j) = \int \phi(u) \prod_{k \ne j} \Phi(V_j - V_k + u) du}.
#'
#' @param object A `choicer_hmnl` or `choicer_hmnp` fit.
#' @param newdata Data frame with the estimation columns (choice column not
#'   required). `NULL` (default) predicts on the estimation data.
#' @param level `"population"` (default) integrates over \eqn{N(b, W)};
#'   `"individual"` uses the respondent-level \eqn{\beta_i} (requires the
#'   prediction rows to belong to estimation respondents; fully
#'   posterior-integrated when the fit kept `keep_beta_i = "draws"`).
#' @param n_draws Number of posterior draws to integrate over (thinned
#'   evenly from the kept draws; default 200).
#' @param aggregate If `TRUE` (default) return a per-alternative posterior
#'   share table (including the outside option); if `FALSE` return the
#'   posterior-mean probability per row of the prediction data.
#' @param ... Ignored.
#' @returns With `aggregate = TRUE`, a `data.table` with columns
#'   `alternative`, `share` (posterior mean), `sd`, `lower`, `upper` (95%
#'   equal-tailed interval); the posterior share draws are attached as
#'   `attr(, "draws")`. With `aggregate = FALSE`, a numeric vector of
#'   posterior-mean choice probabilities, one per prediction row.
#' @examples
#' \donttest{
#' sim <- simulate_hmnl_data(N = 100, T = 3, J = 4, seed = 42)
#' fit <- run_hmnlogit(sim$data, "task", "alt", "choice", c("x1", "x2"),
#'                     person_col = "pid", alt_covariate_cols = "z1",
#'                     mcmc = list(R = 500, burn = 200))
#' predict(fit)                       # posterior shares, estimation data
#' cf <- sim$data
#' cf$x1 <- cf$x1 + 0.5               # a counterfactual attribute change
#' predict(fit, newdata = cf)
#' }
#' @export
predict.choicer_hb <- function(object, newdata = NULL,
                               level = c("population", "individual"),
                               n_draws = 200L, aggregate = TRUE, ...) {
  level <- match.arg(level)
  rd <- .hb_resolve_newdata(object, newdata)
  idx <- .hb_draw_index(object, n_draws)
  core <- .hb_predict_core(object, rd, idx, level)

  if (!aggregate) return(core$row_prob_mean)

  sd_ <- apply(core$share_draws, 2, stats::sd)
  qs <- t(apply(core$share_draws, 2, stats::quantile,
                probs = c(0.025, 0.975)))
  out <- data.table::data.table(
    alternative = colnames(core$share_draws),
    share = colMeans(core$share_draws),
    sd = sd_,
    lower = qs[, 1L],
    upper = qs[, 2L]
  )
  data.table::setattr(out, "draws", core$share_draws)
  out[]
}

# --- wtp ----------------------------------------------------------------------

#' @describeIn wtp Posterior willingness-to-pay for hierarchical Bayes fits:
#'   the per-draw ratio of population-mean utility coefficients,
#'   \eqn{-\bar\gamma_{attr} / \bar\gamma_{price}} (for log-normal
#'   coordinates \eqn{\bar\gamma = \exp(b + W_{kk}/2)}). Ratio posteriors
#'   are heavy-tailed, so the point estimate is the posterior **median**
#'   with equal-tailed quantile intervals — never a posterior mean or a
#'   delta-method SE. A warning is raised when the price coefficient's sign
#'   is not resolved by the posterior. If the price variable was flagged as
#'   endogenous-without-a-control-function at prep time, WTP inherits that
#'   caveat (see `cf_residual_col` in [prepare_hmnl_data()]).
#' @export
wtp.choicer_hb <- function(object, price_var, attr_vars = NULL,
                           level = 0.95, ...) {
  x_names <- colnames(object$draws$b)
  if (!price_var %in% x_names) {
    stop("`price_var` must be one of: ", paste(x_names, collapse = ", "))
  }
  if (is.null(attr_vars)) attr_vars <- setdiff(x_names, price_var)
  if (!all(attr_vars %in% x_names)) {
    stop("Unknown `attr_vars`: ",
         paste(setdiff(attr_vars, x_names), collapse = ", "))
  }
  rc <- .hb_rc_dist(object)
  K <- object$K_struct
  diag_idx <- vapply(seq_len(K), function(k) (k * (k + 1L)) %/% 2L,
                     integer(1L))
  # per-draw population-mean utility coefficients
  gbar <- object$draws$b
  for (k in which(rc == 1L)) {
    gbar[, k] <- exp(object$draws$b[, k] +
                       0.5 * object$draws$w_vech[, diag_idx[k]])
  }
  price_draws <- gbar[, price_var]
  p_pos <- mean(price_draws >= 0)
  if (min(p_pos, 1 - p_pos) > 0.01) {
    warning(sprintf(
      "The sign of `%s` is not resolved by the posterior (P(>= 0) = %.2f): ",
      price_var, p_pos),
      "WTP ratios are ill-defined; consider a log-normal price coefficient ",
      "(rc_dist) or more data.", call. = FALSE)
  }

  alpha <- (1 - level) / 2
  rows <- lapply(attr_vars, function(av) {
    ratio <- -gbar[, av] / price_draws
    qs <- stats::quantile(ratio, probs = c(alpha, 0.5, 1 - alpha))
    data.table::data.table(
      attribute = av, wtp = qs[[2L]], lower = qs[[1L]], upper = qs[[3L]]
    )
  })
  out <- data.table::rbindlist(rows)
  data.table::setattr(out, "price_var", price_var)
  out[]
}

# --- logsum / consumer surplus (HMNL only in v1) -------------------------------

#' @param n_draws Number of posterior draws to integrate over (hierarchical
#'   Bayes methods; thinned evenly from the kept draws).
#' @describeIn logsum Posterior expected logsum for the hierarchical logit:
#'   per choice situation, \eqn{\log(1 + \sum_j \exp V_j)} against the
#'   outside-option anchor, averaged over posterior draws with one
#'   \eqn{\beta \sim N(b_r, W_r)} draw each. Returns the per-task posterior
#'   mean vector.
#' @export
logsum.choicer_hmnl <- function(object, newdata = NULL, n_draws = 200L, ...) {
  rd <- .hb_resolve_newdata(object, newdata)
  idx <- .hb_draw_index(object, n_draws)
  rc <- .hb_rc_dist(object)
  K <- object$K_struct
  delta_new <- .hb_delta_new_draws(object, rd, idx)

  ls_draws <- matrix(0, length(idx), rd$n_tasks)
  for (s in seq_along(idx)) {
    r <- idx[s]
    W <- matrix(0, K, K)
    W[lower.tri(W, diag = TRUE)] <- object$draws$w_vech[r, ]
    W <- W + t(W) - diag(diag(W))
    L <- t(chol(W + diag(1e-10, K)))
    beta <- object$draws$b[r, ] + drop(L %*% stats::rnorm(K))
    V <- drop(rd$X %*% .hb_gamma(beta, rc)) +
      .hb_delta_row(object, rd, r, delta_new)
    denom <- 1 + as.numeric(rowsum(exp(V), rd$task_of_row, reorder = TRUE))
    ls_draws[s, ] <- log(denom)
  }
  out <- colMeans(ls_draws)
  attr(out, "draws") <- ls_draws
  out
}

#' @describeIn logsum The probit expected maximum has no closed form;
#'   simulated-Emax surplus for the HMNP is on the roadmap.
#' @export
logsum.choicer_hmnp <- function(object, newdata = NULL, ...) {
  stop("logsum()/consumer_surplus() are logit-only in this version: the ",
       "EV1 logsum formula does not apply to probit shocks, and the ",
       "simulated-Emax variant is on the roadmap. Use the HMNL for welfare ",
       "analysis.")
}

#' @param n_draws Number of posterior draws to integrate over (hierarchical
#'   Bayes methods).
#' @describeIn consumer_surplus Posterior consumer surplus for the
#'   hierarchical logit: per-task logsum divided by the (positive) marginal
#'   utility of income \eqn{-\bar\gamma_{price}}, per posterior draw. With
#'   `newdata`, the return also carries the compensating variation against
#'   the estimation data (`attr(, "cv")`), i.e. the posterior of
#'   \eqn{(\mathrm{logsum}_{new} - \mathrm{logsum}_{base}) /
#'   (-\bar\gamma_{price})} summed over tasks. Requires a fixed-sign price
#'   coefficient; the posterior-median ratio discipline of
#'   [wtp.choicer_hb()] applies.
#' @export
consumer_surplus.choicer_hmnl <- function(object, price_var, newdata = NULL,
                                          level = 0.95, weights = NULL,
                                          n_draws = 200L, ...) {
  x_names <- colnames(object$draws$b)
  if (!price_var %in% x_names) {
    stop("`price_var` must be one of: ", paste(x_names, collapse = ", "))
  }
  rc <- .hb_rc_dist(object)
  K <- object$K_struct
  diag_idx <- vapply(seq_len(K), function(k) (k * (k + 1L)) %/% 2L,
                     integer(1L))
  gprice <- object$draws$b[, price_var]
  if (rc[match(price_var, x_names)] == 1L) {
    gprice <- exp(gprice + 0.5 * object$draws$w_vech[, diag_idx[
      match(price_var, x_names)]])
  }
  p_pos <- mean(gprice >= 0)
  if (min(p_pos, 1 - p_pos) > 0.01) {
    warning("The sign of the price coefficient is not resolved by the ",
            "posterior; surplus ratios are ill-defined.", call. = FALSE)
  }

  idx <- .hb_draw_index(object, n_draws)
  if (!exists(".Random.seed", envir = globalenv())) stats::runif(1L)
  set.seed_state <- .Random.seed          # reuse one beta path for both
  ls_base <- logsum.choicer_hmnl(object, newdata = NULL, n_draws = n_draws)
  cs_target <- if (is.null(newdata)) {
    ls_base
  } else {
    assign(".Random.seed", set.seed_state, envir = globalenv())
    logsum.choicer_hmnl(object, newdata = newdata, n_draws = n_draws)
  }

  alpha_q <- (1 - level) / 2
  a_draws <- -gprice[idx]                 # marginal utility of income
  cs_draws <- attr(cs_target, "draws") / a_draws
  cs_task <- apply(cs_draws, 2, stats::median)
  out <- data.table::data.table(
    task = seq_along(cs_task),
    cs = cs_task,
    lower = apply(cs_draws, 2, stats::quantile, probs = alpha_q),
    upper = apply(cs_draws, 2, stats::quantile, probs = 1 - alpha_q)
  )
  if (!is.null(newdata)) {
    cv_draws <- rowSums(attr(cs_target, "draws") - attr(ls_base, "draws")) /
      a_draws
    data.table::setattr(out, "cv", stats::quantile(
      cv_draws, probs = c(alpha_q, 0.5, 1 - alpha_q)))
  }
  out[]
}

#' @describeIn consumer_surplus Not available for the probit (see
#'   [logsum.choicer_hmnp()]); roadmapped via simulated Emax.
#' @export
consumer_surplus.choicer_hmnp <- function(object, price_var, newdata = NULL,
                                          level = 0.95, weights = NULL, ...) {
  logsum.choicer_hmnp(object)
}

# --- elasticities / diversion --------------------------------------------------

#' Shared perturbation engine: baseline and perturbed shares per alternative
#' @noRd
.hb_perturb_shares <- function(object, elast_var, eps, n_draws) {
  if (is.null(object$data)) {
    stop("Elasticities need the estimation data; refit with keep_data = TRUE.")
  }
  rd <- .hb_resolve_newdata(object, NULL)
  idx <- .hb_draw_index(object, n_draws)
  x_names <- colnames(object$draws$b)
  if (!elast_var %in% x_names) {
    stop("`elast_var` must be one of: ", paste(x_names, collapse = ", "))
  }
  if (!exists(".Random.seed", envir = globalenv())) stats::runif(1L)
  seed_state <- .Random.seed

  base <- .hb_predict_core(object, rd, idx, "population")
  labels <- colnames(base$share_draws)
  inside <- setdiff(labels, "(outside)")
  J <- length(inside)

  # Perturb elast_var by (1 + eps) for one inside alternative at a time,
  # replaying the same beta path (same RNG state) so the contrast is exact.
  pert <- array(0, dim = c(length(idx), length(labels), J),
                dimnames = list(NULL, labels, inside))
  xk <- match(elast_var, x_names)
  for (jj in seq_len(J)) {
    rd_p <- rd
    rows_j <- rd$alt_label == inside[jj]
    rd_p$X[rows_j, xk] <- rd$X[rows_j, xk] * (1 + eps)
    assign(".Random.seed", seed_state, envir = globalenv())
    pert[, , jj] <- .hb_predict_core(object, rd_p, idx, "population")$share_draws
  }
  list(base = base$share_draws, pert = pert, inside = inside,
       labels = labels, idx = idx)
}

#' @param elast_var Structural covariate to perturb (hierarchical Bayes
#'   methods).
#' @param eps Relative perturbation size (default 0.01).
#' @param n_draws Number of posterior draws to integrate over.
#' @describeIn elasticities Posterior-mean aggregate arc elasticities for
#'   hierarchical Bayes fits: each inside alternative's `elast_var` is
#'   perturbed by `eps` (default 1%) and shares are recomputed per posterior
#'   draw with a common random-coefficient path, giving
#'   \eqn{E_{jk} = (\Delta s_j / s_j) / \epsilon}. Rows are responding
#'   alternatives (including the outside option), columns the perturbed
#'   alternative.
#' @export
elasticities.choicer_hb <- function(object, elast_var, eps = 0.01,
                                    n_draws = 100L, ...) {
  ps <- .hb_perturb_shares(object, elast_var, eps, n_draws)
  base_mean <- colMeans(ps$base)
  out <- matrix(
    0, length(ps$labels), length(ps$inside),
    dimnames = list(ps$labels, ps$inside)
  )
  for (jj in seq_along(ps$inside)) {
    d_share <- colMeans(ps$pert[, , jj]) - base_mean
    out[, jj] <- (d_share / base_mean) / eps
  }
  out
}

#' @param elast_var Structural covariate to perturb (hierarchical Bayes
#'   methods).
#' @param eps Relative perturbation size (default 0.01).
#' @param n_draws Number of posterior draws to integrate over.
#' @describeIn diversion_ratios Posterior-mean diversion ratios for
#'   hierarchical Bayes fits, from the same perturbation engine as
#'   [elasticities.choicer_hb()]: \eqn{DR(j \to k)} is the fraction of the
#'   share alternative j loses (when its `elast_var` worsens) that flows to
#'   k — including the outside option. Columns are the perturbed
#'   alternative j, rows the receiving alternative k; the diagonal is 0.
#' @export
diversion_ratios.choicer_hb <- function(object, elast_var, eps = 0.01,
                                        n_draws = 100L, ...) {
  ps <- .hb_perturb_shares(object, elast_var, eps, n_draws)
  base_mean <- colMeans(ps$base)
  out <- matrix(
    0, length(ps$labels), length(ps$inside),
    dimnames = list(ps$labels, ps$inside)
  )
  for (jj in seq_along(ps$inside)) {
    d_share <- colMeans(ps$pert[, , jj]) - base_mean
    own <- d_share[ps$inside[jj]]
    others <- setdiff(ps$labels, ps$inside[jj])
    out[others, jj] <- d_share[others] / (-own)
    out[ps$inside[jj], jj] <- 0
  }
  out
}
