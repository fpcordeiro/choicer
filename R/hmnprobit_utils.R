# Hierarchical multinomial probit (HMNP) utilities: the data prep
# (prepare_hmnp_data) and the run_hmnprobit() wrapper around the hmnp_gibbs
# kernel. See vignettes/articles/hierarchical_mnp_math.Rmd for the math.

#' Prepare inputs for hierarchical multinomial probit estimation
#'
#' Prepares and validates panel (or cross-sectional) choice data for the
#' hierarchical Bayesian multinomial probit with iid \eqn{N(0, \sigma^2)}
#' utility shocks. The model shares its two-level random-effect structure
#' with [prepare_hmnl_data()]: respondent-level structural tastes
#' \eqn{\beta_i \sim N(b, W)} over the `covariate_cols` (normal only — the
#' probit keeps full conjugacy), and a global alternative-level effect
#' \eqn{\delta_j = z_j'\theta + \xi_j}, \eqn{\xi_j \sim N(0, \sigma_d^2)}.
#'
#' The returned structure is identical to [prepare_hmnl_data()] (both preps
#' share one internal engine), except there is no `rc_dist` field. Unlike
#' [prepare_mnp_data()], utilities are NOT differenced against a base
#' alternative: the iid-shock model works in un-differenced utility space,
#' so unbalanced choice sets are supported and the outside option is
#' implicit (its latent utility is a stochastic \eqn{N(0, \sigma^2)} draw in
#' the kernel, systematic utility 0).
#'
#' @inheritParams prepare_hmnl_data
#' @returns A list of class `c("choicer_data_hmnp", "list")` with the same
#'   components as [prepare_hmnl_data()] (minus `rc_dist`).
#' @seealso [prepare_hmnl_data()] for the component-by-component description.
#' @examples
#' library(data.table)
#' set.seed(42)
#' N <- 20; T <- 3; J <- 4
#' dt <- data.table(
#'   pid  = rep(1:N, each = T * J),
#'   task = rep(seq_len(N * T), each = J),
#'   alt  = rep(1:J, N * T)
#' )
#' dt[, `:=`(x1 = rnorm(.N), x2 = runif(.N, -1, 1))]
#' dt[, choice := 0L]
#' dt[, choice := if (runif(1) < 0.8) sample(c(1L, rep(0L, J - 1))) else 0L,
#'    by = task]
#' input <- prepare_hmnp_data(dt, "task", "alt", "choice", c("x1", "x2"),
#'                            person_col = "pid")
#' input$Ti[1:5]
#' input$alt_mapping
#' @export
prepare_hmnp_data <- function(
    data,
    id_col,
    alt_col,
    choice_col,
    covariate_cols,
    person_col = NULL,
    alt_covariate_cols = NULL,
    outside_opt_label = NULL,
    cf_residual_col = NULL,
    include_outside_option = TRUE
) {
  panel <- .prepare_hb_panel(
    data = data, id_col = id_col, alt_col = alt_col,
    choice_col = choice_col, covariate_cols = covariate_cols,
    person_col = person_col, alt_covariate_cols = alt_covariate_cols,
    outside_opt_label = outside_opt_label, cf_residual_col = cf_residual_col,
    include_outside_option = include_outside_option
  )
  structure(panel, class = c("choicer_data_hmnp", "list"))
}

#' Fit a hierarchical Bayesian multinomial probit (HMNP)
#'
#' Runs the fully conjugate Albert-Chib Gibbs sampler for the hierarchical
#' multinomial probit with iid normal utility shocks in un-differenced
#' utility space:
#' \deqn{U_{ijt} = x_{ijt}'\beta_i + \delta_j + \epsilon_{ijt}, \qquad
#'       U_{iot} = \epsilon_{iot}, \qquad \epsilon \sim N(0, \sigma^2),}
#' choice by argmax within the task including the stochastic implicit
#' outside option, \eqn{\beta_i \sim N(b, W)} (normal coordinates only —
#' log-normal would break conjugacy), and
#' \eqn{\delta_j = z_j'\theta + \xi_j}, \eqn{\xi_j \sim N(0, \sigma_d^2)}.
#'
#' \strong{Identification.} The probit likelihood is invariant to a common
#' rescaling of utilities and \eqn{\sigma}, so the chain runs on the
#' non-identified parameterization (free \eqn{\sigma^2}, better mixing via
#' parameter expansion) and every kept draw is normalized by the matching
#' power of the CURRENT \eqn{\sigma}: reported \eqn{b/\sigma},
#' \eqn{W/\sigma^2}, \eqn{\delta/\sigma}, \eqn{\theta/\sigma},
#' \eqn{\sigma_d^2/\sigma^2}. Raw chains are kept in `draws$*_raw`. The
#' outside option anchors the location of \eqn{\delta} exactly as in
#' [run_hmnlogit()].
#'
#' @inheritParams run_hmnlogit
#' @param input_data A `choicer_data_hmnp` object from [prepare_hmnp_data()].
#' @param prior As in [run_hmnlogit()], plus `a0` (3) and `s0` (3), the
#'   inverse-gamma shape/scale on the non-identified \eqn{\sigma^2}.
#' @param mcmc Named list overriding MCMC defaults: `R` (10000), `burn`
#'   (R %/% 5), `thin` (1), `seed` (drawn via `sample.int()` so `set.seed()`
#'   governs), `trace` (0). No proposal-scale settings: the sampler is fully
#'   conjugate, with no Metropolis steps.
#' @returns A `choicer_hmnp` object (classed `c("choicer_hmnp",
#'   "choicer_hb")`); the same layout as [run_hmnlogit()]'s return, with all
#'   reported summaries on the identified scale, raw chains in
#'   `draws$*_raw`, and the non-identified `draws$sigma2` trace.
#' @seealso [prepare_hmnp_data()], [simulate_hmnp_data()], [run_hmnlogit()]
#' @examples
#' \donttest{
#' sim <- simulate_hmnp_data(N = 100, T = 3, J = 4, seed = 42)
#' fit <- run_hmnprobit(sim$data, "task", "alt", "choice", c("x1", "x2"),
#'                      person_col = "pid", alt_covariate_cols = "z1",
#'                      mcmc = list(R = 500, burn = 200))
#' summary(fit)
#' }
#' @export
run_hmnprobit <- function(
    data = NULL,
    id_col = NULL,
    alt_col = NULL,
    choice_col = NULL,
    covariate_cols = NULL,
    person_col = NULL,
    alt_covariate_cols = NULL,
    outside_opt_label = NULL,
    cf_residual_col = NULL,
    input_data = NULL,
    include_outside_option = TRUE,
    prior = list(),
    mcmc = list(),
    chains = 1,
    keep_beta_i = c("means", "draws", "none"),
    keep_data = TRUE
) {
  cl <- match.call()
  keep_beta_i <- match.arg(keep_beta_i)

  # --- Resolve input pathway --------------------------------------------------
  has_data <- !is.null(data)
  has_input <- !is.null(input_data)
  if (has_data && has_input) {
    stop("Supply either 'data' (convenience) or 'input_data' (advanced), not both.")
  }
  if (!has_data && !has_input) {
    stop("Supply either 'data' (convenience) or 'input_data' (advanced).")
  }
  if (has_data) {
    if (is.null(id_col) || is.null(alt_col) || is.null(choice_col) ||
        is.null(covariate_cols)) {
      stop("Convenience workflow requires: id_col, alt_col, choice_col, ",
           "and covariate_cols.")
    }
    input_list <- prepare_hmnp_data(
      data, id_col, alt_col, choice_col, covariate_cols,
      person_col = person_col, alt_covariate_cols = alt_covariate_cols,
      outside_opt_label = outside_opt_label,
      cf_residual_col = cf_residual_col,
      include_outside_option = include_outside_option
    )
  } else {
    if (!inherits(input_data, "choicer_data_hmnp")) {
      stop("'input_data' must be a choicer_data_hmnp object from ",
           "prepare_hmnp_data().")
    }
    input_list <- input_data
  }

  if (!isTRUE(input_list$include_outside_option)) {
    stop("run_hmnprobit requires include_outside_option = TRUE: the implicit ",
         "outside option anchors the location of delta. A no-outside-good ",
         "mode is on the roadmap.")
  }

  X <- input_list$X
  Z <- input_list$Z
  K <- input_list$K_struct
  P <- input_list$P
  J <- input_list$J
  N <- input_list$N_persons
  n_tasks <- input_list$n_tasks

  # --- Priors -------------------------------------------------------------------
  b_bar <- prior$b_bar %||% rep(0, K)
  A <- prior$A %||% (0.01 * diag(K))
  nu <- prior$nu %||% (K + 3)
  V <- prior$V %||% (nu * diag(K))
  theta_bar <- prior$theta_bar %||% rep(0, P)
  A_theta <- prior$A_theta %||% (0.01 * diag(P))
  sd_prior <- utils::modifyList(
    list(half_cauchy = TRUE, s_d = 1, c0 = 3, d0 = 3),
    prior$sd_prior %||% list()
  )
  a0 <- prior$a0 %||% 3
  s0 <- prior$s0 %||% 3
  if (length(b_bar) != K) stop("prior$b_bar must have length ", K, ".")
  if (!is.matrix(A) || any(dim(A) != K)) {
    stop("prior$A must be a ", K, " x ", K, " matrix.")
  }
  if (nu < K) stop("prior$nu must be >= K = ", K, ".")
  if (!is.matrix(V) || any(dim(V) != K)) {
    stop("prior$V must be a ", K, " x ", K, " matrix.")
  }
  if (length(theta_bar) != P) stop("prior$theta_bar must have length ", P, ".")
  if (!is.matrix(A_theta) || any(dim(A_theta) != P)) {
    stop("prior$A_theta must be a ", P, " x ", P, " matrix.")
  }

  # --- MCMC settings --------------------------------------------------------------
  R <- as.integer(mcmc$R %||% 10000L)
  burn <- as.integer(mcmc$burn %||% (R %/% 5L))
  thin <- as.integer(mcmc$thin %||% 1L)
  trace <- as.integer(mcmc$trace %||% 0L)
  seed <- as.numeric(mcmc$seed %||% sample.int(.Machine$integer.max, 1L))
  if (R <= burn) stop("mcmc$R must be greater than mcmc$burn.")
  if (thin < 1L) stop("mcmc$thin must be >= 1.")
  chains <- as.integer(chains)
  if (chains < 1L) stop("chains must be >= 1.")

  keep_code <- switch(keep_beta_i, none = 0L, means = 1L, draws = 2L)
  if (keep_code == 2L) {
    R_keep_est <- (R - burn + thin - 1L) %/% thin
    bytes <- 8 * as.numeric(K) * N * R_keep_est
    if (bytes > 4e9) {
      stop("keep_beta_i = \"draws\" would allocate ",
           sprintf("%.1f GB", bytes / 1e9),
           " for the beta_i cube; reduce R or use keep_beta_i = \"means\".")
    }
    if (bytes > 1e9) {
      warning("keep_beta_i = \"draws\" allocates ",
              sprintf("%.1f GB", bytes / 1e9), " for the beta_i cube.",
              call. = FALSE)
    }
  }

  # No optimizer init: the sampler is fully conjugate and burn-in absorbs a
  # zero start. delta/theta start at 0 on the raw scale.
  delta_init <- rep(0, J)
  theta_init <- rep(0, P)

  # --- Run the Gibbs sampler ------------------------------------------------------
  run_chain <- function(chain_seed) {
    hmnp_gibbs(
      X = X, Z = Z, M = input_list$M, choice_pos = input_list$choice_pos,
      include_outside_option = TRUE, alt_of_row = input_list$alt_of_row,
      Ti = input_list$Ti, delta_init = delta_init, theta_init = theta_init,
      b_bar = b_bar, A = A, nu = nu, V = V,
      theta_bar = theta_bar, A_theta = A_theta, sd_prior = sd_prior,
      a0 = a0, s0 = s0, R = R, burn = burn, thin = thin, seed = chain_seed,
      keep_beta_i = keep_code, trace = trace
    )
  }
  elapsed <- system.time({
    out <- run_chain(seed)
    extra_chains <- if (chains > 1L) {
      lapply(seq_len(chains - 1L), function(cc) run_chain(seed + cc))
    } else {
      list()
    }
  })
  message("MCMC run time ", convertTime(elapsed))

  # --- Per-draw scale normalization (identified quantities) -------------------
  # E[b / sigma] != E[b] / E[sigma]: every identified summary is computed on
  # the per-draw normalized chains (the run_mnprobit discipline).
  normalize_chain <- function(ch) {
    s <- sqrt(as.numeric(ch$sigma2draw))
    list(
      b = ch$bdraw / s,
      w_vech = ch$wdraw / as.numeric(ch$sigma2draw),
      delta = ch$deltadraw / s,
      theta = ch$thetadraw / s,
      sigma_d2 = as.numeric(ch$sigma_d2draw) / as.numeric(ch$sigma2draw)
    )
  }
  idn <- normalize_chain(out)

  x_names <- colnames(X)
  z_names <- colnames(Z)
  am <- input_list$alt_mapping
  alt_labels <- as.character(am[am$alt_int > 0, ][[input_list$data_spec$alt_col]])
  w_names <- character(K * (K + 1L) / 2L)
  idx <- 1L
  for (a1 in seq_len(K)) {
    for (a2 in seq_len(a1)) {
      w_names[idx] <- sprintf("W[%s,%s]", x_names[a1], x_names[a2])
      idx <- idx + 1L
    }
  }
  colnames(idn$b) <- x_names
  colnames(idn$w_vech) <- w_names
  colnames(idn$delta) <- alt_labels
  colnames(idn$theta) <- z_names
  colnames(out$bdraw) <- x_names
  colnames(out$wdraw) <- w_names
  colnames(out$deltadraw) <- alt_labels
  colnames(out$thetadraw) <- z_names

  chain_blocks <- lapply(c(list(out), extra_chains), function(ch) {
    nn <- normalize_chain(ch)
    colnames(nn$b) <- x_names
    colnames(nn$theta) <- z_names
    list(b = nn$b, theta = nn$theta, sigma_d2 = nn$sigma_d2)
  })
  rhat_tab <- tryCatch(.hb_rhat_table(chain_blocks), error = function(e) NULL)

  coefficients <- colMeans(idn$b)
  se <- apply(idn$b, 2, stats::sd)
  vcov <- stats::cov(idn$b)

  W_mean <- matrix(0, K, K, dimnames = list(x_names, x_names))
  w_means <- colMeans(idn$w_vech)
  idx <- 1L
  for (a1 in seq_len(K)) {
    for (a2 in seq_len(a1)) {
      W_mean[a1, a2] <- W_mean[a2, a1] <- w_means[idx]
      idx <- idx + 1L
    }
  }

  delta_tab <- data.table::data.table(
    alternative = alt_labels,
    mean = as.numeric(out$delta_mean),
    sd = as.numeric(out$delta_sd)
  )
  xi_tab <- data.table::data.table(
    alternative = alt_labels,
    mean = as.numeric(out$xi_mean),
    sd = as.numeric(out$xi_sd)
  )

  beta_i <- NULL
  if (keep_code >= 1L) {
    beta_i <- list(
      mean = matrix(out$beta_i_mean, K, N,
                    dimnames = list(x_names, as.character(input_list$person_ids))),
      sd = matrix(out$beta_i_sd, K, N,
                  dimnames = list(x_names, as.character(input_list$person_ids)))
    )
    if (keep_code == 2L) beta_i$draws <- out$beta_i_draws
  }

  if (N < K) {
    warning("Fewer respondents than structural covariates (N = ", N,
            " < K = ", K, "): W is largely prior-driven.", call. = FALSE)
  }
  if (J < 3 * P) {
    warning("Few alternatives relative to the delta mean function (J = ", J,
            ", P = ", P, "): theta and sigma_d^2 lean on the prior.",
            call. = FALSE)
  }
  if (!is.null(rhat_tab) && any(rhat_tab > 1.05, na.rm = TRUE)) {
    warning("split-R-hat exceeds 1.05 for: ",
            paste(names(rhat_tab)[which(rhat_tab > 1.05)], collapse = ", "),
            ". Consider a longer run.", call. = FALSE)
  }

  new_choicer_hmnp(
    call = cl,
    coefficients = coefficients,
    se = se,
    vcov = vcov,
    theta_summary = build_bayes_coef_table(idn$theta),
    sigma_d2_summary = build_bayes_coef_table(
      matrix(idn$sigma_d2, ncol = 1, dimnames = list(NULL, "sigma_d^2"))
    ),
    W_mean = W_mean,
    delta = delta_tab,
    xi = xi_tab,
    beta_i = beta_i,
    draws = list(
      b = idn$b, w_vech = idn$w_vech, delta = idn$delta, theta = idn$theta,
      sigma_d2 = idn$sigma_d2, sigma2 = as.numeric(out$sigma2draw),
      b_raw = out$bdraw, w_vech_raw = out$wdraw, delta_raw = out$deltadraw,
      theta_raw = out$thetadraw, sigma_d2_raw = as.numeric(out$sigma_d2draw)
    ),
    rhat = rhat_tab,
    prior = list(b_bar = b_bar, A = A, nu = nu, V = V, theta_bar = theta_bar,
                 A_theta = A_theta, sd_prior = sd_prior, a0 = a0, s0 = s0),
    mcmc = list(R = R, burn = burn, thin = thin, seed = seed, trace = trace,
                chains = chains, R_keep = out$R_keep),
    nobs = n_tasks,
    n_persons = N,
    J = J,
    P = P,
    K_struct = K,
    param_map = input_list$param_map,
    alt_mapping = am,
    data_spec = input_list$data_spec,
    cf_active = !is.null(input_list$data_spec$cf_residual_col),
    sampler = list(name = "hmnp_gibbs",
                   elapsed_time = convertTime(elapsed)),
    data = if (keep_data) input_list else NULL
  )
}
