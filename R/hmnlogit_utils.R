# Hierarchical multinomial logit (HMNL) utilities. Phase 0 ships the data
# prep; run_hmnlogit() and the Gibbs kernel arrive with the HMNL phase of
# _plans/hierarchical_bayes_plan.md.

#' Prepare inputs for hierarchical multinomial logit estimation
#'
#' Prepares and validates panel (or cross-sectional) choice data for the
#' hierarchical Bayesian multinomial logit. The model has two random-effect
#' levels: respondent-level structural tastes \eqn{\beta_i \sim N(b, W)}
#' over the `covariate_cols`, and a global alternative-level effect
#' \eqn{\delta_j = z_j'\theta + \xi_j}, \eqn{\xi_j \sim N(0, \sigma_d^2)},
#' with mean-function design \eqn{z_j} built from `alt_covariate_cols`.
#'
#' \strong{Structure.} The design matrix `X` carries structural covariates
#' only — no alternative-specific-constant dummies. The alternative effect
#' \eqn{\delta_j} is indexed by `alt_of_row` (integer codes `1..J`), so
#' memory and compute scale with the number of rows, not with `J` extra
#' design columns.
#'
#' \strong{Outside option.} With `include_outside_option = TRUE` (the
#' default) the outside good is modelled implicitly, following the
#' [prepare_mnl_data()] convention: physical outside rows (identified by
#' `outside_opt_label`) are removed, the estimation kernels add the outside
#' term (systematic utility 0), and a choice situation whose inside rows are
#' all `0` in `choice_col` is coded as "outside chosen" (`choice_pos = 0`).
#' The outside option anchors the location of \eqn{\delta} (mean utility
#' relative to the outside good).
#'
#' \strong{Cross-section vs panel.} `person_col` groups choice situations
#' into respondents sharing one \eqn{\beta_i}. With `person_col = NULL`
#' (default) every choice situation is its own respondent (`Ti` all 1) —
#' the cross-sectional random-coefficients mode.
#'
#' \strong{Control function.} `cf_residual_col` (a user-supplied first-stage
#' residual, Petrin & Train 2010) is appended to `X` as an ordinary
#' covariate; its provenance is recorded in `data_spec`. The first stage is
#' NOT run here — supplying a valid residual is the user's responsibility.
#'
#' @param data Data frame containing choice data.
#' @param id_col Name of the column identifying choice situations (tasks).
#'   Task ids only need to be unique within a respondent.
#' @param alt_col Name of the column identifying alternatives.
#' @param choice_col Name of the column indicating the chosen alternative
#'   (1 = chosen, 0 = not chosen).
#' @param covariate_cols Vector of names of structural covariate columns
#'   (the random-coefficient dimensions).
#' @param person_col Name of the respondent column grouping choice
#'   situations. `NULL` (default) makes each choice situation its own
#'   respondent.
#' @param alt_covariate_cols Names of alternative-level covariate columns
#'   (constant within each alternative) forming the \eqn{\delta} mean
#'   function. `NULL` (default) gives an intercept-only design (P = 1).
#' @param outside_opt_label Label of physical outside-option rows, removed
#'   when `include_outside_option = TRUE` (the outside good is implicit).
#' @param cf_residual_col Name of a first-stage residual column (control
#'   function for an endogenous covariate), appended to `X`. Default `NULL`.
#' @param include_outside_option Logical; if `TRUE` (default) an implicit
#'   outside option with systematic utility 0 is part of every choice set.
#' @param rc_dist Integer vector, one entry per column of `covariate_cols`:
#'   `0` for a normal random coefficient, `1` for log-normal (the
#'   coefficient enters utility as `exp(beta_ik)`; hierarchy normal on the
#'   log scale). Default `NULL` is all-normal. Automatically aligned through
#'   dropped columns; a `cf_residual_col` coordinate is always normal.
#' @returns A list of class `c("choicer_data_hmnl", "list")` containing:
#'   \itemize{
#'     \item `X`: Structural design matrix (total_rows x K_struct), no ASC
#'       columns; `cf_residual_col` last when supplied.
#'     \item `alt_of_row`: Integer alternative code per row (`1..J`).
#'     \item `alt_idx`: Alias of `alt_of_row` for the pooled-MLE init.
#'     \item `Z`: Alternative-level design (J x P), intercept first.
#'     \item `M`: Inside alternatives per choice situation.
#'     \item `choice_pos`: 1-based within-task position of the chosen row;
#'       `0` = outside option chosen.
#'     \item `Ti`: Choice situations per respondent.
#'     \item `person_ids`, `N_persons`, `n_tasks`, `J`, `K_struct`, `P`.
#'     \item `include_outside_option`: Logical flag.
#'     \item `alt_mapping`: Data.table mapping alternatives to summary
#'       statistics (outside option is `alt_int = 0`).
#'     \item `param_map`: Named list of index vectors (`beta`, `theta`),
#'       robust to collinearity drops.
#'     \item `rc_dist`: Integer vector aligned with the columns of `X`.
#'     \item `dropped_cols`, `dropped_z_cols`: Dropped column names, if any.
#'     \item `data_spec`: Column-name metadata (incl. `person_col`,
#'       `outside_opt_label`, `cf_residual_col`, `alt_covariate_cols`).
#'   }
#' @seealso [prepare_hmnp_data()] for the hierarchical probit counterpart.
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
#' dt[, quality := 0.1 * alt]                # alternative-level covariate
#' dt[, choice := 0L]
#' # leave some tasks all-zero: outside option chosen
#' dt[, choice := if (runif(1) < 0.8) sample(c(1L, rep(0L, J - 1))) else 0L,
#'    by = task]
#' input <- prepare_hmnl_data(dt, "task", "alt", "choice", c("x1", "x2"),
#'                            person_col = "pid",
#'                            alt_covariate_cols = "quality")
#' str(input$Z)
#' input$alt_mapping
#' @export
prepare_hmnl_data <- function(
    data,
    id_col,
    alt_col,
    choice_col,
    covariate_cols,
    person_col = NULL,
    alt_covariate_cols = NULL,
    outside_opt_label = NULL,
    cf_residual_col = NULL,
    include_outside_option = TRUE,
    rc_dist = NULL
) {
  ## rc_dist is validated against the user-supplied covariates BEFORE the
  ## shared prep so error messages refer to the input, then aligned to the
  ## surviving X columns afterwards (collinearity / task-constant drops).
  K_user <- length(covariate_cols)
  if (is.null(rc_dist)) rc_dist <- rep(0L, K_user)
  rc_dist <- as.integer(rc_dist)
  if (length(rc_dist) != K_user) {
    stop("`rc_dist` must have length ", K_user,
         " (one entry per column of `covariate_cols`).")
  }
  if (!all(rc_dist %in% c(0L, 1L))) {
    stop("`rc_dist` entries must be 0 (normal) or 1 (log-normal).")
  }

  panel <- .prepare_hb_panel(
    data = data, id_col = id_col, alt_col = alt_col,
    choice_col = choice_col, covariate_cols = covariate_cols,
    person_col = person_col, alt_covariate_cols = alt_covariate_cols,
    outside_opt_label = outside_opt_label, cf_residual_col = cf_residual_col,
    include_outside_option = include_outside_option
  )

  ## Align rc_dist through drops: X columns are the kept structural
  ## covariates (in covariate_cols order) plus, last, the cf residual —
  ## whose coefficient is always a normal coordinate.
  rc_map <- stats::setNames(rc_dist, covariate_cols)
  kept <- colnames(panel$X)
  rc_aligned <- integer(length(kept))
  names(rc_aligned) <- kept
  is_user_cov <- kept %in% covariate_cols
  rc_aligned[is_user_cov] <- rc_map[kept[is_user_cov]]
  panel$rc_dist <- rc_aligned

  structure(panel, class = c("choicer_data_hmnl", "list"))
}

#' Fit a hierarchical Bayesian multinomial logit (HMNL)
#'
#' Runs the adaptive RW-Metropolis-within-Gibbs sampler for the hierarchical
#' (random-coefficients, panel or cross-sectional) multinomial logit with a
#' BLP-style alternative-level random effect:
#' \deqn{U_{ijt} = x_{ijt}'\gamma_i + \delta_j + \epsilon_{ijt}, \qquad
#'       U_{iot} = \epsilon_{iot},}
#' with i.i.d. Gumbel shocks (including on the implicit outside option, whose
#' systematic utility is 0), \eqn{\beta_i \sim N(b, W)} over the structural
#' covariates (\eqn{\gamma_{ik} = \beta_{ik}} or \eqn{\exp(\beta_{ik})} per
#' `rc_dist`), and \eqn{\delta_j = z_j'\theta + \xi_j},
#' \eqn{\xi_j \sim N(0, \sigma_d^2)}. Partial pooling shrinks each
#' \eqn{\delta_j} toward its characteristics-based mean \eqn{z_j'\theta};
#' the outside option anchors the level of \eqn{\delta} (mean utility
#' relative to the outside good), so no base alternative or sum-to-zero
#' constraint is needed.
#'
#' \strong{Initialization.} \eqn{\beta_i} start at the pooled MNL maximum
#' likelihood estimate over the structural covariates (log-normal
#' coordinates transformed to the chain scale with a warn-and-clamp at
#' 0.05); \eqn{\delta} starts at shrunk log choice-share contrasts against
#' the outside option; \eqn{\theta} at the OLS regression of the initial
#' \eqn{\delta} on `Z`.
#'
#' \strong{Priors.} \eqn{b \sim N(b\_bar, A^{-1})}, \eqn{W \sim IW(\nu, V)},
#' \eqn{\theta \sim N(\theta\_bar, A_\theta^{-1})}, and
#' \eqn{\sigma_d \sim} half-Cauchy\eqn{(0, s_d)} via the Makalic-Schmidt
#' scale mixture (set `sd_prior$half_cauchy = FALSE` for a plain
#' \eqn{IG(c_0, d_0)} on \eqn{\sigma_d^2}).
#'
#' \strong{Endogeneity.} If a price-like covariate is endogenous (correlated
#' with \eqn{\xi_j}), supply a first-stage residual via `cf_residual_col`
#' (Petrin & Train 2010); posterior uncertainty does NOT propagate
#' first-stage estimation error.
#'
#' @inheritParams prepare_hmnl_data
#' @param data Data frame (convenience pathway). Supply either `data` (with
#'   the column names) or `input_data`, not both.
#' @param input_data A `choicer_data_hmnl` object from [prepare_hmnl_data()]
#'   (advanced pathway).
#' @param prior Named list overriding prior defaults: `b_bar` (0), `A`
#'   (0.01 I), `nu` (K + 3), `V` (nu I), `theta_bar` (0), `A_theta`
#'   (0.01 I), `sd_prior` (list(half_cauchy = TRUE, s_d = 1, c0 = 3,
#'   d0 = 3)).
#' @param mcmc Named list overriding MCMC defaults: `R` (10000), `burn`
#'   (R %/% 5), `thin` (1), `seed` (drawn via `sample.int()` so
#'   `set.seed()` governs), `trace` (0), `s_init` (2.38 / sqrt(K)),
#'   `accept_target` (0.234).
#' @param chains Number of independent chains (seeds offset by 1). Chain 1
#'   provides the reported draws; all chains feed the split-R-hat table.
#' @param keep_beta_i `"means"` (default) stores posterior means/SDs of the
#'   individual-level \eqn{\beta_i}; `"draws"` additionally stores the full
#'   (K, N, R_keep) draw cube (memory-guarded); `"none"` stores neither.
#' @param keep_data Logical; keep the prepared data on the fit (default
#'   `TRUE`, needed by post-estimation methods).
#' @returns A `choicer_hmnl` object (classed `c("choicer_hmnl",
#'   "choicer_hb")`) with posterior summaries (`coefficients`, `se`,
#'   `vcov` for \eqn{b}; `theta_summary`; `sigma_d2_summary`; `W_mean`;
#'   `delta` and `xi` quality-ladder tables; `beta_i`), the raw thinned
#'   `draws`, acceptance diagnostics in `accept`, the split-R-hat table in
#'   `rhat` (when `chains > 1`), and sampler metadata.
#' @seealso [prepare_hmnl_data()], [simulate_hmnl_data()],
#'   [recovery_table()], [rhat()]
#' @examples
#' \donttest{
#' sim <- simulate_hmnl_data(N = 100, T = 3, J = 4, seed = 42)
#' fit <- run_hmnlogit(sim$data, "task", "alt", "choice", c("x1", "x2"),
#'                     person_col = "pid", alt_covariate_cols = "z1",
#'                     mcmc = list(R = 500, burn = 200))
#' summary(fit)
#' coef(fit, component = "delta")
#' }
#' @export
run_hmnlogit <- function(
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
    rc_dist = NULL,
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
    input_list <- prepare_hmnl_data(
      data, id_col, alt_col, choice_col, covariate_cols,
      person_col = person_col, alt_covariate_cols = alt_covariate_cols,
      outside_opt_label = outside_opt_label,
      cf_residual_col = cf_residual_col,
      include_outside_option = include_outside_option, rc_dist = rc_dist
    )
  } else {
    if (!inherits(input_data, "choicer_data_hmnl")) {
      stop("'input_data' must be a choicer_data_hmnl object from ",
           "prepare_hmnl_data().")
    }
    input_list <- input_data
  }

  # v1 requires the outside option: it is the location anchor identifying the
  # level of delta (and theta_0). A no-outside mode (sum-to-zero xi, no
  # theta_0) is roadmapped, not built.
  if (!isTRUE(input_list$include_outside_option)) {
    stop("run_hmnlogit requires include_outside_option = TRUE: the implicit ",
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
  rc <- as.integer(input_list$rc_dist)

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
  s_init <- as.numeric(mcmc$s_init %||% (2.38 / sqrt(K)))
  accept_target <- as.numeric(mcmc$accept_target %||% 0.234)
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

  # --- Pooled MNL MLE init (structural coordinates) ---------------------------
  # The frequentist kernel shares the prep's encodings exactly: choice_pos is
  # 1-based within task with 0 = outside, and the implicit outside option is
  # switched on. use_asc = FALSE because delta is initialized separately.
  weights <- rep(1, n_tasks)
  eval_f <- function(theta) {
    mnl_loglik_gradient_parallel(
      theta, X, input_list$alt_idx, input_list$choice_pos, input_list$M,
      weights, use_asc = FALSE, include_outside_option = TRUE
    )
  }
  opt <- tryCatch(
    run_optimizer("nloptr", rep(0, K), eval_f),
    error = function(e) NULL
  )
  if (is.null(opt) || !all(is.finite(opt$par))) {
    warning("Pooled MNL MLE initialization failed; starting beta at 0.",
            call. = FALSE)
    gamma_pooled <- rep(0, K)
  } else {
    gamma_pooled <- opt$par
  }
  beta_pooled <- gamma_pooled
  if (any(rc == 1L)) {
    ln <- which(rc == 1L)
    clamped <- gamma_pooled[ln] < 0.05
    if (any(clamped)) {
      warning("Pooled MLE for log-normal coordinate(s) ",
              paste(colnames(X)[ln[clamped]], collapse = ", "),
              " clamped at 0.05 before the log transform.", call. = FALSE)
    }
    beta_pooled[ln] <- log(pmax(gamma_pooled[ln], 0.05))
  }

  # delta init: shrunk log choice-share contrasts vs the outside option
  # (add-half so never-chosen alternatives stay finite); theta init: OLS of
  # delta_init on Z.
  am <- input_list$alt_mapping
  n_choices_in <- am[am$alt_int > 0, ][["N_CHOICES"]]
  n_out <- am[am$alt_int == 0, ][["N_CHOICES"]]
  delta_init <- log((n_choices_in + 0.5) / (n_out + 0.5))
  theta_init <- as.numeric(qr.solve(Z, delta_init))

  # --- Run the Gibbs sampler ------------------------------------------------------
  run_chain <- function(chain_seed) {
    hmnl_gibbs(
      X = X, Z = Z, M = input_list$M, choice_pos = input_list$choice_pos,
      include_outside_option = TRUE, alt_of_row = input_list$alt_of_row,
      Ti = input_list$Ti, rc_dist = rc, beta_pooled = beta_pooled,
      delta_init = delta_init, theta_init = theta_init,
      b_bar = b_bar, A = A, nu = nu, V = V,
      theta_bar = theta_bar, A_theta = A_theta, sd_prior = sd_prior,
      R = R, burn = burn, thin = thin, seed = chain_seed,
      keep_beta_i = keep_code, s_init = s_init,
      accept_target = accept_target, trace = trace
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

  # --- Assemble the posterior object -------------------------------------------
  x_names <- colnames(X)
  z_names <- colnames(Z)
  alt_labels <- as.character(am[am$alt_int > 0, ][[input_list$data_spec$alt_col]])
  w_names <- character(K * (K + 1L) / 2L)
  idx <- 1L
  for (a1 in seq_len(K)) {
    for (a2 in seq_len(a1)) {
      w_names[idx] <- sprintf("W[%s,%s]", x_names[a1], x_names[a2])
      idx <- idx + 1L
    }
  }
  colnames(out$bdraw) <- x_names
  colnames(out$wdraw) <- w_names
  colnames(out$deltadraw) <- alt_labels
  colnames(out$thetadraw) <- z_names

  # split-R-hat over b, theta, sigma_d2 (all chains; single chain uses the
  # within-chain split)
  chain_blocks <- lapply(c(list(out), extra_chains), function(ch) {
    # sigma_d2draw arrives as an n x 1 matrix (arma::vec); flatten so the
    # rhat helper names it rather than inheriting NULL colnames.
    list(b = ch$bdraw, theta = ch$thetadraw,
         sigma_d2 = as.numeric(ch$sigma_d2draw))
  })
  rhat_tab <- tryCatch(.hb_rhat_table(chain_blocks), error = function(e) NULL)

  coefficients <- colMeans(out$bdraw)
  se <- apply(out$bdraw, 2, stats::sd)
  vcov <- stats::cov(out$bdraw)

  W_mean <- matrix(0, K, K, dimnames = list(x_names, x_names))
  w_means <- colMeans(out$wdraw)
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

  mean_acc_beta <- mean(out$accept_rate_beta)
  if (mean_acc_beta < 0.10 || mean_acc_beta > 0.60) {
    warning(sprintf(
      "Mean beta acceptance rate %.2f is outside (0.10, 0.60); consider a ",
      mean_acc_beta), "longer burn-in or adjusting mcmc$s_init.",
      call. = FALSE)
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

  new_choicer_hmnl(
    call = cl,
    coefficients = coefficients,
    se = se,
    vcov = vcov,
    theta_summary = build_bayes_coef_table(out$thetadraw),
    sigma_d2_summary = build_bayes_coef_table(
      matrix(out$sigma_d2draw, ncol = 1, dimnames = list(NULL, "sigma_d^2"))
    ),
    W_mean = W_mean,
    delta = delta_tab,
    xi = xi_tab,
    beta_i = beta_i,
    draws = list(
      b = out$bdraw, w_vech = out$wdraw, delta = out$deltadraw,
      theta = out$thetadraw, sigma_d2 = as.numeric(out$sigma_d2draw),
      loglik = as.numeric(out$loglik_trace)
    ),
    accept = list(
      beta = as.numeric(out$accept_rate_beta),
      delta = as.numeric(out$accept_rate_delta),
      s_final = as.numeric(out$s_final),
      s_delta_final = as.numeric(out$s_delta_final),
      mean_beta = mean_acc_beta,
      mean_delta = mean(out$accept_rate_delta)
    ),
    rhat = rhat_tab,
    prior = list(b_bar = b_bar, A = A, nu = nu, V = V, theta_bar = theta_bar,
                 A_theta = A_theta, sd_prior = sd_prior),
    mcmc = list(R = R, burn = burn, thin = thin, seed = seed, trace = trace,
                s_init = s_init, accept_target = accept_target,
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
    sampler = list(name = "hmnl_gibbs",
                   elapsed_time = convertTime(elapsed)),
    data = if (keep_data) input_list else NULL
  )
}
