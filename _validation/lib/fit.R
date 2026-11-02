# ============================================================================
# Shared fit helpers: the rep record, the estimates table (estimate, truth,
# SE by variance route), the common block of the covariance, joint Wald and
# LR statistics, choicer's post-estimation at the truth, fit diagnostics,
# multi-start checks of the global maximum, and run_unit(), the one place a
# unit (family, cell, method, rep) runs. Definitions only.
# Shared file: one owner at a time; the lead serializes changes.
# ============================================================================

# The variance routes of vcov(fit, type = ) ("robust" is vcov's name for the
# sandwich) and the est table's columns.
VS_ROUTES <- c("hessian", "bhhh", "robust", "cluster")
VS_EST_COLS <- c(name = "character", block = "character", est = "numeric", truth = "numeric",
                 se_hessian = "numeric", se_bhhh = "numeric", se_robust = "numeric",
                 se_cluster = "numeric", gated = "logical")
VS_DERIVED_COLS <- c(name = "character", est = "numeric", se = "numeric", truth = "numeric")

# --- The rep record ---------------------------------------------------------

.vs_check_cols <- function(df, cols, what) {
  if (!is.data.frame(df)) stop("rep record: `", what, "` must be a data.frame or NULL")
  miss <- setdiff(names(cols), names(df))
  if (length(miss)) stop("rep record: `", what, "` lacks column(s) ", paste(miss, collapse = ", "))
  for (k in names(cols)) {
    ok <- switch(cols[[k]], character = is.character(df[[k]]), numeric = is.numeric(df[[k]]),
                 logical = is.logical(df[[k]]) && !anyNA(df[[k]]))
    if (!ok) stop("rep record: `", what, "$", k, "` must be ", cols[[k]])
  }
  if (anyDuplicated(df$name)) stop("rep record: `", what, "$name` must be unique")
  df
}

.vs_check_test <- function(x, stat, what) {
  is_test <- function(t) is.list(t) && all(c(stat, "df", "p") %in% names(t))
  if (is.null(x) || is_test(x)) return(x)
  if (is.list(x) && !is.null(names(x)) && all(vapply(x, is_test, logical(1)))) return(x)
  stop("rep record: `", what, "` must be list(", stat, ", df, p), a named list of them, or NULL")
}

# The rep record (plan §6), every field present (NULL where absent):
#   list(family, cell, method, rep, seed, ok, error, warnings,
#        est     = data.frame(name, block, est, truth, se_hessian, se_bhhh,
#                             se_robust, se_cluster, gated),
#        vcov_c  = list(route, common, asc_var),     vcov_common()
#        joint_wald = list(W, df, p), lr = list(LR, df, p) (or named lists of them),
#        derived = data.frame(name, est, se, truth),
#        diag    = list(code, converged, iterations, decrement, time_sec,
#                       peak_mem_mb = c(heap_mb, rss_mb), rss_kind, multistart,
#                       sim_sec, unit_sec, messages),
#        meta    = list(choicer_version, threads, ref))
# `warnings` is always a data.frame(stage, class, message), with no rows when
# there were none; diag$messages is a named integer (a message's trimmed
# text -> its count; integer(0) when none); the other diag fields not given
# are NA (multistart NULL).
new_rep_record <- function(family, cell, method, rep, seed, ok = TRUE, error = NULL,
                           warnings = NULL, est = NULL, vcov_c = NULL, joint_wald = NULL,
                           lr = NULL, derived = NULL, diag = NULL, meta = NULL) {
  ids <- list(family = family, cell = cell, method = method)
  for (k in names(ids)) {
    if (!is.character(ids[[k]]) || length(ids[[k]]) != 1L || is.na(ids[[k]])) {
      stop("rep record: `", k, "` must be a string")
    }
  }
  if (!is.numeric(rep) || length(rep) != 1L || is.na(rep) || rep < 1 || rep != floor(rep)) {
    stop("rep record: `rep` must be a positive whole number")
  }
  if (!is.numeric(seed) || length(seed) != 1L) {
    stop("rep record: `seed` must be a number (NA allowed)")
  }
  if (!is.logical(ok) || length(ok) != 1L || is.na(ok)) {
    stop("rep record: `ok` must be TRUE or FALSE")
  }
  if (!is.null(error) && (!is.character(error) || length(error) != 1L)) {
    stop("rep record: `error` must be a string or NULL")
  }
  if (is.null(warnings)) warnings <- .vs_no_warnings()
  if (!is.data.frame(warnings) || !all(c("stage", "class", "message") %in% names(warnings))) {
    stop("rep record: `warnings` must be a data.frame(stage, class, message)")
  }
  if (!is.null(est)) est <- .vs_check_cols(est, VS_EST_COLS, "est")
  if (!is.null(derived)) derived <- .vs_check_cols(derived, VS_DERIVED_COLS, "derived")
  if (!is.null(vcov_c) &&
      !(is.list(vcov_c) && all(c("route", "common", "asc_var") %in% names(vcov_c)))) {
    stop("rep record: `vcov_c` must be list(route, common, asc_var), from vcov_common(), or NULL")
  }
  joint_wald <- .vs_check_test(joint_wald, "W", "joint_wald")
  lr <- .vs_check_test(lr, "LR", "lr")
  if (!is.null(diag) && !is.list(diag)) stop("rep record: `diag` must be a list or NULL")
  d <- list(code = NA_integer_, converged = NA, iterations = NA_integer_, decrement = NA_real_,
            time_sec = NA_real_, peak_mem_mb = c(heap_mb = NA_real_, rss_mb = NA_real_),
            rss_kind = NA_character_, multistart = NULL, sim_sec = NA_real_, unit_sec = NA_real_,
            messages = integer(0))
  for (k in names(diag)) d[k] <- list(diag[[k]])
  list(family = family, cell = cell, method = method, rep = as.integer(rep), seed = seed,
       ok = ok, error = error, warnings = warnings, est = est, vcov_c = vcov_c,
       joint_wald = joint_wald, lr = lr, derived = derived, diag = d,
       meta = meta %||% list(choicer_version = NA_character_, threads = NA_integer_, ref = NULL))
}

.vs_no_warnings <- function() {
  data.frame(stage = character(), class = character(), message = character(),
             stringsAsFactors = FALSE)
}

# --- Estimates and variances ------------------------------------------------

# Aligns a named numeric with the names `nm` (the fit's coefficients), by
# name; errors list the missing names and, unless extra_ok, the extra ones.
.vs_align <- function(x, nm, what, extra_ok = FALSE) {
  if (!is.numeric(x) || is.null(names(x))) stop(what, " must be a named numeric")
  miss <- setdiff(nm, names(x))
  extra <- if (extra_ok) character() else setdiff(names(x), nm)
  if (length(miss) || length(extra)) {
    stop(what, " does not match the fit's coefficients by name",
         if (length(miss)) paste0("; missing: ", paste(miss, collapse = ", ")),
         if (length(extra)) paste0("; extra: ", paste(extra, collapse = ", ")), call. = FALSE)
  }
  x[nm]
}

# The fit's covariance by route, a named list (NULL where a route failed,
# with a warning of class "vs_route_failed"). A route that is the fit's own
# se_method (the sandwich's being "robust") takes the stored fit$vcov, which
# is vcov(fit, type = se_method)'s; the others are recomputed.
fit_vcovs <- function(fit, routes = "hessian", cluster = NULL) {
  bad <- setdiff(routes, VS_ROUTES)
  if (length(bad)) {
    stop("unknown variance route(s) ", paste(bad, collapse = ", "), "; routes: ",
         paste(VS_ROUTES, collapse = ", "), " (\"robust\" is the sandwich)")
  }
  own <- c(hessian = "hessian", bhhh = "bhhh", sandwich = "robust",
           cluster = "cluster")[fit$se_method %||% "hessian"]
  nm <- names(coef(fit))
  out <- stats::setNames(vector("list", length(routes)), routes)
  for (r in routes) {
    stored <- !is.na(own) && r == own && !is.null(fit$vcov) && (r != "cluster" || is.null(cluster))
    out[r] <- list(tryCatch({
      V <- if (stored) fit$vcov
           else stats::vcov(fit, type = r, cluster = if (r == "cluster") cluster)
      if (is.null(V)) stop("the information is singular (vcov is NULL)")
      dimnames(V) <- list(nm, nm)
      V
    }, error = function(e) {
      warning(warningCondition(sprintf("vcov route '%s' failed: %s", r, conditionMessage(e)),
                               class = "vs_route_failed"))
      NULL
    }))
  }
  out
}

.vs_se <- function(V) {
  if (is.null(V)) return(NULL)
  v <- diag(V)
  if (any(v < 0, na.rm = TRUE)) {
    warning(warningCondition("negative variance(s) on the diagonal of the covariance",
                             class = "vs_negative_variance"))
  }
  sqrt(ifelse(v >= 0, v, NA_real_))
}

# The est table: one row per coefficient (coef(fit) order) with its block
# from fit$param_map, the estimate, the truth (a named numeric matched by
# name with coef(fit)), the SEs of every route (NA for routes not asked or
# failed) and the gate flag. `V` takes covariances already computed by
# fit_vcovs() (then `routes` and `cluster` are not used). `gated` is TRUE,
# FALSE, or block and parameter names.
est_table <- function(fit, truth, routes = "hessian", cluster = NULL, gated = TRUE, V = NULL) {
  cf <- coef(fit)
  nm <- names(cf)
  truth <- .vs_align(truth, nm, "truth")
  block <- rep("other", length(cf))
  for (b in names(fit$param_map)) block[fit$param_map[[b]]] <- b
  if (is.null(V)) V <- fit_vcovs(fit, routes, cluster)
  se <- lapply(VS_ROUTES, function(r) .vs_se(V[[r]]) %||% rep(NA_real_, length(cf)))
  names(se) <- paste0("se_", VS_ROUTES)
  out <- data.frame(name = nm, block = block, est = unname(cf), truth = unname(truth),
                    lapply(se, unname), gated = .vs_gated(gated, nm, block),
                    stringsAsFactors = FALSE)
  rownames(out) <- NULL
  out
}

.vs_gated <- function(gated, nm, block) {
  if (isTRUE(gated)) return(rep(TRUE, length(nm)))
  if (isFALSE(gated)) return(rep(FALSE, length(nm)))
  if (!is.character(gated)) stop("`gated` must be TRUE, FALSE, or block or parameter names")
  unknown <- setdiff(gated, c(block, nm))
  if (length(unknown)) {
    stop("`gated` names no block or parameter: ", paste(unknown, collapse = ", "))
  }
  block %in% gated | nm %in% gated
}

# The common block of the covariance (plan §6): the full covariance of every
# parameter outside fit$param_map$asc, and only the constants' variances, so
# that rep files stay small with thousands of constants. `V` is a covariance
# already computed; NULL computes the route's. A failed route gives NULLs.
vcov_common <- function(fit, route = "hessian", V = NULL) {
  if (is.null(V)) V <- fit_vcovs(fit, route)[[route]]
  if (is.null(V)) return(list(route = route, common = NULL, asc_var = NULL))
  nm <- names(coef(fit))
  dimnames(V) <- list(nm, nm)
  asc <- as.integer(fit$param_map$asc %||% integer())
  com <- setdiff(seq_along(nm), asc)
  list(route = route, common = V[com, com, drop = FALSE],
       asc_var = stats::setNames(diag(V)[asc], nm[asc]))
}

# Joint Wald statistic W = (est - truth)' V^-1 (est - truth) ~ chi^2_df, on
# the Jacobi-equilibrated V (scale-free). With a named `est`, `truth` (which
# may hold more names) and a named V are taken at est's names. NA (with a
# "vs_wald_failed" warning) when V is missing, not finite or singular.
joint_wald <- function(est, truth, V) {
  if (!is.null(names(est))) {
    truth <- .vs_align(truth, names(est), "truth", extra_ok = TRUE)
    if (!is.null(V) && !is.null(rownames(V))) V <- V[names(est), names(est), drop = FALSE]
  }
  d <- unname(est - truth)
  df <- length(d)
  W <- tryCatch({
    if (is.null(V) || anyNA(d) || !all(is.finite(V))) stop("no finite covariance")
    s <- 1 / sqrt(diag(V))
    if (!all(is.finite(s))) stop("a variance is not positive")
    sd <- s * d
    sum(sd * solve(V * outer(s, s), sd))
  }, error = function(e) {
    warning(warningCondition(paste("joint Wald:", conditionMessage(e)), class = "vs_wald_failed"))
    NA_real_
  })
  list(W = W, df = df, p = stats::pchisq(W, df, lower.tail = FALSE))
}

# Likelihood-ratio statistic LR = 2 (ll_u - ll_r) ~ chi^2_df under the
# restriction (Wilks). A negative LR (the restricted fit higher) is kept: it
# flags an optimizer failure.
lr_test <- function(ll_u, ll_r, df) {
  if (!is.numeric(df) || length(df) != 1L || df < 1) stop("lr_test(): `df` must be >= 1")
  LR <- 2 * (ll_u - ll_r)
  list(LR = LR, df = df, p = stats::pchisq(LR, df, lower.tail = FALSE))
}

# --- Derived estimands --------------------------------------------------------

# vech of Sigma = L L' (row by row: Sigma_11, Sigma_21, Sigma_22, ...; the
# diagonal only for an uncorrelated L), the covariance of the MXL's g = L
# eta, from choicer's sigma block (log diagonal), with delta-method SEs from
# V, the block's covariance (NA without one). A data.frame(name, est, se).
sigma_vech <- function(sigma, K_w, correlated, V = NULL) {
  L <- L_from_theta(unname(sigma), K_w, correlated)
  pij <- if (correlated) cbind(rep(seq_len(K_w), seq_len(K_w)), sequence(seq_len(K_w)))
         else cbind(seq_len(K_w), seq_len(K_w))
  Jac <- matrix(0, nrow(pij), nrow(pij))
  for (p in seq_len(nrow(pij))) {
    dL <- matrix(0, K_w, K_w)
    dL[pij[p, 1L], pij[p, 2L]] <- if (pij[p, 1L] == pij[p, 2L]) L[pij[p, 1L], pij[p, 1L]] else 1
    Jac[, p] <- (dL %*% t(L) + L %*% t(dL))[pij]
  }
  se <- if (is.null(V)) rep(NA_real_, nrow(pij))
        else sqrt(pmax(diag(Jac %*% V %*% t(Jac)), 0))
  data.frame(name = sprintf("Sigma_%d%d", pij[, 1L], pij[, 2L]), est = (L %*% t(L))[pij],
             se = se, stringsAsFactors = FALSE)
}

# --- Post-estimation at the truth -------------------------------------------

# A copy of `fit` at theta0 (matched by name with coef(fit)), for derived
# estimands' truth: choicer's post-estimation at theta0 (plan §6), with the
# mixed logit's draws multiplied, S_mult times the fit's S (the prediction
# sites read draws_info$S in both draw modes). Its variance is removed and
# marked failed (`variance_error`), so that wtp() and summary() do not
# compute a Hessian at theta0 (their SEs are NA); loglik is NA; the NL's
# `lambda` and the MXL's `sigma` follow theta0.
at_truth <- function(fit, theta0, S_mult = 4) {
  if (!is.numeric(S_mult) || length(S_mult) != 1L || is.na(S_mult) || S_mult < 1 ||
      S_mult != floor(S_mult)) {
    stop("at_truth(): `S_mult` must be a positive whole number")
  }
  nm <- names(fit$coefficients)
  fit$coefficients[] <- .vs_align(theta0, nm, "theta0")
  p <- length(nm)
  fit$vcov <- NULL
  fit$se <- stats::setNames(rep(NA_real_, p), nm)
  fit$variance_error <- "at_truth(): the truth carries no variance"
  fit$loglik <- NA_real_
  pm <- fit$param_map
  if (inherits(fit, "choicer_nl") && !is.null(pm$lambda)) {
    fit$lambda <- fit$coefficients[pm$lambda]
  }
  if (inherits(fit, "choicer_mxl")) {
    if (!is.null(pm$sigma) && !is.null(fit$sigma)) {
      build_var_mat <- utils::getFromNamespace("build_var_mat", "choicer")
      sig <- build_var_mat(unname(fit$coefficients[pm$sigma]), fit$draws_info$K_w,
                           fit$rc_correlation)
      dimnames(sig) <- dimnames(fit$sigma)
      fit$sigma <- sig
    }
    fit$draws_info$S <- as.integer(S_mult * fit$draws_info$S)
  }
  fit$at_truth <- list(theta0 = fit$coefficients, S_mult = S_mult)
  fit
}

# --- Diagnostics ------------------------------------------------------------

# The optimizer's record: code (fit$convergence), converged (choicer's
# reading of the code for its optimizer), iterations, the polish's decrement
# (over the mean weight; NA without a polish) and time_sec (the optimizer
# and the polish, fit$optimizer$elapsed_time).
fit_diag <- function(fit) {
  is_converged <- utils::getFromNamespace(".is_converged", "choicer")
  list(code = fit$convergence %||% NA_integer_,
       converged = is_converged(fit),
       iterations = fit$optimizer$iterations %||% NA_integer_,
       decrement = fit$optimizer$polish$decrement %||% NA_real_,
       time_sec = fit$optimizer$elapsed_time %||% NA_real_)
}

# Every family's fits: evaluates `expr`, a choicer fit (polish = "newton"
# where the family polishes), and returns it, or an error unless it
# converged: with a polish, the polish certified the optimum (its decrement
# below tolerance); without, choicer's reading of the optimizer's code
# (fit_diag()). The error makes run_unit() record a failed rep. nloptr's
# code alone cannot tell: NLOPT_FAILURE (-1) returns the start values
# without a warning (choicer's NL from its default start on D2's WESML
# samples, about half the fits), and it also comes at the optimum (an MNL in
# F8, decrement 1e-14 there), which the polish certifies. The polish's
# message, whose text carries the decrement, is muffled (fit_diag() keeps
# the decrement); other messages pass.
require_converged <- function(expr, what = "the fit") {
  fit <- withCallingHandlers(expr, message = function(m) {
    if (startsWith(conditionMessage(m), "Newton polish:")) invokeRestart("muffleMessage")
  })
  pol <- fit$optimizer$polish
  ok <- if (is.null(pol)) fit_diag(fit)$converged else pol$converged
  if (!isTRUE(ok)) {
    stop(what, " did not converge: code ", fit$convergence, ", ", fit$message, call. = FALSE)
  }
  fit
}

# Every mixed logit fit of the suite (S1, 2026-11-02): fit_fun(theta_init)
# runs one fit (theta_init = NULL: choicer's default start); here it runs
# from the default start and from theta0, each through require_converged(),
# and the certified fit with the higher log-likelihood is returned, an error
# (the rep fails) when neither is certified. Why: from choicer's default
# start, the mixed logit fits of D1 and D4 stopped on plateaus where an L
# diagonal goes to -Inf (a random coefficient's SD to 0, or Sigma singular),
# which the polish certifies (on D1 up to 14 log-likelihood units below the
# maximum), and some did not converge at all; from theta0 none did. The record
# attr(fit, "start") = list(winner = "default" | "theta0", dlogL =
# loglik_default - loglik_theta0 (NA unless both are certified),
# default_error, theta0_error (NULL or the message)) measures how often the
# default start loses: a package finding. Methods put it in diag$start.
fit_two_starts <- function(fit_fun, theta0, what = "the fit") {
  run <- function(start, label) {
    tryCatch(require_converged(fit_fun(start), sprintf("%s (%s start)", what, label)),
             error = function(e) conditionMessage(e))
  }
  fits <- list(default = run(NULL, "default"), theta0 = run(theta0, "theta0"))
  failed <- vapply(fits, is.character, logical(1))
  if (all(failed)) {
    stop(what, " did not converge from either start: ", fits$default, "; ", fits$theta0,
         call. = FALSE)
  }
  ll <- vapply(fits, function(f) if (is.character(f)) NA_real_ else f$loglik, numeric(1))
  winner <- names(fits)[which.max(ll)]
  fit <- fits[[winner]]
  attr(fit, "start") <- list(winner = winner, dlogL = unname(ll[["default"]] - ll[["theta0"]]),
                             default_error = if (failed[["default"]]) fits$default,
                             theta0_error = if (failed[["theta0"]]) fits$theta0)
  fit
}

# Refits from n starts theta_hat + spread * SE * U(-1, 1) through the
# caller's refit(theta_init), which returns a fit (MXL and NL; the MNL takes
# no theta_init). SEs that are not finite are replaced by 0.1, and SEs above
# 1 capped at 1: a fit on a zero-variance plateau (an MXL's log SD driven to
# -Inf) has SEs of order 1e4 there, which sent both starts out of range and
# made both refits fail, so the check missed exactly the fits it exists to
# catch. The NL's dissimilarities move on the log scale (lambda_hat *
# exp(spread * SE / lambda_hat * U)), so they stay positive. Returns dlogL =
# max over starts of loglik_alt - loglik_hat (> 0: the fit missed a higher
# maximum), max_dtheta_se = max |theta_alt - theta_hat| / SE (the capped
# SE), n_failed (errors or a loglik that is not finite) and n_warned (refits
# that warned; their warnings are muffled so that they do not count as the
# rep's).
multistart <- function(fit, refit, n = 2, seed, spread = 2) {
  th <- coef(fit)
  p <- length(th)
  se <- fit$se %||% rep(NA_real_, p)
  se[!is.finite(se) | se <= 0] <- 0.1
  se <- pmin(se, 1)
  U <- vs_with_seed(seed, matrix(stats::runif(n * p, -1, 1), n, p))
  lam <- as.integer(fit$param_map$lambda %||% integer())
  res <- lapply(seq_len(n), function(k) {
    start <- th + spread * se * U[k, ]
    start[lam] <- th[lam] * exp(spread * se[lam] / th[lam] * U[k, lam])
    warned <- FALSE
    f <- tryCatch(withCallingHandlers(refit(start), warning = function(w) {
      warned <<- TRUE
      invokeRestart("muffleWarning")
    }), error = function(e) NULL)
    ll <- if (is.null(f)) NA_real_ else f$loglik
    if (!is.finite(ll)) return(c(dlogL = NA, dtheta = NA, warned = warned))
    c(dlogL = ll - fit$loglik, dtheta = max(abs(coef(f)[names(th)] - th) / se),
      warned = warned)
  })
  res <- do.call(rbind, res)
  ok <- !is.na(res[, "dlogL"])
  list(dlogL = if (any(ok)) max(res[ok, "dlogL"]) else NA_real_,
       max_dtheta_se = if (any(ok)) max(res[ok, "dtheta"]) else NA_real_,
       n_failed = sum(!ok), n_warned = sum(res[, "warned"] > 0))
}

# --- Running a unit ---------------------------------------------------------

# Runs one unit (family, cell, method, rep) and returns its rep record; it
# never throws. The seed is vs_seed(family$id, cell, rep), the same for every
# method of the rep; `sim` is the rep's simulated data when the runner caches
# it across methods (else simulate() runs here, after vs_set_seed(seed)).
# The method runs after vs_set_seed(seed) too, so its own draws do not depend
# on whether `sim` was cached. Every warning is recorded (stage, class,
# message) and muffled; every message is muffled too, choicer's timing lines
# ("Optimization run time", "MCMC run time") dropped and the others counted
# by trimmed text in diag$messages. An error, in simulate() or the method,
# gives a record with ok = FALSE and the message. A method returns a list of
# est, vcov_c, joint_wald, lr, derived and diag (missing ones NULL); the
# identity, the times (sim_sec, unit_sec: the method's wall time), the peak
# memory, the messages and meta are filled here.
run_unit <- function(family, cell, method, rep, sim = NULL, ref = NULL) {
  log <- new.env()
  log$w <- list()
  log$m <- character()
  record_conditions <- function(stage, expr) {
    withCallingHandlers(expr, warning = function(w) {
      log$w[[length(log$w) + 1L]] <- data.frame(
        stage = stage, class = class(w)[1], message = conditionMessage(w),
        stringsAsFactors = FALSE)
      invokeRestart("muffleWarning")
    }, message = function(m) {
      txt <- trimws(conditionMessage(m))
      if (!grepl("^(Optimization|MCMC) run time", txt)) log$m <- c(log$m, txt)
      invokeRestart("muffleMessage")
    })
  }
  # Distinct message texts with their counts, in order of first appearance.
  message_counts <- function() {
    u <- unique(log$m)
    if (!length(u)) return(integer(0))
    stats::setNames(tabulate(match(log$m, u), length(u)), u)
  }
  now <- function() proc.time()[["elapsed"]]
  fid <- if (is.list(family) && is.character(family$id)) family$id[1] else NA_character_
  seed <- NA_integer_
  diag <- list()
  loaded <- isNamespaceLoaded("choicer")
  meta <- list(
    choicer_version = if (loaded) getNamespaceVersion("choicer")[["version"]] else NA_character_,
    threads = if (loaded) choicer::thread_info()$omp_get_max_threads else NA_integer_,
    ref = ref)
  record <- function(ok, error = NULL, out = list()) {
    w <- if (length(log$w)) do.call(rbind, log$w) else .vs_no_warnings()
    d <- c(diag, list(messages = message_counts()))
    new_rep_record(family = fid, cell = cell, method = method, rep = rep, seed = seed, ok = ok,
                   error = error, warnings = w, est = out$est, vcov_c = out$vcov_c,
                   joint_wald = out$joint_wald, lr = out$lr, derived = out$derived,
                   diag = c(out$diag[setdiff(names(out$diag), names(d))], d), meta = meta)
  }
  tryCatch({
    row <- family_cell(family, cell)
    fn <- family$methods[[method]]
    if (!is.function(fn)) {
      stop("family ", fid, " has no method '", method, "'; methods: ",
           paste(names(family$methods), collapse = ", "))
    }
    seed <- vs_seed(fid, cell, rep)
    if (is.null(sim)) {
      t0 <- now()
      vs_set_seed(seed)
      sim <- tryCatch(record_conditions("simulate", family$simulate(row, seed)),
                      error = function(e) structure(list(msg = conditionMessage(e)),
                                                    class = "vs_sim_failed"))
      diag$sim_sec <- now() - t0
      if (inherits(sim, "vs_sim_failed")) return(record(FALSE, paste("simulate:", sim$msg)))
    }
    vs_set_seed(seed)
    mem0 <- vs_mem_reset()
    t0 <- now()
    out <- tryCatch(record_conditions("method", fn(sim, row, seed)), error = function(e) e)
    diag$unit_sec <- now() - t0
    mem <- vs_mem_peak(mem0)
    diag$peak_mem_mb <- c(heap_mb = mem$heap_mb, rss_mb = mem$rss_mb)
    diag$rss_kind <- mem$rss_kind
    if (inherits(out, "error")) return(record(FALSE, conditionMessage(out)))
    if (is.null(out)) out <- list()
    extra <- setdiff(names(out), c("est", "vcov_c", "joint_wald", "lr", "derived", "diag"))
    if (!is.list(out) || length(extra)) {
      stop("method '", method, "' must return list(est, vcov_c, joint_wald, lr, derived, diag)",
           if (length(extra)) paste0("; unknown field(s): ", paste(extra, collapse = ", ")))
    }
    tryCatch(record(TRUE, out = out), error = function(e) {
      record(FALSE, paste0("invalid result of method '", method, "': ", conditionMessage(e)))
    })
  }, error = function(e) {
    tryCatch(record(FALSE, conditionMessage(e)), error = function(e2) {
      list(family = fid, cell = as.character(cell)[1], method = as.character(method)[1],
           rep = rep, seed = seed, ok = FALSE,
           error = paste(conditionMessage(e), "|", conditionMessage(e2)),
           warnings = .vs_no_warnings(), est = NULL, vcov_c = NULL, joint_wald = NULL,
           lr = NULL, derived = NULL, diag = diag, meta = meta)
    })
  })
}

# A record's message counts with the simulation's messages added (texts in
# order of first appearance, the simulation's first).
.vs_message_counts <- function(sim_txt, counts) {
  u <- unique(c(sim_txt, names(counts)))
  if (!length(u)) return(integer(0))
  n <- tabulate(match(sim_txt, u), length(u))
  k <- match(names(counts), u)
  n[k] <- n[k] + as.integer(counts)
  stats::setNames(n, u)
}

# The rep's records, a list named by method. The rep is simulated once,
# after vs_set_seed(seed), and the same data go to every method through
# run_unit(), which seeds again before each method, so that the records
# equal those run_unit() gives when it simulates itself. The simulation's
# warnings (stage "simulate"), messages and time join each record; a failed
# simulation gives every method a failed record ("simulate: <message>").
# Every method receives the same `sim`: a method must not modify sim$data by
# reference (`:=`, setorder(), setkey()), or the later methods fit other data.
run_rep <- function(fam, cell, rep, methods, ref = NULL) {
  seed <- vs_seed(fam$id, cell, rep)
  w <- list()
  txt <- character()
  t0 <- proc.time()[["elapsed"]]
  vs_set_seed(seed)
  sim <- tryCatch(withCallingHandlers(
    fam$simulate(family_cell(fam, cell), seed),
    warning = function(c) {
      w[[length(w) + 1L]] <<- data.frame(stage = "simulate", class = class(c)[1],
                                         message = conditionMessage(c), stringsAsFactors = FALSE)
      invokeRestart("muffleWarning")
    }, message = function(c) {
      m <- trimws(conditionMessage(c))
      if (!grepl("^(Optimization|MCMC) run time", m)) txt <<- c(txt, m)
      invokeRestart("muffleMessage")
    }), error = function(e) structure(list(msg = conditionMessage(e)), class = "vs_sim_failed"))
  sim_sec <- proc.time()[["elapsed"]] - t0
  sim_w <- if (length(w)) do.call(rbind, w) else NULL
  failed <- inherits(sim, "vs_sim_failed")
  meta <- list(choicer_version = getNamespaceVersion("choicer")[["version"]],
               threads = choicer::thread_info()$omp_get_max_threads, ref = ref)
  recs <- lapply(methods, function(m) {
    if (failed) {
      return(new_rep_record(fam$id, cell, m, rep, seed, ok = FALSE,
                            error = paste("simulate:", sim$msg), warnings = sim_w,
                            diag = list(sim_sec = sim_sec,
                                        messages = .vs_message_counts(txt, integer(0))),
                            meta = meta))
    }
    rec <- run_unit(fam, cell, m, rep, sim = sim, ref = ref)
    rec$warnings <- rbind(sim_w, rec$warnings)
    rec$diag$sim_sec <- sim_sec
    rec$diag$messages <- .vs_message_counts(txt, rec$diag$messages)
    rec
  })
  names(recs) <- methods
  recs
}
