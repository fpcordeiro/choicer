# ============================================================================
# Phase 8C batteries (the Newton phase of run_mnlogit() and run_nestlogit(),
# polish =, and the MNL/NL gradient kernels' overflow flag), sourced by
# refactor_identity_check.R after p8b_batteries.R.
#   * p8c_grad  : mnl_/nl_loglik_gradient_parallel() at ordinary and
#                 non-finite points (a chosen utility at -Inf, an unchosen one
#                 at +Inf, an unchosen nest at +Inf, a NaN and a tiny
#                 dissimilarity), +- outside option, +- opg_diag. From the
#                 commit that adds it, the lists carry `overflow`: compare
#                 with the tooling's p8c/compare_golden.R, which admits that
#                 leaf only.
#   * p8c_polish: run_mnlogit() and run_nestlogit() with polish = "newton" and
#                 "bhhh" from a loose stop (maxeval 6): outside option and
#                 base, scale_vars none/sd/bhhh, weights, an alternative never
#                 chosen (its constant held; for NL also in a nest whose
#                 dissimilarity exceeds 1), the reference never chosen with a
#                 chosen constant held, every se_method, NL dissimilarities
#                 near and far from 1 and at the bound, singleton nests, a
#                 nest no situation offers twice, optim (and BFGS, which drops
#                 the bound), a custom optimizer, keep_data = FALSE,
#                 max_iter = 0, a polish that fails on the memory check; plus
#                 polish = "none". Class D against main (unused argument).
#   * p8c_mxl   : run_mxlogit()'s polish on the cases that move into shared
#                 helpers: the constants diverging together (outside option
#                 and base, never and always chosen, with and without a bound
#                 or a held constant that stops them), every se_method x
#                 polish pair, a custom optimizer, keep_data = FALSE, and a
#                 polish that fails on the memory check. Class I.
# ============================================================================

# Long data with choices drawn from the model at a known parameter, as
# p8_frame(), with the dissimilarities given (one per nest of two or more
# alternatives, or one value for all) and, with `alone`, alternatives no
# situation offers together.
p8c_frame <- function(seed, N, J, lo, hi, ioo, nests = NULL, lambda = 0.6,
                      alone = NULL) {
  set.seed(seed)
  M <- sample(lo:hi, N, replace = TRUE)
  pool <- setdiff(seq_len(J), alone)
  alt <- unlist(lapply(M, function(m) {
    if (length(alone)) sort(c(sample(pool, m - 1L), sample(alone, 1L)))
    else sort(sample.int(J, m))
  }))
  id <- rep(seq_len(N), M)
  X <- matrix(rnorm(sum(M) * 2L), ncol = 2L)
  delta <- c(if (!ioo) 0, rnorm(J - (!ioo), sd = 0.5))
  beta <- c(0.8, -0.5)
  start <- c(0L, cumsum(M))
  choice <- integer(length(alt))
  if (is.null(nests)) {
    V <- drop(X %*% beta) + delta[alt]
    for (i in seq_len(N)) {
      v <- c(if (ioo) 0, V[(start[i] + 1L):start[i + 1L]])
      p <- exp(v - max(v)); p <- p / sum(p)
      k <- sample.int(length(p), 1L, prob = p) - ioo
      if (k >= 1L) choice[start[i] + k] <- 1L
    }
  } else {
    K_l <- sum(tabulate(nests) > 1L)
    theta <- c(beta, rep_len(lambda, K_l), if (ioo) delta else delta[-1L])
    pr <- choicer:::nl_predict(theta, X, as.integer(alt), as.integer(M),
                               as.integer(nests), TRUE, ioo)$choice_prob
    for (i in seq_len(N)) {
      p <- pr[(start[i] + 1L):start[i + 1L]]
      p <- c(if (ioo) max(0, 1 - sum(p)), p)
      k <- sample.int(length(p), 1L, prob = p) - ioo
      if (k >= 1L) choice[start[i] + k] <- 1L
    }
  }
  w_i <- runif(N, 0.5, 2)
  dt <- data.table::data.table(id = id, alt = alt, choice = choice,
                               x1 = X[, 1L], x2 = X[, 2L], w = w_i[id],
                               clus = (id - 1L) %/% 15L)
  if (!is.null(nests)) dt[, nest := nests[alt]]
  dt[]
}

# A memory limit below the polish's curvature and at or above what the
# gradient and the variance's BHHH matrix need (one thread), by bisection on
# the curvature kernel: so that the polish fails on its first curvature and
# the variance does not. NA where no such limit exists.
p8c_memory_limit <- function(curv, others) {
  needs <- function(f, lim) {
    old <- options(choicer.max_memory = lim)
    on.exit(options(old))
    inherits(tryCatch(f(), error = function(e) e), "error")
  }
  lo <- 1; hi <- 2^40
  if (!needs(curv, lo) || needs(curv, hi)) return(NA_real_)
  for (i in 1:80) {
    mid <- floor((lo + hi) / 2)
    if (needs(curv, mid)) lo <- mid else hi <- mid
    if (hi - lo <= 1) break
  }
  lim <- lo  # the largest limit the curvature kernel refuses
  if (any(vapply(others, needs, NA, lim))) NA_real_ else lim
}

# A fit's record without its call, for fits whose call holds a closure (a
# custom optimizer passed by value), whose environment no other process's
# copy is identical() to.
p8c_rec_nocall <- function(expr) {
  r <- p5a_rec(expr)
  if (is.list(r$value)) r$value$call <- NULL
  r
}

p8c_with_limit <- function(lim, expr) {
  if (is.na(lim)) return(list(value = NULL, warnings = character(),
                              messages = character(), error = "no limit"))
  old <- options(choicer.max_memory = lim)
  on.exit(options(old))
  p5a_rec(expr)
}

# ---------------------------------------------------------------------------
# p8c_grad
# ---------------------------------------------------------------------------

run_p8c_grad <- function() {
  out <- list()
  K <- 2L; J <- 30L
  nest_idx <- as.integer(c(rep(1:8, each = 3L), 9:14))
  for (model in c("mnl", "nl")) {
    for (ioo in c(FALSE, TRUE)) {
      seed <- 15000L + 10L * (model == "nl") + ioo
      d0 <- p4_design(seed, 40L, J, 3L, 8L, K, ioo, FALSE)
      th0 <- if (model == "mnl") {
        p8_theta(seed + 100L, K, n_delta_free(J, TRUE, ioo))
      } else {
        p8_nl_theta(seed + 100L, K, nest_idx, n_delta_free(J, TRUE, ioo))
      }
      t <- 3L
      if (d0$choice_idx[t] == 0L) d0$choice_idx[t] <- 1L
      row0 <- c(0L, cumsum(d0$M))
      rows_t <- (row0[t] + 1L):row0[t + 1L]
      chosen_row <- row0[t] + d0$choice_idx[t]
      other_row <- row0[t] + (if (d0$choice_idx[t] == 1L) 2L else 1L)
      # an unchosen row in another nest than the chosen one (NL), if any
      nest_row <- nest_idx[d0$alt_idx[rows_t]]
      far <- rows_t[nest_row != nest_idx[d0$alt_idx[chosen_row]]]
      cells <- list(
        ordinary = function(d, th) list(d, th),
        chosen_minus_inf = function(d, th) {
          d$X[chosen_row, 1L] <- -1e308; th[1L] <- 2; list(d, th)
        },
        unchosen_plus_inf = function(d, th) {
          d$X[other_row, 1L] <- 1e308; th[1L] <- 2; list(d, th)
        },
        # finite utilities whose scores overflow when summed: an unchosen
        # row of each of the first ten situations at 1e308 times 1e-308
        huge_score = function(d, th) {
          for (s in 1:10) {
            r <- row0[s] + (if (d$choice_idx[s] == 1L) 2L else 1L)
            d$X[r, 2L] <- 1e308
          }
          th[2L] <- 1e-308
          list(d, th)
        })
      if (model == "nl") {
        if (length(far)) {
          cells$unchosen_nest_plus_inf <- function(d, th) {
            d$X[far[1L], 1L] <- 1e308; th[1L] <- 2; list(d, th)
          }
        }
        cells$nan_lambda <- function(d, th) { th[K + 1L] <- NaN; list(d, th) }
        cells$tiny_lambda <- function(d, th) { th[K + 1L] <- 1e-300; list(d, th) }
      }
      for (nm in names(cells)) {
        e <- cells[[nm]](d0, th0)
        d <- e[[1L]]; th <- e[[2L]]
        for (od in c(FALSE, TRUE)) {
          key <- sprintf("grad|%s|%s|oo=%d|opg=%d", model, nm, ioo, od)
          out[[key]] <- p2_val(if (model == "mnl") {
            choicer:::mnl_loglik_gradient_parallel(th, d$X, d$alt_idx, d$choice_idx,
                                                   d$M, d$weights, TRUE, ioo, od)
          } else {
            choicer:::nl_loglik_gradient_parallel(th, d$X, d$alt_idx, d$choice_idx,
                                                  nest_idx, d$M, d$weights, TRUE,
                                                  ioo, od)
          })
        }
      }
    }
  }
  out
}

# ---------------------------------------------------------------------------
# p8c_polish
# ---------------------------------------------------------------------------

run_p8c_polish <- function() {
  out <- list()
  loose <- list(maxeval = 6L)
  custom <- function(theta_init, eval_f, lower, upper, control) {
    r <- nloptr::nloptr(theta_init, eval_f, lb = lower, ub = upper,
                        opts = list(algorithm = "NLOPT_LD_LBFGS", maxeval = 5L))
    list(par = r$solution, value = r$objective, convergence = r$status)
  }

  # --- MNL -----------------------------------------------------------------
  mnl_fit <- function(dt, ...) {
    run_mnlogit(dt, "id", "alt", "choice", c("x1", "x2"), control = loose, ...)
  }
  dm <- p8c_frame(16001L, 150L, 8L, 3L, 6L, FALSE)
  dm_o <- p8c_frame(16002L, 150L, 8L, 3L, 6L, TRUE)
  # an alternative never chosen
  dm_n <- data.table::copy(dm)
  dm_n[alt == 8L, choice := 0L]
  dm_n <- dm_n[id %in% dm_n[choice == 1L, id]]
  # the outside option never chosen, and an alternative chosen wherever it is
  # offered, whose constant held stops the others' rise
  dm_r <- dm_o[id %in% dm_o[choice == 1L, id]]
  dm_r <- dm_r[!(id %in% dm_r[alt == 8L, id] & !id %in% dm_r[alt == 8L & choice == 1L, id])]
  mnl_designs <- list(
    cross = list(dm), outside = list(dm_o, include_outside_option = TRUE),
    none = list(dm, scale_vars = "none"), sd = list(dm, scale_vars = "sd"),
    weights = list(dm, weights_col = "w"),
    held = list(dm_n), ref_held = list(dm_r, include_outside_option = TRUE),
    ref_free = list(dm_o[id %in% dm_o[choice == 1L, id]],
                    include_outside_option = TRUE),
    se_bhhh = list(dm, se_method = "bhhh"),
    se_sandwich = list(dm, weights_col = "w", se_method = "sandwich"),
    se_cluster = list(dm, cluster_col = "clus"),
    noasc = list(dm, use_asc = FALSE),
    optim = list(dm, optimizer = "optim", control = list(maxit = 3L)),
    custom = list(dm, optimizer = custom),
    nodata = list(dm, keep_data = FALSE),
    cert0 = list(dm, polish_control = list(max_iter = 0L)))
  for (nm in names(mnl_designs)) {
    for (pm in c("newton", "bhhh")) {
      a <- mnl_designs[[nm]]
      if (!is.null(a$control)) {
        out[[paste("polish|mnl", nm, pm, sep = "|")]] <- p5a_rec(do.call(
          run_mnlogit, c(list(a[[1L]], "id", "alt", "choice", c("x1", "x2")),
                         a[-1L], list(polish = pm))))
      } else {
        rec_f <- if (nm == "custom") p8c_rec_nocall else p5a_rec
        out[[paste("polish|mnl", nm, pm, sep = "|")]] <-
          rec_f(do.call(mnl_fit, c(a, list(polish = pm))))
      }
    }
  }
  out[["polish|mnl|cross|none"]] <- p5a_rec(mnl_fit(dm, polish = "none"))
  # a polish that fails on its curvature's memory check
  fm <- suppressWarnings(suppressMessages(mnl_fit(dm, se_method = "bhhh")))
  d <- fm$data; th <- unname(coef(fm))
  lim <- p8c_memory_limit(
    function() choicer:::mnl_loglik_hessian_parallel(th, d$X, d$alt_idx, d$choice_idx,
                                                     d$M, d$weights, TRUE, FALSE),
    list(function() choicer:::mnl_loglik_gradient_parallel(th, d$X, d$alt_idx,
                                                           d$choice_idx, d$M, d$weights,
                                                           TRUE, FALSE, TRUE),
         function() choicer:::mnl_bhhh_parallel(th, d$X, d$alt_idx, d$choice_idx,
                                                d$M, d$weights, TRUE, FALSE)))
  out[["polish|mnl|memory|newton"]] <- c(list(limit = lim), p8c_with_limit(
    lim, mnl_fit(dm, se_method = "bhhh", polish = "newton")))

  # --- NL ------------------------------------------------------------------
  nests <- c(1L, 1L, 1L, 2L, 2L, 3L, 4L, 4L)
  nl_fit <- function(dt, ...) {
    run_nestlogit(dt, "id", "alt", "choice", c("x1", "x2"), nest_col = "nest",
                  control = loose, ...)
  }
  dn <- p8c_frame(16101L, 200L, 8L, 3L, 6L, FALSE, nests)
  dn_o <- p8c_frame(16102L, 200L, 8L, 3L, 6L, TRUE, nests)
  dn_far <- p8c_frame(16103L, 300L, 8L, 4L, 8L, FALSE, nests, lambda = 0.25)
  dn_tiny <- p8c_frame(16104L, 30L, 8L, 3L, 6L, FALSE, nests, lambda = 0.05)
  # six singleton nests and one pair
  dn_single <- data.table::copy(dn)[, nest := c(1:6, 7L, 7L)[alt]]
  dn_alone <- p8c_frame(16105L, 200L, 8L, 3L, 6L, FALSE, nests, alone = 7:8)
  dn_n <- data.table::copy(dn)
  dn_n[alt == 5L, choice := 0L]
  dn_n <- dn_n[id %in% dn_n[choice == 1L, id]]
  # an alternative never chosen in the nest whose dissimilarity exceeds 1
  data(mode_choice, package = "choicer", envir = environment())
  mc <- data.table::as.data.table(mode_choice)
  b2 <- mc[mode == "bus"]
  b2[, `:=`(mode = "bus2", choice = 0L, gcost = gcost + 5)]
  mc2 <- rbind(mc, b2)[order(id, mode)]
  mc2[, nest := ifelse(mode == "car", "private", "public")]
  nl_designs <- list(
    cross = list(dn), outside = list(dn_o, include_outside_option = TRUE),
    bhhh_scale = list(dn, scale_vars = "bhhh"), sd = list(dn, scale_vars = "sd"),
    weights = list(dn, weights_col = "w"),
    far = list(dn_far), tiny = list(dn_tiny), single = list(dn_single),
    alone = list(dn_alone), held = list(dn_n),
    se_numeric = list(dn, se_method = "numeric"),
    se_bhhh = list(dn, se_method = "bhhh"),
    se_sandwich = list(dn, weights_col = "w", se_method = "sandwich"),
    se_cluster = list(dn, cluster_col = "clus"),
    optim = list(dn, optimizer = "optim", control = list(maxit = 3L)),
    bfgs = list(dn, optimizer = "optim", control = list(method = "BFGS", maxit = 3L)),
    custom = list(dn, optimizer = custom),
    nodata = list(dn, keep_data = FALSE),
    cert0 = list(dn, polish_control = list(max_iter = 0L)))
  for (nm in names(nl_designs)) {
    for (pm in c("newton", "bhhh")) {
      a <- nl_designs[[nm]]
      if (!is.null(a$control)) {
        out[[paste("polish|nl", nm, pm, sep = "|")]] <- p5a_rec(do.call(
          run_nestlogit, c(list(a[[1L]], "id", "alt", "choice", c("x1", "x2"),
                                nest_col = "nest"), a[-1L], list(polish = pm))))
      } else {
        rec_f <- if (nm == "custom") p8c_rec_nocall else p5a_rec
        out[[paste("polish|nl", nm, pm, sep = "|")]] <-
          rec_f(do.call(nl_fit, c(a, list(polish = pm))))
      }
    }
  }
  for (pm in c("newton", "bhhh")) {
    out[[paste("polish|nl|lambda_above_one", pm, sep = "|")]] <- p5a_rec(run_nestlogit(
      mc2, "id", "mode", "choice", c("gcost", "wait"), nest_col = "nest",
      control = list(maxeval = 60L), polish = pm,
      polish_control = list(max_iter = 50L)))
  }
  out[["polish|nl|cross|none"]] <- p5a_rec(nl_fit(dn, polish = "none"))
  fn <- suppressWarnings(suppressMessages(nl_fit(dn, se_method = "bhhh")))
  d <- fn$data; th <- unname(coef(fn))
  lim <- p8c_memory_limit(
    function() choicer:::nl_loglik_hessian_parallel(th, d$X, d$alt_idx, d$choice_idx,
                                                    d$nest_idx, d$M, d$weights, TRUE, FALSE),
    list(function() choicer:::nl_loglik_gradient_parallel(th, d$X, d$alt_idx,
                                                          d$choice_idx, d$nest_idx, d$M,
                                                          d$weights, TRUE, FALSE, TRUE),
         function() choicer:::nl_bhhh_parallel(th, d$X, d$alt_idx, d$choice_idx,
                                               d$nest_idx, d$M, d$weights, TRUE, FALSE)))
  out[["polish|nl|memory|newton"]] <- c(list(limit = lim), p8c_with_limit(
    lim, nl_fit(dn, se_method = "bhhh", polish = "newton")))
  out
}

# ---------------------------------------------------------------------------
# p8c_mxl
# ---------------------------------------------------------------------------

run_p8c_mxl <- function() {
  out <- list()
  mxl <- function(dt, ...) {
    run_mxlogit(data = dt, id_col = "id", alt_col = "alt", choice_col = "choice",
                covariate_cols = c("x1", "x2"), random_var_cols = c("w1", "w2"),
                S = 20L, ...)
  }
  rec <- function(key, expr) out[[paste("p8c_mxl", key, sep = "|")]] <<- p5a_rec(expr)
  # the outside option never chosen: every inside constant rises together
  sim <- simulate_mxl_data(N = 600L, J = 30L, beta = c(0.8, -0.6),
                           delta = seq(-1, 1, length.out = 30L),
                           Sigma = diag(c(1, 0.5)), outside_option = TRUE,
                           seed = 11L)
  dt <- data.table::as.data.table(sim$data)
  dt <- dt[id %in% dt[alt != 0L & choice == 1L, id] & alt != 0L]
  oo <- list(include_outside_option = TRUE, outside_opt_label = 0L)
  fit_o <- function(...) do.call(mxl, c(list(dt), oo, list(control = list(maxeval = 200L),
                                                          polish = "newton"), list(...)))
  rec("outside_never|cert0", fit_o(polish_control = list(max_iter = 0L, tol_decrement = 1e3)))
  rec("outside_never|steps", fit_o(polish_control = list(max_iter = 100L)))
  rec("outside_never|upper", fit_o(upper = c(ASC_15 = 2), polish_control = list(max_iter = 100L)))
  dt_held <- dt[!(alt == 30L & !id %in% dt[alt == 30L & choice == 1L, id])]
  rec("outside_never|held", do.call(mxl, c(list(dt_held), oo, list(
    control = list(maxeval = 6L), polish = "newton",
    polish_control = list(max_iter = 0L, tol_decrement = 1e3)))))
  dt_b <- dt[!id %in% dt[alt == 30L & choice == 1L, id]]
  rec("outside_never|pinned_never", do.call(mxl, c(list(dt_b), oo, list(
    control = list(maxeval = 200L), polish = "newton",
    lower = c(ASC_30 = -5), upper = c(ASC_30 = -5),
    polish_control = list(max_iter = 0L, tol_decrement = 1e3)))))
  # the base alternative never chosen, and chosen wherever it is offered
  sim <- simulate_mxl_data(N = 600L, J = 6L, beta = c(0.8, -0.6),
                           Sigma = diag(c(1, 0.5)), outside_option = FALSE, seed = 12L)
  dtb <- data.table::as.data.table(sim$data)
  dtb <- dtb[!id %in% dtb[alt == 1L & choice == 1L, id]]
  rec("base_never|cert0", mxl(dtb, control = list(maxeval = 100L), polish = "newton",
                              polish_control = list(max_iter = 0L, tol_decrement = 1e3)))
  sim <- simulate_mxl_data(N = 900L, J = 6L, beta = c(0.8, -0.6),
                           Sigma = diag(c(1, 0.5)), outside_option = FALSE, seed = 13L)
  dta <- data.table::as.data.table(sim$data)
  dta <- dta[!id %in% dta[alt == 1L, id] | id %in% dta[alt == 1L & choice == 1L, id]]
  rec("base_always|upper", mxl(dta, control = list(maxeval = 6L), upper = c(ASC_2 = 5),
                               polish = "newton",
                               polish_control = list(max_iter = 0L, tol_decrement = 1e3)))
  # every se_method with each polish, a custom optimizer, keep_data = FALSE
  sim <- simulate_mxl_data(N = 300L, J = 6L, beta = c(0.8, -0.6),
                           delta = seq(-0.5, 0.5, length.out = 6L),
                           Sigma = diag(c(1, 0.5)), outside_option = TRUE, seed = 71L)
  dx <- data.table::as.data.table(sim$data)[alt != 0L]
  dx[, clus := (id - 1L) %/% 10L]
  set.seed(74)
  wt <- stats::runif(300L, 0.5, 2)
  dx[, w := wt[id]]
  loose <- c(oo, list(control = list(maxeval = 6L)))
  for (se in c("hessian", "bhhh", "sandwich", "cluster")) {
    extra <- switch(se, cluster = list(cluster_col = "clus"),
                    sandwich = list(weights_col = "w"), list())
    for (pm in c("newton", "bhhh")) {
      rec(paste("se", se, pm, sep = "|"),
          do.call(mxl, c(list(dx), loose, extra, list(se_method = se, polish = pm))))
    }
  }
  custom <- function(theta_init, eval_f, lower, upper, control) {
    r <- nloptr::nloptr(theta_init, eval_f, lb = lower, ub = upper,
                        opts = list(algorithm = "NLOPT_LD_LBFGS", maxeval = 5L))
    list(par = r$solution, value = r$objective, convergence = r$status)
  }
  for (pm in c("newton", "bhhh")) {
    out[[paste("p8c_mxl|custom", pm, sep = "|")]] <- p8c_rec_nocall(
      do.call(mxl, c(list(dx), oo, list(optimizer = custom, polish = pm))))
    rec(paste("nodata", pm, sep = "|"),
        do.call(mxl, c(list(dx), loose, list(keep_data = FALSE, polish = pm))))
  }
  # a polish that fails on its curvature's memory check (the variance's BHHH
  # matrix fits)
  f0 <- suppressWarnings(suppressMessages(do.call(mxl, c(list(dx), loose,
                                                         list(se_method = "bhhh")))))
  d <- f0$data; th <- unname(coef(f0))
  gp <- choicer:::.mxl_gen_params(f0$draws_info, N = length(choicer:::.unit_first(d)))
  ka <- list(theta = th, X = d$X, W = d$W, alt_idx = d$alt_idx, choice_idx = d$choice_idx,
             M = d$M, weights = d$weights, eta_draws = gp$eta_draws, rc_dist = f0$rc_dist,
             rc_correlation = f0$rc_correlation, rc_mean = f0$rc_mean, use_asc = TRUE,
             include_outside_option = TRUE)
  lim <- p8c_memory_limit(
    function() do.call(choicer:::mxl_hessian_parallel, ka),
    list(function() do.call(choicer:::mxl_loglik_gradient_parallel, c(ka, list(opg_diag = TRUE))),
         function() do.call(choicer:::mxl_bhhh_parallel, ka)))
  out[["p8c_mxl|memory|newton"]] <- c(list(limit = lim), p8c_with_limit(
    lim, do.call(mxl, c(list(dx), loose, list(se_method = "bhhh", polish = "newton")))))
  out
}
