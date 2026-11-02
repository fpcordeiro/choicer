# ============================================================================
# Monte Carlo performance measures of a cell's replications, their Monte Carlo
# SEs (Morris, White & Crowther 2019, Table 6), and the gate rule (plan §2):
#   metrics_cell()  one (family, cell, method)'s reps -> a table of metrics;
#   gate_family()   Holm over a family's gated rows, then the tolerance rule;
#   pooled_tails()  coverage pooled over the constants within each rep;
#   rate_fit()      a ladder's log-log slope across cells.
# Definitions only. Records are lib/fit.R's rep records; a failed rep enters
# the failure rate only. Every metric needs two reps: below, it is NA and so
# is its MCSE.
# ============================================================================

# The pre-registered tolerances of the gate rule (plan §2; quick: plan Q3);
# `failure` is the tolerance of a variance route's or the joint Wald's
# failure rate.
tolerances <- function(set = c("mc", "quick")) {
  set <- match.arg(set)
  switch(set,
         mc    = c(bias_std = 0.10, se_ratio = 0.05, tail = 0.0075, coverage = 0.015, size = 0.015,
                   failure = 0.005),
         quick = c(bias_std = 0.30, se_ratio = 0.15, tail = 0.025,  coverage = 0.05,  size = 0.025,
                   failure = 0.02))
}

# The metrics: their kind (an estimand's; an estimand's by variance route; a
# cell's test; a cell's by variance route; an SBC quantity's; the cell's),
# whether a family may gate them, and the name of their tolerance in
# tolerances(). ese, lr_mean_q (E LR = q (1 + c/N), Bartlett 1937),
# failure_rate and global_max_rate are findings only; sbc_ecdf is gated on
# its p-value alone (plan §2).
VS_METRICS <- data.frame(
  metric   = c("bias_std", "ese", "se_ratio", "tail_lower", "tail_upper", "coverage",
               "joint_wald_size", "joint_wald_failure_rate", "lr_size", "lr_mean_q",
               "route_failure_rate", "sbc_ecdf", "failure_rate", "global_max_rate"),
  kind     = c("estimand", "estimand", "route", "route", "route", "route",
               "test", "test", "test", "test", "cellroute", "sbc", "cell", "cell"),
  gateable = c(TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE,
               FALSE),
  tol      = c("bias_std", NA, "se_ratio", "tail", "tail", "coverage", "size", "failure", "size",
               NA, "failure", NA, NA, NA),
  stringsAsFactors = FALSE)

VS_METRIC_COLS <- c("family", "cell", "method", "metric", "estimand", "route", "value", "mcse",
                    "target", "tol", "p", "gated", "verdict", "R", "n_failed", "block")

# The metrics that are a share of R independent reps (a binomial count over
# R): their p-value is the exact binomial test's (.gate_p()).
VS_BINOMIAL <- c("tail_lower", "tail_upper", "coverage", "joint_wald_size", "lr_size",
                 "joint_wald_failure_rate", "route_failure_rate")

# The metrics' targets at confidence level `level`.
.metric_targets <- function(level) {
  t <- (1 - level) / 2
  c(bias_std = 0, ese = NA, se_ratio = 1, tail_lower = t, tail_upper = t, coverage = level,
    joint_wald_size = 0.05, joint_wald_failure_rate = 0, lr_size = 0.05, lr_mean_q = 1,
    route_failure_rate = 0, sbc_ecdf = NA, failure_rate = NA, global_max_rate = 0)
}

# --- Proportions and the cell's tests ---------------------------------------

# MCSE of a proportion `value` = x / n over n reps, about the estimate (Morris
# Table 6, "Coverage" and "Power" rows: sqrt(p (1 - p) / n)), at the smoothed
# p~ = (x + 1) / (n + 2): never 0 for a count of 0 or n (with R = 100 a 2.5%
# tail has no exceedance in 8% of runs), Morris's value at large n. The gate
# does not read it: binomial metrics are tested exactly (.gate_p()).
# Vectorized over `value` and `n`; NA where n < 2.
.prop_mcse <- function(value, n) {
  p <- (round(value * n) + 1) / (n + 2)
  out <- sqrt(p * (1 - p) / pmax(n, 2))
  out[rep_len(n < 2, length(out))] <- NA_real_
  out
}

# Size of a chi-square test over reps: the rate of p < alpha (W ~ chi^2_df),
# with its MCSE (Morris Table 6, "Power" row, smoothed as in .prop_mcse()),
# over the reps with a finite statistic only (an infinite W has p = 0, but
# it is a failed statistic, counted in joint_wald_failure_rate); R is their
# count.
joint_wald_size <- function(W, df, alpha = 0.05) {
  ok <- is.finite(W) & is.finite(df)
  p <- stats::pchisq(W[ok], df[ok], lower.tail = FALSE)
  n <- length(p)
  value <- if (n >= 2) mean(p < alpha) else NA_real_
  c(value = value, mcse = .prop_mcse(value, n), R = n)
}

# A likelihood-ratio test over reps (LR ~ chi^2_q, Wilks): its size, as in
# joint_wald_size(), and mean(LR / q) (target 1) with the MCSE of a mean,
# sd / sqrt(R) (Morris Table 6, "Bias" row). A negative LR is kept.
lr_summary <- function(LR, q, alpha = 0.05) {
  ok <- is.finite(LR) & is.finite(q)
  r <- LR[ok] / q[ok]
  n <- length(r)
  s <- joint_wald_size(LR[ok], q[ok], alpha)
  c(size = s[["value"]], size_mcse = s[["mcse"]],
    mean_q = if (n >= 2) mean(r) else NA_real_,
    mean_q_mcse = if (n >= 2) stats::sd(r) / sqrt(n) else NA_real_, R = n)
}

# The reps' tests in a record field (joint_wald or lr: one test, or a named
# list of them), collected by test: list(name, stat, df) per test, name NA
# for an unnamed test.
.collect_tests <- function(recs, field, stat) {
  per <- Filter(Negate(is.null), lapply(recs, `[[`, field))
  single <- vapply(per, function(x) all(c(stat, "df", "p") %in% names(x)), logical(1))
  tests <- unlist(lapply(seq_along(per), function(i) if (single[i]) list(per[[i]]) else per[[i]]),
                  recursive = FALSE)
  if (!length(tests)) return(list())
  nm <- unlist(lapply(seq_along(per), function(i) {
    if (single[i]) NA_character_ else names(per[[i]]) %||% rep(NA_character_, length(per[[i]]))
  }))
  st <- unname(vapply(tests, function(v) as.numeric(v[[stat]])[1], 1))
  df <- unname(vapply(tests, function(v) as.numeric(v$df)[1], 1))
  lapply(unique(nm), function(k) {
    i <- nm %in% k
    list(name = k, stat = st[i], df = df[i])
  })
}

# --- Rows -------------------------------------------------------------------

# Rows of the metrics table from columns (recycled). `p` is set only where it
# does not follow from value and MCSE (sbc_ecdf: the exact ECDF test's).
.metric_rows <- function(metric, estimand = NA_character_, block = NA_character_,
                         route = NA_character_, value, mcse, R, src = "cell", egate = FALSE,
                         p = NA_real_) {
  data.frame(metric = metric, estimand = as.character(estimand), block = as.character(block),
             route = as.character(route), value = as.numeric(value), mcse = as.numeric(mcse),
             R = as.integer(R), src = src, egate = egate, p = as.numeric(p),
             stringsAsFactors = FALSE)
}

# One estimand's rows from its errors d = est - truth over reps and its SEs
# by route (a matrix, one column per route with an SE): standardized bias
# and ESE, then per route the SE ratio, the tails of the Wald interval at
# `level` and its coverage.
.estimand_rows <- function(name, block, d, se, level, src, egate) {
  t <- (1 - level) / 2
  z <- stats::qnorm(1 - t)
  f <- is.finite(d)
  n <- sum(f)
  ese <- if (n >= 2) stats::sd(d[f]) else NA_real_
  ok <- isTRUE(ese > 0)
  out <- list(.metric_rows(
    c("bias_std", "ese"), name, block,
    value = c(if (ok) mean(d[f]) / ese else NA, ese),
    # "Bias" row over ESE: sqrt(sum (est - mean)^2 / (n (n - 1))) / ESE = 1 / sqrt(n);
    # "Empirical SE" row: ESE / sqrt(2 (n - 1)).
    mcse = c(if (ok) 1 / sqrt(n) else NA, if (n >= 2) ese / sqrt(2 * (n - 1)) else NA),
    R = n, src = src, egate = egate))
  for (r in colnames(se)) {
    k <- f & is.finite(se[, r])
    dr <- d[k]
    sr <- se[k, r]
    nr <- length(dr)
    e <- if (nr >= 2) stats::sd(dr) else NA_real_
    v <- m <- rep(NA_real_, 4)
    if (isTRUE(e > 0)) {
      ms2 <- mean(sr^2)
      v[1] <- sqrt(ms2) / e
      # "Relative % error in ModSE" row, without the 100:
      # (ModSE / EmpSE) sqrt(Var(se^2) / (4 n ModSE^4) + 1 / (2 (n - 1))).
      m[1] <- v[1] * sqrt(stats::var(sr^2) / (4 * nr * ms2^2) + 1 / (2 * (nr - 1)))
    }
    if (nr >= 2) {
      v[2] <- mean(dr > z * sr)    # truth below est - z se
      v[3] <- mean(dr < -z * sr)   # truth above est + z se
      v[4] <- 1 - v[2] - v[3]
      m[2:4] <- .prop_mcse(v[2:4], nr)
    }
    out[[length(out) + 1L]] <- .metric_rows(
      c("se_ratio", "tail_lower", "tail_upper", "coverage"), name, block, r,
      value = v, mcse = m, R = nr, src = src, egate = egate)
  }
  do.call(rbind, out)
}

# The OK reps' estimands stacked: src ("est" or "derived"), name, block
# ("derived" for derived estimands), the error d = est - truth, the SEs by
# route (a derived estimand's `se` as the hessian route's) and est$gated.
.stack_estimands <- function(recs) {
  se_cols <- paste0("se_", VS_ROUTES)
  pick <- function(field, src) {
    x <- as.data.frame(data.table::rbindlist(lapply(recs, `[[`, field), use.names = TRUE,
                                             fill = TRUE))
    if (!nrow(x)) return(NULL)
    if (src == "derived") {
      x$block <- "derived"
      x$gated <- NA
      x$se_hessian <- x$se
    }
    out <- data.frame(src = src, name = x$name, block = x$block, d = x$est - x$truth,
                      stringsAsFactors = FALSE)
    for (k in se_cols) out[[k]] <- if (is.null(x[[k]])) NA_real_ else as.numeric(x[[k]])
    out$gated <- as.logical(x$gated)
    out
  }
  rbind(pick("est", "est"), pick("derived", "derived"))
}

# --- Gates ------------------------------------------------------------------

# The gate entries for `method`, parsed: metric[:route][@filter]. The part
# after ':' is a variance route for SE metrics and route_failure_rate and a
# test's name for the cell's tests; '@' filters estimand rows by block or
# estimand name, and sbc_ecdf rows by quantity. Unknown metrics, metrics
# that are never gated and misplaced suffixes are errors, so that a typo
# cannot silently ungate a check.
.parse_gates <- function(gates, method) {
  e <- if (is.null(gates)) NULL
       else if (is.list(gates$metrics)) gates$metrics[[method]]
       else gates$metrics
  e <- trimws(as.character(e))
  at <- regexpr("@", e, fixed = TRUE)
  left <- ifelse(at > 0, substr(e, 1L, at - 1L), e)
  filter <- ifelse(at > 0, substring(e, at + 1L), NA_character_)
  col <- regexpr(":", left, fixed = TRUE)
  metric <- ifelse(col > 0, substr(left, 1L, col - 1L), left)
  route <- ifelse(col > 0, substring(left, col + 1L), NA_character_)
  m <- match(metric, VS_METRICS$metric)
  kind <- VS_METRICS$kind[m]
  bad <- function(i, why) stop("gate entry '", e[i], "' (method '", method, "'): ", why,
                               call. = FALSE)
  for (i in seq_along(e)) {
    if (is.na(m[i])) {
      bad(i, paste0("unknown metric; gateable: ",
                    paste(VS_METRICS$metric[VS_METRICS$gateable], collapse = ", ")))
    }
    if (!VS_METRICS$gateable[m[i]]) bad(i, paste0("'", metric[i], "' is a finding, never gated"))
    if (any(!nzchar(c(route[i], filter[i])), na.rm = TRUE)) bad(i, "empty suffix")
    if (kind[i] %in% c("route", "cellroute") && !is.na(route[i]) && !route[i] %in% VS_ROUTES) {
      bad(i, paste0("unknown route; routes: ", paste(VS_ROUTES, collapse = ", ")))
    }
    if (kind[i] %in% c("estimand", "sbc") && !is.na(route[i])) bad(i, "this metric takes no route")
    if (kind[i] %in% c("test", "cellroute") && !is.na(filter[i])) {
      bad(i, "a cell's metric takes no '@' filter")
    }
  }
  data.frame(metric = metric, route = route, filter = filter, kind = kind,
             stringsAsFactors = FALSE)
}

# Which rows a parsed entry gates. An SE metric without a route matches
# every route; a test's suffix names the test. Without '@', derived
# estimands are not matched; with it, estimand rows whose block or name is
# the filter. Rows of the est table also need est$gated TRUE in every rep.
.entry_matches <- function(en, tab) {
  hit <- tab$metric == en$metric
  if (!is.na(en$route)) {
    hit <- hit & if (en$kind == "test") tab$estimand %in% en$route else tab$route %in% en$route
  }
  hit <- hit & if (is.na(en$filter)) tab$src != "derived"
               else tab$src != "cell" & (tab$block %in% en$filter | tab$estimand %in% en$filter)
  hit & (tab$src != "est" | tab$egate)
}

# --- The cell's table -------------------------------------------------------

# The metrics of one (family, cell, method)'s reps (plan §2), a data.frame
# with columns family, cell, method, metric, estimand, route, value, mcse,
# target, tol, p, gated, verdict (NA until gate_family()), R (the reps
# entering the row), n_failed, and block (the estimand's block, "derived"
# for derived estimands, NA for the cell's rows). Rows, in order:
#   per estimand (est names, then derived names, in order of appearance;
#   the truth is each rep's own): bias_std, ese; per route whose SEs are
#   not all NA: se_ratio, tail_lower, tail_upper, coverage;
#   per SBC quantity, when records carry diag$sbc = list(L, ranks =
#   c(<quantity> = rank, ...)): sbc_ecdf, value the ECDF test's statistic,
#   p its exact p-value (sbc_ecdf_test()), R the reps with a rank;
#   per route with an SE in some rep or named by a gate entry:
#   route_failure_rate (the share of OK reps whose est table lacks an SE of
#   that route; target 0);
#   per cell: joint_wald_size and joint_wald_failure_rate (W not finite;
#   target 0), lr_size and lr_mean_q when records carry the tests (a row per
#   test, `estimand` the test's name when a record holds a named list of
#   them), failure_rate, and global_max_rate when some rep ran the
#   multi-start check.
# `gates` is the family's (validate_family()): a row is gated when its cell
# is among gates$cell and an entry for its method matches it. An entry that
# matches no row of a gated cell adds a row with value NA, so that a gated
# check that cannot be computed is "na" in gate_family(), never "pass".
# `tol_set` is "mc", "quick" or a named numeric like tolerances()'s.
metrics_cell <- function(records, gates = NULL, tol_set = "mc", level = 0.95) {
  tol <- if (is.numeric(tol_set)) tol_set else tolerances(tol_set)
  need <- setdiff(names(tolerances("mc")), "failure")    # without it, failure gates are "na"
  if (!all(need %in% names(tol))) {
    stop("metrics_cell(): `tol_set` must name ", paste(need, collapse = ", "))
  }
  recs <- Filter(is.list, records)
  if (!length(recs)) {
    tab <- .metric_rows("failure_rate", value = NA, mcse = NA, R = 0L)[0, , drop = FALSE]
    return(.finish_rows(tab, rep(NA_character_, 3), tol, level, 0L))
  }
  id <- unique(t(vapply(recs, function(r) as.character(c(r$family, r$cell, r$method))[1:3],
                        character(3))))
  if (nrow(id) != 1L) {
    stop("metrics_cell(): records span several (family, cell, method): ",
         paste(apply(id, 1, paste, collapse = "/"), collapse = ", "))
  }
  ok <- vapply(recs, function(r) isTRUE(r$ok), logical(1))
  good <- recs[ok]
  rows <- list()

  # Estimands.
  long <- .stack_estimands(good)
  if (!is.null(long) && nrow(long)) {
    key <- factor(paste(long$src, long$name, sep = "\r"),
                  levels = unique(paste(long$src, long$name, sep = "\r")))
    for (i in split(seq_len(nrow(long)), key)) {
      se <- as.matrix(long[i, paste0("se_", VS_ROUTES), drop = FALSE])
      colnames(se) <- VS_ROUTES
      se <- se[, colSums(!is.na(se)) > 0, drop = FALSE]
      rows[[length(rows) + 1L]] <- .estimand_rows(long$name[i[1]], long$block[i[1]], long$d[i],
                                                  se, level, long$src[i[1]],
                                                  isTRUE(all(long$gated[i])))
    }
  }

  # SBC quantities: diag$sbc = list(L, ranks = c(<quantity> = rank, ...)).
  sbc <- Filter(Negate(is.null), lapply(good, function(r) r$diag$sbc))
  if (length(sbc)) {
    Ls <- unique(vapply(sbc, function(s) as.numeric(s$L %||% NA)[1], 1))
    if (length(Ls) != 1L || is.na(Ls)) stop("metrics_cell(): SBC records disagree on L or lack it")
    for (k in unique(unlist(lapply(sbc, function(s) names(s$ranks))))) {
      rk <- vapply(sbc, function(s) as.numeric(s$ranks[k])[1], 1)
      rk <- rk[!is.na(rk)]
      t <- if (length(rk) >= 2) sbc_ecdf_test(rk, Ls) else list(stat = NA, p = NA)
      rows[[length(rows) + 1L]] <- .metric_rows("sbc_ecdf", k, value = t$stat, mcse = NA,
                                                R = length(rk), src = "sbc", p = t$p)
    }
  }

  # Variance routes: per rep with an est table, whether the route's SEs are
  # missing for any estimand (a failed or singular route drops the rep from
  # that route's SE metrics). A row per route with an SE in some rep or named
  # by a gate entry.
  en <- .parse_gates(gates, id[1, 3])
  ests <- Filter(function(e) NROW(e) > 0, lapply(good, `[[`, "est"))
  if (length(ests)) {
    se_ok <- vapply(ests, function(e) {
      vapply(VS_ROUTES, function(r) {
        s <- is.finite(e[[paste0("se_", r)]])
        c(all(s), any(s))
      }, logical(2))
    }, matrix(NA, 2, length(VS_ROUTES)))                  # 2 x routes x reps
    named <- en$route[en$kind %in% c("route", "cellroute")]
    for (j in seq_along(VS_ROUTES)) {
      if (!any(se_ok[2, j, ]) && !VS_ROUTES[j] %in% named) next
      n <- length(ests)
      v <- if (n >= 2) mean(!se_ok[1, j, ]) else NA_real_
      rows[[length(rows) + 1L]] <- .metric_rows("route_failure_rate", route = VS_ROUTES[j],
                                                value = v, mcse = .prop_mcse(v, n), R = n)
    }
  }

  # The cell's tests, failures and multi-start checks.
  for (w in .collect_tests(good, "joint_wald", "W")) {
    s <- joint_wald_size(w$stat, w$df)
    n <- length(w$stat)
    v <- if (n >= 2) mean(!is.finite(w$stat)) else NA_real_
    rows[[length(rows) + 1L]] <- .metric_rows(c("joint_wald_size", "joint_wald_failure_rate"),
                                              w$name, value = c(s[["value"]], v),
                                              mcse = c(s[["mcse"]], .prop_mcse(v, n)),
                                              R = c(s[["R"]], n))
  }
  for (w in .collect_tests(good, "lr", "LR")) {
    s <- lr_summary(w$stat, w$df)
    rows[[length(rows) + 1L]] <- .metric_rows(c("lr_size", "lr_mean_q"), w$name,
                                              value = s[c("size", "mean_q")],
                                              mcse = s[c("size_mcse", "mean_q_mcse")], R = s[["R"]])
  }
  n_all <- length(recs)
  fr <- if (n_all >= 2) mean(!ok) else NA_real_
  rows[[length(rows) + 1L]] <- .metric_rows("failure_rate", value = fr,
                                            mcse = .prop_mcse(fr, n_all), R = n_all)
  # A missed global maximum: a refit from a perturbed start reached a
  # meaningfully higher likelihood (dlogL > 0.01) at a different point (max
  # |dtheta| / SE > 0.1); smaller gains are the optimizer's tolerance (with
  # dlogL > 1e-6 one optimum read as missed in 98% of reps).
  ran <- Filter(Negate(is.null), lapply(good, function(r) r$diag$multistart))
  if (length(ran)) {
    dl <- vapply(ran, function(m) as.numeric(m$dlogL %||% NA)[1], 1)
    dt <- vapply(ran, function(m) as.numeric(m$max_dtheta_se %||% NA)[1], 1)
    k <- !is.na(dl)
    missed <- dl[k] > 0.01 & !is.na(dt[k]) & dt[k] > 0.1
    gm <- if (sum(k) >= 2) mean(missed) else NA_real_
    rows[[length(rows) + 1L]] <- .metric_rows("global_max_rate", value = gm,
                                              mcse = .prop_mcse(gm, sum(k)), R = sum(k))
  }
  tab <- do.call(rbind, rows)

  # Gates: parsed in every cell (typos fail early), applied in gated cells.
  tab$gated <- FALSE
  if (!is.null(gates) && id[1, 2] %in% gates$cell) {
    tab <- .apply_gates(tab, en)
  }
  .finish_rows(tab, id[1, ], tol, level, sum(!ok))
}

# Marks the rows the parsed entries gate, and adds a gated row with value NA
# (src "missing") for each check an entry asks for that has no row: per
# estimand the entry covers (est$gated TRUE, or derived and named by '@'),
# its row for the entry's metric and route (any route when the entry names
# none), absent when the route's SEs are NA in every rep; and the entry
# itself when it covers no estimand and matches no row (a test the records
# do not carry, a block without estimands, every rep failed). A gated check
# is therefore never silently absent: gate_family() makes it "na". Gating
# an SE metric of a route also gates that route's route_failure_rate, and
# gating a joint Wald test's size its joint_wald_failure_rate, so that reps
# dropped by a failed variance or Wald statistic count against the gate.
.apply_gates <- function(tab, en) {
  n0 <- nrow(tab)
  est <- unique(tab[tab$src %in% c("est", "derived"), c("src", "estimand", "block", "egate")])
  miss <- list()
  add <- function(metric, estimand, block, route) {
    m <- .metric_rows(metric, estimand, block, route, value = NA, mcse = NA, R = 0L,
                      src = "missing")
    m$gated <- TRUE
    miss[[paste(metric, estimand, route, sep = "\r")]] <<- m
  }
  for (j in seq_len(nrow(en))) {
    e <- en[j, ]
    hit <- .entry_matches(e, tab[seq_len(n0), , drop = FALSE])
    tab$gated[seq_len(n0)] <- tab$gated[seq_len(n0)] | hit
    covered <- NULL
    if (e$kind %in% c("estimand", "route") && nrow(est)) {
      probe <- est
      probe$metric <- e$metric
      probe$route <- e$route
      covered <- est[.entry_matches(e, probe), , drop = FALSE]
      for (k in seq_len(nrow(covered))) {
        has <- hit & tab$src[seq_len(n0)] == covered$src[k] &
          tab$estimand[seq_len(n0)] %in% covered$estimand[k]
        if (!any(has)) add(e$metric, covered$estimand[k], covered$block[k], e$route)
      }
    }
    if (!any(hit) && !NROW(covered)) {
      add(e$metric, if (e$kind == "test") e$route else e$filter, NA,
          if (e$kind %in% c("route", "cellroute")) e$route else NA)
    }
  }
  i0 <- seq_len(n0)
  g <- tab$gated[i0]
  se_metric <- tab$metric[i0] %in% VS_METRICS$metric[VS_METRICS$kind == "route"]
  routes <- unique(c(tab$route[i0][g & se_metric], en$route[en$kind == "route"]))
  for (r in routes[!is.na(routes)]) {
    i <- which(tab$metric[i0] == "route_failure_rate" & tab$route[i0] %in% r)
    if (length(i)) tab$gated[i] <- TRUE else add("route_failure_rate", NA, NA, r)
  }
  tests <- tab$estimand[i0][g & tab$metric[i0] == "joint_wald_size"]
  i <- which(tab$metric[i0] == "joint_wald_failure_rate" & tab$estimand[i0] %in% tests)
  tab$gated[i] <- TRUE
  if (length(miss)) tab <- rbind(tab, do.call(rbind, unname(miss)))
  tab
}

# Adds the identity, targets, tolerances, p and the bookkeeping columns, in
# the table's column order.
.finish_rows <- function(tab, id, tol, level, n_failed) {
  if (is.null(tab$gated)) tab$gated <- rep(FALSE, nrow(tab))
  tm <- match(tab$metric, VS_METRICS$metric)
  tab$target <- unname(.metric_targets(level)[tab$metric])
  tab$tol <- unname(tol[VS_METRICS$tol[tm]])
  tab$p <- .gate_p(tab)
  tab$family <- rep(id[1], nrow(tab))
  tab$cell <- rep(id[2], nrow(tab))
  tab$method <- rep(id[3], nrow(tab))
  tab$verdict <- rep(NA_character_, nrow(tab))
  tab$n_failed <- rep(as.integer(n_failed), nrow(tab))
  rownames(tab) <- NULL
  tab[VS_METRIC_COLS]
}

# Two-sided p of each row's deviation from its target: 2 pnorm(-|value -
# target| / mcse), a zero MCSE giving 1 without a deviation and 0 with one,
# with two small-R references that are exact where the normal one is
# anti-conservative (at the mc tier's R they agree with it):
#   bias_std: sqrt(R) bias_std = mean / (sd / sqrt(R)) is the one-sample t
#   statistic, t with R - 1 df under normal errors (normal reference: 34% of
#   null families fail at R = 2, 10% at R = 4);
#   binomial metrics (VS_BINOMIAL): the exact test of the count x = value R
#   against Binomial(R, target), 2 min(P(X <= x), P(X >= x)) capped at 1
#   (normal reference: 5-10% of tail rows fail at R = 2-4, ~1e-3 per tail row
#   at R = 100 among 30 gated rows; exact: below 1e-4).
# rule = "band" gives gate_family()'s one-sided p against the nearer edge of
# the tolerance band instead (see there). Under either rule an sbc_ecdf row
# keeps its own p, the exact ECDF test's.
.gate_p <- function(tab, rule = "plan") {
  dev <- tab$value - tab$target
  if (rule == "band") {
    p <- .band_p(tab, dev)
  } else {
    zz <- ifelse(tab$mcse > 0, abs(dev) / tab$mcse, ifelse(dev == 0, 0, Inf))
    p <- 2 * stats::pnorm(-zz)
    s <- which(tab$metric == "bias_std" & !is.na(zz) & tab$R >= 2)
    p[s] <- 2 * stats::pt(-zz[s], tab$R[s] - 1)
    b <- which(tab$metric %in% VS_BINOMIAL & !is.na(tab$value) & !is.na(tab$target) &
                 tab$R >= 1)
    x <- round(tab$value[b] * tab$R[b])
    n <- tab$R[b]
    t <- tab$target[b]
    p[b] <- pmin(1, 2 * pmin(stats::pbinom(x, n, t),
                             stats::pbinom(x - 1, n, t, lower.tail = FALSE)))
  }
  s <- which(tab$metric == "sbc_ecdf")
  p[s] <- tab$p[s] %||% NA_real_
  p
}

# One-sided p that the true value lies inside the band [target - tol,
# target + tol], against the band's edge nearer the estimate: 1 inside the
# band (rounding-safe); else pnorm(-(|dev| - tol) / mcse); bias_std
# pt(-(|bias_std| - tol) sqrt(R), R - 1) (a central-t approximation); and for
# the binomial metrics the exact tail at the edge in the direction of the
# deviation, P(X >= x | target + tol) above, P(X <= x | target - tol) below.
.band_p <- function(tab, dev) {
  out <- !is.na(dev) & !is.na(tab$tol) & abs(dev) > tab$tol * (1 + 1e-9)
  p <- ifelse(is.na(dev) | is.na(tab$tol), NA_real_, 1)
  zz <- ifelse(tab$mcse > 0, (abs(dev) - tab$tol) / tab$mcse, Inf)
  i <- which(out & tab$metric != "bias_std" & !tab$metric %in% VS_BINOMIAL)
  p[i] <- stats::pnorm(-zz[i])
  i <- which(out & tab$metric == "bias_std" & tab$R >= 2)
  p[i] <- stats::pt(-zz[i], tab$R[i] - 1)
  i <- which(out & tab$metric %in% VS_BINOMIAL & tab$R >= 1)
  x <- round(tab$value[i] * tab$R[i])
  up <- dev[i] > 0
  edge <- tab$target[i] + ifelse(up, tab$tol[i], -tab$tol[i])
  p[i] <- ifelse(up, stats::pbinom(x - 1, tab$R[i], pmin(edge, 1), lower.tail = FALSE),
                 stats::pbinom(x, tab$R[i], pmax(edge, 0)))
  p
}

# The gate over metrics_cell() rows: p by the rule, p_holm the Holm
# adjustment over each family's gated rows with a p-value, and verdict
# "fail" or "pass" where gated, "na" where gated but the value, its MCSE,
# the target or the tolerance is NA, "finding" where not gated.
#   rule = "plan" (plan §2, the default): p tests value = target (two-sided,
#   .gate_p()); "fail" where p_holm < alpha and |value - target| > tol. A
#   true value inside the band but off target (a tail of 0.0278 against
#   0.025 +- 0.0075 at R = 1000) fails at a rate that grows with R.
#   rule = "band": p tests that the true value lies inside [target - tol,
#   target + tol] (one-sided, against the nearer edge, .band_p()); "fail"
#   where p_holm < alpha, so a true value anywhere inside the band fails a
#   family with probability at most alpha.
# Under both rules an sbc_ecdf row fails where p_holm < alpha, its exact p
# Holm-adjusted with the family's other gated rows (plan §2: "SBC: Holm p <
# 0.01 on the ECDF test"); it is "na" only without a p. The tolerance is
# compared as tol (1 + 1e-9), so that rounding cannot decide a verdict (the
# 95% tail target (1 - 0.95) / 2 is 0.025000000000000022).
gate_family <- function(tab, alpha = 0.01, rule = c("plan", "band")) {
  rule <- match.arg(rule)
  dev <- tab$value - tab$target
  tab$p <- .gate_p(tab, rule)
  tab$p_holm <- rep(NA_real_, nrow(tab))
  for (f in unique(tab$family)) {
    i <- which(tab$family %in% f & tab$gated & !is.na(tab$p))
    tab$p_holm[i] <- stats::p.adjust(tab$p[i], "holm")
  }
  sbc <- tab$metric == "sbc_ecdf"
  na <- is.na(tab$p_holm) | (!sbc & (is.na(tab$value) | is.na(tab$mcse) | is.na(tab$target) |
                                       is.na(tab$tol)))
  fail <- tab$p_holm < alpha & (sbc | rule == "band" | abs(dev) > tab$tol * (1 + 1e-9))
  tab$verdict <- ifelse(!tab$gated, "finding", ifelse(na, "na", ifelse(fail, "fail", "pass")))
  tab
}

# --- Pooled constants -------------------------------------------------------

# Coverage pooled over the estimands of `block` within each rep (plan §2: for
# thousands of constants): per rep, the shares of constants whose truth lies
# below and above the Wald interval of `route`, averaged over reps. A rep's
# constants are correlated (they share its data), so the MCSE is the
# between-rep SD of the per-rep shares over sqrt(R). `groups` maps constant
# names to groups (e.g. expected-count terciles) for one set of rows per
# group (constants it does not name are left out); without it the estimand
# column is "all". Rows bind with metrics_cell()'s: pooled_tail_lower,
# pooled_tail_upper, pooled_coverage, never gated, R the reps entering.
pooled_tails <- function(records, block = "asc", groups = NULL, level = 0.95, route = "hessian") {
  if (!route %in% VS_ROUTES) stop("pooled_tails(): unknown route '", route, "'")
  all_recs <- Filter(is.list, records)
  recs <- Filter(function(r) isTRUE(r$ok) && !is.null(r$est), all_recs)
  lev <- if (is.null(groups)) "all" else if (is.factor(groups)) levels(groups)
         else unique(as.character(groups))
  if (!is.null(groups)) groups <- stats::setNames(as.character(groups), names(groups))
  t <- (1 - level) / 2
  z <- stats::qnorm(1 - t)
  # Per rep and group: the shares below and above the interval.
  shares <- lapply(recs, function(r) {
    e <- r$est[r$est$block == block, , drop = FALSE]
    d <- e$est - e$truth
    s <- e[[paste0("se_", route)]]
    grp <- if (is.null(groups)) rep("all", nrow(e)) else unname(groups[e$name])
    k <- is.finite(d) & is.finite(s) & !is.na(grp)
    lo <- tapply(d[k] > z * s[k], factor(grp[k], levels = lev), mean)
    up <- tapply(d[k] < -z * s[k], factor(grp[k], levels = lev), mean)
    cbind(lo = as.numeric(lo), up = as.numeric(up))
  })
  tab <- do.call(rbind, lapply(seq_along(lev), function(g) {
    sh <- do.call(rbind, lapply(shares, function(x) x[g, ]))
    sh <- if (is.null(sh)) matrix(numeric(), 0, 2) else sh[!is.na(sh[, 1]), , drop = FALSE]
    n <- nrow(sh)
    cv <- 1 - sh[, 1] - sh[, 2]
    v <- if (n >= 2) c(mean(sh[, 1]), mean(sh[, 2]), mean(cv)) else rep(NA_real_, 3)
    m <- if (n >= 2) c(stats::sd(sh[, 1]), stats::sd(sh[, 2]), stats::sd(cv)) / sqrt(n)
         else rep(NA_real_, 3)
    .metric_rows(c("pooled_tail_lower", "pooled_tail_upper", "pooled_coverage"), lev[g], block,
                 route, value = v, mcse = m, R = n, src = "pooled")
  }))
  tab$target <- rep(c((1 - level) / 2, (1 - level) / 2, level), length(lev))
  r1 <- if (length(all_recs)) all_recs[[1]] else list()
  tab$family <- r1$family %||% NA_character_
  tab$cell <- r1$cell %||% NA_character_
  tab$method <- r1$method %||% NA_character_
  tab$tol <- NA_real_
  tab$p <- .gate_p(tab)
  tab$gated <- FALSE
  tab$verdict <- NA_character_
  tab$n_failed <- as.integer(length(all_recs) - sum(vapply(all_recs, function(r) isTRUE(r$ok),
                                                          logical(1))))
  rownames(tab) <- NULL
  tab[VS_METRIC_COLS]
}

# --- Ladders ----------------------------------------------------------------

# A ladder's rate: weighted least squares of log|value| on log x across the
# cells of `tab` (one method's rows: one row per cell after the filters),
# weights 1 / (mcse / |value|)^2, the delta-method variance of log|value|.
# The slope's SE is the known-variance one, sqrt([(X'WX)^-1]_22), not
# rescaled by the residuals (three cells leave one degree of freedom). `x`
# is a named numeric cell -> x (N, S, N/J ...); estimand NULL takes the
# metric's only estimand (rows without one for the cell's metrics); route
# NA takes rows without a route. Returns c(slope, se, intercept, n).
rate_fit <- function(tab, x, metric, estimand = NULL, route = NA) {
  if (!is.numeric(x) || is.null(names(x))) stop("rate_fit(): `x` must be a numeric named by cell")
  sel <- tab$metric == metric & tab$cell %in% names(x)
  sel <- sel & if (is.na(route)) is.na(tab$route) else tab$route %in% route
  if (!is.null(estimand)) {
    sel <- sel & tab$estimand %in% estimand
  } else if (length(unique(tab$estimand[sel])) > 1L) {
    stop("rate_fit(): metric '", metric, "' has several estimands; name one")
  }
  r <- tab[sel, , drop = FALSE]
  if (anyDuplicated(r$cell)) {
    stop("rate_fit(): several rows per cell; pass one family's and one method's rows")
  }
  r <- r[is.finite(r$value) & r$value != 0 & is.finite(r$mcse) & r$mcse > 0, , drop = FALSE]
  n <- nrow(r)
  na <- c(slope = NA_real_, se = NA_real_, intercept = NA_real_, n = n)
  lx <- log(x[r$cell])
  if (n < 2 || !all(is.finite(lx)) || length(unique(lx)) < 2) return(na)
  X <- cbind(1, lx)
  w <- (abs(r$value) / r$mcse)^2
  A <- crossprod(X, w * X)
  Ai <- solve(A)
  b <- Ai %*% crossprod(X, w * log(abs(r$value)))
  c(slope = b[2], se = sqrt(Ai[2, 2]), intercept = b[1], n = n)
}
