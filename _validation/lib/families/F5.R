# ============================================================================
# F5: simulation error of maximum simulated likelihood (plan §5; claim C4),
# paired. Simulated-likelihood bias is about b / S (E log P_S = log P -
# Var(P_S) / (2 P^2) + O(S^-2)), so bias / SD = O(sqrt(N) / S), and choicer's
# permuted Halton draws carry no randomized-QMC guarantee (MXL article §2.5,
# §5.4): the simulation error is measured here, on the gated regimes' own
# data, to tell whether the gates' condition sqrt(N) / S <= 0.25 is the
# right one. Findings only.
# Cells are regimes: another family's cell, whose simulate() gives the data
# (load_family() by name; dataset d is that family's rep d, on its seed):
#   F3         F3's top rung, 16k students (D4-40, cross-section)
#   F4         F4's top rung, 3.2k persons x 8 tasks (D1, the panel
#              likelihood)
#   F6         F6's middle rung, 20k patients (D3: 200 hospitals with
#              constants, no outside option, cross-section; source S = 200)
#   F7.claims  F7's claims, 625k patients, about 10^6 admissions (D3-claims,
#              2,000 hospitals, the panel likelihood; source S = 100)
#   F7.census  F7's census, 10^6 students (D4-census, 3,000 schools in 150
#              markets, an outside option; source S = 100)
# 30 datasets each, F7's regimes 1 (plan §5), whose reference is S = 800
# and ladder ends at 400. Store mode at S = 200 builds the K_w x S x U cube
# of draws: 3.2 GB for the census (10^6 students), 2.0 GB for the claims
# (625k patients), beside the 10-15 GB an F7 fit holds.
# One method, `msl`, per dataset (the ladder's fits depend on the
# reference's, so they run in one unit):
#   1. the source family's own fit of the rep, at the source cell's S and
#      draws seed (the rep's seed + 1, as F3, F4, F6 and F7 seed theirs):
#      fit_two_starts() on F3, F4 and F6 (choicer's default start and theta0,
#      the better certified fit kept: the S1 rule for estimation from
#      scratch), F7's f7_estimate() on its regimes (Newton steps from theta0,
#      S2, with the default start on the source's reps 10, 20, ...), so that
#      it stays the source's own estimate;
#   2. the reference, theta_ref, at S_ref generate-mode draws, warm-started
#      from it (step 1 has chosen the maximum's basin; at S_ref only the
#      simulation error moves it), which costs well under the two cold starts
#      at S_ref the S1 rule would take;
#   3. warm-started from theta_ref, fits at every S of S_gen in generate
#      mode with `seeds` draw seeds, and in store mode (plain Halton, no
#      seed) at every S of S_store. They measure simulation error, not
#      optimization: one start each, next to the maximum.
# The warm fits (2, 3) take Newton steps alone (fit_newton(), lib/fit.R: the
# optimizer stopped after one evaluation, the polish up to 30 steps), falling
# back to the optimizer where they are not certified
# (diag$fits$newton_only): from a start within
# the simulation error of the maximum the polish certifies it in 2-3 steps,
# the same maximum as the optimizer's (|difference| / SE 2e-5 on F7's census
# at N / 25, 3e-9 on its claims at N / 100) at 1/3 to 1/13 of its time:
# L-BFGS took 75 and 847 iterations there.
# Every fit is certified by the Newton polish (require_converged()). Draw
# seeds derive from the source rep's seed s: s + 1 (step 1, the source
# family's own), then streams k = 1 (the reference) and 1 + j (seed j of
# the ladder, shared across S) at s + k * 1000003 (.f5_subseed(), the
# stride of dgp_bayes.R's db_subseed()): the source's seeds are consecutive
# over reps, so a small offset would give one rep's ladder the next rep's
# reference seed, and the datasets' simulation errors would not be
# independent.
# The est table: the reference's rows "<parameter>/ref" (truth theta0) and
# the ladder's "<parameter>/S<S>/g<j>", "<parameter>/S<S>/store" and, for
# step 1's fit, "<parameter>/S<S>/src" (a generate-mode point at the
# source's own S; truth theta_ref). No row carries an SE: the reference's
# Hessian SE, the sampling SE that standardizes the findings, is in
# diag$se_ref (an SE beside a ladder row would make the generic SE ratio
# and route-failure rows meaningless). So the generic metrics are bias_std
# and ese only: for "/ref", the reference fit's against theta0 over the
# datasets; for a ladder row, bias_std = mean(theta_S - theta_ref) /
# sd(theta_S - theta_ref) over datasets and ese = sd(theta_S - theta_ref),
# the simulation error's spread across datasets in the parameter's units.
# The plan's quantities are f5_findings(records): simulation bias and SD in
# units of the reference's sampling SE against S, the bias's rate in S and
# the smallest S with |bias| < 0.1 SE, against S = 4 sqrt(U).
# Sourced by load_family() into a child of the library's environment.
# ============================================================================

# The regimes: the source family and cell, the cell's N (the source
# cell's: students, persons or patients; checks/wp7_sim.R --only=load
# compares), the reference's S, the generate ladder, the draw seeds per S
# and the store ladder (comma-separated S).
F5_REGIMES <- data.frame(
  cell = c("F3", "F4", "F6", "F7.claims", "F7.census"),
  src = c("F3", "F4", "F6", "F7", "F7"),
  src_cell = c("16k", "3.2k", "20k", "claims", "census"),
  N = c(16000L, 3200L, 20000L, 625000L, 1000000L),
  S_ref = c(3200L, 3200L, 3200L, 800L, 800L),
  S_gen = c(rep("50,100,200,400,800", 3L), rep("50,100,200,400", 2L)),
  seeds = 2L,
  S_store = "100,200",
  R = c(30L, 30L, 30L, 1L, 1L),
  stringsAsFactors = FALSE)

# Each regime's mixed logit, by regime: the source family's own fit call
# (F3's f3_run(), F4's f4_run(), F6's f6_run(), F7's f7_run() by design),
# its model arguments copied and checked against the family's fit by
# checks/wp7_sim.R --only=spec; polished; `...` carries S, draws, seed and
# theta_init.
F5_MODELS <- list(
  F3 = function(sim, ...) {
    choicer::run_mxlogit(sim$data, "id", "alt", "choice", D4_X, D4_W, rc_dist = c(1L, 0L),
                         rc_mean = FALSE, include_outside_option = TRUE, polish = "newton", ...)
  },
  F4 = function(sim, ...) {
    choicer::run_mxlogit(sim$data, "id", "alt", "choice", SP_X, SP_W, rc_dist = c(1L, 0L),
                         rc_mean = TRUE, rc_correlation = TRUE, include_outside_option = TRUE,
                         person_col = "person", polish = "newton", ...)
  },
  F6 = function(sim, ...) {
    choicer::run_mxlogit(sim$data, "id", "alt", "choice", D3_X$default, D3_W,
                         rc_dist = c(1L, 0L), rc_mean = TRUE, polish = "newton", ...)
  },
  F7.claims = function(sim, ...) {
    choicer::run_mxlogit(sim$data, "id", "alt", "choice", D3_X$default, D3_W,
                         rc_dist = c(1L, 0L), rc_mean = TRUE, person_col = "patient",
                         polish = "newton", ...)
  },
  F7.census = function(sim, ...) {
    choicer::run_mxlogit(sim$data, "id", "alt", "choice", D4_X, D4_W, rc_dist = c(1L, 0L),
                         rc_mean = FALSE, include_outside_option = TRUE, polish = "newton", ...)
  })

# The regimes that run: a model and a family file.
.f5_lib <- parent.env(environment())
.f5_cells <- F5_REGIMES[F5_REGIMES$cell %in% names(F5_MODELS) &
                          file.exists(file.path(.f5_lib$.vs_lib, "families",
                                                paste0(F5_REGIMES$src, ".R"))), ]
rownames(.f5_cells) <- NULL

.f5_loaded <- new.env(parent = emptyenv())
.f5_source <- function(id) {
  if (!exists(id, envir = .f5_loaded, inherits = FALSE)) {
    assign(id, load_family(id, .f5_lib), envir = .f5_loaded)
  }
  get(id, envir = .f5_loaded, inherits = FALSE)
}

# Stream k's seed from s, k * 1000003 steps on, in 1..2^31 - 2 as vs_seed()
# gives them: two reps' streams meet only 10^6 reps apart.
.f5_subseed <- function(s, k) as.integer(1 + (as.numeric(s) + k * 1000003) %% (2^31 - 2))

.f5_S <- function(x) if (is.na(x) || !nzchar(x)) integer() else as.integer(strsplit(x, ",")[[1]])

# Dataset d of a regime: the source family's rep d, simulated as its runner
# does (vs_set_seed(), then simulate() on the source cell's row), with
# meta$f5 naming it.
f5_simulate <- function(cell, seed) {
  src <- .f5_source(cell$src)
  d <- vs_seed_rep("F5", cell$cell, seed)
  s <- vs_seed(cell$src, cell$src_cell, d)
  row <- family_cell(src, cell$src_cell)
  vs_set_seed(s)
  sim <- src$simulate(row, s)
  sim$meta$f5 <- list(src = cell$src, src_cell = cell$src_cell, rep = d, seed = s,
                      S_src = row$S)
  sim
}

f5_msl <- function(sim, cell, seed) {
  model <- F5_MODELS[[cell$cell]]
  s0 <- sim$meta$f5$seed
  fits <- list()
  # newton_steps is the kept fit's polish, after the optimizer where
  # newton_only is FALSE; newton_error and newton_failed_sec are fit_newton()'s
  # record of a failed attempt, which time_sec leaves out.
  log_fit <- function(f, variant, S, draws, dseed) {
    fits[[length(fits) + 1L]] <<- data.frame(
      variant = variant, S = S, draws = draws, seed = dseed, code = f$convergence,
      iterations = f$optimizer$iterations %||% NA_integer_,
      newton_steps = f$optimizer$polish$iterations %||% NA_integer_,
      newton_only = isTRUE(attr(f, "newton_only")),
      newton_error = attr(f, "newton_error") %||% NA_character_,
      newton_failed_sec = attr(f, "newton_failed_sec") %||% NA_real_,
      decrement = f$optimizer$polish$decrement %||% NA_real_,
      time_sec = f$optimizer$elapsed_time %||% NA_real_, loglik = f$loglik,
      stringsAsFactors = FALSE)
  }
  # A warm fit (step 2 and 3).
  warm <- function(S, draws, dseed, start, what) {
    fit_newton(function(...) model(sim, S = S, draws = draws, seed = dseed, ...), start, what)
  }
  rows <- function(est, suffix) {
    est$name <- paste0(est$name, "/", suffix)
    est
  }
  S_src <- sim$meta$f5$S_src
  dsrc <- s0 + 1L
  src_fit <- function(...) model(sim, S = S_src, draws = "generate", seed = dsrc, ...)
  src_what <- sprintf("the fit at the source's S = %d", S_src)
  start <- if (cell$src == "F7") {
    environment(.f5_source("F7")$simulate)$f7_estimate(src_fit, sim$truth, sim$meta$f5$rep,
                                                       src_what)
  } else {
    fit_two_starts(function(theta_init) src_fit(theta_init = theta_init), sim$truth, src_what)
  }
  log_fit(start, "src", S_src, "generate", dsrc)
  dref <- .f5_subseed(s0, 1L)
  ref <- warm(cell$S_ref, "generate", dref, coef(start), "the reference fit")
  log_fit(ref, "ref", cell$S_ref, "generate", dref)
  th_ref <- coef(ref)
  V <- fit_vcovs(ref, "hessian")
  if (is.null(V$hessian)) stop("the reference fit has no Hessian variance", call. = FALSE)
  se_ref <- sqrt(diag(V$hessian))
  est <- list(rows(est_table(ref, sim$truth, V = list(), gated = FALSE), "ref"),
              rows(est_table(start, th_ref, V = list(), gated = FALSE),
                   sprintf("S%d/src", S_src)))
  ladder <- function(S, draws, j) {
    gen <- draws == "generate"
    tag <- if (gen) sprintf("g%d", j) else "store"
    dseed <- if (gen) .f5_subseed(s0, 1L + j) else NA_integer_
    f <- warm(S, draws, if (gen) dseed, th_ref, sprintf("the %s fit at S = %d (%s)", draws, S, tag))
    log_fit(f, tag, S, draws, dseed)
    rows(est_table(f, th_ref, V = list(), gated = FALSE), sprintf("S%d/%s", S, tag))
  }
  for (S in .f5_S(cell$S_gen)) for (j in seq_len(cell$seeds)) {
    est[[length(est) + 1L]] <- ladder(S, "generate", j)
  }
  for (S in .f5_S(cell$S_store)) est[[length(est) + 1L]] <- ladder(S, "store", 0L)
  list(est = do.call(rbind, est),
       vcov_c = vcov_common(ref, V = V$hessian),
       diag = c(fit_diag(ref), list(
         start = attr(start, "start"), S_ref = cell$S_ref, se_ref = se_ref,
         units = ref$n_persons %||% ref$nobs, src = unlist(sim$meta$f5),
         fits = do.call(rbind, fits))))
}

f5_estimands <- function(fit, sim) NULL

# --- Findings -------------------------------------------------------------------

# A record's ladder rows in the reference's SE units: z = (theta_S -
# theta_ref) / SE_ref per parameter, S, kind (generate, store, or source:
# step 1's fit, kept apart so that a source S on the ladder is not pooled
# with the ladder's seeds) and seed.
.f5_long <- function(r) {
  e <- r$est
  m <- regmatches(e$name, regexec("^(.*)/(ref|S([0-9]+)/(g([0-9]+)|store|src))$", e$name))
  part <- function(k) vapply(m, function(x) if (length(x)) x[k] else NA_character_, "")
  param <- part(2L)
  ref <- part(3L) %in% "ref"
  tag <- part(5L)
  out <- data.frame(regime = r$cell, rep = r$rep, param = param, block = e$block,
                    kind = ifelse(tag %in% "store", "store",
                                  ifelse(tag %in% "src", "source", "generate")),
                    S = as.integer(part(4L)), seed = part(5L),
                    z = (e$est - e$truth) / unname(r$diag$se_ref[param]),
                    stringsAsFactors = FALSE)
  out[!ref & !is.na(param), , drop = FALSE]
}

# Bias and SD of z at one (regime, parameter, kind, S): bias = the mean over
# datasets of each dataset's mean over seeds, its MCSE their SD / sqrt(D)
# (with one dataset, the seeds' SD / sqrt(seeds)); SD = the pooled
# within-dataset SD over seeds (nu = sum of seeds - 1 degrees of freedom),
# its MCSE SD / sqrt(2 nu) (normal theory). Store mode has one fit per
# dataset: no SD.
.f5_cell <- function(z, rep) {
  zd <- tapply(z, rep, mean)
  nd <- tapply(z, rep, length)
  D <- length(zd)
  ss <- sum(tapply(z, rep, function(v) sum((v - mean(v))^2)))
  nu <- sum(nd - 1L)
  sd <- if (nu > 0) sqrt(ss / nu) else NA_real_
  c(datasets = D, bias_se = mean(zd),
    bias_se_mcse = if (D >= 2) stats::sd(zd) / sqrt(D) else if (nu > 0) sd / sqrt(length(z))
                   else NA_real_,
    sd_se = sd, sd_se_mcse = if (nu > 0) sd / sqrt(2 * nu) else NA_real_)
}

# One ladder's z by independent unit (rows) and S (columns): a dataset's mean
# over seeds; with one dataset, one row per seed. A unit keeps its S points
# together, because they share a dataset, a reference fit and draw seeds, so a
# rate's uncertainty comes from the spread of its value over units and not from
# the S-level MCSEs, which ignore the covariance between S points. NA where a
# unit has no fit at an S.
.f5_units <- function(long, d) {
  z <- long[long$regime == d$regime[1] & long$param == d$param[1] & long$kind == d$kind[1] &
              long$S %in% d$S, , drop = FALSE]
  by <- if (length(unique(z$rep)) >= 2L) z$rep else z$seed
  tapply(z$z, list(by, z$S), mean)[, as.character(d$S), drop = FALSE]
}

# The MCSE of a statistic linear in the S-level means, sum_S w_S y_S: its SD over
# the units complete at the S points it uses / sqrt(units), so that every
# covariance among those points counts (the delta method for a statistic that is
# not linear; weights treated as fixed). NA below two units.
.f5_rate_se <- function(u, w) {
  on <- w != 0
  v <- drop(u[, on, drop = FALSE] %*% w[on])
  v <- v[!is.na(v)]
  if (length(v) >= 2L) stats::sd(v) / sqrt(length(v)) else NA_real_
}

# F5's findings from its records (every regime): a list of
#   by_S    per regime, parameter (the common ones; the constants pooled as
#           "(constants, RMS)": the root mean square over constants), kind
#           (generate, store, source) and S: bias / SE and SD / SE with
#           MCSEs (SE the reference fit's Hessian SE, per dataset);
#   rates   per regime, common parameter and ladder (generate, store): b, the bias's
#           coefficient in bias = b (1/S - 1/S_ref) (the reference's own
#           b / S_ref subtracted; weighted least squares on the MCSEs), so
#           that b / S is the bias against infinite S in SE units; S_01 =
#           |b| / 0.1, the S at which it is 0.1 SE; the log-log slope of
#           |bias| on S (C4 predicts -1 well below S_ref); and S_ladder, the
#           smallest S of the ladder from which every |bias| < 0.1 SE; b_se and
#           slope_se from complete units (.f5_units(), .f5_rate_se());
#   regimes per regime and ladder: the likelihood units U, the gates' S =
#           4 sqrt(U) (sqrt(U) / S = 0.25), the source family's own S and
#           the largest |bias| / SE of its fits there over the common
#           parameters (kind source), and the largest S_01 and S_ladder
#           over them, with the parameter.
f5_findings <- function(records) {
  recs <- Filter(function(r) is.list(r) && identical(r$family, "F5") && isTRUE(r$ok) &&
                   !is.null(r$est), records)
  if (!length(recs)) return(list())
  long <- do.call(rbind, lapply(recs, .f5_long))
  info <- do.call(rbind, lapply(recs, function(r) {
    data.frame(regime = r$cell, S_ref = r$diag$S_ref, units = r$diag$units,
               S_src = as.integer(r$diag$src[["S_src"]]))
  }))
  info <- info[!duplicated(info$regime), , drop = FALSE]
  key <- with(long, paste(regime, param, kind, S, sep = "\r"))
  by_S <- do.call(rbind, lapply(split(seq_len(nrow(long)), factor(key, unique(key))), function(i) {
    d <- long[i[1], ]
    data.frame(regime = d$regime, param = d$param, block = d$block, kind = d$kind, S = d$S,
               as.list(.f5_cell(long$z[i], long$rep[i])), stringsAsFactors = FALSE)
  }))
  com <- by_S[by_S$block != "asc", , drop = FALSE]
  asc <- by_S[by_S$block == "asc", , drop = FALSE]
  if (nrow(asc)) {
    k <- with(asc, paste(regime, kind, S, sep = "\r"))
    asc <- do.call(rbind, lapply(split(asc, factor(k, unique(k))), function(d) {
      rms <- function(x) sqrt(mean(x^2))
      data.frame(regime = d$regime[1], param = "(constants, RMS)", block = "asc", kind = d$kind[1],
                 S = d$S[1], datasets = max(d$datasets), bias_se = rms(d$bias_se),
                 bias_se_mcse = NA_real_, sd_se = rms(d$sd_se), sd_se_mcse = NA_real_)
    }))
  }
  by_S <- rbind(com, asc)
  by_S <- by_S[order(match(by_S$regime, F5_REGIMES$cell), by_S$kind, by_S$param, by_S$S), ]
  rownames(by_S) <- NULL
  lad <- com[com$kind != "source", , drop = FALSE]
  key <- with(lad, paste(regime, param, kind, sep = "\r"))
  rates <- do.call(rbind, lapply(split(lad, factor(key, unique(key))), function(d) {
    d <- d[order(d$S), , drop = FALSE]
    S_ref <- info$S_ref[match(d$regime[1], info$regime)]
    y <- d$bias_se
    m <- d$bias_se_mcse
    x <- 1 / d$S - 1 / S_ref
    wt <- is.finite(m) & m > 0
    u <- .f5_units(long, d)
    if (sum(wt) >= 2) {
      b <- sum(x[wt] * y[wt] / m[wt]^2) / sum(x[wt]^2 / m[wt]^2)
      cb <- numeric(length(x))
      cb[wt] <- (x[wt] / m[wt]^2) / sum(x[wt]^2 / m[wt]^2)
      b_se <- .f5_rate_se(u, cb)
    } else {
      b <- sum(x * y) / sum(x^2)
      b_se <- NA_real_
    }
    lw <- (y / m)^2
    k <- is.finite(log(abs(y))) & is.finite(lw)
    sl <- if (sum(k) >= 2) {
      X <- cbind(1, log(d$S[k]))
      f <- stats::lm.wfit(X, log(abs(y[k])), lw[k])
      a <- numeric(length(y))
      a[k] <- (solve(crossprod(X * sqrt(lw[k]))) %*% t(X * lw[k]))[2, ] / y[k]
      c(f$coefficients[[2]], .f5_rate_se(u, a))
    } else c(NA_real_, NA_real_)
    below <- abs(y) < 0.1
    ok_from <- rev(cumprod(rev(below))) == 1
    data.frame(regime = d$regime[1], param = d$param[1], kind = d$kind[1],
               S = paste(d$S, collapse = ","), b = b, b_se = b_se, S_01 = abs(b) / 0.1,
               slope = sl[1], slope_se = sl[2],
               S_ladder = if (any(ok_from)) min(d$S[ok_from]) else NA_integer_,
               stringsAsFactors = FALSE)
  }))
  rownames(rates) <- NULL
  rk <- paste(rates$regime, rates$kind)
  regimes <- do.call(rbind, lapply(split(rates, factor(rk, unique(rk))), function(d) {
    k <- match(d$regime[1], info$regime)
    U <- info$units[k]
    i <- if (any(is.finite(d$S_01))) which.max(d$S_01) else NA_integer_
    sl <- if (anyNA(d$S_ladder)) NA_integer_ else max(d$S_ladder)
    at <- com$bias_se[com$regime == d$regime[1] & com$kind == "source"]
    data.frame(regime = d$regime[1], kind = d$kind[1], units = U, S_gates = ceiling(4 * sqrt(U)),
               S_src = info$S_src[k], max_bias_at_S_src = if (length(at)) max(abs(at)) else NA,
               S_01_max = d$S_01[i], S_01_param = d$param[i],
               S_ladder_max = sl, S_ladder_param = if (is.na(sl)) d$param[is.na(d$S_ladder)][1]
                                                   else d$param[which.max(d$S_ladder)],
               stringsAsFactors = FALSE)
  }))
  rownames(regimes) <- NULL
  list(by_S = by_S, rates = rates, regimes = regimes)
}

family <- list(
  id = "F5",
  claims = "C4",
  dgp = "the source families' (F3: D4-40, F4: D1, F6: D3, F7: D3-claims, D4-census)",
  cells = .f5_cells[c("cell", "N", "src", "src_cell", "S_ref", "S_gen", "seeds", "S_store")],
  methods = list(msl = f5_msl),
  methods_by_cell = NULL,
  R = list(mc = stats::setNames(.f5_cells$R, .f5_cells$cell)),
  simulate = f5_simulate,
  estimands = f5_estimands,
  gates = NULL,
  findings = f5_findings
)
