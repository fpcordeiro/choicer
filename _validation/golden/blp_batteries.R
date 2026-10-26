# The share inversion's golden batteries and the solver tier (BLP engine plan,
# Phase A).
#
# The inversion's results may move within its tolerance (user decision U1,
# 2026-10-25), so its leaves are not checked for bit identity. Every BLP call
# of the harness goes through golden_blp(), which runs the call and records one
# leaf (class "golden_blp", kept whole by flatten()) with what the solver tier
# needs to judge the result on its own:
#   * the call's value, warnings (with the classed warning's $convergence),
#     messages, error and printed output;
#   * the residual at the returned delta, recomputed through the untouched
#     *_predict_shares() kernel: max |log S_j - log s_j| over the inside
#     alternatives, re-anchored on each component's reference without an
#     outside option, and each component's outside share's with one;
#   * a reference solution: the delta the targets were computed at, where a
#     call site knows it, or else Newton's on the same kernel from the call's
#     start, to the rounding floor; with kappa = ||B^-1||_inf, B the
#     normalized residual's Jacobian on the free constants (central
#     differences), and the error of the returned delta against it,
#     re-anchored per component;
#   * the engine's input checks that main did not make, read off the call's
#     arguments by the harness itself (gb_expect_error()): where one fails,
#     the call must stop with its message.
# gb_tier() judges a leaf that is not identical to main's (the comparison loop
# of refactor_identity_check.R calls it); gb_tier_summary() prints the tier's
# report. run_blp() is the battery of designs that make the inversion iterate.
# Sourced before the inline batteries, which use golden_blp().

`%||%` <- function(a, b) if (is.null(a)) b else a

# ---------------------------------------------------------------------------
# Running a call
# ---------------------------------------------------------------------------

# The function a call goes to: an export, or blp()'s method for the fit.
gb_fun <- function(fn, args) {
  if (identical(fn, "blp")) {
    return(utils::getS3method("blp", class(args[[1L]])[1L]))
  }
  get(fn, envir = asNamespace("choicer"))
}

# The call's arguments by name, with the function's defaults for the others
# (evaluated where the function would evaluate them).
gb_args <- function(fn, args) {
  f <- gb_fun(fn, args)
  mc <- match.call(f, as.call(c(list(as.name(fn)), args)))
  a <- as.list(mc)[-1L]
  fm <- formals(f)
  env <- list2env(a, parent = environment(f))
  for (nm in setdiff(names(fm), c(names(a), "..."))) {
    if (!identical(fm[[nm]], quote(expr = ))) {
      a[[nm]] <- eval(fm[[nm]], env)
      assign(nm, a[[nm]], envir = env)
    }
  }
  a
}

# Value, warnings ("class | message", and the classed warning's
# $convergence), messages, error ("class | message | call=TRUE/FALSE") and
# printed lines of a call, made by name (do.call() with a string, so that a
# fit's stored call never holds a closure), and its arguments: they are
# forced here, once, so that what forming them raises is the call's (and NULL
# comes back for them when that fails).
gb_run <- function(fn, args) {
  warns <- character()
  conv_warn <- NULL
  msgs <- character()
  err <- NA_character_
  val <- NULL
  forced <- NULL
  printed <- utils::capture.output(val <- tryCatch(
    withCallingHandlers({
      forced <- args
      do.call(fn, forced)
    },
      warning = function(w) {
        warns <<- c(warns, paste(class(w)[1L], conditionMessage(w), sep = " | "))
        if (inherits(w, "choicer_blp_unconverged")) conv_warn <<- w$convergence
        invokeRestart("muffleWarning")
      },
      message = function(m) {
        msgs <<- c(msgs, conditionMessage(m))
        invokeRestart("muffleMessage")
      }),
    error = function(e) {
      err <<- sprintf("%s | %s | call=%s", class(e)[1L], conditionMessage(e),
                      !is.null(conditionCall(e)))
      NULL
    }))
  list(value = val, warnings = warns, conv_warn = conv_warn, messages = msgs,
       error = err, stdout = printed, args = forced)
}

# ---------------------------------------------------------------------------
# The problem a call solves
# ---------------------------------------------------------------------------

# Components of the situation-alternative graph over the situations of
# positive weight: each alternative's label is its component's smallest code
# (the component's reference), found by propagating the smallest label; with
# each component's weight (mass, named by label) and each alternative's
# offer and captive weight (the situations offering it, and those offering
# only it; a code listed twice counts its situation once), over the total
# weight.
gb_components <- function(alt, M, w, J) {
  sit <- rep.int(seq_along(M), M)
  keep <- w[sit] > 0
  a <- alt[keep]
  s <- sit[keep]
  offered <- tabulate(a, J) > 0
  lab <- seq_len(J)
  if (length(a)) {
    repeat {
      sm <- tapply(lab[a], s, min)
      am <- tapply(as.numeric(sm[as.character(s)]), a, min)
      idx <- as.integer(names(am))
      new <- lab
      new[idx] <- pmin(lab[idx], as.integer(am))
      new <- new[new]
      if (identical(new, lab)) break
      lab <- new
    }
  }
  W <- sum(w)
  offer <- captive <- numeric(J)
  mass <- numeric()
  if (length(a)) {
    sets <- lapply(split(a, factor(s, levels = unique(s))), unique)
    ws <- w[as.integer(names(sets))]
    o <- tapply(rep(ws, lengths(sets)), unlist(sets), sum)
    offer[as.integer(names(o))] <- o
    one <- lengths(sets) == 1L
    if (any(one)) {
      cc <- tapply(ws[one], unlist(sets[one]), sum)
      captive[as.integer(names(cc))] <- cc
    }
    offer <- offer / W
    captive <- captive / W
    mass <- tapply(ws, lab[vapply(sets, `[`, 1L, 1L)], sum) / W
  }
  refs <- sort(unique(lab[offered]))
  list(lab = lab, offered = offered, ncomp = length(refs), refs = refs,
       mass = mass, offer = offer, captive = captive,
       size = tabulate(lab, J)[lab])
}

# The free constants of a full delta: without an outside option, the shares
# depend on differences only, so the kernels' theta takes delta - delta_1
# without its first entry.
gb_free <- function(d, outside) if (outside) d else (d - d[1L])[-1L]

# The full delta (J entries) of a theta whose constants follow its first n_head
# parameters: zeros without constants, the baseline's zero prepended without an
# outside option. Where a call site's targets are the shares at theta, this is
# the exact reference.
gb_delta_star <- function(theta, n_head, J, use_asc, outside) {
  if (!use_asc) return(rep(0, J))
  if (outside) theta[n_head + seq_len(J)] else c(0, theta[n_head + seq_len(J - 1L)])
}

# A fit's full delta, the reference of blp() on its own predicted shares.
gb_fit_delta <- function(fit) {
  outside <- isTRUE(fit$include_outside_option)
  J <- nrow(fit$alt_mapping) - outside
  asc <- fit$param_map$asc
  if (is.null(asc)) return(rep(0, J))
  d <- as.numeric(fit$coefficients[asc])
  if (outside) d else c(0, d)
}

# A design without covariates as one covariate of zeros (the prediction
# kernels want one; the utilities are the same, bit for bit).
gb_x0 <- function(X, beta) {
  if (NCOL(X) == 0L) list(X = matrix(0, NROW(X), 1L), beta = 0) else
    list(X = X, beta = as.numeric(beta))
}

# The shares function of a call (the outside option's share first when it is
# modelled) through the untouched prediction kernel, its targets rescaled to
# sum to 1 (as the engine rescales them within its tolerance), and the
# constants' layout. NULL where the call's inputs make no such problem.
gb_problem <- function(fn, a) {
  tryCatch(gb_problem_(fn, a), warning = function(w) NULL,
           error = function(e) NULL)
}

gb_problem_ <- function(fn, a) {
  if (identical(fn, "blp")) return(gb_problem_fit(a))
  outside <- isTRUE(as.logical(a$include_outside_option))
  alt <- a$alt_idx
  M <- a$M
  w <- as.numeric(a$weights)
  start <- as.numeric(a$delta)
  J <- if (identical(fn, "mxl_blp_contraction")) max(as.integer(alt)) else length(start)
  if (identical(fn, "mxl_blp_contraction") && !outside && length(start) == J - 1L) {
    start <- c(0, start)
  }
  xb <- gb_x0(a$X, a$beta)
  X <- xb$X
  beta <- xb$beta
  sfun <- switch(fn,
    blp_contraction = {
      function(d) as.numeric(choicer:::mnl_predict_shares(
        c(beta, gb_free(d, outside)), X, alt, M, w, TRUE, outside))
    },
    nl_blp_contraction = {
      nest <- a$nest_idx
      counts <- tabulate(as.integer(nest), max(as.integer(nest)))
      lam <- as.numeric(a$lambda)[counts != 1L]
      function(d) as.numeric(choicer:::nl_predict_shares(
        c(beta, lam, gb_free(d, outside)), X, alt, M, w, nest, TRUE, outside))
    },
    mxl_blp_contraction = {
      rc_mean <- isTRUE(a$rc_mean)
      # (with rc_mean = FALSE the solver takes mu as given: for normal
      # coefficients, the kernel's means; a log-normal one's has no match)
      if (!rc_mean && any(as.numeric(a$mu) != 0)) {
        if (any(a$rc_dist != 0L)) return(NULL)
        rc_mean <- TRUE
      }
      head <- c(beta, if (rc_mean) as.numeric(a$mu), as.numeric(a$L_params))
      W <- a$W; eta <- a$eta_draws
      function(d) as.numeric(choicer:::mxl_predict_shares(
        c(head, gb_free(d, outside)), X, W, alt, M, w, eta, a$rc_dist,
        isTRUE(a$rc_correlation), rc_mean, TRUE, outside, a$gen_seed,
        a$gen_scramble, a$gen_S))
    })
  gb_problem_common(sfun, a$target_shares, outside, J, alt, M, w, start,
                    a$tol, a$max_iter)
}

# blp() on a fit: the shares predict(type = "shares") evaluates, with the
# fit's draws, at the fit's parameters and a given full delta.
gb_problem_fit <- function(a) {
  fit <- a[[1L]]
  d <- fit[["data"]]
  if (is.null(d)) return(NULL)
  pm <- fit$param_map
  co <- fit$coefficients
  outside <- isTRUE(fit$include_outside_option)
  J <- nrow(fit$alt_mapping) - outside
  start <- a$delta_init
  if (is.null(start)) {
    start <- if (is.null(pm$asc)) rep(0, J) else if (outside) co[pm$asc] else c(0, co[pm$asc])
  }
  start <- as.numeric(start)
  w <- d$weights
  xb <- gb_x0(d$X, co[pm$beta])
  if (inherits(fit, "choicer_mxl")) {
    gp <- choicer:::.mxl_pred_draws(fit$draws_info, length(d$M), "blp()")
    head <- c(xb$beta, if (!is.null(pm$mu)) co[pm$mu], co[pm$sigma])
    sfun <- function(dl) as.numeric(choicer:::mxl_predict_shares(
      c(head, gb_free(dl, outside)), xb$X, d$W, d$alt_idx, d$M, w, gp$eta_draws,
      fit$rc_dist, fit$rc_correlation, fit$rc_mean, TRUE, outside, gp$gen_seed,
      gp$gen_scramble, gp$gen_S))
  } else if (inherits(fit, "choicer_nl")) {
    head <- c(xb$beta, co[pm$lambda])
    sfun <- function(dl) as.numeric(choicer:::nl_predict_shares(
      c(head, gb_free(dl, outside)), xb$X, d$alt_idx, d$M, w, d$nest_idx, TRUE,
      outside))
  } else {
    sfun <- function(dl) as.numeric(choicer:::mnl_predict_shares(
      c(xb$beta, gb_free(dl, outside)), xb$X, d$alt_idx, d$M, w, TRUE, outside))
  }
  gb_problem_common(sfun, a$target_shares, outside, J, d$alt_idx, d$M, w,
                    start, a$tol, a$max_iter)
}

gb_problem_common <- function(sfun, target, outside, J, alt, M, w, start, tol,
                              max_iter) {
  S <- as.numeric(target)
  if (length(S) != J + outside || any(!is.finite(S)) || any(S <= 0)) return(NULL)
  if (length(start) != J || any(!is.finite(start))) return(NULL)
  cp <- gb_components(as.integer(alt), as.integer(M), w, J)
  if (!all(cp$offered)) return(NULL)
  S <- S / sum(S)
  if (!outside) start <- start - start[1L]
  free <- if (outside) seq_len(J) else which(cp$lab != seq_len(J))
  S_in <- if (outside) S[-1L] else S
  list(sfun = sfun, S_in = S_in, S0 = if (outside) S[1L] else NA_real_,
       outside = outside, J = J, lab = cp$lab, ncomp = cp$ncomp, refs = cp$refs,
       mass = cp$mass, free = free, start = start,
       resid_scale = max(1, abs(log(S))),
       tol = as.numeric(tol), max_iter = as.numeric(max_iter))
}

# The normalized log-share residuals at a full delta (J entries; NA where the
# kernel's output does not fit the problem) and the stopping residual's two
# parts: the largest absolute inside residual (res_in) and, with an outside
# option, each component's outside share's (res_out); res is the larger.
# With two or more components a component's outside share is its weight less
# its inside shares, a difference that loses about eps * weight / share,
# relatively: slack bounds that loss in res_out alone, and an outside share
# below 1e3 eps of its component's weight leaves res_out unverifiable (NA).
gb_resid <- function(P, d) {
  # Without an outside option and without free constants (every component a
  # single alternative) each share is its component's weight: nothing to
  # solve, and the residual is 0 whatever the kernel would say
  if (!P$outside && !length(P$free)) {
    return(list(rn = rep(0, P$J), res = 0, res_in = 0, res_out = NA_real_, slack = 0))
  }
  s <- tryCatch(P$sfun(d), error = function(e) NULL)
  if (length(s) != P$J + P$outside) {
    return(list(rn = rep(NA_real_, P$J), res = NA_real_, res_in = NA_real_,
                res_out = NA_real_, slack = 0))
  }
  s_in <- if (P$outside) s[-1L] else s
  r <- log(P$S_in) - log(s_in)
  rn <- if (P$outside) r else r - r[P$lab]
  # (shares that are not finite give residuals that count as Inf)
  inf_if <- function(x) if (is.finite(x)) x else Inf
  res_in <- inf_if(max(abs(rn)))
  res_out <- NA_real_
  slack <- 0
  if (P$outside) {
    if (P$ncomp == 1L) {
      res_out <- inf_if(abs(log(P$S0) - log(s[1L])))
    } else {
      res_out <- 0
      for (cl in names(P$mass)) {
        inc <- P$lab == as.integer(cl)
        m <- P$mass[[cl]]
        S0c <- m - sum(P$S_in[inc])
        s0c <- m - sum(s_in[inc])
        if (isTRUE(min(S0c, s0c) <= 1e3 * .Machine$double.eps * m)) {
          res_out <- NA_real_  # lost to cancellation: not verifiable
          break
        }
        res_out <- max(res_out, inf_if(abs(log(S0c) - log(s0c))))
        slack <- max(slack, 16 * .Machine$double.eps * m * (1 / S0c + 1 / s0c))
      }
    }
  }
  if (!is.finite(slack)) slack <- Inf
  list(rn = rn, res = max(res_in, res_out, na.rm = TRUE), res_in = res_in,
       res_out = res_out, slack = slack)
}

# The normalized residuals' Jacobian on the free constants, by central
# differences.
gb_jacobian <- function(P, d) {
  k <- length(P$free)
  B <- matrix(NA_real_, k, k)
  for (l in seq_len(k)) {
    j <- P$free[l]
    h <- 1e-6 * max(1, abs(d[j]))
    up <- d; up[j] <- d[j] + h
    dn <- d; dn[j] <- d[j] - h
    B[, l] <- (gb_resid(P, up)$rn[P$free] - gb_resid(P, dn)$rn[P$free]) / (2 * h)
  }
  B
}

# Linear solves by LINPACK's QR (serial, so the references do not depend on
# how many threads a BLAS uses).
gb_solve <- function(B, b) qr.solve(B, b, tol = 1e-14)

# Newton on the free constants from the call's start, the components'
# references held at it, its steps at most 5 log units long (weakly linked
# markets otherwise send the first steps far out, where the residual can
# still fall while the constants diverge), with backtracking on the largest
# residual; stops at its rounding floor.
gb_newton <- function(P) {
  d <- P$start
  f <- gb_resid(P, d)
  best <- max(abs(f$rn[P$free]))
  if (!is.finite(best)) return(NULL)
  for (it in seq_len(200L)) {
    if (best <= 1e-15) break
    B <- gb_jacobian(P, d)
    dx <- tryCatch(gb_solve(B, -f$rn[P$free]), error = function(e) NULL)
    if (is.null(dx) || any(!is.finite(dx))) break
    if (max(abs(dx)) > 5) dx <- dx * (5 / max(abs(dx)))
    t <- 1
    moved <- FALSE
    for (ls in seq_len(40L)) {
      dn <- d
      dn[P$free] <- d[P$free] + t * dx
      fn <- gb_resid(P, dn)
      bn <- max(abs(fn$rn[P$free]))
      if (is.finite(bn) && bn < best) { moved <- TRUE; break }
      t <- t / 2
    }
    if (!moved) break
    d <- dn; f <- fn; best <- bn
  }
  d
}

# The reference: the given delta, or Newton's; its residual floor and kappa.
# NULL where Newton stalls above 1e-11 (or a given delta misses its targets
# by that much: golden_blp() notes it).
gb_reference <- function(P, ref) {
  d <- if (is.null(ref)) gb_newton(P) else as.numeric(ref)
  if (is.null(d) || length(d) != P$J) return(NULL)
  if (!P$outside) d <- d - d[1L]
  fl <- gb_resid(P, d)$res
  if (!is.finite(fl) || fl > 1e-11) return(NULL)
  if (!length(P$free)) return(list(ref = d, floor = fl, kappa = 0))
  B <- gb_jacobian(P, d)
  Bi <- tryCatch(gb_solve(B, diag(nrow(B))), error = function(e) NULL)
  kappa <- if (is.null(Bi) || any(!is.finite(Bi))) Inf else max(rowSums(abs(Bi)))
  list(ref = d, floor = fl, kappa = kappa)
}

# Error of a full delta against the reference, re-anchored per component
# without an outside option.
gb_error <- function(P, d, ref) {
  if (P$outside) d - ref else (d - d[P$lab]) - (ref - ref[P$lab])
}

# The engine's input checks that main did not make, in the engine's order,
# read off a call's arguments: the regular expression of the message of the
# first that fails, NA if none does (or the arguments are not the kind the
# checks read). Main's own checks come first in the engine, so where they
# fail the error is main's; these are the necessary conditions only.
gb_expect_error <- function(fn, a) {
  tryCatch(gb_expect_error_(fn, a), warning = function(w) NA_character_,
           error = function(e) NA_character_)
}

gb_expect_error_ <- function(fn, a) {
  if (is.null(a)) return(NA_character_)
  mxl <- identical(fn, "mxl_blp_contraction")
  if (identical(fn, "blp")) {
    fit <- a[[1L]]
    d <- fit[["data"]]
    if (is.null(d)) return(NA_character_)
    outside <- isTRUE(fit$include_outside_option)
    alt <- as.integer(d$alt_idx)
    M <- as.integer(d$M)
    w <- as.numeric(d$weights)
    J <- nrow(fit$alt_mapping) - outside
    start <- as.numeric(a$delta_init %||% rep(0, J))
    lam <- NULL
    # the solver's message (the mixed logit's or the others'), or blp()'s own
    # check in R, which comes first from Phase A3 on
    positive <- paste(if (inherits(fit, "choicer_mxl")) "non-positive entries" else
      "strictly positive", "must be finite and positive", sep = "|")
  } else {
    outside <- isTRUE(as.logical(a$include_outside_option))
    alt <- as.integer(a$alt_idx)
    M <- as.integer(a$M)
    w <- as.numeric(a$weights)
    start <- as.numeric(a$delta)
    J <- if (mxl) max(alt) else length(start)
    lam <- if (identical(fn, "nl_blp_contraction")) as.numeric(a$lambda)
    positive <- if (mxl) "Target shares contains non-positive entries" else
      "all target_shares must be strictly positive"
  }
  S <- as.numeric(a$target_shares)
  if (length(S) != J + outside) return(NA_character_)  # main's length check
  # Every entry positive, NaN and NA failing (main's checks let NaN through,
  # and the mixed logit's looked at the inside shares only)
  if (!isTRUE(all(S > 0))) return(positive)
  if (!is.null(lam) && !isTRUE(all(lam > 0 & is.finite(lam)))) {
    return("lambda must be positive and finite")
  }
  if (!all(is.finite(start))) return("delta must be finite")
  if (!isTRUE(all(w >= 0 & is.finite(w)))) return("weights must be finite and non-negative")
  cp <- gb_components(alt, M, w, J)
  if (!all(cp$offered)) {
    nowhere <- which(!cp$offered)
    shown <- paste(utils::head(nowhere, 10L), collapse = ", ")
    if (length(nowhere) > 10L) shown <- paste0(shown, ", \\.\\.\\.")
    return(sprintf("Alternatives? %s (is|are) offered in no choice situation of positive weight",
                   shown))
  }
  tol <- 1.4901161193847656e-08
  S_in <- if (outside) S[-1L] else S
  if (outside) {
    if (!(abs(sum(S) - 1) <= tol)) return("must sum to 1, the outside option's share included")
    S_in <- S_in / sum(S)
    if (cp$ncomp > 1L) {
      inside <- tapply(S_in, cp$lab, sum)[names(cp$mass)]
      if (any(!(cp$mass - inside > 0))) return("leave no share to the outside option")
    }
  } else {
    sums <- tapply(S_in, cp$lab, sum)[names(cp$mass)]
    if (any(!(abs(sums - cp$mass) <= tol * cp$mass))) {
      return(if (cp$ncomp == 1L) "must sum to 1; they sum to" else "situations weigh")
    }
    S_in <- S_in * as.numeric(cp$mass[as.character(cp$lab)] / sums[as.character(cp$lab)])
  }
  for (j in seq_len(J)) {
    if (!outside && cp$size[j] == 1L) next
    if (!(S_in[j] < cp$offer[j])) {
      return(sprintf(paste0("Alternative %d's target share \\(.*\\) is not below the ",
                            "weight of the choice situations that offer it \\("), j))
    }
    if (!outside && !(S_in[j] > cp$captive[j])) {
      return(sprintf(paste0("Alternative %d's target share \\(.*\\) is not above the ",
                            "weight of the choice situations that offer it alone"), j))
    }
  }
  NA_character_
}

# ---------------------------------------------------------------------------
# The recorder
# ---------------------------------------------------------------------------

# One BLP call's leaf. fn: "blp_contraction", "nl_blp_contraction",
# "mxl_blp_contraction" or "blp"; args: its arguments (formed inside the
# recorder); ref: the delta the targets were computed at (full length), NULL
# for a Newton reference, or FALSE for none; expect: "any" (judged by its
# outcome; a call that main solved to its tolerance must still converge),
# "converged" (must converge), "capped" (max_iter binds by design) or "error"
# (must stop; also wherever gb_expect_error() finds a failing check);
# expect_error: a regular expression a new error may match (beyond the
# failing check's message).
golden_blp <- function(fn, args, ref = NULL, expect = "any", expect_error = NULL) {
  run <- gb_run(fn, args)
  a <- if (is.null(run$args)) NULL else tryCatch(gb_args(fn, run$args), error = function(e) NULL)
  P <- if (is.null(a)) NULL else gb_problem(fn, a)
  auto <- gb_expect_error(fn, a)
  if (!is.na(auto)) expect <- "error"
  pattern <- c(expect_error, if (!is.na(auto)) auto)
  delta <- if (is.null(run$value)) NULL else as.numeric(run$value)
  printed_max <- any(grepl("Maximum iterations", run$stdout, fixed = TRUE))
  conv_attr <- if (is.null(run$value)) NULL else attr(run$value, "convergence")
  converged <- is.na(run$error) && !is.null(delta) && all(is.finite(delta)) &&
    !printed_max && !any(startsWith(run$warnings, "choicer_blp_unconverged")) &&
    !isFALSE(conv_attr$converged)
  # (degenerate designs make the kernels' logs warn; the call's own warnings
  # were recorded above)
  # (no reference where the inputs fail a check: the call must stop)
  R <- NULL
  ref_note <- NA_character_
  if (!is.null(P) && !isFALSE(ref) && is.na(auto)) {
    R <- suppressWarnings(gb_reference(P, ref))
    if (is.null(R) && !is.null(ref)) ref_note <- "the given reference misses its targets"
  }
  rr <- list(res = NA_real_, res_in = NA_real_, res_out = NA_real_, slack = 0)
  err <- NULL
  if (!is.null(P) && !is.null(delta) && length(delta) == P$J) {
    dd <- if (P$outside) delta else delta - delta[1L]
    rr <- suppressWarnings(gb_resid(P, dd))
    if (!is.null(R)) err <- gb_error(P, dd, R$ref)
  }
  leaf <- list(
    fn = fn, expect = expect,
    expect_error = if (length(pattern)) paste(pattern, collapse = "|") else NA_character_,
    value = run$value, warnings = run$warnings, conv_warn = run$conv_warn,
    messages = run$messages, stdout = run$stdout, error = run$error,
    has_problem = !is.null(P),
    tol = if (is.null(P)) NA_real_ else P$tol,
    max_iter = if (is.null(P)) NA_real_ else P$max_iter,
    J = if (is.null(P)) NA_integer_ else P$J,
    outside = if (is.null(P)) NA else P$outside,
    refs = if (is.null(P) || P$outside) NULL else P$refs,
    converged = converged, residual = rr$res, residual_in = rr$res_in,
    residual_out = rr$res_out, slack = rr$slack,
    resid_scale = if (is.null(P)) NA_real_ else P$resid_scale,
    start = if (is.null(P)) NULL else P$start,
    ref = R$ref, ref_floor = R$floor, kappa = R$kappa, ref_note = ref_note,
    err = err, err_max = if (is.null(err)) NA_real_ else max(abs(err)))
  leaf <- unserialize(serialize(leaf, NULL))
  attr(leaf, "golden_leaf") <- TRUE
  class(leaf) <- "golden_blp"
  leaf
}

# ---------------------------------------------------------------------------
# The solver tier
# ---------------------------------------------------------------------------

gb_unconverged_warning <- function(x) startsWith(x$warnings, "choicer_blp_unconverged")

# A leaf without its call site's expectations, for the identity test: a
# battery may refine what it expects of the engine without recording a new
# baseline.
gb_core <- function(x) {
  if (inherits(x, "golden_blp")) x[c("expect", "expect_error")] <- NULL
  x
}

# A value's attributes other than those blp() adds (row names and the
# convergence record, Phase A3).
gb_attrs <- function(v) {
  at <- attributes(v)
  at[setdiff(names(at), c("dimnames", "convergence"))]
}

# The residual a leaf's tol admits: tol and the rounding of log shares (the
# outside terms also admit their slack, gb_resid()).
gb_admits <- function(x) {
  x$tol * (1 + 1e-9) + 16 * .Machine$double.eps * x$resid_scale
}

# Whether a leaf's recomputed residual meets its tol: the inside terms at the
# rounding margin, the outside terms with their slack too; NA where neither
# part could be recomputed.
gb_within_tol <- function(x) {
  if (is.na(x$residual_in)) return(NA)
  ok <- is.finite(x$residual_in) && x$residual_in <= gb_admits(x)
  if (!is.na(x$residual_out)) {
    ok <- ok && is.finite(x$residual_out) && x$residual_out <= gb_admits(x) + x$slack
  }
  ok
}

# Judge a leaf n that is not identical to main's leaf b. Returns pass (TRUE or
# FALSE), why (the failed rules, or the outcome) and the metrics the summary
# reports.
gb_tier <- function(b, n) {
  out <- list(pass = TRUE, why = character(), res_ratio = NA_real_,
              bound_ratio = NA_real_, vs_main = NA_real_, vs_main_res = NA_real_)
  fail <- function(msg) {
    out$pass <<- FALSE
    out$why <<- c(out$why, msg)
  }
  if (!inherits(b, "golden_blp") || !inherits(n, "golden_blp")) {
    fail("not a golden_blp leaf on both sides")
    return(out)
  }
  # The predict kernels are untouched and the calls the same: references,
  # starts and settings agree bit for bit
  if (!identical(b$ref, n$ref) || !identical(b$kappa, n$kappa) ||
      !identical(b$ref_floor, n$ref_floor) || !identical(b$start, n$start) ||
      !identical(b$tol, n$tol) || !identical(b$max_iter, n$max_iter)) {
    fail("reference, start or settings changed")
  }
  if (length(n$stdout)) fail("printed output")
  if (!identical(b$messages, n$messages)) fail("messages changed")
  if (!identical(b$warnings[!gb_unconverged_warning(b)],
                 n$warnings[!gb_unconverged_warning(n)])) {
    fail("other warnings changed")
  }
  matches <- function(msg) !is.na(n$expect_error) && grepl(n$expect_error, msg)
  if (!is.na(n$error)) {
    if (identical(n$expect, "error")) {
      if (!matches(n$error)) fail("error, but not the failing check's")
    } else if (!is.na(b$error)) {
      if (!identical(b$error, n$error) && !matches(n$error)) fail("error message changed")
    } else if (!matches(n$error)) {
      fail("new error")
    }
    out$why <- c(out$why, "error")
    return(out)
  }
  if (identical(n$expect, "error")) {
    fail(sprintf("no error, though the inputs fail a check (%s)", n$expect_error))
    return(out)
  }
  if (!isTRUE(n$has_problem)) {
    fail("no problem recorded for the call (the harness cannot judge it)")
    return(out)
  }
  v <- n$value
  if (is.null(v) || !is.double(v) || !identical(dim(v), c(as.integer(n$J), 1L))) {
    fail("value is not a J x 1 double matrix")
    return(out)
  }
  if (!is.null(b$value) && !identical(gb_attrs(b$value), gb_attrs(v))) {
    fail("value attributes changed")
  }
  delta <- as.numeric(v)
  # Without an outside option each component's reference is held at its
  # start (the first constant at 0)
  if (!isTRUE(n$outside) && length(n$refs) &&
      !identical(delta[n$refs], n$start[n$refs])) {
    fail("a component's reference moved")
  }
  # NA: no residual could be recomputed (the prediction kernels refuse the
  # call's inputs), and the leaf is reported as unverified; the outside terms
  # alone can be unverifiable (an outside share lost to cancellation)
  checkable <- !is.na(n$residual_in)
  if (!checkable) out$why <- c(out$why, "residual not recomputable")
  if (checkable && isTRUE(n$outside) && is.na(n$residual_out)) {
    out$why <- c(out$why, "outside residual not verifiable")
  }
  if (n$converged) {
    out$res_ratio <- n$residual / n$tol
    if (checkable && !isTRUE(gb_within_tol(n))) fail("residual above tol")
    if (!is.null(n$ref)) {
      bound <- 2 * n$kappa * (n$tol + n$ref_floor) +
        10 * .Machine$double.eps * pmax(1, abs(n$ref))
      out$bound_ratio <- max(abs(n$err) / bound)
      if (!is.finite(out$bound_ratio) || out$bound_ratio > 1) fail("outside the reference bound")
      if (isTRUE(b$converged) && is.finite(b$err_max) && b$err_max > 0) {
        out$vs_main <- n$err_max / b$err_max
      }
    }
    out$why <- c(out$why, "converged")
  } else {
    # main's own flag trusts its step criterion (the damped step under tol,
    # then one more step); it counts here where main's residual is within a
    # hundred times tol
    main_reached <- isTRUE(b$converged) && is.finite(b$residual) &&
      b$residual <= 100 * gb_admits(b)
    if (main_reached && !identical(n$expect, "capped")) fail("regression: main converged")
    if (identical(n$expect, "converged")) fail("expected to converge")
    if (!any(gb_unconverged_warning(n))) fail("no choicer_blp_unconverged warning")
    cv <- n$conv_warn
    ok_cv <- is.list(cv) && is.numeric(cv$evaluations) && length(cv$evaluations) == 1L &&
      !is.na(cv$evaluations) && is.numeric(cv$residual) && length(cv$residual) == 1L
    if (!ok_cv) {
      fail("the warning carries no $convergence record")
    } else {
      if (is.finite(n$max_iter) && cv$evaluations > n$max_iter) fail("evaluations above max_iter")
      if (identical(cv$status, "max_iter") && is.finite(n$max_iter) &&
          cv$evaluations != n$max_iter) {
        fail("stopped for max_iter before using it")
      }
      if (cv$evaluations > 0) {
        if (checkable) {
          # the engine's residual against the recomputed one: equal to
          # rounding, or, where the outside terms are not verifiable, at
          # least the inside part
          r1 <- cv$residual
          r2 <- n$residual
          same <- if (isTRUE(n$outside) && is.na(n$residual_out)) {
            is.finite(r1) && r1 >= n$residual_in - 1e-12 * max(1, n$residual_in)
          } else {
            (!is.finite(r1) && !is.finite(r2)) ||
              (is.finite(r1) && is.finite(r2) &&
               abs(r1 - r2) <= 1e-12 * max(1, r2) + n$slack)
          }
          if (!same) fail("reported residual differs from the recomputed one")
        }
      } else if (!identical(delta, n$start)) {
        # no evaluation: the start, normalized, comes back
        fail("no evaluation, yet the start was not returned")
      }
    }
    if (isFALSE(b$converged) && is.finite(b$residual) && b$residual > 0 &&
        is.finite(n$residual)) {
      out$vs_main_res <- n$residual / b$residual
    }
    out$why <- c(out$why, "unconverged")
  }
  out
}

# A leaf's field, NULL for a leaf that is not a list.
gb_field <- function(x, nm) if (is.list(x)) x[[nm]] else NULL

# The tier's report over the leaves (rows of the comparison table with their
# baseline and new leaves).
gb_tier_summary <- function(keys, base_flat, new_flat, results) {
  rows <- lapply(seq_along(keys), function(i) {
    b <- base_flat[[keys[i]]]
    n <- new_flat[[keys[i]]]
    r <- results[[i]]
    has_ref <- !is.null(gb_field(n, "ref"))
    tol <- gb_field(n, "tol") %||% NA_real_
    data.frame(entry = keys[i], identical = isTRUE(r$identical), pass = isTRUE(r$pass),
               why = paste(r$why, collapse = "; "),
               main_conv = isTRUE(gb_field(b, "converged")),
               new_conv = isTRUE(gb_field(n, "converged")),
               main_error = !is.na(gb_field(b, "error") %||% NA),
               new_error = !is.na(gb_field(n, "error") %||% NA),
               has_ref = has_ref,
               ref_note = !is.na(gb_field(n, "ref_note") %||% NA),
               unverified = any(grepl("not recomputable", r$why, fixed = TRUE)),
               illcond = has_ref && is.finite(tol) && 2 * n$kappa * tol > 1e-4,
               res_ratio = r$res_ratio %||% NA_real_,
               bound_ratio = r$bound_ratio %||% NA_real_,
               vs_main = r$vs_main %||% NA_real_,
               vs_main_res = r$vs_main_res %||% NA_real_,
               main_only = isTRUE(gb_field(b, "converged")) &&
                 !isTRUE(gb_field(n, "converged")) && is.na(gb_field(n, "error") %||% NA),
               stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

# ---------------------------------------------------------------------------
# Designs
# ---------------------------------------------------------------------------

# Choice sets over J alternatives in K regions: situation i belongs to region
# ((i - 1) mod K) + 1 and takes m of its alternatives (all of them for m at
# least the region's size), plus one from another region with probability
# `border`. The first situation of each region offers the whole region, so
# every alternative is offered.
gb_sets <- function(N, J, m, K = 1L, border = 0) {
  reg <- ((seq_len(J) - 1L) %% K) + 1L
  by_reg <- split(seq_len(J), reg)
  sets <- vector("list", N)
  for (i in seq_len(N)) {
    k <- ((i - 1L) %% K) + 1L
    own <- by_reg[[k]]
    pick <- if (i <= K || length(own) <= m) own else sort(own[sample.int(length(own), m)])
    if (border > 0 && K > 1L && stats::runif(1L) < border) {
      others <- seq_len(J)[reg != k]
      pick <- sort(unique(c(pick, others[sample.int(length(others), 1L)])))
    }
    sets[[i]] <- pick
  }
  list(alt_idx = as.integer(unlist(sets)), M = lengths(sets), reg = reg)
}

gb_design <- function(seed, N, J, m, K = 1L, border = 0, Kx = 2L, Kw = 0L) {
  set.seed(seed)
  s <- gb_sets(N, J, m, K, border)
  n <- length(s$alt_idx)
  list(X = matrix(stats::rnorm(n * Kx), n, Kx), alt_idx = s$alt_idx, M = s$M,
       weights = stats::runif(N, 0.5, 2),
       W = if (Kw > 0L) matrix(stats::rnorm(n * Kw), n, Kw) else NULL,
       J = J, N = N, reg = s$reg)
}

# ---------------------------------------------------------------------------
# The battery
# ---------------------------------------------------------------------------

run_blp <- function() {
  out <- list()
  beta <- c(0.8, -0.5)
  # Messages of the engine's input checks (the leaves match them only where
  # main did not stop)
  E_SUM <- "must sum to 1|situations weigh"
  E_OFFER <- "offered in no choice situation"
  E_BOUND <- "situations that offer it"
  E_DELTA <- "delta must be finite"
  E_LAMBDA <- "lambda must be positive"
  E_DAMP <- "damping"
  E_SCALAR <- "`tol`|`max_iter`"

  # -- Multinomial logit: the export ------------------------------------------
  mnl <- function(d, delta, target, ioo, ...) {
    c(list(delta = delta, target_shares = target, X = d$X, beta = beta,
           alt_idx = d$alt_idx, M = d$M, weights = d$weights,
           include_outside_option = ioo), list(...))
  }
  mnl_shares <- function(d, ds, ioo, w = d$weights) {
    as.numeric(choicer:::mnl_predict_shares(c(beta, gb_free(ds, ioo)), d$X,
                                            d$alt_idx, d$M, w, TRUE, ioo))
  }
  designs <- list(
    full = gb_design(51001L, 300L, 8L, 8L),
    regional = gb_design(51002L, 600L, 12L, 3L, K = 3L, border = 0.05))
  for (dn in names(designs)) {
    d <- designs[[dn]]
    for (ioo in c(FALSE, TRUE)) {
      set.seed(51100L + ioo)
      ds <- stats::rnorm(d$J)
      if (!ioo) ds <- ds - ds[1L]
      S <- mnl_shares(d, ds, ioo)
      tag <- sprintf("mnl|%s|oo=%d", dn, ioo)
      out[[paste0(tag, "|known")]] <- golden_blp("blp_contraction",
        mnl(d, rep(0, d$J), S, ioo), ref = ds, expect = "converged")
      out[[paste0(tag, "|known|tol12")]] <- golden_blp("blp_contraction",
        mnl(d, rep(0, d$J), S, ioo, tol = 1e-12), ref = ds, expect = "converged")
      out[[paste0(tag, "|known|far")]] <- golden_blp("blp_contraction",
        mnl(d, rep(c(4, -4), length.out = d$J), S, ioo), ref = ds, expect = "converged")
      out[[paste0(tag, "|max_iter=0")]] <- golden_blp("blp_contraction",
        mnl(d, rep(0, d$J), S, ioo, max_iter = 0L), ref = FALSE, expect = "capped")
      out[[paste0(tag, "|max_iter=3")]] <- golden_blp("blp_contraction",
        mnl(d, rep(0, d$J), S, ioo, max_iter = 3L), ref = FALSE, expect = "capped")
      if (dn == "full") {
        out[[paste0(tag, "|symmetrized")]] <- golden_blp("blp_contraction",
          mnl(d, rep(0, d$J), (S + rev(S)) / 2, ioo))
        out[[paste0(tag, "|uniform")]] <- golden_blp("blp_contraction",
          mnl(d, rep(0, d$J), rep(1 / (d$J + ioo), d$J + ioo), ioo))
        out[[paste0(tag, "|tiny")]] <- {
          dt <- ds; dt[d$J] <- -22
          golden_blp("blp_contraction", mnl(d, rep(0, d$J), mnl_shares(d, dt, ioo), ioo),
                     ref = dt, expect = "converged")
        }
        out[[paste0(tag, "|sum_off")]] <- golden_blp("blp_contraction",
          mnl(d, rep(0, d$J), S * 1.1, ioo), ref = FALSE, expect = "error",
          expect_error = E_SUM)
        out[[paste0(tag, "|na_delta")]] <- golden_blp("blp_contraction",
          mnl(d, replace(rep(0, d$J), 2L, NA), S, ioo), ref = FALSE, expect = "error",
          expect_error = E_DELTA)
        out[[paste0(tag, "|underflow_start")]] <- golden_blp("blp_contraction",
          mnl(d, replace(rep(0, d$J), d$J, -1000), S, ioo), ref = FALSE)
        out[[paste0(tag, "|tol_neg")]] <- golden_blp("blp_contraction",
          mnl(d, rep(0, d$J), S, ioo, tol = -1), ref = FALSE, expect = "error",
          expect_error = E_SCALAR)
        out[[paste0(tag, "|max_iter_na")]] <- golden_blp("blp_contraction",
          mnl(d, rep(0, d$J), S, ioo, max_iter = NA_integer_), ref = FALSE,
          expect = "error", expect_error = E_SCALAR)
      }
    }
  }
  # A zero-weight situation among others, and an alternative only a
  # zero-weight situation offers
  d <- designs$full
  w0 <- replace(d$weights, 5L, 0)
  ds <- c(0, stats::rnorm(d$J - 1L))
  a0 <- mnl(d, rep(0, d$J), mnl_shares(d, ds, FALSE, w0), FALSE)
  a0$weights <- w0
  out[["mnl|zero_weight"]] <- golden_blp("blp_contraction", a0, ref = ds,
                                         expect = "converged")
  d9 <- d
  d9$alt_idx <- c(d$alt_idx[seq_len(sum(d$M[1:4]))], 9L,
                  d$alt_idx[-seq_len(sum(d$M[1:4]))])
  d9$M <- d$M
  d9$M[5L] <- d9$M[5L] + 1L
  d9$X <- rbind(d$X[seq_len(sum(d$M[1:4])), , drop = FALSE], c(0.1, 0.2),
                d$X[-seq_len(sum(d$M[1:4])), , drop = FALSE])
  out[["mnl|offered_by_zero_weight"]] <- golden_blp("blp_contraction",
    list(delta = rep(0, 9L), target_shares = rep(1 / 9, 9L), X = d9$X, beta = beta,
         alt_idx = d9$alt_idx, M = d9$M, weights = w0, include_outside_option = FALSE),
    ref = FALSE, expect = "error", expect_error = E_OFFER)
  # Two components (regions without a border): feasible targets, and targets
  # whose components' sums miss their situations' weight; with an outside
  # option, each component's outside share enters the residual
  d2 <- gb_design(51003L, 400L, 8L, 4L, K = 2L, border = 0)
  set.seed(51104L)
  ds <- c(0, stats::rnorm(7L))
  S2 <- mnl_shares(d2, ds, FALSE)
  out[["mnl|two_components"]] <- golden_blp("blp_contraction",
    mnl(d2, rep(0, 8L), S2, FALSE), ref = ds, expect = "converged")
  out[["mnl|two_components|start"]] <- golden_blp("blp_contraction",
    mnl(d2, c(0, 1, 0.5, -0.5, 2, 1, 0, -1), S2, FALSE), ref = ds, expect = "converged")
  # region 1's targets up by a fifth, region 2's down to keep the total at 1
  in1 <- d2$reg == 1L
  bad <- S2
  bad[in1] <- S2[in1] * 1.2
  bad[!in1] <- S2[!in1] * (1 - sum(bad[in1])) / sum(S2[!in1])
  out[["mnl|two_components|infeasible"]] <- golden_blp("blp_contraction",
    mnl(d2, rep(0, 8L), bad, FALSE), ref = FALSE, expect = "error", expect_error = E_SUM)
  S2o <- mnl_shares(d2, ds - 1, TRUE)
  out[["mnl|two_components|oo=1"]] <- golden_blp("blp_contraction",
    mnl(d2, rep(0, 8L), S2o, TRUE), ref = ds - 1, expect = "converged")
  out[["mnl|two_components|oo=1|max_iter=3"]] <- golden_blp("blp_contraction",
    mnl(d2, rep(0, 8L), S2o, TRUE, max_iter = 3L), ref = FALSE, expect = "capped")
  # Repeated codes (a situation listing an alternative twice), both layouts
  dr <- designs$full
  r0 <- c(0L, cumsum(dr$M))
  for (t in seq(2L, 60L, by = 3L)) dr$alt_idx[r0[t] + 1L] <- dr$alt_idx[r0[t] + 2L]
  for (ioo in c(FALSE, TRUE)) {
    set.seed(51200L + ioo)
    ds <- stats::rnorm(dr$J)
    if (!ioo) ds <- ds - ds[1L]
    out[[sprintf("mnl|repeated|oo=%d", ioo)]] <- golden_blp("blp_contraction",
      mnl(dr, rep(0, dr$J), mnl_shares(dr, ds, ioo), ioo), ref = ds, expect = "converged")
  }
  # An alternative's target above the weight of the situations offering it
  # (regional sets with a border: alternative 12 is offered by its region)
  d <- designs$regional
  offer12 <- sum(d$weights[vapply(split(d$alt_idx, rep.int(seq_along(d$M), d$M)),
                                  function(a) 12L %in% a, TRUE)]) / sum(d$weights)
  tg <- rep((1 - 1.05 * offer12) / 11, 12L)
  tg[12L] <- 1.05 * offer12
  out[["mnl|offer_bound"]] <- golden_blp("blp_contraction",
    mnl(d, rep(0, 12L), tg, FALSE), ref = FALSE, expect = "error", expect_error = E_BOUND)

  # -- Nested logit: the export -----------------------------------------------
  nest_idx <- rep(1:3, each = 3L)
  dn <- gb_design(52001L, 500L, 9L, 9L)
  nl <- function(delta, target, lambda, ioo, ...) {
    c(list(delta = delta, target_shares = target, X = dn$X, beta = beta,
           lambda = lambda, alt_idx = dn$alt_idx, nest_idx = nest_idx, M = dn$M,
           weights = dn$weights, include_outside_option = ioo), list(...))
  }
  for (ioo in c(FALSE, TRUE)) {
    for (lam in c(0.3, 0.9, 1.4)) {
      set.seed(52100L + 10L * ioo + round(10 * lam))
      ds <- stats::rnorm(9L)
      if (!ioo) ds <- ds - ds[1L]
      S <- as.numeric(choicer:::nl_predict_shares(c(beta, rep(lam, 3L), gb_free(ds, ioo)),
                                                  dn$X, dn$alt_idx, dn$M, dn$weights,
                                                  nest_idx, TRUE, ioo))
      for (damp in c(1, 0.6, 0.3)) {
        out[[sprintf("nl|oo=%d|lambda=%g|damping=%g", ioo, lam, damp)]] <-
          golden_blp("nl_blp_contraction", nl(rep(0, 9L), S, rep(lam, 3L), ioo,
                                              damping = damp), ref = ds,
                     expect = "converged")
      }
      if (lam == 0.3) {
        out[[sprintf("nl|oo=%d|max_iter=3", ioo)]] <- golden_blp("nl_blp_contraction",
          nl(rep(0, 9L), S, rep(lam, 3L), ioo, max_iter = 3L), ref = FALSE,
          expect = "capped")
      }
    }
  }
  # Nests of six with lambda = 3 (where the exact diagonal overshoots a nest's
  # common shift), and repeated codes, both layouts
  dn6 <- gb_design(52002L, 600L, 12L, 12L)
  nest6 <- rep(1:2, each = 6L)
  dn6r <- dn6
  r0 <- c(0L, cumsum(dn6$M))
  for (t in seq(2L, 90L, by = 4L)) dn6r$alt_idx[r0[t] + 1L] <- dn6r$alt_idx[r0[t] + 2L]
  for (ioo in c(FALSE, TRUE)) {
    for (dd in list(list(tag = "big_lambda", d = dn6, lam = c(3, 3)),
                    list(tag = "repeated", d = dn6r, lam = c(0.4, 1.5)))) {
      set.seed(52200L + ioo)
      ds <- stats::rnorm(12L)
      if (!ioo) ds <- ds - ds[1L]
      d <- dd$d
      S <- as.numeric(choicer:::nl_predict_shares(c(beta, dd$lam, gb_free(ds, ioo)),
                                                  d$X, d$alt_idx, d$M, d$weights,
                                                  nest6, TRUE, ioo))
      out[[sprintf("nl|%s|oo=%d", dd$tag, ioo)]] <- golden_blp("nl_blp_contraction",
        list(delta = rep(0, 12L), target_shares = S, X = d$X, beta = beta,
             lambda = dd$lam, alt_idx = d$alt_idx, nest_idx = nest6, M = d$M,
             weights = d$weights, include_outside_option = ioo),
        ref = ds, expect = "converged")
    }
  }
  S <- rep(1 / 9, 9L)
  out[["nl|lambda_zero"]] <- golden_blp("nl_blp_contraction",
    nl(rep(0, 9L), S, c(0.5, 0, 0.5), FALSE, max_iter = 10L), ref = FALSE,
    expect = "error", expect_error = E_LAMBDA)
  out[["nl|damping_zero"]] <- golden_blp("nl_blp_contraction",
    nl(rep(0, 9L), S, rep(0.5, 3L), FALSE, damping = 0, max_iter = 10L), ref = FALSE,
    expect = "error", expect_error = E_DAMP)
  out[["nl|damping_big"]] <- golden_blp("nl_blp_contraction",
    nl(rep(0, 9L), S, rep(0.5, 3L), FALSE, damping = 1.5, max_iter = 10L), ref = FALSE,
    expect = "error", expect_error = E_DAMP)

  # -- Mixed logit: the export ------------------------------------------------
  dm <- gb_design(53001L, 240L, 6L, 6L, Kw = 2L)
  Sd <- 20L
  cube <- get_halton_normals(Sd, length(dm$M), 2L)
  draws <- list(
    cube = list(eta_draws = cube, gen_seed = -1L, gen_scramble = 1L, gen_S = 0L),
    store = list(eta_draws = array(0, c(2L, 0L, 0L)), gen_seed = 0L, gen_scramble = 2L,
                 gen_S = Sd),
    generate = list(eta_draws = array(0, c(2L, 0L, 0L)), gen_seed = 7L,
                    gen_scramble = 1L, gen_S = Sd))
  mu <- c(0.4, -0.3)
  Lp <- c(log(0.8), 0.3, log(0.6))
  for (rc in list(normal = c(0L, 0L), lognormal = c(0L, 1L))) {
    for (dr in names(draws)) {
      for (ioo in c(FALSE, TRUE)) {
        set.seed(53100L + 2L * ioo + rc[2L])
        ds <- stats::rnorm(6L)
        if (!ioo) ds <- ds - ds[1L]
        g <- draws[[dr]]
        S <- as.numeric(choicer:::mxl_predict_shares(
          c(beta, mu, Lp, gb_free(ds, ioo)), dm$X, dm$W, dm$alt_idx, dm$M,
          dm$weights, g$eta_draws, rc, TRUE, TRUE, TRUE, ioo, g$gen_seed,
          g$gen_scramble, g$gen_S))
        args <- function(delta, ...) c(list(
          delta = delta, target_shares = S, X = dm$X, W = dm$W, beta = beta, mu = mu,
          L_params = Lp, alt_idx = dm$alt_idx, M = dm$M, weights = dm$weights,
          eta_draws = g$eta_draws, rc_dist = rc, rc_correlation = TRUE,
          rc_mean = TRUE, include_outside_option = ioo), list(...),
          list(gen_seed = g$gen_seed, gen_scramble = g$gen_scramble, gen_S = g$gen_S))
        tag <- sprintf("mxl|%s|%s|oo=%d", if (rc[2L]) "lognormal" else "normal", dr, ioo)
        out[[tag]] <- golden_blp("mxl_blp_contraction", args(rep(0, 6L)), ref = ds,
                                 expect = "converged")
        if (!ioo) {
          out[[paste0(tag, "|free")]] <- golden_blp("mxl_blp_contraction",
            args(rep(0, 5L)), ref = ds, expect = "converged")
          # the free start is the full one with its leading zero
          out[[paste0(tag, "|same_free")]] <- identical(
            as.numeric(out[[tag]]$value), as.numeric(out[[paste0(tag, "|free")]]$value))
        }
        if (dr == "store" && rc[2L] == 0L) {
          out[[paste0(tag, "|max_iter=0")]] <- golden_blp("mxl_blp_contraction",
            args(rep(0, 6L), max_iter = 0L), ref = FALSE, expect = "capped")
          out[[paste0(tag, "|max_iter=3")]] <- golden_blp("mxl_blp_contraction",
            args(rep(0, 6L), max_iter = 3L), ref = FALSE, expect = "capped")
        }
      }
    }
  }

  # -- blp() on fits ----------------------------------------------------------
  for (ioo in c(FALSE, TRUE)) {
    s <- data.table::as.data.table(simulate_mnl_data(
      N = 300L, J = 5L, beta = c(0.8, -0.6), delta = c(0.4, -0.3, 0.2, 0.1, -0.2),
      outside_option = TRUE, seed = 54001L)$data)
    if (ioo) s <- s[alt != 0L]
    s[, nest := ifelse(alt <= 2L, 1L, 2L)]
    oo_args <- if (ioo) {
      list(include_outside_option = TRUE, outside_opt_label = 0L)
    } else {
      list(outside_opt_label = 0L)
    }
    for (ua in c(TRUE, FALSE)) {
      common <- c(list(data = s, id_col = "id", alt_col = "alt", choice_col = "choice",
                       covariate_cols = c("x1", "x2"), use_asc = ua), oo_args)
      fits <- list(
        mnl = suppressWarnings(suppressMessages(do.call("run_mnlogit", common))),
        nl = suppressWarnings(suppressMessages(do.call("run_nestlogit",
                                                        c(common, nest_col = "nest")))))
      for (m in names(fits)) {
        fit <- fits[[m]]
        J <- nrow(fit$alt_mapping) - ioo
        pm <- fit$param_map
        d0 <- if (is.null(pm$asc)) rep(0, J) else if (ioo) fit$coefficients[pm$asc] else
          c(0, fit$coefficients[pm$asc])
        set.seed(54100L + ioo + 2L * ua)
        ds <- as.numeric(d0) + stats::rnorm(J, 0, 0.5)
        if (!ioo) ds <- ds - ds[1L]
        P <- gb_problem("blp", list(fit, target_shares = rep(1 / (J + ioo), J + ioo),
                                    tol = 1e-8, max_iter = 1000))
        tag <- sprintf("fit|%s|oo=%d|asc=%d", m, ioo, ua)
        if (is.null(P)) {
          # (a harness fault: noted, and the run goes on to its report)
          warning("run_blp(): no shares function for ", tag, call. = FALSE)
          out[[paste0(tag, "|no_problem")]] <- TRUE
          next
        }
        S <- P$sfun(ds)
        out[[paste0(tag, "|known")]] <- golden_blp("blp", list(fit, target_shares = S),
                                                   ref = ds, expect = "converged")
        tg <- as.numeric(predict(fit, type = "shares"))
        out[[paste0(tag, "|symmetrized")]] <- golden_blp("blp",
          list(fit, target_shares = (tg + rev(tg)) / 2))
        out[[paste0(tag, "|uniform")]] <- golden_blp("blp",
          list(fit, target_shares = rep(1 / (J + ioo), J + ioo)))
        if (ua) {
          # (without constants the fit's start is zeros already)
          out[[paste0(tag, "|start0")]] <- golden_blp("blp",
            list(fit, target_shares = S, delta_init = rep(0, J)), ref = ds,
            expect = "converged")
          out[[paste0(tag, "|max_iter=0")]] <- golden_blp("blp",
            list(fit, target_shares = S, max_iter = 0), ref = FALSE, expect = "capped")
          out[[paste0(tag, "|max_iter=3")]] <- golden_blp("blp",
            list(fit, target_shares = S, max_iter = 3), ref = FALSE, expect = "capped")
        }
        if (m == "nl") {
          out[[paste0(tag, "|damping")]] <- golden_blp("blp",
            list(fit, target_shares = S, damping = 0.6), ref = ds)
        }
      }
    }
  }
  # Mixed logit fits, stored and generated draws, the draws kept or not
  for (mode in c("store", "generate")) {
    for (ioo in c(FALSE, TRUE)) {
      s <- data.table::as.data.table(simulate_mxl_data(
        N = 200L, J = 4L, beta = c(0.8, -0.6), delta = c(0.4, -0.3, 0.2, 0.1),
        outside_option = TRUE, seed = 55001L)$data)
      if (ioo) s <- s[alt != 0L]
      oo_args <- if (ioo) {
        list(include_outside_option = TRUE, outside_opt_label = 0L)
      } else {
        list(outside_opt_label = 0L)
      }
      fit <- suppressWarnings(suppressMessages(do.call("run_mxlogit", c(list(
        data = s, id_col = "id", alt_col = "alt", choice_col = "choice",
        covariate_cols = c("x1", "x2"), random_var_cols = c("w1", "w2"), S = 20L,
        draws = mode, seed = 55002L, control = list(maxeval = 40L)), oo_args))))
      J <- nrow(fit$alt_mapping) - ioo
      pm <- fit$param_map
      d0 <- if (ioo) fit$coefficients[pm$asc] else c(0, fit$coefficients[pm$asc])
      set.seed(55100L + ioo)
      ds <- as.numeric(d0) + stats::rnorm(J, 0, 0.5)
      if (!ioo) ds <- ds - ds[1L]
      P <- gb_problem("blp", list(fit, target_shares = rep(1 / (J + ioo), J + ioo),
                                  tol = 1e-8, max_iter = 1000))
      tag <- sprintf("fit|mxl|%s|oo=%d", mode, ioo)
      if (is.null(P)) {
        warning("run_blp(): no shares function for ", tag, call. = FALSE)
        out[[paste0(tag, "|no_problem")]] <- TRUE
        next
      }
      S <- P$sfun(ds)
      for (kb in c(0, Inf)) {
        out[[sprintf("%s|keep=%g", tag, kb)]] <- golden_blp("blp",
          list(fit, target_shares = S, keep_draws_bytes = kb), ref = ds,
          expect = "converged")
      }
      # kept draws are the values each evaluation would form
      out[[paste0(tag, "|same_keep")]] <- identical(
        as.numeric(out[[paste0(tag, "|keep=0")]]$value),
        as.numeric(out[[paste0(tag, "|keep=Inf")]]$value))
      out[[paste0(tag, "|max_iter=3")]] <- golden_blp("blp",
        list(fit, target_shares = S, max_iter = 3), ref = FALSE, expect = "capped")
    }
  }
  out
}
