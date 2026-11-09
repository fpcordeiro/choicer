# --- Share-matching constants (concentrate_out_asc = "blp") ------------------

#' The share-matching constants' variance
#'
#' The variance of `run_mxlogit(concentrate_out_asc = "blp")`'s estimates.
#' The constants delta solve the share equations sum_u w_u c_u(theta_c,
#' delta) = 0, and theta_c maximizes the log-likelihood along them, so the
#' estimate solves the stacked equations
#'   sum_u w_u m_u = 0,  m_u = s_c,u + D' s_d,u   (the projected score)
#'   sum_u w_u c_u = 0                             (the share residuals)
#' over the likelihood units u, with D = d delta / d theta_c' = -C_d^{-1} C_c
#' (C the Jacobian of the weighted predicted counts, C_d on the free
#' constants, C_c on theta_c). Linearized at the truth, with I the
#' information (the negated Hessian of the weighted log-likelihood, or the
#' BHHH matrix for `"bhhh"`), It = [I D'] I [I; D] (the curvature along the
#' constraint) and R = (I_cd + D' I_dd) C_d^{-1},
#'   theta_c - theta0 = It^{-1} sum_u w_u (m_u + R c_u)
#'   delta - delta0   = D (theta_c - theta0) - C_d^{-1} sum_u w_u c_u,
#' so that, with Omega the variance of the stacked sums (blocks 11: m, cc: c,
#' 12: their covariance),
#'   V_c  = It^{-1} (Omega_11 + R Omega_21 + Omega_12 R' + R Omega_cc R') It^{-1}
#'   X    = It^{-1} (Omega_12 + R Omega_cc) C_d^{-1}
#'   V_cd = V_c D' - X
#'   V_d  = D V_c D' - D X - X' D' + C_d^{-1} Omega_cc C_d^{-1}.
#' D's own sampling error does not enter to first order: it multiplies the
#' constants' score, whose mean is zero at the truth. Omega by type:
#' * `"hessian"` and `"bhhh"`, under the model with the weights as
#'   frequencies: Omega_11 = It (the information equality), Omega_12 = 0, and
#'   Omega_cc = C_d + K for the sample's own shares (K = sum_u w_u
#'   Cov_s(Q_us), the share Jacobian kernel's variance mode: the model's
#'   variance of each unit's chosen counts) or sum_u w_u c_u c_u' for
#'   supplied targets (the sample's covariate composition; the meat kernel
#'   with the weights);
#' * `"sandwich"`: sum_u w_u^2 v_u v_u', v_u = (m_u, c_u) (the meat kernel
#'   with the squared weights);
#' * `"cluster"`: sum_g v_g v_g', v_g = sum_{u in g} w_u v_u.
#' Supplied targets estimated from an independent sample add
#' diag(W) V_S diag(W) to Omega_cc, W the weight of each alternative's
#' market. A supplied target's residual is dense, -T_u s*_m on the market's
#' alternatives the unit never faced; the meat kernel forms the part on the
#' unit's own alternatives, and `.blp_meat_correct()` adds the rest.
#'
#' Every type assumes correctly specified choice probabilities, as the
#' plan's estimator does: `"sandwich"` and `"cluster"` relax the information
#' equality (sampling weights such as WESML's, dependence within clusters),
#' not the model. Under misspecification the constants' expected score is not
#' zero, and two terms this omits become first order: the curvature of D in
#' the bread (sum_j E[s_d,j] d^2 delta_j / d theta_c d theta_c'), and D's own
#' sampling error in the meat.
#'
#' @param info The information at the estimate, n x n (the BHHH matrix for
#'   `"bhhh"`), theta_c first.
#' @param jac `mxl_blp_jacobian_parallel()` at the estimate; with `K` (its
#'   variance mode) for the sample's shares under `"hessian"` and `"bhhh"`.
#' @param type One of `"hessian"`, `"bhhh"`, `"sandwich"`, `"cluster"`.
#' @param tg The targets: `source` (`"sample"` or `"supplied"`); for supplied
#'   targets, `s` (each free constant's target, a share within its market),
#'   `market` (each free constant's market, 1..M), `b` (per market, from
#'   `.blp_kappa()`), and, for targets with sampling error, `vcov` (their
#'   covariance over the free constants) and `W` (each free constant's
#'   market weight).
#' @param meat `mxl_blp_meat_parallel()`'s result: for `"sandwich"` and
#'   `"cluster"`, and for supplied targets under `"hessian"` and `"bhhh"`.
#' @param n_c The number of continuous parameters (theta_c).
#' @param D The implicit derivative the meat was formed with, or NULL to
#'   form it here.
#' @returns List with `vcov` (n x n, or NULL when a factor is singular or
#'   not positive definite) and `se`, as `invert_hessian()` returns them, and
#'   `It`, the curvature along the constraint (for the outer certificate;
#'   NULL when the solve with C_d failed before it was formed).
#' @noRd
.blp_vcov <- function(info, jac, type, tg, meat = NULL, n_c, D = NULL) {
  stopifnot(length(type) == 1L, type %in% c("hessian", "bhhh", "sandwich", "cluster"),
            length(tg[["source"]]) == 1L, tg[["source"]] %in% c("sample", "supplied"))
  n <- nrow(info)
  ci <- seq_len(n_c)
  di <- n_c + seq_len(n - n_c)
  fail <- function(why, It = NULL) {
    message("Standard errors set to NA: ", why)
    list(vcov = NULL, se = rep(NA_real_, n), It = It)
  }
  # Finite inputs, so that every solve below has finite right-hand sides
  # (spd_solve() returns zeros for one that is not)
  parts <- list(info, jac$jacobian, if (identical(tg[["source"]], "sample")) jac$K,
                meat$meat, meat$A_c, meat$A_d, tg[["s"]], tg[["b"]], tg[["vcov"]],
                tg[["W"]], D)
  if (!all(vapply(parts, function(x) all(is.finite(x)), TRUE))) {
    return(fail("the information, the share Jacobian, the meat or the targets are not finite."))
  }
  # Every solve with C_d reads it in place in the Jacobian (its trailing
  # block), one factorization each
  sc <- .jacobi_scale_of(diag(jac$jacobian)[di])
  cd_solve <- function(B) spd_solve(jac$jacobian, B, sc, n_c)
  s1 <- cd_solve(jac$jacobian[di, ci, drop = FALSE])
  if (s1$info != 0L) {
    return(fail(paste0("the share Jacobian on the constants is not positive ",
                       "definite, or is ill-conditioned.")))
  }
  if (is.null(D)) D <- -s1$x
  rm(s1)
  I_dc <- info[di, ci, drop = FALSE]
  I_dd_D <- info[di, di, drop = FALSE] %*% D
  DI <- crossprod(D, I_dc)
  It <- info[ci, ci, drop = FALSE] + DI + t(DI) + crossprod(D, I_dd_D)
  It <- (It + t(It)) / 2
  Rt <- cd_solve(I_dc + I_dd_D)$x
  rm(I_dd_D, DI, I_dc)

  # Omega's blocks
  if (type %in% c("hessian", "bhhh")) {
    O11 <- It
    O12 <- matrix(0, n_c, n - n_c)
    Occ <- if (identical(tg[["source"]], "sample")) {
      jac$jacobian[di, di, drop = FALSE] + jac$K
    } else {
      .blp_meat_correct(meat, tg, n_c)[di, di, drop = FALSE]
    }
  } else {
    O <- .blp_meat_correct(meat, tg, n_c)
    O11 <- O[ci, ci, drop = FALSE]
    O12 <- O[ci, di, drop = FALSE]
    Occ <- O[di, di, drop = FALSE]
    rm(O)
  }
  if (!is.null(tg[["vcov"]])) Occ <- Occ + tg[["vcov"]] * tcrossprod(tg[["W"]])

  Iti <- tryCatch(equilibrated_inverse(It, .jacobi_scale(It)), error = function(e) {
    if (!grepl("singular", conditionMessage(e), fixed = TRUE)) stop(e)
    NULL
  })
  if (is.null(Iti)) {
    return(fail("the information along the share constraint is singular.", It))
  }
  RO <- crossprod(Rt, Occ)                      # R Omega_cc
  M11 <- O11 + crossprod(Rt, t(O12)) + O12 %*% Rt + RO %*% Rt
  Vc <- Iti %*% M11 %*% Iti
  Vc <- (Vc + t(Vc)) / 2
  # X' = C_d^{-1} (Omega_12 + R Omega_cc)' It^{-1}
  Xt <- cd_solve(t(O12 + RO))$x %*% Iti
  rm(RO, O12, M11, Rt)
  # C_d^{-1} Omega_cc C_d^{-1}, then D V_c D' - D X - X' D', the last two
  # as -2 D X before the symmetrization
  Z <- cd_solve(Occ)$x
  rm(Occ)
  Zt <- t(Z)
  rm(Z)
  Vd <- cd_solve(Zt)$x
  rm(Zt)
  Vcd <- tcrossprod(Vc, D) - t(Xt)              # V_c D' - X
  Vd <- Vd + D %*% (Vcd - t(Xt))
  Vd <- (Vd + t(Vd)) / 2
  V <- matrix(0, n, n)
  V[ci, ci] <- Vc
  V[ci, di] <- Vcd
  V[di, ci] <- t(Vcd)
  V[di, di] <- Vd
  rm(Vd)

  dv <- diag(V)
  neg <- !is.na(dv) & dv < 0
  if (any(neg)) {
    message("The variance is not positive definite; ", sum(neg),
            " variance(s) negative. Standard errors set to NA for those ",
            "parameters.")
    dv[neg] <- NA_real_
  }
  list(vcov = V, se = sqrt(dv), It = It)
}

#' A meat with the supplied targets' dense residuals
#'
#' The meat kernel's sum over v^s_u, v_u's part on the unit's own
#' alternatives, plus the rest of v_u = v^s_u - T_u e_m(u), e_m = (0,
#' s*_m on market m's free constants): with A_m the kernel's per-market sums
#' sum_u alpha_u v^s_u and b_m `.blp_kappa()`'s,
#'   Omega = meat - sum_m (A_m e_m' + e_m A_m') + sum_m b_m e_m e_m'.
#' A unit lies within one market and markets share no alternatives, so the
#' constants' part of every A_m is one vector (`A_d`).
#'
#' @param meat `mxl_blp_meat_parallel()`'s result: `meat`, `A_c` (n_c x M)
#'   and `A_d` (one per free constant).
#' @param tg The targets, as `.blp_vcov()` takes them.
#' @param n_c The number of continuous parameters.
#' @returns The n x n meat.
#' @noRd
.blp_meat_correct <- function(meat, tg, n_c) {
  O <- meat$meat
  if (!identical(tg[["source"]], "supplied")) return(O)
  ci <- seq_len(n_c)
  di <- n_c + seq_along(tg[["s"]])
  s <- tg[["s"]]
  mk <- tg[["market"]]
  oc <- -meat$A_c[, mk, drop = FALSE] * rep(s, each = n_c)
  O[ci, di] <- O[ci, di] + oc
  O[di, ci] <- O[di, ci] + t(oc)
  for (m in unique(mk)) {
    j <- which(mk == m)
    sj <- s[j]
    aj <- meat$A_d[j]
    O[di[j], di[j]] <- O[di[j], di[j]] - (aj %o% sj + sj %o% aj) +
      tg[["b"]][m] * (sj %o% sj)
  }
  O
}

#' The per-unit and per-market sums for supplied targets' dense residuals
#'
#' `mxl_blp_meat_parallel()` sums alpha_u v^s_u per market, and
#' `.blp_meat_correct()` scales e_m e_m' by b_m. Robust meats (sum_u w_u^p
#' v_u v_u', p = 1 for Omega_cc under `"hessian"` and `"bhhh"`, 2 for the
#' sandwich): alpha_u = w_u^p T_u, b_m = sum_{u in m} w_u^p T_u^2. Cluster
#' meats (clusters within markets): kappa_g = sum_{u in g} w_u T_u, alpha_u =
#' kappa_g(u) w_u, b_m = sum_{g in m} kappa_g^2.
#'
#' @param w Each unit's weight.
#' @param T_u Each unit's number of choice situations.
#' @param market Each unit's market, 1..M.
#' @param n_markets M.
#' @param cluster Each unit's cluster code, or NULL for a robust meat.
#' @param power p, for a robust meat.
#' @returns List with `alpha` (per unit) and `b` (per market).
#' @noRd
.blp_kappa <- function(w, T_u, market, n_markets, cluster = NULL, power = 2) {
  per_market <- function(x, m) {
    as.numeric(tapply(x, factor(m, seq_len(n_markets)), sum, default = 0))
  }
  if (is.null(cluster)) {
    wp <- w^power
    return(list(alpha = wp * T_u, b = per_market(wp * T_u^2, market)))
  }
  kappa <- rowsum(w * T_u, cluster, reorder = FALSE)
  kg <- kappa[match(cluster, rownames(kappa))]
  g_first <- !duplicated(cluster)
  list(alpha = kg * w, b = per_market(kg[g_first]^2, market[g_first]))
}

#' The share-matching constants' requirements, checked before the draws
#'
#' `concentrate_out_asc = "blp"` solves the share equations for the constants
#' at every value of the other parameters, without bounds, and factors their
#' share Jacobian densely: it needs constants, none of them bounded (or
#' pinned), no polish (a Newton phase on the likelihood would move the
#' estimate to the maximum-likelihood estimate), and at most `max_dense` of
#' them. With the sample's own shares as targets, every alternative must be
#' chosen somewhere and not wherever it is offered (`.check_profile_inputs()`
#' reads the same cases): such a target is no share the model attains.
#'
#' @param n_asc Number of free constants.
#' @param asc_lower,asc_upper Their bounds.
#' @param asc_names Their names.
#' @param alt_mapping The preparation's `alt_mapping`.
#' @param ctl From `.concentrate_control()`.
#' @param d The prepared `input_data`.
#' @param polish The fit's `polish`.
#' @param sample_targets Whether the targets are the sample's own shares.
#' @returns NULL, invisibly.
#' @noRd
.check_blp_inputs <- function(n_asc, asc_lower, asc_upper, asc_names,
                              alt_mapping, ctl, d, polish, sample_targets) {
  what <- "concentrate_out_asc = \"blp\""
  if (n_asc == 0L) {
    stop(what, " solves the share equations for the alternative-specific ",
         "constants, so it needs them: use_asc = TRUE and more than one ",
         "alternative.", call. = FALSE)
  }
  b <- is.finite(asc_lower) | is.finite(asc_upper)
  if (any(b)) {
    stop(what, " solves the share equations for the constants without ",
         "bounds, but `lower`/`upper` bound ",
         paste(utils::head(asc_names[b], 5L), collapse = ", "),
         if (sum(b) > 5L) paste0(" and ", sum(b) - 5L, " more"),
         ". Remove those bounds, or fit without concentrate_out_asc.",
         call. = FALSE)
  }
  if (!identical(polish, "none")) {
    stop(what, " estimates the constants from the share equations, and ",
         "polish = \"", polish, "\" would move them, with the other ",
         "parameters, to the maximum-likelihood estimate. Fit with ",
         "polish = \"none\".", call. = FALSE)
  }
  if (sample_targets) {
    am <- alt_mapping
    if (all(c("alt_int", "N_OBS", "N_CHOICES") %in% names(am))) {
      am <- as.list(am)
      code <- am$alt_int
      inside <- code >= 1L & !is.na(code)
      un <- .unchosen_offers(d, max(c(0L, code[inside])))
      am$N_OBS[inside] <- am$N_CHOICES[inside] + un[code[inside]]
    }
    u <- .unbounded_constants(am)
    if (!is.null(u)) {
      msg <- .unbounded_constants_text(u, alt_mapping)
      stop(errorCondition(paste(c(paste0(
        what, " with the sample's own shares as targets needs every ",
        "alternative chosen in some choice situation, and not in every ",
        "situation that offers it: no constants attain such a share. Supply ",
        "target_shares, or:"), msg[-1L]), collapse = "\n"),
        class = "choicer_unbounded_constants", call = NULL))
    }
  }
  if (n_asc > ctl$max_dense) {
    stop(what, " factors the constants' share Jacobian densely, up to ",
         "concentrate_control$max_dense = ", ctl$max_dense, " constants; this ",
         "model has ", n_asc, ". Each evaluation holds two J x J matrices ",
         "(8 J^2 bytes each) and factors one (about J^3 / 3 operations): ",
         "raise max_dense, or fit without concentrate_out_asc.", call. = FALSE)
  }
  invisible(NULL)
}

#' The choice situations that offer one alternative alone
#'
#' Whether every row of each situation carries the same alternative (the
#' share inversion's captive situations, without an outside option, however
#' often a situation lists it), in blocks of about 2^22 rows.
#'
#' @param d A prepared `input_data`.
#' @returns Logical vector, one per situation.
#' @noRd
.situations_alone <- function(d) {
  M <- d$M
  n <- length(M)
  out <- logical(n)
  if (!n) return(out)
  off <- c(0, cumsum(as.numeric(M)))
  step <- max(1, floor(2^22 / max(1, off[n + 1L] / n)))
  lo <- 1
  while (lo <= n) {
    hi <- min(n, lo + step - 1)
    code <- d$alt_idx[seq.int(off[lo] + 1, off[hi + 1L])]
    other <- code != rep.int(d$alt_idx[off[lo:hi] + 1], M[lo:hi])
    out[lo:hi] <- tabulate(rep.int(seq_len(hi - lo + 1), M[lo:hi])[other],
                           hi - lo + 1) == 0L
    lo <- hi + 1
  }
  out
}

#' The sample's own shares as targets: the checks market by market
#'
#' `.check_blp_inputs()` checks each alternative over the whole sample
#' (chosen somewhere, and not wherever it is offered). Here, on the exact
#' chosen counts, what no constants attain within a market: with an outside
#' option, a market where no one chooses it (its outside share 0); without
#' one, an alternative chosen only in the situations that offer it alone
#' (its share their weight, which it has at any constant).
#'
#' @param chosen Each situation's chosen inside alternative, 0 for the
#'   outside option.
#' @param sit_mk Each situation's market.
#' @param d The prepared `input_data`.
#' @param outside Whether there is an outside option.
#' @param n_mk,J The markets and the inside alternatives.
#' @param market Each inside alternative's market.
#' @param inside_lab Their labels.
#' @returns NULL, invisibly; an error of class `choicer_unbounded_constants`.
#' @noRd
.check_blp_sample_targets <- function(chosen, sit_mk, d, outside, n_mk, J,
                                      market, inside_lab) {
  fmt <- function(x) {
    paste0(paste(utils::head(x, 5L), collapse = ", "),
           if (length(x) > 5L) paste0(" and ", length(x) - 5L, " more"))
  }
  what <- paste0("concentrate_out_asc = \"blp\" with the sample's own shares ",
                 "as targets")
  if (outside) {
    none <- which(tabulate(sit_mk[chosen == 0L], n_mk) == 0L)
    if (length(none)) {
      stop(errorCondition(paste0(
        what, " needs the outside option chosen in every market (a group of ",
        "alternatives that share choice situations): no constants attain an ",
        "outside share of 0, and no one chooses it in the market of ",
        fmt(inside_lab[market == none[1L]]),
        if (length(none) > 1L) {
          paste0(" (and in ", length(none) - 1L, " other market",
                 if (length(none) > 2L) "s", ")")
        },
        ". Supply target_shares, or drop that market's choice situations."),
        class = "choicer_unbounded_constants", call = NULL))
    }
  } else {
    alone <- which(tabulate(chosen[chosen > 0L & !.situations_alone(d)], J) == 0L)
    if (length(alone)) {
      stop(errorCondition(paste0(
        what, " needs every alternative chosen in some choice situation that ",
        "offers another: ", fmt(inside_lab[alone]), " ",
        if (length(alone) > 1L) "are" else "is",
        " chosen only where offered alone, a share no constant attains (such ",
        "situations say nothing of a constant). Supply target_shares, or drop ",
        "or pool these alternatives."),
        class = "choicer_unbounded_constants", call = NULL))
    }
  }
  invisible(NULL)
}

#' Is a target covariance positive semi-definite, to rounding?
#'
#' A row of zero variance must be zero, and the others, on the correlation
#' scale, must leave a positive definite matrix when 1e-6 is added to the
#' diagonal (every eigenvalue above -1e-6): a covariance given to a few
#' digits passes, a matrix that is not a covariance fails. One Cholesky
#' factorization. The caller passes the block the variance takes, the free
#' constants': without an outside option, shares within a market from n
#' choices have a covariance singular by construction, but not on that
#' block.
#'
#' @param V A symmetric matrix.
#' @returns TRUE or FALSE.
#' @noRd
.blp_vcov_psd <- function(V) {
  d <- diag(V)
  if (any(d < 0)) return(FALSE)
  z <- d == 0
  if (any(V[z, ] != 0)) return(FALSE)
  if (all(z)) return(TRUE)
  s <- 1 / sqrt(d[!z])
  R <- V[!z, !z, drop = FALSE] * (s %o% s)
  diag(R) <- 1 + 1e-6
  !inherits(tryCatch(chol(R), error = identity), "error")
}

#' The share inversion's checks of the targets, in the fit's terms
#'
#' `blp_check_targets()` with the targets (the offer and captive bounds, each
#' market's outside share), its error restated: which targets they are, the
#' labels of the inside alternatives its message numbers, and, over several
#' markets, that its shares are of the whole sample's weight.
#'
#' @param alt_idx,M,w,J,outside As `blp_check_targets()` takes them.
#' @param engine The targets in the share inversion's convention.
#' @param source `"sample"`, `"supplied"` or `"wesml"`.
#' @param inside_lab The inside alternatives' labels, by code.
#' @param n_mk The markets.
#' @returns NULL, invisibly.
#' @noRd
.blp_engine_check <- function(alt_idx, M, w, J, outside, engine, source,
                              inside_lab, n_mk) {
  r <- tryCatch(blp_check_targets(alt_idx, M, w, J, outside, target_shares = engine),
                error = identity)
  if (!inherits(r, "error")) return(invisible(NULL))
  msg <- conditionMessage(r)
  # the codes after "alternative(s) " (no lookbehind: PCRE2 before 10.43
  # takes none of variable length); the restatement never masks the message
  codes <- tryCatch({
    hit <- regmatches(msg, gregexpr("lternatives? [0-9]+(, [0-9]+)*", msg))[[1L]]
    x <- unique(as.integer(unlist(strsplit(sub("^lternatives? ", "", hit), ", ",
                                           fixed = TRUE))))
    x[!is.na(x) & x >= 1L & x <= J]
  }, error = function(e) integer())
  stop("concentrate_out_asc = \"blp\": ",
       switch(source, sample = "the sample's own shares",
              wesml = "the choice-based sample's population shares Q",
              "target_shares"),
       ", as targets, fail the share inversion's checks. ", msg,
       if (length(codes)) {
         paste0(" (In its numbering, ",
                paste0(codes, " is ", encodeString(inside_lab[codes], quote = '"'),
                       collapse = ", "), ".)")
       },
       if (n_mk > 1L) {
         paste0(" Its shares are of the whole sample's weight: a target within ",
                "a market times the market's share of the weight.")
       }, call. = FALSE)
}

#' The share-matching constants' targets
#'
#' The targets in the share inversion's convention (shares of the total
#' weight, the outside option's first), from `target_shares` (shares within
#' each market: a market is a group of alternatives that share choice
#' situations, `blp_check_targets()`'s components), a choice-based sample's
#' population shares Q (its `choice_sampling` provenance), or the sample's
#' own weighted shares; checked by the share inversion's rules, with
#' `target_vcov`, before the draws. A market's within-market targets s*
#' become the total weight's W_m / W s*, so that its share equations are
#' sum_{t in m} w_t Pbar_tj = W_m s*_j.
#'
#' @param target_shares,target_vcov The user's.
#' @param d The prepared `input_data`.
#' @param n_asc Number of free constants.
#' @param choice_sampling The resolved provenance (`.resolve_choice_sampling()`).
#' @returns List: `source` (`"sample"`, `"supplied"` or `"wesml"`), `engine`
#'   (the share inversion's targets), `counts` (each free constant's target
#'   count, W times its share of the total), `shares` (each inside
#'   alternative's target within its market, named), `market` (each inside
#'   alternative's market), `market_weight` (each market's weight), `vcov`
#'   (`target_vcov` in the alternatives' order, or NULL) and `unit_market`
#'   (each likelihood unit's market).
#' @noRd
.blp_resolve_targets <- function(target_shares, target_vcov, d, n_asc,
                                 choice_sampling) {
  what <- "concentrate_out_asc = \"blp\""
  outside <- isTRUE(d$include_outside_option)
  J <- if (outside) n_asc else n_asc + 1L
  w <- as.numeric(d$weights)
  mk <- blp_check_targets(d$alt_idx, d$M, w, J, outside)
  n_mk <- max(c(1L, mk$market), na.rm = TRUE)
  if (!outside && n_mk > 1L) {
    stop(what, " needs, without an outside option, alternatives that all ",
         "share choice situations, directly or through others: these data ",
         "split them into ", n_mk, " groups that never meet, whose constants ",
         "have no common level. Fit with an outside option, or each group ",
         "on its own.", call. = FALSE)
  }
  am <- d$alt_mapping
  code <- am[["alt_int"]]
  labels <- as.character(am[[2L]])
  inside_lab <- labels[match(seq_len(J), code)]
  if (anyNA(mk$market)) {
    stop(what, ": ", paste(utils::head(inside_lab[is.na(mk$market)], 5L),
                           collapse = ", "),
         " offered in no choice situation of positive weight, so no target ",
         "share can be met there and the constant is not identified.",
         call. = FALSE)
  }
  W_tot <- mk$weight
  first <- c(0, cumsum(as.numeric(d$M)))[seq_along(d$M)]
  sit_mk <- mk$market[d$alt_idx[first + 1]]   # each situation's market
  wesml <- is.list(choice_sampling) && identical(choice_sampling$scheme, "wesml") &&
    !is.null(choice_sampling$Q)
  source <- if (!is.null(target_shares)) "supplied" else if (wesml) "wesml" else "sample"
  if (!is.null(target_vcov) && !identical(source, "supplied")) {
    stop("`target_vcov` is the sampling covariance of supplied `target_shares`; ",
         "without them it has nothing to describe.", call. = FALSE)
  }
  perm <- NULL   # named target_shares: their entries in the alternatives' order
  S_in <- switch(source,
    sample = {
      ch <- d$choice_idx > 0L
      chosen <- integer(length(d$M))   # 0: the outside option
      chosen[ch] <- d$alt_idx[first[ch] + d$choice_idx[ch]]
      .check_blp_sample_targets(chosen, sit_mk, d, outside, n_mk, J,
                                mk$market, inside_lab)
      cnt <- numeric(J)
      if (any(ch)) {
        s <- rowsum(w[ch], chosen[ch])
        cnt[as.integer(rownames(s))] <- s[, 1L]
      }
      cnt / W_tot
    },
    wesml = {
      Q <- choice_sampling$Q
      miss <- setdiff(inside_lab, names(Q))
      if (length(miss)) {
        stop(what, " takes a choice-based sample's population shares Q as ",
             "its targets, but Q has none for ",
             paste(utils::head(miss, 5L), collapse = ", "),
             ". Supply target_shares.", call. = FALSE)
      }
      as.numeric(Q[inside_lab])
    },
    supplied = {
      ts <- .int64_to_double(target_shares, "`target_shares`")
      if (!is.numeric(ts) || !is.null(dim(ts)) && length(dim(ts)) > 1L &&
          ncol(ts) != 1L) {
        stop("`target_shares` must be a numeric vector of shares.", call. = FALSE)
      }
      nm <- names(ts)
      ts <- as.numeric(ts)
      bad <- which(!is.finite(ts) | !(ts > 0))
      if (length(bad)) {
        stop(sprintf("`target_shares` must be finite and positive (entry %d%s is %s).",
                     bad[1L],
                     if (!is.null(nm)) paste0(", ", encodeString(nm[bad[1L]], quote = '"'), ",")
                     else "",
                     format(ts[bad[1L]])), call. = FALSE)
      }
      if (!is.null(nm)) {
        if (length(ts) != J || anyNA(nm) || anyDuplicated(nm) ||
            !setequal(nm, inside_lab)) {
          stop("`target_shares` names must be the inside alternatives' labels, ",
               "each once (", J, if (outside) {
                 "; the outside option's share follows from theirs"
               }, "): ", paste(utils::head(inside_lab, 5L), collapse = ", "),
               if (J > 5L) ", ...", ".", call. = FALSE)
        }
        perm <- match(inside_lab, nm)
        ts <- ts[perm]
      } else if (outside && length(ts) == J + 1L) {
        # blp()'s convention, one market: the outside option's first
        if (n_mk > 1L) {
          stop("`target_shares` with the outside option's share first is ",
               "taken for one market only; these data have ", n_mk, ". Give ",
               "each inside alternative's share within its market.",
               call. = FALSE)
        }
        if (!(abs(sum(ts) - 1) <= sqrt(.Machine$double.eps))) {
          stop(sprintf(paste0("`target_shares` with the outside option's share ",
                              "first must sum to 1; they sum to %.15g."), sum(ts)),
               call. = FALSE)
        }
        ts <- ts[-1L] / sum(ts)
      } else if (length(ts) != J) {
        stop(sprintf(paste0("`target_shares` must have one share per inside ",
                            "alternative (%d)%s; it has %d."), J,
                     if (outside) ", or one more, the outside option's first" else "",
                     length(ts)), call. = FALSE)
      }
      mk$mass[mk$market] * ts
    })
  engine <- if (outside) c(1 - sum(S_in), S_in) else S_in
  .blp_engine_check(d$alt_idx, d$M, w, J, outside, engine, source, inside_lab,
                    n_mk)
  free <- if (outside) seq_len(J) else seq_len(J)[-1L]
  if (!outside) {
    # one market, whose targets the share inversion rescales to sum to 1
    # (within its tolerance): the counts and the variance take them so
    S_in <- S_in / sum(S_in)
    engine <- S_in
  }
  if (!is.null(target_vcov)) {
    V <- .int64_to_double(target_vcov, "`target_vcov`")
    if (!is.numeric(V) || !is.matrix(V) || !all(dim(V) == c(J, J)) ||
        !all(is.finite(V)) || max(abs(V - t(V))) > 1e-8 * max(1e-300, abs(V))) {
      stop(sprintf(paste0("`target_vcov` must be a finite symmetric %d x %d ",
                          "matrix: the targets' sampling covariance, over the ",
                          "inside alternatives in target_shares' order."), J, J),
           call. = FALSE)
    }
    # into the alternatives' order: by its dimnames, when it has them (the
    # labels), else as target_shares are
    dn <- dimnames(V)
    if (!is.null(dn[[1L]]) || !is.null(dn[[2L]])) {
      if (!identical(dn[[1L]], dn[[2L]]) || anyNA(dn[[1L]]) ||
          anyDuplicated(dn[[1L]]) || !setequal(dn[[1L]], inside_lab)) {
        stop("`target_vcov`'s row and column names, when it has them, must ",
             "both be the inside alternatives' labels, in one order.",
             call. = FALSE)
      }
      perm <- match(inside_lab, dn[[1L]])
    }
    if (!is.null(perm)) V <- V[perm, perm, drop = FALSE]
    V <- unname((V + t(V)) / 2)
    if (!.blp_vcov_psd(V[free, free, drop = FALSE])) {
      stop("`target_vcov` must be a covariance, positive semi-definite to ",
           "within rounding (checked on ",
           if (outside) "the inside alternatives" else {
             "the alternatives other than the base, whose row and column the variance drops"
           },
           ", on the correlation scale).", call. = FALSE)
    }
    target_vcov <- V
  }
  # each likelihood unit's market: its situations' (one market per unit for
  # targets other than the sample's, which the variance takes)
  first_u <- .unit_first(d)
  unit_market <- sit_mk[first_u]
  if (!identical(source, "sample") && n_mk > 1L && !is.null(d$Ti) &&
      any(sit_mk != rep(unit_market, d$Ti))) {
    stop(what, " with targets other than the sample's own shares needs each ",
         "decision maker's choice situations in one market (a group of ",
         "alternatives that share situations); these data have a decision ",
         "maker in two.", call. = FALSE)
  }
  list(source = source, engine = engine,
       counts = W_tot * S_in[free],
       shares = stats::setNames(S_in / mk$mass[mk$market], inside_lab),
       market = mk$market, market_weight = W_tot * mk$mass,
       vcov = target_vcov, unit_market = unit_market)
}

#' The inner solve of the share-matching constants
#'
#' At fixed continuous parameters theta_c: the constants delta(theta_c) that
#' solve the share equations on the estimation draws (the share inversion,
#' `engine_f`), started from `delta0` advanced along the best point's
#' implicit derivative (delta_best + D_best (theta_c - theta_c,best)), then
#' one pass of the share Jacobian kernel there (`jac_f`). With its
#' objective f, gradient g and Jacobian C, and lambda = C_d^{-1} g_d (one
#' factorization of C_d, which also gives D = -C_d^{-1} C_c):
#' * the value is f - lambda' (counts - target counts), the objective at the
#'   exact solution to second order in the inner solve's error (the
#'   envelope theorem does not hold here: the constants' score is not zero);
#' * the gradient's theta_c entries are g_c - C_c' lambda, the total
#'   derivative along delta(theta_c) (the implicit-function term).
#' A solve fails (`converged` FALSE, the overflow sentinel to the outer
#' optimizer) when the share inversion does not converge, the likelihood
#' overflows, or C_d is not positive definite.
#'
#' @param engine_f Function(theta_c, delta0, max_iter) returning the share
#'   inversion's record with `delta_free`.
#' @param jac_f Function(theta) returning `mxl_blp_jacobian_parallel()`'s.
#' @param n_c,n The continuous parameters (first) and all.
#' @param counts The free constants' target counts.
#' @param free The free constants' inside alternatives.
#' @returns Function(theta_c, delta0, max_iter, best = NULL) returning the
#'   profile solver's record (`.mxl_profile_solver()`), with `extra` (`D`,
#'   `lambda`, and `objective`, the negated log-likelihood at the solved
#'   constants, uncorrected).
#' @noRd
.mxl_blp_solver <- function(engine_f, jac_f, n_c, n, counts, free) {
  force(engine_f)
  force(jac_f)
  ci <- seq_len(n_c)
  di <- n_c + seq_len(n - n_c)
  function(theta_c, delta0, max_iter, best = NULL) {
    if (!is.null(best$extra$D)) {
      delta0 <- delta0 + as.numeric(best$extra$D %*% (theta_c - best$theta_c))
    }
    out <- list(converged = FALSE, value = NA_real_, delta = delta0,
                gradient = NULL, iterations = 0L, passes = 0L,
                residual = NA_real_, decrement = NA_real_, message = NULL,
                error = NULL)
    r <- tryCatch(engine_f(theta_c, delta0, max_iter), error = identity)
    if (inherits(r, "error")) {
      out$error <- conditionMessage(r)
      return(out)
    }
    out$iterations <- r$evaluations
    out$passes <- r$evaluations
    out$delta <- r$delta_free
    out$residual <- r$residual
    if (!identical(r$code, 0L)) {
      out$message <- paste0("the share inversion stopped (", .BLP_STATUS[r$code + 1L],
                            ") with a log-share residual of ",
                            format(r$residual, digits = 3))
      return(out)
    }
    th <- numeric(n)
    th[ci] <- theta_c
    th[di] <- r$delta_free
    j <- tryCatch(jac_f(th), error = identity)
    if (inherits(j, "error")) {
      out$error <- conditionMessage(j)
      return(out)
    }
    out$passes <- out$passes + 1L
    if (isTRUE(j$overflow) || isTRUE(j$nonfinite)) {
      out$message <- "the likelihood overflows at the solved constants"
      return(out)
    }
    Cc <- j$jacobian[di, ci, drop = FALSE]
    sol <- spd_solve(j$jacobian, cbind(j$gradient[di], Cc),
                     .jacobi_scale_of(diag(j$jacobian)[di]), n_c)
    if (sol$info != 0L) {
      out$message <- "the share Jacobian on the constants is not positive definite"
      return(out)
    }
    lambda <- sol$x[, 1L]
    g <- j$gradient
    g[ci] <- g[ci] - as.numeric(crossprod(Cc, lambda))
    out$value <- j$objective - sum(lambda * (as.numeric(j$counts)[free] - counts))
    out$gradient <- g
    out$converged <- TRUE
    out$extra <- list(D = -sol$x[, -1L, drop = FALSE], lambda = lambda,
                      objective = j$objective)
    out
  }
}

#' The variance of a share-matching fit, from its stored data
#'
#' `.blp_vcov()` with the kernels at the estimate on the fit's draws: the
#' information (the Hessian, or the BHHH matrix for `"bhhh"`), the share
#' Jacobian (with K for the sample's own shares under `"hessian"` and
#' `"bhhh"`), and the meat where Omega needs it (`"robust"` and `"cluster"`,
#' and supplied or WESML targets under `"hessian"` and `"bhhh"`).
#'
#' @param object A fitted `choicer_mxl` with `share_targets` and its data.
#' @param type One of `"hessian"`, `"bhhh"`, `"robust"`, `"cluster"`.
#' @param cluster,what As `.assemble_score_vcov()` takes them.
#' @param gp Draw arguments from `.mxl_gen_params()`, or NULL.
#' @returns `.blp_vcov()`'s list.
#' @noRd
.blp_fit_vcov <- function(object, type, cluster = NULL, gp = NULL,
                          what = "`cluster`") {
  d <- object[["data"]]
  tg <- object$share_targets
  theta <- object$coefficients
  n <- length(theta)
  J_d <- length(object$param_map$asc)
  n_c <- n - J_d
  di <- n_c + seq_len(J_d)
  first_u <- .unit_first(d)
  U <- length(first_u)
  sample <- identical(tg[["source"]], "sample")
  cl <- NULL
  if (identical(type, "cluster")) {
    cluster <- .unit_clusters(object, cluster, what)
    .check_cluster_labels(cluster, U)
    cl <- .cluster_codes(cluster)
    if (!sample && length(unique(tg$unit_market)) > 1L &&
        any(tapply(tg$unit_market, cl, function(m) length(unique(m))) > 1L)) {
      stop("Cluster-robust standard errors of concentrate_out_asc = \"blp\" with ",
           "targets other than the sample's own shares need clusters within ",
           "markets (groups of alternatives that share choice situations); a ",
           "cluster here spans two.", call. = FALSE)
    }
  }
  if (is.null(gp)) gp <- .mxl_gen_params(object$draws_info, N = U)
  args <- list(theta = theta, X = d$X, W = d$W, alt_idx = d$alt_idx,
               choice_idx = d$choice_idx, M = d$M, weights = d$weights,
               eta_draws = gp$eta_draws, rc_dist = object$rc_dist,
               rc_correlation = object$rc_correlation, rc_mean = object$rc_mean,
               use_asc = TRUE,
               include_outside_option = isTRUE(object$include_outside_option),
               gen_seed = gp$gen_seed, gen_scramble = gp$gen_scramble,
               gen_S = gp$gen_S, Ti = d$Ti)
  model_based <- type %in% c("hessian", "bhhh")
  jac <- do.call(mxl_blp_jacobian_parallel,
                 c(args, list(variance = sample && model_based)))
  info <- do.call(if (identical(type, "bhhh")) mxl_bhhh_parallel else mxl_hessian_parallel,
                  args)
  outside <- isTRUE(object$include_outside_option)
  free <- if (outside) seq_len(J_d) else seq_len(J_d + 1L)[-1L]
  tv <- list(source = if (sample) "sample" else "supplied")
  meat <- NULL
  if (!(sample && model_based)) {
    s1 <- spd_solve(jac$jacobian, jac$jacobian[di, seq_len(n_c), drop = FALSE],
                    .jacobi_scale_of(diag(jac$jacobian)[di]), n_c)
    if (s1$info != 0L) {
      message("Standard errors set to NA: the share Jacobian on the constants ",
              "is not positive definite.")
      return(list(vcov = NULL, se = rep(NA_real_, n)))
    }
    ka <- NULL
    if (!sample) {
      T_u <- diff(c(first_u, length(d$M) + 1L))
      ka <- .blp_kappa(as.numeric(d$weights)[first_u], T_u, tg$unit_market,
                       length(tg$market_weight), cluster = cl,
                       power = if (identical(type, "robust")) 2 else 1)
      mf <- tg[["market"]][free]
      tv <- list(source = "supplied", s = as.numeric(tg$shares)[free], market = mf,
                 b = ka$b,
                 vcov = if (!is.null(tg[["vcov"]])) tg[["vcov"]][free, free, drop = FALSE],
                 W = tg$market_weight[mf])
    }
    margs <- args
    if (identical(type, "robust")) margs$weights <- d$weights^2
    D <- -s1$x
    rm(s1)
    meat <- do.call(mxl_blp_meat_parallel,
                    c(margs, list(D = D, supplied = !sample, cluster = cl,
                                  alpha = ka$alpha,
                                  market = if (!sample) tg$unit_market)))
  }
  .blp_vcov(info, jac, switch(type, robust = "sandwich", type), tv, meat, n_c,
            D = if (!(sample && model_based)) D)
}

#' The outer objective's Hessian of a share-matching fit, by differences
#'
#' Central differences of the outer gradient (the total derivative of the
#' negated log-likelihood along delta(theta_c), `.mxl_blp_solver()`'s), one
#' pair of solves of the constants per continuous parameter, each started
#' along the estimate's implicit derivative: the curvature of the objective
#' the optimizer minimized. The curvature along the share constraint that the
#' variance takes, It = T' I T with T = (identity; D), omits the constants'
#' score times the curvature of delta(theta_c), and that score is not zero
#' where the shares match, whatever the targets: the share equations average
#' each situation's choice probabilities over the draws, while the score
#' weights each draw by its posterior given the situation's choices. Under
#' `se_method = "bhhh"` It is also the BHHH matrix's, positive semidefinite
#' at any point.
#'
#' @param solve The inner solve, `.mxl_blp_solver()`'s function.
#' @param theta_c The continuous parameters at the estimate.
#' @param at The solve's record there (its `delta` and `extra`).
#' @param h The steps, one per continuous parameter.
#' @param max_iter The solves' evaluations of the shares.
#' @returns List: `hessian`, the differences' symmetrized matrix (n_c x n_c),
#'   or NULL with `why`, the reason, where a solve failed.
#' @noRd
.blp_outer_hessian <- function(solve, theta_c, at, h, max_iter) {
  n_c <- length(theta_c)
  best <- list(theta_c = theta_c, extra = at$extra)
  H <- matrix(0, n_c, n_c)
  for (k in seq_len(n_c)) {
    th <- list(theta_c, theta_c)
    th[[1L]][k] <- theta_c[k] + h[k]
    th[[2L]][k] <- theta_c[k] - h[k]
    g <- vector("list", 2L)
    for (i in 1:2) {
      r <- solve(th[[i]], at$delta, max_iter, best)
      if (!isTRUE(r$converged)) {
        return(list(hessian = NULL,
                    why = r$error %||% r$message %||% "the solve did not converge"))
      }
      g[[i]] <- r$gradient[seq_len(n_c)]
    }
    H[, k] <- (g[[1L]] - g[[2L]]) / (th[[1L]][k] - th[[2L]][k])
  }
  list(hessian = (H + t(H)) / 2, why = NULL)
}

#' The outer decrement of a share-matching fit
#'
#' Half the total gradient's squared length in a curvature of the outer
#' objective (the variance's along the share constraint, or the differences'
#' Hessian), over the mean weight: NA where that curvature is not positive
#' definite (no certificate then, as for the polish's decrement).
#'
#' @param It The curvature.
#' @param g_c The total gradient in theta_c (of the negated log-likelihood).
#' @param w_bar The mean weight.
#' @returns A number, or NA.
#' @noRd
.blp_outer_decrement <- function(It, g_c, w_bar) {
  tryCatch({
    chol(It)
    0.5 * sum(g_c * solve(It, g_c)) / w_bar
  }, error = function(e) NA_real_)
}

#' The share-matching constants' inner solve for a model and its draws
#'
#' The closures of `.mxl_blp_solver()` over the prepared data and the
#' estimation draws: the share inversion with Ti (each situation integrates
#' over its likelihood unit's draws), Anderson acceleration, exponentials
#' kept and compensated sums at scale (blp()'s passes), to `tol_shares`,
#' within 2^32 bytes of kept draws and exponentials; and the share Jacobian
#' kernel. The constants (the free ones, after the continuous parameters in
#' theta) start from the targets' log ratios to the reference's.
#'
#' @param d The prepared `input_data` (or a fit's stored data).
#' @param outside Whether the model has an outside option.
#' @param eta_draws The draws' cube (store mode), or an empty one.
#' @param gen List of `seed`, `scramble` and `S`, the generator's (seed -1
#'   in store mode).
#' @param rc_dist,rc_correlation,rc_mean The model's.
#' @param param_map,n_params The parameters' map and count.
#' @param tg `.blp_resolve_targets()`'s list.
#' @param ctl `.concentrate_control()`'s, for `"blp"`.
#' @returns List: `solve` (`.mxl_blp_solver()`'s), `start` (the free
#'   constants' start), `engine_f` and `jac_f`.
#' @noRd
.mxl_blp_setup <- function(d, outside, eta_draws, gen, rc_dist, rc_correlation,
                           rc_mean, param_map, n_params, tg, ctl) {
  c_idx <- setdiff(seq_len(n_params), param_map$asc)
  n_asc <- length(param_map$asc)
  # (values compared: the map's positions are doubles with rc_correlation)
  stopifnot(all(param_map$asc == length(c_idx) + seq_len(n_asc)))
  K_w <- ncol(d$W)
  S_eng <- tg$engine
  S_in <- if (outside) S_eng[-1L] else S_eng
  free <- if (outside) seq_len(n_asc) else seq_len(n_asc + 1L)[-1L]
  engine_f <- function(theta_c, delta0, max_iter) {
    co <- numeric(n_params)
    co[c_idx] <- theta_c
    r <- mxl_blp_solve(
      delta = if (outside) delta0 else c(0, delta0), target_shares = S_eng,
      X = d$X, W = d$W, beta = co[param_map$beta],
      mu = if (rc_mean) co[param_map$mu] else rep(0, K_w),
      L_params = co[param_map$sigma], alt_idx = d$alt_idx, M = d$M,
      weights = d$weights, eta_draws = eta_draws, rc_dist = rc_dist,
      rc_correlation = rc_correlation, rc_mean = rc_mean,
      include_outside_option = outside,
      control = .blp_control(ctl$tol_shares, max_iter, method = "anderson",
                             pass = .BLP_PASS, compensate = .BLP_COMPENSATE),
      gen_seed = gen$seed, gen_scramble = gen$scramble, gen_S = gen$S,
      cache_bytes = 2^32, Ti = d$Ti)
    r$delta_free <- if (outside) as.numeric(r$delta) else as.numeric(r$delta)[-1L]
    r
  }
  jac_f <- function(theta) {
    mxl_blp_jacobian_parallel(
      theta = theta, X = d$X, W = d$W, alt_idx = d$alt_idx,
      choice_idx = d$choice_idx, M = d$M, weights = d$weights,
      eta_draws = eta_draws, rc_dist = rc_dist, rc_correlation = rc_correlation,
      rc_mean = rc_mean, use_asc = TRUE, include_outside_option = outside,
      gen_seed = gen$seed, gen_scramble = gen$scramble, gen_S = gen$S,
      Ti = d$Ti)
  }
  list(solve = .mxl_blp_solver(engine_f, jac_f, length(c_idx), n_params,
                               tg$counts, free),
       start = log(S_in[free] / S_eng[1L]), engine_f = engine_f, jac_f = jac_f)
}
