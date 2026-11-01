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
