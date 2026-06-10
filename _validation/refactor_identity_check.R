#!/usr/bin/env Rscript
# ============================================================================
# Bit-identity harness for the C++ deduplication refactor.
#
# Calls every exported C++ entry point of mnlogit.cpp / mxlogit.cpp /
# nestlogit.cpp over a grid of flag combinations and several random theta
# points, then either saves the outputs as the reference baseline or compares
# the current build's outputs against that baseline with identical()
# (bit-exact). Also asserts error messages for invalid inputs.
#
# Usage:
#   Rscript _validation/refactor_identity_check.R baseline
#   Rscript _validation/refactor_identity_check.R check
#
# Must be run single-threaded: with >1 OpenMP thread the per-thread partial
# accumulators are combined under `omp critical` in nondeterministic order,
# so results are not bit-stable even without any code change.
# ============================================================================

suppressMessages(library(choicer))

args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args) >= 1) args[1] else ""
if (!mode %in% c("baseline", "check")) {
  stop("Usage: Rscript _validation/refactor_identity_check.R baseline|check")
}
baseline_path <- file.path("_validation", "refactor_baseline.rds")

getFromNamespace("set_num_threads", "choicer")(1)

get_halton_normals <- getFromNamespace("get_halton_normals", "choicer")

# ---------------------------------------------------------------------------
# Data generators (all RNG under fixed seeds)
# ---------------------------------------------------------------------------

# Unbalanced choice sets over J global alternatives. Individual 1 always sees
# the full set so every alternative has positive predicted share (needed for
# the BLP contraction targets).
make_choice_data <- function(seed, N, J, K, include_outside_option) {
  set.seed(seed)
  M <- integer(N)
  alt_list <- vector("list", N)
  for (i in seq_len(N)) {
    if (i == 1) {
      alts <- seq_len(J)
    } else {
      alts <- sort(sample.int(J, sample(2:J, 1)))
    }
    M[i] <- length(alts)
    alt_list[[i]] <- alts
  }
  alt_idx <- unlist(alt_list)
  total <- sum(M)
  X <- matrix(rnorm(total * K), total, K)
  choice_idx <- integer(N)
  for (i in seq_len(N)) {
    lo <- if (include_outside_option) 0L else 1L
    choice_idx[i] <- sample(lo:M[i], 1L)
  }
  weights <- runif(N, 0.5, 2)
  list(X = X, alt_idx = alt_idx, choice_idx = choice_idx, M = M,
       weights = weights, J = J, K = K, N = N, total = total)
}

n_delta_free <- function(J, use_asc, ioo) {
  if (!use_asc) 0L else if (ioo) J else J - 1L
}

# ---------------------------------------------------------------------------
# MNL battery
# ---------------------------------------------------------------------------

run_mnl <- function() {
  out <- list()
  d <- make_choice_data(seed = 101, N = 40, J = 5, K = 3,
                        include_outside_option = FALSE)
  d_oo <- make_choice_data(seed = 102, N = 40, J = 5, K = 3,
                           include_outside_option = TRUE)
  for (use_asc in c(TRUE, FALSE)) {
    for (ioo in c(TRUE, FALSE)) {
      dat <- if (ioo) d_oo else d
      n_par <- dat$K + n_delta_free(dat$J, use_asc, ioo)
      for (rep in 1:2) {
        set.seed(1000 + 10 * rep + ioo)
        theta <- rnorm(n_par, sd = 0.4)
        key <- sprintf("mnl|asc=%d|oo=%d|rep=%d", use_asc, ioo, rep)
        res <- list()
        res$grad <- mnl_loglik_gradient_parallel(
          theta, dat$X, dat$alt_idx, dat$choice_idx, dat$M, dat$weights,
          use_asc, ioo)
        res$pred <- mnl_predict(theta, dat$X, dat$alt_idx, dat$M, use_asc, ioo)
        res$shares <- mnl_predict_shares(
          theta, dat$X, dat$alt_idx, dat$M, dat$weights, use_asc, ioo)
        res$hess <- mnl_loglik_hessian_parallel(
          theta, dat$X, dat$alt_idx, dat$choice_idx, dat$M, dat$weights,
          use_asc, ioo)
        res$elas <- mnl_elasticities_parallel(
          theta, dat$X, dat$alt_idx, dat$choice_idx, dat$M, dat$weights,
          elast_var_idx = 1L, use_asc, ioo)
        res$dr <- mnl_diversion_ratios_parallel(
          theta, dat$X, dat$alt_idx, dat$M, dat$weights, use_asc, ioo)
        if (use_asc && rep == 1) {
          beta <- theta[seq_len(dat$K)]
          res$blp <- blp_contraction(
            rep(0, dat$J), res$shares, dat$X, beta, dat$alt_idx, dat$M,
            dat$weights, ioo)
        }
        out[[key]] <- res
      }
    }
  }
  out
}

# ---------------------------------------------------------------------------
# MXL battery
# ---------------------------------------------------------------------------

run_mxl <- function() {
  out <- list()
  N <- 30; J <- 4; K_x <- 2; K_w <- 2; S <- 8
  d <- make_choice_data(seed = 201, N = N, J = J, K = K_x,
                        include_outside_option = FALSE)
  d_oo <- make_choice_data(seed = 202, N = N, J = J, K = K_x,
                           include_outside_option = TRUE)
  set.seed(203)
  eta <- get_halton_normals(S, N, K_w)
  rc_dist_mixes <- list(normal = c(0L, 0L), mixed = c(0L, 1L),
                        lognormal = c(1L, 1L))

  for (use_asc in c(TRUE, FALSE)) {
    for (ioo in c(TRUE, FALSE)) {
      dat <- if (ioo) d_oo else d
      for (rc_corr in c(TRUE, FALSE)) {
        for (rc_mean in c(TRUE, FALSE)) {
          for (dist_name in names(rc_dist_mixes)) {
            rc_dist <- rc_dist_mixes[[dist_name]]
            for (w_type in c("row", "alt")) {
              set.seed(300 + 7 * ioo + 13 * rc_corr)
              W <- if (w_type == "row") {
                matrix(rnorm(dat$total * K_w), dat$total, K_w)
              } else {
                matrix(rnorm(J * K_w), J, K_w)
              }
              L_size <- if (rc_corr) K_w * (K_w + 1) / 2 else K_w
              n_par <- K_x + (if (rc_mean) K_w else 0L) + L_size +
                n_delta_free(J, use_asc, ioo)
              set.seed(400 + 3 * use_asc + 5 * rc_mean)
              theta <- rnorm(n_par, sd = 0.3)
              key <- sprintf(
                "mxl|asc=%d|oo=%d|corr=%d|mean=%d|dist=%s|W=%s",
                use_asc, ioo, rc_corr, rc_mean, dist_name, w_type)
              res <- list()
              res$grad <- mxl_loglik_gradient_parallel(
                theta, dat$X, W, dat$alt_idx, dat$choice_idx, dat$M,
                dat$weights, eta, rc_dist, rc_corr, rc_mean, use_asc, ioo)
              res$hess <- mxl_hessian_parallel(
                theta, dat$X, W, dat$alt_idx, dat$choice_idx, dat$M,
                dat$weights, eta, rc_dist, rc_corr, rc_mean, use_asc, ioo)
              res$bhhh <- mxl_bhhh_parallel(
                theta, dat$X, W, dat$alt_idx, dat$choice_idx, dat$M,
                dat$weights, eta, rc_dist, rc_corr, rc_mean, use_asc, ioo)
              res$pred <- mxl_predict(
                theta, dat$X, W, dat$alt_idx, dat$M, eta, rc_dist,
                rc_corr, rc_mean, use_asc, ioo)
              res$logsum <- mxl_logsum(
                theta, dat$X, W, dat$alt_idx, dat$M, eta, rc_dist,
                rc_corr, rc_mean, use_asc, ioo)
              res$shares <- mxl_predict_shares(
                theta, dat$X, W, dat$alt_idx, dat$M, dat$weights, eta,
                rc_dist, rc_corr, rc_mean, use_asc, ioo)
              res$elas_fixed <- mxl_elasticities_parallel(
                theta, dat$X, W, dat$alt_idx, dat$choice_idx, dat$M,
                dat$weights, eta, rc_dist, elast_var_idx = 1L,
                is_random_coef = FALSE, rc_corr, rc_mean, use_asc, ioo)
              res$elas_random <- mxl_elasticities_parallel(
                theta, dat$X, W, dat$alt_idx, dat$choice_idx, dat$M,
                dat$weights, eta, rc_dist, elast_var_idx = 1L,
                is_random_coef = TRUE, rc_corr, rc_mean, use_asc, ioo)
              res$dr_fixed <- mxl_diversion_ratios_parallel(
                theta, dat$X, W, dat$alt_idx, dat$M, dat$weights, eta,
                rc_dist, elast_var_idx = 1L, is_random_coef = FALSE,
                rc_corr, rc_mean, use_asc, ioo)
              res$dr_random <- mxl_diversion_ratios_parallel(
                theta, dat$X, W, dat$alt_idx, dat$M, dat$weights, eta,
                rc_dist, elast_var_idx = 2L, is_random_coef = TRUE,
                rc_corr, rc_mean, use_asc, ioo)
              if (dist_name == "normal") {
                beta <- theta[seq_len(K_x)]
                mu <- if (rc_mean) theta[K_x + seq_len(K_w)] else rep(0, K_w)
                L_params <- theta[K_x + (if (rc_mean) K_w else 0L) +
                                    seq_len(L_size)]
                res$blp <- mxl_blp_contraction(
                  rep(0, J), res$shares, dat$X, W, beta, mu, L_params,
                  dat$alt_idx, dat$M, dat$weights, eta, rc_dist,
                  rc_corr, rc_mean, ioo)
              }
              out[[key]] <- res
            }
          }
        }
      }
    }
  }

  # Cholesky utilities
  set.seed(500)
  for (rc_corr in c(TRUE, FALSE)) {
    L_size <- if (rc_corr) 6L else 3L
    L_params <- rnorm(L_size, sd = 0.5)
    key <- sprintf("mxl_chol|corr=%d", rc_corr)
    out[[key]] <- list(
      L = getFromNamespace("build_L_mat", "choicer")(L_params, 3L, rc_corr),
      Sigma = build_var_mat(L_params, 3L, rc_corr),
      jac = jacobian_vech_Sigma(L_params, 3L, rc_corr)
    )
  }
  out
}

# ---------------------------------------------------------------------------
# NL battery
# ---------------------------------------------------------------------------

run_nl <- function() {
  out <- list()
  N <- 40; J <- 5; K <- 3
  # nests: {1,2} -> 1, {3,4} -> 2, {5} -> 3 (singleton)
  nest_idx <- c(1L, 1L, 2L, 2L, 3L)
  n_nonsing <- 2L
  d <- make_choice_data(seed = 301, N = N, J = J, K = K,
                        include_outside_option = FALSE)
  d_oo <- make_choice_data(seed = 302, N = N, J = J, K = K,
                           include_outside_option = TRUE)
  for (use_asc in c(TRUE, FALSE)) {
    for (ioo in c(TRUE, FALSE)) {
      dat <- if (ioo) d_oo else d
      n_par <- K + n_nonsing + n_delta_free(J, use_asc, ioo)
      for (rep in 1:2) {
        set.seed(2000 + 10 * rep + ioo)
        theta <- rnorm(n_par, sd = 0.4)
        theta[K + seq_len(n_nonsing)] <- runif(n_nonsing, 0.6, 1.0)
        key <- sprintf("nl|asc=%d|oo=%d|rep=%d", use_asc, ioo, rep)
        res <- list()
        res$grad <- nl_loglik_gradient_parallel(
          theta, dat$X, dat$alt_idx, dat$choice_idx, nest_idx, dat$M,
          dat$weights, use_asc, ioo)
        res$pred <- nl_predict(
          theta, dat$X, dat$alt_idx, dat$M, nest_idx, use_asc, ioo)
        res$shares <- nl_predict_shares(
          theta, dat$X, dat$alt_idx, dat$M, dat$weights, nest_idx,
          use_asc, ioo)
        res$elas <- nl_elasticities_parallel(
          theta, dat$X, dat$alt_idx, dat$choice_idx, nest_idx, dat$M,
          dat$weights, elast_var_idx = 1L, use_asc, ioo)
        res$dr <- nl_diversion_ratios_parallel(
          theta, dat$X, dat$alt_idx, nest_idx, dat$M, dat$weights,
          use_asc, ioo)
        if (rep == 1) {
          res$num_hess <- nl_loglik_numeric_hessian(
            theta, dat$X, dat$alt_idx, dat$choice_idx, nest_idx, dat$M,
            dat$weights, use_asc, ioo)
        }
        if (use_asc && rep == 1) {
          beta <- theta[seq_len(K)]
          lambda_full <- rep(1, 3)
          lambda_full[1:2] <- theta[K + 1:2]
          res$blp <- nl_blp_contraction(
            rep(0, J), res$shares, dat$X, beta, lambda_full, dat$alt_idx,
            nest_idx, dat$M, dat$weights, ioo)
        }
        out[[key]] <- res
      }
    }
  }
  out
}

# ---------------------------------------------------------------------------
# Error-message subharness
# ---------------------------------------------------------------------------

run_errors <- function() {
  grab <- function(expr) tryCatch({ expr; "<<no error>>" },
                                  error = function(e) conditionMessage(e))
  out <- list()
  d <- make_choice_data(seed = 901, N = 10, J = 4, K = 2,
                        include_outside_option = FALSE)
  nest_idx <- c(1L, 1L, 2L, 2L)
  set.seed(902)
  eta <- get_halton_normals(4, d$N, 2)
  W <- matrix(rnorm(d$total * 2), d$total, 2)

  # MNL: ASC expected but theta has only beta
  out$mnl_missing_asc <- grab(mnl_loglik_gradient_parallel(
    rnorm(2), d$X, d$alt_idx, d$choice_idx, d$M, d$weights, TRUE, FALSE))
  # MXL: theta too short for L block (gradient path)
  out$mxl_short_L <- grab(mxl_loglik_gradient_parallel(
    rnorm(2), d$X, W, d$alt_idx, d$choice_idx, d$M, d$weights, eta,
    c(0L, 0L), FALSE, FALSE, FALSE, FALSE))
  # MXL: theta too short for delta block (gradient path, use_asc)
  out$mxl_short_delta <- grab(mxl_loglik_gradient_parallel(
    rnorm(4), d$X, W, d$alt_idx, d$choice_idx, d$M, d$weights, eta,
    c(0L, 0L), FALSE, FALSE, TRUE, FALSE))
  # MXL: delta stop on a predict-family path
  out$mxl_predict_short_delta <- grab(mxl_predict(
    rnorm(4), d$X, W, d$alt_idx, d$M, eta, c(0L, 0L),
    FALSE, FALSE, TRUE, FALSE))
  # MXL: bad rc_dist length
  out$mxl_bad_rcdist <- grab(mxl_loglik_gradient_parallel(
    rnorm(7), d$X, W, d$alt_idx, d$choice_idx, d$M, d$weights, eta,
    c(0L), FALSE, FALSE, TRUE, FALSE))
  # MXL: mismatched eta draws (wrong N)
  set.seed(903)
  eta_bad <- get_halton_normals(4, d$N + 1, 2)
  out$mxl_bad_eta <- grab(mxl_loglik_gradient_parallel(
    rnorm(7), d$X, W, d$alt_idx, d$choice_idx, d$M, d$weights, eta_bad,
    c(0L, 0L), FALSE, FALSE, TRUE, FALSE))
  # NL: non-positive lambda
  out$nl_bad_lambda <- grab(nl_loglik_gradient_parallel(
    c(rnorm(2), -0.5, 0.8, rnorm(3)), d$X, d$alt_idx, d$choice_idx,
    nest_idx, d$M, d$weights, TRUE, FALSE))
  # NL: theta too short for lambda block
  out$nl_short <- grab(nl_loglik_gradient_parallel(
    rnorm(3), d$X, d$alt_idx, d$choice_idx, nest_idx, d$M, d$weights,
    FALSE, FALSE))
  out
}

# ---------------------------------------------------------------------------
# Run, save or compare
# ---------------------------------------------------------------------------

cat("Running MNL battery...\n");   res_mnl <- run_mnl()
cat("Running MXL battery...\n");   res_mxl <- run_mxl()
cat("Running NL battery...\n");    res_nl <- run_nl()
cat("Running error battery...\n"); res_err <- run_errors()

results <- list(mnl = res_mnl, mxl = res_mxl, nl = res_nl, errors = res_err)

if (mode == "baseline") {
  saveRDS(results, baseline_path)
  cat(sprintf("Baseline saved to %s (%d MNL, %d MXL, %d NL configs, %d error cases)\n",
              baseline_path, length(res_mnl), length(res_mxl), length(res_nl),
              length(res_err)))
  quit(status = 0)
}

# mode == "check"
if (!file.exists(baseline_path)) {
  stop("No baseline found at ", baseline_path,
       " — run with 'baseline' on the pre-refactor build first.")
}
baseline <- readRDS(baseline_path)

flatten <- function(x, prefix = "") {
  if (is.list(x)) {
    nms <- names(x)
    if (is.null(nms)) nms <- as.character(seq_along(x))
    out <- list()
    for (i in seq_along(x)) {
      out <- c(out, flatten(x[[i]], paste0(prefix, "$", nms[i])))
    }
    out
  } else {
    stats::setNames(list(x), prefix)
  }
}

base_flat <- flatten(baseline)
new_flat  <- flatten(results)

fail <- FALSE
if (!identical(sort(names(base_flat)), sort(names(new_flat)))) {
  fail <- TRUE
  cat("MISMATCH in output structure:\n")
  cat("  only in baseline:", setdiff(names(base_flat), names(new_flat)), "\n")
  cat("  only in new:     ", setdiff(names(new_flat), names(base_flat)), "\n")
}
n_checked <- 0L
for (key in intersect(names(base_flat), names(new_flat))) {
  n_checked <- n_checked + 1L
  if (!identical(base_flat[[key]], new_flat[[key]])) {
    fail <- TRUE
    b <- base_flat[[key]]; n <- new_flat[[key]]
    if (is.numeric(b) && is.numeric(n) && length(b) == length(n)) {
      dev <- max(abs(b - n), na.rm = TRUE)
      cat(sprintf("NOT IDENTICAL: %s  (max abs dev: %.3e)\n", key, dev))
    } else {
      cat(sprintf("NOT IDENTICAL: %s  (type/length mismatch)\n", key))
      cat("  all.equal:", paste(all.equal(b, n), collapse = "; "), "\n")
    }
  }
}

if (fail) {
  cat("\nFAILED: outputs are not bit-identical to baseline.\n")
  quit(status = 1)
}
cat(sprintf("\nOK: all %d outputs bit-identical to baseline.\n", n_checked))
quit(status = 0)
