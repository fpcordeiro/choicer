# Kernel tests for the panel mixed logit likelihood (Revelt & Train 1998).
#
# The four MXL estimation kernels and mxl_conditional_tastes_parallel() take a
# decision-maker structure `Ti` (the number of consecutive choice situations of
# each decision maker). Every check compares the kernels with a computation
# that shares no code with the C++:
#   * the brute-force oracle in helper-mxl-panel.R, derived from the estimand
#     (objective, per-decision-maker log-likelihood, conditional tastes);
#   * numerical derivatives of the kernel objective and of the oracle;
#   * exact identities: BHHH and the gradient as weighted sums of the unit
#     score rows, Hessian symmetry and linearity in the weights, Ti = rep(1, N)
#     as the cross-section, generate mode as store mode, and the collapse of
#     the panel likelihood to a cross-section on shared draws when S = 1.
# The cells of mxlp_cell_configs() are small unbalanced panels (T_u in 1..4,
# 1 to 4 inside alternatives per situation) that jointly cover every model
# option; each failure message names its cell.

mxlp_cells <- mxlp_build_cells()

# --- Objective, gradient and Hessian -----------------------------------------

test_that("panel objective equals the brute-force oracle", {
  for (fx in mxlp_cells) {
    res <- mxlp_call("gradient", fx)
    orc <- mxlp_oracle(fx$theta, fx)
    expect_length(res$gradient, length(fx$theta))
    mxlp_expect_close(res$objective, -orc$loglik, 1e-10,
                      sprintf("[%s] objective vs oracle", fx$name))
  }
})

test_that("panel gradient matches numDeriv of kernel objective and oracle", {
  skip_if_not_installed("numDeriv")
  for (fx in mxlp_cells) {
    grad <- mxlp_call("gradient", fx)$gradient
    obj_kernel <- function(th) mxlp_call("gradient", fx, theta = th)$objective
    obj_oracle <- function(th) -mxlp_oracle(th, fx)$loglik
    mxlp_expect_close(grad, numDeriv::grad(obj_kernel, fx$theta), TOL_GRAD,
                      sprintf("[%s] gradient vs numDeriv(kernel objective)",
                              fx$name))
    mxlp_expect_close(grad, numDeriv::grad(obj_oracle, fx$theta), TOL_GRAD,
                      sprintf("[%s] gradient vs numDeriv(oracle)", fx$name))
  }
})

test_that("panel Hessian matches numDeriv of the kernel objective", {
  skip_if_not_installed("numDeriv")
  for (fx in mxlp_cells) {
    H <- mxlp_call("hessian", fx)
    p <- length(fx$theta)
    expect_equal(dim(H), c(p, p), label = sprintf("[%s] dim(H)", fx$name))
    obj <- function(th) mxlp_call("gradient", fx, theta = th)$objective
    mxlp_expect_close(H, numDeriv::hessian(obj, fx$theta), TOL_HESS,
                      sprintf("[%s] Hessian vs numDeriv", fx$name))
  }
})

test_that("panel Hessian is symmetric and linear in the weights", {
  for (fx in mxlp_cells) {
    H <- mxlp_call("hessian", fx)
    mxlp_expect_close(H, t(H), 1e-10,
                      sprintf("[%s] Hessian symmetry", fx$name))
    H2 <- mxlp_call("hessian", fx, weights = 2 * fx$weights)
    mxlp_expect_close(H2, 2 * H, 1e-10,
                      sprintf("[%s] H(2w) vs 2 H(w)", fx$name))
  }
})

# --- Unit scores, BHHH and conditional tastes --------------------------------

test_that("BHHH and the gradient are weighted sums over the unit score rows", {
  for (fx in mxlp_cells) {
    S_u <- mxlp_call("scores", fx)
    w_u <- mxlp_unit_weights(fx)
    expect_equal(dim(S_u), c(fx$U, length(fx$theta)),
                 label = sprintf("[%s] dim(scores)", fx$name))
    mxlp_expect_close(mxlp_call("bhhh", fx), crossprod(sqrt(w_u) * S_u),
                      1e-10, sprintf("[%s] BHHH vs sum_u w_u s_u s_u'",
                                     fx$name))
    # The kernel returns the gradient of the negated log-likelihood.
    mxlp_expect_close(colSums(w_u * S_u), -mxlp_call("gradient", fx)$gradient,
                      1e-10, sprintf("[%s] sum_u w_u s_u vs -gradient",
                                     fx$name))
  }
})

test_that("score rows match numDeriv of each decision maker's log-likelihood", {
  skip_if_not_installed("numDeriv")
  for (fx in mxlp_cells) {
    S_u <- mxlp_call("scores", fx)
    unit_of <- rep(seq_len(fx$U), fx$Ti)
    for (u in fx$probe_units) {
      what <- sprintf("[%s] score row %d (T_u = %d)", fx$name, u, fx$Ti[u])
      g_oracle <- numDeriv::grad(mxlp_oracle_unit_loglik, fx$theta,
                                 fx = fx, u = u)
      mxlp_expect_close(S_u[u, ], g_oracle, 1e-6,
                        paste(what, "vs numDeriv(oracle l_u)"))
      # Indicator weights, constant within every unit, select decision
      # maker u in the kernel's own objective.
      w_ind <- as.numeric(unit_of == u)
      loglik_u <- function(th) {
        -mxlp_call("gradient", fx, theta = th, weights = w_ind)$objective
      }
      mxlp_expect_close(S_u[u, ], numDeriv::grad(loglik_u, fx$theta), 1e-6,
                        paste(what, "vs numDeriv(indicator-weighted kernel)"))
    }
  }
})

test_that("conditional tastes equal the oracle's posterior draw moments", {
  for (fx in mxlp_cells) {
    ct <- mxlp_call("tastes", fx)
    orc <- mxlp_oracle(fx$theta, fx)
    expect_named(ct, c("mean", "sd"), ignore.order = TRUE)
    expect_equal(dim(ct$mean), c(fx$K_w, fx$U),
                 label = sprintf("[%s] dim(tastes$mean)", fx$name))
    expect_equal(dim(ct$sd), c(fx$K_w, fx$U),
                 label = sprintf("[%s] dim(tastes$sd)", fx$name))
    mxlp_expect_close(ct$mean, orc$mean, 1e-10,
                      sprintf("[%s] conditional mean vs oracle", fx$name))
    mxlp_expect_close(ct$sd, orc$sd, 1e-10,
                      sprintf("[%s] conditional SD vs oracle", fx$name))
  }
})

# --- Equivalences -----------------------------------------------------------

test_that("Ti = rep(1, N) reproduces the cross-sectional kernels", {
  for (fx in mxlp_cells) {
    xs <- mxlp_cross_section(fx)
    ones <- rep(1L, xs$N)
    for (k in mxlp_kernels) {
      res_null <- mxlp_call(k, xs, Ti = NULL)
      res_ones <- mxlp_call(k, xs, Ti = ones)
      mxlp_expect_close(res_ones, res_null, 1e-12,
                        sprintf("[%s] %s, Ti = rep(1, N) vs NULL",
                                fx$name, k))
    }
    orc <- mxlp_oracle(xs$theta, xs)
    mxlp_expect_close(mxlp_call("gradient", xs)$objective, -orc$loglik, 1e-10,
                      sprintf("[%s] cross-section objective vs oracle",
                              fx$name))
    ct <- mxlp_call("tastes", xs)
    expect_equal(dim(ct$mean), c(xs$K_w, xs$N),
                 label = sprintf("[%s] cross-section dim(tastes$mean)",
                                 fx$name))
    mxlp_expect_close(ct$mean, orc$mean, 1e-10,
                      sprintf("[%s] cross-section conditional mean vs oracle",
                              fx$name))
    mxlp_expect_close(ct$sd, orc$sd, 1e-10,
                      sprintf("[%s] cross-section conditional SD vs oracle",
                              fx$name))
  }
})

test_that("generate mode reproduces store-mode panel draws", {
  # gen_scramble = 0 generates plain Halton normals: decision maker u gets
  # the same block as slice u of get_halton_normals(S, U, K_w).
  for (fx in mxlp_cells) {
    for (k in mxlp_kernels) {
      mxlp_expect_close(mxlp_call(k, fx, generate = TRUE), mxlp_call(k, fx),
                        1e-12, sprintf("[%s] %s, generate vs store",
                                       fx$name, k))
    }
  }
})

test_that("with S = 1 the panel kernels are a cross-section on shared draws", {
  # With S = 1, l_u = sum_t log P_t(c_u), so the panel objective, gradient and
  # Hessian equal the cross-sectional ones when each situation uses its
  # decision maker's draw; unit scores are within-unit sums of situation
  # scores; and the conditional tastes are that draw, with zero SD.
  for (fx in mxlp_cells) {
    f1 <- fx
    f1$S <- 1L
    f1$eta <- get_halton_normals(1L, fx$U, fx$K_w)
    unit_of <- rep(seq_len(fx$U), fx$Ti)
    xs <- f1
    xs[["Ti"]] <- NULL
    xs$eta <- f1$eta[, , unit_of, drop = FALSE]

    pan <- mxlp_call("gradient", f1)
    crs <- mxlp_call("gradient", xs)
    mxlp_expect_close(pan$objective, crs$objective, 1e-10,
                      sprintf("[%s] S = 1 objective", fx$name))
    mxlp_expect_close(pan$gradient, crs$gradient, 1e-10,
                      sprintf("[%s] S = 1 gradient", fx$name))
    mxlp_expect_close(mxlp_call("hessian", f1), mxlp_call("hessian", xs),
                      1e-10, sprintf("[%s] S = 1 Hessian", fx$name))
    mxlp_expect_close(mxlp_call("scores", f1),
                      rowsum(mxlp_call("scores", xs), unit_of), 1e-10,
                      sprintf("[%s] S = 1 unit scores", fx$name))
    ct <- mxlp_call("tastes", f1)
    mxlp_expect_close(ct$mean, mxlp_oracle(f1$theta, f1)$mean, 1e-10,
                      sprintf("[%s] S = 1 conditional mean", fx$name))
    mxlp_expect_close(ct$sd, matrix(0, fx$K_w, fx$U), 1e-10,
                      sprintf("[%s] S = 1 conditional SD", fx$name))
  }
})

# --- Numerical stability -----------------------------------------------------

test_that("panel kernels stay accurate when a unit's likelihood underflows", {
  # Decision maker 2 makes 150 low-probability choices, so exp(lambda_us) is
  # exactly 0 in double precision for every draw. The panel likelihood and
  # the posterior draw weights must be computed in log space.
  skip_if_not_installed("numDeriv")
  fx <- mxlp_long_panel_fixture()
  orc <- mxlp_oracle(fx$theta, fx)
  expect_true(all(exp(orc$lambda[2L, ]) == 0))  # the fixture does underflow

  res <- mxlp_call("gradient", fx)
  mxlp_expect_close(res$objective, -orc$loglik, 1e-10,
                    "[long panel] objective vs oracle")
  obj <- function(th) mxlp_call("gradient", fx, theta = th)$objective
  mxlp_expect_close(res$gradient, numDeriv::grad(obj, fx$theta), TOL_GRAD,
                    "[long panel] gradient vs numDeriv")
  mxlp_expect_close(mxlp_call("hessian", fx), numDeriv::hessian(obj, fx$theta),
                    TOL_HESS, "[long panel] Hessian vs numDeriv")

  S_u <- mxlp_call("scores", fx)
  mxlp_expect_close(colSums(S_u), -res$gradient, 1e-10,
                    "[long panel] sum_u s_u vs -gradient")
  mxlp_expect_close(mxlp_call("bhhh", fx), crossprod(S_u), 1e-10,
                    "[long panel] BHHH vs sum_u s_u s_u'")
  ct <- mxlp_call("tastes", fx)
  mxlp_expect_close(ct$mean, orc$mean, 1e-10,
                    "[long panel] conditional mean vs oracle")
  mxlp_expect_close(ct$sd, orc$sd, 1e-10,
                    "[long panel] conditional SD vs oracle")
})

# The edge fixtures put one chosen alternative `gap` utils below its
# competitors, once in a cross-section and once inside a multi-situation
# decision maker (see mxlp_edge_fixture()).

test_that("Hessian keeps a unit whose simulated likelihood is below 1e-12", {
  # A gap of -36 puts P(chosen) near 1e-16 at every draw: tiny but finite, so
  # the unit belongs in the Hessian. The pre-panel kernel dropped every unit
  # whose simulated probability was below 1e-12.
  skip_if_not_installed("numDeriv")
  for (panel in c(FALSE, TRUE)) {
    fx <- mxlp_edge_fixture(-36, panel)
    lik <- sum(exp(mxlp_oracle(fx$theta, fx)$lambda[fx$edge_unit, ]))
    expect_true(lik > 0 && lik < 1e-12,
                label = sprintf("[%s] summed unit likelihood in (0, 1e-12)",
                                fx$name))
    grad <- function(th) drop(mxlp_call("gradient", fx, theta = th)$gradient)
    mxlp_expect_close(mxlp_call("hessian", fx),
                      numDeriv::jacobian(grad, fx$theta), TOL_HESS,
                      sprintf("[%s] Hessian vs numDeriv::jacobian(gradient)",
                              fx$name))
  }
})

test_that("subnormal choice probabilities keep the objective and scores exact", {
  # P(chosen) is subnormal (positive, below .Machine$double.xmin) at every
  # draw. The kernels form log P(chosen) as the chosen utility minus the
  # log-sum-exp denominator, never as the log of the subnormal probability,
  # so the objective and the unit's score stay exact even at a gap of -740,
  # where P(chosen) itself keeps only about 6 significant bits. (Earlier
  # kernels took log(P(chosen)) and lost that precision; before the panel
  # work they also returned a zero gradient here.)
  skip_if_not_installed("numDeriv")
  for (gap in c(-720, -740)) {
    for (panel in c(FALSE, TRUE)) {
      fx <- mxlp_edge_fixture(gap, panel)
      u <- fx$edge_unit
      p <- exp(mxlp_situation_logp(fx$theta, fx, fx$edge_situation))
      expect_true(all(p > 0 & p < .Machine$double.xmin),
                  label = sprintf("[%s] P(chosen) subnormal at every draw",
                                  fx$name))
      res <- mxlp_call("gradient", fx)
      S_u <- mxlp_call("scores", fx)
      B <- mxlp_call("bhhh", fx)
      mxlp_expect_close(res$objective, -mxlp_oracle(fx$theta, fx)$loglik,
                        1e-10, sprintf("[%s] objective vs oracle", fx$name))
      expect_true(all(is.finite(S_u)) && all(is.finite(B)),
                  label = sprintf("[%s] finite scores and BHHH", fx$name))
      w_u <- mxlp_unit_weights(fx)
      mxlp_expect_close(colSums(w_u * S_u), -res$gradient, 1e-10,
                        sprintf("[%s] sum_u w_u s_u vs -gradient", fx$name))
      mxlp_expect_close(B, crossprod(sqrt(w_u) * S_u), 1e-10,
                        sprintf("[%s] BHHH vs sum_u w_u s_u s_u'", fx$name))
      g_oracle <- numDeriv::grad(mxlp_oracle_unit_loglik, fx$theta,
                                 fx = fx, u = u)
      mxlp_expect_close(S_u[u, ], g_oracle, 1e-6,
                        sprintf("[%s] score row %d vs numDeriv(oracle l_u)",
                                fx$name, u))
    }
  }
})

test_that("a choice probability that underflows to zero keeps an exact log-likelihood", {
  # Gaps of -1000 and -1e5 make P(chosen) exactly 0 in double precision at
  # every draw. In log space the choice keeps a finite, exact log
  # probability, so the unit's log-likelihood and its derivatives stay
  # exact: no optimizer sentinel, the Hessian keeps the unit, and its
  # conditional tastes are defined. With 10^7 choice situations some sit this
  # far in the tail, and one such unit used to turn the whole objective into
  # the sentinel. Checked for an inside choice with and without the outside
  # option, and for the outside option chosen against inside alternatives
  # far above it.
  skip_if_not_installed("numDeriv")
  cases <- expand.grid(layout = c("inside", "no_outside", "outside"),
                       panel = c(FALSE, TRUE), stringsAsFactors = FALSE)
  for (k in seq_len(nrow(cases))) {
    G <- H <- list()
    for (gap in c(-1000, -1e5)) {
      fx <- mxlp_edge_fixture(gap, cases$panel[k], cases$layout[k])
      p <- exp(mxlp_situation_logp(fx$theta, fx, fx$edge_situation))
      expect_true(all(p == 0),
                  label = sprintf("[%s] P(chosen) = 0 at every draw", fx$name))

      orc <- mxlp_oracle(fx$theta, fx)
      res <- mxlp_call("gradient", fx)
      mxlp_expect_close(res$objective, -orc$loglik, 1e-10,
                        sprintf("[%s] objective vs oracle", fx$name))
      mxlp_expect_close(colSums(mxlp_unit_weights(fx) *
                                  mxlp_call("scores", fx)),
                        -res$gradient, 1e-10,
                        sprintf("[%s] sum_u w_u s_u vs -gradient", fx$name))
      ct <- mxlp_call("tastes", fx)
      mxlp_expect_close(ct$mean, orc$mean, 1e-10,
                        sprintf("[%s] conditional mean vs oracle", fx$name))
      mxlp_expect_close(ct$sd, orc$sd, 1e-10,
                        sprintf("[%s] conditional SD vs oracle", fx$name))
      G[[length(G) + 1L]] <- res$gradient
      H[[length(H) + 1L]] <- mxlp_call("hessian", fx)
      if (gap == -1000) {
        obj <- function(th) -mxlp_oracle(th, fx)$loglik
        mxlp_expect_close(res$gradient, numDeriv::grad(obj, fx$theta),
                          TOL_GRAD,
                          sprintf("[%s] gradient vs numDeriv(oracle)", fx$name))
        grad <- function(th) drop(mxlp_call("gradient", fx, theta = th)$gradient)
        J <- numDeriv::jacobian(grad, fx$theta,
                                method.args = list(d = 1e-3, r = 6))
        mxlp_expect_close(H[[1L]], J, TOL_HESS,
                          sprintf("[%s] Hessian vs numDeriv::jacobian(gradient)",
                                  fx$name))
      }
    }
    # Once P(chosen) is 0 at every draw, the gap enters the unit's
    # log-likelihood only linearly, through beta_1 * x_1, so moving it from
    # -1000 to -1e5 shifts the gradient only along beta_1, by
    # w_u (1e5 - 1000) / beta_1, and leaves the Hessian unchanged. (At -1e5
    # the objective is about 1e5, too large for finite differences at these
    # tolerances. An uncentered Louis identity would miss the Hessian check by
    # about 0.15: the draw weights sum to one only to ulp(lambda), about
    # 1e-11, times g^2 = 4e10.)
    w_u <- mxlp_unit_weights(fx)[fx$edge_unit]
    shift <- replace(numeric(length(fx$theta)), 1L,
                     w_u * (1e5 - 1000) / fx$theta[1L])
    mxlp_expect_close(G[[2L]] - G[[1L]], shift, 1e-9,
                      sprintf("[%s] gradient shift from gap -1000 to -1e5",
                              fx$name))
    mxlp_expect_close(H[[2L]], H[[1L]], 1e-4,
                      sprintf("[%s] Hessian at gap -1e5 vs gap -1000",
                              fx$name))
  }
})

test_that("overflowing utilities trigger the sentinel, the Hessian skip and NA tastes", {
  # With log-space probabilities a unit's log-likelihood is non-finite only
  # when its utilities are. Here the chosen alternative's fixed utility
  # overflows to -Inf (x = -1e308, beta_1 = 2), so every draw gives the
  # unit's choices zero probability. The objective returns the optimizer
  # sentinel, the Hessian leaves the unit out, and the unit's conditional
  # tastes are NA while every other unit's stay exact.
  for (panel in c(FALSE, TRUE)) {
    fx <- mxlp_edge_fixture(-36, panel)
    u <- fx$edge_unit
    t <- fx$edge_situation
    fx$theta[1L] <- 2
    fx$X[sum(fx$M[seq_len(t - 1L)]) + fx$choice_idx[t], 1L] <- -1e308
    expect_true(all(mxlp_situation_logp(fx$theta, fx, t) == -Inf),
                label = sprintf("[%s] chosen utility is -Inf at every draw",
                                fx$name))

    res <- mxlp_call("gradient", fx)
    expect_equal(res$objective, 1e10,
                 label = sprintf("[%s] sentinel objective", fx$name))
    expect_true(all(res$gradient == 0),
                label = sprintf("[%s] sentinel gradient is zero", fx$name))

    mxlp_expect_close(mxlp_call("hessian", fx),
                      mxlp_call("hessian", mxlp_drop_unit(fx, u)), 1e-10,
                      sprintf("[%s] Hessian vs Hessian without unit %d",
                              fx$name, u))

    ct <- mxlp_call("tastes", fx)
    orc <- mxlp_oracle(fx$theta, fx)
    expect_true(all(is.na(ct$mean[, u])) && all(is.na(ct$sd[, u])),
                label = sprintf("[%s] conditional tastes of unit %d are NA",
                                fx$name, u))
    mxlp_expect_close(ct$mean[, -u], orc$mean[, -u], 1e-10,
                      sprintf("[%s] other units' conditional mean vs oracle",
                              fx$name))
    mxlp_expect_close(ct$sd[, -u], orc$sd[, -u], 1e-10,
                      sprintf("[%s] other units' conditional SD vs oracle",
                              fx$name))
  }
})

# --- Input validation --------------------------------------------------------

test_that("panel kernels reject non-positive or missing Ti entries", {
  fx <- mxlp_cells[[2L]]
  Ti <- fx$Ti
  # Same length and sum as Ti, so only positivity is violated. Uniform
  # weights keep the regrouped units weight-constant.
  Ti_zero <- Ti
  Ti_zero[2L] <- Ti[1L] + Ti[2L]
  Ti_zero[1L] <- 0L
  Ti_neg <- Ti
  Ti_neg[2L] <- Ti[1L] + Ti[2L] + 1L
  Ti_neg[1L] <- -1L
  Ti_na <- Ti
  Ti_na[2L] <- NA_integer_
  w1 <- rep(1, fx$N)
  for (k in mxlp_kernels) {
    expect_error(mxlp_call(k, fx, Ti = Ti_zero, weights = w1),
                 "Ti must be positive", label = paste(k, "with Ti[1] = 0"))
    expect_error(mxlp_call(k, fx, Ti = Ti_neg, weights = w1),
                 "Ti must be positive", label = paste(k, "with Ti[1] = -1"))
    expect_error(mxlp_call(k, fx, Ti = Ti_na, weights = w1),
                 "Ti must be positive", label = paste(k, "with Ti[2] = NA"))
  }
})

test_that("panel kernels reject an empty Ti (needs the empty-Ti guard)", {
  # Ti = integer(0) describes no decision maker at all; it must not fall
  # back to the cross-section. The draw cube is empty too, so no other check
  # can fire first. Each kernel's message is captured so that every kernel
  # reports its own result.
  fx <- mxlp_cells[[2L]]
  eta_empty <- fx$eta[, , integer(0), drop = FALSE]
  for (k in mxlp_kernels) {
    msg <- tryCatch({
      mxlp_call(k, fx, Ti = integer(0), weights = rep(1, fx$N),
                eta = eta_empty)
      "<no error>"
    }, error = conditionMessage)
    expect_match(msg, "at least one respondent",
                 label = paste(k, "error message with Ti = integer(0)"))
  }
})

test_that("panel kernels reject Ti not summing to the number of situations", {
  fx <- mxlp_cells[[2L]]
  Ti_long <- fx$Ti
  Ti_long[fx$U] <- Ti_long[fx$U] + 1L
  u <- fx$probe_units[2L]  # an interior unit with T_u >= 2
  Ti_short <- fx$Ti
  Ti_short[u] <- Ti_short[u] - 1L
  w1 <- rep(1, fx$N)
  for (k in mxlp_kernels) {
    expect_error(mxlp_call(k, fx, Ti = Ti_long, weights = w1),
                 "does not match the number of choice situations",
                 label = paste(k, "with sum(Ti) = N + 1"))
    expect_error(mxlp_call(k, fx, Ti = Ti_short, weights = w1),
                 "does not match the number of choice situations",
                 label = paste(k, "with sum(Ti) = N - 1"))
  }
})

test_that("panel kernels reject weights that vary within a decision maker", {
  fx <- mxlp_cells[[2L]]
  u <- fx$probe_units[2L]  # an interior unit with T_u >= 2
  first <- sum(fx$Ti[seq_len(u - 1L)]) + 1L
  w_bad <- fx$weights
  w_bad[first + 1L] <- 1.5 * w_bad[first + 1L]
  for (k in c("gradient", "hessian", "bhhh")) {
    expect_error(mxlp_call(k, fx, weights = w_bad),
                 "constant within each decision maker",
                 label = paste(k, "with a within-unit weight change"))
  }
})

test_that("store-mode draws need one slice per decision maker", {
  fx <- mxlp_cells[[2L]]
  eta_situations <- get_halton_normals(fx$S, fx$N, fx$K_w)  # N != U slices
  eta_short <- fx$eta[, , -1L, drop = FALSE]                  # U - 1 slices
  for (k in mxlp_kernels) {
    expect_error(mxlp_call(k, fx, eta = eta_situations),
                 "eta_draws 3rd dimension",
                 label = paste(k, "with one slice per situation"))
    expect_error(mxlp_call(k, fx, eta = eta_short),
                 "eta_draws 3rd dimension",
                 label = paste(k, "with U - 1 slices"))
  }
  # The cross-section keeps its check: one slice per situation.
  expect_error(mxlp_call("gradient", fx, Ti = NULL, eta = fx$eta),
               "eta_draws 3rd dimension")
})
