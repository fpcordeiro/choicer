# Predicted shares, one per alternative of the fit (after 6D follow-up 6):
# predict(type = "shares") of MNL, NL and MXL fits with and without constants,
# with and without the outside option, on newdata offering every alternative
# or lacking the last-coded, the last two, a middle or the first one, through
# the data.frame and list paths; probabilities on the same newdata; and
# in-sample shares and probabilities.
#
# Deliberate against the parent commit (40bbe86), which numbered the
# alternatives up to the largest code present: keys "shares|...|exit_3|..."
# and "shares|...|exit_2_3|..." for MNL and NL fits (with or without
# constants) and MXL fits without constants ("...|mxl|noasc|..."). Every
# other key is class I, in-sample shares included.

shr_rec <- function(expr) {
  r <- p2_val(expr)
  r$messages <- r$messages[!startsWith(r$messages, "Optimization run time") &
                            !startsWith(r$messages, "Timing stopped at")]
  r
}

run_pred_shares <- function() {
  out <- list()
  for (oo in c(FALSE, TRUE)) {
    sim <- simulate_mxl_data(N = 150L, J = 3L, beta = c(0.8, -0.6),
                             Sigma = diag(c(1, 0.5)),
                             delta = c(0.3, -0.2, 0.1), outside_option = oo,
                             vary_choice_set = FALSE, seed = 41L)
    dt <- data.table::as.data.table(sim$data)
    if (oo) dt <- dt[alt != 0L]
    dt[, nest := ifelse(alt <= 2L, 1L, 2L)]
    args <- list(data = dt, id_col = "id", alt_col = "alt",
                 choice_col = "choice", covariate_cols = c("x1", "x2"))
    if (oo) {
      args <- c(args, list(include_outside_option = TRUE,
                           outside_opt_label = 0L))
    }
    o <- if (oo) "oo" else "nooo"
    offered <- list(all = 1:3, exit_3 = 1:2, exit_2_3 = 1L,
                    exit_2 = c(1L, 3L), exit_1 = 2:3)
    for (ua in c(TRUE, FALSE)) {
      a <- c(args, list(use_asc = ua))
      u <- if (ua) "asc" else "noasc"
      fits <- list(
        mnl = suppressMessages(do.call("run_mnlogit", a)),
        nl = suppressMessages(do.call("run_nestlogit",
                                      c(a, list(nest_col = "nest")))),
        mxl = suppressMessages(do.call("run_mxlogit", c(a, list(
          random_var_cols = c("w1", "w2"), S = 10L)))),
        mxl_gen = suppressMessages(do.call("run_mxlogit", c(a, list(
          random_var_cols = c("w1", "w2"), S = 10L, draws = "generate",
          seed = 9L)))))
      for (m in names(fits)) {
        fit <- fits[[m]]
        key <- function(...) paste("shares", o, m, u, ..., sep = "|")
        out[[key("insample", "shares")]] <- shr_rec(predict(fit, type = "shares"))
        out[[key("insample", "probs")]] <- shr_rec(predict(fit))
        for (case in names(offered)) {
          nd <- dt[alt %in% offered[[case]]]
          out[[key(case, "df")]] <- shr_rec(
            predict(fit, newdata = nd, type = "shares"))
          out[[key(case, "probs")]] <- shr_rec(predict(fit, newdata = nd))
          d <- choicer:::prepare_newdata(fit, nd)
          out[[key(case, "list")]] <- shr_rec(predict(
            fit, type = "shares",
            newdata = list(X = d$X, W = d$W, alt_idx = d$alt_idx, M = d$M)))
        }
      }
    }
  }
  out
}
