# Phase 5A batteries.
# polish: run_mxlogit(polish = "newton"/"bhhh") from a loose optimizer stop
#   (maxeval 6) over a J = 6 cross-section with the outside option, scale_vars
#   "sd" and "none", a binding upper bound, a never-chosen alternative (its
#   constant held), a weighted panel with generated draws and the sandwich,
#   optim, and max_iter = 0; plus one explicit polish = "none". Class D
#   against main, which has no polish argument; pinned for Phase 5B.
# ucon: run_mnlogit(), run_mxlogit() and run_nestlogit() on data whose
#   constants have no finite maximum-likelihood estimate (an alternative never
#   chosen, or chosen wherever offered; the outside option or the base
#   alternative so), with and without constants, plus a control design. From
#   the commit that adds the warning these entries record it (class W against
#   main, where they have none); every other slot is identical.

p5a_rec <- function(expr) {
  r <- p2_val(expr)
  r$messages <- r$messages[!startsWith(r$messages, "Optimization run time") &
                            !startsWith(r$messages, "Timing stopped at")]
  # the polish's message ends with its run time
  r$messages <- sub(", [0-9]+h:[0-9]+m:[0-9.]+s\n?$", "", r$messages)
  r
}

# N situations over J alternatives; the chosen one cycles through `chosen`.
p5a_ucon_data <- function(N = 40L, J = 4L, chosen = seq_len(J),
                          labels = seq_len(J), seed = 1L) {
  set.seed(seed)
  pick <- chosen[(seq_len(N) - 1L) %% length(chosen) + 1L]
  dt <- data.table::data.table(id = rep(seq_len(N), each = J),
                               code = rep(seq_len(J), N))
  dt[, choice := as.integer(code == pick[id])]
  dt[, `:=`(alt = labels[code], x1 = rnorm(.N), w1 = rnorm(.N),
            nest = ifelse(code <= J %/% 2L, "A", "B"))]
  dt[, code := NULL]
  dt[]
}

run_ucon <- function() {
  out <- list()
  d_never <- p5a_ucon_data(chosen = 1:3)
  d_ok <- p5a_ucon_data(seed = 2L)
  d_always <- p5a_ucon_data(seed = 3L)
  d_always[id <= 3L, choice := 0L]
  d_always <- rbind(d_always, data.table::data.table(
    id = 1:3, choice = 1L, alt = 5L, x1 = c(0.1, -0.2, 0.3),
    w1 = c(0.4, 0.5, -0.6), nest = "B"))
  data.table::setorderv(d_always, c("id", "alt"))
  d_base_never <- p5a_ucon_data(chosen = 2:4, seed = 4L)
  d_base_always <- p5a_ucon_data(chosen = 2:4, seed = 5L)
  d_base_always <- d_base_always[!(alt == 1L & id > 5L)]
  d_base_always[id <= 5L, choice := as.integer(alt == 1L)]
  d_many <- p5a_ucon_data(N = 60L, J = 9L, chosen = 1:2, seed = 6L)
  d_text <- p5a_ucon_data(chosen = 1:3, labels = c("bus", "car", "rail", "walk"),
                          seed = 7L)
  d_out_none <- p5a_ucon_data(seed = 8L)
  d_out_none[, choice := 0L]
  designs <- list(never = list(d_never), ok = list(d_ok),
                  always = list(d_always), base_never = list(d_base_never),
                  base_always = list(d_base_always), many = list(d_many),
                  text = list(d_text),
                  outside_never = list(p5a_ucon_data(seed = 9L),
                                       include_outside_option = TRUE),
                  outside_always = list(d_out_none, include_outside_option = TRUE))
  fit_mnl <- function(dt, ...) run_mnlogit(dt, "id", "alt", "choice", "x1", ...)
  fit_mxl <- function(dt, ...) {
    run_mxlogit(data = dt, id_col = "id", alt_col = "alt", choice_col = "choice",
                covariate_cols = "x1", random_var_cols = "w1", S = 10L, ...)
  }
  fit_nl <- function(dt, ...) {
    run_nestlogit(data = dt, id_col = "id", alt_col = "alt", choice_col = "choice",
                  covariate_cols = "x1", nest_col = "nest", ...)
  }
  models <- list(mnl = fit_mnl, mxl = fit_mxl, nl = fit_nl)
  for (nm in names(designs)) {
    dz <- designs[[nm]]
    for (md in names(models)) {
      # the nested logit with the outside option is not in this battery
      if (md == "nl" && !is.null(dz$include_outside_option)) next
      for (ua in c(TRUE, FALSE)) {
        args <- c(dz, list(use_asc = ua))
        if (md == "mxl" && isTRUE(dz$include_outside_option)) {
          args$outside_opt_label <- 0L
        }
        out[[paste("ucon", md, nm, if (ua) "asc1" else "asc0", sep = "|")]] <-
          p5a_rec(do.call(models[[md]], args))
      }
    }
  }
  out
}

run_polish_bat <- function() {
  out <- list()
  sim <- function(seed, N = 300L, J = 6L, T = 1L) {
    d <- data.table::as.data.table(simulate_mxl_data(
      N = N, J = J, T = T, beta = c(0.8, -0.6),
      delta = seq(-0.5, 0.5, length.out = J), Sigma = diag(c(1, 0.5)),
      outside_option = TRUE, seed = seed)$data)
    d[alt != 0L]
  }
  base <- function(dt) {
    list(data = dt, id_col = "id", alt_col = "alt", choice_col = "choice",
         covariate_cols = c("x1", "x2"), random_var_cols = c("w1", "w2"),
         S = 20L, include_outside_option = TRUE, outside_opt_label = 0L,
         control = list(maxeval = 6L))
  }
  dt_x <- sim(71L)
  dt_n <- data.table::copy(dt_x)
  dt_n[alt == 6L, choice := 0L]
  dt_p <- sim(72L, N = 100L, T = 3L)
  set.seed(73)
  wt <- stats::runif(100L, 0.5, 2)
  dt_p[, w := wt[pid]]
  b_x <- base(dt_x)
  designs <- list(
    cross = b_x,
    sd = c(b_x, list(scale_vars = "sd")),
    none = c(b_x, list(scale_vars = "none")),
    bound = c(b_x, list(upper = c(x1 = 0.5))),
    held = base(dt_n),
    panel_gen_w = c(base(dt_p), list(person_col = "pid", draws = "generate",
                                     seed = 5L, weights_col = "w",
                                     se_method = "sandwich")),
    optim = c(b_x[names(b_x) != "control"],
              list(optimizer = "optim", control = list(maxit = 3L))),
    cert0 = c(b_x, list(polish_control = list(max_iter = 0L)))
  )
  for (nm in names(designs)) {
    for (pm in c("newton", "bhhh")) {
      out[[paste("polish", nm, pm, sep = "|")]] <-
        p5a_rec(do.call("run_mxlogit", c(designs[[nm]], list(polish = pm))))
    }
  }
  out[["polish|cross|none"]] <-
    p5a_rec(do.call("run_mxlogit", c(b_x, list(polish = "none"))))
  out
}
