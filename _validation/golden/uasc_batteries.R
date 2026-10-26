# 6D follow-up 6 (run_mxlogit(use_asc = FALSE)): mixed logit fits without
# constants and their twins with constants, over the outside option's two
# layouts, panels, generate mode, rc_mean/correlation/log-normal, the advanced
# workflow, se_method, bounds and scale_vars; the post-estimation methods on
# each; and blp() on multinomial and nested logit fits with and without
# constants and the outside option.
#
# Deliberate differences (class D) against main, where every use_asc = FALSE
# mixed logit fit stops ("Theta vector too long") and blp() on a multinomial
# or nested logit fit without constants and with the outside option stops:
#   * every "uasc" key whose config has the field "asc0" and starts with
#     "fit|" or "post|" (the mixed logit without constants);
#   * "blp|mnl|oo|asc0" and "blp|nl|oo|asc0".
# Every other entry (the asc1 twins, MNL/NL fits, blp() with constants or
# without the outside option) is class I.

uasc_rec <- function(expr) {
  r <- p2_val(expr)
  r$messages <- r$messages[!startsWith(r$messages, "Optimization run time") &
                            !startsWith(r$messages, "Timing stopped at")]
  r
}

# A fit's predicted shares, symmetrized: the targets of the BLP leaves.
uasc_sym_shares <- function(fit) {
  tg <- predict(fit, type = "shares")
  (tg + rev(tg)) / 2
}

uasc_sim <- function(seed, oo = FALSE, T = 1L, N = 150L, ...) {
  data.table::as.data.table(simulate_mxl_data(
    N = N, J = 3L, T = T, beta = c(0.8, -0.6), delta = c(0.4, -0.3, 0.2),
    outside_option = oo, seed = seed, ...)$data)
}

run_uasc <- function() {
  out <- list()
  dt_x <- uasc_sim(61L, Sigma = diag(c(1, 0.5)))
  dt_o <- uasc_sim(62L, oo = TRUE, Sigma = diag(c(1, 0.5)))
  dt_i <- dt_o[alt != 0L]
  dt_p <- uasc_sim(63L, T = 3L, N = 60L, Sigma = diag(c(1, 0.5)))
  dt_l <- uasc_sim(64L, Sigma = matrix(c(1, 0.3, 0.3, 0.5), 2),
                   mu = c(0.5, -1), rc_dist = c(0L, 1L))
  dt_x[, cl := (id - 1L) %/% 5L]
  base <- function(dt) {
    list(data = dt, id_col = "id", alt_col = "alt", choice_col = "choice",
         covariate_cols = c("x1", "x2"), random_var_cols = c("w1", "w2"),
         S = 20L)
  }
  oo <- list(include_outside_option = TRUE, outside_opt_label = 0L)
  designs <- list(
    cross = list(args = base(dt_x), sv = c("bhhh", "none", "sd")),
    outside = list(args = c(base(dt_i), oo), sv = c("bhhh", "none")),
    base_row = list(args = c(base(dt_o), list(outside_opt_label = 0L)),
                    sv = c("bhhh", "none")),
    panel = list(args = c(base(dt_p), list(person_col = "pid")),
                 sv = c("bhhh", "none")),
    generate = list(args = c(base(dt_i), oo,
                             list(draws = "generate", seed = 5L)),
                    sv = c("bhhh", "none")),
    corr_ln = list(args = c(base(dt_l), list(rc_mean = TRUE,
                                             rc_correlation = TRUE,
                                             rc_dist = c(0L, 1L))),
                   sv = c("bhhh", "none")),
    se_bhhh = list(args = c(base(dt_x), list(se_method = "bhhh")),
                   sv = "bhhh"),
    se_sandwich = list(args = c(base(dt_i), oo, list(se_method = "sandwich")),
                       sv = "bhhh"),
    se_cluster = list(args = c(base(dt_x), list(cluster_col = "cl")),
                      sv = "bhhh"),
    bound_L = list(args = c(base(dt_x), list(lower = c(L_11 = -0.5),
                                             upper = c(x1 = 2))),
                   sv = c("bhhh", "none")),
    bound_asc = list(args = c(base(dt_x), list(lower = c(ASC_2 = -1))),
                     sv = "bhhh"),
    theta_init = list(args = c(base(dt_x),
                               list(theta_init = c(0.5, -0.5, 0, 0))),
                      sv = c("bhhh", "none"))
  )
  d_adv <- prepare_mxl_data(dt_i, "id", "alt", "choice", c("x1", "x2"),
                            c("w1", "w2"), include_outside_option = TRUE,
                            outside_opt_label = 0L)
  designs$advanced <- list(
    args = list(input_data = d_adv,
                eta_draws = get_halton_normals(20L, d_adv$N, 2L)),
    sv = c("bhhh", "none"))

  nd_of <- function(a) {
    if (is.null(a$data)) return(NULL)
    nd <- data.table::copy(a$data)
    nd[, x1 := x1 + 0.5]
    nd
  }
  for (nm in names(designs)) {
    dz <- designs[[nm]]
    for (sv in dz$sv) {
      for (ua in c(FALSE, TRUE)) {
        args <- c(dz$args, list(use_asc = ua, scale_vars = sv))
        if (nm == "theta_init" && ua) {
          args$theta_init <- c(args$theta_init, 0, 0)
        }
        tag <- paste(nm, sv, if (ua) "asc1" else "asc0", sep = "|")
        fit <- NULL
        rec <- uasc_rec(fit <- do.call("run_mxlogit", args))
        out[[paste0("fit|", tag)]] <- rec
        if (sv != "bhhh") next
        nd <- nd_of(dz$args)
        f <- function(expr) {
          if (is.null(fit)) return(uasc_rec(stop("no fit")))
          uasc_rec(expr)
        }
        post <- list(
          coef_names = f(names(coef(fit))),
          param_map = f(fit$param_map),
          pred = f(predict(fit)),
          shares = f(predict(fit, type = "shares")),
          shares_nd = f(if (!is.null(nd)) predict(fit, newdata = nd,
                                                  type = "shares")),
          elas_x1 = f(elasticities(fit, "x1")),
          elas_w1 = f(elasticities(fit, "w1", is_random_coef = TRUE)),
          div_x1 = f(diversion_ratios(fit, "x1")),
          div_w2 = f(diversion_ratios(fit, "w2", is_random_coef = TRUE)),
          logsum = f(logsum(fit)),
          cs = f(consumer_surplus(fit, price_var = "x2")),
          wtp = f(wtp(fit, price_var = "x2")),
          tastes = f(conditional_tastes(fit)[c("mean", "sd")]),
          gof = f(gof(fit)),
          summary = f(summary(fit)$coefficients),
          v_hessian = f(vcov(fit, type = "hessian")),
          v_bhhh = f(vcov(fit, type = "bhhh")),
          v_robust = f(vcov(fit, type = "robust")),
          wesml = f(wesml_vcov(fit, type = "se")),
          # (the targets formed inside the recorder, which records their
          # errors and warnings as the call's)
          blp = golden_blp("blp", list(fit, target_shares = uasc_sym_shares(fit)))
        )
        for (p in names(post)) out[[paste0("post|", tag, "|", p)]] <- post[[p]]
      }
    }
  }

  # blp() on multinomial and nested logit fits
  s <- data.table::as.data.table(simulate_mnl_data(
    N = 200L, J = 3L, beta = c(0.8, -0.6), delta = c(0.4, -0.3, 0.2),
    outside_option = TRUE, seed = 65L)$data)
  data_mn <- list(oo = s[alt != 0L], nooo = s)
  for (o in names(data_mn)) {
    d <- data_mn[[o]]
    d[, nest := ifelse(alt <= 1L, 1L, 2L)]
    oo_args <- if (o == "oo") {
      list(include_outside_option = TRUE, outside_opt_label = 0L)
    } else {
      list(outside_opt_label = 0L)
    }
    for (ua in c(FALSE, TRUE)) {
      a <- if (ua) "asc1" else "asc0"
      mnl <- suppressMessages(do.call("run_mnlogit", c(
        list(data = d, id_col = "id", alt_col = "alt", choice_col = "choice",
             covariate_cols = c("x1", "x2"), use_asc = ua), oo_args)))
      nl <- suppressMessages(do.call("run_nestlogit", c(
        list(data = d, id_col = "id", alt_col = "alt", choice_col = "choice",
             covariate_cols = c("x1", "x2"), nest_col = "nest", use_asc = ua),
        oo_args)))
      for (m in c("mnl", "nl")) {
        fit <- if (m == "mnl") mnl else nl
        out[[paste("blp", m, o, a, sep = "|")]] <-
          golden_blp("blp", list(fit, target_shares = uasc_sym_shares(fit)))
        out[[paste("blp", m, o, a, "init", sep = "|")]] <-
          golden_blp("blp", list(fit, target_shares = uasc_sym_shares(fit),
                                 delta_init = rep(0, nrow(fit$alt_mapping) - (o == "oo"))))
      }
    }
  }
  out
}
