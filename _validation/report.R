# ============================================================================
# The report of a run: <run>/REPORT.md, <run>/tables/*.csv, <run>/plots/*.png.
# For each family in the run, against the tier's design (its units, reps and
# gate cells): the units' rep counts, units absent from the run and missing
# reps flagged; the metrics of every (cell, method) (metrics_cell(), then
# gate_family() over the family; an absent unit of a gate cell gets its gate
# entries as "na" rows), a table per cell and method with the gated rows (H)
# first and the findings below; failures (distinct errors with counts),
# warnings by class and messages; plots in base graphics: a QQ plot of the
# Hessian z at the gate cells, ladders (ESE*sqrt(N), N*bias, tail imbalance)
# against every numeric cell column that varies, with rate_fit() slopes,
# and SBC ECDF differences where the records carry ranks; then the run's
# metadata. Reads one (cell, method)'s records at a time. Runs with R = 2
# per cell (NA metrics) and with families whose reps all failed. A section
# that errors is reported in REPORT.md and makes the exit status 1, as does
# a run without rep files.
#
# Usage (the build under test is the FIRST entry of R_LIBS):
#   R_LIBS=<lib> OMP_NUM_THREADS=1 Rscript _validation/report.R --run=DIR [--tier=mc|quick]
#   --run   a run directory, or a run tag under _validation/output/
#   --tier  mc: tolerances "mc", the gate cells gates$cell, the reps of
#           R_by_method or R$mc; quick (Q3's rep files): tolerances "quick",
#           the units quick.R ran (q3_plan.rds; else the cells of R$quick at
#           every method), all of them gated, as quick.R gates them
# ============================================================================

root <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE),
                                                       value = TRUE)[1])))
sys.source(file.path(root, "lib", "common.R"), envir = (boot <- new.env()), keep.source = FALSE)
b <- boot$vs_boot(root)
build <- b$build
L <- b$L

args <- vs_args(defaults = list(run = NULL, tier = "mc"))
if (is.null(args$run)) stop("--run=DIR is required", call. = FALSE)
if (!args$tier %in% c("mc", "quick")) stop("--tier must be mc or quick", call. = FALSE)
run_dir <- if (dir.exists(args$run)) normalizePath(args$run) else
  file.path(vs_root(), "output", args$run)
if (!dir.exists(run_dir)) stop("no run directory ", args$run, call. = FALSE)
tab_dir <- vs_out_dir(out = file.path(run_dir, "tables"))
plot_dir <- vs_out_dir(out = file.path(run_dir, "plots"))
# The units quick.R ran (Q3's rep files): data.frame(family, cell, method, R).
q3_plan <- if (args$tier == "quick" && file.exists(file.path(run_dir, "q3_plan.rds"))) {
  readRDS(file.path(run_dir, "q3_plan.rds"))
}

MAX_BLOCK <- 10L    # blocks with more estimands are summarized in REPORT.md
MAX_PANELS <- 12L   # estimands per QQ or ladder figure
errors <- character()

# --- Helpers ------------------------------------------------------------------

# Evaluates `expr`; an error is recorded (and returns NULL) so that one
# section cannot stop the report.
rp_try <- function(what, expr) {
  tryCatch(expr, error = function(e) {
    errors <<- c(errors, sprintf("%s: %s", what, conditionMessage(e)))
    NULL
  })
}

.fmt <- function(x, digits = 4L) {
  out <- formatC(as.numeric(x), digits = digits, format = "g")
  out[is.na(x)] <- "NA"
  trimws(out)
}
.pm <- function(v, s) ifelse(is.na(v), "NA", paste0(.fmt(v), " ± ", .fmt(s, 2L)))
.short <- function(x, n = 160L) {
  x <- gsub("[\r\n\t]+", " ", x)
  ifelse(nchar(x) > n, paste0(substr(x, 1L, n - 3L), "..."), x)
}

# A data.frame as a Markdown table (cells as text, pipes escaped).
md_table <- function(df) {
  if (is.null(df) || !nrow(df)) return("(none)")
  cells <- vapply(df, function(x) gsub("|", "\\|", ifelse(is.na(x), "NA", as.character(x)),
                                       fixed = TRUE), character(nrow(df)))
  if (!is.matrix(cells)) cells <- matrix(cells, nrow = 1L)
  c(paste0("| ", paste(names(df), collapse = " | "), " |"),
    paste0("|", paste(rep("---", ncol(df)), collapse = "|"), "|"),
    paste0("| ", apply(cells, 1L, paste, collapse = " | "), " |"))
}

# rbind of data.frames whose columns differ (missing columns NA).
.rbind_fill <- function(dfs) {
  dfs <- Filter(function(d) !is.null(d) && nrow(d), dfs)
  if (!length(dfs)) return(NULL)
  cols <- unique(unlist(lapply(dfs, names)))
  do.call(rbind, lapply(dfs, function(d) {
    for (k in setdiff(cols, names(d))) d[[k]] <- NA
    d[cols]
  }))
}

# Reps as runs ("1-40,45,47-50"); more than six runs are cut to the first
# three and the last two, with the count of reps.
.ranges <- function(r) {
  if (!length(r)) return("")
  r <- sort(unique(r))
  br <- c(0L, which(diff(r) != 1L), length(r))
  s <- vapply(seq_len(length(br) - 1L), function(k) {
    a <- r[br[k] + 1L]; b <- r[br[k + 1L]]
    if (a == b) as.character(a) else sprintf("%d-%d", a, b)
  }, "")
  if (length(s) <= 6L) return(paste(s, collapse = ","))
  sprintf("%s,...,%s (%d reps)", paste(s[1:3], collapse = ","),
          paste(utils::tail(s, 2L), collapse = ","), length(r))
}

.png <- function(file, width, height) {
  dev <- list(filename = file, width = width, height = height, res = 96)
  if (isTRUE(capabilities("cairo"))) dev$type <- "cairo"
  do.call(grDevices::png, dev)
}

# --- Reading a (cell, method) -------------------------------------------------

# One (cell, method)'s summary from its records: counts, failures,
# warnings, messages, the metrics (with pooled_tails() rows for blocks of
# more than MAX_BLOCK estimands), the Hessian z of each estimand (gate
# cells only) and SBC ranks.
#
# SBC ranks are the records' diag$sbc, list(L = <number of posterior
# draws>, ranks = c(<test quantity> = <rank in 0..L>)), as metrics_cell()
# reads them; records without them have no SBC section.
rp_unit <- function(fam_id, cell, method, gates, want_z) {
  recs <- withCallingHandlers(read_reps(run_dir, fam_id, cell, method),
                              vs_unreadable_rep = function(w) invokeRestart("muffleWarning"))
  ok <- vapply(recs, function(r) isTRUE(r$ok), logical(1))
  reps <- vapply(recs, function(r) as.integer(r$rep), integer(1))
  err <- vapply(recs[!ok], function(r) as.character(r$error %||% "(no message)")[1], "")
  fails <- if (length(err)) {
    t <- table(err)
    data.frame(cell = cell, method = method, error = .short(names(t), 300L),
               count = as.integer(t), stringsAsFactors = FALSE)
  }
  w <- .rbind_fill(lapply(recs, function(r) {
    if (!is.data.frame(r$warnings) || !nrow(r$warnings)) return(NULL)
    data.frame(rep = r$rep, r$warnings[c("stage", "class", "message")], stringsAsFactors = FALSE)
  }))
  warns <- if (!is.null(w)) {
    key <- paste(w$stage, w$class, sep = "\r")
    do.call(rbind, lapply(split(seq_len(nrow(w)), factor(key, unique(key))), function(i) {
      data.frame(cell = cell, method = method, stage = w$stage[i[1]], class = w$class[i[1]],
                 count = length(i), reps = length(unique(w$rep[i])),
                 example = .short(w$message[i[1]]), stringsAsFactors = FALSE)
    }))
  }
  m <- unlist(lapply(recs, function(r) r$diag$messages))
  msgs <- if (length(m)) {
    s <- tapply(as.integer(m), factor(names(m), unique(names(m))), sum)
    nr <- tapply(names(m), factor(names(m), unique(names(m))), length)
    data.frame(cell = cell, method = method, message = .short(names(s)), count = as.integer(s),
               reps = as.integer(nr), stringsAsFactors = FALSE)
  }
  blocks <- character()
  z <- list()
  sbc <- list()
  for (r in recs[ok]) {
    e <- r$est
    if (is.data.frame(e) && nrow(e)) {
      blocks[e$name] <- e$block
      if (want_z) z <- c(z, list(data.frame(name = e$name, block = e$block,
                                            z = (e$est - e$truth) / e$se_hessian)))
    }
    d <- r$derived
    if (is.data.frame(d) && nrow(d)) {
      blocks[d$name] <- "derived"
      if (want_z) z <- c(z, list(data.frame(name = d$name, block = "derived",
                                            z = (d$est - d$truth) / d$se)))
    }
    s <- r$diag$sbc
    if (!is.null(s$ranks)) sbc[[length(sbc) + 1L]] <- s
  }
  met <- if (length(recs)) {
    rp_try(sprintf("metrics_cell(%s, %s, %s)", fam_id, cell, method),
           metrics_cell(recs, gates = gates, tol_set = args$tier))
  }
  # Blocks of more than MAX_BLOCK estimands (thousands of constants): their
  # coverage pooled within reps, rows that bind with metrics_cell()'s.
  eb <- table(blocks[blocks != "derived"])
  for (bk in names(eb)[eb > MAX_BLOCK]) {
    met <- .rbind_fill(list(met, rp_try(sprintf("pooled_tails(%s, %s, %s, %s)", fam_id, cell,
                                                method, bk), pooled_tails(recs, block = bk))))
  }
  bad <- length(attr(recs, "unreadable"))
  list(files = length(recs) + bad, ok = sum(ok), failed = sum(!ok), unreadable = bad,
       reps = reps, fails = fails, warns = warns, msgs = msgs, met = met,
       z = if (length(z)) do.call(rbind, z), sbc = sbc)
}

# A unit of the design without rep files: nothing read, and at a gate cell
# its gate entries as "na" rows. metrics_cell() adds such a row for every
# gate entry that matches no row of a gated cell, so they are taken from it
# on one placeholder failed record, keeping the gated rows (R = 0).
rp_absent <- function(fam_id, cell, method, gates) {
  met <- if (cell %in% gates$cell) {
    ph <- new_rep_record(fam_id, cell, method, 1L, NA_integer_, ok = FALSE,
                         error = "absent from the run")
    tab <- rp_try(sprintf("metrics_cell(%s, %s, %s) (absent)", fam_id, cell, method),
                  metrics_cell(list(ph), gates = gates, tol_set = args$tier))
    tab <- tab[tab$gated, , drop = FALSE]
    tab$R <- rep(0L, nrow(tab))
    tab$n_failed <- rep(0L, nrow(tab))
    tab
  }
  list(files = 0L, ok = 0L, failed = 0L, unreadable = 0L, reps = integer(), met = met,
       sbc = list())
}

# --- Tables -------------------------------------------------------------------

# The blocks of a metrics table with more than MAX_BLOCK estimands (rows of
# the cell's metrics and tests have no block; pooled_tails() rows are not
# estimands).
rp_big_blocks <- function(tab) {
  eb <- unique(tab[!is.na(tab$block) & !rp_pooled(tab), c("estimand", "block"), drop = FALSE])
  k <- table(eb$block)
  names(k)[k > MAX_BLOCK]
}

rp_pooled <- function(tab) startsWith(tab$metric, "pooled_")

.blank <- function(x) ifelse(is.na(x), "", x)

# Rows in long form: metric, estimand, route, value ± MCSE, R (the reps
# entering the row), target, tol, Holm p and verdict.
rp_long <- function(d) {
  data.frame(metric = d$metric, estimand = .blank(d$estimand), route = .blank(d$route),
             value = .pm(d$value, d$mcse), R = d$R, target = .fmt(d$target), tol = .fmt(d$tol),
             p_holm = .fmt(d$p_holm, 2L), verdict = d$verdict, stringsAsFactors = FALSE)
}
RP_FINDING_COLS <- c("metric", "estimand", "route", "value", "R", "target", "tol")

# A summarized block: one row per (metric, route) with the mean value over
# its estimands, their range, the range of R and the verdict counts.
rp_block_rows <- function(d, b) {
  key <- paste(d$metric, d$route, sep = "\r")
  do.call(rbind, lapply(split(seq_len(nrow(d)), factor(key, unique(key))), function(i) {
    v <- table(d$verdict[i], useNA = "no")
    rng <- suppressWarnings(range(d$value[i], na.rm = TRUE))
    rr <- range(d$R[i])
    data.frame(metric = d$metric[i[1]], estimand = sprintf("%s (%d)", b, length(unique(d$estimand[i]))),
               route = .blank(d$route[i[1]]),
               value = sprintf("mean %s [%s, %s]", .fmt(mean(d$value[i], na.rm = TRUE)),
                               .fmt(rng[1]), .fmt(rng[2])),
               R = if (rr[1] == rr[2]) as.character(rr[1]) else paste(rr, collapse = "-"),
               target = .fmt(d$target[i[1]]), tol = .fmt(d$tol[i[1]]), p_holm = "",
               verdict = paste(names(v), v, collapse = ", "), stringsAsFactors = FALSE)
  }))
}

# The metrics of the estimands' pivot.
RP_PIVOT <- c("bias_std", "ese", "se_ratio", "tail_lower", "tail_upper", "coverage")

# The estimands' findings pivoted: a row per (estimand, route), a column per
# metric (bias_std and ese on the estimand's first row), "H" where the entry
# is gated (its value is in the gated table).
rp_pivot <- function(d) {
  ms <- RP_PIVOT
  out <- list()
  for (e in unique(d$estimand)) {
    de <- d[d$estimand == e, , drop = FALSE]
    routes <- unique(de$route[!is.na(de$route)])
    if (!length(routes)) routes <- NA_character_
    for (k in seq_along(routes)) {
      cell <- function(m) {
        x <- de[de$metric == m & (if (m %in% ms[1:2]) is.na(de$route) & k == 1L
                                  else de$route %in% routes[k]), , drop = FALSE]
        if (!nrow(x)) "" else if (x$gated[1]) "H" else .pm(x$value[1], x$mcse[1])
      }
      out[[length(out) + 1L]] <- data.frame(
        estimand = if (k == 1L) e else "", route = .blank(routes[k]),
        t(vapply(ms, cell, "")), stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, out)
}

# A block's worst rows: per metric with a target, the k rows farthest from
# it in MCSE units.
rp_worst <- function(d, k = 3L) {
  d <- d[!is.na(d$target) & is.finite(d$value) & is.finite(d$mcse) & d$mcse > 0, , drop = FALSE]
  if (!nrow(d)) return(NULL)
  d$z <- abs(d$value - d$target) / d$mcse
  d <- do.call(rbind, lapply(split(d, factor(d$metric, unique(d$metric))), function(x) {
    utils::head(x[order(-x$z), , drop = FALSE], k)
  }))
  data.frame(rp_long(d)[c("metric", "estimand", "route", "value", "R", "target")],
             `|dev|/MCSE` = .fmt(d$z, 3L), verdict = d$verdict, check.names = FALSE,
             stringsAsFactors = FALSE)
}

# Gated rows for the family's gate section (at most `max` of them).
rp_gate_rows <- function(b, max = 50L) {
  out <- md_table(utils::head(data.frame(
    cell = b$cell, method = b$method, metric = b$metric, estimand = .blank(b$estimand),
    route = .blank(b$route), value = .pm(b$value, b$mcse), R = b$R, target = .fmt(b$target),
    tol = .fmt(b$tol), p_holm = .fmt(b$p_holm, 2L), verdict = b$verdict), max))
  if (nrow(b) > max) out <- c(out, "", sprintf("... and %d more in the metrics table.", nrow(b) - max))
  out
}

# One (cell, method)'s metrics as Markdown: the gated rows (H) with target,
# tolerance, Holm p and verdict; then the findings: the estimands pivoted
# (targets and tolerances in the line above), the cell's rows in long form;
# then each block of more than MAX_BLOCK estimands (thousands of constants)
# summarized: its coverage pooled within reps (pooled_tails()), one row per
# (metric, route) and its worst rows. Every row is in tables/<family>_metrics.csv.
rp_metrics_md <- function(tab) {
  if (is.null(tab) || !nrow(tab)) return("(no metrics)")
  big <- rp_big_blocks(tab)
  pooled <- rp_pooled(tab)
  inbig <- !is.na(tab$block) & tab$block %in% big & !pooled
  out <- character()
  g <- tab[tab$gated & !inbig, , drop = FALSE]
  gb <- lapply(big, function(b) {
    d <- tab[tab$gated & inbig & tab$block %in% b, , drop = FALSE]
    if (nrow(d)) rp_block_rows(d, b)
  })
  if (nrow(g) || length(unlist(gb))) {
    out <- c(out, "Gated (H):", "", md_table(rbind(rp_long(g), do.call(rbind, gb))), "")
  }
  bad <- tab[tab$gated & inbig & tab$verdict %in% c("fail", "na"), , drop = FALSE]
  if (nrow(bad)) out <- c(out, "**Failing or NA gated rows of summarized blocks:**", "",
                          rp_gate_rows(bad), "")
  est <- tab[!is.na(tab$block) & !inbig & !pooled & tab$metric %in% RP_PIVOT, , drop = FALSE]
  if (nrow(est) && !all(est$gated)) {
    tg <- unique(est[!est$metric %in% "ese", c("metric", "target", "tol"), drop = FALSE])
    out <- c(out, sprintf("Findings by estimand (target, tol: %s; R in the metrics table):",
                          paste(sprintf("%s %s, %s", tg$metric, .fmt(tg$target), .fmt(tg$tol)),
                                collapse = "; ")), "", md_table(rp_pivot(est)), "")
  }
  # The cell's rows, and estimand rows of metrics the pivot does not show
  # (the SBC section lists the sbc_ecdf rows).
  cl <- tab[!tab$gated & tab$metric != "sbc_ecdf" &
              (is.na(tab$block) | (!inbig & !pooled & !tab$metric %in% RP_PIVOT)), , drop = FALSE]
  if (nrow(cl)) out <- c(out, "Findings of the cell:", "", md_table(rp_long(cl)[RP_FINDING_COLS]), "")
  for (b in big) {
    d <- tab[inbig & tab$block %in% b, , drop = FALSE]
    out <- c(out, sprintf("Block %s, %d estimands (every row in the metrics table):", b,
                          length(unique(d$estimand))), "")
    p <- tab[pooled & tab$block %in% b, , drop = FALSE]
    if (nrow(p)) out <- c(out, md_table(rp_long(p)[RP_FINDING_COLS]), "")
    f <- d[!d$gated, , drop = FALSE]
    if (nrow(f)) out <- c(out, md_table(rp_block_rows(f, b)), "")
    w <- rp_worst(d)
    if (!is.null(w)) out <- c(out, "Its worst rows:", "", md_table(w), "")
  }
  out
}

# --- Plots --------------------------------------------------------------------

# QQ plots of z = (estimate - truth) / se_hessian, one panel per estimand
# of the small blocks (at most MAX_PANELS) and one per summarized block,
# pooled over its estimands.
rp_qq <- function(z, file, title) {
  if (is.null(z) || !any(is.finite(z$z))) return(NULL)
  k <- tapply(z$name, z$block, function(x) length(unique(x)))
  big <- names(k)[k > MAX_BLOCK]
  small <- unique(z$name[!z$block %in% big])
  panels <- c(stats::setNames(lapply(utils::head(small, MAX_PANELS),
                                     function(n) z$z[z$name == n]), utils::head(small, MAX_PANELS)),
              stats::setNames(lapply(big, function(b) z$z[z$block == b]),
                              sprintf("%s (pooled, %d)", big, k[big])))
  panels <- Filter(function(v) sum(is.finite(v)) >= 2L, panels)
  if (!length(panels)) return(NULL)
  nc <- min(4L, length(panels))
  nr <- ceiling(length(panels) / nc)
  .png(file, 260L * nc, 260L * nr + 30L)
  on.exit(grDevices::dev.off())
  graphics::par(mfrow = c(nr, nc), mar = c(4, 4, 2.5, 1), oma = c(0, 0, 2, 0))
  for (n in names(panels)) {
    v <- panels[[n]][is.finite(panels[[n]])]
    lim <- range(v, stats::qnorm(stats::ppoints(length(v))), -2, 2)
    stats::qqnorm(v, main = n, cex.main = 0.9, pch = 20, col = grDevices::adjustcolor(1, 0.6),
                  xlim = lim, ylim = lim)
    graphics::abline(0, 1, col = 2)
  }
  graphics::mtext(title, outer = TRUE, cex = 0.9)
  basename(file)
}

# A method's ladder rows: for each cell and estimand, ESE*sqrt(N), N*bias
# (bias = bias_std * ESE; MCSE N * ESE / sqrt(R)) and the Hessian tail
# imbalance upper - lower, its MCSE from the multinomial, Var = (pu + pl -
# (pu - pl)^2) / R, at the smoothed shares p = (x + 1) / (R + 2) of x reps
# out of R, as metrics_cell() evaluates a proportion's MCSE.
rp_ladder_rows <- function(tab, N) {
  rows_of <- function(metric, route = NA) {
    d <- tab[tab$metric == metric & (if (is.na(route)) is.na(tab$route) else tab$route %in% route),
             c("cell", "estimand", "value", "mcse", "R", "target"), drop = FALSE]
    d$key <- paste(d$cell, d$estimand, sep = "\r")
    d
  }
  ladder <- function(d, quantity, value, mcse) {
    if (!nrow(d)) return(NULL)
    data.frame(cell = d$cell, estimand = d$estimand, quantity = quantity, value = value,
               mcse = mcse, stringsAsFactors = FALSE)
  }
  ese <- rows_of("ese")
  b <- rows_of("bias_std")[match(ese$key, rows_of("bias_std")$key), , drop = FALSE]
  lo <- rows_of("tail_lower", "hessian")
  up <- rows_of("tail_upper", "hessian")[match(lo$key, rows_of("tail_upper", "hessian")$key), ,
                                         drop = FALSE]
  Nc <- N[ese$cell]
  R <- lo$R
  pl <- (lo$value * R + 1) / (R + 2)
  pu <- (up$value * R + 1) / (R + 2)
  .rbind_fill(list(
    ladder(ese, "ESE*sqrt(N)", ese$value * sqrt(Nc), ese$mcse * sqrt(Nc)),
    ladder(ese, "N*bias", Nc * b$value * ese$value, Nc * ese$value / sqrt(b$R)),
    ladder(lo, "tail upper-lower", up$value - lo$value,
           ifelse(R >= 2, sqrt((pu + pl - (pu - pl)^2) / pmax(R, 1)), NA_real_))))
}

# The ladder figure of one method against one cell column: a row of three
# panels per estimand, log x (and log y for ESE*sqrt(N)), MCSE bars.
rp_ladder_plot <- function(lad, x, xname, file, title) {
  ests <- utils::head(unique(lad$estimand[is.finite(lad$value)]), MAX_PANELS)
  if (!length(ests)) return(NULL)
  qs <- c("ESE*sqrt(N)", "N*bias", "tail upper-lower")
  .png(file, 3L * 280L, 200L * length(ests) + 40L)
  on.exit(grDevices::dev.off())
  graphics::par(mfrow = c(length(ests), 3L), mar = c(4, 4, 2, 1), oma = c(0, 0, 2, 0))
  for (e in ests) for (q in qs) {
    d <- lad[lad$estimand == e & lad$quantity == q & is.finite(lad$value), , drop = FALSE]
    d$x <- x[d$cell]
    d <- d[is.finite(d$x) & d$x > 0, , drop = FALSE]
    if (!nrow(d)) { graphics::plot.new(); graphics::title(paste(e, q), cex.main = 0.9); next }
    m <- ifelse(is.finite(d$mcse), d$mcse, 0)
    lo <- d$value - m; hi <- d$value + m
    logy <- q == "ESE*sqrt(N)" && all(lo > 0)
    graphics::plot(d$x, d$value, log = if (logy) "xy" else "x", pch = 19, xlab = xname,
                   ylab = q, main = paste(e, "-", q), cex.main = 0.9,
                   ylim = range(c(lo, hi, if (!logy) 0)))
    graphics::arrows(d$x, lo, d$x, hi, angle = 90, code = 3, length = 0.03)
    if (!logy) graphics::abline(h = 0, col = "grey50", lty = 2)
  }
  graphics::mtext(title, outer = TRUE, cex = 0.9)
  basename(file)
}

# rate_fit() slopes of a method against one cell column: ese, bias_std and
# the tail imbalance (as metric rows), per estimand. The predicted slope is
# given against N only: -1/2 for all three (C1, C5).
rp_rates <- function(tab, lad, x, xname, method) {
  imb <- lad[lad$quantity == "tail upper-lower", ]
  imb_tab <- data.frame(cell = imb$cell, metric = "tail_imbalance", estimand = imb$estimand,
                        route = "hessian", value = imb$value, mcse = imb$mcse,
                        stringsAsFactors = FALSE)
  both <- .rbind_fill(list(tab[c("cell", "metric", "estimand", "route", "value", "mcse")],
                           imb_tab))
  out <- list()
  for (e in utils::head(unique(lad$estimand), MAX_PANELS)) {
    for (m in c("ese", "bias_std", "tail_imbalance")) {
      route <- if (m == "tail_imbalance") "hessian" else NA
      f <- tryCatch(rate_fit(both, x, m, estimand = e, route = route), error = function(e) NULL)
      if (is.null(f)) next
      out[[length(out) + 1L]] <- data.frame(
        method = method, x = xname, estimand = e, metric = m,
        slope = as.numeric(f[["slope"]]), se = as.numeric(f[["se"]]), n = as.integer(f[["n"]]),
        predicted = if (xname == "N") -0.5 else NA_real_, stringsAsFactors = FALSE)
    }
  }
  .rbind_fill(out)
}

# SBC (ranks from diag$sbc, see rp_unit()): one ECDF-difference panel per
# test quantity, with the panels' sbc_ecdf_test() p-values. The verdicts
# are the metrics table's sbc_ecdf rows (gate_family()); the p-values are
# shown, ungated, only for units whose metrics carry no such rows.
rp_sbc <- function(sbc, file, title) {
  if (!length(sbc)) return(NULL)
  Lr <- sbc[[1]]$L
  qn <- unique(unlist(lapply(sbc, function(s) names(s$ranks))))
  ranks <- lapply(qn, function(q) unlist(lapply(sbc, function(s) s$ranks[[q]])))
  names(ranks) <- qn
  tests <- lapply(ranks, function(r) sbc_ecdf_test(r, Lr))
  p <- vapply(tests, function(t) as.numeric(t$p), numeric(1))
  nc <- min(4L, length(qn))
  nr <- ceiling(length(qn) / nc)
  .png(file, 340L * nc, 290L * nr + 30L)
  on.exit(grDevices::dev.off())
  graphics::par(mfrow = c(nr, nc), mar = c(4, 4, 3.5, 1), oma = c(0, 0, 2, 0), cex.main = 0.9)
  for (q in qn) sbc_ecdf_plot(ranks[[q]], Lr, test = tests[[q]], main = q)
  graphics::mtext(title, outer = TRUE, cex = 0.9)
  list(file = basename(file),
       tab = data.frame(quantity = qn, R = lengths(ranks), L = Lr, p = p, stringsAsFactors = FALSE))
}

# --- A family ------------------------------------------------------------------

# A unit's expected reps: R_by_method[[cell]][method] where the family gives
# it (mc tier), else R[[tier]][cell]; NA when the family names neither.
rp_expected <- function(fam, cell, method) {
  n <- if (args$tier == "mc") fam$R_by_method[[cell]][method]
  if (is.null(n) || is.na(n)) n <- fam$R[[args$tier]][cell]
  if (length(n) != 1L) NA_integer_ else as.integer(n)
}

# The tier's design of a family: its units data.frame(cell, method, R) and
# its gates. mc: every cell's methods at rp_expected() reps, gated at
# gates$cell. quick: the units quick.R ran (q3_plan.rds), else the cells of
# R$quick at every method of the cell, all of them gated, as quick.R gates
# them. The gate metrics are the family's.
rp_design <- function(fam) {
  none <- data.frame(cell = character(), method = character(), R = integer())
  if (is.null(fam)) return(list(units = none, gates = NULL))
  p <- if (!is.null(q3_plan)) q3_plan[q3_plan$family == fam$id, c("cell", "method", "R")]
  if (is.null(p) || !nrow(p)) {
    cells <- if (args$tier == "mc") fam$cells$cell else names(fam$R$quick)[fam$R$quick > 0]
    p <- .rbind_fill(lapply(cells, function(cl) {
      ms <- cell_methods(fam, cl)
      if (length(ms)) {
        data.frame(cell = cl, method = ms, R = vapply(ms, function(m) rp_expected(fam, cl, m), 1L),
                   stringsAsFactors = FALSE, row.names = NULL)
      }
    })) %||% none
    p <- p[!is.na(p$R) & p$R > 0, , drop = FALSE]
  }
  gate_cells <- if (args$tier == "mc") fam$gates$cell else unique(p$cell)
  list(units = p, gates = list(cell = gate_cells, metrics = fam$gates$metrics))
}

# A family's section of REPORT.md (its units, gates, a section per cell and
# method, plots, rates, SBC), with its tables and plots written; returns
# list(md, summary). `idx` lists the family's rep files.
rp_family <- function(fam_id, idx) {
  md <- c(sprintf("## %s", fam_id), "")
  fam <- rp_try(sprintf("load_family(%s)", fam_id), load_family(fam_id, L))
  design <- rp_design(fam)
  gates <- design$gates
  if (is.null(fam)) {
    md <- c(md, "Family file not loaded: no design, gates or cell columns.", "")
  } else {
    gm <- gates$metrics
    gm <- if (is.list(gm)) {
      paste(sprintf("%s: %s", names(gm), vapply(gm, paste, "", collapse = ", ")), collapse = "; ")
    } else {
      paste(gm, collapse = ", ")
    }
    md <- c(md, sprintf("Claims %s; DGP %s; tier %s: gate cell(s) %s; gated metrics %s.",
                        paste(fam$claims, collapse = ", "), fam$dgp, args$tier,
                        paste(gates$cell, collapse = ", "), gm), "")
  }

  # The units: the design's and those with rep files, in the family's cell order.
  cell_order <- unique(c(fam$cells$cell, sort(unique(idx$cell))))
  keys <- unique(rbind(idx[c("cell", "method")], design$units[c("cell", "method")]))
  keys <- keys[order(match(keys$cell, cell_order), keys$method), , drop = FALSE]
  present <- paste(keys$cell, keys$method) %in% paste(idx$cell, idx$method)
  units <- list()
  for (k in seq_len(nrow(keys))) {
    cl <- keys$cell[k]; m <- keys$method[k]
    u <- if (present[k]) {
      rp_try(sprintf("reading %s/%s/%s", fam_id, cl, m),
             rp_unit(fam_id, cl, m, gates, want_z = cl %in% gates$cell))
    } else {
      rp_absent(fam_id, cl, m, gates)
    }
    if (!is.null(u)) units[[paste(cl, m, sep = "\r")]] <- c(list(cell = cl, method = m), u)
  }

  # Rep counts: the reps missing against the design (absent units, gaps and
  # reps missing at the end of the range alike) and the gaps in the files.
  cnt <- .rbind_fill(lapply(units, function(u) {
    i <- which(design$units$cell == u$cell & design$units$method == u$method)
    exp <- if (length(i)) as.integer(design$units$R[i[1]]) else NA_integer_
    miss <- if (is.na(exp)) 0L else length(setdiff(seq_len(exp), u$reps))
    gaps <- if (length(u$reps)) setdiff(seq_len(max(u$reps)), u$reps) else integer()
    data.frame(cell = u$cell, method = u$method,
               status = if (!u$files) "**absent**" else if (miss) "**incomplete**" else "",
               files = u$files, ok = u$ok, failed = u$failed, unreadable = u$unreadable,
               reps = .ranges(u$reps), expected = exp, missing = miss,
               gaps = if (length(gaps)) .ranges(gaps) else "", stringsAsFactors = FALSE)
  }))
  md <- c(md, sprintf("### Units (tier %s)", args$tier), "", md_table(cnt), "")
  n_miss <- sum(cnt$missing %||% 0L)
  if (n_miss || any(cnt$unreadable > 0)) {
    md <- c(md, sprintf(paste("**%d reps missing against the design (%d units absent), %d unreadable",
                              "files: the run is incomplete.**"),
                        n_miss, sum(cnt$files == 0L), sum(cnt$unreadable)), "")
  }

  # Metrics, gated over the family.
  tab <- .rbind_fill(lapply(units, `[[`, "met"))
  if (!is.null(tab)) tab <- rp_try(sprintf("gate_family(%s)", fam_id), gate_family(tab))
  if (!is.null(tab)) {
    utils::write.csv(tab, file.path(tab_dir, paste0(fam_id, "_metrics.csv")), row.names = FALSE)
  }
  g <- if (is.null(tab)) NULL else tab[tab$gated, , drop = FALSE]
  vc <- if (is.null(g)) integer() else table(factor(g$verdict, c("pass", "fail", "na")))
  md <- c(md, "### Gates", "",
          if (is.null(g) || !nrow(g)) "No gated rows."
          else sprintf("%d gated rows: %d pass, %d fail, %d na (Holm over the family, tolerances %s).",
                       nrow(g), vc[["pass"]], vc[["fail"]], vc[["na"]], args$tier), "")
  if (!is.null(g) && any(g$verdict %in% "na")) {
    md <- c(md, sprintf(paste("**%d gated check(s) without a value (verdict na, never a pass): a gate",
                              "entry that matches no row (R = 0), or a value or MCSE that is NA.**"),
                        sum(g$verdict %in% "na")), "", rp_gate_rows(g[g$verdict %in% "na", ]), "")
  }
  if (!is.null(g) && any(g$verdict %in% "fail")) {
    md <- c(md, sprintf("**%d gated check(s) fail.**", sum(g$verdict %in% "fail")), "",
            rp_gate_rows(g[g$verdict %in% "fail", ]), "")
  }

  # One section per (cell, method).
  for (u in units) {
    md <- c(md, sprintf("### %s / cell %s / method %s", fam_id, u$cell, u$method), "",
            if (!u$files) "**Absent from the run: no rep file.**"
            else sprintf("%d reps: %d ok, %d failed.", u$files, u$ok, u$failed), "")
    ut <- if (is.null(tab)) NULL else tab[tab$cell == u$cell & tab$method == u$method, , drop = FALSE]
    md <- c(md, rp_metrics_md(ut), "")
    if (!is.null(u$fails)) md <- c(md, "Failures:", "", md_table(u$fails[c("error", "count")]), "")
    if (!is.null(u$warns)) {
      md <- c(md, "Warnings by class:", "",
              md_table(u$warns[c("stage", "class", "count", "reps", "example")]), "")
    }
    if (!is.null(u$msgs)) {
      md <- c(md, "Messages:", "", md_table(utils::head(u$msgs[c("message", "count", "reps")], 10L)),
              "")
    }
  }
  for (what in c("fails", "warns", "msgs")) {
    d <- .rbind_fill(lapply(units, `[[`, what))
    if (!is.null(d)) {
      utils::write.csv(d, file.path(tab_dir, sprintf("%s_%s.csv", fam_id,
                                                     c(fails = "failures", warns = "warnings",
                                                       msgs = "messages")[[what]])),
                       row.names = FALSE)
    }
  }
  if (!is.null(cnt)) utils::write.csv(cnt, file.path(tab_dir, paste0(fam_id, "_units.csv")),
                                      row.names = FALSE)

  # QQ plots at the gate cells.
  qq <- character()
  for (u in units) {
    if (!u$cell %in% gates$cell) next
    f <- rp_try(sprintf("QQ %s/%s/%s", fam_id, u$cell, u$method),
                rp_qq(u$z, file.path(plot_dir, sprintf("%s_qq_%s_%s.png", fam_id, u$cell, u$method)),
                      sprintf("%s cell %s, %s: z = (est - truth) / se_hessian", fam_id, u$cell,
                              u$method)))
    if (!is.null(f)) qq <- c(qq, sprintf("![QQ %s %s](plots/%s)", u$cell, u$method, f))
  }
  if (length(qq)) md <- c(md, "### QQ plots of z at the gate cells", "", qq, "")

  # Ladders against every numeric cell column that varies across a method's
  # cells; the rates of the estimands whose fits have two cells or more.
  rates <- list()
  lad_md <- character()
  if (!is.null(fam) && !is.null(tab)) {
    cols <- setdiff(names(fam$cells)[vapply(fam$cells, is.numeric, logical(1))], "cell")
    N <- stats::setNames(as.numeric(fam$cells$N), fam$cells$cell)
    big <- rp_big_blocks(tab)
    for (m in unique(tab$method)) {
      mt <- tab[tab$method == m & !is.na(tab$block) & !tab$block %in% big, , drop = FALSE]
      mcells <- unique(mt$cell)
      for (xc in cols) {
        x <- stats::setNames(as.numeric(fam$cells[[xc]]), fam$cells$cell)[mcells]
        if (length(unique(x[is.finite(x)])) < 2L) next
        lad <- rp_try(sprintf("ladder %s/%s/%s", fam_id, m, xc), rp_ladder_rows(mt, N))
        if (is.null(lad)) next
        f <- rp_try(sprintf("ladder plot %s/%s/%s", fam_id, m, xc),
                    rp_ladder_plot(lad, x, xc, file.path(plot_dir, sprintf("%s_ladder_%s_%s.png",
                                                                          fam_id, m, xc)),
                                   sprintf("%s, method %s, against %s", fam_id, m, xc)))
        if (!is.null(f)) lad_md <- c(lad_md, sprintf("![ladder %s %s](plots/%s)", m, xc, f), "")
        rates[[length(rates) + 1L]] <- rp_try(sprintf("rates %s/%s/%s", fam_id, m, xc),
                                              rp_rates(mt, lad, x, xc, m))
      }
    }
  }
  if (length(lad_md)) md <- c(md, "### Ladders (MCSE bars)", "", lad_md)
  rates <- .rbind_fill(rates)
  if (!is.null(rates)) {
    utils::write.csv(rates, file.path(tab_dir, paste0(fam_id, "_rates.csv")), row.names = FALSE)
    r <- rates[is.finite(rates$slope), , drop = FALSE]
    md <- c(md, "### Rates (rate_fit: slope of log|value| on log x)", "",
            if (!nrow(r)) "No slope: fewer than two cells with a finite, non-zero value."
            else md_table(data.frame(method = r$method, x = r$x, estimand = r$estimand,
                                     metric = r$metric, slope = .pm(r$slope, r$se), n = r$n,
                                     predicted = .fmt(r$predicted))), "")
  }

  # SBC: the plots; the verdicts are the metrics table's sbc_ecdf rows
  # (statistic, exact p, Holm p and verdict from gate_family()); the plots'
  # p-values are listed, ungated, only for units whose metrics carry none.
  sb <- list()
  sb_md <- character()
  for (u in units) {
    s <- rp_try(sprintf("SBC %s/%s/%s", fam_id, u$cell, u$method),
                rp_sbc(u$sbc, file.path(plot_dir, sprintf("%s_sbc_%s_%s.png", fam_id, u$cell,
                                                          u$method)),
                       sprintf("%s cell %s, %s: SBC ECDF differences", fam_id, u$cell, u$method)))
    if (is.null(s)) next
    sb_md <- c(sb_md, sprintf("![SBC %s %s](plots/%s)", u$cell, u$method, s$file), "")
    rows <- if (is.null(tab)) NULL
            else tab[tab$cell == u$cell & tab$method == u$method & tab$metric == "sbc_ecdf", ]
    if (!is.null(rows) && nrow(rows)) {
      sb_md <- c(sb_md, md_table(data.frame(
        cell = rows$cell, method = rows$method, quantity = rows$estimand, R = rows$R,
        statistic = .fmt(rows$value), p = .fmt(rows$p, 3L), p_holm = .fmt(rows$p_holm, 3L),
        H = ifelse(rows$gated, "H", ""), verdict = rows$verdict)), "")
    } else {
      sb[[length(sb) + 1L]] <- data.frame(cell = u$cell, method = u$method, s$tab)
    }
  }
  sb <- .rbind_fill(sb)
  if (length(sb_md)) md <- c(md, "### SBC ECDF tests", "", sb_md)
  if (!is.null(sb)) {
    utils::write.csv(sb, file.path(tab_dir, paste0(fam_id, "_sbc.csv")), row.names = FALSE)
    md <- c(md, "Ungated (the metrics table has no sbc_ecdf rows for these units):", "",
            md_table(data.frame(sb[c("cell", "method", "quantity", "R", "L")], p = .fmt(sb$p, 3L))),
            "")
  }

  summ <- data.frame(family = fam_id, units = nrow(keys), absent = sum(cnt$files == 0L),
                     missing = n_miss, files = sum(cnt$files %||% 0L),
                     ok = sum(cnt$ok %||% 0L), failed = sum(cnt$failed %||% 0L),
                     gated = if (is.null(g)) 0L else nrow(g),
                     pass = if (length(vc)) vc[["pass"]] else 0L,
                     fail = if (length(vc)) vc[["fail"]] else 0L,
                     na = if (length(vc)) vc[["na"]] else 0L,
                     findings = if (is.null(tab)) 0L else sum(!tab$gated),
                     stringsAsFactors = FALSE)
  list(md = md, summary = summ)
}

# --- Run metadata ----------------------------------------------------------------

rp_meta <- function() {
  files <- list.files(run_dir, pattern = "^meta_.*\\.rds$", full.names = TRUE)
  if (!length(files)) return(c("## Run metadata", "", "No meta_*.rds file in the run.", ""))
  ms <- lapply(files, readRDS)
  thr <- function(t) if (is.list(t)) as.character(t$omp_get_max_threads %||% NA) else NA
  tab <- data.frame(
    file = basename(files), host = vapply(ms, function(m) m$host %||% NA_character_, ""),
    time = vapply(ms, function(m) m$time %||% NA_character_, ""),
    choicer = vapply(ms, function(m) m$choicer_version %||% NA_character_, ""),
    path = vapply(ms, function(m) m$choicer_path %||% NA_character_, ""),
    ref = vapply(ms, function(m) as.character(m$ref %||% ""), ""),
    threads = vapply(ms, function(m) thr(m$threads), ""),
    git = vapply(ms, function(m) substr(m$git_sha %||% NA_character_, 1L, 10L), ""),
    dirty = vapply(ms, function(m) as.character(m$git_dirty %||% NA), ""),
    stringsAsFactors = FALSE)
  utils::write.csv(tab, file.path(tab_dir, "meta.csv"), row.names = FALSE)
  m <- ms[[1]]
  s <- m$session
  pk <- c(s$otherPkgs, s$loadedOnly)
  pk <- pk[intersect(c("choicer", "data.table", "Rcpp", "RcppArmadillo", "nloptr", "randtoolbox"),
                     names(pk))]
  out <- c("## Run metadata", "", md_table(tab), "",
           sprintf("- R: %s; platform %s; %s", s$R.version$version.string %||% NA,
                   s$platform %||% NA, s$running %||% NA),
           sprintf("- Packages: %s", paste(names(pk), vapply(pk, function(p) p$Version %||% "?", ""),
                                           collapse = ", ")),
           sprintf("- BLAS: %s; LAPACK: %s (%s)", m$blas %||% NA, m$lapack %||% NA,
                   m$lapack_version %||% NA),
           sprintf("- CPU: %s", m$cpu %||% NA), "")
  mixed <- vapply(c("choicer", "path", "ref", "git", "threads"),
                  function(k) length(unique(tab[[k]])) > 1L, logical(1))
  if (any(mixed)) {
    out <- c(out, sprintf("**The run's processes differ in: %s.**",
                          paste(names(mixed)[mixed], collapse = ", ")), "")
  }
  dp <- unique(unlist(lapply(ms, `[[`, "git_dirty_paths")))
  if (length(dp) && !all(is.na(dp))) {
    out <- c(out, sprintf("Uncommitted paths: %s", paste(dp[!is.na(dp)], collapse = ", ")), "")
  }
  out
}

# --- Main -------------------------------------------------------------------------

idx <- .vs_rep_files(run_dir)
fams <- sort(unique(idx$family), method = "radix")
tl <- tolerances(args$tier)
parts <- lapply(fams, function(f) {
  r <- rp_try(sprintf("family %s", f), rp_family(f, idx[idx$family == f, , drop = FALSE]))
  r %||% list(md = c(sprintf("## %s", f), "", "(the family's section failed; see Report errors)",
                     ""), summary = NULL)
})
summ <- .rbind_fill(lapply(parts, `[[`, "summary"))
if (!is.null(summ)) utils::write.csv(summ, file.path(tab_dir, "summary.csv"), row.names = FALSE)
meta_md <- rp_try("run metadata", rp_meta())
if (!length(fams)) errors <- c(errors, "no rep files in the run")

# The summary, a family's missing reps in bold when there are any.
summ_md <- summ
if (!is.null(summ_md)) {
  summ_md$missing <- ifelse(summ$missing > 0, sprintf("**%d**", summ$missing), "0")
}
md <- c(sprintf("# Validation report: %s", basename(run_dir)), "",
        sprintf("Run directory `%s`; %s; tier %s, tolerances %s; report by choicer %s.", run_dir,
                format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), args$tier,
                paste(names(tl), tl, sep = " ", collapse = ", "), build$version), "",
        "## Summary", "", md_table(summ_md), "")
for (v in c("missing", "fail", "na")) {
  k <- if (is.null(summ)) integer() else which(summ[[v]] > 0)
  if (length(k)) {
    what <- c(missing = "Reps missing against the design (an incomplete run)",
              fail = "Gated rows failing", na = "Gated rows without a value (na)")[[v]]
    md <- c(md, sprintf("**%s: %s.**", what,
                        paste(sprintf("%s (%d)", summ$family[k], summ[[v]][k]), collapse = ", ")), "")
  }
}
if (length(errors)) {
  md <- c(md, "## Report errors", "", paste0("- ", .short(errors, 400L)), "")
}
md <- c(md, unlist(lapply(parts, `[[`, "md")), meta_md)
writeLines(md, file.path(run_dir, "REPORT.md"))
cat(sprintf("report: %s (%d families, %d rep files)\n", file.path(run_dir, "REPORT.md"),
            length(fams), nrow(idx)))
if (!is.null(summ)) {
  cat(sprintf("  %s: %d gated rows, %d pass, %d fail, %d na; %d reps missing (%d units absent)\n",
              summ$family, summ$gated, summ$pass, summ$fail, summ$na, summ$missing, summ$absent),
      sep = "")
}
if (length(errors)) {
  cat("report errors:\n", paste0("  ", errors, "\n"), sep = "")
  quit(save = "no", status = 1L)
}
