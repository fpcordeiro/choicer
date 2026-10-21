#!/usr/bin/env Rscript
# Time and memory of one data preparation, prepare_{mxl,mnl,nl}_data(), on a
# synthetic long data frame: one build (R_LIBS), one configuration, one process.
#
#   R_LIBS=<lib> Rscript _benchmarks/scale/prep_mem.R regime=claims U=15000 model=mxl \
#     [tag=<label>] [csv=<file.csv>]
#   (under the malloc shim, output/Rmc with MC_SHIM set, it also reports the
#   heap blocks and bytes moved, the peak live heap and the sampled peak RSS)
#
# The long data hold person, id (choice situation), alt, choice, x1..x10 and
# w1..w3, plus nest (alt %% 10) for the nested logit, built from gen_regime().
# The MXL preparation takes person as its decision-maker column.
source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE),
  value = TRUE)))), "scale_common.R"))
lib <- load_choicer()
opt <- scale_args()
U <- as.integer(opt$U %||% "15000"); regime <- opt$regime %||% "claims"
model <- opt$model %||% "mxl"; tag <- opt$tag %||% basename(lib)
d <- gen_regime(regime, U, seed = 1)
sit <- rep(seq_along(d$M), d$M)
person <- rep(rep(seq_along(d$Ti), d$Ti), d$M)
chosen <- integer(sum(d$M)); chosen[cumsum(d$M) - d$M + d$choice_idx] <- 1L
df <- data.frame(person = person, id = sit, alt = d$alt_idx, choice = chosen)
Xc <- paste0("x", seq_len(ncol(d$X))); Wc <- paste0("w", seq_len(ncol(d$W)))
for (k in seq_along(Xc)) df[[Xc[k]]] <- d$X[, k]
for (k in seq_along(Wc)) df[[Wc[k]]] <- d$W[, k]
if (model == "nl") df$nest <- d$alt_idx %% 10L
rm(d, sit, person, chosen); invisible(gc())
sz <- as.numeric(object.size(df))
call_prep <- function() {
  switch(model,
    mxl = prepare_mxl_data(df, id_col = "id", alt_col = "alt", choice_col = "choice",
                           covariate_cols = Xc, random_var_cols = Wc, person_col = "person"),
    mnl = prepare_mnl_data(df, id_col = "id", alt_col = "alt", choice_col = "choice",
                           covariate_cols = c(Xc, Wc)),
    nl = prepare_nl_data(df, id_col = "id", alt_col = "alt", choice_col = "choice",
                         covariate_cols = c(Xc, Wc), nest_col = "nest"),
    stop("unknown model ", model))
}
invisible(gc()); invisible(gc())
rss_base <- rss_now()[["rss"]]
mc_reset(); rss_start(5)
t0 <- Sys.time()
prep <- call_prep()
t_prep <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
m <- rss_stop(); a <- mc_get(); mc_disable()
res <- data.frame(tag = tag, model = model, regime = regime, U = U, rows = nrow(df),
                  input_gb = sz / 1e9, time_s = t_prep, shim = shim_loaded(),
                  peak_rss_gb = (m[["peak_rss"]] - rss_base) / 1e9,
                  peak_live_gb = a[["peak_live"]] / 1e9,
                  alloc_gb = a[["bytes_alloc"]] / 1e9, blocks = a[["n_alloc"]],
                  result_gb = as.numeric(object.size(prep)) / 1e9, load = loadavg())
res$peak_rss_x <- res$peak_rss_gb / res$input_gb
res$peak_live_x <- res$peak_live_gb / res$input_gb
cat(sprintf("[%s] %s %s U=%d rows=%.3g: %.2fs | peak RSS +%.2f GB (%.2fx input %.2f GB), peak live +%.2f GB (%.2fx) | %.1f GB in %.0f blocks | result %.2f GB\n",
            tag, model, regime, U, res$rows, t_prep, res$peak_rss_gb, res$peak_rss_x, res$input_gb,
            res$peak_live_gb, res$peak_live_x, res$alloc_gb, res$blocks, res$result_gb))
if (!is.null(opt$csv)) {
  utils::write.table(res, opt$csv, sep = ",", row.names = FALSE,
                     col.names = !file.exists(opt$csv), append = file.exists(opt$csv))
}
