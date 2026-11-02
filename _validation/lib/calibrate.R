# ============================================================================
# Calibration of a registered DGP: simulates it at a large N and prints its
# calibration moments against their targets (plan §4), then its truth theta0.
# The DGPs never call choicer's kernels; their theta0 is frozen in the
# registry, and this script shows how the moments sit at it.
#
# Usage (the build under test is the FIRST entry of R_LIBS; one thread):
#   R_LIBS=<lib> OMP_NUM_THREADS=1 Rscript _validation/lib/calibrate.R --dgp=<id> \
#     [--variant=default] [--N=200000] [--seed=1]
# --N is the DGP's size argument (situations, or decision makers in a panel).
# Exits 1 with a message for an unknown DGP or variant.
# ============================================================================

main <- function() {
  lib_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE),
                                                            value = TRUE)[1])))
  source(file.path(lib_dir, "families.R"), local = (boot <- new.env()))
  L <- boot$vs_source_lib(lib_dir)
  a <- L$vs_args(defaults = list(dgp = NULL, variant = "default", N = 2e5, seed = 1L))
  if (is.null(a$dgp)) stop("name the DGP: --dgp=<id> (registered: ",
                           paste(L$list_dgps(), collapse = ", "), ")", call. = FALSE)
  L$vs_load_choicer()
  dgp <- L$get_dgp(a$dgp)
  if (!a$variant %in% dgp$variants) {
    stop("DGP '", dgp$id, "' has no variant '", a$variant, "'; variants: ",
         paste(dgp$variants, collapse = ", "), call. = FALSE)
  }
  L$vs_threads(1L)
  t0 <- proc.time()[["elapsed"]]
  sim <- dgp$simulate(N = a$N, seed = a$seed, variant = a$variant)
  m <- dgp$moments(sim)
  sec <- proc.time()[["elapsed"]] - t0

  tg <- dgp$targets
  mom <- union(names(m), tg$moment)
  i <- match(mom, tg$moment)
  tab <- data.frame(moment = mom, value = unname(m[mom]), lo = tg$lo[i], hi = tg$hi[i],
                    stringsAsFactors = FALSE)
  tab$within <- ifelse(is.na(tab$lo) & is.na(tab$hi), NA,
                       tab$value >= tab$lo & tab$value <= tab$hi)
  cat(sprintf("DGP %s (%s), variant %s, N = %s, seed %d: simulated in %.1f s\n\n", dgp$id,
              dgp$name, a$variant, format(a$N, scientific = FALSE), a$seed, sec))
  print(tab, row.names = FALSE, digits = 4)
  n_t <- sum(!is.na(tab$within))
  cat(sprintf("\n%d of %d targeted moments within their targets\n\ntheta0:\n",
              sum(tab$within, na.rm = TRUE), n_t))
  print(L$dgp_theta0(dgp, a$variant))
  invisible(tab)
}

status <- tryCatch({ main(); 0L }, error = function(e) {
  message("calibrate: ", conditionMessage(e))
  1L
})
quit(save = "no", status = status)
