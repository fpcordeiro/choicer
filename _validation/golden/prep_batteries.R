# ============================================================================
# Phase 6A batteries for the scratch golden harness (sourced by
# refactor_identity_check.R). They run on both main and the candidate build,
# so they only call functions that exist in both:
#   * prep: prepare_mnl_data(), prepare_mxl_data(), prepare_nl_data() over a
#     grid of inputs (data.frame, keyed data.table, tibble, list and matrix
#     inputs, extra columns of every type, NA rows, outside option, weights,
#     clusters, panels, integer and mixed covariates, collinear columns, error
#     cases) plus two designs with n * p >= 1e6 and a dependent column, which
#     take the chunked-QR path in the Phase 6A build.
#   * cube: get_halton_normals() over (S, N, K_w), including sizes that span
#     several blocks of the Phase 6A builder (digests for the large ones).
# Every prepared object is stored after a serialize round trip (a live
# data.table is never identical() to a deserialized one) with X and W in
# double storage (Phase 6A always builds double; main returns integer for
# all-integer designs; the kernels convert either way), and is marked as a
# leaf so flatten() does not recurse into it.
# ============================================================================

rt <- function(x) unserialize(serialize(x, NULL))

as_leaf <- function(x) {
  x <- rt(x)
  if (is.list(x)) attr(x, "golden_leaf") <- TRUE
  x
}

normalize_prep <- function(x) {
  if (is.list(x)) {
    for (nm in c("X", "W")) {
      if (is.matrix(x[[nm]]) && is.integer(x[[nm]])) storage.mode(x[[nm]]) <- "double"
    }
  }
  x
}

# Run f(...) capturing value, warnings, messages, error (message and the
# function name of its call) and whether the data argument was left unchanged.
capture_call <- function(fname, args) {
  before <- rt(args$data)
  warns <- character(); msgs <- character()
  err_msg <- NA_character_; err_fun <- NA_character_
  value <- tryCatch(
    withCallingHandlers(
      do.call(fname, args, envir = asNamespace("choicer")),
      warning = function(w) {
        warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning")
      },
      message = function(m) {
        msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage")
      }),
    error = function(e) {
      err_msg <<- conditionMessage(e)
      cl <- conditionCall(e)
      err_fun <<- if (is.null(cl)) "<none>" else paste(deparse(cl[[1]]), collapse = "")
      NULL
    })
  list(value = if (is.null(value)) NULL else as_leaf(normalize_prep(value)),
       warnings = warns, messages = msgs, error = err_msg, error_fun = err_fun,
       unchanged = identical(rt(args$data), before))
}

# Long-format data: N choice situations nested in persons (T per person),
# unbalanced choice sets over J alternatives, every column type a user might
# carry along.
make_long <- function(seed, N, J, T = 1L, min_alts = 2L) {
  set.seed(seed)
  sets <- lapply(seq_len(N), function(i) {
    if (i == 1L) seq_len(J) else sort(sample.int(J, sample(min_alts:J, 1L)))
  })
  M <- lengths(sets)
  n <- sum(M)
  id <- rep(seq_len(N), M)
  person <- rep(rep(seq_len(ceiling(N / T)), each = T)[seq_len(N)], M)
  alt <- unlist(sets)
  choice <- integer(n)
  first <- cumsum(M) - M
  choice[first + vapply(M, function(m) sample.int(m, 1L), 1L)] <- 1L
  x1 <- rnorm(n); x2 <- rnorm(n); w1 <- rnorm(n); w2 <- rnorm(n); w3 <- rnorm(n)
  df <- data.frame(
    person = person, id = id, alt = alt, choice = choice,
    x1 = x1, x2 = x2, x3 = 2 * x1 - x2, w1 = w1, w2 = w2, w3 = w3,
    xi1 = sample(0:3, n, TRUE), xi2 = sample(-2:2, n, TRUE),
    wt = rep(runif(N, 0.5, 2), M), cl = rep(sample.int(max(2L, N %/% 4L), N, TRUE), M),
    nest = ifelse(alt <= ceiling(J / 2), "A", "B"),
    note = sample(c("a", "b", NA), n, TRUE),
    stringsAsFactors = FALSE)
  df$pw <- rep(runif(ceiling(N / T), 0.5, 2), each = T)[df$person]
  df$pcl <- df$person %% 3L
  df$lst <- as.list(seq_len(n))
  df
}

drop_outside_choice <- function(df, frac = 0.25) {
  # Make a share of situations choose nobody (implicit outside option).
  ids <- unique(df$id)
  none <- ids[seq(2L, length(ids), by = max(2L, round(1 / frac)))]
  df$choice[df$id %in% none] <- 0L
  df
}

run_prep <- function() {
  out <- list()
  base <- make_long(601, N = 60, J = 6, T = 3)
  keyed <- data.table::as.data.table(base)
  data.table::setkeyv(keyed, "id")
  data.table::setindexv(keyed, "alt")
  na <- base
  na$x1[na$id %in% c(5, 17)] <- NA
  na$w2[na$id == 23 & na$alt == max(na$alt[na$id == 23])] <- NA
  na$note[1:5] <- NA   # unneeded column: must not drop anything
  oo <- drop_outside_choice(base)
  b0 <- base[, setdiff(names(base), "lst")]
  labeled <- rbind(b0, transform(b0[!duplicated(b0$id), ], alt = 0L, choice = 0L, nest = "O"))
  labeled <- labeled[order(labeled$id, labeled$alt), ]
  fac <- base; fac$alt <- factor(paste0("alt", fac$alt))
  chr <- base; chr$id <- paste0("s", chr$id)
  attr_alt <- base; attr(attr_alt$alt, "label") <- "alternative"
  tib <- tibble::as_tibble(base[, setdiff(names(base), "lst")])
  lst <- as.list(base[, c("id", "alt", "choice", "x1", "x2", "w1", "w2", "nest")])
  mat <- as.matrix(base[, c("id", "alt", "choice", "x1", "x2", "w1", "w2")])
  dup <- base[, c("id", "alt", "choice", "x1", "x2", "w1", "w2", "nest")]
  dup <- cbind(dup, x1 = rnorm(nrow(dup)))           # duplicated name "x1"
  matcol <- base[, c("id", "alt", "choice", "x1", "x2", "w1", "w2", "nest")]
  matcol$m <- matrix(rnorm(2 * nrow(matcol)), ncol = 2)
  bad_choice <- base; bad_choice$choice[3] <- 2L
  non_num <- base; non_num$x2 <- as.character(non_num$x2)
  inf_cov <- base; inf_cov$x2[7] <- Inf
  big <- make_long(602, N = 100000, J = 4, min_alts = 4L)   # 4e5 rows
  big <- big[, c("id", "alt", "choice", "x1", "x2", "x3", "w1", "w2", "w3", "nest")]
  big$w3 <- big$w1 + big$w2                                  # dependent W column

  mnl <- function(data, cov = c("x1", "x2"), ...) {
    capture_call("prepare_mnl_data", list(data = data, id_col = "id", alt_col = "alt",
                                          choice_col = "choice", covariate_cols = cov, ...))
  }
  mxl <- function(data, cov = c("x1", "x2"), rv = c("w1", "w2"), ...) {
    capture_call("prepare_mxl_data", list(data = data, id_col = "id", alt_col = "alt",
                                          choice_col = "choice", covariate_cols = cov,
                                          random_var_cols = rv, ...))
  }
  nl <- function(data, cov = c("x1", "x2"), nest = "nest", ...) {
    capture_call("prepare_nl_data", list(data = data, id_col = "id", alt_col = "alt",
                                         choice_col = "choice", covariate_cols = cov,
                                         nest_col = nest, ...))
  }

  for (m in c("mnl", "mxl", "nl")) {
    f <- get(m)
    out[[paste0(m, "_df")]] <- f(base)
    out[[paste0(m, "_keyed_dt")]] <- f(keyed)
    out[[paste0(m, "_tibble")]] <- f(tib)
    out[[paste0(m, "_na")]] <- f(na)
    out[[paste0(m, "_oo")]] <- f(oo, include_outside_option = TRUE)
    out[[paste0(m, "_oo_label")]] <- f(labeled, include_outside_option = TRUE,
                                       outside_opt_label = 0L)
    out[[paste0(m, "_label_first")]] <- f(labeled, outside_opt_label = 0L,
                                          include_outside_option = FALSE)
    out[[paste0(m, "_wcol")]] <- f(base, weights_col = "wt")
    out[[paste0(m, "_wvec")]] <- f(base, weights = runif(60, 0.5, 2))
    out[[paste0(m, "_cl")]] <- f(base, cluster_col = "cl")
    out[[paste0(m, "_int")]] <- if (m == "mxl") f(base, cov = c("xi1", "xi2"), rv = c("xi1", "xi2"))
                                else f(base, cov = c("xi1", "xi2"))
    out[[paste0(m, "_mixed")]] <- f(base, cov = c("xi1", "x2"))
    out[[paste0(m, "_collinear")]] <- f(base, cov = c("x1", "x2", "x3"))
    out[[paste0(m, "_factor_alt")]] <- f(fac)
    out[[paste0(m, "_char_id")]] <- f(chr)
    out[[paste0(m, "_attr_alt")]] <- f(attr_alt)
    out[[paste0(m, "_list_input")]] <- f(lst)
    out[[paste0(m, "_dup_names")]] <- f(dup)
    out[[paste0(m, "_matrix_col")]] <- f(matcol)
    out[[paste0(m, "_err_missing")]] <- f(base, cov = c("x1", "nope"))
    out[[paste0(m, "_err_both_w")]] <- f(base, weights = rep(1, 60), weights_col = "wt")
    out[[paste0(m, "_err_choice")]] <- f(bad_choice)
    out[[paste0(m, "_err_non_num")]] <- f(non_num)
    out[[paste0(m, "_err_inf")]] <- f(inf_cov)
    out[[paste0(m, "_err_wcol_varies")]] <- f(base, weights_col = "x1")
    out[[paste0(m, "_big")]] <- f(big, cov = c("x1", "x2", "x3"))
  }
  out$mnl_matrix_input <- mnl(mat)
  out$mxl_matrix_input <- mxl(mat)
  out$mxl_person <- mxl(base, person_col = "person")
  out$mxl_person_w <- mxl(base, person_col = "person", weights_col = "pw", cluster_col = "pcl")
  out$mxl_person_err_w <- mxl(base, person_col = "person", weights_col = "wt")
  out$mxl_rc_corr <- mxl(base, rc_correlation = TRUE)
  out$mxl_rv3 <- mxl(base, rv = c("w1", "w2", "w3"))
  out$mxl_big_w <- mxl(big, cov = c("x1", "x2"), rv = c("w1", "w2", "w3"))
  out$mxl_overlap <- mxl(base, cov = c("x1", "w1"), rv = c("w1", "w2"))
  out$nl_err_no_nest_col <- nl(base, nest = "nope")
  multi <- base; multi$nest[multi$alt == 1][1] <- "B"
  out$nl_err_multi_nest <- nl(multi)
  one <- base; one$nest <- "A"
  out$nl_err_one_nest <- nl(one)
  na_nest <- base; na_nest$nest[na_nest$alt == 2] <- NA
  out$nl_err_na_nest <- nl(na_nest)
  miss_alt <- base; names(miss_alt)[names(miss_alt) == "alt"] <- "alt0"
  out$nl_err_missing_alt <- nl(miss_alt)
  out$mnl_err_missing_alt <- mnl(miss_alt)
  # Named covariate vectors: names must not reach colnames(X) / colnames(W).
  out$mnl_named_cols <- mnl(base, cov = c(a = "x1", b = "x2"))
  out$mxl_named_cols <- mxl(base, cov = c(a = "x1"), rv = c(w = "w1", v = "w2"))
  # A data.frame subclass whose as.data.frame() method changes its columns
  # (as plm's pdata.frame does): data.table >= 1.17 calls that method.
  sub <- base[, setdiff(names(base), "lst")]
  class(sub) <- c("golden_sub", "data.frame")
  assign("as.data.frame.golden_sub", function(x, ...) {
    x <- unclass(x)
    x <- lapply(x, function(v) structure(v, golden = TRUE))
    as.data.frame.list(x, stringsAsFactors = FALSE)
  }, envir = globalenv())
  registerS3method("as.data.frame", "golden_sub",
                   get("as.data.frame.golden_sub", envir = globalenv()))
  out$mnl_df_subclass <- mnl(sub)
  # Duplicated name of an unused column whose second copy has NAs: main's
  # `:=` drop removes only the first copy, and the survivor's NAs drop rows.
  dupz <- base[, c("id", "alt", "choice", "x1", "x2", "w1", "w2", "nest")]
  dupz <- cbind(dupz, z = 1, z = ifelse(dupz$id %in% c(4, 9), NA, 2))
  out$mnl_dup_unused_na <- mnl(dupz)
  out$mxl_dup_unused_na <- mxl(dupz)
  out$nl_dup_unused_na <- nl(dupz)
  out$mxl_df_subclass <- mxl(sub)
  out$nl_df_subclass <- nl(sub)
  out
}

cube_digest <- function(x) {
  f <- tempfile(fileext = ".bin")
  on.exit(unlink(f))
  con <- file(f, "wb"); writeBin(as.vector(x), con); close(con)
  c(md5 = unname(tools::md5sum(f)), dim = paste(dim(x), collapse = "x"),
    type = typeof(x))
}

run_cube <- function() {
  small <- list(c(1, 1, 1), c(7, 3, 1), c(50, 20, 2), c(100, 30, 3), c(13, 1000, 5),
                c(1, 1000, 13), c(1000, 7, 1))
  out <- list()
  for (p in small) {
    out[[sprintf("S%d_N%d_K%d", p[1], p[2], p[3])]] <-
      get_halton_normals(S = p[1], N = p[2], K_w = p[3])
  }
  # Several blocks of 2^22 draws, the last one partial.
  large <- list(c(100, 30001, 3), c(1000, 5003, 1), c(7, 400009, 2))
  for (p in large) {
    out[[sprintf("S%d_N%d_K%d_md5", p[1], p[2], p[3])]] <-
      cube_digest(get_halton_normals(S = p[1], N = p[2], K_w = p[3]))
  }
  out
}
