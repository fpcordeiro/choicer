# The cluster meat kernels take labels coded 1..G in order of first
# appearance, the grouping of rowsum() over as.character(labels).
# .cluster_codes() forms them from the labels' own values where those group
# as their strings do (factor codes, plain integer, logical and character
# vectors) and through as.character() otherwise, so the codes, and the
# meats, are those of the string route.

cc_string_codes <- function(x) {
  s <- as.character(x)
  match(s, unique(s))
}

cc_label_sets <- function() {
  set.seed(6101)
  g <- sample.int(7L, 60L, replace = TRUE)
  lev <- c("b", "a", "unused", "c", "e", "d", "g", "f")
  sets <- list(
    integer = g,
    logical = g %% 2L == 0L,
    character = letters[g],
    factor = factor(letters[g]),
    factor_reordered = factor(letters[g], levels = rev(letters[1:7])),
    factor_unused = factor(letters[g], levels = c(letters[1:7], "zz")),
    factor_na_level = addNA(factor(ifelse(g == 3L, NA, letters[g]))),
    # two levels with one label: as.character() merges them
    factor_dup_levels = structure(g, levels = c("a", "b", "a", "c", "d", "c", "e"),
                                  class = "factor"),
    # 0.1 + 0.2 and 0.3 print alike: as.character() merges them
    double = c(0.1, 0.3, 0.1 + 0.2, 1, 2)[(g %% 5L) + 1L],
    date = as.Date("2020-01-01") + g,
    named_integer = setNames(g, paste0("u", seq_along(g)))
  )
  if (requireNamespace("bit64", quietly = TRUE)) {
    sets$integer64 <- bit64::as.integer64(g * 1e12)
  }
  sets
}

test_that(".cluster_codes() numbers clusters as the string route does", {
  for (nm in names(cc_label_sets())) {
    x <- cc_label_sets()[[nm]]
    expect_identical(.cluster_codes(x), cc_string_codes(x), label = nm)
  }
  expect_identical(.cluster_codes(c(0.1 + 0.2, 0.3, 0.5)), c(1L, 1L, 2L))
  expect_identical(.cluster_codes(structure(c(1L, 2L, 3L),
                                            levels = c("x", "x", "y"),
                                            class = "factor")),
                   c(1L, 1L, 2L))
})

test_that(".cluster_codes() takes the string route wherever the codes might not group alike", {
  e_latin1 <- "\xe9"
  Encoding(e_latin1) <- "latin1"
  z_bytes <- "z\xff"
  Encoding(z_bytes) <- "bytes"
  sets <- list(
    # an unused level marked "bytes" makes the levels compare bytewise, but
    # as.character() of the labels translates the latin1 and UTF-8 ones alike
    factor_bytes_level = structure(c(1L, 2L, 1L),
                                   levels = c(e_latin1, enc2utf8(e_latin1), z_bytes),
                                   class = "factor"),
    # an NA level and NA codes, both NA as strings
    factor_na_level_and_code = structure(c(1L, NA, 2L, 1L), levels = c("a", NA),
                                         class = "factor"),
    ordered = factor(c("lo", "hi", "lo"), levels = c("lo", "hi"), ordered = TRUE),
    # a matrix's unique() is unique.matrix()'s rows; its strings are a vector
    matrix = matrix(c(1L, 1L, 2L, 3L), ncol = 2L),
    character_matrix = matrix(c("a", "a", "b", "c"), ncol = 2L)
  )
  for (nm in names(sets)) {
    expect_identical(.cluster_codes(sets[[nm]]), cc_string_codes(sets[[nm]]), label = nm)
  }
  # malformed factors stop as as.character() does
  for (x in list(structure(c(1L, 5L, 1L), levels = c("a", "b"), class = "factor"),
                 structure(c(1L, 2L, 1L), class = "factor"),
                 structure(c(1L, 2L, 1L), levels = c(10, 20), class = "factor"))) {
    expect_identical(tryCatch(.cluster_codes(x), error = conditionMessage),
                     tryCatch(cc_string_codes(x), error = conditionMessage))
  }
  # a factor subclass keeps its own as.character() method
  registerS3method("as.character", "cc_caseless",
                   function(x, ...) toupper(levels(x))[unclass(x)])
  caseless <- structure(c(1L, 2L, 3L, 1L), levels = c("a", "A", "b"),
                        class = c("cc_caseless", "factor"))
  expect_identical(.cluster_codes(caseless), c(1L, 1L, 2L, 1L))
})

test_that("cluster meats do not depend on the labels' type for one partition", {
  on.exit(set_omp_threads(2L), add = TRUE)
  set_omp_threads(1L)
  set.seed(6111)
  N <- 80L; J <- 4L
  dt <- data.table(id = rep(seq_len(N), each = J), alt = rep(seq_len(J), N),
                   x1 = rnorm(N * J), x2 = rnorm(N * J))
  dt[, choice := as.integer(seq_len(.N) == sample.int(.N, 1L)), by = id]
  dt[, nest := ifelse(alt <= 2L, "A", "B")]
  fits <- list(mnl = run_mnlogit(dt, "id", "alt", "choice", c("x1", "x2")),
               nl = run_nestlogit(dt, "id", "alt", "choice", c("x1", "x2"), "nest"))
  g <- sample.int(9L, N, replace = TRUE)
  labels <- list(integer = g, character = sprintf("g%02d", g),
                 factor = factor(sprintf("g%02d", g)),
                 factor_reordered = factor(g, levels = rev(seq_len(9L))),
                 double = g + 0.5)
  for (m in names(fits)) {
    # labels named by situation id, as vcov() asks
    ids <- fits[[m]]$data$situation_ids
    named <- lapply(labels, function(x) { names(x) <- ids; x })
    ref <- vcov(fits[[m]], type = "cluster", cluster = named$integer)
    for (nm in names(named)[-1L]) {
      expect_identical(vcov(fits[[m]], type = "cluster", cluster = named[[nm]]),
                       ref, label = paste(m, nm))
    }
  }
})
