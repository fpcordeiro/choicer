# Warnings that most tests' fits raise and do not concern. Arguments pass through
# `...`, so the fit still sees which ones are missing. A test about one of these
# calls choicer::run_*().
#
# - Many tests stop a fit at a cap on purpose, for speed, and the fit warns
#   (choicer_optimizer_capped): the subject of test-optimizer-failure.R's last
#   section.
# - Most mixed logit tests fit small samples without real heterogeneity, where a
#   variance ends at zero and the fit warns (choicer_variance_plateau):
#   test-mxl-plateau.R's subject.
run_mnlogit <- function(...) {
  suppressWarnings(choicer::run_mnlogit(...), classes = "choicer_optimizer_capped")
}
run_nestlogit <- function(...) {
  suppressWarnings(choicer::run_nestlogit(...), classes = "choicer_optimizer_capped")
}
run_mxlogit <- function(...) {
  suppressWarnings(choicer::run_mxlogit(...),
                   classes = c("choicer_optimizer_capped", "choicer_variance_plateau"))
}
