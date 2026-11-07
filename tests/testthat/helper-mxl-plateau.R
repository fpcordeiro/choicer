# Most mixed logit tests fit small samples without real heterogeneity, where a
# variance ends at zero and the fit warns (choicer_variance_plateau): that
# warning is test-mxl-plateau.R's subject, which calls choicer::run_mxlogit().
# Arguments pass through `...`, so the fit still sees which ones are missing.
run_mxlogit <- function(...) {
  suppressWarnings(choicer::run_mxlogit(...), classes = "choicer_variance_plateau")
}
