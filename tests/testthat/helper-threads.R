# Set choicer's threads where OpenMP runs. A build without OpenMP runs every
# kernel on one thread, and set_num_threads() there only warns, so tests set
# threads through this helper; test-utilities.R tests set_num_threads() itself.
set_omp_threads <- function(n) {
  if (isTRUE(thread_info()$openmp_enabled)) set_num_threads(n)
  invisible(NULL)
}
