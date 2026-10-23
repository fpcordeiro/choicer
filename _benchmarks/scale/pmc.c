/* Instructions retired and cycles of this process so far (all threads),
   from the kernel's per-process counters (Apple silicon, macOS 12+; no root
   needed). bench_pmc.R compiles it into output/pmc.so on first use. */
#include <libproc.h>
#include <sys/resource.h>
#include <unistd.h>
#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>

SEXP pmc_read(void) {
  struct rusage_info_v4 ri;
  SEXP out = PROTECT(Rf_allocVector(REALSXP, 2));
  if (proc_pid_rusage(getpid(), RUSAGE_INFO_V4, (rusage_info_t *) &ri) == 0) {
    REAL(out)[0] = (double) ri.ri_instructions;
    REAL(out)[1] = (double) ri.ri_cycles;
  } else {
    REAL(out)[0] = NA_REAL;
    REAL(out)[1] = NA_REAL;
  }
  UNPROTECT(1);
  return out;
}
